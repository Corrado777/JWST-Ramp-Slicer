#!/usr/bin/env python
"""
Produce group difference images from JWST NIRCam ramps.

Usage:
    python ramp_slicer.py jw..._uncal.fits -o outdir
    python ramp_slicer.py jw..._ramp.fits  -o outdir     # skips Detector1

Outputs one calibrated ``*_int{N}_grp{K}_cal.fits`` per integration and group pair.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import astropy.io.fits as fits

# CRDS must be configured before any jwst import touches a reference file.
os.environ.setdefault("CRDS_PATH", str(Path("~/crds_cache").expanduser()))
os.environ.setdefault("CRDS_SERVER_URL", "https://jwst-crds.stsci.edu")

# --------------------------------------------------------------------------- #
# Stage 1 — UNCAL to ramp                                                      #
# --------------------------------------------------------------------------- #

def uncal_to_ramp(uncal_path: str | Path, out_dir: str | Path) -> Path:
    """Run Detector1Pipeline on an UNCAL file, stopping before the ramp fit.

    Returns the path to the ``*_ramp.fits`` file (4-D SCI array in DN).
    """
    from jwst.pipeline import Detector1Pipeline

    uncal_path = Path(uncal_path).resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    det1 = Detector1Pipeline()
    det1.jump.skip = True 
    det1.ramp_fit.skip = True
    det1.save_results = True
    det1.output_dir = str(out_dir)
    det1.run(str(uncal_path))

    stem = uncal_path.stem.replace("_uncal", "")
    ramp_path = out_dir / f"{stem}_ramp.fits"
    if not ramp_path.exists():
        raise FileNotFoundError(
            f"Detector1Pipeline did not write {ramp_path.name} into {out_dir}"
        )
    return ramp_path


# --------------------------------------------------------------------------- #
# Stage 2 — ramp to group-difference images                                    #
# --------------------------------------------------------------------------- #

def slice_ramp(
    ramp_path: str | Path,
    out_dir: str | Path,
    keep_rate: bool = False,
) -> list[Path]:
    """Slice a ramp into one calibrated image per consecutive group pair.

    Returns the calibrated file paths, ordered by group pair then integration.
    """
    ramp_path = Path(ramp_path).resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(ramp_path, memmap=False) as hdul:
        sci_shape = hdul["SCI"].shape
        header = hdul[0].header.copy()
        if "INT_TIMES" not in hdul:
            raise KeyError(
                f"{ramp_path.name} has no INT_TIMES table, so the observation "
                f"time of each group-difference frame cannot be determined."
            )
        int_times = hdul["INT_TIMES"].data.copy()

    # SCI is (n_ints, n_groups, ny, nx), or (n_groups, ny, nx) when n_ints == 1.
    if len(sci_shape) == 4:
        n_ints, n_groups = sci_shape[0], sci_shape[1]
    else:
        n_ints, n_groups = 1, sci_shape[0]

    # We need the timing information to assign the correct MJD to each group-difference frame.
    tgroup = _require(header, "TGROUP", ramp_path, float)    # s between group starts
    tframe = _require(header, "TFRAME", ramp_path, float)    # s per frame read
    nframes = _require(header, "NFRAMES", ramp_path, int)    # frames averaged per group

    gain_ref, readnoise_ref = _crds_refs(ramp_path)
    stem = ramp_path.stem.replace("_ramp", "")

    cal_paths: list[Path] = []

    for k in range(n_groups - 1):
        # Standard processing runs RampFitStep once over all groups.  We run it
        # per 2-group window instead.  
        rate_models = _rampfit_pair(ramp_path, k, gain_ref, readnoise_ref)
        # =================================================================== #

        for int_idx, rate_model in enumerate(rate_models):
            if int_idx >= len(int_times):
                raise IndexError(
                    f"{ramp_path.name} has {n_ints} integrations but its "
                    f"INT_TIMES table has only {len(int_times)} rows."
                )
            int_start_mjd = float(int_times[int_idx]["int_start_MJD_UTC"])
            t_start, t_mid, t_end = _group_diff_timing(
                int_start_mjd, tframe, nframes, tgroup, k
            )
            rate_path = _write_rate(
                rate_model, stem, int_idx, k, t_start, t_mid, t_end, tgroup, out_dir
            )
            rate_model.close()

            # Image2Pipeline is unmodified
            cal_paths.append(_run_image2(rate_path, out_dir))

            if not keep_rate:
                rate_path.unlink()

    return cal_paths


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #

def _require(header, key: str, path: Path, cast):
    """Return ``header[key]``, raising if absent."""
    if key not in header:
        raise KeyError(
            f"Required keyword {key!r} is missing from {path.name}. "
            f"This is probably not a NIRCam ramp product."
        )
    return cast(header[key])


def _crds_refs(path: Path) -> tuple[str, str]:
    """Return the CRDS (gain, readnoise) reference paths matching this exposure.

    These are handed to RampFitStep as overrides purely so the n_groups - 1
    calls below do not each re-query the CRDS server
    """
    import crds
    from jwst import datamodels

    with datamodels.open(str(path)) as model:
        params = model.get_crds_parameters()

    refs = crds.getreferences(
        params, reftypes=["gain", "readnoise"], observatory="jwst"
    )
    return refs["gain"], refs["readnoise"]


def _group_diff_timing(
    int_start_mjd: float,
    tframe: float,
    nframes: int,
    tgroup: float,
    k: int,
) -> tuple[float, float, float]:
    """Return (t_start, t_mid, t_end) in MJD for the frame D_k = G_{k+1} - G_k.

    t_eff(k) = int_start + TFRAME * (NFRAMES + 1) / 2 + k * TGROUP
    """
    offset_s = tframe * (nframes + 1) / 2.0

    t_start = int_start_mjd + (offset_s + k * tgroup) / 86400.0
    t_end = int_start_mjd + (offset_s + (k + 1) * tgroup) / 86400.0
    t_mid = (t_start + t_end) / 2.0

    return t_start, t_mid, t_end


def _rampfit_pair(ramp_path: Path, k: int, gain_ref: str, readnoise_ref: str) -> list:
    """Fit groups (k, k+1) only; return one 2-D rate model per integration.

    ``firstgroup``/``lastgroup`` make the step flag every other group DO_NOT_USE.
    """
    from jwst.ramp_fitting import RampFitStep

    result = RampFitStep.call(
        str(ramp_path),
        firstgroup=k,
        lastgroup=k + 1,
        maximum_cores="1",
        override_gain=gain_ref,
        override_readnoise=readnoise_ref,
    )
    rate_model, rateints_model = result if isinstance(result, tuple) else (result, None)

    has_cube = (
        rateints_model is not None
        and getattr(rateints_model, "data", None) is not None
        and rateints_model.data.ndim == 3
    )
    if has_cube and rateints_model.data.shape[0] > 1:
        models = [
            _slice_rateints(rateints_model, i)
            for i in range(rateints_model.data.shape[0])
        ]
        rate_model.close()
        rateints_model.close()
    else:
        # NINTS == 1: the 2-D rate model already is the single integration.
        models = [rate_model]
        if rateints_model is not None:
            rateints_model.close()

    return models


def _slice_rateints(rateints_model, int_idx: int):
    """Extract one integration of a rateints cube as a standalone ImageModel."""
    from jwst import datamodels

    im = datamodels.ImageModel(data=rateints_model.data[int_idx].copy())
    im.update(rateints_model)
    for attr in ("err", "dq", "var_poisson", "var_rnoise"):
        plane = getattr(rateints_model, attr, None)
        if plane is not None and getattr(plane, "ndim", 0) == 3:
            setattr(im, attr, plane[int_idx].copy())
    return im


def _write_rate(
    rate_model,
    stem: str,
    int_idx: int,
    k: int,
    t_start: float,
    t_mid: float,
    t_end: float,
    tgroup: float,
    out_dir: Path,
) -> Path:
    """Stamp the group-difference window's timing on a rate model and save it.
    """
    rate_model.meta.exposure.start_time = t_start
    rate_model.meta.exposure.mid_time = t_mid
    rate_model.meta.exposure.end_time = t_end
    rate_model.meta.exposure.exposure_time = tgroup
    rate_model.meta.exposure.effective_exposure_time = tgroup

    out_path = out_dir / f"{stem}_int{int_idx:03d}_grp{k:03d}_rate.fits"
    rate_model.save(str(out_path))
    return out_path


def _run_image2(rate_path: Path, out_dir: Path) -> Path:
    """Run Image2Pipeline on one group-difference rate file; return the cal path.

    Unmodified stock processing
    """
    from jwst.pipeline import Image2Pipeline

    Image2Pipeline.call(str(rate_path), output_dir=str(out_dir), save_results=True)

    cal_path = out_dir / rate_path.name.replace("_rate.fits", "_cal.fits")
    if not cal_path.exists():
        raise FileNotFoundError(
            f"Image2Pipeline did not write {cal_path.name} into {out_dir}"
        )
    return cal_path


# --------------------------------------------------------------------------- #
# Command line                                                                 #
# --------------------------------------------------------------------------- #

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Produce group-difference images from a JWST NIRCam "
                    "*_uncal.fits or *_ramp.fits file.",
    )
    parser.add_argument(
        "input", type=Path, help="a *_uncal.fits or *_ramp.fits file"
    )
    parser.add_argument(
        "-o", "--out-dir", type=Path, required=True,
        help="directory to write the group-difference images into",
    )
    parser.add_argument(
        "--keep-rate", action="store_true",
        help="also keep the intermediate DN/s *_rate.fits files",
    )
    args = parser.parse_args(argv)

    name = args.input.name
    if name.endswith("_uncal.fits"):
        ramp_path = uncal_to_ramp(args.input, args.out_dir)
    elif name.endswith("_ramp.fits"):
        ramp_path = args.input
    else:
        parser.error(f"expected a *_uncal.fits or *_ramp.fits file, got {name}")

    cal_paths = slice_ramp(ramp_path, args.out_dir, keep_rate=args.keep_rate)

    print(f"\nWrote {len(cal_paths)} group-difference images to {args.out_dir}")
    for path in cal_paths:
        print(f"  {path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())