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