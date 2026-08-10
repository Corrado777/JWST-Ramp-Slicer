# JWST Ramp Slicer
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21875188-blue)](https://zenodo.org/badge/latestdoi/1317534449)

Author: Anthony Girmenia

Please contact agirmen@uwo.ca if you have any difficulty installing, questions about the code, or general inquiries. 
## Citing
If you use this code, please cite both the paper describing the method and the software itself:

**Paper:**
> Girmenia, A., & Metchev, S. (2026). Increasing Sensitivity to Trailed Solar System Objects in Archival JWST NIRCam Imaging with Group Differencing. Manuscript submitted for peer review.

**Software:**
> Girmenia, A. (2026). JWST-Ramp-Slicer (Version v1.0.2) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.21879108

See [`CITATION.cff`](./CITATION.cff) for full citation metadata.

## Description

This code recalibrates an $N$ group JWST NIRCam ramp into an $N-1$ sequence of group-difference images. 

Each NIRCam integration is a series of non-destructive reads (groups), each one measuring the cumulative charge since the last reset. The standard pipeline fits a slope over all the groups to estimate the count rate. While this works well for static sources, moving sources will only remain in each pixel briefly, causing the per-pixel count rate to be lower as the object trails over the detector.

Group-differencing mitigates this by only fitting the ramp between consecutive pairs of groups, minimizing the amount a moving object can trail across the detector. A proper noise aware co-addition of the group difference images can significantly increase sensitivity to moving objects, and equals the sensitivity of the original ramp-fit at the $v=0$ limiting case.

## Install
```bash
pip install -r requirements.txt
```

The `jwst` package needs CRDS reference files. The script defaults to

```bash
export CRDS_PATH=$HOME/crds_cache
export CRDS_SERVER_URL=https://jwst-crds.stsci.edu
```

and respects these variables if you have already set them. The first run downloads a few
hundred MB of reference files into `CRDS_PATH`; later runs read from that cache.

## Example Use

From an uncalibrated exposure:

```bash
python ramp_slicer.py jw02211023001_02201_00002_nrcalong_uncal.fits -o outdir
```

From a ramp file you already have (skips the detector-level stage):

```bash
python ramp_slicer.py jw02211023001_02201_00002_nrcalong_ramp.fits -o outdir
```

Or from Python:

```python
from ramp_slicer import uncal_to_ramp, slice_ramp

ramp_path = uncal_to_ramp("jw..._uncal.fits", "outdir")
cal_paths = slice_ramp(ramp_path, "outdir")
```

Add `--keep-rate` (or `keep_rate=True`) to also keep the intermediate DN/s `*_rate.fits`
files.

## Output

One file per integration `N` and group pair `K`:

```
outdir/jw02211023001_02201_00002_nrcalong_int000_grp000_cal.fits
outdir/jw02211023001_02201_00002_nrcalong_int000_grp001_cal.fits
...
```

For an exposure with `NGROUPS` groups and `NINTS` integrations that is
`(NGROUPS - 1) x NINTS` images. Each one is a standard JWST `cal` product.

## Notes
A two-point least-squares slope is exact, so running `RampFitStep` on a 2-group window is not an approximation of a difference, but exactly `(G_{k+1} - G_k) / TGROUP`

