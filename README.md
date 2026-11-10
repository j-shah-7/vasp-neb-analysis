# VASP NEB analysis

Extract image energies, relative energy profiles and atomic path coordinates
from **existing** VASP nudged elastic band calculations. The script discovers
the image folders and migrating metal, preserves energy-source provenance, and
writes fresh CSV/JSON/figure exports for each analysis.

Developed for Jameela Fatheema's metal diffusion studies on pristine and
defective two-dimensional materials, including Au, Ag and Cu on h-BN.

## Files

| File | Purpose |
|---|---|
| `extract_neb.py` | Standalone command-line extractor; copy into an individual NEB parent folder or supply `--folder`. |
| `neb_extract_and_review.ipynb` | Self-contained notebook with the same extractor, settings, table previews and interpretation. |
| `requirements.txt` | Script and plotting dependencies. |
| `requirements-notebook.txt` | Additional local Jupyter environment. |
| `tests/test_neb.py` | Synthetic profile and parser validation; no research files required. |

## Install and run

Python 3.10 or newer is recommended; validation used Python 3.12.

```powershell
python -m pip install -r requirements.txt
python extract_neb.py --folder "D:\research\my-neb"
```

Alternatively, copy `extract_neb.py` into the NEB parent containing `00`, `01`,
`02`, etc. and run the file in VS Code. With no `--folder`, it checks beside the
script first, then the current working directory. The number of images is
discovered automatically; the `IMAGES` tag, when available, provides an extra
completeness check.

For the notebook, install `requirements-notebook.txt`, select that Python
kernel, set `NEB_ROOT`, and run the cells in order. The notebook embeds its
functions and does not require the `.py` file beside it.

## Raw input

Each numbered image can supply `OUTCAR`, `OSZICAR`, `CONTCAR` and/or `POSCAR`.
Files compressed as `.gz`, `.bz2` or `.xz` are supported. Energy and geometry
from the same complete OUTCAR ionic record are preferred. Completed OSZICAR
summaries and structure files provide documented fallbacks.

Both fixed endpoints need matching energies. If their output is stored in
separate relaxation folders:

```powershell
python extract_neb.py --folder "D:\research\my-neb" --initial-outcar "D:\research\initial\OUTCAR" --final-outcar "D:\research\final\OUTCAR"
```

The endpoint outputs must represent the same structures, cell, composition and
electronic-structure settings as the NEB pathway.

## Options

| Option | Meaning |
|---|---|
| `--energy toten` | Default: electronic free energy TOTEN consistently across images. |
| `--energy sigma0` | Use the extrapolated E0 consistently instead. |
| `--metal Au` | Override automatic metal assignment. |
| `--atom-index 129` | Select a migrating atom by its **1-based** structure index. |
| `--pbc xy` | Slab default; alternatives `xyz` or `none`. |
| `--no-frequencies` | Extract the profile without exploratory local curvature fitting. |
| `--no-plot` | Write numerical results only. |
| `--output-base PATH` | Place new timestamped reports under another parent folder. |

If multiple atoms of one metal are present, a clearly dominant moving atom may
be identified. Ambiguous assignments require an explicit atom index. Periodic
path reconstruction includes neighboring translations for hexagonal cells and
removes host center-of-mass translation where suitable substrate atoms exist.

## Output

Each invocation creates `NEB_extracted/<timestamp>/`:

| File | Contents |
|---|---|
| `neb_images.csv` | TOTEN/E0, selected and relative energies, image labels, atomic/path coordinates, source filenames and diagnostics. |
| `neb_summary.json` | Composition, selected atom, completeness, sampled forward/reverse energy rises and warnings. |
| `neb_atoms_path.csv` | Aligned atomic coordinates along the path, when usable geometry exists. |
| `neb_local_frequencies.csv` | Exploratory endpoint curvature fits and review flags, when sufficient data exist. |
| `neb_profile.png/.pdf/.svg` | Energy profile. |
| `local_frequency_fits.png` | Exploratory curvature fits, when available. |

Missing image indices or selected energies produce a diagnostic export with
unavailable profile-rise fields. The extractor does not silently mix energy
conventions or invent missing image energies. Exit codes: `0` complete energy
profile, `2` diagnostic/incomplete profile, `1` input or processing error.

## Interpretation

The sampled rises are `max(E_images) - E_initial` and `max(E_images) - E_final`.
They become meaningful migration-barrier estimates only after the path and
NEB convergence have been assessed. Reading the final energy does not establish
SCF or NEB convergence. Raw atomic forces are not NEB projected-force criteria.

The optional local fit uses a cumulative mass-weighted path coordinate Q and
an endpoint model `E(Q)-E(0) ≈ K_Q Q²/2`. It reports a **one-dimensional local
frequency** with window-sensitivity diagnostics. Sparse NEB images may not
resolve a harmonic basin. These estimates are distinct from full Vineyard
attempt-frequency prefactors, which need appropriate minimum/saddle vibrational
modes [3]. A diffusion coefficient or measured switching time is not inferred
from the static NEB profile.

## Validation

```powershell
python -m unittest discover -s tests -v
```

Synthetic tests check numerical image ordering, complete-profile extraction,
consistent energy selection, missing-image behavior and rejection of a later
incomplete OUTCAR ionic record. User hard-drive calculations require local
inspection of the exported provenance and convergence diagnostics.

The companion MD workflow is at
[j-shah-7/vasp-aimd-analysis](https://github.com/j-shah-7/vasp-aimd-analysis).
Repository notebooks have cleared execution outputs. No raw calculation results
or VASP pseudopotentials are included. A redistribution license has not yet
been selected.

## References — ready to copy

1. VASP Wiki. Nudged elastic bands. https://vasp.at/wiki/Nudged_elastic_bands.
2. VASP Wiki. OSZICAR. https://vasp.at/wiki/OSZICAR.
3. Vineyard, G. H. Frequency factors and isotope effects in solid state rate processes. *Journal of Physics and Chemistry of Solids* **3**, 121–127 (1957). https://doi.org/10.1016/0022-3697(57)90059-8.
