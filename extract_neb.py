#!/usr/bin/env python3
"""Read one existing VASP NEB directory. No VASP jobs are submitted.

Run: python extract_neb.py
Or:  python extract_neb.py --folder "D:/research/my-neb"

Requires numpy. Matplotlib is optional for the command-line program.
Numbered image directories, elements, and the moving metal are discovered.
TOTEN is selected by default; --energy sigma0 selects E0 consistently instead.
Every invocation writes a new NEB_extracted/<timestamp> directory.
Local frequency fits are exploratory 1-D path-curvature estimates, not full
Vineyard prefactors. No diffusion coefficient is inferred from an NEB profile.
"""

from pathlib import Path
from itertools import product
from datetime import datetime
import argparse
import bz2
import csv
import gzip
import json
import lzma
import math
import re
import sys

try:
    import numpy as np
except ImportError:
    raise SystemExit("Install NumPy in this Python environment: python -m pip install numpy")

VERSION = "1.0"
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][-+]?\d+)?"
NUM_RE = re.compile(NUMBER)
ITER_RE = re.compile(r"Iteration\s+(\d+)\s*\(\s*\d+\s*\)")
TOTEN_RE = re.compile(r"free\s+energy\s+TOTEN\s*=\s*(" + NUMBER + ")", re.I)
E0_RE = re.compile(r"energy\s*\(\s*sigma\s*->\s*0\s*\)\s*=\s*(" + NUMBER + ")", re.I)
WO_RE = re.compile(r"energy\s+without\s+entropy\s*=\s*(" + NUMBER + ")", re.I)
MASSES = {
    "H": 1.00794, "Li": 6.94, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.9984, "Na": 22.98977, "Mg": 24.305,
    "Al": 26.98154, "Si": 28.085, "P": 30.97376, "S": 32.06,
    "Cl": 35.45, "K": 39.0983, "Ca": 40.078, "Ti": 47.867,
    "V": 50.9415, "Cr": 51.9961, "Mn": 54.93804, "Fe": 55.845,
    "Co": 58.93319, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.92160, "Se": 78.971,
    "Zr": 91.224, "Nb": 92.90637, "Mo": 95.95, "Ru": 101.07,
    "Rh": 102.9055, "Pd": 106.42, "Ag": 107.8682, "Cd": 112.414,
    "In": 114.818, "Sn": 118.710, "Sb": 121.760, "Te": 127.60,
    "I": 126.90447, "Hf": 178.49, "Ta": 180.9479, "W": 183.84,
    "Re": 186.207, "Os": 190.23, "Ir": 192.217, "Pt": 195.084,
    "Au": 196.96657, "Hg": 200.592, "Pb": 207.2, "Bi": 208.9804,
}
HOST_ELEMENTS = {"B", "N", "C", "Mo", "W", "S", "Se", "Te"}
FREQUENCY_CONVERSION = math.sqrt(1.602176634e-19 / 1e-20 / 1.66053906660e-27) / (2 * math.pi * 1e12)


def number(value):
    return float(value.replace("D", "E").replace("d", "e"))


def open_text(path):
    path = Path(path)
    opener = {".gz": gzip.open, ".bz2": bz2.open, ".xz": lzma.open}.get(path.suffix, open)
    return opener(path, "rt", encoding="utf-8-sig", errors="replace")


def existing_file(folder, stem):
    return next((Path(folder) / name for name in
                 (stem, stem + ".gz", stem + ".bz2", stem + ".xz")
                 if (Path(folder) / name).is_file()), None)


def read_incar(path):
    if not Path(path).is_file():
        return {}
    result = {}
    with open_text(path) as handle:
        for line in handle:
            line = re.split(r"[#!]", line, maxsplit=1)[0]
            for item in line.split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    result[key.strip().upper()] = value.strip()
    return result


def discover_images(folder):
    images = sorted((p for p in Path(folder).iterdir()
                     if p.is_dir() and p.name.isdecimal()), key=lambda p: int(p.name))
    if len(images) < 3:
        raise ValueError(f"Expected at least three numbered image folders in {folder}. "
                         "Choose the NEB parent containing 00, 01, ..., final image.")
    indices = [int(p.name) for p in images]
    if len(set(indices)) != len(indices):
        raise ValueError("Duplicate numeric image labels, such as 01 and 1, are ambiguous.")
    return images


def read_oszicar(path):
    last = None
    with open_text(path) as handle:
        for line in handle:
            step = re.match(r"\s*(\d+)\s+", line)
            free = re.search(r"(?<!\w)F\s*=\s*(" + NUMBER + ")", line)
            zero = re.search(r"(?<!\w)E0\s*=\s*(" + NUMBER + ")", line)
            if step and free and zero:
                last = {"toten": number(free.group(1)), "sigma0": number(zero.group(1)),
                        "without_entropy": None, "ionic_step": int(step.group(1)),
                        "positions": None, "forces": None, "cell": None,
                        "electronic_converged": None, "record_complete": True}
    return last


def read_outcar(path):
    """Stream complete ionic records; discard a later incomplete energy block.

    Energies in incomplete electronic iterations are not allowed to overwrite
    a prior complete ionic record. Geometry from that same record is retained.
    """
    metadata = {"species": [], "counts": [], "masses": [], "warnings": []}
    titles, vrhfin, pomass = [], [], []
    lattice, natoms, ionic_step = None, None, None
    record, last_complete, last_partial = {}, None, None

    def snapshot(rec):
        return {"toten": rec.get("toten"), "sigma0": rec.get("sigma0"),
                "without_entropy": rec.get("without_entropy"),
                "ionic_step": rec.get("ionic_step"),
                "positions": rec.get("positions"), "forces": rec.get("forces"),
                "cell": rec.get("cell"),
                "electronic_converged": rec.get("electronic_converged"),
                "record_complete": (rec.get("sigma0") is not None
                                    and rec.get("positions") is not None)}

    def remember(rec):
        nonlocal last_complete, last_partial
        if rec.get("toten") is not None or rec.get("sigma0") is not None:
            last_partial = snapshot(rec)
            if rec.get("sigma0") is not None and rec.get("positions") is not None:
                last_complete = snapshot(rec)

    with open_text(path) as handle:
        stream = iter(handle)
        for line in stream:
            if "TITEL" in line and "=" in line:
                tokens = line.split("=", 1)[1].split()
                if len(tokens) > 1:
                    match = re.match(r"([A-Z][a-z]?)", tokens[1])
                    if match:
                        titles.append(match.group(1))
            if "VRHFIN" in line:
                match = re.search(r"VRHFIN\s*=\s*([A-Z][a-z]?)", line)
                if match:
                    vrhfin.append(match.group(1))
            if "ions per type" in line and "=" in line:
                metadata["counts"] = [int(v) for v in line.split("=", 1)[1].split()]
                natoms = sum(metadata["counts"])
            if "NIONS" in line and natoms is None:
                match = re.search(r"NIONS\s*=\s*(\d+)", line)
                if match:
                    natoms = int(match.group(1))
            if "POMASS" in line:
                match = re.search(r"POMASS\s*=\s*(" + NUMBER + ")", line)
                if match:
                    pomass.append(number(match.group(1)))
            if "direct lattice vectors" in line.lower():
                vectors = []
                for _ in range(3):
                    values = NUM_RE.findall(next(stream, ""))
                    if len(values) >= 3:
                        vectors.append([number(v) for v in values[:3]])
                if len(vectors) == 3:
                    lattice = np.asarray(vectors, float)
            match = ITER_RE.search(line)
            if match:
                current = int(match.group(1))
                if ionic_step is not None and current != ionic_step:
                    remember(record)
                    record = {}
                ionic_step = current
                record["ionic_step"] = current
            if "EDIFF is reached" in line:
                record["electronic_converged"] = True
            if "POSITION" in line and "TOTAL-FORCE" in line and natoms is not None:
                if record.get("positions") is not None:
                    remember(record)
                    record = {"ionic_step": ionic_step}
                values = []
                for body in stream:
                    parts = body.split()
                    if not parts or set(body.strip()) <= {"-"}:
                        if not values:
                            continue
                        break
                    try:
                        row = [number(v) for v in parts[:6]]
                    except ValueError:
                        break
                    if len(row) != 6:
                        break
                    values.append(row)
                    if len(values) == natoms:
                        break
                if len(values) == natoms:
                    data = np.asarray(values, float)
                    record["positions"], record["forces"] = data[:, :3], data[:, 3:]
                    record["cell"] = None if lattice is None else lattice.copy()
                else:
                    metadata["warnings"].append("An incomplete POSITION block was ignored.")
            match = TOTEN_RE.search(line)
            if match:
                record["toten"] = number(match.group(1))
            match = WO_RE.search(line)
            if match:
                record["without_entropy"] = number(match.group(1))
            match = E0_RE.search(line)
            if match:
                record["sigma0"] = number(match.group(1))
        remember(record)

    ntypes = len(metadata["counts"])
    candidates = titles if len(titles) >= ntypes and ntypes else vrhfin
    metadata["species"] = candidates[:ntypes] if ntypes else []
    metadata["masses"] = pomass[:ntypes] if len(pomass) >= ntypes else []
    selected = last_complete or last_partial
    if selected is not None and not selected["record_complete"]:
        metadata["warnings"].append("Only a partial energy record is available; inspect this image.")
    if (last_complete is not None and last_partial is not None
            and last_partial["ionic_step"] != last_complete["ionic_step"]):
        metadata["warnings"].append("A later incomplete ionic record was ignored.")
    return selected, metadata


def read_poscar(path, species_hint=None):
    with open_text(path) as handle:
        lines = handle.readlines()
    if len(lines) < 8:
        raise ValueError(f"Empty or incomplete structure: {path}")
    scales = np.asarray([number(v) for v in lines[1].split()], float)
    base = np.asarray([[number(v) for v in line.split()[:3]] for line in lines[2:5]])
    if base.shape != (3, 3) or abs(np.linalg.det(base)) < 1e-12:
        raise ValueError(f"Invalid lattice: {path}")
    if len(scales) == 1:
        factor = scales[0]
        if factor < 0:
            factor = (abs(factor) / abs(np.linalg.det(base))) ** (1 / 3)
        if factor <= 0:
            raise ValueError(f"Invalid scale factor: {path}")
        scale_xyz = np.repeat(factor, 3)
    elif len(scales) == 3 and np.all(scales > 0):
        scale_xyz = scales
    else:
        raise ValueError(f"Unsupported POSCAR scale factors: {path}")
    cell = base * scale_xyz[None, :]
    tokens = lines[5].split()
    if tokens and all(v.isdigit() for v in tokens):
        counts = [int(v) for v in tokens]
        species = list(species_hint or [])
        cursor = 6
    else:
        species = tokens
        counts = [int(v) for v in lines[6].split()]
        cursor = 7
    if len(species) != len(counts):
        raise ValueError(f"VASP-4 structure needs element names from an OUTCAR: {path}")
    if lines[cursor].strip().lower().startswith("s"):
        cursor += 1
    mode = lines[cursor].strip().lower()
    cursor += 1
    if not mode.startswith(("d", "c", "k")):
        raise ValueError(f"Expected Direct or Cartesian coordinates: {path}")
    natoms = sum(counts)
    coordinates = np.asarray([[number(v) for v in line.split()[:3]]
                              for line in lines[cursor:cursor + natoms]], float)
    if coordinates.shape != (natoms, 3):
        raise ValueError(f"Incomplete atom positions: {path}")
    positions = coordinates @ cell if mode.startswith("d") else coordinates * scale_xyz
    if not np.isfinite(positions).all():
        raise ValueError(f"Nonfinite atom positions: {path}")
    return {"cell": cell, "positions": positions, "species": species, "counts": counts,
            "symbols": np.repeat(species, counts), "source": str(path), "quality": "structure_file"}


def minimum_image(displacements, cell, pbc=(True, True, False)):
    """Nearest image, including neighboring translations for hexagonal cells."""
    shape = np.asarray(displacements).shape
    fractional = np.asarray(displacements).reshape(-1, 3) @ np.linalg.inv(cell)
    periodic = np.asarray(pbc, bool)
    fractional[:, periodic] -= np.rint(fractional[:, periodic])
    best = fractional @ cell
    distance2 = np.sum(best ** 2, axis=1)
    for shift in product(*[(-1, 0, 1) if flag else (0,) for flag in pbc]):
        candidate = (fractional + shift) @ cell
        trial = np.sum(candidate ** 2, axis=1)
        improve = trial < distance2
        best[improve], distance2[improve] = candidate[improve], trial[improve]
    return best.reshape(shape)


def load_image(folder, energy_override=None, fallback_species=None):
    warnings, record = [], None
    metadata = {"species": [], "counts": [], "masses": []}
    outcar = Path(energy_override).expanduser().resolve() if energy_override else existing_file(folder, "OUTCAR")
    if energy_override and not outcar.is_file():
        raise FileNotFoundError(f"Endpoint OUTCAR override does not exist: {outcar}")
    if outcar is not None:
        record, metadata = read_outcar(outcar)
        warnings.extend(metadata.get("warnings", []))
    source = str(outcar) if record is not None else ""
    if (record is None or not record.get("record_complete")) and not energy_override:
        oszicar = existing_file(folder, "OSZICAR")
        if oszicar is not None:
            ionic_record = read_oszicar(oszicar)
            if ionic_record is not None:
                record = ionic_record
                source = str(oszicar)
    structure = None
    # A complete OUTCAR record supplies geometry and energy at the same ionic step.
    if (record is not None and record.get("positions") is not None
            and record.get("cell") is not None and metadata.get("species")
            and sum(metadata.get("counts", [])) == len(record["positions"])):
        structure = {"cell": record["cell"], "positions": record["positions"],
                     "symbols": np.repeat(metadata["species"], metadata["counts"]),
                     "species": metadata["species"], "counts": metadata["counts"],
                     "source": source, "quality": "OUTCAR_same_ionic_record"}
    if structure is None:
        hint = metadata.get("species") or fallback_species
        for name in ("CONTCAR", "POSCAR"):
            path = existing_file(folder, name)
            if path is None:
                continue
            try:
                structure = read_poscar(path, hint)
                structure["quality"] = name
                break
            except (ValueError, IndexError, StopIteration) as exc:
                warnings.append(str(exc))
    if record is None:
        warnings.append("No complete ionic energy was found in OUTCAR or OSZICAR.")
    return {"label": folder.name, "index": int(folder.name), "folder": str(folder),
            "record": record, "energy_source": source, "structure": structure,
            "metadata": metadata, "warnings": warnings}


def choose_metal(symbols, aligned_positions=None, requested=None, atom_index=None):
    unique = set(str(v) for v in symbols)
    if atom_index is not None:
        index = int(atom_index) - 1
        if not 0 <= index < len(symbols):
            raise ValueError("--atom-index is 1-based and must identify an existing atom.")
        metal = str(symbols[index])
        if requested and metal != requested:
            raise ValueError("--metal and --atom-index identify different elements.")
        return metal, index, "explicit 1-based atom index"
    if requested:
        if requested not in unique:
            raise ValueError(f"Requested metal {requested} is absent; elements are {sorted(unique)}")
        candidates = [requested]
    else:
        candidates = sorted(unique & {"Au", "Ag", "Cu", "Li"})
        if not candidates:
            candidates = sorted(unique - HOST_ELEMENTS)
    if len(candidates) != 1:
        return None, None, "ambiguous elements; use --metal and, if needed, --atom-index"
    metal = candidates[0]
    indices = np.flatnonzero(np.asarray(symbols) == metal)
    if len(indices) == 1:
        return metal, int(indices[0]), "only atom of this metal"
    if aligned_positions is None:
        return metal, None, "several metal atoms; geometry is needed to identify motion"
    excursion = np.max(np.linalg.norm(aligned_positions[:, indices] - aligned_positions[0, indices], axis=2), axis=0)
    order = np.argsort(excursion)[::-1]
    if len(order) > 1 and excursion[order[0]] <= 1.3 * max(excursion[order[1]], 1e-8):
        return metal, None, "several similarly moving metal atoms; use --atom-index"
    return metal, int(indices[order[0]]), "largest excursion among metal atoms; verify selection"


def build_path(images, pbc, requested_metal=None, atom_index=None):
    warnings = []
    structures = [image["structure"] for image in images]
    available = [s for s in structures if s is not None]
    if not available:
        return {"warnings": ["No usable structures: geometry and frequency estimates unavailable."]}
    reference = available[0]
    symbols = reference["symbols"]
    composition = {str(s): int(np.sum(symbols == s)) for s in dict.fromkeys(symbols)}
    metal, selected, detection = choose_metal(symbols, requested=requested_metal, atom_index=atom_index)
    result = {"symbols": symbols, "composition": composition, "metal": metal,
              "metal_index0": selected, "metal_detection": detection, "warnings": warnings}
    if len(available) != len(images):
        warnings.append("Some images lack geometry; continuous reaction coordinates unavailable.")
        return result
    if any(not np.array_equal(s["symbols"], symbols) for s in structures):
        warnings.append("Atom counts/order differ between images; geometric analysis disabled.")
        return result
    cell = reference["cell"]
    if any(not np.allclose(s["cell"], cell, rtol=1e-6, atol=1e-5) for s in structures):
        warnings.append("Cell changes between images; fixed-cell path analysis disabled.")
        return result
    wrapped = np.asarray([s["positions"] for s in structures])
    unwrapped = np.empty_like(wrapped)
    unwrapped[0] = wrapped[0]
    for i in range(1, len(images)):
        unwrapped[i] = unwrapped[i - 1] + minimum_image(wrapped[i] - wrapped[i - 1], cell, pbc)
    host = np.asarray([str(s) in HOST_ELEMENTS for s in symbols])
    if not np.any(host):
        host = np.ones(len(symbols), bool)
    weights = np.asarray([MASSES.get(str(s), 1.0) for s in symbols])
    # Remove only a common translation, using the host center of mass.
    centers = np.average(unwrapped[:, host], axis=1, weights=weights[host])
    aligned = unwrapped - (centers - centers[0])[:, None, :]
    metal, selected, detection = choose_metal(symbols, aligned, requested_metal, atom_index)
    if metal is None or selected is None:
        warnings.append("Moving-metal identification: " + detection)
    result.update(metal=metal, metal_index0=selected, metal_detection=detection,
                  cell=cell, positions=aligned,
                  path_A=np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(aligned, axis=0).reshape(len(images) - 1, -1), axis=1))])
    masses = np.asarray([MASSES.get(str(s), np.nan) for s in symbols])
    # Actual POTCAR masses, when present, take precedence over standard masses.
    for image in images:
        meta = image["metadata"]
        if len(meta.get("masses", [])) == len(meta.get("species", [])):
            for symbol, mass in zip(meta.get("species", []), meta.get("masses", [])):
                masses[symbols == symbol] = mass
    result["masses_amu"] = masses
    if np.isfinite(masses).all():
        steps = np.diff(aligned, axis=0)
        result["q_sqrtamu_A"] = np.r_[0.0, np.cumsum(np.sqrt(np.sum(steps ** 2 * masses[None, :, None], axis=(1, 2))))]
    else:
        warnings.append("Unknown atomic masses; mass-weighted coordinate unavailable.")
    if selected is not None:
        metal_xyz = aligned[:, selected]
        result["metal_path_A"] = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(metal_xyz, axis=0), axis=1))]
        normal = np.cross(cell[0], cell[1])
        normal /= np.linalg.norm(normal)
        result["metal_height_A"] = (metal_xyz - np.mean(aligned[:, host], axis=1)) @ normal
    intermediate_poscar = [image["label"] for image in images[1:-1]
                           if image["structure"]["quality"] == "POSCAR"]
    result["frequency_geometry_usable"] = not intermediate_poscar
    if intermediate_poscar:
        warnings.append("Only input POSCAR geometry is available for intermediate images "
                        + ", ".join(intermediate_poscar) + "; frequency fitting disabled.")
    return result


def curvature_fits(energies, q, tolerance_eV=0.001, windows=(3, 4, 5)):
    """Exploratory endpoint fits E(q)-E(0)=0.5*K*q^2, with fit diagnostics.

    q is cumulative mass-weighted displacement in sqrt(amu)*Angstrom.
    Frequencies describe only this constrained path direction. All estimates
    require human review; no full transition-state prefactor is claimed.
    """
    fits = []
    energy, q = np.asarray(energies, float), np.asarray(q, float)
    for endpoint in ("initial", "final"):
        e = energy if endpoint == "initial" else energy[::-1]
        x = q - q[0] if endpoint == "initial" else q[-1] - q[::-1]
        for n in sorted(set(min(int(w), len(e)) for w in windows if w >= 3)):
            row = {"endpoint": endpoint, "n_images": n,
                   "frequency_THz": None, "curvature_eV_per_amu_A2": None,
                   "fit_rmse_eV": None, "relative_rmse": None,
                   "q_span_fraction": None, "free_fit_endpoint_gradient": None,
                   "status": "unavailable", "flags": ""}
            xx, yy = x[:n], e[:n] - e[0]
            if not np.isfinite(xx).all() or not np.isfinite(yy).all() or np.any(np.diff(xx) <= 1e-10):
                row["status"] = "missing_or_repeated_coordinate"
            elif yy[1] < -tolerance_eV:
                row["status"] = "endpoint_not_minimum_along_path"
            elif yy[1] <= tolerance_eV:
                row["status"] = "near_minimum_energy_change_unresolved"
            elif np.any(yy < -tolerance_eV):
                row["status"] = "window_extends_below_endpoint_energy"
            else:
                squared = xx ** 2
                coefficient = float(np.dot(squared, yy) / np.dot(squared, squared))
                curvature = 2 * coefficient
                row["curvature_eV_per_amu_A2"] = curvature
                if curvature > 0:
                    predicted = coefficient * squared
                    rmse = float(np.sqrt(np.mean((predicted - yy) ** 2)))
                    span = float(xx[-1] / x[-1]) if x[-1] > 0 else None
                    free = np.polynomial.polynomial.polyfit(xx, yy, 2)
                    row.update(frequency_THz=FREQUENCY_CONVERSION * math.sqrt(curvature),
                               fit_rmse_eV=rmse, relative_rmse=rmse / max(float(np.ptp(yy)), 1e-15),
                               q_span_fraction=span, free_fit_endpoint_gradient=float(free[1]),
                               status="exploratory_1d_only")
                    flags = ["one_sided_fit", "not_a_full_Vineyard_prefactor"]
                    if span is not None and span > 0.35:
                        flags.append("broad_fraction_of_path")
                    if row["relative_rmse"] > 0.10:
                        flags.append("nonparabolic_window")
                    if abs(float(free[1]) * xx[-1]) > 0.25 * max(float(np.ptp(yy)), tolerance_eV):
                        flags.append("free_fit_has_endpoint_slope")
                    row["flags"] = "; ".join(flags)
                else:
                    row["status"] = "nonpositive_curvature"
            fits.append(row)
        valid = [f for f in fits if f["endpoint"] == endpoint and f["frequency_THz"] is not None]
        if len(valid) >= 2:
            values = [f["frequency_THz"] for f in valid]
            if max(values) / min(values) > 1.25:
                for row in valid:
                    row["flags"] += "; fit_window_sensitive"
    return fits


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, Path):
        return str(value)
    return value


def write_csv(path, rows, headers=None):
    if not rows and not headers:
        return
    fields = headers or list(rows[0])
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{key: json_safe(row.get(key)) for key in fields} for row in rows])


def new_output_folder(root, output_base=None):
    base = Path(output_base).expanduser().resolve() if output_base else root / "NEB_extracted"
    base.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    for index in range(1000):
        path = base / (stamp if index == 0 else f"{stamp}_{index:03d}")
        try:
            path.mkdir()
            return path
        except FileExistsError:
            continue
    raise FileExistsError(f"Could not allocate a new output folder in {base}")


def plot_profile(result, show=False):
    import matplotlib.pyplot as plt
    rows, summary = result["rows"], result["summary"]
    path = result["path"]
    if "path_A" in path and summary["image_sequence_complete"]:
        x, xlabel = path["path_A"], "Reaction coordinate, all-atom path length (Å)"
    else:
        x, xlabel = [r["image_index"] for r in rows], "Image index"
    y = [np.nan if r["relative_to_initial_eV"] is None else r["relative_to_initial_eV"] for r in rows]
    fig, ax = plt.subplots(figsize=(8.0, 5.4))
    x_plot, y_plot = [], []
    for i, (xx, yy) in enumerate(zip(x, y)):
        if i and rows[i]["image_index"] - rows[i - 1]["image_index"] > 1:
            x_plot.append(np.nan)
            y_plot.append(np.nan)
        x_plot.append(xx)
        y_plot.append(yy)
    ax.plot(x_plot, y_plot, "o-", lw=2.5, ms=6, color="#B18416")
    if not np.isfinite(y).any():
        ax.text(0.5, 0.5, "Initial energy unavailable; inspect neb_images.csv", ha="center", va="center", transform=ax.transAxes)
    for xx, yy, row in zip(x, y, rows):
        if np.isfinite(yy):
            ax.annotate(row["image"], (xx, yy), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=10)
    ax.set(xlabel=xlabel, ylabel=f"Energy relative to initial image (eV), {summary['energy_label']}",
           title=f"{summary.get('metal') or 'Metal unassigned'} — {Path(summary['neb_folder']).name}")
    ax.tick_params(labelsize=13)
    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14)
    ax.set_title(ax.get_title(), fontsize=16)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(result["output_dir"] / f"neb_profile.{ext}", dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def plot_frequency_fits(result, show=False):
    import matplotlib.pyplot as plt
    q = result["path"].get("q_sqrtamu_A")
    if q is None or not result["frequency_fits"]:
        return None
    energy = np.asarray([r["selected_energy_eV"] for r in result["rows"]], float)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, endpoint in zip(axes, ("initial", "final")):
        xx = q - q[0] if endpoint == "initial" else q[-1] - q[::-1]
        yy = energy - energy[0] if endpoint == "initial" else energy[::-1] - energy[-1]
        ax.plot(xx, yy, "ko-", label="NEB images", lw=1.5)
        for row in result["frequency_fits"]:
            if row["endpoint"] == endpoint and row["frequency_THz"] is not None:
                grid = np.linspace(0, xx[row["n_images"] - 1], 100)
                ax.plot(grid, 0.5 * row["curvature_eV_per_amu_A2"] * grid ** 2,
                        label=f"{row['n_images']} images: {row['frequency_THz']:.3g} THz")
        ax.set(xlabel="Mass-weighted path distance (√amu Å)", ylabel="Energy above endpoint (eV)",
               title=f"{endpoint.capitalize()} endpoint: exploratory curvature")
        ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(result["output_dir"] / "local_frequency_fits.png", dpi=250, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return fig


def analyze_neb(folder, energy="toten", metal=None, atom_index=None, pbc="xy",
                output_base=None, initial_outcar=None, final_outcar=None,
                estimate_frequencies=True, energy_tolerance_eV=0.001, make_plots=True):
    root = Path(folder).expanduser().resolve()
    if not root.is_dir():
        raise FileNotFoundError(f"NEB folder does not exist: {root}")
    if energy not in ("toten", "sigma0"):
        raise ValueError("energy must be 'toten' or 'sigma0'.")
    if any(character not in "xyz" for character in pbc):
        raise ValueError("pbc must contain only x, y, z, or be an empty string.")
    if not math.isfinite(energy_tolerance_eV) or energy_tolerance_eV < 0:
        raise ValueError("energy_tolerance_eV must be nonnegative and finite.")
    folders = discover_images(root)
    incar = read_incar(root / "INCAR")
    indices = [int(p.name) for p in folders]
    sequence_ok = indices == list(range(0, max(indices) + 1))
    warnings = []
    if not sequence_ok:
        warnings.append("Image indices must start at zero and contain no gaps for a complete NEB profile.")
    if "IMAGES" in incar:
        expected = int(incar["IMAGES"].split()[0]) + 2
        if len(folders) != expected:
            sequence_ok = False
            warnings.append(f"INCAR expects {expected} folders including endpoints; found {len(folders)}.")
    images, hint = [], None
    for i, folder_path in enumerate(folders):
        override = initial_outcar if i == 0 else final_outcar if i == len(folders) - 1 else None
        image = load_image(folder_path, override, hint)
        if image["structure"] is not None:
            hint = image["structure"]["species"]
        images.append(image)
    # A later VASP-5 structure can also identify species for earlier VASP-4 images.
    if hint:
        for image in images:
            if image["structure"] is None:
                for stem in ("CONTCAR", "POSCAR"):
                    candidate = existing_file(image["folder"], stem)
                    if candidate:
                        try:
                            image["structure"] = read_poscar(candidate, hint)
                            image["structure"]["quality"] = stem
                            break
                        except (ValueError, IndexError):
                            pass
    path = build_path(images, tuple(c in pbc for c in "xyz"), metal, atom_index)
    warnings.extend(path.get("warnings", []))
    selected = [image["record"].get(energy) if image["record"] else None for image in images]
    complete = (sequence_ok and all(v is not None and math.isfinite(v) for v in selected)
                and all(image["record"] and image["record"].get("record_complete") for image in images))
    intermediate_steps = {image["record"].get("ionic_step") for image in images[1:-1]
                          if image["record"] is not None and image["record"].get("ionic_step") is not None}
    if len(intermediate_steps) > 1:
        warnings.append("Intermediate images end at different ionic iterations. This is not a confirmed synchronized NEB iteration; local frequency fitting is disabled.")
        path["frequency_geometry_usable"] = False
    initial, final = selected[0], selected[-1]
    finite = [v for v in selected if v is not None and math.isfinite(v)]
    minimum = min(finite) if finite else None
    rows = []
    for i, image in enumerate(images):
        rec = image["record"] or {}
        structure = image["structure"]
        warnings.extend(f"Image {image['label']}: {note}" for note in image["warnings"])
        value = selected[i]
        if value is None:
            warnings.append(f"Image {image['label']}: selected {energy} energy is missing; it was not replaced by a different energy type.")
        row = {"image": image["label"], "image_index": image["index"],
               "selected_energy_eV": value, "energy_type": energy,
               "TOTEN_eV": rec.get("toten"), "E0_sigma0_eV": rec.get("sigma0"),
               "energy_without_entropy_eV": rec.get("without_entropy"),
               "relative_to_initial_eV": value - initial if value is not None and initial is not None else None,
               "relative_to_final_eV": value - final if value is not None and final is not None else None,
               "relative_to_lowest_image_eV": value - minimum if value is not None and minimum is not None else None,
               "ionic_step": rec.get("ionic_step"), "complete_ionic_record": rec.get("record_complete"),
               "energy_source": image["energy_source"],
               "geometry_source": structure["source"] if structure else "",
               "geometry_quality": structure["quality"] if structure else "",
               "electronic_convergence_marker": rec.get("electronic_converged"),
               "max_raw_force_eV_per_A": float(np.max(np.linalg.norm(rec["forces"], axis=1))) if rec.get("forces") is not None else None,
               "metal": path.get("metal"),
               "metal_atom_index_1based": path.get("metal_index0") + 1 if path.get("metal_index0") is not None else None,
               "path_all_atoms_A": None, "path_metal_A": None, "q_sqrtamu_A": None,
               "metal_x_A": None, "metal_y_A": None, "metal_z_A": None, "metal_height_A": None}
        if "positions" in path:
            row["path_all_atoms_A"] = float(path["path_A"][i])
            if "q_sqrtamu_A" in path:
                row["q_sqrtamu_A"] = float(path["q_sqrtamu_A"][i])
            if path.get("metal_index0") is not None:
                xyz = path["positions"][i, path["metal_index0"]]
                row.update(metal_x_A=float(xyz[0]), metal_y_A=float(xyz[1]), metal_z_A=float(xyz[2]),
                           path_metal_A=float(path["metal_path_A"][i]), metal_height_A=float(path["metal_height_A"][i]))
        rows.append(row)
    fits = []
    if (estimate_frequencies and complete and path.get("frequency_geometry_usable")
            and "q_sqrtamu_A" in path):
        fits = curvature_fits(selected, path["q_sqrtamu_A"], energy_tolerance_eV)
    maximum_index = int(np.argmax(selected)) if complete else None
    output_dir = new_output_folder(root, output_base)
    summary = {"tool_version": VERSION, "neb_folder": str(root), "output_folder": str(output_dir),
               "image_labels": [p.name for p in folders], "number_of_folders": len(folders),
               "number_of_intermediate_images": len(folders) - 2,
               "image_sequence_complete": sequence_ok, "complete_energy_profile": complete,
               "energy_type": energy, "energy_label": "TOTEN" if energy == "toten" else "E0 (sigma→0)",
               "composition": path.get("composition"), "metal": path.get("metal"),
               "metal_atom_index_1based": path.get("metal_index0") + 1 if path.get("metal_index0") is not None else None,
               "metal_detection": path.get("metal_detection"),
               "geometry_complete": "positions" in path,
               "pbc": pbc, "translation_reference": "host center of mass; all atoms if no recognized host",
               "initial_energy_eV": initial, "final_energy_eV": final,
               "endpoint_energy_change_eV": final - initial if initial is not None and final is not None else None,
               "sampled_forward_energy_rise_eV": max(selected) - initial if complete else None,
               "sampled_reverse_energy_rise_eV": max(selected) - final if complete else None,
               "maximum_energy_image": images[maximum_index]["label"] if complete else None,
               "maximum_is_interior": 0 < maximum_index < len(images) - 1 if complete else None,
               "energy_tolerance_for_curvature_eV": energy_tolerance_eV,
               "frequency_estimates_are": "exploratory one-dimensional local frequencies, not full harmonic transition-state prefactors",
               "convergence_note": "Reading energies does not establish NEB/SCF convergence. Raw forces are not NEB projected-force convergence measures.",
               "incar_parameters": incar, "warnings": warnings}
    result = {"summary": summary, "rows": rows, "path": path, "frequency_fits": fits, "output_dir": output_dir}
    write_csv(output_dir / "neb_images.csv", rows)
    write_csv(output_dir / "neb_local_frequencies.csv", fits,
              ["endpoint", "n_images", "frequency_THz", "curvature_eV_per_amu_A2", "fit_rmse_eV",
               "relative_rmse", "q_span_fraction", "free_fit_endpoint_gradient", "status", "flags"])
    if "positions" in path:
        atom_rows = []
        for i, image in enumerate(images):
            for atom, (symbol, xyz) in enumerate(zip(path["symbols"], path["positions"][i])):
                atom_rows.append({"image": image["label"], "atom_index_1based": atom + 1,
                                  "element": str(symbol), "mass_amu": path["masses_amu"][atom],
                                  "x_A": xyz[0], "y_A": xyz[1], "z_A": xyz[2]})
        write_csv(output_dir / "neb_atoms_path.csv", atom_rows)
    with (output_dir / "neb_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(json_safe(summary), handle, indent=2, ensure_ascii=False, allow_nan=False)
    if make_plots:
        try:
            plot_profile(result)
            plot_frequency_fits(result)
        except ImportError:
            print("Plotting skipped. Optional install: python -m pip install matplotlib")
    return result


def print_summary(result):
    summary = result["summary"]
    print(f"NEB folder: {summary['neb_folder']}")
    print(f"Detected images: {', '.join(summary['image_labels'])}")
    print(f"Metal: {summary.get('metal') or 'unassigned'}; atom index (1-based): {summary.get('metal_atom_index_1based')}")
    print(f"Metal selection: {summary.get('metal_detection') or 'unavailable'}")
    print(f"Energy convention: {summary['energy_label']}")
    for row in result["rows"]:
        absolute = "MISSING" if row["selected_energy_eV"] is None else f"{row['selected_energy_eV']:.9f}"
        relative = "" if row["relative_to_initial_eV"] is None else f"{row['relative_to_initial_eV']:+.6f}"
        print(f"  {row['image']:>4}  E = {absolute:>17} eV   E-E(initial) = {relative} eV")
    if summary["complete_energy_profile"]:
        print(f"Sampled forward energy rise: {summary['sampled_forward_energy_rise_eV']:.6f} eV")
        print(f"Sampled reverse energy rise: {summary['sampled_reverse_energy_rise_eV']:.6f} eV")
        print(f"Final - initial energy: {summary['endpoint_energy_change_eV']:+.6f} eV")
        if not summary["maximum_is_interior"]:
            print("Maximum is an endpoint: an interior saddle was not resolved by these energies.")
    else:
        print("INCOMPLETE energy profile: barriers are not reported. See image/source diagnostics.")
    print("Local frequency candidates are exploratory path-curvature estimates; review fit-window flags.")
    for warning in summary["warnings"]:
        print("NOTE:", warning)
    print(f"Saved outputs: {result['output_dir']}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folder", type=Path, help="NEB parent directory; otherwise auto-detect beside this script, then current folder.")
    parser.add_argument("--energy", choices=("toten", "sigma0"), default="toten")
    parser.add_argument("--metal", help="Optional element override, e.g. Au, Ag, Cu, Li.")
    parser.add_argument("--atom-index", type=int, help="Optional migrating atom index, 1-based.")
    parser.add_argument("--pbc", default="xy", help="Periodic axes: xy (slab default), xyz, or none.")
    parser.add_argument("--output-base", type=Path, help="Parent of the new timestamped output folder.")
    parser.add_argument("--initial-outcar", type=Path, help="Matching relaxed initial-state OUTCAR if absent from image 00.")
    parser.add_argument("--final-outcar", type=Path, help="Matching relaxed final-state OUTCAR if absent from the final image.")
    parser.add_argument("--energy-tolerance", type=float, default=0.001, help="Near-minimum energy resolution for exploratory fits, eV; default 0.001.")
    parser.add_argument("--no-frequencies", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args(argv)
    folder = args.folder
    if folder is None:
        script_folder = Path(__file__).resolve().parent
        for candidate in (script_folder, Path.cwd()):
            try:
                discover_images(candidate)
                folder = candidate
                break
            except (ValueError, OSError):
                pass
    if folder is None:
        parser.error("No NEB images detected. Copy this script into the NEB parent or use --folder PATH.")
    try:
        result = analyze_neb(folder, energy=args.energy, metal=args.metal, atom_index=args.atom_index,
                             pbc="" if args.pbc == "none" else args.pbc,
                             output_base=args.output_base, initial_outcar=args.initial_outcar,
                             final_outcar=args.final_outcar, estimate_frequencies=not args.no_frequencies,
                             energy_tolerance_eV=args.energy_tolerance, make_plots=not args.no_plot)
        print_summary(result)
        return 0 if result["summary"]["complete_energy_profile"] else 2
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
