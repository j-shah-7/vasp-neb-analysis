from pathlib import Path
import importlib.util
import json
import shutil
import tempfile
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("extract_neb", ROOT / "extract_neb.py")
neb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(neb)


class NebExtractionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        (self.root / "INCAR").write_text("IMAGES = 3\n")
        self.energies = [-10.0, -9.9, -9.5, -9.8, -10.2]
        for index, energy in enumerate(self.energies):
            folder = self.root / f"{index:02d}"
            folder.mkdir()
            (folder / "OSZICAR").write_text(f" 5 F= {energy:.8f} E0= {energy+.02:.8f} d E = 0.0\n")
            (folder / "CONTCAR").write_text(
                "Synthetic NEB\n1\n10 0 0\n0 10 0\n0 0 20\nB N Au\n1 1 1\nCartesian\n"
                f"1 1 8\n2.4 1 8\n{3+index*.2} 3 10\n")

    def tearDown(self):
        self.temp.cleanup()

    def test_complete_profile_and_energy_convention(self):
        result = neb.analyze_neb(self.root, make_plots=False, estimate_frequencies=False)
        summary = result["summary"]
        self.assertTrue(summary["complete_energy_profile"])
        self.assertEqual(summary["image_labels"], ["00", "01", "02", "03", "04"])
        self.assertEqual(summary["metal"], "Au")
        self.assertAlmostEqual(summary["sampled_forward_energy_rise_eV"], .5)
        self.assertAlmostEqual(summary["sampled_reverse_energy_rise_eV"], .7)
        alternate = neb.analyze_neb(self.root, energy="sigma0", make_plots=False, estimate_frequencies=False)
        np.testing.assert_allclose([row["selected_energy_eV"] for row in alternate["rows"]],
                                   np.array(self.energies) + .02)
        self.assertTrue((result["output_dir"] / "neb_images.csv").is_file())

    def test_missing_image_does_not_invent_barrier(self):
        shutil.rmtree(self.root / "02")
        result = neb.analyze_neb(self.root, make_plots=False, estimate_frequencies=False)
        self.assertFalse(result["summary"]["complete_energy_profile"])
        self.assertIsNone(result["summary"]["sampled_forward_energy_rise_eV"])
        self.assertTrue((result["output_dir"] / "neb_summary.json").is_file())

    def test_later_incomplete_ionic_record_is_ignored(self):
        outcar = self.root / "OUTCAR_test"
        outcar.write_text(
            " ions per type = 1 1 1\n"
            " direct lattice vectors reciprocal lattice vectors\n"
            " 10 0 0 0 0 0\n 0 10 0 0 0 0\n 0 0 20 0 0 0\n"
            " Iteration 1( 1)\n POSITION TOTAL-FORCE (eV/Angst)\n ------\n"
            " 1 1 8 0 0 0\n 2.4 1 8 0 0 0\n 3 3 10 0 0 0\n"
            " free energy TOTEN = -10.000 eV\n"
            " energy without entropy = -9.990 energy(sigma->0) = -9.980\n"
            " Iteration 2( 1)\n free energy TOTEN = -99.000 eV\n")
        record, metadata = neb.read_outcar(outcar)
        self.assertEqual(record["ionic_step"], 1)
        self.assertAlmostEqual(record["toten"], -10)
        self.assertTrue(any("incomplete" in message for message in metadata["warnings"]))

    def test_notebook_cells_are_valid_and_outputs_cleared(self):
        notebook = json.loads((ROOT / "neb_extract_and_review.ipynb").read_text())
        for index, cell in enumerate(notebook["cells"]):
            if cell["cell_type"] == "code":
                compile("".join(cell["source"]), f"cell_{index}", "exec")
                self.assertEqual(cell["outputs"], [])
                self.assertIsNone(cell["execution_count"])


if __name__ == "__main__":
    unittest.main()
