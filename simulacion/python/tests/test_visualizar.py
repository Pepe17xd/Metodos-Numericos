"""Pruebas de lectura y validación de exportaciones MATLAB."""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import visualizar


class ResultsReaderTests(unittest.TestCase):
    """Comprueba el contrato de lectura sin ejecutar el modelo ni graficar."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        for scenario in visualizar.SCENARIOS:
            (self.root / scenario).mkdir()

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_reads_csv_declared_by_manifest(self) -> None:
        """Carga una tabla declarada y conserva sus valores exportados."""
        folder = self.root / "escenario_base"
        pd.DataFrame({"x_m": [0.0], "y_m": [0.0], "T_C": [25.0]}).to_csv(
            folder / "campo.csv", index=False
        )
        (folder / "manifest.json").write_text(
            json.dumps({"datasets": {"field": "campo.csv"}}), encoding="utf-8"
        )

        results = visualizar.Results(self.root)
        field = results.table("escenario_base", "field")

        self.assertIsNotNone(field)
        self.assertEqual(list(field.columns), ["x_m", "y_m", "T_C"])
        self.assertEqual(float(field.loc[0, "T_C"]), 25.0)
        self.assertFalse(results.warnings)

    def test_rejects_manifest_path_outside_scenario(self) -> None:
        """Evita que un manifiesto lea archivos ajenos a su escenario."""
        folder = self.root / "escenario_base"
        (folder / "manifest.json").write_text(
            json.dumps({"datasets": {"field": "../fuera.csv"}}), encoding="utf-8"
        )

        with self.assertRaises(visualizar.InputError):
            visualizar.Results(self.root)

    def test_reports_incomplete_field_schema(self) -> None:
        """Informa cuando una tabla no permite identificar x, y y temperatura."""
        folder = self.root / "escenario_base"
        pd.DataFrame({"x_m": [0.0], "valor": [25.0]}).to_csv(
            folder / "campo.csv", index=False
        )
        (folder / "manifest.json").write_text(
            json.dumps({"datasets": {"field": "campo.csv"}}), encoding="utf-8"
        )
        results = visualizar.Results(self.root)

        self.assertIsNone(visualizar.field_data(results))
        self.assertTrue(any("field_grid" in warning for warning in results.warnings))


if __name__ == "__main__":
    unittest.main()
