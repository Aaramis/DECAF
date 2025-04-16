import os
import sys
import tempfile
from decaf.cli import main


def test_entry_point(monkeypatch):

    with tempfile.NamedTemporaryFile(suffix=".fa", delete=False) as tmp_input:
        tmp_input.write(b">seq1\nAGCTAGCTAG\n")
        tmp_input_path = tmp_input.name

    with tempfile.TemporaryDirectory() as tmp_output_dir:
        monkeypatch.setattr(
            sys,
            "argv",
            [
                "decaf",
                "--barcode",
                "ITS",
                "--taxa",
                "plants",
                "--input_fastq",
                tmp_input_path,
                "--output_folder",
                tmp_output_dir,
            ],
        )

        try:
            main()
        except SystemExit as e:
            assert e.code == 0

    os.unlink(tmp_input_path)
