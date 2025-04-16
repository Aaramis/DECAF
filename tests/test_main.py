import os
import sys
import tempfile
import torch

from decaf.cli import main
from unittest.mock import MagicMock, patch


@patch("transformers.AutoModel.from_pretrained")
@patch("transformers.AutoTokenizer.from_pretrained")
def test_entry_point(mock_tokenizer, mock_model, monkeypatch):
    mock_tokenizer_instance = MagicMock()
    mock_tokenizer_instance.return_value = {
        "input_ids": torch.tensor([[1, 2, 3, 0]]),
        "attention_mask": torch.tensor([[1, 1, 1, 0]])
    }
    mock_tokenizer.return_value = mock_tokenizer_instance

    mock_model.return_value = MagicMock()

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
