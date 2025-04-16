import os
import tempfile
import unittest
from unittest.mock import patch

from click.testing import CliRunner

from decaf.cli import main


class TestCLI(unittest.TestCase):
    """Test cases for the DECAF CLI functionality."""

    def setUp(self):
        """Set up temporary files and directories for testing."""
        # Create a temporary FASTA file
        self.temp_dir = tempfile.TemporaryDirectory()
        self.fake_fasta = os.path.join(self.temp_dir.name, "input.fa")
        with open(self.fake_fasta, "w") as f:
            f.write(">seq1\nAGCTAGCTAG\n>seq2\nTGCATGCATG\n")

        # Create a temporary output directory
        self.output_dir = os.path.join(self.temp_dir.name, "output")

    def tearDown(self):
        """Clean up temporary files and directories."""
        self.temp_dir.cleanup()

    def test_cli_exits_on_model_config_error(self):
        """Test that CLI exits when no model configuration is found."""
        runner = CliRunner()
        with patch("decaf.cli.get_model_config", return_value=None):
            result = runner.invoke(
                main,
                [
                    "--barcode",
                    "NONEXISTENT",
                    "--taxa",
                    "NONEXISTENT",
                    "--input_fastq",
                    self.fake_fasta,
                    "--output_folder",
                    self.output_dir,
                ],
            )
            self.assertEqual(result.exit_code, 1)
            self.assertIn(
                "Starting DECAF processing with NONEXISTENT barcode for NONEXISTENT",
                result.output,
            )
            self.assertIn(
                "No model configuration found for barcode=NONEXISTENT",
                result.output,
            )

    def test_cli_full_success_flow(self):
        """Test the full CLI flow with all dependencies mocked."""
        runner = CliRunner()
        mock_model_config = {"model_name": "MockModel", "preprocessing": {}}
        mock_model = unittest.mock.MagicMock()
        mock_model.tokenizer = "fake_tokenizer"

        with (
            patch("decaf.cli.load_config", return_value={"mock": "config"}),
            patch("decaf.cli.get_model_config", return_value=mock_model_config),
            patch("decaf.cli.load_model", return_value=mock_model),
            patch(
                "decaf.cli.read_sequences", return_value=["AGCTAGCTAG", "TGCATGCATG"]
            ),
            patch("decaf.cli.AmpliconDataset"),
            patch("decaf.cli.DataLoader"),
            patch("decaf.cli.pl.Trainer.predict", return_value=["pred1", "pred2"]),
            patch("decaf.cli.process_predictions") as mock_process_predictions,
        ):

            result = runner.invoke(
                main,
                [
                    "--barcode",
                    "ITS",
                    "--taxa",
                    "plants",
                    "--input_fastq",
                    self.fake_fasta,
                    "--output_folder",
                    self.output_dir,
                    "--batch_size",
                    "16",
                    "--cpus",
                    "2",
                    "--gpus",
                    "0",
                    "--threshold",
                    "0.8",
                    "--verbose",
                ],
            )

            self.assertEqual(result.exit_code, 0)
            mock_process_predictions.assert_called_once()
            self.assertIn(
                "Starting DECAF processing with ITS barcode for plants", result.output
            )
            self.assertIn("Created output directory", result.output)

    def test_output_directory_exists(self):
        os.makedirs(self.output_dir)  # Crée le dossier avant exécution

        runner = CliRunner()
        with (patch("decaf.cli.get_model_config", return_value=None),):
            result = runner.invoke(
                main,
                [
                    "--barcode",
                    "BAD",
                    "--taxa",
                    "BAD",
                    "--input_fastq",
                    self.fake_fasta,
                    "--output_folder",
                    self.output_dir,
                ],
            )
            # Vérifie que le message "Created output directory" n’apparaît pas
            self.assertNotIn("Created output directory", result.output)

    def test_verbose_logging(self):
        runner = CliRunner()
        with (patch("decaf.cli.get_model_config", return_value=None),):
            result = runner.invoke(
                main,
                [
                    "--barcode",
                    "BAD",
                    "--taxa",
                    "BAD",
                    "--input_fastq",
                    self.fake_fasta,
                    "--output_folder",
                    self.output_dir,
                    "--verbose",
                ],
            )
            # Niveau DEBUG devrait s'afficher
            self.assertIn("Starting DECAF processing", result.output)
