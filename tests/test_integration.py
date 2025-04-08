import unittest
import os
import tempfile
from unittest.mock import MagicMock, patch
import torch
import json
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Import modules
from decaf.data import AmpliconDataset
from decaf.models import AmpliconClassifier, load_model
from decaf.utils import read_sequences, write_sequences
from decaf.config import load_config, get_model_config


class TestDECAFIntegration(unittest.TestCase):
    def setUp(self):
        # Create temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()

        # Create mock sequences
        self.test_sequences = [
            SeqRecord(
                Seq("ATGCTAGCTAGCTAGCT"), id="seq1", description="test sequence 1"
            ),
            SeqRecord(
                Seq("GCATCGATCGATCGATCG"), id="seq2", description="test sequence 2"
            ),
            SeqRecord(
                Seq("TACGTACGATCGATCGAT"), id="seq3", description="test sequence 3"
            ),
        ]

        # Create test FASTA file
        self.test_fasta = os.path.join(self.temp_dir.name, "test_sequences.fasta")
        write_sequences(self.test_sequences, self.test_fasta)

        # Create mock config file
        self.config_file = os.path.join(self.temp_dir.name, "config.json")
        self.config = {
            "models": {
                "ITS": {
                    "plants": {
                        "model_name": "ITS_Plant",
                        "barcode": "ITS",
                        "model_path": "mock_path",
                        "categories": ["target", "contaminant"],
                        "preprocessing": {"max_length": 128},
                    }
                }
            }
        }

        with open(self.config_file, "w") as f:
            json.dump(self.config, f)

    def tearDown(self):
        self.temp_dir.cleanup()

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_end_to_end_flow(
        self, mock_model_from_pretrained, mock_tokenizer_from_pretrained
    ):
        """Test end-to-end flow from loading data to classification"""
        # Set up mocks
        mock_tokenizer = MagicMock()
        mock_tokenizer.return_value = {
            "input_ids": torch.tensor([[1, 2, 3, 4, 0, 0]]),
            "attention_mask": torch.tensor([[1, 1, 1, 1, 0, 0]]),
        }
        mock_tokenizer_from_pretrained.return_value = mock_tokenizer

        mock_model = MagicMock()
        mock_model.classifier = MagicMock()
        mock_model.classifier.out_features = 2

        # Mock model outputs
        mock_outputs = MagicMock()
        mock_outputs.logits = torch.tensor([[0.8, 0.2], [0.3, 0.7]])
        mock_model.return_value = mock_outputs

        mock_model_from_pretrained.return_value = mock_model

        # Load configuration
        configs = load_config()
        self.assertIsNotNone(configs)

        # Get model config
        model_config = get_model_config(configs, "ITS", "plants")
        self.assertIsNotNone(model_config)
        self.assertEqual(model_config["model_name"], "ITS_Plant")

        # Load model
        classifier = load_model(model_config)
        self.assertIsInstance(classifier, AmpliconClassifier)

        # Read sequences
        sequences = read_sequences(self.test_fasta)
        self.assertEqual(len(sequences), 3)

        # Create dataset
        dataset = AmpliconDataset(
            sequences=sequences,
            tokenizer=classifier.tokenizer,
            config=model_config.get("preprocessing", {}),
        )
        self.assertEqual(len(dataset), 3)

        # Test processing batch
        batch_size = 2
        dataloader = torch.utils.data.DataLoader(dataset, batch_size=batch_size)

        # Process first batch
        batch = next(iter(dataloader))
        self.assertIn("input_ids", batch)
        self.assertIn("attention_mask", batch)

        # Ensure model can process the batch
        with torch.no_grad():
            inputs = {
                "input_ids": batch["input_ids"],
                "attention_mask": batch["attention_mask"],
            }
            outputs = classifier(**inputs)

            # Check outputs format
            logits = outputs.logits
            self.assertEqual(logits.shape[0], batch_size)  # Batch dimension
            self.assertEqual(
                logits.shape[1], len(classifier.categories)
            )  # Number of classes

            # Compute predictions
            probs = torch.softmax(logits, dim=-1)
            preds = torch.argmax(probs, dim=-1)

            # Ensure predictions have the right shape
            self.assertEqual(preds.shape[0], batch_size)


if __name__ == "__main__":
    unittest.main()
