import unittest
from unittest.mock import MagicMock
import torch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from decaf.data import AmpliconDataset


class TestAmpliconDataset(unittest.TestCase):
    def setUp(self):
        # Create mock tokenizer
        self.mock_tokenizer = MagicMock()
        self.mock_tokenizer.return_value = {
            "input_ids": torch.tensor([[1, 2, 3, 4, 0, 0]]),
            "attention_mask": torch.tensor([[1, 1, 1, 1, 0, 0]]),
        }

        # Create test sequences
        self.test_sequences = [
            SeqRecord(Seq("ATGCTAGCTAGCT"), id="seq1"),
            SeqRecord(Seq("GCATCGATCGATCG"), id="seq2"),
            SeqRecord(Seq("TACGTACGAT"), id="seq3"),
        ]

        # Config for dataset
        self.config = {"max_length": 10}

    def test_dataset_initialization(self):
        """Test dataset initialization"""
        dataset = AmpliconDataset(
            sequences=self.test_sequences,
            tokenizer=self.mock_tokenizer,
            config=self.config,
        )

        self.assertEqual(len(dataset), 3)
        self.assertEqual(dataset.max_length, 10)
        self.assertEqual(dataset.sequences, self.test_sequences)
        self.assertEqual(dataset.tokenizer, self.mock_tokenizer)

    def test_dataset_getitem(self):
        """Test dataset __getitem__ method"""
        dataset = AmpliconDataset(
            sequences=self.test_sequences,
            tokenizer=self.mock_tokenizer,
            config=self.config,
        )

        item = dataset[0]

        # Check if tokenizer was called with the right sequence
        self.mock_tokenizer.assert_called_with(
            str(self.test_sequences[0].seq),
            padding="max_length",
            truncation=True,
            max_length=10,
            return_tensors="pt",
        )

        # Check returned item
        self.assertIn("input_ids", item)
        self.assertIn("attention_mask", item)
        self.assertIn("sequence_id", item)
        self.assertEqual(item["sequence_id"], "seq1")


if __name__ == "__main__":
    unittest.main()
