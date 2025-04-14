import os
import tempfile
import unittest

import torch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from decaf.utils import (
    classify_sequences,
    create_classified_sequence,
    determine_file_format,
    get_categories_from_config,
    process_prediction_batch,
    process_predictions,
    write_classified_sequences,
    write_sequences,
)


class TestNewUtilsFunctions(unittest.TestCase):
    def setUp(self):
        # Create temporary files for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.fasta_file = os.path.join(self.temp_dir.name, "test.fasta")
        self.fastq_file = os.path.join(self.temp_dir.name, "test.fastq")
        self.txt_file = os.path.join(self.temp_dir.name, "test.txt")

        # Sample sequences for testing
        self.test_sequences = [
            SeqRecord(Seq("ATGCTAGCTAGCT"), id="seq1", description="test sequence 1"),
            SeqRecord(Seq("GCATCGATCGATCG"), id="seq2", description="test sequence 2"),
        ]

        # Quality score for fastq
        for record in self.test_sequences:
            record.letter_annotations["phred_quality"] = [40] * len(record.seq)

        # Write test files
        write_sequences(self.test_sequences, self.fasta_file)
        write_sequences(self.test_sequences, self.fastq_file)

        # Create a text file (invalid format)
        with open(self.txt_file, "w") as f:
            f.write("This is not a sequence file")

        # Sample model config for testing
        self.model_config = {"categories": {"0": "bacteria", "1": "virus", "2": "host"}}

        # Sample predictions for testing
        self.test_predictions = [
            {
                "sequence_ids": ["seq1", "seq2"],
                "sequence_strs": ["ATGCTAGCTAGCT", "GCATCGATCGATCG"],
                "predictions": torch.tensor([0, 1]),
                "confidences": torch.tensor([0.95, 0.87]),
            }
        ]

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_determine_file_format(self):
        """Test file format determination"""
        self.assertEqual(determine_file_format(self.fasta_file), "fasta")
        self.assertEqual(determine_file_format(self.fastq_file), "fastq")

        # Test with file extensions variations
        fna_file = os.path.join(self.temp_dir.name, "test.fna")
        with open(fna_file, "w") as f:
            f.write(">seq1\nATGC\n")
        self.assertEqual(determine_file_format(fna_file), "fasta")

        fq_file = os.path.join(self.temp_dir.name, "test.fq")
        with open(fq_file, "w") as f:
            f.write("@seq1\nATGC\n+\nIIII\n")
        self.assertEqual(determine_file_format(fq_file), "fastq")

        # Test unsupported format
        with self.assertRaises(ValueError):
            determine_file_format(self.txt_file)

    def test_get_categories_from_config(self):
        """Test category extraction from config"""
        categories = get_categories_from_config(self.model_config)
        self.assertEqual(categories, {"0": "bacteria", "1": "virus", "2": "host"})

        # Test with empty config
        categories = get_categories_from_config({})
        self.assertEqual(categories, {})

    def test_create_classified_sequence(self):
        """Test creating classified sequence record"""
        # Test with fasta format
        seq = create_classified_sequence(
            "ATGCTAGCTAGCT", "seq1", "bacteria", 0.95, "fasta"
        )

        self.assertEqual(seq.id, "seq1")
        self.assertEqual(str(seq.seq), "ATGCTAGCTAGCT")
        self.assertIn("DECAF_pred=bacteria", seq.description)
        self.assertIn("DECAF_confidence=0.9500", seq.description)

        # Test with fastq format using SeqRecord
        original_seq = self.test_sequences[0]
        seq = create_classified_sequence(
            original_seq.seq, "seq1", "virus", 0.87, "fastq"
        )

        self.assertEqual(seq.id, "seq1")
        self.assertEqual(str(seq.seq), "ATGCTAGCTAGCT")
        self.assertIn("DECAF_pred=virus", seq.description)
        self.assertIn("DECAF_confidence=0.8700", seq.description)

    def test_classify_sequences(self):
        """Test sequence classification directly without mocking"""
        categories = {"0": "bacteria", "1": "virus", "2": "host"}

        # Test the function with actual data and default threshold
        classified_seqs = classify_sequences(
            self.test_predictions, categories, "fasta", threshold=0.5
        )

        # Check that the function returned the expected dictionary structure
        self.assertIn("bacteria", classified_seqs)
        self.assertIn("virus", classified_seqs)
        self.assertIn("host", classified_seqs)
        self.assertIn("unclassified", classified_seqs)

        # Check that sequences were classified correctly
        self.assertEqual(len(classified_seqs["bacteria"]), 1)
        self.assertEqual(len(classified_seqs["virus"]), 1)
        self.assertEqual(len(classified_seqs["host"]), 0)
        self.assertEqual(
            len(classified_seqs["unclassified"]), 0
        )  # None below threshold

        # Check bacteria sequence
        bact_seq = classified_seqs["bacteria"][0]
        self.assertEqual(bact_seq.id, "seq1")
        self.assertEqual(str(bact_seq.seq), "ATGCTAGCTAGCT")

        # Check virus sequence
        virus_seq = classified_seqs["virus"][0]
        self.assertEqual(virus_seq.id, "seq2")
        self.assertEqual(str(virus_seq.seq), "GCATCGATCGATCG")

    def test_classify_sequences_with_threshold(self):
        """Test sequence classification with confidence threshold"""
        categories = {"0": "bacteria", "1": "virus", "2": "host"}

        # Test with normal threshold (should classify as usual)
        classified_seqs = classify_sequences(
            self.test_predictions, categories, "fasta", threshold=0.5
        )

        # Check that "unclassified" category exists
        self.assertIn("unclassified", classified_seqs)

        # Check that sequences were still classified (both confidences are > 0.5)
        self.assertEqual(len(classified_seqs["bacteria"]), 1)
        self.assertEqual(len(classified_seqs["virus"]), 1)
        self.assertEqual(len(classified_seqs["unclassified"]), 0)

        # Test with high threshold (should classify both as unclassified)
        classified_seqs = classify_sequences(
            self.test_predictions, categories, "fasta", threshold=0.99
        )

        self.assertEqual(len(classified_seqs["bacteria"]), 0)
        self.assertEqual(len(classified_seqs["virus"]), 0)
        self.assertEqual(len(classified_seqs["unclassified"]), 2)

        # Test with medium threshold (should classify only one)
        classified_seqs = classify_sequences(
            self.test_predictions, categories, "fasta", threshold=0.9
        )

        self.assertEqual(len(classified_seqs["bacteria"]), 1)  # 0.95 > 0.9
        self.assertEqual(len(classified_seqs["virus"]), 0)  # 0.87 < 0.9
        self.assertEqual(len(classified_seqs["unclassified"]), 1)

    def test_process_prediction_batch(self):
        """Test processing prediction batch"""
        categories = {"0": "bacteria", "1": "virus", "2": "host"}
        classified_seqs = {category: [] for category in categories.values()}

        batch = {
            "sequence_ids": ["seq1", "seq2"],
            "sequence_strs": ["ATGCTAGCTAGCT", "GCATCGATCGATCG"],
            "predictions": torch.tensor([0, 1]),
            "confidences": torch.tensor([0.95, 0.87]),
        }

        process_prediction_batch(batch, categories, "fasta", classified_seqs)

        # Check that sequences were classified correctly
        self.assertEqual(len(classified_seqs["bacteria"]), 1)
        self.assertEqual(len(classified_seqs["virus"]), 1)
        self.assertEqual(len(classified_seqs["host"]), 0)

        # Check bacteria sequence
        bact_seq = classified_seqs["bacteria"][0]
        self.assertEqual(bact_seq.id, "seq1")
        self.assertEqual(str(bact_seq.seq), "ATGCTAGCTAGCT")

        # Check virus sequence
        virus_seq = classified_seqs["virus"][0]
        self.assertEqual(virus_seq.id, "seq2")
        self.assertEqual(str(virus_seq.seq), "GCATCGATCGATCG")

        # Test with unknown category
        batch = {
            "sequence_ids": ["seq3"],
            "sequence_strs": ["ATGCATGC"],
            "predictions": torch.tensor([3]),  # Category 3 is not in our mapping
            "confidences": torch.tensor([0.75]),
        }

        # Add "unknown" category to our classified_seqs dict
        classified_seqs["unknown"] = []

        process_prediction_batch(batch, categories, "fasta", classified_seqs)

        # Check that unknown sequence was classified correctly
        self.assertEqual(len(classified_seqs["unknown"]), 1)
        unknown_seq = classified_seqs["unknown"][0]
        self.assertEqual(unknown_seq.id, "seq3")

    def test_process_prediction_batch_with_threshold(self):
        """Test processing prediction batch with threshold"""
        categories = {"0": "bacteria", "1": "virus", "2": "host"}
        classified_seqs = {category: [] for category in categories.values()}
        classified_seqs["unclassified"] = []  # Add unclassified category

        batch = {
            "sequence_ids": ["seq1", "seq2", "seq3"],
            "sequence_strs": ["ATGCTAGCTAGCT", "GCATCGATCGATCG", "ATATAGCGCGCTAT"],
            "predictions": torch.tensor([0, 1, 2]),
            "confidences": torch.tensor([0.95, 0.6, 0.3]),
        }

        # Test with threshold 0.7
        process_prediction_batch(
            batch, categories, "fasta", classified_seqs, threshold=0.7
        )

        # Check that sequences were classified correctly with threshold
        self.assertEqual(len(classified_seqs["bacteria"]), 1)  # 0.95 > 0.7
        self.assertEqual(len(classified_seqs["virus"]), 0)  # 0.6 < 0.7
        self.assertEqual(len(classified_seqs["host"]), 0)  # 0.3 < 0.7
        self.assertEqual(
            len(classified_seqs["unclassified"]), 2
        )  # 2 sequences below threshold

    def test_process_predictions(self):
        """Test the process_predictions function directly without mocking"""
        # Create output folder
        output_folder = os.path.join(self.temp_dir.name, "output")
        os.makedirs(output_folder, exist_ok=True)

        # Call the function with real data
        process_predictions(
            self.test_predictions, self.model_config, output_folder, self.fasta_file
        )

        # Check that output files were created
        bacteria_file = os.path.join(output_folder, "classified_bacteria.fasta")
        virus_file = os.path.join(output_folder, "classified_virus.fasta")

        self.assertTrue(
            os.path.exists(bacteria_file), f"File {bacteria_file} does not exist"
        )
        self.assertTrue(os.path.exists(virus_file), f"File {virus_file} does not exist")

        # Check file contents
        with open(bacteria_file, "r") as f:
            content = f.read()
            self.assertIn(">seq1", content)
            self.assertIn("ATGCTAGCTAGCT", content)

        with open(virus_file, "r") as f:
            content = f.read()
            self.assertIn(">seq2", content)
            self.assertIn("GCATCGATCGATCG", content)

    def test_process_predictions_with_threshold(self):
        """Test the process_predictions function with threshold"""
        # Create output folder
        output_folder = os.path.join(self.temp_dir.name, "output_threshold")
        os.makedirs(output_folder, exist_ok=True)

        # Add a low confidence prediction
        test_predictions_with_low_conf = [
            {
                "sequence_ids": ["seq1", "seq2", "seq3"],
                "sequence_strs": ["ATGCTAGCTAGCT", "GCATCGATCGATCG", "ATATAGCGCGCTAT"],
                "predictions": torch.tensor([0, 1, 2]),
                "confidences": torch.tensor([0.95, 0.6, 0.3]),
            }
        ]

        # Call the function with threshold
        process_predictions(
            test_predictions_with_low_conf,
            self.model_config,
            output_folder,
            self.fasta_file,
            threshold=0.7,  # Set threshold to 0.7
        )

        # Check that output files were created
        bacteria_file = os.path.join(output_folder, "classified_bacteria.fasta")
        unclassified_file = os.path.join(output_folder, "classified_unclassified.fasta")

        self.assertTrue(os.path.exists(bacteria_file))
        self.assertTrue(os.path.exists(unclassified_file))

        # Check file contents
        with open(bacteria_file, "r") as f:
            content = f.read()
            self.assertIn(">seq1", content)  # Only seq1 has high enough confidence

        with open(unclassified_file, "r") as f:
            content = f.read()
            self.assertIn(">seq2", content)  # seq2 and seq3 should be unclassified
            self.assertIn(">seq3", content)

    def test_write_classified_sequences(self):
        """Test writing classified sequences to files"""
        # Create test classified sequences
        bacteria_seq = create_classified_sequence(
            "ATGCTAGCT", "seq1", "bacteria", 0.95, "fasta"
        )
        virus_seq = create_classified_sequence(
            "GCATCGATC", "seq2", "virus", 0.87, "fasta"
        )

        classified_seqs = {
            "bacteria": [bacteria_seq],
            "virus": [virus_seq],
            "host": [],  # Empty category should be skipped
        }

        # Write sequences
        output_folder = self.temp_dir.name
        write_classified_sequences(
            classified_seqs, output_folder, self.fasta_file, "fasta"
        )

        # Check that files were created
        bacteria_file = os.path.join(output_folder, "classified_bacteria.fasta")
        virus_file = os.path.join(output_folder, "classified_virus.fasta")
        host_file = os.path.join(output_folder, "classified_host.fasta")

        self.assertTrue(os.path.exists(bacteria_file))
        self.assertTrue(os.path.exists(virus_file))
        self.assertFalse(
            os.path.exists(host_file)
        )  # Should not create file for empty category

        # Check file contents
        with open(bacteria_file, "r") as f:
            content = f.read()
            self.assertIn(">seq1", content)
            self.assertIn("ATGCTAGCT", content)
            self.assertIn("DECAF_pred=bacteria", content)

        # Test with fastq format
        fastq_seq = self.test_sequences[0]
        fastq_seq.description = " DECAF_pred=bacteria DECAF_confidence=0.9500"
        classified_fastq = {"bacteria": [fastq_seq]}

        write_classified_sequences(
            classified_fastq, output_folder, self.fastq_file, "fastq"
        )
        fastq_output = os.path.join(output_folder, "classified_bacteria.fastq")
        self.assertTrue(os.path.exists(fastq_output))
