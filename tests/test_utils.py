import logging
import os
import tempfile
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from decaf.utils import (
    get_sequence_stats,
    read_sequences,
    setup_logger,
    write_sequences,
)


class TestUtils(unittest.TestCase):
    def setUp(self):
        # Create temporary files for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.fasta_file = os.path.join(self.temp_dir.name, "test.fasta")
        self.fastq_file = os.path.join(self.temp_dir.name, "test.fastq")

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

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_setup_logger(self):
        """Test logger setup"""
        setup_logger(logging.DEBUG)
        logger = logging.getLogger("decaf")
        self.assertEqual(logger.level, logging.DEBUG)
        self.assertGreater(len(logger.handlers), 0)
        self.assertIsInstance(logger.handlers[0], logging.StreamHandler)

    def test_read_sequences(self):
        """Test reading sequences from file"""
        fasta_sequences = read_sequences(self.fasta_file)

        self.assertEqual(len(fasta_sequences), 2)
        self.assertEqual(str(fasta_sequences[0].seq), "ATGCTAGCTAGCT")
        self.assertEqual(fasta_sequences[0].id, "seq1")

        fastq_sequences = read_sequences(self.fastq_file)
        self.assertEqual(len(fastq_sequences), 2)

    def test_get_sequence_stats(self):
        """Test sequence statistics calculation"""
        stats = get_sequence_stats(self.test_sequences)

        self.assertEqual(stats["count"], 2)
        self.assertEqual(stats["min_length"], 13)  # "ATGCTAGCTAGCT" length
        self.assertEqual(stats["max_length"], 14)  # "GCATCGATCGATCG" length
        self.assertEqual(stats["avg_length"], 13.5)

    def test_empty_sequence_stats(self):
        """Test sequence statistics with empty list"""
        stats = get_sequence_stats([])

        self.assertEqual(stats["count"], 0)
        self.assertEqual(stats["min_length"], 0)
        self.assertEqual(stats["max_length"], 0)
        self.assertEqual(stats["avg_length"], 0)

    def test_unsupported_format(self):
        """Test error on unsupported format"""
        invalid_file = os.path.join(self.temp_dir.name, "invalid.txt")
        with open(invalid_file, "w") as f:
            f.write("Invalid content")

        with self.assertRaises(ValueError):
            read_sequences(invalid_file)

    def test_write_sequences_invalid_output(self):
        """Test writing sequences with an unsupported output file format"""
        # Prepare the sequences
        test_sequences = [
            SeqRecord(Seq("ATGCTAGCTAGCT"), id="seq1", description="test sequence 1"),
            SeqRecord(Seq("GCATCGATCGATCG"), id="seq2", description="test sequence 2"),
        ]

        # Define an invalid output file format
        invalid_output_file = os.path.join(self.temp_dir.name, "output.txt")

        # Test that it raises a ValueError for unsupported output file format
        with self.assertRaises(ValueError):
            write_sequences(test_sequences, invalid_output_file)
