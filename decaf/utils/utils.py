"""
Utility functions for DECAF.
"""

import os
import logging
import sys
from typing import List, Dict, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def setup_logger(level: int = logging.INFO) -> logging.Logger:
    """
    Set up the logger for the application.

    Args:
        level: Logging level

    Returns:
        Configured logger
    """
    logger = logging.getLogger("decaf")
    logger.setLevel(level)

    # Create console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)

    # Create formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console_handler.setFormatter(formatter)

    # Add handler to logger
    logger.addHandler(console_handler)

    return logger


def read_sequences(input_file: str) -> List[SeqRecord]:
    """
    Read sequences from a FASTQ or FASTA file.

    Args:
        input_file: Path to the input file

    Returns:
        List of SeqRecord objects

    Raises:
        ValueError: If the file format is not supported
    """
    file_ext = os.path.splitext(input_file)[1].lower()

    if file_ext in [".fastq", ".fq"]:
        format_type = "fastq"
    elif file_ext in [".fasta", ".fa", ".fna"]:
        format_type = "fasta"
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    with open(input_file, "r") as handle:
        sequences = list(SeqIO.parse(handle, format_type))

    return sequences


def write_sequences(sequences: List[SeqRecord], output_file: str) -> None:
    """
    Write sequences to a FASTQ or FASTA file.

    Args:
        sequences: List of SeqRecord objects
        output_file: Path to the output file

    Raises:
        ValueError: If the file format is not supported
    """
    file_ext = os.path.splitext(output_file)[1].lower()

    if file_ext in [".fastq", ".fq"]:
        format_type = "fastq"
    elif file_ext in [".fasta", ".fa", ".fna"]:
        format_type = "fasta"
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    with open(output_file, "w") as handle:
        SeqIO.write(sequences, handle, format_type)


def get_sequence_stats(sequences: List[SeqRecord]) -> Dict[str, Any]:
    """
    Calculate statistics for a list of sequences.

    Args:
        sequences: List of SeqRecord objects

    Returns:
        Dictionary with sequence statistics
    """
    if not sequences:
        return {"count": 0, "min_length": 0, "max_length": 0, "avg_length": 0}

    lengths = [len(seq) for seq in sequences]

    return {
        "count": len(sequences),
        "min_length": min(lengths),
        "max_length": max(lengths),
        "avg_length": sum(lengths) / len(lengths),
    }


def convert_file_format(input_file: str, output_file: str) -> None:
    """
    Convert between FASTQ and FASTA formats.

    Args:
        input_file: Path to the input file
        output_file: Path to the output file

    Raises:
        ValueError: If the file format is not supported
    """
    in_ext = os.path.splitext(input_file)[1].lower()
    out_ext = os.path.splitext(output_file)[1].lower()

    # Determine input format
    if in_ext in [".fastq", ".fq"]:
        in_format = "fastq"
    elif in_ext in [".fasta", ".fa", ".fna"]:
        in_format = "fasta"
    else:
        raise ValueError(f"Unsupported input file format: {in_ext}")

    # Determine output format
    if out_ext in [".fastq", ".fq"]:
        out_format = "fastq"
    elif out_ext in [".fasta", ".fa", ".fna"]:
        out_format = "fasta"
    else:
        raise ValueError(f"Unsupported output file format: {out_ext}")

    # Convert
    with open(input_file, "r") as in_handle:
        sequences = list(SeqIO.parse(in_handle, in_format))

    with open(output_file, "w") as out_handle:
        SeqIO.write(sequences, out_handle, out_format)
