"""
Utility functions for DECAF.
"""

import logging
import os
import sys
from typing import Any, Dict, List, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)


def setup_logger(level: int = logging.INFO) -> None:
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


def process_predictions(
    predictions: Union[List[Dict[str, Any]], List[List[Dict[str, Any]]], None],
    model_config: Dict[str, Any],
    output_folder: str,
    input_file: str,
    threshold: float = 0.5,
) -> None:
    """
    Process model predictions and write sequences to classified files.

    Args:
        predictions: List of prediction batches from predict_step
        model_config: Model configuration dictionary
        output_folder: Output directory path
        input_file: Path to the input file (for determining format)
        threshold: Confidence threshold for classification (default: 0.5)
    """
    processed_predictions: List[Dict[str, Any]] = []

    if predictions is not None:
        # Flatten predictions if they're nested
        for batch in predictions:
            if isinstance(batch, list):
                processed_predictions.extend(batch)
            else:
                processed_predictions.append(batch)

    format_type = determine_file_format(input_file)
    categories = get_categories_from_config(model_config)
    classified_seqs = classify_sequences(
        processed_predictions, categories, format_type, threshold
    )
    write_classified_sequences(classified_seqs, output_folder, input_file, format_type)


def determine_file_format(input_file: str) -> str:
    """
    Determine the file format based on the input file extension.

    Args:
        input_file: Path to the input file

    Returns:
        The file format ('fasta' or 'fastq')

    Raises:
        ValueError: If the file format is not supported
    """
    file_ext = os.path.splitext(input_file)[1].lower()
    if file_ext in [".fastq", ".fq"]:
        return "fastq"
    elif file_ext in [".fasta", ".fa", ".fna"]:
        return "fasta"
    raise ValueError(f"Unsupported file format: {file_ext}")


def get_categories_from_config(model_config: Dict[str, Any]) -> Dict[str, str]:
    """
    Extract category mapping from model configuration.

    Args:
        model_config: Model configuration dictionary

    Returns:
        Dictionary mapping prediction IDs to category names
    """
    return model_config.get("categories", {})


def classify_sequences(
    predictions: List[Dict[str, Any]],
    categories: Dict[str, str],
    format_type: str,
    threshold: float = 0.5,  # Add threshold parameter
) -> Dict[str, List[SeqRecord]]:
    """
    Classify sequences based on model predictions.

    Args:
        predictions: List of prediction batches
        categories: Category mapping
        format_type: Input file format
        threshold: Confidence threshold for classification

    Returns:
        Dictionary mapping categories to classified sequences
    """
    # Initialize with all categories plus "unclassified"
    classified_seqs: Dict[str, List[SeqRecord]] = {
        category: [] for category in categories.values()
    }
    classified_seqs["unclassified"] = []  # Add "unclassified" category

    for batch in predictions:
        process_prediction_batch(
            batch, categories, format_type, classified_seqs, threshold
        )

    return classified_seqs


def process_prediction_batch(
    batch: Dict[str, Any],
    categories: Dict[str, str],
    format_type: str,
    classified_seqs: Dict[str, List[SeqRecord]],
    threshold: float = 0.5,  # Add threshold parameter
) -> None:
    """
    Process a single batch of predictions and classify sequences.

    Args:
        batch: Prediction batch
        categories: Category mapping
        format_type: Input file format
        classified_seqs: Dictionary to store classified sequences
        threshold: Confidence threshold for classification
    """
    sequence_ids = batch["sequence_ids"]
    sequence_strs = batch["sequence_strs"]
    batch_preds = batch["predictions"].cpu().numpy()
    batch_confs = batch["confidences"].cpu().numpy()

    for seq_id, seq_str, pred, conf in zip(
        sequence_ids, sequence_strs, batch_preds, batch_confs
    ):
        # Check if confidence is below threshold
        if conf < threshold:
            category = "unclassified"
        else:
            category = categories.get(str(pred.item()), "unknown")

        classified_seq = create_classified_sequence(
            seq_str,
            seq_id,
            category,
            conf,
            format_type,
        )
        classified_seqs[category].append(classified_seq)


def create_classified_sequence(
    original_seq: str,
    original_id: str,
    category: str,
    confidence: float,
    format_type: str,
) -> SeqRecord:
    """
    Create a new SeqRecord with classification metadata.

    Args:
        original_seq: Original sequence string
        original_id: Original suequence id
        category: Predicted category
        confidence: Prediction confidence
        format_type: Input file format

    Returns:
        New SeqRecord with classification info
    """
    new_seq = SeqRecord(
        seq=Seq(original_seq),
        id=original_id,
        description=f" DECAF_pred={category} DECAF_confidence={confidence:.4f}",
    )

    if format_type == "fastq" and hasattr(original_seq, "letter_annotations"):
        new_seq.letter_annotations = original_seq.letter_annotations

    return new_seq


def write_classified_sequences(
    classified_seqs: Dict[str, List[SeqRecord]],
    output_folder: str,
    input_file: str,
    format_type: str,
) -> None:
    """
    Write classified sequences to separate files in single-line FASTA format.

    Args:
        classified_seqs: Dictionary of classified sequences
        output_folder: Output directory path
        input_file: Original input file path (for extension)
        format_type: File format ('fasta' or 'fastq')
    """
    file_ext = os.path.splitext(input_file)[1].lower()

    for category, seqs in classified_seqs.items():
        if not seqs:
            continue

        output_file = os.path.join(output_folder, f"classified_{category}{file_ext}")

        if format_type == "fasta":
            with open(output_file, "w") as handle:
                for record in seqs:
                    handle.write(f">{record.id} {record.description}\n")
                    handle.write(f"{str(record.seq)}\n")
        else:
            with open(output_file, "w") as handle:
                SeqIO.write(seqs, handle, format_type)

        logger.info(f"Wrote {len(seqs)} sequences to {output_file}")
