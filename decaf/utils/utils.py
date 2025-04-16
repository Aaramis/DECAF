"""
Utility functions for DECAF.

This module provides utility functions for file operations, sequence handling,
logging setup, and results processing for the DECAF project.

Functions
---------
setup_logger
    Set up the logger for the application.
read_sequences
    Read biological sequences from FASTQ or FASTA files.
write_sequences
    Write biological sequences to FASTQ or FASTA files.
get_sequence_stats
    Calculate statistics for a list of sequences.
process_predictions
    Process model predictions and write sequences to classified files.
determine_file_format
    Determine the file format based on the input file extension.
get_categories_from_config
    Extract category mapping from model configuration.
classify_sequences
    Classify sequences based on model predictions.
process_prediction_batch
    Process a single batch of predictions and classify sequences.
create_classified_sequence
    Create a new SeqRecord with classification metadata.
write_classified_sequences
    Write classified sequences to separate files.
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

    Parameters
    ----------
    level : int, default=logging.INFO
        Logging level

    Examples
    --------
    >>> setup_logger(logging.DEBUG)
    >>> logger = logging.getLogger("decaf")
    >>> logger.debug("Debug message")
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

    Parameters
    ----------
    input_file : str
        Path to the input file

    Returns
    -------
    list of SeqRecord
        List of BioPython SeqRecord objects

    Raises
    ------
    ValueError
        If the file format is not supported

    Examples
    --------
    >>> sequences = read_sequences("sample.fastq")
    >>> len(sequences)
    42
    >>> sequences[0].id
    'seq1'
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

    Parameters
    ----------
    sequences : list of SeqRecord
        List of BioPython SeqRecord objects
    output_file : str
        Path to the output file

    Raises
    ------
    ValueError
        If the file format is not supported

    Examples
    --------
    >>> from Bio.Seq import Seq
    >>> from Bio.SeqRecord import SeqRecord
    >>> seq1 = SeqRecord(Seq("ACGT"), id="seq1")
    >>> seq2 = SeqRecord(Seq("TGCA"), id="seq2")
    >>> write_sequences([seq1, seq2], "output.fasta")
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

    Parameters
    ----------
    sequences : list of SeqRecord
        List of BioPython SeqRecord objects

    Returns
    -------
    dict
        Dictionary with sequence statistics including:
        - count: Total number of sequences
        - min_length: Length of shortest sequence
        - max_length: Length of longest sequence
        - avg_length: Average sequence length

    Examples
    --------
    >>> sequences = read_sequences("sample.fasta")
    >>> stats = get_sequence_stats(sequences)
    >>> stats["count"]
    100
    >>> stats["avg_length"]
    250.5
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


def process_predictions(
    predictions: Union[List[Dict[str, Any]], List[List[Dict[str, Any]]], None],
    model_config: Dict[str, Any],
    output_folder: str,
    input_file: str,
    threshold: float = 0.5,
) -> None:
    """
    Process model predictions and write sequences to classified files.

    Parameters
    ----------
    predictions : list of dict or list of list of dict or None
        List of prediction batches from predict_step
    model_config : dict
        Model configuration dictionary
    output_folder : str
        Output directory path
    input_file : str
        Path to the input file (for determining format)
    threshold : float, default=0.5
        Confidence threshold for classification

    Examples
    --------
    >>> process_predictions(predictions, model_config, "output_dir", "input.fasta", 0.7)
    """
    processed_predictions: List[Dict[str, Any]] = []

    if predictions is not None:
        # Flatten predictions if they're nested
        for batch in predictions:
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

    Parameters
    ----------
    input_file : str
        Path to the input file

    Returns
    -------
    str
        The file format ('fasta' or 'fastq')

    Raises
    ------
    ValueError
        If the file format is not supported

    Examples
    --------
    >>> format_type = determine_file_format("sequences.fasta")
    >>> format_type
    'fasta'
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

    Parameters
    ----------
    model_config : dict
        Model configuration dictionary

    Returns
    -------
    dict
        Dictionary mapping prediction IDs to category names

    Examples
    --------
    >>> config = {"categories": {"0": "fungi", "1": "plants"}}
    >>> categories = get_categories_from_config(config)
    >>> categories["0"]
    'fungi'
    """
    return model_config.get("categories", {})


def classify_sequences(
    predictions: List[Dict[str, Any]],
    categories: Dict[str, str],
    format_type: str,
    threshold: float = 0.5,
) -> Dict[str, List[SeqRecord]]:
    """
    Classify sequences based on model predictions.

    Parameters
    ----------
    predictions : list of dict
        List of prediction batches
    categories : dict
        Category mapping
    format_type : str
        Input file format
    threshold : float, default=0.5
        Confidence threshold for classification

    Returns
    -------
    dict
        Dictionary mapping categories to classified sequences

    Examples
    --------
    >>> classified = classify_sequences(predictions, categories, "fasta", 0.6)
    >>> len(classified["fungi"])
    23
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
    threshold: float = 0.5,
) -> None:
    """
    Process a single batch of predictions and classify sequences.

    Parameters
    ----------
    batch : dict
        Prediction batch
    categories : dict
        Category mapping
    format_type : str
        Input file format
    classified_seqs : dict
        Dictionary to store classified sequences
    threshold : float, default=0.5
        Confidence threshold for classification

    Examples
    --------
    >>> process_prediction_batch(batch, categories, "fasta", classified_seqs, 0.7)
    """
    sequence_ids = batch["sequence_ids"]
    sequence_strs = batch["sequence_strs"]
    seq_qualities = batch["seq_qualities"]
    batch_preds = batch["predictions"].cpu().numpy()
    batch_confs = batch["confidences"].cpu().numpy()

    for seq_id, seq_str, seq_quality, pred, conf in zip(
        sequence_ids, sequence_strs, seq_qualities, batch_preds, batch_confs
    ):
        # Check if confidence is below threshold
        if conf < threshold:
            category = "unclassified"
        else:
            category = categories.get(str(pred.item()), "unknown")

        classified_seq = create_classified_sequence(
            seq_str,
            seq_id,
            seq_quality,
            category,
            conf,
            format_type,
        )
        classified_seqs[category].append(classified_seq)


def create_classified_sequence(
    original_seq: str,
    original_id: str,
    seq_quality: str,
    category: str,
    confidence: float,
    format_type: str,
) -> SeqRecord:
    """
    Create a new SeqRecord with classification metadata.

    Parameters
    ----------
    original_seq : str
        Original sequence string
    original_id : str
        Original sequence id
    seq_quality : str
        Original sequence quality
    category : str
        Predicted category
    confidence : float
        Prediction confidence
    format_type : str
        Input file format

    Returns
    -------
    SeqRecord
        New SeqRecord with classification info

    Examples
    --------
    >>> seq = create_classified_sequence("ACGT", "seq1", "fungi", 0.95, "fasta")
    >>> seq.description
    ' DECAF_pred=fungi DECAF_confidence=0.9500'
    """
    new_seq = SeqRecord(
        seq=Seq(original_seq),
        id=original_id,
        description=f" DECAF_pred={category} DECAF_confidence={confidence:.4f}",
    )

    if format_type == "fastq":
        phred_quality = [ord(char) - 33 for char in seq_quality]
        new_seq.letter_annotations = {"phred_quality": phred_quality}

    return new_seq


def write_classified_sequences(
    classified_seqs: Dict[str, List[SeqRecord]],
    output_folder: str,
    input_file: str,
    format_type: str,
) -> None:
    """
    Write classified sequences to separate files in single-line FASTA format.

    Parameters
    ----------
    classified_seqs : dict
        Dictionary of classified sequences
    output_folder : str
        Output directory path
    input_file : str
        Original input file path (for extension)
    format_type : str
        File format ('fasta' or 'fastq')

    Examples
    --------
    >>> write_classified_sequences(classified_seqs, "results", "input.fasta", "fasta")
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
