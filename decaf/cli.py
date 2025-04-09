"""
Command Line Interface for DECAF.
"""

import logging
import os
import sys

import click
import pytorch_lightning as pl
from torch.utils.data import DataLoader

from .config import get_model_config, load_config
from .data import AmpliconDataset
from .models import load_model
from .utils import process_predictions, read_sequences, setup_logger

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "-b",
    "--barcode",
    required=True,
    help="Molecular barcode type (e.g., ITS)",
)
@click.option(
    "-t",
    "--taxa",
    required=True,
    help="Target taxonomic group (e.g., plants)",
)
@click.option(
    "-i",
    "--input_fastq",
    required=True,
    type=click.Path(exists=True),
    help="Input FASTQ or FASTA file with sequences to classify",
)
@click.option(
    "-o",
    "--output_folder",
    required=True,
    type=click.Path(),
    help="Output folder for results",
)
@click.option(
    "--batch_size",
    default=32,
    type=int,
    help="Batch size for processing (default: 32)",
)
@click.option(
    "--cpus",
    default=1,
    type=int,
    help="Number of CPU cores to use (default: 1)",
)
@click.option(
    "--gpus",
    default=0,
    type=int,
    help="Number of GPUs to use (default: 0)",
)
@click.option(
    "--verbose",
    is_flag=True,
    help="Enable verbose output",
)
def main(
    barcode: str,
    taxa: str,
    input_fastq: str,
    output_folder: str,
    batch_size: int = 32,
    cpus: int = 1,
    gpus: int = 0,
    verbose: bool = False,
) -> None:
    """
    DECAF: DEcontamination and Classification of Amplicon Fragment

    Classify and decontaminate amplicon sequencing data.
    """
    # Setup logging
    log_level = logging.DEBUG if verbose else logging.INFO
    setup_logger(log_level)

    logger.info(f"Starting DECAF processing with {barcode} barcode for {taxa}")

    # Create output directory if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        logger.info(f"Created output directory: {output_folder}")

    # Load configuration
    logger.info("Loading model configuration")
    configs = load_config()
    model_config = get_model_config(configs, barcode, taxa)

    if not model_config:
        logger.error(f"No model configuration found for barcode={barcode}, taxa={taxa}")
        sys.exit(1)

    # Load model
    logger.info(f"Loading model: {model_config['model_name']}")
    model = load_model(model_config)

    # Read input sequences
    logger.info(f"Reading sequences from: {input_fastq}")
    sequences = read_sequences(input_fastq)
    logger.info(f"Read {len(sequences)} sequences")

    # Create Dataset
    dataset = AmpliconDataset(
        sequences=sequences,
        tokenizer=model.tokenizer,
        config=model_config.get("preprocessing", {}),
    )

    # Create DataLoader
    dataloader = DataLoader(
        dataset,
        batch_size=batch_size,
        num_workers=cpus,
        shuffle=False,
    )

    # Initialize trainer
    trainer = pl.Trainer(
        accelerator="gpu" if gpus > 0 else "cpu",
        devices=gpus if gpus > 0 else "auto",
        logger=False,
    )

    # Run predictions
    logger.info("Starting predictions...")
    predictions = trainer.predict(model, dataloader)

    # Process predictions and write output files
    logger.info("Processing predictions and writing output files...")
    process_predictions(
        predictions=predictions,
        model_config=model_config,
        output_folder=output_folder,
        input_file=input_fastq,
    )
    logger.info("DECAF processing completed successfully")


if __name__ == "__main__":
    main()
