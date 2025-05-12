# Quick Start

## Installation

```bash
pip install decaf
```

## Basic Analysis

To analyze your DNA sequences:

```bash
decaf --input_fastq data/test.fasta --output_folder output/ --taxa plants --barcode ITS --cpus 4
```

### Basic Options

| Option | Description | Example |
|--------|-------------|---------|
| `--input_fastq` | Input file (FASTA/FASTQ) | `--input_fastq sequences.fasta` |
| `--output_folder` | Output directory | `--output_folder results/` |
| `--taxa` | Taxonomic group | `--taxa plants` |
| `--barcode` | Barcode | `--barcode ITS` |
| `--batch-size` | Batch size | `--batch-size 32` |
| `--gpus` | Use GPU (0/1) | `--gpus 0` |
| `--cpus` | Number of CPU | `--cpus 1` |

## Complete Example

```bash
decaf \
  --input_fastq data/test.fasta \
  --output_folder output/ \
  --taxa plants \
  --barcode ITS \
  --batch_size 64 \
  --gpus 0 \
  --cpus 1
```
