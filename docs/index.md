# DECAF

*DECAF -  DEcontamination and Classification of Amplicon Fragment*

Welcome to the DECAF documentation, a bioinformatics tool for the classification and decontamination of ITS (Internal Transcribed Spacer) sequences from plants.

## Main Features

- Classification of ITS sequences for plants
- Detection and filtering of contaminants
- Support of FASTQ/FASTA formats
- Command-line interface
- Modular architecture for extension to other markers and taxa

## Quick Start

```bash
# Installation
pip install decaf

# Classification of sequences
decaf --input_fastq data/test.fasta --output_folder output/ --taxa plants --barcode ITS --cpus 4
```


## Documentation Structure

- [ ] Installation - How to install DECAF
- [ ] Getting Started - First steps with DECAF
- [ ] Guides - Detailed guides for specific tasks
- [ ] API Reference - Technical API documentation
- [ ] Contributing - How to contribute to the project