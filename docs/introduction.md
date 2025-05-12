# Introduction to DECAF

DECAF (DEcontamination and Classification of Amplicon Fragment) is a bioinformatics framework designed for the analysis and decontamination of environmental DNA sequences. It uses deep learning models to enhance the reliability of genomic analyses.

## Objectives

- **Sequence Classification**
  - ITS sequence identification
  - Support for FASTQ and FASTA formats

- **Contaminant Management**
  - Detection of non-target sequences
  - Automatic contaminant filtering
  - Generation of analysis reports

- **Performance and Scalability**
  - Optimized for batch processing
  - GPU support via PyTorch
  - Modular extensible architecture

## Architecture

DECAF is built around several key components :

```
DECAF/
├── models/           # Implementation of deep learning models
├── data/            # Data management and preprocessing
├── utils/           # Utility functions
├── tests/           # Unit and integration tests
└── docs/           # Documentation
```

## Technologies Used

- **Deep Learning Framework**
  - PyTorch
  - PyTorch Lightning
  - Transformers

- **Data Management**
  - Polars
  - Pandas
  - Biopython

- **Tests and Quality**
  - pytest
  - pytest-cov
  - black
  - flake8

## Use Cases

DECAF is particularly useful for :

1. **Environmental Research**
   - Analysis of environmental DNA sequences
   - Biodiversity studies
   - Ecological monitoring

2. **Bioinformatics**
   - Processing of large datasets
   - Taxonomic classification
   - Sequence decontamination

3. **Academic Research**
   - Validation of experimental results
   - Comparative analysis
   - Phylogenetic studies

## Advantages

- **Precision**
  - Optimized deep learning models
  - Robust architecture
  - Continuous validation

- **Flexibility**
  - Support of multiple formats
  - Extensibility
  - Customization

- **Performance**
  - GPU support via PyTorch
  - Batch processing
  - Memory optimization

## Next Steps

To get started with DECAF, see the [Installation](installation.md) section.

For more information on specific features, see the [API Documentation](api/index.md).
