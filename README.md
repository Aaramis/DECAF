<!-- ![DECAF](docs/source/images/decaf_logo.png) -->

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![PyPI version](https://badge.fury.io/py/decaf.svg)](https://pypi.org/project/decaf/)
[![Documentation Status](https://readthedocs.org/projects/decaf/badge/?version=latest)](https://decaf.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# DECAF: Decontamination and Classification of Amplicon Fragments

DECAF (Decontamination and Classification of Amplicon Fragments) is a bioinformatics framework designed for the analysis and decontamination of environmental DNA sequences. It uses deep learning models to enhance the reliability of environmental genomic analyses.

## 📋 Description

DECAF provides a complete solution for:
- Classification of ITS (Internal Transcribed Spacer) DNA sequences
- Detection and filtering of contaminants
- Large-scale environmental sequence analysis
- Integration into existing bioinformatics pipelines

## 🚀 Main Features

- **Advanced Classification**
  - Deep learning models optimized for DNA
  - Support for FASTQ and FASTA formats
  - Intuitive command-line interface

- **Contaminant Management**
  - Precise detection of non-target sequences
  - Automatic contaminant filtering
  - Detailed analysis reports

- **Performance and Scalability**
  - Optimized for batch processing
  - GPU support via PyTorch
  - Extensible modular architecture

- **Complete Documentation**
  - Detailed user guide
  - Use case examples
  - API documentation

## 📦 Installation

### Prerequisites

- Python 3.8 or higher
- Git
- NVIDIA graphics card (recommended for fast processing)

### Installing via PyPI

```bash
pip install decaf
```

### Installing from source

1. Clone the repository:
```bash
git clone https://github.com/Aaramis/DECAF.git
cd DECAF
```

2. Create a virtual environment:
```bash
python -m venv decaf-env
source decaf-env/bin/activate  # On Linux/Mac
# decaf-env\Scripts\activate  # On Windows
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## 🏃 Quick Start

```bash
decaf --input_fastq data/test.fasta --output_folder output/ --taxa plants --barcode ITS --cpus 4
```

For more options and examples, consult the complete documentation.

## 📚 Documentation

The complete documentation is available at:
[https://decaf.readthedocs.io](https://decaf.readthedocs.io)

To generate the documentation locally:

1. Install development dependencies:
```bash
pip install -r requirements.txt
```

2. Start the documentation server:
```bash
mkdocs serve
```

Then open your browser at: http://127.0.0.1:8000

## 🤝 Contributing

We welcome contributions to DECAF!

1. Open an issue to report bugs or suggest features
2. Create a pull request to contribute code
3. Follow the code style guidelines

To check your code style:
```bash
pip install black
black --check .
```

To format your code:
```bash
black .
```

## 🏗️ Project Structure

```
DECAF/
├── decaf/                 # Main code
│   ├── models/           # Model implementation
│   ├── data/             # Data management
│   └── utils/            # Utility functions
├── tests/                # Unit and integration tests
├── docs/                 # Documentation
├── config/               # Configuration files
└── data/                 # Example data
```

## 📜 License

DECAF is under the MIT license. See the [LICENSE](LICENSE) file for more details.

## 🙏 Acknowledgments

- [Auguste_GARDETTE](https://github.com/Aaramis) - Lead Developer
- [Contributors](https://github.com/Aaramis/DECAF/graphs/contributors) - All contributors

## 📞 Support

For any questions or issues, please open an issue on GitHub or contact the development team.
