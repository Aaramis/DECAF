# DECAF: DNA Environmental Contaminant Analysis Framework

<!-- Optional Badges -->
<!-- ![GitHub issues](https://img.shields.io/github/issues/Aaramis/DECAF) -->
<!-- ![GitHub forks](https://img.shields.io/github/forks/Aaramis/DECAF) -->
<!-- ![GitHub stars](https://img.shields.io/github/stars/Aaramis/DECAF) -->
<!-- ![GitHub license](https://img.shields.io/github/license/Aaramis/DECAF) -->
<!-- ![Documentation Status](https://readthedocs.org/projects/decaf/badge/?version=latest) -->

DECAF is a bioinformatics tool designed for the classification and decontamination of plant Internal Transcribed Spacer (ITS) DNA sequences. Leveraging deep learning models, DECAF accurately identifies and filters out contaminants from environmental samples, improving the reliability of downstream analyses.

## Key Features

- **Plant ITS Sequence Classification:** Accurately classifies ITS sequences belonging to plants.
- **Contaminant Detection & Filtering:** Identifies and removes non-target sequences.
- **Supports Common Formats:** Works with FASTQ and FASTA file formats.
- **Command-Line Interface:** Easy to integrate into bioinformatics pipelines.
- **Modular Design:** Built with a modular architecture, allowing for potential future extensions to other genetic markers and taxa.

## Prerequisites

- Python 3.x

## Installation

1. **Clone the Repository (Optional):**
   ```bash
   git clone https://github.com/Aaramis/DECAF.git
   cd DECAF
   ```

2. **Create a Virtual Environment:**
   It's highly recommended to use a virtual environment to manage dependencies.

   - On Windows:
     ```bash
     python -m venv decaf-env
     decaf-env\Scripts\activate
     ```
   - On macOS and Linux:
     ```bash
     python -m venv decaf-env
     source decaf-env/bin/activate
     ```

3. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

## Quick Start

Basic usage instructions will be added here soon. This will include examples like:

```bash
# Example command structure (TO BE UPDATED)
# decaf process --input your_sequences.fastq --output cleaned_sequences.fastq --model its_plant_model
```

For detailed usage instructions, please refer to the full documentation.

## Documentation

Comprehensive documentation for DECAF is available at:
https://Aaramis.github.io/DECAF/

To build and view the documentation locally:

1. Ensure you have installed the development dependencies:
   ```bash
   pip install -r requirements.txt
   ```
   (You might want a separate requirements-dev.txt for docs and testing tools)

2. Serve the documentation:
   ```bash
   mkdocs serve
   ```
   Then open your browser to http://127.0.0.1:8000.

## Development Status

This project is currently in the initial development phase. Feedback and contributions are welcome!

## Contributing

We welcome contributions to DECAF! Please feel free to:

- Report bugs or suggest features by opening an issue.
- Contribute code by submitting a pull request.

Before submitting a pull request, please ensure your code adheres to the project's coding style. We use black for code formatting. You can check your code with:

```bash
pip install black  # If you haven't already
black --check .
```

And format it with:

```bash
black .
```

More detailed contribution guidelines will be added soon.

## Contributors

- [Your Name/GitHub Handle](Your GitHub Link)
- Feel free to add yourself here when you contribute!

## License

This project is licensed under the [Spécifiez votre licence, e.g., MIT License]. See the LICENSE file for details.

(Action: You should choose a license and add a LICENSE file to your repository. Common choices include MIT, Apache 2.0, and GPL.)