# Installation of DECAF

## Installation via PyPI

The simplest way to install DECAF is via PyPI :

```bash
pip install decaf
```

## Installation from Source

### Prerequisites

Before installing DECAF, make sure you have the following dependencies :

- Python 3.8 or higher
- Git
- A NVIDIA graphics card (recommended for fast processing)
- pip (Python package manager)

### Installation Steps

1. **Clone the Repository**
```bash
git clone https://github.com/Aaramis/DECAF.git
cd DECAF
```

2. **Create a Virtual Environment**

It is strongly recommended to use a virtual environment to manage dependencies.

On Linux/Mac :
```bash
python -m venv decaf-env
source decaf-env/bin/activate
```

Sur Windows :
```bash
python -m venv decaf-env
decaf-env\Scripts\activate
```

3. **Install Dependencies**

Install the required dependencies :
```bash
pip install -r requirements.txt
```

For developers, install the development dependencies :
```bash
pip install -r requirements-dev.txt
```

4. **Install DECAF in Development Mode**
```bash
pip install -e .
```

## Verification of Installation

To verify that DECAF is correctly installed, run :
```bash
decaf --version
```

## Environment Configuration

### Environment Variables

DECAF uses the following environment variables :

- `DECAF_CONFIG_PATH` : Path to the configuration file
- `DECAF_LOG_LEVEL` : Log level (DEBUG, INFO, WARNING, ERROR)
- `DECAF_GPU` : Use GPU (0 to disable)

### Log Configuration

The log configuration can be customized in the `config/logging.yaml` file.

## Troubleshooting

### Common Issues

1. **Python Version Error**
   - Ensure you are using Python 3.8 or higher
   - Check the version with : `python --version`

2. **Permission Issues**
   - Use a virtual environment
   - Run commands with appropriate permissions

3. **Missing Dependencies**
   - Verify that all dependencies are installed
   - Reinstall dependencies if necessary

## Update

To update DECAF :

```bash
pip install --upgrade decaf
```

Or from the source code :
```bash
git pull
pip install --upgrade -e .
```

## Uninstall

To uninstall DECAF :
```bash
pip uninstall decaf
```

## Support

For any questions or issues, please :

1. Open an issue on GitHub
2. Consult the FAQ
3. Contact the development team
