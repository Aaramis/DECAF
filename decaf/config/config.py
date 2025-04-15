"""
Configuration management for DECAF.

This module provides functions to load, retrieve, and save model configurations
for different molecular barcodes and taxonomic groups.

Functions
---------
load_config
    Load model configurations from a JSON file.
get_model_config
    Get configuration for a specific barcode and taxonomic group.
save_config
    Save model configurations to a JSON file.
"""

import json
import logging
from typing import Any, Dict, List, Optional

logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config_models.json") -> List[Dict[str, Any]]:
    """
    Load the model configuration from a JSON file.

    Parameters
    ----------
    config_path : str, default="config/config_models.json"
        Path to the configuration file

    Returns
    -------
    list of dict
        List of model configurations

    Raises
    ------
    FileNotFoundError
        If the configuration file is not found
    json.JSONDecodeError
        If the configuration file is not valid JSON

    Examples
    --------
    >>> configs = load_config("config/my_config.json")
    >>> len(configs)
    3
    """
    try:
        with open(config_path, "r") as f:
            config = json.load(f)
        return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found: {config_path}")
        raise
    except json.JSONDecodeError:
        logger.error(f"Invalid JSON in configuration file: {config_path}")
        raise


def get_model_config(
    configs: List[Dict[str, Any]], barcode: str, taxa: str
) -> Optional[Dict[str, Any]]:
    """
    Get the configuration for a specific barcode and taxonomic group.

    Parameters
    ----------
    configs : list of dict
        List of model configurations
    barcode : str
        Molecular barcode type
    taxa : str
        Taxonomic group

    Returns
    -------
    dict or None
        Model configuration if found, None otherwise

    Examples
    --------
    >>> configs = load_config()
    >>> config = get_model_config(configs, "ITS", "fungi")
    >>> config["model_type"]
    'random_forest'
    """
    for config in configs:
        if (
            config["barcode"].lower() == barcode.lower()
            and config["taxa"].lower() == taxa.lower()
        ):
            return config
    return None


def save_config(
    configs: List[Dict[str, Any]], config_path: str = "config.json"
) -> None:
    """
    Save the model configurations to a JSON file.

    Parameters
    ----------
    configs : list of dict
        List of model configurations
    config_path : str, default="config.json"
        Path to the configuration file

    Examples
    --------
    >>> configs = load_config()
    >>> configs[0]["threshold"] = 0.95
    >>> save_config(configs, "updated_config.json")
    """
    with open(config_path, "w") as f:
        f.write(json.dumps(configs, indent=2))
    logger.info(f"Configuration saved to: {config_path}")
