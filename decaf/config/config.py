"""
Configuration management for DECAF.
"""

import json
import logging
from typing import Dict, List, Any, Optional

logger = logging.getLogger(__name__)


def load_config(config_path: str = "config/config_models.json") -> List[Dict[str, Any]]:
    """
    Load the model configuration from a JSON file.

    Args:
        config_path: Path to the configuration file

    Returns:
        List of model configurations

    Raises:
        FileNotFoundError: If the configuration file is not found
        json.JSONDecodeError: If the configuration file is not valid JSON
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

    Args:
        configs: List of model configurations
        barcode: Molecular barcode type
        taxa: Taxonomic group

    Returns:
        Model configuration if found, None otherwise
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

    Args:
        configs: List of model configurations
        config_path: Path to the configuration file
    """
    with open(config_path, "w") as f:
        json.dump(configs, f, indent=2)
    logger.info(f"Configuration saved to: {config_path}")
