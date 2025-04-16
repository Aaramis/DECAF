"""
Model management for DECAF using PyTorch Lightning.

This module provides classes and functions for loading and using
pre-trained transformer models for amplicon classification.

Functions
---------
load_model
    Load a model based on the configuration.

Classes
-------
AmpliconClassifier
    Amplicon classifier using pre-trained transformer models.
"""

import logging
from typing import Dict, Any

import torch
import pytorch_lightning as pl
from transformers import AutoModelForSequenceClassification, AutoTokenizer


logger = logging.getLogger(__name__)


class AmpliconClassifier(pl.LightningModule):
    """
    Complete amplicon classifier for inference using a pre-trained Lightning model.

    This class wraps a HuggingFace transformer model for sequence classification
    and provides methods for prediction.

    Parameters
    ----------
    model_config : dict
        Configuration dictionary containing model parameters

    Attributes
    ----------
    model_path : str
        Path to the pre-trained model
    categories : list
        List of output categories/classes
    preprocessing : dict
        Preprocessing configuration options
    tokenizer : AutoTokenizer
        HuggingFace tokenizer for the model
    model : AutoModelForSequenceClassification
        Pre-trained transformer model

    Raises
    ------
    ValueError
        If the model output dimensions don't match the expected number of classes
    """

    def __init__(self, model_config: Dict[str, Any]):
        super().__init__()
        self.model_config = model_config
        self.model_path = model_config["model_path"]
        self.categories = model_config["categories"]
        self.preprocessing = model_config.get("preprocessing", {})

        num_classes = len(self.categories)

        # Load tokenizer and model
        self.tokenizer = AutoTokenizer.from_pretrained(
            self.model_path, trust_remote_code=True
        )
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_path)

        # Check classifier output size
        actual_out_features = self.model.classifier.out_features
        if actual_out_features != num_classes:
            raise ValueError(
                f"Model classifier output features ({actual_out_features}) "
                f"do not match expected number of classes ({num_classes})"
            )

        self.model.eval()

    def forward(self, input_ids, attention_mask):
        """
        Forward pass through the model.

        Parameters
        ----------
        input_ids : torch.Tensor
            Token IDs of input sequences
        attention_mask : torch.Tensor
            Attention mask for input sequences

        Returns
        -------
        transformers.modeling_outputs.SequenceClassifierOutput
            Model outputs including logits

        Examples
        --------
        >>> classifier = AmpliconClassifier(model_config)
        >>> outputs = classifier(input_ids, attention_mask)
        >>> outputs.logits.shape
        torch.Size([batch_size, num_classes])
        """
        return self.model(input_ids=input_ids, attention_mask=attention_mask)

    def predict_step(self, batch, batch_idx):
        """
        PyTorch Lightning predict step.

        Parameters
        ----------
        batch : dict
            Batch of data containing input_ids, attention_mask, etc.
        batch_idx : int
            Index of the batch

        Returns
        -------
        dict
            Dictionary containing predictions and metadata:
            - sequence_ids: Original sequence identifiers
            - sequence_strs: Original sequence strings
            - seq_qualities: Original sequence qualities
            - predictions: Class predictions
            - confidences: Confidence scores for predictions
            - probabilities: Full probability distribution over classes

        Examples
        --------
        >>> predictions = classifier.predict_step(batch, 0)
        >>> predictions["predictions"].shape
        torch.Size([batch_size])
        """
        input_ids = batch["input_ids"]
        attention_mask = batch["attention_mask"]

        with torch.no_grad():
            outputs = self(input_ids, attention_mask)
            logits = outputs.logits
            probs = torch.softmax(logits, dim=1)

            # Get predictions and confidence scores
            confidences, preds = torch.max(probs, dim=1)

            return {
                "sequence_ids": batch["sequence_id"],
                "sequence_strs": batch["sequence_str"],
                "seq_qualities": batch["seq_quality"],
                "predictions": preds,
                "confidences": confidences,
                "probabilities": probs,
            }


def load_model(model_config: Dict[str, Any]) -> AmpliconClassifier:
    """
    Load a model based on the configuration.

    Parameters
    ----------
    model_config : dict
        Model configuration from the config file

    Returns
    -------
    AmpliconClassifier
        Initialized AmpliconClassifier

    Examples
    --------
    >>> config = {"model_path": "path/to/model", "categories": ["classA", "classB"]}
    >>> model = load_model(config)
    """
    return AmpliconClassifier(model_config)
