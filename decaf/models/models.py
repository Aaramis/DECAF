"""
Model management for DECAF using PyTorch Lightning.
"""

import torch
import logging
from typing import Dict, Any
from transformers import AutoModelForSequenceClassification, AutoTokenizer
import pytorch_lightning as pl

logger = logging.getLogger(__name__)


class AmpliconClassifier(pl.LightningModule):
    """
    Complete amplicon classifier for inference using a pre-trained Lightning model.
    """

    def __init__(self, model_config: Dict[str, Any]):
        """
        Initialize the classifier from config dict.

        Args:
            model_config: Configuration dictionary
        """
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
        return self.model(input_ids=input_ids, attention_mask=attention_mask)

    def predict_step(self, batch, batch_idx):
        """PyTorch Lightning predict step."""
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
            "predictions": preds,
            "confidences": confidences,
            "probabilities": probs,
        }


def load_model(model_config: Dict[str, Any]) -> AmpliconClassifier:
    """
    Load a model based on the configuration.

    Args:
        model_config: Model configuration from the config file

    Returns:
        Initialized AmpliconClassifier
    """
    return AmpliconClassifier(model_config)
