import unittest
from unittest.mock import MagicMock, patch

import torch

from decaf.models import AmpliconClassifier


class TestAmpliconClassifier(unittest.TestCase):
    def setUp(self):
        # Mock model configuration
        self.model_config = {
            "model_path": "mock_model_path",
            "categories": ["target", "contaminant"],
            "preprocessing": {"max_length": 128},
        }

        # Create patches for AutoModelForSequenceClassification and AutoTokenizer
        self.mock_tokenizer = MagicMock()
        self.mock_model = MagicMock()
        self.mock_model.classifier = MagicMock()
        self.mock_model.classifier.out_features = 2

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_classifier_initialization(
        self, mock_model_from_pretrained, mock_tokenizer_from_pretrained
    ):
        """Test classifier initialization"""
        # Set up mocks
        mock_tokenizer_from_pretrained.return_value = self.mock_tokenizer
        mock_model_from_pretrained.return_value = self.mock_model

        # Initialize classifier
        classifier = AmpliconClassifier(self.model_config)

        # Check if tokenizer and model were loaded correctly
        mock_tokenizer_from_pretrained.assert_called_once_with(
            self.model_config["model_path"], trust_remote_code=True
        )
        mock_model_from_pretrained.assert_called_once_with(
            self.model_config["model_path"]
        )

        # Check model properties
        self.assertEqual(classifier.tokenizer, self.mock_tokenizer)
        self.assertEqual(classifier.model, self.mock_model)
        self.assertEqual(classifier.categories, ["target", "contaminant"])

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_classifier_mismatch_classes(
        self, mock_model_from_pretrained, mock_tokenizer_from_pretrained
    ):
        """Test classifier raises error when classes mismatch"""
        # Set up mocks
        mock_tokenizer_from_pretrained.return_value = self.mock_tokenizer

        # Set up model with 3 output features (mismatched with config)
        mock_model = MagicMock()
        mock_model.classifier = MagicMock()
        mock_model.classifier.out_features = 3
        mock_model_from_pretrained.return_value = mock_model

        # Test that initialization raises ValueError
        with self.assertRaises(ValueError):
            AmpliconClassifier(self.model_config)

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_forward(self, mock_model_from_pretrained, mock_tokenizer_from_pretrained):
        """Test forward pass"""
        # Set up mocks
        mock_tokenizer_from_pretrained.return_value = self.mock_tokenizer
        mock_model_from_pretrained.return_value = self.mock_model

        # Initialize classifier
        classifier = AmpliconClassifier(self.model_config)

        # Create test input
        input_ids = torch.tensor([[1, 2, 3, 4]])
        attention_mask = torch.tensor([[1, 1, 1, 1]])

        # Call forward
        classifier.forward(input_ids, attention_mask)

        # Check if model's forward was called with correct arguments
        self.mock_model.assert_called_once_with(
            input_ids=input_ids, attention_mask=attention_mask
        )
