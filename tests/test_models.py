import unittest
from unittest.mock import MagicMock, patch

import torch

from decaf.models import AmpliconClassifier, load_model


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

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_predict_step(
        self, mock_model_from_pretrained, mock_tokenizer_from_pretrained
    ):
        """Test predict_step method"""
        # Set up mocks
        mock_tokenizer_from_pretrained.return_value = self.mock_tokenizer
        mock_model_from_pretrained.return_value = self.mock_model

        # Mock output from the model
        mock_outputs = MagicMock()
        # Create a 2x2 tensor with values [0.3, 0.7] and [0.6, 0.4]
        mock_logits = torch.tensor([[0.3, 0.7], [0.6, 0.4]])
        mock_outputs.logits = mock_logits
        self.mock_model.return_value = mock_outputs

        # Initialize classifier
        classifier = AmpliconClassifier(self.model_config)

        # Create test batch
        batch = {
            "input_ids": torch.tensor([[1, 2, 3], [4, 5, 6]]),
            "attention_mask": torch.tensor([[1, 1, 1], [1, 1, 1]]),
            "sequence_id": ["seq1", "seq2"],
            "sequence_str": ["ATCG", "GCTA"],
        }

        # Call predict_step
        result = classifier.predict_step(batch, batch_idx=0)

        # Expected probabilities after softmax on mock_logits
        expected_probs = torch.softmax(mock_logits, dim=1)
        # Expected predictions and confidences
        expected_preds = torch.tensor([1, 0])  # Class indices with max probability
        expected_confidences = torch.tensor(
            [0.6682, 0.6457]
        )  # Max probabilities (softmax values)

        # Check results
        self.assertEqual(result["sequence_ids"], batch["sequence_id"])
        self.assertEqual(result["sequence_strs"], batch["sequence_str"])
        expected_confidences, _ = torch.max(expected_probs, dim=1)
        self.assertTrue(torch.allclose(result["predictions"], expected_preds))
        self.assertTrue(
            torch.allclose(result["confidences"], expected_confidences, rtol=1e-4)
        )
        self.assertTrue(
            torch.allclose(result["probabilities"], expected_probs, rtol=1e-4)
        )

    @patch("decaf.models.models.AmpliconClassifier")
    def test_load_model(self, mock_amplicon_classifier):
        """Test load_model function"""
        # Set up mock
        mock_classifier_instance = MagicMock()
        mock_amplicon_classifier.return_value = mock_classifier_instance

        # Call load_model
        result = load_model(self.model_config)

        # Check if AmpliconClassifier was initialized correctly
        mock_amplicon_classifier.assert_called_once_with(self.model_config)
        self.assertEqual(result, mock_classifier_instance)

    @patch("decaf.models.models.AutoTokenizer.from_pretrained")
    @patch("decaf.models.models.AutoModelForSequenceClassification.from_pretrained")
    def test_classifier_without_preprocessing(
        self, mock_model_from_pretrained, mock_tokenizer_from_pretrained
    ):
        """Test initialization without preprocessing parameter"""
        # Set up mocks
        mock_tokenizer_from_pretrained.return_value = self.mock_tokenizer
        mock_model_from_pretrained.return_value = self.mock_model

        # Create config without preprocessing
        config_without_preprocessing = {
            "model_path": "mock_model_path",
            "categories": ["target", "contaminant"],
        }

        # Initialize classifier
        classifier = AmpliconClassifier(config_without_preprocessing)

        # Check that preprocessing is an empty dict
        self.assertEqual(classifier.preprocessing, {})
