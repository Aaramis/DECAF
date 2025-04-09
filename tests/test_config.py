import json
import unittest
from unittest.mock import mock_open, patch

from decaf.config import get_model_config, load_config


class TestConfigManagement(unittest.TestCase):

    @patch("logging.getLogger")
    def setUp(self, mock_logger):
        # Exemple de configuration pour les tests
        self.sample_config = [
            {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"},
            {"barcode": "ITS", "taxa": "Fungi", "model_path": "path/to/model2"},
            {"barcode": "COI", "taxa": "Birds", "model_path": "path/to/model3"},
        ]
        self.config_path = "config/config_models.json"

        # Maintenant, on utilise directement mock_logger.return_value
        self.logger_mock = mock_logger.return_value

    def tearDown(self):
        self.logger_mock.handlers.clear()

    def test_load_config_success(self):
        """Test le chargement d'un fichier de configuration valide"""
        # Simuler un fichier JSON valide
        mock_file_content = json.dumps(self.sample_config)
        mock_file = mock_open(read_data=mock_file_content)

        with patch("builtins.open", mock_file):
            config = load_config(self.config_path)

        # Vérifier que le fichier a été ouvert correctement
        mock_file.assert_called_once_with(self.config_path, "r")

        # Vérifier que la configuration est correctement chargée
        self.assertEqual(config, self.sample_config)

    def test_get_model_config_found(self):
        """Test la récupération d'une configuration existante"""
        result = get_model_config(self.sample_config, "COI", "Fish")
        expected = {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"}
        self.assertEqual(result, expected)

    def test_get_model_config_case_insensitive(self):
        result = get_model_config(self.sample_config, "coi", "fish")
        expected = {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"}
        self.assertEqual(result, expected)

    def test_get_model_config_not_found(self):
        result = get_model_config(self.sample_config, "COI", "Mammals")
        self.assertIsNone(result)
