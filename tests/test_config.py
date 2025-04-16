import json
import unittest
from unittest.mock import mock_open, patch
from decaf.config import get_model_config, load_config, save_config


class TestConfigManagement(unittest.TestCase):
    @patch("logging.getLogger")
    def setUp(self, mock_logger):
        # Sample configuration for tests
        self.sample_config = [
            {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"},
            {"barcode": "ITS", "taxa": "Fungi", "model_path": "path/to/model2"},
            {"barcode": "COI", "taxa": "Birds", "model_path": "path/to/model3"},
        ]
        self.config_path = "config/config_models.json"
        # Now, directly use mock_logger.return_value
        self.logger_mock = mock_logger.return_value

    def tearDown(self):
        self.logger_mock.handlers.clear()

    def test_load_config_success(self):
        """Test loading a valid configuration file"""
        # Simulate a valid JSON file
        mock_file_content = json.dumps(self.sample_config)
        mock_file = mock_open(read_data=mock_file_content)
        with patch("builtins.open", mock_file):
            config = load_config(self.config_path)
            # Verify the file was opened correctly
            mock_file.assert_called_once_with(self.config_path, "r")
            # Verify the configuration is correctly loaded
            self.assertEqual(config, self.sample_config)

    def test_get_model_config_found(self):
        """Test retrieving an existing configuration"""
        result = get_model_config(self.sample_config, "COI", "Fish")
        expected = {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"}
        self.assertEqual(result, expected)

    def test_get_model_config_case_insensitive(self):
        """Test case-insensitive configuration retrieval"""
        result = get_model_config(self.sample_config, "coi", "fish")
        expected = {"barcode": "COI", "taxa": "Fish", "model_path": "path/to/model1"}
        self.assertEqual(result, expected)

    def test_get_model_config_not_found(self):
        """Test handling when configuration is not found"""
        result = get_model_config(self.sample_config, "COI", "Mammals")
        self.assertIsNone(result)

    def test_load_config_file_not_found(self):
        """Test error handling when configuration file doesn't exist"""
        # Mock the logger directly in the module where it's used
        with patch("decaf.config.config.logger") as mock_logger:
            with patch("builtins.open", side_effect=FileNotFoundError()):
                with self.assertRaises(FileNotFoundError):
                    load_config("non_existent_file.json")
                # Verify the error is properly logged
                mock_logger.error.assert_called_once_with(
                    "Configuration file not found: non_existent_file.json"
                )

    def test_load_config_invalid_json(self):
        """Test error handling when file contains invalid JSON"""
        # Simulate an invalid JSON file
        mock_file = mock_open(read_data="{invalid json")
        # Mock the logger directly in the module where it's used
        with patch("decaf.config.config.logger") as mock_logger:
            with patch("builtins.open", mock_file):
                with self.assertRaises(json.JSONDecodeError):
                    load_config(self.config_path)
                # Verify the error is properly logged
                mock_logger.error.assert_called_once_with(
                    f"Invalid JSON in configuration file: {self.config_path}"
                )

    def test_save_config_success(self):
        """Test successfully saving configuration to a file"""
        # Create a mock for the open function
        mock_file = mock_open()

        # Mock the logger directly in the module where it's used
        with patch("decaf.config.config.logger") as mock_logger:
            with patch("builtins.open", mock_file):
                # Call the function
                save_config(self.sample_config, "test_config.json")

                # Verify the file was opened correctly for writing
                mock_file.assert_called_once_with("test_config.json", "w")

                # Verify the JSON was written to the file
                handle = mock_file()
                expected_json = json.dumps(self.sample_config, indent=2)
                handle.write.assert_called_once_with(expected_json)

                # Verify the success message was logged
                mock_logger.info.assert_called_once_with(
                    "Configuration saved to: test_config.json"
                )
