# Basic Configuration

## Configuration File Structure

DECAF configuration files use the JSON format. Here is the basic structure :

```json
[
    {
        "model_name": "ITS_Plant",
        "barcode": "ITS",
        "taxa": "plants",
        "categories": {
            "1": "plants",
            "0": "contaminations"
        },
        "model_path": "models/ITS_Plant",
        "preprocessing": {
            "max_length": 150,
            "truncation": true,
            "padding": "max_length"
        },
        "training": {
            "batch_size": 32,
            "learning_rate": 2e-5,
            "epochs": 5
        }
    }
]
```

## Configuration Fields

### Model Configuration

| Field | Description | Example Value |
|-------|-------------|---------------|
| `model_name` | Name of the model | `ITS_Plant` |
| `barcode` | Barcode marker used | `ITS` |
| `taxa` | Target taxa for classification | `plants` |
| `categories` | Mapping of prediction categories | `{"1": "plants", "0": "contaminations"}` |
| `model_path` | Path to the model files | `models/ITS_Plant` |

### Preprocessing Parameters

| Field | Description | Default Value |
|-------|-------------|---------------|
| `max_length` | Maximum sequence length | `150` |
| `truncation` | Whether to truncate sequences | `true` |
| `padding` | Padding strategy | `max_length` |

### Training Parameters

| Field | Description | Default Value |
|-------|-------------|---------------|
| `batch_size` | Batch size | `32` |
| `learning_rate` | Learning rate | `2e-5` |
| `epochs` | Number of epochs | `5` |

## System Configuration

The system configuration is handled separately and is not part of the model configuration JSON. It should be configured through environment variables or command line arguments.

## Complete Configuration Example

```json
[
    {
        "model_name": "ITS_Plant",
        "barcode": "ITS",
        "taxa": "plants",
        "categories": {
            "1": "plants",
            "0": "contaminations"
        },
        "model_path": "models/ITS_Plant",
        "preprocessing": {
            "max_length": 200,
            "truncation": true,
            "padding": "max_length"
        },
        "training": {
            "batch_size": 64,
            "learning_rate": 2e-5,
            "epochs": 10
        }
    }
]
```

## Best Practices

1. **Organisation des Fichiers**
   - Separate configurations by analysis type
   - Maintain a reference configuration
   - Document major changes

2. **Version Management**
   - Version the configuration files
   - Maintain a history of changes
   - Test each new configuration

3. **Performance Optimization**
   - Adjust `batch_size` based on GPU memory
   - Adjust `num_workers` based on the number of cores
   - Adjust `max_length` based on input data

## Troubleshooting

### Common Issues

1. **Insufficient Memory**
   - Reduce `batch_size`
   - Disable GPU usage through environment variables
   - Increase `num_workers`

2. **Slow Performance**
   - Increase `num_workers`
   - Optimize `max_length`
   - Check GPU configuration

3. **Format Issues**
   - Validate sequence input
   - Check maximum length
   - Verify category mapping in `categories` field
