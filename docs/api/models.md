# DECAF Models

## Overview

DECAF uses deep learning models optimized for DNA analysis. The main models available are :

- `its_plant`: Model specialized for ITS classification in plants
- Other models will be added as needed

## Models Architecture

```python
class DNAClassificationModel(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.encoder = TransformerEncoder(config)
        self.classifier = nn.Linear(config.hidden_size, config.num_labels)
        self.dropout = nn.Dropout(config.dropout_rate)

    def forward(self, input_ids, attention_mask):
        # Sequence encoding
        embeddings = self.encoder(input_ids, attention_mask)
        # Classification
        logits = self.classifier(self.dropout(embeddings))
        return logits
```

## Models Configuration

Models can be configured via a YAML file :

```yaml
model:
  name: its_plant
  hidden_size: 768
  num_labels: 2
  dropout_rate: 0.1
  max_length: 150
  learning_rate: 5e-5
```

## Models Usage

```python
from decaf.models import load_model

# Load a model
model = load_model("its_plant")

# Prediction on a sequence
sequence = "ATCG..."
result = model.predict(sequence)
```