# Modèles DECAF

## Vue d'Ensemble

DECAF utilise des modèles d'apprentissage profond optimisés pour l'analyse d'ADN. Les principaux modèles disponibles sont :

- `its_plant`: Modèle spécialisé pour la classification d'ITS chez les plantes
- D'autres modèles seront ajoutés en fonction des besoins

## Architecture des Modèles

```python
class DNAClassificationModel(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.encoder = TransformerEncoder(config)
        self.classifier = nn.Linear(config.hidden_size, config.num_labels)
        self.dropout = nn.Dropout(config.dropout_rate)

    def forward(self, input_ids, attention_mask):
        # Encodage des séquences
        embeddings = self.encoder(input_ids, attention_mask)
        # Classification
        logits = self.classifier(self.dropout(embeddings))
        return logits
```

## Configuration des Modèles

Les modèles peuvent être configurés via un fichier YAML :

```yaml
model:
  name: its_plant
  hidden_size: 768
  num_labels: 2
  dropout_rate: 0.1
  max_length: 150
  learning_rate: 5e-5
```

## Utilisation des Modèles

```python
from decaf.models import load_model

# Charger un modèle
model = load_model("its_plant")

# Prédiction sur une séquence
sequence = "ATCG..."
result = model.predict(sequence)
```