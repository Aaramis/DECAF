# Configuration de Base

## Structure des Fichiers de Configuration

Les fichiers de configuration de DECAF utilisent le format YAML. Voici la structure de base :

```yaml
# Configuration du modèle
model:
  name: its_plant
  hidden_size: 768
  num_labels: 2
  dropout_rate: 0.1

# Paramètres d'entraînement
training:
  batch_size: 32
  learning_rate: 5e-5
  epochs: 10
  patience: 3

# Paramètres de données
data:
  max_length: 150
  num_workers: 4
  input_format: fasta

# Paramètres système
system:
  use_gpu: true
  num_gpus: 1
  num_nodes: 1
  seed: 42
```

## Paramètres du Modèle

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| `name` | Nom du modèle | `its_plant` |
| `hidden_size` | Taille des vecteurs cachés | `768` |
| `num_labels` | Nombre de classes | `2` |
| `dropout_rate` | Taux de dropout | `0.1` |

## Paramètres d'Entraînement

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| `batch_size` | Taille des lots | `32` |
| `learning_rate` | Taux d'apprentissage | `5e-5` |
| `epochs` | Nombre d'époques | `10` |
| `patience` | Patience pour early stopping | `3` |

## Paramètres de Données

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| `max_length` | Longueur maximale des séquences | `150` |
| `num_workers` | Nombre de workers pour le chargement | `4` |
| `input_format` | Format des données d'entrée | `fasta` |

## Paramètres Système

| Paramètre | Description | Valeur par défaut |
|-----------|-------------|-------------------|
| `use_gpu` | Utiliser le GPU | `true` |
| `num_gpus` | Nombre de GPUs | `1` |
| `num_nodes` | Nombre de nœuds | `1` |
| `seed` | Seed pour la reproductibilité | `42` |

## Exemple de Configuration Complète

```yaml
model:
  name: its_plant
  hidden_size: 768
  num_labels: 2
  dropout_rate: 0.1

training:
  batch_size: 64
  learning_rate: 5e-5
  epochs: 20
  patience: 5





data:
  max_length: 200
  num_workers: 8
  input_format: fasta

system:
  use_gpu: true
  num_gpus: 2
  num_nodes: 1
  seed: 42
```

## Bonnes Pratiques

1. **Organisation des Fichiers**
   - Séparer les configurations par type d'analyse
   - Maintenir une configuration de référence
   - Documenter les changements majeurs

2. **Gestion des Versions**
   - Versionner les fichiers de configuration
   - Maintenir un historique des modifications
   - Tester chaque nouvelle configuration

3. **Optimisation des Performances**
   - Adapter `batch_size` selon la mémoire GPU
   - Optimiser `num_workers` selon le nombre de cœurs
   - Ajuster `max_length` selon les données d'entrée

## Dépannage

### Problèmes Communs

1. **Mémoire Insuffisante**
   - Réduire `batch_size`
   - Désactiver `use_gpu`
   - Augmenter `num_workers`

2. **Performance Lente**
   - Augmenter `num_workers`
   - Optimiser `max_length`
   - Vérifier la configuration GPU

3. **Problèmes de Format**
   - Vérifier `input_format`
   - Valider les séquences d'entrée
   - Vérifier la longueur maximale
