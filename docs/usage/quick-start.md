# Démarrage Rapide

## Installation

```bash
pip install decaf
```

## Analyse de Base

Pour analyser vos séquences d'ADN :

```bash
decaf analyze --input data/sequences.fasta --output results/ --model its_plant
```

### Options de Base

| Option | Description | Exemple |
|--------|-------------|---------|
| `--input` | Fichier d'entrée (FASTA/FASTQ) | `--input sequences.fasta` |
| `--output` | Dossier de sortie | `--output results/` |
| `--model` | Modèle à utiliser | `--model its_plant` |
| `--batch-size` | Taille des lots | `--batch-size 32` |
| `--gpu` | Utiliser le GPU (0/1) | `--gpu 1` |

## Exemple Complet

```bash
decaf analyze \
  --input data/sequences.fasta \
  --output results/ \
  --model its_plant \
  --batch-size 64 \
  --gpu 1
```
