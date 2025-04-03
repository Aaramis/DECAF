# DECAF

*DECAF -  DEcontamination and Classification of Amplicon Fragment*

Bienvenue dans la documentation de DECAF, un outil bioinformatique pour la classification et la décontamination de séquences ITS (Internal Transcribed Spacer) de plantes.

## Fonctionnalités principales

- Classification des séquences ITS pour les plantes
- Détection et filtrage des contaminants
- Support des formats FASTQ/FASTA
- Interface en ligne de commande
- Architecture modulaire pour extension à d'autres marqueurs et taxons

## Démarrage rapide

```bash
# Installation
pip install decaf

# Classification de séquences
decaf -i sequences.fasta -o results.csv --barcode ITS --taxa plantes --cpus 10 --gpus 1
```


## Structure de la documentation

- [ ] Installation - Comment installer DECAF
- [ ] Guide de démarrage - Premiers pas avec DECAF
- [ ] Guides - Guides détaillés pour les tâches spécifiques
- [ ] Référence API - Documentation technique de l'API
- [ ] Contribuer - Comment contribuer au projet