![DECAF](https://raw.githubusercontent.com/Aaramis/DECAF/main/docs/source/images/decaf_logo.png)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![PyPI version](https://badge.fury.io/py/decaf.svg)](https://pypi.org/project/decaf/)
[![Documentation Status](https://readthedocs.org/projects/decaf/badge/?version=latest)](https://decaf.readthedocs.io/en/latest/?badge=latest)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# DECAF: Deep Learning Framework for Environmental Contaminant Analysis in DNA Sequences

DECAF (Deep Learning Framework for Environmental Contaminant Analysis in DNA Sequences) est un framework bioinformatique moderne conçu pour l'analyse et la décontamination de séquences d'ADN environnementales. Il utilise des modèles d'apprentissage profond pour améliorer la fiabilité des analyses génomiques environnementales.

## 📋 Description

DECAF offre une solution complète pour :
- La classification de séquences d'ADN ITS (Internal Transcribed Spacer)
- La détection et la filtration des contaminants
- L'analyse de séquences environnementales à grande échelle
- L'intégration dans des pipelines bioinformatiques existants

## 🚀 Caractéristiques Principales

- **Classification Avancée**
  - Modèles de deep learning optimisés pour l'ADN
  - Support des formats FASTQ et FASTA
  - Interface en ligne de commande intuitive

- **Gestion des Contaminants**
  - Détection précise des séquences non-cibles
  - Filtrage automatique des contaminants
  - Rapports détaillés d'analyse

- **Performance et Scalabilité**
  - Optimisé pour le traitement par lots
  - Support GPU via PyTorch
  - Architecture modulaire extensible

- **Documentation Complète**
  - Guide d'utilisation détaillé
  - Exemples de cas d'utilisation
  - API documentation

## 📦 Installation

### Prérequis

- Python 3.8 ou supérieur
- Git
- Une carte graphique NVIDIA (recommandé pour le traitement rapide)

### Installation via PyPI

```bash
pip install decaf
```

### Installation depuis le code source

1. Cloner le dépôt :
```bash
git clone https://github.com/Aaramis/DECAF.git
cd DECAF
```

2. Créer un environnement virtuel :
```bash
python -m venv decaf-env
source decaf-env/bin/activate  # Sur Linux/Mac
# decaf-env\Scripts\activate  # Sur Windows
```

3. Installer les dépendances :
```bash
pip install -r requirements.txt
```

## 🏃‍♂️ Utilisation Rapide

```bash
decaf analyze --input data/sequences.fasta --output results/ --model its_plant
```

Pour plus d'options et d'exemples, consultez la documentation complète.

## 📚 Documentation

La documentation complète est disponible sur :
[https://decaf.readthedocs.io](https://decaf.readthedocs.io)

Pour générer la documentation localement :

1. Installer les dépendances de développement :
```bash
pip install -r requirements.txt
```

2. Lancer le serveur de documentation :
```bash
mkdocs serve
```

Puis ouvrir votre navigateur à l'adresse : http://127.0.0.1:8000

## 🤝 Contribuer

Nous accueillons avec plaisir les contributions à DECAF !

1. Ouvrez une issue pour signaler des bugs ou proposer des fonctionnalités
2. Créez une pull request pour contribuer du code
3. Suivez les directives de style de code

Pour vérifier le style de votre code :
```bash
pip install black
black --check .
```

Pour formater votre code :
```bash
black .
```

## 🏗️ Structure du Projet

```
DECAF/
├── decaf/                 # Code source principal
│   ├── models/           # Implémentation des modèles
│   ├── data/             # Gestion des données
│   └── utils/            # Fonctions utilitaires
├── tests/                # Tests unitaires et d'intégration
├── docs/                 # Documentation
├── config/               # Fichiers de configuration
└── data/                 # Données d'exemple
```

## 📜 Licence

DECAF est sous licence MIT. Voir le fichier [LICENSE](LICENSE) pour plus de détails.

## 🙏 Remerciements

- [Auguste_GARDETTE](https://github.com/Aaramis) - Développeur principal
- [Contributors](https://github.com/Aaramis/DECAF/graphs/contributors) - Tous les contributeurs

## 📞 Support

Pour toute question ou problème, veuillez ouvrir une issue sur GitHub ou contacter l'équipe de développement.