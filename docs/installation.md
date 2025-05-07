# Installation de DECAF

## Installation via PyPI

La manière la plus simple d'installer DECAF est via PyPI :

```bash
pip install decaf
```

## Installation depuis le Code Source

### Prérequis

Avant d'installer DECAF, assurez-vous d'avoir les dépendances suivantes :

- Python 3.8 ou supérieur
- Git
- Une carte graphique NVIDIA (recommandé pour le traitement rapide)
- pip (gestionnaire de paquets Python)

### Étapes d'Installation

1. **Cloner le Dépôt**
```bash
git clone https://github.com/Aaramis/DECAF.git
cd DECAF
```

2. **Créer un Environnement Virtuel**

Il est fortement recommandé d'utiliser un environnement virtuel pour gérer les dépendances.

Sur Linux/Mac :
```bash
python -m venv decaf-env
source decaf-env/bin/activate
```

Sur Windows :
```bash
python -m venv decaf-env
decaf-env\Scripts\activate
```

3. **Installer les Dépendances**

Installez les dépendances nécessaires :
```bash
pip install -r requirements.txt
```

Pour les développeurs, installez également les dépendances de développement :
```bash
pip install -r requirements-dev.txt
```

4. **Installer DECAF en Mode Développement**
```bash
pip install -e .
```

## Vérification de l'Installation

Pour vérifier que DECAF est correctement installé, exécutez :
```bash
decaf --version
```

## Configuration de l'Environnement

### Variables d'Environnement

DECAF utilise les variables d'environnement suivantes :

- `DECAF_CONFIG_PATH` : Chemin vers le fichier de configuration
- `DECAF_LOG_LEVEL` : Niveau de log (DEBUG, INFO, WARNING, ERROR)
- `DECAF_GPU` : Utiliser le GPU (0 pour désactiver)

### Configuration du Log

La configuration du log peut être personnalisée dans le fichier `config/logging.yaml`.

## Dépannage

### Problèmes Communs

1. **Erreur de Version Python**
   - Assurez-vous d'utiliser Python 3.8 ou supérieur
   - Vérifiez la version avec : `python --version`

2. **Problèmes de Permissions**
   - Utilisez un environnement virtuel
   - Exécutez les commandes avec les droits appropriés

3. **Dépendances Manquantes**
   - Vérifiez que toutes les dépendances sont installées
   - Réinstallez les dépendances si nécessaire

## Mise à Jour

Pour mettre à jour DECAF :

```bash
pip install --upgrade decaf
```

Ou depuis le code source :
```bash
git pull
pip install --upgrade -e .
```

## Suppression

Pour désinstaller DECAF :
```bash
pip uninstall decaf
```

## Support

Pour toute question ou problème, veuillez :

1. Ouvrir une issue sur GitHub
2. Consulter la FAQ
3. Contacter l'équipe de développement
