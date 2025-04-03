# Standards de codage pour DECAF

Ce document définit les normes de codage à suivre pour le projet DECAF.

## Style de code général

- Suivre la norme [PEP 8](https://www.python.org/dev/peps/pep-0008/) pour le style de code Python.
- Utiliser [Black](https://black.readthedocs.io/) avec les paramètres par défaut pour le formatage automatique du code.
- Utiliser [isort](https://pycqa.github.io/isort/) pour organiser les imports.
- Utiliser [flake8](https://flake8.pycqa.org/) pour la vérification statique du code.

## Structure des fichiers

- Chaque module doit avoir un fichier `__init__.py`.
- Les fichiers de test doivent être nommés `test_*.py` et placés dans le répertoire `tests/`.
- Les scripts exécutables doivent être placés dans le répertoire `decaf/`.

## Longueur des lignes et indentation

- Limite de 88 caractères par ligne (norme Black).
- Utiliser 4 espaces pour l'indentation (pas de tabulations).
- Les commentaires de ligne doivent être limités à 72 caractères.

## Nommage

- Classes: Utiliser le format CamelCase (`class ModelTrainer:`).
- Fonctions et variables: Utiliser le format snake_case (`def process_sequence():`, `sequence_data = ...`).
- Constantes: Utiliser des majuscules avec underscore (`MAX_SEQUENCE_LENGTH = 1000`).
- Modules: Utiliser des noms courts, tout en minuscules (`utils.py`, `preprocessing.py`).
- Éviter les noms de variables d'une seule lettre, sauf pour les compteurs ou les itérateurs.

## Docstrings

- Toutes les fonctions, classes et méthodes publiques doivent avoir des docstrings.
- Utiliser le format [NumPy/SciPy](https://numpydoc.readthedocs.io/en/latest/format.html) pour les docstrings.
- Exemple de docstring pour une fonction:

```python
from typing import Tuple

def align_sequences(seq1: str, seq2: str, method: str = 'global') -> Tuple[str, str, float]:
    """
    Aligne deux séquences biologiques.
    
    Parameters
    ----------
    seq1 : str
        Première séquence à aligner
    seq2 : str
        Deuxième séquence à aligner
    method : {'global', 'local', 'semi-global'}, default='global'
        Méthode d'alignement à utiliser
        
    Returns
    -------
    tuple
        Un tuple contenant les séquences alignées (str, str) et le score (float)
        
    Examples
    --------
    >>> align_sequences("ACGT", "ACT")
    ('ACGT', 'AC-T', 2.5)
    """

```

## Commentaires

- Les commentaires doivent expliquer le "pourquoi", pas le "comment".
- Pour les sections complexes, ajouter des commentaires expliquant la logique.
- Les commentaires TODO doivent inclure un identifiant de la personne responsable:
  ```python
  # TODO(@username): Implémenter la fonction de normalisation
  ```

## Gestion des erreurs

- Utiliser des exceptions spécifiques, pas des exceptions génériques.
- Définir des exceptions personnalisées dans un module `exceptions.py`.
- Les messages d'erreur doivent être informatifs et suggérer des solutions.

## Tests unitaires

- Chaque fonction doit avoir au moins un test associé.
- Utiliser pytest comme framework de test.
- Les tests doivent être indépendants les uns des autres.
- Les tests doivent vérifier les cas normaux et les cas limites.

## Validation du code

Avant chaque commit, assurez-vous que:
1. Tous les tests passent (`pytest`)
2. Le code est correctement formaté (`black .` et `isort .`)
3. Aucun avertissement n'est signalé par flake8 (`flake8`)