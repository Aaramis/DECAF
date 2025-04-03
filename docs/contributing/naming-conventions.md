# Conventions de nommage pour DECAF

Ce document définit les conventions de nommage à suivre pour le projet DECAF.

## Principes généraux

- Les noms doivent être descriptifs et explicites.
- Préférer des noms longs et clairs à des noms courts et ambigus.
- Éviter les abréviations sauf si elles sont standard dans le domaine (ex: DNA, RNA, ITS).
- Utiliser l'anglais pour tous les noms.

## Modules et packages

- Noms en minuscules.
- Mots séparés par des underscores si nécessaire.
- Noms courts mais descriptifs.
- Éviter les caractères spéciaux et les chiffres.

**Exemples:**
- `preprocessing`
- `model_training`
- `data_loaders`
- `visualization`

## Classes

- Utiliser le format CamelCase (chaque mot commence par une majuscule).
- Les noms doivent être des substantifs ou des groupes nominaux.
- Suffixe spécifique pour certains types de classes:
  - `*Error` ou `*Exception` pour les exceptions
  - `*Manager` pour les classes de gestion
  - `*Factory` pour les classes factory

**Exemples:**
- `SequenceProcessor`
- `ModelTrainer`
- `DataLoader`
- `ConfigManager`
- `SequenceQualityError`

## Fonctions et méthodes

- Utiliser le format snake_case (minuscules avec underscores).
- Utiliser des verbes ou phrases verbales.
- Nom descriptif du comportement.
- Méthodes privées préfixées par un underscore (`_`).

**Exemples:**
- `clean_sequence()`
- `train_model()`
- `calculate_similarity_score()`
- `_validate_input()`

## Variables et attributs

- Utiliser le format snake_case.
- Nom précis et descriptif du contenu.
- Attributs privés préfixés par un underscore (`_`).
- Éviter les noms d'une seule lettre sauf pour les compteurs.

**Exemples:**
- `sequence_data`
- `model_parameters`
- `gc_content`
- `_internal_state`

## Constantes

- Tout en MAJUSCULES avec underscores.
- Définies au niveau du module.

**Exemples:**
- `MAX_SEQUENCE_LENGTH`
- `DEFAULT_QUALITY_THRESHOLD`
- `DNA_BASES`

## Arguments de fonctions

- Utiliser le format snake_case.
- Être cohérent avec les noms des variables.
- Pour les booléens, utiliser des préfixes tels que `is_`, `has_`, `should_`.

**Exemples:**
- `def process_sequence(sequence_data, is_quality_filtered=True):`

## Itérateurs et boucles

- Pour les itérateurs courts ou les indices, des noms courts sont acceptables.
- Si l'objet a une signification spécifique, utiliser un nom plus descriptif.

**Exemples:**
```python
# Court et simple
for i in range(10):
    print(i)

# Plus descriptif
for sequence in sequences:
    process_sequence(sequence)

# Avec index et valeur
for idx, base in enumerate(sequence):
    if is_valid_base(base):
        valid_bases[idx] = base
```

## Fichiers

- Noms en minuscules.
- Mots séparés par des underscores.
- Extensions appropriées (.py, .md, .txt, etc.).
- Les noms de fichiers de test doivent commencer par `test_`.

**Exemples:**
- `sequence_processor.py`
- `model_training.py`
- `test_sequence_processor.py`

## Tables de base de données et champs

- Tables au pluriel, en snake_case.
- Champs en snake_case.
- Clé primaire nommée `id`.
- Clés étrangères nommées `*_id`.

**Exemples:**
- Table: `sequences`, Champs: `id`, `sequence_data`, `quality_score`
- Table: `taxonomies`, Champs: `id`, `sequence_id`, `taxon_name`

## Types génériques et TypeHints

- Utiliser CamelCase pour les types génériques personnalisés.
- Suivre les conventions de la bibliothèque standard pour les annotations de type.

**Exemples:**
```python
from typing import List, Dict, Optional, Union

SequenceData = Dict[str, Union[str, float]]
TaxonomyTree = Dict[str, List[str]]

def process_sequences(sequences: List[str]) -> List[SequenceData]:
    pass
```