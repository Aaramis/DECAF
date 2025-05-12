# Naming Conventions for DECAF

This document defines the naming conventions to be followed for the DECAF project.

## General Principles

- Names must be descriptive and explicit.
- Prefer long and clear names to short and ambiguous names.
- Avoid abbreviations unless they are standard in the domain (e.g., DNA, RNA, ITS).
- Use English for all names.

## Modules and packages

- Use lowercase.
- Words separated by underscores if necessary.
- Short but descriptive names.
- Avoid special characters and numbers.

**Examples:**
- `preprocessing`
- `model_training`
- `data_loaders`
- `visualization`

## Classes

- Use the CamelCase format (each word starts with a capital letter).
- Names must be nouns or nominal groups.
- Specific suffixes for certain types of classes:
  - `*Error` or `*Exception` for exceptions
  - `*Manager` for management classes
  - `*Factory` for factory classes

**Examples:**
- `SequenceProcessor`
- `ModelTrainer`
- `DataLoader`
- `ConfigManager`
- `SequenceQualityError`

## Functions and methods

- Use the snake_case format (lowercase with underscores).
- Use verbs or verb phrases.
- Descriptive name of the behavior.
- Private methods prefixed with an underscore (`_`).

**Examples:**
- `clean_sequence()`
- `train_model()`
- `calculate_similarity_score()`
- `_validate_input()`

## Variables and attributes

- Use the snake_case format.
- Precise and descriptive name of the content.
- Private attributes prefixed with an underscore (`_`).
- Avoid single-letter names except for counters.

**Examples:**
- `sequence_data`
- `model_parameters`
- `gc_content`
- `_internal_state`

## Constants

- All uppercase with underscores.
- Defined at the module level.

**Examples:**
- `MAX_SEQUENCE_LENGTH`
- `DEFAULT_QUALITY_THRESHOLD`
- `DNA_BASES`

## Function arguments

- Use the snake_case format.
- Be consistent with variable names.
- For booleans, use prefixes such as `is_`, `has_`, `should_`.

**Examples:**
- `def process_sequence(sequence_data, is_quality_filtered=True):`

## Iterators and loops

- For short iterators or indices, short names are acceptable.
- If the object has a specific meaning, use a more descriptive name.

**Examples:**
```python
# Short and simple
for i in range(10):
    print(i)

# More descriptive
for sequence in sequences:
    process_sequence(sequence)

# With index and value
for idx, base in enumerate(sequence):
    if is_valid_base(base):
        valid_bases[idx] = base
```

## Files

- Use lowercase.
- Words separated by underscores.
- Appropriate extensions (.py, .md, .txt, etc.).
- Test file names must start with `test_`.

**Examples:**
- `sequence_processor.py`
- `model_training.py`
- `test_sequence_processor.py`

## Database tables and fields

- Tables in plural, snake_case.
- Fields in snake_case.
- Primary key named `id`.
- Foreign key named `*_id`.

**Examples:**
- Table: `sequences`, Fields: `id`, `sequence_data`, `quality_score`
- Table: `taxonomies`, Fields: `id`, `sequence_id`, `taxon_name`

## Generic types and TypeHints

- Use CamelCase for custom generic types.
- Follow the conventions of the standard library for type annotations.

**Examples:**
```python
from typing import List, Dict, Optional, Union

SequenceData = Dict[str, Union[str, float]]
TaxonomyTree = Dict[str, List[str]]

def process_sequences(sequences: List[str]) -> List[SequenceData]:
    pass
```