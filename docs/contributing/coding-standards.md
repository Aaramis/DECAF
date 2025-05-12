# DECAF Coding Standards

This document defines the coding standards to be followed for the DECAF project.

## General Code Style

- Follow the [PEP 8](https://www.python.org/dev/peps/pep-0008/) Python code style.
- Use [Black](https://black.readthedocs.io/) with default parameters for automatic code formatting.
- Use [isort](https://pycqa.github.io/isort/) to organize imports.
- Use [flake8](https://flake8.pycqa.org/) for static code analysis.

## File Structure

- Each module must have an `__init__.py` file.
- Test files must be named `test_*.py` and placed in the `tests/` directory.
- Executable scripts must be placed in the `decaf/` directory.

## Line Length and Indentation

- Limit of 88 characters per line (Black norm).
- Use 4 spaces for indentation (no tabs).
- Line comments must be limited to 72 characters.

## Naming

- Classes: Use the CamelCase format (`class ModelTrainer:`).
- Functions and variables: Use the snake_case format (`def process_sequence():`, `sequence_data = ...`).
- Constants: Use uppercase with underscores (`MAX_SEQUENCE_LENGTH = 1000`).
- Modules: Use short names in lowercase (`utils.py`, `preprocessing.py`).
- Avoid single-letter variable names, except for counters or iterators.

## Docstrings

- All public functions, classes, and methods must have docstrings.
- Use the [NumPy/SciPy](https://numpydoc.readthedocs.io/en/latest/format.html) format for docstrings.
- Example docstring for a function:

```python
from typing import Tuple

def align_sequences(seq1: str, seq2: str, method: str = 'global') -> Tuple[str, str, float]:
    """
    Align two biological sequences.
    
    Parameters
    ----------
    seq1 : str
        First sequence to align
    seq2 : str
        Second sequence to align
    method : {'global', 'local', 'semi-global'}, default='global'
        Alignment method to use
        
    Returns
    -------
    tuple
        A tuple containing the aligned sequences (str, str) and the score (float)
        
    Examples
    --------
    >>> align_sequences("ACGT", "ACT")
    ('ACGT', 'AC-T', 2.5)
    """

```

## Comments

- Comments should explain the "why", not the "how".
- For complex sections, add comments explaining the logic.
- TODO comments should include the responsible person's identifier:
  ```python
  # TODO(@username): Implement the normalization function
  ```

## Error Handling

- Use specific exceptions, not generic exceptions.
- Define custom exceptions in an `exceptions.py` module.
- Error messages should be informative and suggest solutions.

## Unit Tests

- Each function must have at least one associated test.
- Use pytest as the test framework.
- Tests must be independent of each other.
- Tests must verify normal and edge cases.

## Code Validation

Before each commit, make sure that:
1. All tests pass (`pytest`)
2. The code is correctly formatted (`black .` and `isort .`)
3. No warnings are reported by flake8 (`flake8`)
4. Type checking passes (`mypy decaf/`)