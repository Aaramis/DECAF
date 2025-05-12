# Docstrings

This file contains templates for different types of docstrings to be used in DECAF.
These examples serve as a reference to maintain consistent documentation.


## Template for a module
```python
"""
Module for ITS sequence processing.

This module provides functions for cleaning, normalizing and preparing
ITS sequences for analysis and classification.

Functions
---------
clean_sequence
    Clean a sequence by removing non-nucleotide characters.
normalize_sequence
    Normalize a sequence for model training.
extract_features
    Extract features from a sequence.
"""
```


## Template for a class

```python
class SequenceProcessor:
    """
    Class for biological sequence processing.
    
    This class implements various methods for cleaning, normalizing and extracting
    information from biological sequences, particularly ITS sequences from plants.
    
    Parameters
    ----------
    min_length : int, default=100
        Minimum length of sequences to process
    max_length : int, default=1000
        Maximum length of sequences to process
    quality_threshold : float, default=0.9
        Minimum quality threshold for sequences
        
    Attributes
    ----------
    processed_sequences : list
        List of processed sequences
    stats : dict
        Statistics on processed sequences
        
    Notes
    -----
    This class is thread-safe and can be used in a multiprocessing context.
    """
    

    def __init__(self, min_length: int = 100, max_length: int = 1000, quality_threshold: float = 0.9) -> None:
        self.min_length = min_length
        self.max_length = max_length
        self.quality_threshold = quality_threshold
        self.processed_sequences: List[str] = []
        self.stats: Dict[str, Union[int, float]] = {}

    def process_batch(self, sequences: List[str]) -> List[Dict[str, Union[str, float]]]:
        """
        Process a batch of sequences.
        
        Parameters
        ----------
        sequences : list of str
            List of sequences to process
            
        Returns
        -------
        list of dict
            List of processed sequences with their metadata
            
        Raises
        ------
        ValueError
            If a sequence contains invalid characters
        
        Examples
        --------
        >>> processor = SequenceProcessor()
        >>> result = processor.process_batch(["ACGT", "TACG"])
        >>> len(result)
        2
        """
        pass
```

## Template for a function

```python
def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content of a sequence.
    
    Parameters
    ----------
    sequence : str
        DNA sequence
        
    Returns
    -------
    float
        Proportion of G and C in the sequence (between 0 and 1)
        
    Raises
    ------
    ValueError
        If the sequence contains invalid characters
        
    Examples
    --------
    >>> calculate_gc_content("ACGT")
    0.5
    >>> calculate_gc_content("AAAA")
    0.0
    """
    pass
```

## Template for a custom exception

```python
class SequenceQualityError(Exception):
    """
    Exception raised when a sequence does not meet quality criteria.
    
    Parameters
    ----------
    sequence_id : str
        ID of the problematic sequence
    issue : str
        Description of the quality issue
    """
    
    def __init__(self, sequence_id: str, issue: str) -> None:
        self.sequence_id = sequence_id
        self.issue = issue
        super().__init__(f"Quality issue for sequence {sequence_id}: {issue}")
```