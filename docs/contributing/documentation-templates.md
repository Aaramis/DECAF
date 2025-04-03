# Docstrings

Ce fichier contient des templates pour les différents types de docstrings à utiliser dans DECAF.
Ces exemples servent de référence pour maintenir une documentation cohérente.


## Template pour un module
```python
"""
Module pour le traitement des séquences ITS.

Ce module fournit des fonctions pour nettoyer, normaliser et préparer
des séquences ITS pour l'analyse et la classification.

Functions
---------
clean_sequence
    Nettoie une séquence en supprimant les caractères non-nucléotidiques.
normalize_sequence
    Normalise une séquence pour l'entraînement du modèle.
extract_features
    Extrait des caractéristiques à partir d'une séquence.
"""
```


## Template pour une classe

```python
class SequenceProcessor:
    """
    Classe pour le traitement des séquences biologiques.
    
    Cette classe implémente diverses méthodes pour nettoyer, normaliser
    et extraire des informations des séquences biologiques, en particulier
    des séquences ITS de plantes.
    
    Parameters
    ----------
    min_length : int, default=100
        Longueur minimale des séquences à traiter
    max_length : int, default=1000
        Longueur maximale des séquences à traiter
    quality_threshold : float, default=0.9
        Seuil de qualité minimum pour les séquences
        
    Attributes
    ----------
    processed_sequences : list
        Liste des séquences traitées
    stats : dict
        Statistiques sur les séquences traitées
        
    Notes
    -----
    Cette classe est thread-safe et peut être utilisée dans un
    contexte multiprocessing.
    """
    

    def __init__(self, min_length: int = 100, max_length: int = 1000, quality_threshold: float = 0.9) -> None:
        self.min_length = min_length
        self.max_length = max_length
        self.quality_threshold = quality_threshold
        self.processed_sequences: List[str] = []
        self.stats: Dict[str, Union[int, float]] = {}

    def process_batch(self, sequences: List[str]) -> List[Dict[str, Union[str, float]]]:
        """
        Traite un lot de séquences.
        
        Parameters
        ----------
        sequences : list of str
            Liste des séquences à traiter
            
        Returns
        -------
        list of dict
            Liste des séquences traitées avec leurs métadonnées
            
        Raises
        ------
        ValueError
            Si une séquence contient des caractères non valides
        
        Examples
        --------
        >>> processor = SequenceProcessor()
        >>> result = processor.process_batch(["ACGT", "TACG"])
        >>> len(result)
        2
        """
        pass
```

## Template pour une fonction

```python
def calculate_gc_content(sequence: str) -> float:
    """
    Calcule le contenu en GC d'une séquence.
    
    Parameters
    ----------
    sequence : str
        Séquence d'ADN
        
    Returns
    -------
    float
        Proportion de G et C dans la séquence (entre 0 et 1)
        
    Raises
    ------
    ValueError
        Si la séquence contient des caractères non-ATGC
        
    Examples
    --------
    >>> calculate_gc_content("ACGT")
    0.5
    >>> calculate_gc_content("AAAA")
    0.0
    """
    pass
```

## Template pour une exception personnalisée

```python
class SequenceQualityError(Exception):
    """
    Exception levée lorsqu'une séquence ne répond pas aux critères de qualité.
    
    Parameters
    ----------
    sequence_id : str
        Identifiant de la séquence problématique
    issue : str
        Description du problème de qualité
    """
    
    def __init__(self, sequence_id: str, issue: str) -> None:
        self.sequence_id = sequence_id
        self.issue = issue
        super().__init__(f"Problème de qualité pour la séquence {sequence_id}: {issue}")
```