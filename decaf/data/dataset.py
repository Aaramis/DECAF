from typing import Any, Dict, List

from Bio.SeqRecord import SeqRecord
from torch.utils.data import Dataset
from transformers import PreTrainedTokenizer


class AmpliconDataset(Dataset):
    """
    Dataset for amplicon sequences used in transformer models.

    This class provides a PyTorch Dataset implementation for biological
    amplicon sequences, handling tokenization and preparation for model input.

    Parameters
    ----------
    sequences : list of SeqRecord
        List of BioPython SeqRecord objects (from FASTQ/FASTA)
    tokenizer : PreTrainedTokenizer
        A HuggingFace tokenizer to convert sequences into tokens
    config : dict
        Configuration dictionary containing parameters like max_length

    Attributes
    ----------
    sequences : list
        The list of sequence records
    tokenizer : PreTrainedTokenizer
        The tokenizer instance
    max_length : int
        Maximum sequence length for padding/truncation
    """

    def __init__(
        self,
        sequences: List[SeqRecord],
        tokenizer: PreTrainedTokenizer,
        config: Dict[str, Any],
    ):
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.max_length = config["max_length"]

    def __len__(self) -> int:
        """
        Get the number of sequences in the dataset.

        Returns
        -------
        int
            Number of sequences

        Examples
        --------
        >>> dataset = AmpliconDataset(sequences, tokenizer, config)
        >>> len(dataset)
        150
        """
        return len(self.sequences)

    def __getitem__(self, idx: int) -> Dict[str, Any]:
        """
        Get a tokenized sequence by index.

        Parameters
        ----------
        idx : int
            Index of the sequence to retrieve

        Returns
        -------
        dict
            Dictionary containing:
            - input_ids: Tensor of token IDs
            - attention_mask: Tensor indicating non-padded tokens
            - sequence_id: Original sequence identifier
            - seq_quality: Original sequence quality
            - sequence_str: Original sequence string

        Examples
        --------
        >>> dataset = AmpliconDataset(sequences, tokenizer, config)
        >>> item = dataset[0]
        >>> item["input_ids"].shape
        torch.Size([512])
        """
        seq_record = self.sequences[idx]
        seq_str = str(seq_record.seq)
        if "phred_quality" in seq_record.letter_annotations:
            phred_scores = seq_record.letter_annotations["phred_quality"]
            seq_quality = "".join(chr(q + 33) for q in phred_scores)
        else:
            seq_quality = ""

        # Tokenization
        encoded = self.tokenizer(
            seq_str,
            padding="max_length",
            truncation=True,
            max_length=self.max_length,
            return_tensors="pt",
        )

        # Remove batch dimension
        input_ids = encoded["input_ids"].squeeze(0)
        attention_mask = encoded["attention_mask"].squeeze(0)

        return {
            "input_ids": input_ids,
            "attention_mask": attention_mask,
            "sequence_id": seq_record.id,
            "seq_quality": seq_quality,
            "sequence_str": seq_str,
        }
