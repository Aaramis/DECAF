from typing import Any, Dict, List

from Bio.SeqRecord import SeqRecord
from torch.utils.data import Dataset
from transformers import PreTrainedTokenizer


class AmpliconDataset(Dataset):
    def __init__(
        self,
        sequences: List[SeqRecord],
        tokenizer: PreTrainedTokenizer,
        config: Dict[str, Any],
    ):
        """
        Args:
            sequences: List of BioPython SeqRecord objects (from FASTQ/FASTA)
            tokenizer: A HuggingFace tokenizer
        """
        self.sequences = sequences
        self.tokenizer = tokenizer
        self.max_length = config["max_length"]

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        seq_record = self.sequences[idx]
        seq_str = str(seq_record.seq)  # Convert Seq to string

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
            "sequence_str": seq_str,
        }
