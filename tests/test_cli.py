import pytest
from click.testing import CliRunner
from decaf.cli import main


@pytest.fixture
def fake_fasta(tmpdir):
    """Crée un fichier FASTA fictif temporaire pour les tests."""
    fake_fasta_file = tmpdir.join("input.fa")
    fake_fasta_file.write(">seq1\nAGCTAGCTAG\n>seq2\nTGCATGCATG\n")
    return str(fake_fasta_file)


def test_cli(fake_fasta):
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            '--barcode', 'ITS',
            '--taxa', 'plants',
            '--input_fastq', fake_fasta,
            '--output_folder', 'output/'
        ]
    )
    
    # Vérification de la sortie et du succès de l'exécution
    assert result.exit_code == 0
    assert "Starting DECAF processing with ITS barcode for plants" in result.output
