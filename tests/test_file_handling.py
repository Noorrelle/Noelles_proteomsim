from pathlib import Path

from Noelles_proteomsim import read_fasta


def test_read_fasta():
    fasta_path = Path(__file__).parent / "dummy.fasta"
    protein_map = read_fasta(str(fasta_path))

    assert protein_map["P12345"] == "MKWVTFISLLFLFSSAYSGSLRLSFSGK"
    assert protein_map["P67890"] == "MALWMRLLPLLALLALWGPGPGSGK"

