from Noelles_proteomsim import (
       enzyme_cleavage_patterns,
       digest_protein_sequence,
       digest_protein_collection,
       compute_sequence_coverage,
       compute_sequence_coverage_dictionary,
   )

def test_digest_protein_sequence():
    dummy_seq_1 = "RQFLRGLDFLHANCIVHRDLKPENILVTSGGTVKLADFGLARIYSYQMALTPVVVTLWYRAPEVLLQSTYATPVDMWSVGCIFAEMFRRKPLFCGNSEADQLGKIFDLIGLPPEDDWPR"
    dummy_seq_2 = "MEDEEKAVEILGNTEAAHPPSPIRCCWLRLRCLAATSIICGCSCLGVMALVFAIKAEERHKAGRSEEAVRWGARARKLILASFAVWLAVLILGPLLLWLLSYAIAQAPEPTIDERE"
    test1 = digest_protein_sequence(dummy_seq_1, "LysN", 1, 20)
    test2 = digest_protein_sequence(dummy_seq_2, "LysN", 1, 20)
    assert len(test1) == 4
    assert len(test2) == 3
    

test_digest_protein_sequence()

def test_digest_protein_collection():
    dummy_proteins = {
        "Protein_1": "MKKLAAARQFLRGLDFLHANCIVHRDLKPENILVTSGGTVKL",
        "Protein_2": "MEDEEKAVEILGNTEAAHPPSPIRCCWLRLRCLAATSIICGCS",
        "Protein_3": "MALWMRLLPLLALLALWGPGPGQGQGQGQKLVFFAEDVGSNK",
        "Protein_4": "MKWVTFISLLFLFSSAYS",
        "Protein_5": "GILGFVFTLTVPSERLEKAKGLL"
    }
    assert digest_protein_collection(dummy_proteins, "LysN") == {'Protein_1': ['KLAAARQFLRGLDFLHANCIVHRDL', 'KPENILVTSGGTV'],
    'Protein_2': ['MEDEE'],
    'Protein_3': ['MALWMRLLPLLALLALWGPGPGQGQGQGQ', 'KLVFFAEDVGSN'],
    'Protein_4': ['KWVTFISLLFLFSSAYS'],
    'Protein_5': ['GILGFVFTLTVPSERLE']}

test_digest_protein_collection()

def test_compute_sequence_coverage():
    dummy_prot_seq = "MKKLAAARQFLRGLDFLHANCIVHRDLKPENILVTSGGTVKL"
    dummy_peps = ['KLAAARQFLRGLDFLHANCIVHRDL', 'KPENILVTSGGTV']
    coverage = 0.9047619047619048
    assert coverage == compute_sequence_coverage(dummy_prot_seq, dummy_peps)

test_compute_sequence_coverage()