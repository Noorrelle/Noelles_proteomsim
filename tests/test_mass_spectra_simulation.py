from Noelles_proteomsim import (
    amino_acid_mass_dalton,
    calculate_mol_mass,
    calculate_mol_mass_collection,
    calculate_mz_collection,
    plot_spectrum,
    fragment_peptide,
)

def test_calculate_mol_mass():
    aa_mass_dict = {"A":1, "B":2, "C":3, "D":4}
    test_peptide = "ABCD"
    expected_mass = {"ABCD":1+2+3+4}
    calculated_mass = calculate_mol_mass(test_peptide, aa_mass_dict)
    assert calculated_mass == expected_mass
    # was passiert wenn etwas nicht in seinem dictionary vorhanden ist 

test_calculate_mol_mass()

def test_calculate_mol_mass_collection():
    aa_mass_dict = {"A":1, "B":2, "C":3, "D":4}
    peptides = ["ABCD", "ABCD"]
    expected = {"ABCD": 10, "ABCD":10}
    actual = calculate_mol_mass_collection(peptides, aa_mass_dict)
    peptides2 = ["ABCD"]
    expected2 = {"ABCD": 10}
    actual2 = calculate_mol_mass_collection(peptides2, aa_mass_dict)

    assert actual == expected
    assert actual2 == expected2

test_calculate_mol_mass_collection()

def test_calculate_mz_collection():
    peptide_mass_map = {"PEPTIDE": 799.36,
        "TEST": 450.25,
        "MASS": 600.50}
    actual = calculate_mz_collection(peptide_mass_map, charge=2)
    #expected = {"PEPTIDE": 400.687, "TEST":226.132, "MASS":301.257}
    expected = {"PEPTIDE":(799.36 + 2 * 1.007) / 2, "TEST":(450.25 + 2 * 1.007) / 2, "MASS":(600.50 + 2 * 1.007) / 2}

    assert actual == expected

test_calculate_mz_collection()

def test_fragment_peptide():
    peptide = "PEPT"
    expected = fragment_peptide(peptide)
    actual = ["P", "PE", "PEP", "T", "PT", "EPT"]

    print(fragment_peptide("PEPT"))
    print(actual)

    assert set(actual) == set(expected)
    assert "PEPT" not in expected

test_fragment_peptide()