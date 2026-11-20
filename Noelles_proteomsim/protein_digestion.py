import re

enzyme_cleavage_patterns = {
    'LysC': r'(?<=K)',
    'LysN': r'(?=K)',
    'ArgC': r'(?<=R)',
    'Trypsin': r'(?<=[KR])(?!P)',
}

def digest_protein_sequence(protein_seq, enzyme, min_pep_len=5, max_pep_len=30):
    """
    function that accepts a protein sequence plus a cleavage pattern and returns the resulting list of peptides and filters paptides within a certain range

    Parameters
    protein_seq : str
        sequence of your protein
    cleave_pattern : str
        Enzyme beaviour that is mimicked
    max_pep_len : float
        the max length of each peptide, default = 30
    min_pep_len
        the min length of each peptide, default = 5
    Returns
    a list with all peptides that are within in the length range you put in
    """
    sequence = protein_seq
    cleave_pattern = enzyme_cleavage_patterns[enzyme]
    peptides = re.split(cleave_pattern, sequence)

    filtered = [
        seq for seq in peptides 
            if len(seq) >= min_pep_len and len(seq) <= max_pep_len
        ] # list comprehension for filtering

    return filtered
    
def digest_protein_collection(protein_map, enzyme, min_pep_len=5, max_pep_len=30):
    """ 
    function that loops over every protein in the dictionary and cleaves each protein by the selected pattern and returns the resulting list of peptides and filters pwptides within a certain range.

    Parameters
    ----------
    protein_seq : str
        sequence of your protein
    cleave_pattern : str
        Enzyme beaviour that is mimicked
    max_pep_len : float
        the max length of each peptide
    min_pep_len
        the min length of each peptide
    Returns
    -------
    a list with all peptides that are within in the length range you put in
    """
    peptide_dictionary = {}
    for protein, sequence in protein_map.items():
        peptide_dictionary[protein] = digest_protein_sequence(sequence, enzyme, min_pep_len, max_pep_len)
        
        
    return peptide_dictionary

def compute_sequence_coverage(protein_seq, peptides):
    """
    Compute the sequence coverage of a protein based on its identified peptides.

    The sequence coverage is calculated as the sum of the lengths of all
    identified peptides divided by the length of the protein sequence.

    Parameters
    ----------
    protein_seq : str
        Amino acid sequence of the protein.
    peptides : list of str
        List of identified peptide sequences belonging to the protein.

    Returns
    -------
    float
        Sequence coverage of the protein, expressed as a fraction between
        0 and 1.
    """
    
    len_pro = len(protein_seq)
    len_pepall = 0
    for pep in peptides:
        len_pep = len(pep)
        len_pepall += len_pep
    
    coverage = len_pepall / len_pro 
    return coverage

def compute_sequence_coverage_dictionary(peptide_dictionary, protein_dictionary):
    """
    Compute the sequence coverage for each protein in a dictionary.

    The sequence coverage is calculated as the sum of the lengths of all
    peptides associated with a protein divided by the length of the
    corresponding protein sequence.

    Parameters
    ----------
    peptide_dictionary : dict
        Dictionary containing protein identifiers as keys and lists of
        peptide sequences as values.
    protein_dictionary : dict
        Dictionary containing protein identifiers as keys and their
        corresponding amino acid sequences as values.

    Returns
    -------
    dict
        Dictionary containing protein identifiers as keys and their
        sequence coverage values as strings.
    """
    
    coverage_dictionary = {}
    
    for protein, sequence in peptide_dictionary.items():
        n_pro = len(protein_dictionary[protein])
        n_allpep = 0
        for seq in sequence:
            
            n_pep = len(seq)
            n_allpep = n_pep + n_allpep
        
        coverage = n_allpep / n_pro
        coverage_dictionary[protein] = str(coverage)
    return coverage_dictionary
