import numpy as np
import matplotlib.pyplot as plt

amino_acid_mass_dalton = {
    'A': 71.08, 'R': 156.19, 'N': 114.10, 'D': 115.09,
    'C': 103.15, 'E': 129.12, 'Q': 128.13, 'G': 57.05,
    'H': 137.14, 'I': 113.16, 'L': 113.16, 'K': 128.17,
    'M': 131.19, 'F': 147.18, 'P': 97.12, 'S': 87.08,
    'T': 101.11, 'W': 186.21, 'Y': 163.18, 'V': 99.13,
}

def calculate_mol_mass(peptide_seq, amino_acid_mass_dict):
    """
    Calculate the molecular mass of a peptide sequence.

    Calculates the sum of the molecular masses of all amino acids
    in the peptide sequence and accounts for the loss of water
    molecules during peptide bond formation.

    Parameters
    ----------
    peptide_seq : str
        Amino acid sequence of the peptide (one-letter code).
    amino_acid_mass_dict : dict
        Dictionary mapping amino acid one-letter codes to their
        molecular masses (in Da).

    Returns
    -------
    float
        Calculated molecular mass of the peptide in Da.
    """
    
    molmass_dic = {}
    molm = sum(amino_acid_mass_dict[aa] for aa in peptide_seq)
    molmass_dic[peptide_seq]=molm
    return molmass_dic

def calculate_mol_mass_collection(peptides, amino_acid_mass_dict):
    """
    Calculate the molecular mass of each peptide based on its amino acid sequence.

    Parameters
    ----------
    peptides : list
        List of peptide sequences for which the molecular masses are calculated.
    amino_acid_mass_dict : dict
        Dictionary containing amino acids as keys and their corresponding masses as values.

    Returns
    -------
    molmass_dic : dict
        Dictionary containing peptide sequences as keys and their calculated molecular masses as values.
    """
    molmass_dic = {}
    for peptide_seq in peptides:
        molm = sum(amino_acid_mass_dict[aa] for aa in peptide_seq)
        molmass_dic[peptide_seq]=molm
    return molmass_dic

def calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007):
    """
    Calculate the mass-to-charge ratio (m/z) for each peptide.

    Parameters
    ----------
    peptide_mass_map : dict
        Dictionary containing peptide sequences and their corresponding molecular masses.
    charge : int
        Charge state of the peptides.
    proton_mass : float
        Mass of a proton.

    Returns
    -------
    peptides_mz_dic : dict
        Dictionary containing peptide sequences as keys and their calculated m/z values as values.
    """
    peptides_mz_dic = {}
    for peptide, mass in peptide_mass_map.items():
        peptides_mz_dic[peptide]= (mass + charge * proton_mass) / charge
    return peptides_mz_dic

def plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42):
    """
    Plot a mass spectrum using m/z values and randomly generated intensities.

    Parameters
    ----------
    mz_values : dict
        Dictionary containing peptide sequences and their corresponding m/z values.
    random_count_range : tuple, optional
        Minimum and maximum values for randomly generated intensities.
        Default is (0, 30000).
    seed : int, optional
        Random seed to ensure reproducibility of the generated intensities.
        Default is 42.

    Returns
    -------
    bar_plot : BarContainer
        Matplotlib bar chart representing the mass spectrum.
    """
    np.random.seed(seed)
    mz_list = []
    intensity = []
    for mz in mz_values.values():
        mz_list.append(mz)
        intensity.append(np.random.randint(random_count_range[0], random_count_range[1]+1))
    
    bar_plot = plt.bar(mz_list, intensity, width = 1)
    plt.title("Mass Spectrum")
    plt.xlabel("m/z")
    plt.ylabel("Intensity")
    return(bar_plot)

def fragment_peptide(peptide):
    """
    Generate b- and y-ion fragments for peptide sequences.
    """

    # Allow both a single peptide string and a list of peptides
    if isinstance(peptide, str):
        peptide = [peptide]

    fragment_list_b = []
    fragment_list_y = []

    for pep in peptide:

        # b-ion fragments: N-terminus to C-terminus
        for i in range(1, len(pep)):
            fragment_list_b.append(pep[:i])

        # y-ion fragments: C-terminus to N-terminus
        for i in range(1, len(pep)):
            fragment_list_y.append(pep[len(pep)-i:])

    return fragment_list_b + fragment_list_y
      