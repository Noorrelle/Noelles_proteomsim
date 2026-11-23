# Noelles_proteomsim

A small Python package that simulates a bottom-up mass spectrometry proteomics experiment: from protein sequences in a FASTA file, through enzymatic digestion and liquid chromatography, to MS1 and MS2 spectra.
It was developed as part of the *Advanced Python* course and is still under active development.

## What it does

| Step | Input | Output |
|---|---|---|
| 1. Read proteins | FASTA file (UniProt-style headers) | dictionary `{protein_id: sequence}` |
| 2. Protein digestion | protein dictionary + protease name (`"Trypsin"`, `"LysC"`, `"LysN"`, `"ArgC"`) | dictionary `{protein_id: [peptides]}`, sequence coverage per protein |
| 3. Liquid chromatography | list of peptides | dictionary `{peptide: retention_time}`, retention time histogram, peptides inside a chosen retention time window |
| 4. Mass spectrometry (MS1) | selected peptides | dictionaries of peptide masses and m/z values, MS1 spectrum plot |
| 5. Fragmentation (MS2) | one peptide | list of b- and y-ion fragments |

## Installation

Requires Python 3.10 or newer.

```bash
# 1. Clone the repository
git clone https://github.com/Noorrelle/Noelles_proteomsim.git
cd Noelles_proteomsim

# 2. Create and activate a fresh environment
conda create -n proteomsim python=3.10
conda activate proteomsim

# 3. Install the exact dependency versions, then the package itself
pip install -r requirements.txt
pip install .
```

Alternatively, install the package directly from GitHub:

```bash
pip install git+https://github.com/Noorrelle/Noelles_proteomsim.git
```

## Quick start

The easiest way to try the package is the example notebook
**[`tutorials/ms_experiment_final.ipynb`](tutorials/ms_experiment_final.ipynb)**, which runs the whole workflow end to end on [`data/sample_proteins.fasta`](data/sample_proteins.fasta).

The same workflow in a few lines (run from the repository root):

```python
import Noelles_proteomsim as ps

# 1. Read proteins
proteins = ps.read_fasta("data/sample_proteins.fasta")

# 2. Digest with trypsin, keep peptides with 5-30 amino acids
peptides = ps.digest_protein_collection(proteins, "Trypsin", 5, 30)
coverage = ps.compute_sequence_coverage_dictionary(peptides, proteins)

# 3. Liquid chromatography
pooled_peptides = [pep for peps in peptides.values() for pep in peps]
rt_map = ps.predict_lc_retention_times(pooled_peptides)
ps.plot_retention_time(rt_map.values())
selected = ps.select_retention_time_window(rt_map, 3, 20)

# 4. MS1: masses, m/z values (charge 2) and spectrum
masses = ps.calculate_mol_mass_collection(selected, ps.amino_acid_mass_dalton)
mz_values = ps.calculate_mz_collection(masses, 2)
ps.plot_spectrum(mz_values)

# 5. MS2: fragment one peptide
fragments = ps.fragment_peptide("MATSR")
```

## Package overview

### `file_handling`
- `read_fasta(filepath)` – reads a FASTA file into a dictionary. The protein ID is taken from the second field of a UniProt header (`>sp|P11802|CDK4_HUMAN ...` → `P11802`).

### `protein_digestion`
- `enzyme_cleavage_patterns` – regular expressions describing where each protease cuts (`LysC`, `LysN`, `ArgC`, `Trypsin`).
- `digest_protein_sequence(protein_seq, enzyme, min_pep_len=5, max_pep_len=30)` – digests one protein sequence and returns the peptides within the length range.
- `digest_protein_collection(protein_map, enzyme, min_pep_len=5, max_pep_len=30)` – digests every protein of a dictionary.
- `compute_sequence_coverage(protein_seq, peptides)` – fraction of one protein sequence covered by its peptides.
- `compute_sequence_coverage_dictionary(peptide_dictionary, protein_dictionary)` – sequence coverage for every protein.

### `liquid_chromatography`
- `predict_lc_retention_times(peptides)` – predicts retention times with pyteomics `achrom` (Guo pH 7.0 coefficients); returns `{peptide: retention_time}`.
- `plot_retention_time(retention_times, resolution=30)` – histogram of the retention times.
- `select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time)` – keeps only the peptides eluting inside the window.

### `mass_spectra_simulation`
- `amino_acid_mass_dalton` – average amino acid residue masses in Dalton.
- `calculate_mol_mass(peptide_seq, amino_acid_mass_dict)` – mass of one peptide, returned as `{peptide: mass}`.
- `calculate_mol_mass_collection(peptides, amino_acid_mass_dict)` – masses of many peptides.
- `calculate_mz_collection(peptide_mass_map, charge, proton_mass=1.007)` – m/z values for a given charge state.
- `plot_spectrum(mz_values, random_count_range=(0, 30000), seed=42)` – plots a spectrum with simulated random intensities.
- `fragment_peptide(peptide)` – returns the b-ion and y-ion fragment sequences of a peptide (or a list of peptides).

## Simplifications

This is a teaching simulation, not an analysis tool:
- Peptide masses are the sum of residue masses; terminal water is not added.
- Peak intensities in `plot_spectrum` are random numbers (reproducible through `seed`).
- Fragments are returned as sequences; their masses can be computed with `calculate_mol_mass_collection`.

## Running the tests

```bash
pip install pytest
pytest -v tests
```

The tests also run automatically on GitHub Actions for every push and pull request to `main`.

## License

MIT, see [LICENSE](LICENSE).
