# GPCR Selectivity Explorer

A command-line tool designed to analyze and compare amino acid conservation and selectivity between two groups of G-protein coupled receptors (GPCRs) based on data from GPCRdb.

This is a command-line version of server version running at https://carlssonlabtools.icm.uu.se/GPCR_Selectivity_Explorer

## Features

*   Compare two user-defined lists of receptors.
*   Focus analysis on specific structural segments (e.g., transmembrane helices, loops).
*   Set a conservation cutoff to identify key residues.
*   Utilize various substitution matrices (e.g., BLOSUM62, PAM250) for analysis.
*   Supports custom, user-modified alignments and residue tables for advanced analysis.


## Requirements
- Operating system: Linux, macOS, Windows
- Python: 3.8+
- Memory: ≥ 4 GB RAM recommended
- Disk: ~500 MB free space
- Conda or Miniconda
- Required Python packages (pandas, biopython, requests, openpyxl, pymol-open-source).

## Install Miniconda

```bash
# --- For Linux ---
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# --- For macOS ---
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
# bash Miniconda3-latest-MacOSX-x86_64.sh

# Initialize Conda
source ~/.bashrc  # or 'source ~/.zshrc' depending on your shell

## Installation
Clone the repository and install dependencies, it may takes a few minutes:
```bash
git clone https://github.com/Carlssonlab/GPCR_Selectivity_Explorer.git
cd GPCR_Selectivity_Explorer

# Create a conda environment named 'venv' with Python 3.10
conda create -n venv python=3.10 -y

# Activate the environment
conda activate venv

# Install dependencies listed in requirements.txt
pip install -r requirements.txt

# Install PyMOL separately (not available on PyPI)
conda install -c conda-forge pymol-open-source -y
## Usage

The script is run from the command line, providing arguments to define the receptor groups, analysis parameters, and output location.

```bash
python GPCR_Selectivity_Explorer.py [-h] [-l1 LIST1 [LIST1 ...]] [-r1 REFERENCE1] [-l2 LIST2 [LIST2 ...]] [-r2 REFERENCE2] [-sg SEGMENTS [SEGMENTS ...]] [-c CONSERVATION_CUTOFF] [-sb SUBSTITUTION_MATRIX_METHOD] [-o OUTPUTFOLDER] [--custom] [-r1f REFERENCE1_RESIDUES_TABLE] [-r1A REFERENCE1_ALIGNMENT] [-r2f REFERENCE2_RESIDUES_TABLE] [-r2A REFERENCE2_ALIGNMENT]
```

### Arguments

| Argument | Description |
| :--- | :--- |
| `-h`, `--help` | Show this help message and exit. |
| `-l1`, `--list1` | Space-separated list for Reference 1 Receptor segments (e.g., `adrb1_human adrb2_human adrb3_human`). |
| `-r1`, `--reference1` | Reference 1 Receptor Name in GPCRdb (e.g., `adrb2_human`). |
| `-l2`, `--list2` | Space-separated list for Reference 2 Receptor segments (e.g., `drd1_human drd5_human`). |
| `-r2`, `--reference2` | Reference 2 Receptor Name in GPCRdb (e.g., `drd1_human`). |
| `-sg`, `--segments` | Segment information for reference receptor residues from GPCRdb (e.g., `TM1 TM2 ECL2 ICL1`). |
| `-c`, `--conservation_cutoff` | Conservation Cutoff value (e.g., `80`). |
| `-sb`, `--substitution_matrix_method` | Substitution matrix method. Options include: `BENNER22`, `BENNER6`, `BENNER74`, `BLOSUM45`, `BLOSUM50`, `BLOSUM62`, `BLOSUM80`, `BLOSUM90`, `DAYHOFF`, `FENG`, `GENETIC`, `GONNET1992`, `HOXD70`, `JOHNSON`, `JONES`, `LEVIN`, `MCLACHLAN`, `MDM78`, `NUC.4.4`, `PAM250`, `PAM30`, `PAM70`, `RAO`, `RISLER`, `SCHNEIDER`, `STR`, `TRANS`. |
| `-o`, `--outputfolder` | The path to the output folder. |
| `--custom` | A flag to indicate the use of custom alignments and residue tables. |
| `-r1f`, `--reference1_residues_table` | Path to the Reference 1 modified residues table Excel file from GPCRdb. |
| `-r1A`, `--reference1_alignment` | Path to the Reference 1 modified alignment Excel file from GPCRdb. |
| `-r2f`, `--reference2_residues_table` | Path to the Reference 2 modified residues table Excel file from GPCRdb. |
| `-r2A`, `--reference2_alignment` | Path to the Reference 2 modified alignment Excel file from GPCRdb. |


### Outputs

- Sequence alignments, mutagenesis files, and residue numbering tables are provided in Excel format (.xlsx).
- 3D receptor structures are delivered as PDB files (.pdb), along with a ready-to-use PyMOL visualization session (.pse).

The PDB files include additional annotations encoded in the B-factor field:
- GPCRdb numbering is stored in the B-factor field of the C atoms in the backbone.
- Consensus conservation values are stored in the B-factor field of the N atoms in the backbone.
- BLOSUM-based substitution classification is encoded in the B-factor field of the Cα atoms in the backbone using the following scheme:
    - Conservative substitution, above conservation cutoff in one receptor → B-factor = 12.5 (cyan in PyMOL)
    - Conservative substitution, above conservation cutoff in both receptors → B-factor = 25 (dark blue in PyMOL)
    - Non-conservative substitution, above conservation cutoff in one receptor → B-factor = 50 (orange in PyMOL)
    - Non-conservative substitution, above conservation cutoff in both receptors → B-factor = 100 (red in PyMOL)


## Examples

### 1. Standard Analysis

This example compares a list of adrenergic receptors against a list of dopamine receptors, focusing on the transmembrane helices and ECL2. It uses a conservation cutoff of 80% and the BLOSUM62 matrix, the run time is about few minutes.

```bash
python GPCR_Selectivity_Explorer.py \
-l1 adrb1_human adrb2_human adrb3_human \
-r1 adrb2_human \
-l2 drd1_human drd5_human \
-r2 drd1_human \
-sg TM1 TM2 TM3 TM4 TM5 TM6 TM7 ECL2 \
-c 80 \
-sb BLOSUM62 \
-o outputs_default
```

### 2. Using Custom Alignments

This example demonstrates how to use manually edited alignment and residue table files as input.

**Workflow for Custom Alignments:**

To use custom alignments, first run the tool in standard mode. Then, edit the output alignment files (`*_GPCRdb_alignment.xlsx`) for List 1 and List 2 by adding or removing sequences. If you add new positions (e.g., the full ECL3), you must also update the corresponding reference residue table (`*_GPCRdb_residues_table_full_seq.xlsx`) by adding the amino acid number and its GPCRdb number (e.g., `E193 67x10`). If a GPCRdb number does not exist for a new position, create a unique one that is not already assigned.

Finally, use the custom alignment and residue tables as inputs for the workflow using the `--custom` flag.

```bash
python GPCR_Selectivity_Explorer.py \
-r1 adrb2_human \
-r2 drd1_human \
-c 80 \
-r1f outputs_default/adrb2_human_GPCRdb_residues_table_full_seq.xlsx \
-r1A outputs_default/adrb2_human_GPCRdb_alignment.xlsx \
-r2f outputs_default/drd1_human_GPCRdb_residues_table_full_seq.xlsx \
-r2A outputs_default/drd1_human_GPCRdb_alignment.xlsx \
--custom \
-o outputs_custom
```

## Reproducibility
To reproduce the results of β adrenergic and dopamine D1-like receptors described in work, you can run:
```bash
cd GPCR_Selectivity_Explorer/ADRBs_D1R-like_analysis_repreducability
bash run_analysis.sh
```

## Citing this work
If you use GPCR Selectivity Explorer in your research, please cite the following publication:

Kahlous NA, Rinne M, Zhang X, Li Y, Chen Y, Motso A, Gao K, Bergqvist C, Shen H, Wang Y, De Vaca I, Díaz-Holguín A, Ullmann P, Bengtsson T, Lauschke VM, Kukkonen JP, Delemotte L, Wright SC, Liu X, Larhammar D, Carlsson J, “[Molecular mechanisms of ligand selectivity by catecholamine G protein-coupled receptors],” Nature communications, 2025, Volume(Issue), Pages.  DOI: [Paper DOI]

## License

This project is licensed under the MIT License.
