# SCN2A: An Analysis of Protein Structure Deviations in SCN2A Mutations

This repository contains an analysis pipeline for investigating the hypothesis that mutations in the SCN2A gene causing larger deviations in protein structure (measured by RMSD and SASA) are more likely to be classified as pathogenic and associated with neurodevelopmental disorders such as epilepsy and autism.

## Table of Contents
- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Setup and Installation](#setup-and-installation)
- [How to Run the Analysis](#how-to-run-the-analysis)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Acknowledgements](#acknowledgements)

## Introduction

SCN2A encodes the sodium channel protein NaV1.2, which is critical for neural signaling. Mutations in this gene have been linked to a range of neurodevelopmental disorders. This project uses structural analysis and machine learning techniques to predict the pathogenicity of SCN2A mutations and their association with clinical outcomes.

## Project Structure

The repository is organized as follows:

```
SCN2A/
│
├── SCN2A__.py                         # Main analysis script
├── 6J8E.pdb                           # Template PDB file
├── wild_type.B99990001.pdb            # Generated wild-type PDB file
├── classification_association_confusion_matrix.png  # Confusion matrix for association analysis
├── classification_confusion_matrix.png              # Confusion matrix for classification analysis
├── regression_residuals.png                          # Regression residuals plot
├── regression_true_vs_predicted.png                 # Regression true vs predicted plot
├── Individual_Reflection.pdf         # Individual reflection document
├── Presentation_Script.pdf           # Presentation script
├── SCN2A_Project_Presentation_with_Detailed_Content.pptx  # Detailed presentation
├── SCN2A_Project_Report.pdf          # Final project report
├── README.md                         # This README file
```

## Requirements

This project requires the following dependencies:

- Python 3.x
- Modeller
- Biopython
- SciPy
- Matplotlib
- Seaborn
- Scikit-learn
- FreeSASA
- FPDF
- python-pptx

You can install the required packages using `pip`:

```bash
pip install biopython matplotlib seaborn scikit-learn freesasa fpdf python-pptx
```

**Note:** Modeller requires a license key. Make sure you have Modeller installed and properly configured. You can obtain a license key from [Modeller's official site](https://salilab.org/modeller/).

## Setup and Installation

1. Clone this repository to your local machine:

```bash
git clone https://github.com/hannahgwimpy/SCN2A.git
```

2. Navigate to the project directory:

```bash
cd SCN2A
```

3. Ensure Modeller is installed and properly configured. Add your Modeller license key to your environment configuration.

4. Download the necessary UniProt data for SCN2A by running the analysis script (details below).

## How to Run the Analysis

To perform the full analysis:

1. Run the main script:

```bash
python SCN2A__.py
```

The script will:
- Fetch SCN2A variant data from UniProt.
- Generate mutant protein sequences.
- Perform structural analysis, including calculating RMSD and SASA.
- Run regression and classification models.
- Generate visualizations, a project report, a presentation, and a reflection document.

The script is designed to be modular and will save intermediate results, plots, and reports in the project directory.

## Output Files

After running the script, the following files will be generated:

- **SCN2A_Project_Presentation_with_Detailed_Content.pptx**: The final presentation with detailed content.
- **SCN2A_Project_Report.pdf**: The final project report summarizing the analysis.
- **Individual_Reflection.pdf**: A reflection document.
- **Presentation_Script.pdf**: A script to accompany the presentation.
- Plots generated during the analysis: `classification_association_confusion_matrix.png`, `classification_confusion_matrix.png`, `regression_residuals.png`, and `regression_true_vs_predicted.png`.

## Troubleshooting

### Common Issues

1. **Modeller License Key Error**: Ensure Modeller is correctly installed and the license key is configured.
2. **Missing Dependencies**: Double-check that all required packages are installed using `pip install -r requirements.txt`.
3. **FileNotFoundError**: Ensure that the necessary PDB files (e.g., `6J8E.pdb`) are in the correct directory.

## Acknowledgements

This project was developed by Hannah Wimpy as part of an academic research initiative. It leverages open-source tools such as UniProt, Modeller, and SciPy for structural bioinformatics and data science.
