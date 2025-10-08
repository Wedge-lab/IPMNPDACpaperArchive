# IPMNPDAC_WGS
Scripts for analyzing IPMN-PDAC multi-region WGS data and reconstructing clonal trees.

## ipynb View

### 1) Install VS Code: 
https://code.visualstudio.com/

### 2) Paste to URL: 
https://github.com/xtgithubhe/IPMNPDACpaperArchive/blob/main/IPMNPDAC_WGS/ipmnpdac_WGS.ipynb
Press Enter. VS Code will open the repository in a temporary (virtual) workspace — no need to clone it locally.

**OR:**
### 3) View the notebooks directly using:
- nbviewer
  
  https://nbviewer.org/github.com/xtgithubhe/IPMNPDACpaperArchive/blob/main/IPMNPDAC_WGS/ipmnpdac_WGS.ipynb
- Google Colab
  
  https://colab.research.google.com/github.com/xtgithubhe/IPMNPDACpaperArchive/blob/main/IPMNPDAC_WGS/ipmnpdac_WGS.ipynb

## Requirements

- **Python**: >= 3.9  
- **R**: >= 4.2  
  ### Python dependencies
  - pandas>=2.2.3
  - numpy>=2.2.6
  - matplotlib>=3.10.0
  ### R dependencies
  - dplyr(>=1.1.4)       
  - ggplot2(>=3.5.2)     
  - tibble(>=3.3.0)
  - tidyr(>=1.3.1)
    
  ## Input Data Formats
  1. CSV file
     - Example: case16_pindelSNVmsDPC.csv
     - Required columns:
       * chr
       * pos
       * most.likely.cluster
  2. Text table
     - Example: *.DPinput.txt
     - Required columns:
       * chr
       * end
       * subclonal.CN
       * nMaj1
       * nMin1
       * frac1
       * nMaj2
       * nMin2
       * frac2
       * no.chrs.bearing.mut
           
  ## Data Use
 The path of all datasets can be adjusted to access the Data within the Git repository structure, but the file names must remain the same.

**Contact:** xiaotong.he@manchester.ac.uk

# BBDPCtiming
Code for the BBDtiming approach to mutation timing and clonal inference using inputs from Battenberg and DPClust.

## Overview
To infer the likely timing and clonal status of mutations within the tumour, we categorized them into six groups (from clonal early to subclone) using the BBD timing approach, 
which integrates cancer cell fraction (CCF) estimates and local copy number information derived from Battenberg and DPClust ( https://github.com/Wedge-lab ). 

Mutations with a CCF ≤ 0 were designated as unspecialized (unSp), indicating likely sequencing artefacts or biologically insignificant variants with minimal representation in the cancer cell population. 
Mutations with a CCF between 0 and 0.95 were classified as subclone, representing variants likely emerging after the initial clonal expansion.
Mutations with a CCF of at least 0.95, indicative of high cancer cell prevalence, were further subdivided based on copy number features and chromosomal distribution. 
If either of the major alleles (nMaj1 or nMaj2) had a copy number greater than 1 and the mutation was present on more than one chromosome (noChrsBearingMut > 1), 
the mutation was classified as clone early (cloneEarly), suggesting early clonal mutations that underwent chromosomal duplication. 
If, instead, the minor allele (nMin1 or nMin2) had a copy number of 1, or the mutation was detected on only one or no chromosome (noChrsBearingMut ≤ 1),
the mutation was classified as possibly clone late (possiblylate), potentially indicating limited clonal expansion or structural constraints in copy number architecture. 
 
 Mutations with high CCF that did not meet the specific criteria for cloneEarly or possiblylate were assigned to the clone late category ((cloneLate), 
 representing clonal mutations with less distinctive chromosomal features. Any remaining mutations, which did not clearly fall into the above categories, were grouped as ‘clonalNA’.

## Getting started
Clone the repository git clone
git clone git@github.com-xtgithubhe:xtgithubhe/BBDPCtiming.git

## Input and output
Example input and output files are located in ./Data.

## Basic usage example (./bbDPCtiming.sh is designed for systems using SLURM.)
python3 ./bbDPCtiming.py ${sampleID} ${ssDPCIfolderpath} ${dpcoutfolderpath} ${outpath4in}

**Contact:** xiaotong.he@manchester.ac.uk


# IPMNPDAC_RNAseq
Scripts for integrative analysis of genomic, transcriptomic, and immune data in IPMN-PDAC.
Key features include:

- Dataset preparation and preprocessing

- Clustering and subtype assignment

- Batch correction using ComBat-Seq

- Cancer hallmark annotation

- Immune deconvolution with EPIC

**Contact:** LeonorPatricia.SchubertSantana@glasgow.ac.uk
  
