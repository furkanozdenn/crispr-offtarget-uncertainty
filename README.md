# Learning to quantify uncertainty in off-target activity for CRISPR guide RNAs

> crispAI is a deep learning-based tool that predicts the off-target cleavage activity for a given single guide RNA (sgRNA) and target DNA sequence pair and quantifies the uncertainty in off-target cleavage activity.
> crispAI-aggregate is an uncertainty aware genome-wide specificity score for sgRNA's.

> <a href="https://en.wikipedia.org/wiki/Deep_learning" target="_blank">**Deep Learning**</a>,  <a href="https://www.nature.com/articles/s41467-018-04252-2" target="_blank">**CRISPR-based genome editing**</a>, <a href="https://link.springer.com/content/pdf/10.1007/978-3-319-54339-0.pdf" target="_blank">**Uncertainty Quantification**</a>


---

## Authors

<a href="https://www.cs.ox.ac.uk/people/furkan.ozden/" target="_blank">**Furkan Ozden**</a> and <a href="https://www.cs.ox.ac.uk/people/peter.minary/" target="_blank">**Peter Minary**</a>

---

## Questions & comments 

[lastauthorname].[lastauthorsurname]@cs.ox.ac.uk

## Table of Contents 

> Warning: Please note that crispAI software is free for academic usage, while it is licenced for commercial usage. Please get in touch for more info.

- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Citations](#citations)
- [License](#license)

---


## Installation

### Step-1: Download and untar USCS chroms required for Cas-OFFinder (if you do not plan to use crispAI-aggregate score you can skip this step).

Download target organism's chromosome FASTA files.

- http://hgdownload.soe.ucsc.edu/downloads.html (UCSC genome sequences library)

Extract all FASTA files in a directory and name it 'ucsc_chroms'.

For example (human chromosomes, in POSIX environment):

    $ wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
    $ mkdir -p /var/chromosome/human_hg19
    $ tar zxf chromFa.tar.gz -C /var/chromosome/human_hg19
    $ ls -al /var/chromosome/human_hg19
      drwxrwxr-x.  2 user group      4096 2013-10-18 11:49 .
      drwxrwxr-x. 16 user group      4096 2013-11-12 12:44 ..
      -rw-rw-r--.  1 user group 254235640 2009-03-21 00:58 chr1.fa
      -rw-rw-r--.  1 user group 138245449 2009-03-21 01:00 chr10.fa
      -rw-rw-r--.  1 user group 137706654 2009-03-21 01:00 chr11.fa
      -rw-rw-r--.  1 user group 136528940 2009-03-21 01:01 chr12.fa
      -rw-rw-r--.  1 user group 117473283 2009-03-21 01:01 chr13.fa
      -rw-rw-r--.  1 user group 109496538 2009-03-21 01:01 chr14.fa

Move obtained directory 'ucsc_chroms' to project folder ./crispAI_score

    $ mv ucsc_chroms ./crispAI_score

### Step-2: Use conda to create environment with required packages.
    
    $ conda env create -f env/crispAI_env.yml
    $ conda activate crispAI_env

### Step-3: Install required R packages for NuPoP library using R_environment.csv (optionally run restore_environment.R).
- Note: crispAI uses R version 4.2 to annotate target sites. You should have R version 4.2 installed on your system.

    $ cd env/
    $ Rscript restore_environment.R


### Step-4: Test installation by running on example input data, this will run with default parameters.

    $ python crispAI.py --input_file example_offt_input.txt --mode offt-score


---

## Usage Example 

### offt-score mode: Off-target cleavage activity prediction for sgRNA-target pairs.

- Input: See example_offt_input.txt file for example input.
  
