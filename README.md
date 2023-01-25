# uASPIre_UTR_CDS

This repository contains the data of _**u**ltradeep **A**cquisition of **S**equence-**P**henotype **I**nter**re**lations_ (uASPIre) and the code for processing raw NGS data as presented in the manuscript "Ultradeep characterisation of translational sequence determinants refutes rare-codon hypothesis and unveils quadruplet base pairing of initiator tRNA and transcript".


## Organization of the repository

1. [**data**](data) contains the data required to reproduce all figures in the main text.


## Hardware requirements
This code package requires a Unix system with >10 GB RAM and optimally multiple cores.


## Software requirements
The code runs on Unix systems only, and was developed and tested on R version 4.2.1 running on a Red Hat Enterprise Linux Server (release 7.9).


## Installation
No installation is needed, the code can be run as is.


## Software dependencies
This code requires the following UNIX tools and packages:

+ AGREP 3.41.5 (e.g. available via https://github.com/Wikinaut/agrep)
+ TRE agrep 0.8.0 (e.g. available via https://wiki.ubuntuusers.de/tre-agrep/)


## Functionality
"The algorithms for processing of NGS data for this project were written in bash and python. Briefly, forward and reverse reads retrieved from fastq files were paired and all reads with more than six consecutive unidentified nucleotides were removed. Afterwards, target fragments were selected by a 10-bp constant region (GAGCTCGCAT, max. 3 mismatches) and sequences from different samples were deconvoluted by their unique combination of two 6-bp indices (Supplementary Tab. 3). Next, the discriminator state was determined by searching for the presence of an attP or attR site corresponding to the sequences GGGTTTGTACCGTACAC or GCCCGGATGATCCTGAC, respectively (max. 3 mismatches). RBS sequences were determined by retrieving the 17 nucleotides upstream of the bxb1 start codon." (from Höllerer _et al_., 2020)


## License
Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)

(c) Simon Höllerer, 2023, ETH Zurich

https://creativecommons.org/licenses/by-nc/4.0/legalcode

This license is valid for the data and the code shared in this repository.


## Date
January 2023

## Contact
simon.hoellerer@bsse.ethz.ch

Simon Höllerer, ETH Zürich, D-BSSE, BPL, Basel, Switzerland 
