# Context
This project directory houses the scripts to infer cell composition of bulk RNAseq from mouse samples treated with a range of exposures to model orbital conditions for astronauts. It uses several publically available tools and follows the general framework established in [Gural et al. 2025](doi.org/10.1371/journal.pgen.1011807). Specifically, it aims to: 
1) Perform standard QC check and exploratory analysis of bulk tissue RNAseq from 120 mouse samples
2) Load, process, and generate marker sets from publicly available single-cell/nucleus RNAseq
3) Join the data from 1) and 2) in a cell type deconvolution approach to produce estimates of the relative proportion of each cell type contibutuing to the observed whole tissue counts
4) Integrate these estimates in downstream analysis to attribute sub-components of the overall expression changes between conditions and samples to specific underlying changes in cell makeup