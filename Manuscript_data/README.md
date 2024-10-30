# Manuscript data
This folder contains data that were used to generated all figures in the the manuscript "The RNA-dependent association of phosphatidylinositol 4,5-bisphosphate with intrinsically disordered proteins contribute to nuclear compartmentalization" by Sztacho M., Cervenka, J., Salovska, B., et al.

## FASTA_db
Contains all FASTA sequence databases that were used for the bioinformatic analyses presented in the manuscript. 

## IDR_KR_motif_results

Contains all data that were used to generate the following main figures: 2G, 2H, 2I, and Table 2.
Contains all data that were used to generate the following supplementary figures: S11, S13, and S15. 

Contains all results from the analysis, which uses information about the presence of K/R motif in the dataset from the ScanProsite tool, and searches UniProtKB accession numbers of such proteins against D2P2. K/R motif were counted as present in IDR and reported only if at least three predictors from the D2P2 predicted IDR with a minimum length of 20 amino acid residues at the site of the K/R motif. The pI and hydrophobicity of the disordered regions were calculated using the R package “peptides”, functions pI, and hydrophobicity. The hydrophobicity of each IDR was determined based on grand average of hydropathy (GRAVY) value, calculated as the sum of hydropathy values of all amino acids in the IDR divided by the length of the IDR.

For the PTM sites mapping, a database of known posttranslational modification (PTM) sites was downloaded (August 6, 2022) from the PhosphoSitePlus database. For the following known PTMs: acetylation, methylation, phosphorylation, sumoylation, and ubiquitination, we explored whether they are located in the IDR containing the K/R motifs identified by the abovementioned analysis of K/R motifs enrichment in IDRs.

## IDR_prediction_results

Contains all data that were used to generate the following main figures: 2D, 2E, and 2F. 
Contains all data that were used to generate the following supplementary figures: S6, S7, S8, S9, S10, and S14. 

Contains all results from IDR prediction using the Database of Disordered Protein Predictions (D2P2).The pI and hydrophobicity of the disordered regions were calculated using the R package “peptides”, functions pI, and hydrophobicity. The hydrophobicity of each IDR was determined based on grand average of hydropathy (GRAVY) value, calculated as the sum of hydropathy values of all amino acids in the IDR divided by the length of the IDR. 

## Main_figures
Contains all data that were used to generate the following main figures: 1C, 4A, 4B, 4C, 4D, 5B, and 5D. 

## Supplementary_figures
Contains all data that were used to generate the following supplementary figures: S1, S2, S18, S19, S20, S21, S22B, S22C, and S23B. 
