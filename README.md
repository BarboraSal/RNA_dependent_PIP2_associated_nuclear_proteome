# RNA_dependent_PIP2_associated_nuclear_proteome
The code and analysis files associated with a manuscript "The RNA-dependent interaction of phosphatidylinositol 4,5-bisphosphate with intrinsically disordered proteins contributes to nuclear compartmentalization" by Sztacho M., Cervenka, J., Salovska, B., et al. 

## R version
The code was built on R version 4.3.1. 

## Input datasets
An example FASTA sequence database, ScanProsite motif table, and an example ouput from Analysis 1 are provided in the Example_data folder. 

## R files description 
### Analysis 1

The first analysis takes the input FASTA sequence database to retrieve intrinsincally disordered regions (IDRs) in proteins from the D2P2 database (https://d2p2.pro). 

The IDRs are then matched to three K/R rich sequence motifs identified by the ScanProsite tool (https://prosite.expasy.org/scanprosite). 

The K/R motif-matched IDRs are then analyzed in terms of their physicochemical properties (pI, GRAVY score) and length.

The Modification site datasets were downloaded from the PhosphoSite Plus database (https://www.phosphosite.org/homeAction) and matched to the IDRs. 

### Analysis 2

The second analysis takes the input FASTA sequence database to retrieve intrinsincally disordered regions (IDRs) in proteins from the D2P2 database (https://d2p2.pro). 

The IDRs are then analyzed in terms of their physicochemical properties (pI, GRAVY score) and length. 

## Further information

The analysis is described in the Methods section of the manuscript. 
