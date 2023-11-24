#' This code was used to generate results presented in Figure 2 and relevant supplementary files
#' First, we retrieved information about IDRs from the D2P2 database: https://d2p2.pro/ 
#' The the IDRs are filtered based on minimum length of 20 and being predicted by at least 3 predictors
#' Second, we match those to motifs identified and imported from the ScanProsite tool: https://prosite.expasy.org/scanprosite/ 
#' And export the results as an excel table for further analysis

########################################################################
# Required packages
library("dplyr")
library("openxlsx")
library("Peptides")
library("purrr")
library('rjson')
library("seqinr")

########################################################################
# Setup path to working directory
setwd('c:/Users/user/Documents/Barbora 2021/Projects_OLD/13_MSz_JCer_Script/_git/analysis_01')
list.files()

########################################################################
# DEFINE ANALYSIS PARAMETERS
########################################################################

# define IDR length that is going to be used for minimum length filtering as length >= idr_length
min_idr_length <- 20

# define the minimum number of predictors for high confident results reporting
min_predictor <- 3

########################################################################
# DEFINE INPUT FILES
########################################################################

# FASTA sequence database containing proteins of interest
input_file_name <-  "example_fasta.fasta"

# define naming for the output files and visualization
name <- "example"
output_file_name <- paste(name, "results.xlsx", sep = "_")

# read the FASTA sequence database from UniProt
fasta <- read.fasta(file = input_file_name, 
                     as.string = FALSE,  
                    forceDNAtolower = TRUE)

# reformat the FASTA headers
names(fasta) <- gsub("sp\\|(\\w.*)\\|.*", "\\1", names(fasta) )

# read motif table, these are tab-delimited predictions from the ScanProsite tool
motif_table <- read.delim(file = "example_motif.txt", sep = "\t", header = T,  stringsAsFactors = F )

# create UniProt accession column from the FASTA header
motif_table$uniprot <- gsub("sp\\|(\\w.*)\\|.*", "\\1", motif_table$id )

########################################################################
# RETRIVE PREDICTIONS FROM D2P2
########################################################################

# function to process the response from D2P2
process_response <- function(response) {
  prediction_retrieved <- unlist(response[[1]][[1]][[3]]$disorder$disranges)
  prediction_retrieved <- matrix(res, ncol = 3, byrow = TRUE)
  prediction_retrieved <- data.frame(prediction_retrieved)
  
  # Define column names and ordering to match the original code's output
  colnames(prediction_retrieved) <- c("predictor", "start_pos", "end_pos")
  
  # Calculate the length and other derived columns
  prediction_retrieved$length <- (as.numeric(prediction_retrieved$end_pos) - as.numeric(prediction_retrieved$start_pos)) + 1
  prediction_retrieved$passed_min_idr_length <- ifelse(prediction_retrieved$length >= min_idr_length, "yes", "no")
  prediction_retrieved$seq_id <- paste(prediction_retrieved$start_pos, prediction_retrieved$end_pos, sep = "_")
  prediction_retrieved$protein_id <- rep(NA, nrow(prediction_retrieved))
  
  # Adjust the order of columns to match the original structure
  prediction_retrieved <- prediction_retrieved[, c("predictor", "start_pos", "end_pos", "length", "passed_min_idr_length", "seq_id", "protein_id")]
  
  return(prediction_retrieved)
}

# make protein ids table that will be used to submit to D2P2
ids <- data.frame(id = names(fasta), 
                  call = sapply(names(fasta), 
                                FUN = function(x) paste('http://d2p2.pro/api/seqid/["',  x,  '"]', sep = "")
                  )
)

# pre-allocate space for 'predictions' and 'match'
predictions <- setNames(vector("list", length(ids$id)), ids$id)
match <- data.frame(uniprotID_submitted = ids$id, d2p2_matched = rep(NA, length(ids$id)))

# retrieve the predictions from D2P2 and process the response
for (i in seq_along(ids$id)) {
  protein_id <- ids$id[i]
  response <- fromJSON(readLines(ids$call[i], warn = FALSE))
  
  if (length(response[[1]]) > 0) {
    res <- process_response(response)
    res$protein_id <- protein_id
    predictions[[i]] <- res
    match$d2p2_matched[i] <- "yes"
  } else {
    predictions[[i]] <- "Not_matched"
    match$d2p2_matched[i] <- "no"
  }
  
  rm(res, response)
}

########################################################################
# Analyze the D2P2 results and motif matching
########################################################################

# count the number of proteins and calculate the percentage retrieved from D2P2
count_submitted <- nrow(match)
count_d2p2_retrieved <- sum(match$d2p2_matched == "yes")
percentage_retrieved <- (count_d2p2_retrieved / count_submitted) * 100

# indicate whether each protein in motif_table matched in D2P2
matched_proteins <- match$uniprotID_submitted[match$d2p2_matched == "yes"]
motif_table$uniprot_in_d2p2_predictions <- ifelse(motif_table$uniprot %in% matched_proteins, "yes", "no")

# indicate in the match table whether the protein id is found in scanprosite results
match$scanprosite_motif_matched <- ifelse(match$uniprotID_submitted %in% motif_table$uniprot, "yes", "no")

########################################################################
# save results to excel sheets
sheet1 <- match
sheet2 <- motif_table

########################################################################
# report the number and percentage of protein IDs with KR motifs
count_scanprosite_match <- sum(match$scanprosite_motif_matched == "yes")
percentage_scanprosite_match <- (count_scanprosite_match / count_submitted) * 100

# filter unique protein IDs that are both in the motif table and D2P2 prediction
uq <- intersect(unique(motif_table$uniprot), matched_proteins)

# report the number and percentage of unique protein IDs both in the motif table and D2P2 prediction
count_d2p2_and_scanprosite_match <- length(uq)
percentage_d2p2_and_scanprosite_match <- count_d2p2_and_scanprosite_match/nrow(match)*100

# filter the original motif table to only include rows where the UniProt ID has a corresponding prediction in D2P2
motif_table_filtered <- motif_table %>%
  
  filter(uniprot_in_d2p2_predictions == "yes") %>%
  
  # initialize new columns with default NA values
  mutate(
    number_of_predictors = NA_integer_,
    idr_max_start = NA_integer_,
    idr_min_end = NA_integer_,
    idr_overlap_length = NA_integer_,
    idr_sequence = NA_character_,
    pI_sequence = NA_real_,
    hydrophobicity_sequence = NA_real_
  )


# convert prediction start and end positions to numeric
predictions <- lapply(predictions, function(pred) {
  if (is.list(pred) || is.data.frame(pred)) {
    pred$start_pos <- as.numeric(pred$start_pos)
    pred$end_pos <- as.numeric(pred$end_pos)
  }
  return(pred)
})

########################################################################
# iterate through the filtered motif table to process IDR and motif matches
for (i in 1:nrow(motif_table_filtered)) {
  
  # extract the motif details and the corresponding protein sequence
  motif_selection <- motif_table_filtered[i, ]
  uniprot_id <- motif_selection$uniprot
  sequence_id <- fasta[uniprot_id][[1]]
  
  # fetch the prediction data for the current UniProt ID
  prediction_selection <- predictions[[uniprot_id]]
  
  # check overlap of motif region with each predicted region using a vectorized operation
  overlap <- (motif_selection$start_pos >= prediction_selection$start_pos) & 
    (motif_selection$end_pos <= prediction_selection$end_pos)
  
  # if overlaps exist
  if (any(overlap)) {
    
    # filter matched predictions
    prediction_selection_match <- prediction_selection[overlap, ]
    
    # extract relevant information from the matched predictions and store in the filtered motif table
    motif_table_filtered$number_of_predictors[i] <- length(unique(prediction_selection_match$predictor))
    motif_table_filtered$idr_max_start[i] <- max(prediction_selection_match$start_pos, na.rm = TRUE)
    motif_table_filtered$idr_min_end[i] <- min(prediction_selection_match$end_pos, na.rm = TRUE)
    motif_table_filtered$idr_overlap_length[i] <- motif_table_filtered$idr_min_end[i] - motif_table_filtered$idr_max_start[i] + 1
    
    # extract sequence information of the overlapping IDR and calculate its properties
    idr_sequence <- sequence_id[motif_table_filtered$idr_max_start[i]:motif_table_filtered$idr_min_end[i]]
    idr_sequence <- paste(idr_sequence, collapse = "")
    motif_table_filtered$idr_sequence[i] <- idr_sequence
    motif_table_filtered$pI_sequence[i] <- Peptides::pI(idr_sequence)
    motif_table_filtered$hydrophobicity_sequence[i] <- Peptides::hydrophobicity(idr_sequence)
    
  } else {
    # if no matches were found, set relevant columns to default values
    motif_table_filtered$number_of_predictors[i] <- 0
    motif_table_filtered$idr_max_start[i] <- NA
    motif_table_filtered$idr_min_end[i] <- NA
    motif_table_filtered$idr_overlap_length[i] <- NA
    motif_table_filtered$idr_sequence[i] <- NA 
    motif_table_filtered$pI_sequence[i] <- NA
    motif_table_filtered$hydrophobicity_sequence[i] <- NA  
  }
  
}

########################################################################
# save results to excel sheets
sheet3 <- motif_table_filtered

# filter the output table based on the minimum length and minimum predictors parameters (as defined above)
motif_table_filtered_2 <- motif_table_filtered %>% 
  filter(number_of_predictors >= min_predictor, idr_overlap_length >= min_idr_length)

# save filtered results to excel sheets
sheet4 <- motif_table_filtered_2

########################################################################
# create analysis summary values
values <- c(
  'Number of protein ids in the id file used for D2P2 predictions (or total)' = count_submitted,
  'Number of protein ids successfully retrieved from D2P2 with some predictions' = count_d2p2_retrieved,
  'Percentage of successfully retrieved protein ids from total (related to above)' = percentage_retrieved,
  'Number of unique protein ids in the Scanprosite output file' = count_scanprosite_match,
  'Percentage of protein ids found in the Scanprosite table from total' = percentage_scanprosite_match,
  'Number of protein ids with both D2P2 prediction and Scanprosite motif' = count_d2p2_and_scanprosite_match,
  'Percentage of protein ids with both D2P2 prediction and Scanprosite motif' = percentage_d2p2_and_scanprosite_match,
  'Number of motifs in the Scanprosite output table' = nrow(motif_table),
  'Number of motifs used for the analysis (some motifs were discarded due to missing information from D2P2)' = nrow(motif_table_filtered),
  'How many of these motifs were mapped to predicted disordered region (at least 3 predictors)' = nrow(motif_table_filtered_2),
  'Percentage of motifs mapped to disordered region with high confidence from total motifs analyzed' = (nrow(motif_table_filtered_2) / nrow(motif_table_filtered)) * 100,
  'Average pI of the IDRs matched to KR motifs' = mean(motif_table_filtered_2$pI_sequence, na.rm = TRUE),
  'Median pI of the IDRs matched to KR motifs' = median(motif_table_filtered_2$pI_sequence, na.rm = TRUE),
  'Average GRAVY of the IDRs matched to KR motifs' = mean(motif_table_filtered_2$hydrophobicity_sequence, na.rm = TRUE),
  'Median GRAVY of the IDRs matched to KR motifs' = median(motif_table_filtered_2$hydrophobicity_sequence, na.rm = TRUE)
)

# create summary dataframe
summary <- tibble(
  parameter = names(values),
  value = unname(values),
  legend = names(values)
)

# save results to excel sheets
sheet5 <- summary

########################################################################
# write the output excel file with all results
excel_sheets <- list("protein_match" = sheet1,
                         "KR_motif_table_all" = sheet2,
                         "KR_motif_table_filt1" = sheet3,
                         "KR_motif_table_filt2" = sheet4,
                         "analysis_summary" = sheet5)

openxlsx::write.xlsx(excel_sheets, file = output_file_name)
