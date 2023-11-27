#' This code was used to generate results presented in Figure 2 and relevant supplementary files
#' Retrieve information about IDRs from the D2P2 database: https://d2p2.pro/ 
#' The IDRs are filtered based on minimum length of 20 and being predicted by at least 3 predictors
#` Calculate their pI using the pI() function from the R package "peptides"
#` Calculate their GRAVY scores using the hydrophobicity() function from the R package "peptides"

########################################################################
# Required packages
library("dplyr")
library("ggplot2")
library("openxlsx")
library("Peptides")
library("purrr")
library('rjson')
library("seqinr")

########################################################################
# Setup path to working directory
setwd('c:/path/to/wd')

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

########################################################################
# RETRIVE PREDICTIONS FROM D2P2
########################################################################

# function to process the response from D2P2
process_response <- function(response) {
  prediction_retrieved <- unlist(response[[1]][[1]][[3]]$disorder$disranges)
  prediction_retrieved <- matrix(prediction_retrieved, ncol = 3, byrow = TRUE)
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
  
}

########################################################################
# Analyze the D2P2 results
########################################################################

# For each protein retrieved from d2p2 count IDR that passed the selected cutoff and  Count the number of predictors
# note, this analysis does not address individual predictors and the overlap between them

# how many IDRs per protein passed minimal length
match$idr_passed_min_length <- sapply(predictions, function(prediction) {
  if (is.data.frame(prediction)) {
    sum(prediction$passed_min_idr_length == "yes")
  } else {
    NA
  }
})

# how many predictors correspond to these IDRs
match$predictor_number <- sapply(predictions, function(prediction) {
  if (is.data.frame(prediction)) {
    prediction <- prediction %>%
      filter(passed_min_idr_length == "yes")
    length(unique(prediction$predictor))
  } else {
    NA
  }
})


match$predictor_passed_cutoff <- match$predictor_number >= min_predictor
match$predictor_passed_cutoff <- ifelse(match$predictor_passed_cutoff, "yes", "no")

# Count the number of proteins and calculate the percentage retrieved from D2P2
count_submitted <- nrow(match)
count_d2p2_retrieved <- sum(match$d2p2_matched == "yes", na.rm = TRUE) # Ensure NA values are not included
percentage_retrieved <- (count_d2p2_retrieved / count_submitted) * 100

count_idr_passed_cutoff <- sum(match$predictor_passed_cutoff == "yes", na.rm = TRUE)
percentage_idr_passed_cutoff <- (count_idr_passed_cutoff / count_d2p2_retrieved) * 100

# Create a summary dataframe with named columns
summary_percentage <- data.frame(
  count_submitted = count_submitted,
  count_d2p2_retrieved = count_d2p2_retrieved,
  percentage_retrieved = percentage_retrieved,
  count_idr_passed_cutoff = count_idr_passed_cutoff,
  percentage_idr_passed_cutoff = percentage_idr_passed_cutoff, 
  min_idr_length = min_idr_length,
  min_predictor = min_predictor
)

# Create a summary data frame
summary_percentage <- t(summary_percentage)
summary_percentage <- as.data.frame(summary_percentage)
names(summary_percentage) <- c("result")
summary_percentage$value <- row.names(summary_percentage)

########################################################################
# Analyze the IDR sequences, their pI, hydrophobicity, and length
########################################################################

# retrieve the IDR sequences
# calculate their pI using the pI() function from the R package "peptides"
# calculate their GRAVY scores using the hydrophobicity() function from the R package "peptides"

# Select unique protein IDs matched to D2P2
unique_ids <- match$uniprotID_submitted[match$d2p2_matched == "yes"]

# Filter FASTA sequences of the unique proteins from UniProt matched to D2P2
fasta_predicted <- fasta[names(fasta) %in% unique_ids]

# Define a function to process IDR sequences
process_idr <- function(prediction, fasta_sequence) {
  prediction %>%
    dplyr::mutate(
      start_pos = as.numeric(start_pos),
      end_pos = as.numeric(end_pos),
      idr_sequence = map2_chr(start_pos, end_pos, ~paste(fasta_sequence[.x:.y], collapse = "")),
      pI_sequence = map_dbl(idr_sequence, ~Peptides::pI(.x)),
      hydrophobicity_sequence = map_dbl(idr_sequence, ~Peptides::hydrophobicity(.x))
    )
}

# Process predictions and combine them into a single data frame
all_results_table <- map_df(unique_ids, ~{
  protein_id <- .x
  fasta_sequence <- fasta_predicted[[protein_id]]
  prediction <- predictions[[protein_id]]
  process_idr(prediction, fasta_sequence)
})

########################################################################
# Summarize the results
########################################################################

# select IDR based on the minimum length
filtered_results <- all_results_table %>% 
  filter(passed_min_idr_length == "yes")

# calculate summary statistics per predictor
summary <- filtered_results %>% 
  group_by(predictor) %>% 
  summarise(
    count_IDR = n(),
    count_protein = n_distinct(protein_id),
    median_pI = median(pI_sequence, na.rm = TRUE),
    median_gravy = median(hydrophobicity_sequence, na.rm = TRUE),
    median_length = median(end_pos - start_pos + 1, na.rm = TRUE),
    mean_pI = mean(pI_sequence, na.rm = TRUE),
    mean_gravy = mean(hydrophobicity_sequence, na.rm = TRUE),
    mean_length = mean(end_pos - start_pos + 1, na.rm = TRUE)
  ) %>%
  ungroup()

########################################################################
# Export the results
########################################################################

# write the output excel file with all results
results_sheets <- list("D2P2_match" = match,
                         "idr_percentage" = summary_percentage,
                         "predictor_summary" = summary)

openxlsx::write.xlsx(results_sheets, file = output_file_name)


# save the R object for the following analysis
results <- list(match, summary_percentage, all_results_table, filtered_results, summary)
names(results) <- c("match.table", "percentage.summary", "idr.all", "idr.filtered", "predictor.summary")
save(results, file = paste(output_file_name, ".RData", sep = ""))


  