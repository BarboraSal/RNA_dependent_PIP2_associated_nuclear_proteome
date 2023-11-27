#' The results were used to perform enrichment analysis shown in Table 2
#' Process modification files downloaded from the phosphosite plus database
#' Match the IDR data to the site table to identify if IDRs containing PTMs
#' Summarize and export the results

########################################################################
# Required packages

library(data.table)
library(tidyverse)
library(openxlsx)

########################################################################
# Setup path to working directory
setwd("c:/path/to/wd")

########################################################################
# function to process the PTM sites files downloaded from PSPdb
process_ptms <- function(ptm_file) {
  
  site_data <- fread(ptm_file)
  
  # clean
  site_data$MOD_RSD_site <- gsub("-.*", "", site_data$MOD_RSD)
  site_data$MOD_RSD_site <- gsub("\\D{1}", "", site_data$MOD_RSD_site)
  site_data$MOD_RSD_site <- as.numeric(site_data$MOD_RSD_site)
  
  return(site_data)
  
}

# function to match the PTM sites to IDRs
match_ptms <- function(idr_file, modification_list) {
  
  prox <- 0
  
  data <- idr_file %>%
    dplyr::rename("UniprotID" = uniprot) %>%
    dplyr::mutate(
      
      PROX_Ace_Prot_IN = UniprotID %in% modification_list$Acetylation$ACC_ID,
      PROX_Ace_Site_IN = rep(0, times = nrow(.)),
      
      PROX_Phos_Prot_IN = UniprotID %in% modification_list$Phosphorylation$ACC_ID,
      PROX_Phos_Site_IN = rep(0, times = nrow(.)),
      
      PROX_Ubi_Prot_IN = UniprotID %in% modification_list$Ubiquitination$ACC_ID,
      PROX_Ubi_Site_IN = rep(0, times = nrow(.)),
      
      PROX_Met_Prot_IN = UniprotID %in% modification_list$Methylation$ACC_ID,
      PROX_Met_Site_IN = rep(0, times = nrow(.)),
      
      PROX_Sumo_Prot_IN = UniprotID %in% modification_list$Sumoylation$ACC_ID,
      PROX_Sumo_Site_IN = rep(0, times = nrow(.))
    )
  
  # Acetylation
  for (i in 1:nrow(data)){
    
    if (data$PROX_Ace_Prot_IN[i] == TRUE){
      
      p <- data$UniprotID[i]
      start <- data$idr_max_start[i]
      end <- data$idr_min_end[i]
      
      min <- start-prox-1
      max <- end+prox+1
      
      mod <- modification_list$Acetylation %>% filter(ACC_ID == p)
      
      for (j in 1:nrow(mod)){
        mod$IN[j] <- mod$MOD_RSD_site[j] < max & mod$MOD_RSD_site[j] > min
      }
      
      result <- as.numeric(sum(mod$IN == TRUE))
      
      data$PROX_Ace_Site_IN[i] <- result
      
    }else{
      data$PROX_Ace_Site_IN[i] <- 0
    }
    
  }
  
  # Phosphorylation
  for (i in 1:nrow(data)){
    
    if (data$PROX_Phos_Prot_IN[i] == TRUE){
      
      p <- data$UniprotID[i]
      start <- data$idr_max_start[i]
      end <- data$idr_min_end[i]
      
      
      min <- start-prox-1
      max <- end+prox+1
      
      mod <- modification_list$Phosphorylation %>% filter(ACC_ID == p)
      
      for (j in 1:nrow(mod)){
        mod$IN[j] <- mod$MOD_RSD_site[j] < max & mod$MOD_RSD_site[j] > min
      }
      
      result <- as.numeric(sum(mod$IN == TRUE))
      
      data$PROX_Phos_Site_IN[i] <- result
      
    }else{
      data$PROX_Phos_Site_IN[i] <- 0
    }
    
  }
  
  # Ubiquitination
  for (i in 1:nrow(data)){
    
    if (data$PROX_Ubi_Prot_IN[i] == TRUE){
      
      p <- data$UniprotID[i]
      start <- data$idr_max_start[i]
      end <- data$idr_min_end[i]
      
      
      min <- start-prox-1
      max <- end+prox+1
      
      mod <- modification_list$Ubiquitination %>% filter(ACC_ID == p)
      
      for (j in 1:nrow(mod)){
        mod$IN[j] <- mod$MOD_RSD_site[j] < max & mod$MOD_RSD_site[j] > min
      }
      
      result <- as.numeric(sum(mod$IN == TRUE))
      
      data$PROX_Ubi_Site_IN[i] <- result
      
    }else{
      data$PROX_Ubi_Site_IN[i] <- 0
    }
    
  }
  
  # Methylation
  for (i in 1:nrow(data)){
    
    if (data$PROX_Met_Prot_IN[i] == TRUE){
      
      p <- data$UniprotID[i]
      start <- data$idr_max_start[i]
      end <- data$idr_min_end[i]
      
      
      min <- start-prox-1
      max <- end+prox+1
      
      mod <- modification_list$Methylation %>% filter(ACC_ID == p)
      
      for (j in 1:nrow(mod)){
        mod$IN[j] <- mod$MOD_RSD_site[j] < max & mod$MOD_RSD_site[j] > min
      }
      
      result <- as.numeric(sum(mod$IN == TRUE))
      
      data$PROX_Met_Site_IN[i] <- result
      
    }else{
      data$PROX_Met_Site_IN[i] <- 0
    }
    
  }
  
  # Sumoylation
  for (i in 1:nrow(data)){
    
    if (data$PROX_Sumo_Prot_IN[i] == TRUE){
      
      p <- data$UniprotID[i]
      p <- data$UniprotID[i]
      start <- data$idr_max_start[i]
      end <- data$idr_min_end[i]
      
      
      min <- start-prox-1
      max <- end+prox+1
      
      mod <- modification_list$Sumoylation %>% filter(ACC_ID == p)
      
      for (j in 1:nrow(mod)){
        mod$IN[j] <- mod$MOD_RSD_site[j] < max & mod$MOD_RSD_site[j] > min
      }
      
      result <- as.numeric(sum(mod$IN == TRUE))
      
      data$PROX_Sumo_Site_IN[i] <- result
      
    }else{
      data$PROX_Sumo_Site_IN[i] <- 0
    }
    
    
}
  
  return(data)
  
}

########################################################################
# Process the PTM site tables 
########################################################################

# list all site tables
filenames <- list.files(pattern = "dataset_mod.txt")
names <- gsub("(.*)_.*_.*_mod.txt", "\\1", filenames)

# load and process all site tables
modification_list <- lapply(filenames, process_ptms)
names(modification_list) <- names

rm(filenames, names)

########################################################################
# Match the PTMsites to IDRs
########################################################################

# load all input IDR files
filenames <- list.files(pattern = "_results.xlsx")
names <- gsub("_results.xlsx", "", filenames)

files <- lapply(filenames, openxlsx::read.xlsx, sheet = "KR_motif_table_filt2")
names(files) <- names

# match to PTMsites
files_processed <- lapply(files, match_ptms, modification_list = modification_list )
names(files_processed) <- names

########################################################################
# Export results
########################################################################
# working directory
setwd("c:/path/to/wd")

# save Rdata with the matched files
save(files_processed, file = "All_files_processed.Rdata")

for (i in 1:length(files_processed)){
  
  data <- files_processed[[i]]
  name <- names(files_processed)[i]

  legend <- "Each sheet contains the number of idrs proximal to a specific modification"
  
  summary_acetylation <- data %>% 
    group_by(PROX_Ace_Site_IN) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::rename("Number_of_proximal_sites" = PROX_Ace_Site_IN)
  
  summary_phosphorylation <- data %>% 
    group_by(PROX_Phos_Site_IN) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::rename("Number_of_proximal_sites" = PROX_Phos_Site_IN)
  
  summary_ubiquitination <- data %>% 
    group_by(PROX_Ubi_Site_IN) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::rename("Number_of_proximal_sites" = PROX_Ubi_Site_IN)
  
  summary_methylation <- data %>% 
    group_by(PROX_Met_Site_IN) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::rename("Number_of_proximal_sites" = PROX_Met_Site_IN)
  
  summary_sumoylation <- data %>% 
    group_by(PROX_Sumo_Site_IN) %>%
    dplyr::summarise(count = n()) %>%
    dplyr::rename("Number_of_proximal_sites" = PROX_Sumo_Site_IN)
  
  # set output filename
  output_file_name <- paste(name, "PTM_in_IDR_summary.xlsx", sep = "_" )
  
  # write the output excel file with all results
  results_sheets <- list("Legend" = legend,
                         "Acetylation" = summary_acetylation,
                         "Phosphorylation" = summary_phosphorylation, 
                         "Ubiquitination" = summary_ubiquitination,
                         "Methylation" = summary_methylation, 
                         "Sumoylation" = summary_sumoylation)
  
  openxlsx::write.xlsx(results_sheets, file = output_file_name)

}

