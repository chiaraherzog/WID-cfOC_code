# Preparation of patterns and summary files for cfDNAme assay
# Author: chiara herzog, chiara.herzog@uibk.ac.at
# Date: 20 Dec 2021
# Outputs:
#    1 detail file with pattern per sample and methylation status. Removing any duplicate calls (not so many)
#    1 summary file with region and methylated reads

count_patterns <- function(input, output, targets){

  dir <- input
  dir1 <- output
  
  suppressPackageStartupMessages(require(dplyr))
  suppressPackageStartupMessages(require(tidyverse))
  suppressPackageStartupMessages(require(stringr))
  
  if(!dir.exists(paste0(dir1, "/processed/"))){
    dir.create(paste0(dir1, "/processed/"))
  }
  
  if(!dir.exists(paste0(dir1, "/summary/"))){
    dir.create(paste0(dir1, "/summary/"))
  }
  
  
  # Define regions and CpGs of interest
  regions <- read.table(paste0(targets, "/efc_5_target_genomic_seq_CpG_positions.fna"))
  EFC_144_genome <- as.numeric(regions$V1[25:35])
  EFC_204_genome <- as.numeric(regions$V1[38:56])
  EFC_228_genome <- as.numeric(regions$V1[59:67])
  
  #----------------------------
  # Read files 
  setwd(dir)
  files <- list.files(dir)
  files <- files[grepl("CpG", files)]
  
  pB <- txtProgressBar(min=1,max=length(files), width =50L, style = 3)
  
  for(i in 1:length(files)){
    
    setTxtProgressBar(pB, i)
      # Define the name of the file
    name <- ifelse(grepl("IP12_5MP|NA_5MP|NTC_5MP", files[i]),
                   paste(stringr::str_split(files[i], "_", simplify = TRUE)[,4],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,5],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,6],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,11],
                         sep = "_"),
                   paste(stringr::str_split(files[i], "_", simplify = TRUE)[,3],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,4],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,5],
                         stringr::str_split(files[i], "_", simplify = TRUE)[,6],
                         sep = "_"))
      
      file <- read.table(files[i], skip = 1)
      colnames(file) <- c("seqID", "methylationstate", "region", "position", "call")

      # Removing any CpGs not of interest
      file <- file %>%
        dplyr::filter(
                 (region == "EFC_144_genome" & position %in% EFC_144_genome) | 
                 (region == "EFC_204_genome" & position %in% EFC_204_genome) | 
                 (region == "EFC_228_genome" & position %in% EFC_228_genome)) %>%
        dplyr::distinct()
      
      # Methylated CpGs per read
      tmp1 <- suppressMessages(file %>%
        dplyr::group_by(seqID, region) %>%
        dplyr::summarise(methylated = sum(call=="Z")/n()))
      
      #----------------------------
      # Patterns resolved by region
      
      # find which regions are available
      
      regions <- unique(file$region)
      
      # EFC_144
      if("EFC_144_genome" %in% regions){
        tmp <- file |> 
          dplyr::filter(region == "EFC_144_genome" & position %in% EFC_144_genome) |> 
          dplyr::ungroup() |> 
          tidyr::pivot_wider(id_cols = c("seqID", "region"),
                      names_from = position,
                      values_from = call) |> 
          as.data.frame() 
        
        tmp <- tmp %>%
          tidyr::unite(pattern, colnames(tmp)[3:ncol(tmp)], sep = "")
        tmp$pattern <- gsub("z", "0", tmp$pattern)
        tmp$pattern <- gsub("Z", "1", tmp$pattern)
        tmp$pattern <- gsub("NULL", "x", tmp$pattern)
        
        efc144 <- tmp[nchar(tmp$pattern)==length(EFC_144_genome),]
      } else {
        efc144 <- data.frame(seqID = NA,
                             region = "EFC_144_genome",
                             pattern = NA)
      }
      
      # EFC_204
      if("EFC_204_genome" %in% regions){
        tmp <- file %>%
          dplyr::filter(region == "EFC_204_genome" & position %in% EFC_204_genome) %>%
          dplyr::ungroup() |> 
          tidyr::pivot_wider(id_cols = c("seqID", "region"),
                      names_from = "position",
                      values_from = "call") %>%
          as.data.frame() 
        
        tmp <- tmp %>%
          tidyr::unite(pattern, colnames(tmp)[3:ncol(tmp)], sep = "")
        tmp$pattern <- gsub("z", "0", tmp$pattern)
        tmp$pattern <- gsub("Z", "1", tmp$pattern)
        tmp$pattern <- gsub("NULL", "x", tmp$pattern)
        
        efc204 <- tmp[nchar(tmp$pattern)==length(EFC_204_genome),]
      } else {
        efc204 <- data.frame(seqID = NA,
                             region = "EFC_204_genome",
                             pattern = NA)
      }
      
      # EFC_228
      if("EFC_228_genome" %in% regions){
        tmp <- file %>%
          dplyr::filter(region == "EFC_228_genome" & position %in% EFC_228_genome) %>%
          dplyr::ungroup() |> 
          tidyr::pivot_wider(id_cols = c("seqID", "region"),
                      names_from = "position",
                      values_from = "call") %>%
          as.data.frame() 
        
        tmp <- tmp %>%
          tidyr::unite(pattern, colnames(tmp)[3:ncol(tmp)], sep = "")
        tmp$pattern <- gsub("z", "0", tmp$pattern)
        tmp$pattern <- gsub("Z", "1", tmp$pattern)
        tmp$pattern <- gsub("NULL", "x", tmp$pattern)
        
        efc228 <- tmp[nchar(tmp$pattern)==length(EFC_228_genome),]
      } else {
        efc228 <- data.frame(seqID = NA,
                             region = "EFC_228_genome",
                             pattern = NA)
      }
      
      
      #----------------------------
      patterns <- rbind(efc144, efc204, efc228)
      patterns <- suppressMessages(patterns %>%
                                     dplyr::filter(!is.na(seqID)) %>%
                                     dplyr::inner_join(tmp1, by = c("seqID", "region")))
      
      patterns <- patterns %>%
        dplyr::mutate(onemissing = ifelse(str_count(pattern, pattern = "1") == nchar(pattern)-1, "TRUE",
                                   "FALSE")) %>%
        dplyr::mutate(sample = name)
      
      summary <- patterns %>%
        dplyr::group_by(region) %>%
        dplyr::summarise(sample = name,
                  reads = n(),
                  reads_onemissingorfull = sum((methylated==1) | (onemissing==TRUE)),
                  onemissingorfull = sum(methylated==1 | onemissing==TRUE)/n()*100,
                  reads_methylated = sum(methylated == 1),
                  fully_methylated = sum(methylated == 1)/n()*100,
                  reads_onemissing = sum(onemissing == TRUE),
                  onemissing = sum(onemissing==TRUE)/n()*100) %>%
        unique()
      
      #----------------------------
      # Save by-read output
      write.table(patterns, file = paste0(dir1, "/processed/", name, "_patterns.csv"),
                sep = ",",
                quote = FALSE,
                row.names = FALSE,
                col.names = TRUE)
      
      write.table(summary, file = paste0(dir1, "/summary/", name, "_summary.csv"),
                  sep = ",",
                  quote = FALSE,
                  row.names = FALSE,
                  col.names = TRUE)
      
      rm(file, tmp, summary, patterns, efc144, efc228, efc204, name);gc()
      
  }
  close(pB)
  
  cat("\ndone")
}