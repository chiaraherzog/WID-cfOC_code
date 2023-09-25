grab_stats <- function(dir, lib = NULL){
  
  if(!exists("lib")){
    stop("select library ukfoccs, ucl, or circadian")
  }
  
  report <- list.files(dir, full.names = T)
  report <- report[!grepl("splitting", report)]
  stats <- data.frame(matrix(nrow = length(report),
                             ncol = 6))
  colnames(stats) <- c("name", "reads", "mapping", "totalr", "meth_cpgs", "bisconv")
  
  for (i in 1:length(report)){
    
    cat(paste0(basename(report[i]), "..."))
    if(lib == "ukfoccs"){
    # fix the name 
      name <- case_when(grepl("NA", report[i]) ~ paste(substr(basename(report[i]), 9, 18), stringr::str_split(basename(report[i]), "_", simplify = TRUE)[,9], sep = "_"),
                        grepl("NTC", report[i]) ~ paste(substr(basename(report[i]), 9, 19), stringr::str_split(basename(report[i]), "_", simplify = TRUE)[,9], sep = "_"),
                        grepl("IP12", report[i]) ~ paste(substr(basename(report[i]), 9, 20), stringr::str_split(basename(report[i]), "_", simplify = TRUE)[,9], sep = "_"),
                        grepl("UKFOCS", report[i]) ~ str_split(basename(report[i]), "_", simplify = TRUE)[,3])
    } 
    
    if(lib == "ucl"){
      name <- case_when(grepl("NA", report[i]) ~ paste(substr(basename(report[i]), 1, 13)),
                        grepl("ST_", report[i]) ~ paste(stringr::str_split(basename(report[i]), "_", simplify = TRUE)[,1], 
                                                        stringr::str_split(basename(report[i]), "_", simplify = TRUE)[,2],
                                                        sep = "_"),
                        grepl("UKFOCS", report[i]) ~ paste(substr(basename(report[i]), 1, 16)),
                        grepl("S_", report[i]) ~ gsub("_", "-", paste(str_split(basename(report[i]), "_", simplify = TRUE)[,1],
                                                                      str_split(basename(report[i]), "_", simplify = TRUE)[,2],
                                                                      str_split(basename(report[i]), "_", simplify = TRUE)[,3],
                                                                      sep = "_")))
    }
    
    if(lib == "circadian"){
      name <- case_when(grepl("NA", report[i]) ~ substr(basename(report[i]), 9, 17),
                        grepl("NTC", report[i]) ~ substr(basename(report[i]), 9, 11),
                        grepl("IP12", report[i]) ~ substr(basename(report[i]), 9, 12),
                        grepl("_ST_", report[i]) ~ substr(basename(report[i]), 9, 21))  
    }
    
    # read in mapping rate
    alignment1 <- suppressMessages(read.table(report[i], skip = 6,
                                              sep = "\t",
                                              nrow = 6))
    meth <- suppressMessages(read.table(report[i], skip = 35,
                                        sep = "\t",
                                        nrow = 3))
    
    stats$name[i] <- name
    stats$reads[i] <- as.numeric(alignment1$V2[2])*2
    stats$totalr[i] <- as.numeric(alignment1$V2[1])*2
    stats$mapping[i] <- as.numeric(substr(alignment1$V2[3], 1, nchar(alignment1$V2[3])-2))
    stats$meth_cpgs[i] <- as.numeric(substr(meth$V2[1], 1, nchar(meth$V2[1])-1))
    stats$bisconv[i] <- 100-as.numeric(gsub("%", "", meth$V2[3]))
    
    cat("done\n")
  }
  
  return(stats)
}