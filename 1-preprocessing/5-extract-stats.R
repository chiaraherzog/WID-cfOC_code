# Extract stats from run
#---------------------------------------
# Compile phenotypic and epigenetic data from all plates
# including:
#   - stats
#   - gDNA
#   - cfDNAme patterns
# Edit 6 Apr: adding CA125 UCL data
# Edit 18 Aug: adding updated gDNA data
#---------------------------------------

library(dplyr)
library(tidyverse)
library(stringr)

dir <- "<path to results>"

# Extract mapping/stats
report <- list.files(paste0(output, "reports/"))
report <- report[!grepl("splitting", report)]
stats <- data.frame(matrix(nrow = length(report),
                           ncol = 4))
colnames(stats) <- c("name", "reads", "mapping", "meth_cpgs")

for (i in 1:length(report)){
  # fix the name 
  name <- case_when(grepl("NA", report[i]) ~ paste(substr(report[i], 9, 18), stringr::str_split(report[i], "_", simplify = TRUE)[,9], sep = "_"),
                    grepl("NTC", report[i]) ~ paste(substr(report[i], 9, 19), stringr::str_split(report[i], "_", simplify = TRUE)[,9], sep = "_"),
                    grepl("IP12", report[i]) ~ paste(substr(report[i], 9, 20), stringr::str_split(report[i], "_", simplify = TRUE)[,9], sep = "_"),
                    grepl("UKFOCS", report[i]) ~ str_split(report[i], "_", simplify = TRUE)[,3])
  
  # read in mapping rate
  alignment1 <- suppressMessages(read.table(report[i], skip = 6,
                                            sep = "\t",
                                            nrow = 6))
  meth <- suppressMessages(read.table(report[i], skip = 35,
                                      sep = "\t",
                                      nrow = 3))
  
  stats$name[i] <- name
  stats$mapping[i] <- as.numeric(substr(alignment1$V2[3], 1, nchar(alignment1$V2[3])-2))
  stats$reads[i] <- as.numeric(alignment1$V2[2])
  stats$meth_cpgs[i] <- as.numeric(substr(meth$V2[1], 1, nchar(meth$V2[1])-1))
  
}
rm(meth, alignment1, name, report);gc()

# Stats is combined with pheno file