cutoff <- function(data, type = "WID-cfOC"){
  if(type == "WID-cfOC"){
    efc144 <- 0.36973
    efc204 <- 0.26203
    efc228 <- 0.00029
  
    tmp <- data %>%
    mutate(cf_ovarian_144 = ifelse(data[,grepl("fully_methylated_EFC_144", colnames(data))] >= efc144, "pos", "neg")) %>%
    mutate(cf_ovarian_204 = ifelse(data[,grepl("fully_methylated_EFC_204", colnames(data))] >= efc204, "pos", "neg")) %>%
    mutate(cf_ovarian_228 = ifelse(data[,grepl("fully_methylated_EFC_228", colnames(data))] >= efc228, "pos", "neg")) %>%
    mutate(cf_ovarian_score = ifelse(cf_ovarian_144 == "pos" | cf_ovarian_204 == "pos" | cf_ovarian_228=="pos", "pos", "neg"))
  } else {
    efc93 <- 0.00008671026
    tmp <- data %>%
      mutate(cf_breast_score = ifelse(data[,grepl("fully_methylated_EFC_93", colnames(data))] >= efc93, "pos", "neg"))
  }
  
  return(tmp)
    
}