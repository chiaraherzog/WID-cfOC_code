
write_sens_spec <- function(type, score, label = "Ovarian cancer", scorename = "cfDNAme score"){
  library(epiR)
  require(gt)
  dat <- data.frame(matrix(nrow = 5, ncol = 3))
  colnames(dat) <- c(label, "Control", "All")
  rownames(dat) <- c(paste0(scorename, " positive"),
                     paste0(scorename, " negative"),
                     "Total",
                     "Sensitivity | Specificity",
                     "95% CI")
  dat[1,1] <- sum(type!="Control" & score == "pos")
  dat[1,2] <- sum(type=="Control" & score == "pos")
  dat[1,3] <- sum(score=="pos")
  dat[2,1] <- sum(type!="Control" & score != "pos")
  dat[2,2] <- sum(type=="Control" & score != "pos")
  dat[2,3] <- sum(score != "pos")
  dat[3,1] <- sum(type!="Control")
  dat[3,2] <- sum(type=="Control")
  dat[3,3] <- length(type)
  
  # SESP
  type <- factor(type, levels = c("Cancer", "Control"))
  score <- factor(score, levels = c("pos", "neg"))
  tab <- table(score, type)
  rval <- summary(epi.tests(tab))
  dat[4,1] <- paste0(round(rval$est[3]*100,1), "%")
  dat[5,1] <- paste0(round(rval$lower[3]*100,1), "-", round(rval$upper[3]*100,1), "%")
  dat[4,2] <- paste0(round(rval$est[4]*100,1), "%")
  dat[5,2] <- paste0(round(rval$lower[4]*100,1), "-", round(rval$upper[4]*100,1), "%")
  dat[4,3] <- ""
  dat[5,3] <- ""
  
  dat
  
  
  
  tab <- dat %>%
    gt::gt(rownames_to_stub=T) %>%
    text_transform(locations = cells_body(),
                   fn = function(x){
                     gsub("[.]", "·", x)
                   }
    ) %>%
    fmt_markdown(columns = everything()) %>%
    tab_style(
      cell_text(v_align = "top"),
      locations = cells_stub()
    )%>%
    tab_style(
      cell_text(align = "left"),
      locations = cells_stub()
    ) %>%
    gt::tab_options(column_labels.font.weight = "bold",
                    row_group.font.weight = "bold",
                    stub.font.weight = "bold",
                    data_row.padding = 2,
                    column_labels.font.size = 15,
                    table.font.size = 15,
                    row_group.padding = 2,
                    table.font.names = "Guardian Sans",
                    table.width = px(400),
                    row_group.border.right.width = px(10),
                    summary_row.padding = 2,
                    table.border.top.color = "white",
                    row_group.border.top.width = px(1),
                    row_group.border.bottom.width = px(1),
                    stub.border.width = px(0),
                    heading.title.font.size = 12)
  
  return(list(data = dat,
              table = tab))
  
}
