plot_sens_spec <- function(type, scorelist, title = "Title", print = T){
  
  require(epiR)
  require(ggsci)
  
  cols <- MetBrewer::met.brewer("Isfahan1", n = 8)
  
  # relevel for tab
  type <- factor(type, levels = c("Cancer", "Control"))
  
  # create df
  df <- data.frame(matrix(nrow = length(scorelist),
                          ncol = 7))
  colnames(df) <- c("Score",
                    "Sensitivity",
                    "Specificity",
                    "Lower_Sens",
                    "Higher_Sens",
                    "Lower_Spec",
                    "Higher_Spec")
  
  
  for (i in 1:length(scorelist)){
    scorelist[[i]] <- factor(scorelist[[i]], levels = c("pos", "neg"))
    
    tab <- table(scorelist[[i]], type)
    rval <- summary(epi.tests(tab))
    
    df$Score[i] <- names(scorelist)[i]
    df$Sensitivity[i] <- rval[3,2]
    df$Lower_Sens[i] <- rval[3,3]
    df$Higher_Sens[i] <- rval[3,4]
    df$Specificity[i] <- rval[4,2]
    df$Lower_Spec[i] <- rval[4,3]
    df$Higher_Spec[i] <- rval[4,4]
  }
  
  
  plot <- df %>%
    ggplot(aes(x = Specificity,
               y = Sensitivity,
               colour = Score)) +
    geom_point()  +
    geom_pointrange(aes(xmin = Lower_Spec,
                      xmax = Higher_Spec)) +
    geom_pointrange(aes(ymin = Lower_Sens,
                      ymax = Higher_Sens)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_bw() +
    theme(legend.position = "top",
          text = element_text(family="Guardian Sans")) +
    scale_colour_manual(values = cols[c(3,5,1,7,2,4,6)],
                        name = "") +
    ggtitle(title)
  
  data <- list(plot = plot,
               df = df)
    
  if(print == T){
    print(plot)
  } else {
    return(data)
  }
  
}
