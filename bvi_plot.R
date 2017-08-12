bvi_plot <- function(bvi_scores){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  taxunits <- colnames(bvi_scores)[1]
  
  bvi_scores %>%
    select(-c(BVI, rBVI)) %>%
    gather(Sample, Score, -1) %>% 
    set_colnames(value = c("Spp", "Sample", "Score")) %>% 
    ggplot(aes(x = Sample, y = Score, fill = Spp)) +
    geom_col(position = "fill", color = "black") +
    theme_bw()
}