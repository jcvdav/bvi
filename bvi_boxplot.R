bvi_boxplot <- function(bvi_scores){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  taxunits <- colnames(bvi_scores)[1]
  
  bvi_scores %>%
    select(-c(BVI, rBVI)) %>%
    gather(Sample, Score, -1) %>% 
    set_colnames(value = c("Spp", "Sample", "Score")) %>% 
    transform(Spp = reorder(Spp, Score)) %>% 
    ggplot(aes(x = Spp, y = Score)) +
    geom_boxplot(outlier.shape = NULL,
                 outlier.colour = "red",
                 outlier.alpha = 0,
                 fill = "gray",
                 alpha = 0.5) +
    geom_point(aes(color = Sample), size = 2) +
    theme_bw() +
    coord_flip() +
    labs(x = taxunits, y = "Score")

}