bvi_col <- function(bvi_scores){
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  taxunits <- colnames(bvi_scores)[1]
  
  maxrBVI <- max(bvi_scores$rBVI)*1.1/100
  
  transform(bvi_scores, Spp = reorder(Spp, rBVI)) %>% 
    ggplot(aes(x = Spp, y = rBVI/100)) +
    geom_col(fill = "gray", color = "black", alpha = 0.5) +
    labs(x = taxunits, y = "relative BVI") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, maxrBVI), labels = scales::percent) +
    coord_flip() +
    theme_bw()
}