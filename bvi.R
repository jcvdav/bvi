#' Biological Value Index
#' 
#' @description Calculates the Biological Value Index (Sanders, 1960) with modifications from Loya-Salina & Escofet (1990). It also provides a relative index, that allows for comparisons.
#' 
#' @param data a data.frame in wide format. One column (usually the first one) must contain taxonomic units (e.g. species). Then, each column represents a sample, and the cells contain observed abundances or standardized densities (i.e. abundance / area or abundance / time)
#' @param taxunits A character vector that contains the name of the column with the taxonomic units. Default value is "Spp".
#' @param p The cutoff percentage to use, expressed relative to 1 (i.e. for a cutoff of 95 percent use p = 0.95)
#' @param other a logical that indicates if species that do not contribute to p should be added together in a category "other". If FALSE, they are displayed individually.
#' 
#' @return results An object of class data.frame containing a column with species name (spp), the scores for each species in each sample, a column (BVI) with the Biological Value Index, and a column (rBVI) with the Relative Biological Value Index. The results can be directly passed to bvi_plot to get a graphical representation.
#' 
#' @seealso bvi_plot
#' 
#' @author Villasenor-Derbez, J.C.
#' 
#' @export

bvi <- function(data, taxunits = "Spp", p = 0.95, others = T){
  
  # Make sure that taxunits is present
  if(!taxunits %in% colnames(data)){
    stop(paste0("Variable '", taxunits, "' not present in supplied data"))
  }
  
  #Test that the cutoff percentage is between 0 and 1
  if (p > 1 | p < 0){
    stop("Your cutoff percentage must be between 0 and 1")
  }
  
  #Test that the data passed is a matrix, else it transforms it
  spp <- data %$%
    get(taxunits) %>%
    as.character()
  
  data %<>%
    select(-matches(taxunits)) %>%  # Delete the column that contains the taxunits
    as.matrix() %>%                 # Conver to a matrix
    set_rownames(value = spp)       # Assign taxunits as rownames to the matrix
  
  data[is.na(data)] <- 0            # Replace all NAs with 0s
  
  n_samples <- dim(data)[2]         # Get number of samples in data
  n_taxunits <- length(spp)         # Get the number of species (or taxonomic units) in data
  
  prop <- prop.table(data, 2)       # Calculate matrix of relative abundances by sample
  
  ni <- rowSums(data)/sum(data)     # Calculate total relative abundances by species
  
  nsp <- rep_len(0, n_samples)      #vector with zeros of length=number of samples
  
  #The folowing cycle calculates the number of species (nsp) needed to reach the cutoff percentage (p)
  for (j in seq(1, n_samples)){
    i <- 1                             #Re-set i at 1
    nij <- sort(prop[,j],decreasing=T) #Sort the column j in a decreasing order
    ni_p <- nij[i]                     #ni_p will be the cumulative relative abundance. Here it is the largest value of the sorted column
    
    #This conditional keeps adding species' relative abundance to ni_p until the cutoff percentage p is reached
    while(!ni_p >= p){
      ni_p <- ni_p + nij[i + 1]        #ni_p increases with the addition of every species [i+1]
      i <- i + 1                       #i keeps track of the number of species that have been added to reach p
    }
    nsp[j] <- i                        #For sample j, the number of species to reach p is i
  }
  
  #Extract the maximum number of species needed and overwrite nsp
  nsp <- max(nsp)
  
  #Pre-demension the matrix of scores so that it is equal to data, but with scores of 0
  scores <- matrix(0, nrow = n_taxunits, ncol = n_samples) %>% 
    set_colnames(colnames(data)) %>% 
    set_rownames(rownames(data))
  
  #Cycle iterates across samples (j)
  for (j in seq(1, n_samples)){
    score <- unname(prop[, j])                #Extract abundances of sample j
    un <- sort(unique(score), decreasing = T) #Extract unique vector of abundances to correct for ties, it has length <= length(score)
    N <- nsp                                  #N has the maximum number of species needed, and is modified later in the correction for ties
    score_j <- matrix(0, n_taxunits)          #Set initial score_j = 0 and with size of score
    
    # Cycle iterates across species (i) either:
    # from 1 to nsp (max number of species needed to reach p) or 
    # length(un) if length(un) < nsp
    
    for (i in seq(1, min(length(un)-1,nsp))){
      score_j[score == max(un)] <- N #The score N is assigned to species i (i.e. all those that have an abundance of max?score)
      N <- N - sum(score == max(un))   #Correction for ties
      un[un == max(un)] <- 0         # Make the previous values of max(un) = 0, so that the next max(un) is updated
    }
    scores[,j] <- score_j            #Assign scores of sample j to the score matrix
  }
  
  BVI <- rowSums(scores)       #Calculates BVI by adding scores for each species across samples
  rBVI <- BVI / sum(BVI) * 100 #Calculate relative BVI that can be used to compare
  
  #Build preeliminar results into a data.frame
  results <- data.frame(Spp = rownames(data),
                       scores,
                       BVI,
                       rBVI,
                       stringsAsFactors = F) %>% 
    arrange(desc(rBVI))

  #If other=TRUE, adds species that do not contribute to p in a single raw
  if(others){
    results_b <- results[1:nsp,]
    scores_b <- colSums(results[(nsp+1):n_taxunits, 1:n_samples+1])
    results_b[nsp+1, 1:n_samples+1] <- scores_b
    results_b$Spp[nsp+1] <- "Others"
    results_b$BVI[nsp+1] <-  sum(scores_b)
    results_b$rBVI[nsp+1] <-  sum(scores_b) / sum(results_b$BVI)
    results <- results_b
  }
  
  colnames(results) <- c(taxunits, colnames(data), "BVI", "rBVI")
  
  #Return the results
  return(results)
}


