###
# Random Peptide Generator
###

computeAADist <- computePriorDist <- function(peptides) {
  #` Compute the distribution of amino acids in a peptide
  #` @param peptides A list of peptides wherein each list member is a character string of any length without spaces
  #` @example \code{computeAADist(list('AAKVCF', 'CWHGGDRQ'))}
  aa <- c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T")
  # freqs <- setNames(vector("list", length(aa)), aa)
  freqs <- rep(0, length(aa))
  totalAAs <- sum(sapply(peptides, function(peptide){ length(strsplit(peptide,'')[[1]]) }))
  for(peptide in peptides){
    for(p in strsplit(peptide,'')[[1]]){
      w <- which(aa == p)    
      freqs[[w]] <- freqs[w] + 1
    }
  }
  return(freqs/totalAAs)
}

randomAA <- function(prob=c(), aas=c()){
  #` @return A single random amino acid
  aa <- c("G","P","A","V","L","I","M","C","F","Y","W","H","K","R","Q","N","E","D","S","T")
  if(length(aas)>0) aa <- aas
  return(sample(aa, 1, prob=prob))
}

randomPeptides <- function(n.peps=1, min.aas=6, max.aas=9, prob=c()){
  #` Generate a list of length \code{n.peps} random peptides where a each peptide is a random length between min.aas and max.aas To generate random peptides of the same length set min.aas == max.aas .
  #` @param n.peps The length of the desired list of peptides
  #` @ param min.aas The minimum length of each peptide to be returned in the desired list of peptides
  #` @param max.aas The maximum length of each peptide to be returned in the desired list of peptides
  #' @param prob A probability distribution from which to sample amino acids
  #` @example \code{randomPeptides(n.peps=100, mon.aas=10, max.aas=10)}
  #` @example \code{randomPeptides(prob=computeAADist(peptides.positive))}
  if(n.peps<1) stop('n.peps must be greater than or equal to 1')
  if(min.aas > max.aas) {
    stop('min.aas cannot be less than max.aas')
  }else if(min.aas == max.aas){
    to.sample <- max.aas
  }else{
    to.sample <- min.aas:max.aas    
  } 
  ## This is extremely slow
  sapply(1:n.peps, function(x){
    n.aas <- sample(to.sample, 1)
    paste(sapply(1:n.aas, function(y){
      randomAA(prob=prob)
    }), collapse='')
  })
}


###
# Peptide descriptor generator
###

makeDescriptors <- function(peptides, descriptors) {
  #` Create a dataframe of descriptors with rows as peptides and columns as descriptors
  #` @param peptides Each peptide for which several descriptors will be calculated
  #` @param descriptors A list of method names (from the Rcpi package) for calculating protein descriptors
  out <- lapply(peptides, function(peptide){
    sapply(descriptors, function(descriptor){
      do.call(descriptor, list(peptide))
    })
  })
  t(sapply(out, unlist))
}

