# Calculate odds of a random 17.5 kb slice of the KSHV genome GK18 will include the 'minimum region' frequently over-covered in whole viral genome sequencing of many KS tumors. 17.5 kb is the average size of the overrepresented regions.

# Initializing variables to be used
docs <- "C:/Users/donju/OneDrive/Documents/"
fasta <- readDNAStringSet(paste0(docs, "GK18.txt"))
GK18 <- DNAString(fasta[[1]]) 
bkpts <- read.csv(paste0(docs,"bkpts5.csv"))

# Mini-fuction to evaluate if a segment contains the K5 position, to be nested in function 'bootstrap' below
have_K5T1 <- function(centerpoint, segment=17500){ #centerpoint of over-coverage region, segment is its length
  dup <- segment/2
  if (((centerpoint-dup) < 25000) & ((centerpoint+dup) > 30000)){  #if read over-covered segment's 5'end is less than 2500 and 3'end more than 3000
    return(TRUE) # then return true
  } else {
    return(FALSE)
  }
}

# Sample random 10 numbers between coordinates [segment] and length(GK18)-[segment], with replacement; then replicate this 10 random number set 100 times; outputs a matrix of 10 columns and 100 rows
bootstrap <- function(t=100, sampling=10, segment=17500, s=GK18){
  sim_sv <- replicate(t, sample(segment:(length(s)-segment), sampling, replace=TRUE)) 
  svK5 <- apply(sim_sv, c(1,2), have_K5T1)
  svp <- apply(svK5, 2, mean)
  # svp <- vector()
  # for (c in 1:ncol(sim_sv)){                    #for each row of matrix
  #   l <- sapply(sim_sv[,c], have_K5T1, segment) #do have_K5T1 on list of 10 random numbers on the row
  #   p_K5T1 <- mean(l)                           #average the 10-results
  #   svp <- append(svp, p_K5T1)                  #add to vector svp
  # }
  return(svp)
}

svp <- bootstrap(10000, sampling=10, segment=176000)   

# What is the probability, that given 10 tumors with random 17.kb windows, all 10 windows contain K5-K6 region? 
pnorm(1, mean=mean(svp), sd=sd(svp), lower.tail=FALSE)