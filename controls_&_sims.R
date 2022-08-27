# CREATING CONTROLS AND SIMULATIONS
# Used in making Table 3, and for 1000 simulations for testing significance of association of observed breakpoints with G4


# Create permuted control sequences
# Function to shuffle or randomize the order a given DNA sequence 
shuffle <- function(bd){
  bd.shuff <- bd
  for (s in 1:length(bd)){
    bd.shuff[[s]] <- sample(bd[[s]], length(bd[[s]]))
  }
  return(bd.shuff)
}

# Shuffle only middle 1 kb of 5 kb windows sliced in plot_nB-DNA.R, then write out fasta files
bkpt5kb <- sliceWindow(bkpts, 5000)
bkpt5kb.500shuff <- bkpt5kb
for (s in 1:length(bkpt5kb)){
  bkpt5kb.500shuff[[s]][2001:3000] <- bkpt1kb.shuff[[s]]
  writeXStringSet(bkpt5kb.500shuff[s], paste0(docs, "windows_shuffled/", names(bkpt5kb.500shuff[s]), "shuffled.fasta"))
}
#The permuted sequences were then input on plot_nB=DNA.R



# Create random control sequences
# Pick random coordinates inside GK18 sequence, in a dataframe of same length as bkpts
bkptr <- data.frame(sample(1000:(length(GK18[1:137168])-1000), nrow(bkpts), replace=TRUE)) 
# Plot non-B DNA probablities and scores in 1 kb windows centering at random coordinates. SIST results for entire GK18 have been created above in 'dfo'
df_listr <- plot_bkpt(bkptr, dfo=dfo, b_window=500, write_plot=FALSE, template=GK18)



# Combine non-B DNA scores or probabilities by non-B DNA across observed breakpoints
# Function comb_nbdna defined in common_functions.R
melt <- comb_nbdna(df_list_all, "P_melt")
zdna <- comb_nbdna(df_list_all, "P_Z")
cruci <- comb_nbdna(df_list_all, "P_cruciform")
g4 <- comb_nbdna(df_list_all, "S_G4")
t3 <- comb_nbdna(df_list_all, "S_triplex")

#Function to sum up the probabilities/scores 200 bp before and after the breakpoint, 400 bp total
sumprob <-  function(df, lower=201, upper=600){
  vp <- vector()
  for (p in 2:length(df)){
    vp <- append(vp, sum(df[upper:lower,p]))
  }
  names(vp) <- colnames(df)[2:length(df)]
  return(vp)
}

sum_melt <- sumprob(melt, 301, 700)
sum_zdna <- sumprob(zdna, 401, 600)
sum_cruci <- sumprob(cruci, 301, 700)
sum_g4 <- sumprob(g4, 301, 700)
sum_t3 <- sumprob(t3, 301, 700)



# Combine across provided breakpoints in df_listr (random control) each non-B DNA probability and score mrand <- comb_nbdna(df_listr, "P_melt")
zrand <- comb_nbdna(df_listr, "P_Z")
crand <- comb_nbdna(df_listr, "P_cruciform")
g4rand <- comb_nbdna(df_listr, "S_G4")
t3rand <- comb_nbdna(df_listr, "S_triplex")

# Sum probabilities or scores the from the middle 400 bp of the 1 kb window 
sum_mrand <- sumprob(mrand, 301, 700)
sum_zrand <- sumprob(zrand, 301, 700)
sum_crand <- sumprob(crand, 301, 700)
sum_g4rand <- sumprob(g4rand, 301, 700)
sum_t3rand <- sumprob(t3rand, 301, 700)

# Test 2 sample sets for significance. Results in Table 3
t.test(sum_melt, sum_mrand)
t.test(sum_zdna, sum_zrand)
t.test(sum_cruci, sum_crand)
t.test(sum_t3, sum_t3rand)
t.test(sum_g4, sum_g4rand)

# Measure distance from breakpoint to nearest position with significant (>0.2) non-B DNA probability or score in either direction
# Function to measure distance
dist <- function(df_list, bdna="S_G4", min=0.2){
  d_list <- vector()
  for (dfw in df_list){
    o <- ifelse(dfw[500, bdna] > min, 0, 501)
    right <-  which(dfw[501:1000, bdna] > min)[1]
    right <- ifelse(is.na(right), 501, right)
    left <- which(dfw[499:1, bdna] > min)[1]
    left <- ifelse(is.na(left), 501, left)
    d_list <- append(d_list, min(o, right, left))
  }
  #names(d_list) <- names(df_list)
  return(d_list)
}

# Measure distance for observed breakpoints
melt_d <- dist(df_list_all, "P_melt")
zdna_d <- dist(df_list_all, "P_Z")
cruci_d <- dist(df_list_all, "P_cruciform")
g4_d <- dist(df_list_all, "S_G4")
t3_d <- dist(df_list_all, "S_triplex")

# Measure distance for random controls
melt_dr <- dist(df_listr1, "P_melt")
zdna_dr <- dist(df_listr1, "P_Z")
cruci_dr <- dist(df_listr1, "P_cruciform")
t3_dr <- dist(df_listr1, "S_triplex")
g4_dr <- dist(df_listr1, "S_G4")

# Test 2 sampel sets for significance. Results in Table 3
t.test(melt_d, melt_dr)
t.test(zdna_d, zdna_dr)
t.test(cruci_d, cruci_dr)
t.test(t3_d, t3_dr)
t.test(g4_d, g4_dr)



# Assess association of observed breakpoints with G4, by comparing G4 association with 1000 simulations of 31 random breakpoints on the GK18 genome

# function to sum up g4 scores in the middle 400 bp of random control breakpoints
simG4_sum <- function(points=31, df=dfog, seq=GK18){
  pts <- list(sample(1000:(length(seq)-1000), points, replace=TRUE))
  g4_list_sim <- lapply(pts[[1]], function(x) dfog[(x-499):(x+500),5])
  g4_sum_sim <- sapply(g4_list_sim, function(x) sum(x[301:700]))
  return(g4_sum_sim)
}
g4.sums.array <- replicate(1000, simG4_sum()) #create 1000 simulations
g4.sums.mean <- apply(g4.sums.array, 2, mean) #mean of 31 scores from each simulation
g4.sums.mean.mean <- mean(g4.sums.mean) #mean of the 31 means
g4.sums.sd <- sd(g4.sums.mean) #standard deviation from the 31 means
g4.sums.pval <- pnorm(mean(sum_g4), mean=g4.sums.mean.mean , sd=g4.sums.sd, lower.tail=FALSE) #sum_g4 is from the observed breakpoints, defined above
g4.sums.pval

# function to measure distances to nearest position with g4 scores above a certain score, in random control breakpoints
simG4_dist <- function(points=31, df=dfog, seq=GK18, min=0.1){
  pts <- list(sample(1000:(length(seq)-1000), points, replace=TRUE))
  g4_list_sim <- lapply(pts[[1]], function(x) dfog[(x-499):(x+500),5])
  d_list <- vector()
  for (window in g4_list_sim){
    o <- ifelse(window[500] > min, 0, 501)
    right <-  which(window[501:1000] > min)[1]
    right <- ifelse(is.na(right), 501, right)
    left <- which(window[499:1] > min)[1]
    left <- ifelse(is.na(left), 501, left)
    d_list <- append(d_list, min(o, right, left))
  }
  return(d_list)
}

g4.dist.array <- replicate(1000, simG4_dist(min=0.2)) #create 1000 simulations
g4.dist.mean <- apply(g4.dist.array, 2, mean) #mean of 31 distances from each simulation
g4.dist.mean.mean <- mean(g4.dist.mean)  #mean of the 31 means
g4.dist.sd <- sd(g4.dist.mean) #standard deviation from the 31 means
g4.dist.pval <- pnorm(mean(g4_d), mean=g4.dist.mean.mean , sd=g4.dist.sd, lower.tail=TRUE)
g4.dist.pval