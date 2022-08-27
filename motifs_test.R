# Counting Motifs and their association with breakpoints
# Tested to see if motifs identified in (Brown, 2014) are present and associated with the observed breakpoints
# Not enough sample size to be robust

docs <- "C:/Users/donju/OneDrive/Documents/"
HIPPOS_B <- readDNAStringSet(paste0(docs, "HIPPOSB.fasta"))

# To create DNA stringset of length 'ws' centered on coordinates (provided in data frame 'bkpts')
sliceWindow <- function(df=bkpts, ws=400){
  bkpt.dna <- DNAStringSet()
  for (x in 1:length(df$Coordinate)){
    s <- df$Sequence[x]
    b <- df$Coordinate[x]
    n <- df$Track[x]
    h <- ws/2
    bkpt.dna <- append(bkpt.dna, DNAStringSet(HIPPOS_B[[s]][(b-h+1):(b+h)])) 
    names(bkpt.dna)[x] <- paste0(n, "_",  b)
  } 
  return(bkpt.dna)
}

# To 
countMotifs <- function(bkpt.dna, mod="total"){
  bkpt.rc <- reverseComplement(bkpt.dna)
  tggtgg <- vcountPattern("tggtgg", bkpt.dna) + vcountPattern("tggtgg", bkpt.rc)
  cctcccct <- vcountPattern("cctcccct", bkpt.dna) + vcountPattern("cctcccct", bkpt.rc)
  cccag <- vcountPattern("cccag", bkpt.dna) + vcountPattern("cccag", bkpt.rc)
  tggag <- vcountPattern("tggag", bkpt.dna) + vcountPattern("tggag", bkpt.rc)
  gggct <- vcountPattern("gggct", bkpt.dna) + vcountPattern("gggct", bkpt.rc)
  aggag <- vcountPattern("aggag", bkpt.dna) + vcountPattern("aggag", bkpt.rc)
  df <- data.frame(tggtgg, cctcccct, cccag, tggag, gggct, aggag)
  rownames(df) <- names(bkpt.dna)
  if (mod=="at_least"){       #if counting only at least one motif present
    df[] <- lapply(df, function(x) ifelse(x>=1, 1, 0))
  }
  return(df)
}

randomCountMotifs <- function(bkpt_df, dnaset=HIPPOS, ws=200, mod="total"){
  bkpt.randseq <- DNAStringSet()
  wh <- ws/2
  for (x in 1:nrow(bkpt_df)){
    s <- bkpt_df$Sequence[x]
    b <- sample(wh:(length(dnaset[[s]])-wh), 1) # pick 1 random number
    dna <- DNAStringSet(dnaset[[s]][(b-wh+1):(b+wh)])
    bkpt.randseq <- append(bkpt.randseq, dna)
    #names(bkpt.randseq)[x] <- paste0(bkpt_df$Sequence[x], "_", b)
  }
  rand.motifs <- countMotifs(bkpt.randseq, mod)
  return(apply(rand.motifs, 2, sum))
}

permutedCountMotifs <- function(bd=bkpt1kb){
  shuff.dna <- shuffle(bd)
  shuff.motifs <- countMotifs(shuff.dna)
  shuff.motifv <- apply(shuff.motifs, 2, sum)
  return(shuff.motifv)
}

bkpt50bp <- sliceWindow(bkpts, 50)
bkpt1kb <- sliceWindow(bkpts, 1000)

bkpt.motifs200 <- countMotifs(sliceWindow(bkpts, 200), 200)
bkpt.motifs200.al <- countMotifs(sliceWindow(bkpts, 200), "at_least")
bkpt.motifs100.al <- countMotifs(sliceWindow(bkpts, 100), "at_least")
bkpt.motifs50 <- countMotifs(bkpt50bp, 50)
bkpt.motifv <- apply(bkpt.motifs100, 2, sum)
bkpt.motift <- sum(bkpt.motifv)

motifv.100full.al <- apply(bkpt.motifs100.al, 2, sum)
motifv.200full <- apply(bkpt.motifs200, 2, sum)
motifv.200full.al <- apply(bkpt.motifs200.al, 2, sum)
motifv.400full <- apply(bkpt.motifs400, 2, sum)

#Randomly sample & count motifs around however many breakpoints, 1000 times
rand.200full.array <- replicate(1000, randomCountMotifs(bkpts, HIPPOS_B, 200))#, "at_least"))
rand.array <- rand.200full.array
rand.100full.array.al
rand.200full.array
rand.200full.array.al
rand.400full.array
rand.400mr.array

rand.mean <- apply(rand.array, 1, mean)
rand.sd <- apply(rand.array, 1, sd)
rand200full.pval <- pnorm(motifv.200full, mean=rand.mean, sd=rand.sd, lower.tail=FALSE)

rand200full.pval.al
rand100full.pval
rand200full.pval
rand400full.pval
rand200mr.pval
rand100mr.pval

#Saving to table file
motif_df <- data.frame(motifv.200full, shuff.mean, shuff.sd, shuff.pval, rand.mean, rand.sd, rand200full.pval)
motif_df
write.csv(motif_df, file=paste0(docs, "breakpoint motifs pvalue.csv"))
write.csv(bkpt.motifs200.al, file=paste0(docs, "motif present within 100-bp.csv"))

bkpt200bp.shuff <- shuffle(bkpt200bp)

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

#Shuffle nucleotides of (28) observed breakpoints, 1000 times   
perm.array <- replicate(1000, permutedCountMotifs(bkpt200bp))
shuff.mean <- apply(perm.array, 1, mean)
shuff.sd <- apply(perm.array, 1, sd)
shuff.pval <- pnorm(motifv.200full, mean=shuff.mean, sd=shuff.sd, lower.tail=FALSE)
shuff.pval

perm.sum <- replicate(1000, sum(permutedCountMotifs(bkpt200bp)))
shuff.total.mean <- mean(perm.sum)
shuff.total.sd <- sd(perm.sum)
shuff.total.pval <- pnorm(bkpt.motift, mean=shuff.total.mean, sd=shuff.total.sd, lower.tail=FALSE)
shuff.total.pval