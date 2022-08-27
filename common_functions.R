# Install and load the following R packages:
# pqsfinder, triplex, Biostrings
# This script initializes many functions used in the study
# SIST (Stress-Induced Structural Transitions) (Zhabinskaya, et al, 2015) for calculation of cruciform, denaturation and Z-DNA formation probabilities is run in Bash. R packages triplex (Hon, et al, 2013) and pqsfinder(Hon, et al, 2017) is used for calculating triplex and G-quadruplex scores, respectively.

library(pqsfinder)
library(triplex)
library(Biostrings)

#Initialize some directories, and variables for testing
docs <- "C:/Users/donju/OneDrive/Documents/"
fasta <- readDNAStringSet(paste0(docs, "GK18.txt"))
GK18 <- DNAString(fasta[[1]])#[1:125000]
HIPPOS <- readDNAStringSet(paste0(docs, "HIPPOS.fasta")) #HIPPOS.fasta contains genome seuquences of isolates to be analyzed
windows <- list.files(paste0(docs, "windows/"))
dir_in <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\input\\"

bkpts <- read.csv(paste0(docs,"bkpts2.csv"))
bkpt <- sscs_bkpts$Minimum[12]
w_size <- 400
bkpt <- 16450
f <- pqsfinder(GK18)

# Adapt the G4 prediction tool pqsfinder for use here; refer to the pqsfinder manual for the modifiers used e.g. 'deep'
# To output only a list of G4 scores for each base in a given sequence window

g4_finder <- function(bkpt, w_size, seq){
  w <- seq[(bkpt - w_size + 1):(bkpt + w_size)]
  pqs <- pqsfinder(w, deep=TRUE) #, min_score=20)
  return(pqs)
}

g4_deep <- function(bkpt, w_size, seq){
  pqs <- g4_finder(bkpt, w_size, seq)
  pqsd <- density(pqs)
  df <- data.frame((-w_size+1):w_size, pqsd)
  names(df) <- c('coordinates','G4_density')
  return(df) 
}

g4_scores <- function(bkpt, w_size, seq){
  pqs <- g4_finder(bkpt, w_size, seq)
  pqsS <- maxScores(pqs)
  df <- data.frame((-w_size+1):w_size-1, pqsS)
  names(df) <- c('coordinates','G4_scores')
  return(df) 
}

#COMBINE ALL G4 DENSITIES/SCORES
comb_g4 <- function(bkpt_df, operation=g4_deep, w_size=400, dna=GK18){
  df <- operation(bkpt_df$Minimum[1], w_size, dna)
  for (breakpoint in bkpt_df$Minimum[1:(nrow(bkpt_df))]){
    df_add <- operation(breakpoint, w_size, dna)
    df <- cbind(df, df_add[,2]) #G4_density or G4_scores
  }
  #print(paste("length(df) is", length(df)))
  #print(paste("nrow(bkpt_df) is", nrow(bkpt_df)))
  colnames(df) <- c("coordinates", sub("\\+hu_map.vcf","",bkpt_df$Track.Name[1:(nrow(bkpt_df))])) # +1 column total
  return(df)
}

comb_nbdna <- function(df_list, bdna){
  dfc <- data.frame(df_list[[1]]$Position, df_list[[1]][bdna])
  names(dfc) <- c("Position", names(df_list)[1])
  for (b in 2:length(df_list)){
    dfc <- cbind(dfc, df_list[[b]][bdna])
    names(dfc)[length(dfc)] <- names(df_list)[b]
  }
  return(dfc)
}

#PLOT AVERAGE OF COMBINED SCORES/DENSITIES
avg <- function(df_binned, FUN=mean){  #max or rowMeans or sum
  avg_v <- apply(df_binned[,-1], 1, FUN)
  plot_depth <- data.frame(df_binned[,1], avg_v)
  #colnames(plot_depth)[1] <- "Position"
  return(plot_depth)
}

#BIN G4 SCORES (Optional)
bin_by <- function(v, bin_size){
  ms_vector <- vector()
  for (i in (1:(length(v)/bin_size))){
    c <- i*bin_size
    ms <- max(v[(c-(bin_size-1)):c])
    ms_vector <- append(ms_vector, ms)
  }
  if ((length(v)%%bin_size != 0)){
    r <- length(v)%%bin_size
    ls <- length(v)
    last <- max(v[(ls-r):ls])
    ms_vector <- append(ms_vector, last)
  }
  return(ms_vector)
}

bin_df <- function(df, bin_size=40){
  df_binned <- seq(df[,1][1], df[,1][nrow(df)], bin_size)
  for (b in 2:ncol(df)){
    full_v <- df[,b]
    binned_v <- bin_by(full_v, bin_size)
    df_binned <- cbind(df_binned, binned_v)
  }
  df_binned <- data.frame(df_binned)
  colnames(df_binned) <- colnames(df)
  return(df_binned)
}

#SAVE G4 PLOTS FOR INDIVIDUAL BREAKPOINT
for (breakpoint in bkpts$Minimum){
  g4 <- g4_deep(breakpoint, w_size, GK18) #g4_scores or g4_deep
  source <- bkpts[which(bkpts$Minimum==breakpoint),"Track.Name"]
  source <- sub("\\+hu_map.vcf", "_", source)
  picture <- paste(source, "position", breakpoint)
  g4_plot <- ggplot(g4, aes(x=coordinates, y=G4_density)) + geom_line() + geom_vline(xintercept=0) + ggtitle((picture)) #y=G4_density or G4_scores
  ggsave(paste0(picture, ".png"), g4_plot, "png", paste0(docs,"g4-depth")) #"g4-depth" or g4-scores
}

#PLOT ONE G4 DENSITY
# df <- g4_deep(bkpts$Minimum[42], w_size, GK18)
# source <- bkpts[42,"Track.Name"]
# source <- sub("\\+hu_map.vcf", " ", source)
# picture <- paste0(source, " position: ", bkpts$Minimum[42])
# ggplot(df, aes(x=coordinates, y=G4_density)) + geom_line() + geom_vline(xintercept=0) + ggtitle((picture))

# TRANSFORM OUTPUT OF triplex.search() INTO A LIST OF TRIPLEX SCORES FOR PLOTTING
# Triplex.search() normally outputs in a text format, refer to its manual
t3form <- function(tlist, window){          #for 2-kb window, recording only inner 800 bp
  tlist2 <- vector(mode="list")
  for (x in 1:length(tlist)){       
    v <- rep(0, (2*window))                 #create a vector of just 0s as column
    if (length(tlist[[x]])==0){             #if TriplexViews is empty or NONE
      tlist2 <- append(tlist2, list(v))     #add the vector of just 0s right away
    }else{                                  #else if there are triplexes
      for (y in (1:length(tlist[[x]]))){    #for each triplex in a window
        c <- start(tlist[[x]])[y]           #c is start of triplex y in window "t[[x]]"
        # if (c > 500 & c < 1501){          #only if triplex is inside the middle 800-bp window
          # shifted <- c-500                #adjust c position for plotting in 800-bp window
          # v[shifted:(shifted+width(tlist[[x]])[y])] <- score(tlist[[x]])[y]/50
          #replace 0s to scores at corresponding positions
        # }
        v[c:(c+width(tlist[[x]])[y])] <- score(tlist[[x]])[y]/50
      } #after looping though all contents of tlist[[x]], only then add vector
      tlist2 <- append(tlist2, list(v))  #add the vector of just 0s with/without change
    }
  }
  tlist2 <- data.frame(tlist2)
  names(tlist2) <- names(tlist)
  return(tlist2)
}

# TO CREATE A TABLE OF SIST RESULTS OF THE WHOLE KSHV GENOME GK18 SEQUENCE, FOR SLICING WINDOWS FOR 'RANDOM CONTROL' DATASET
# SIST is run on Bash. Because SIST is only able to efficiently process a maximum of 5 kb, the 137-kb KSHV genome GK18 was first sliced (in Bash) into windows of 5 kb, overlapping by 2.5 kb. Then, the function all_sist below combines the resulting table of scores per base back into one table, to list the scores at each position in the 137 kb genome. The higher of the scores will be listed where the windows overlap.

#Mini-function to evaluate the higher score in overlaps of windows, to be nested inside function 'all_sist' below
higher <- function(pos, dfn, dfo){
  # print("pos is")
  # print(pos)
  dfo[dfo$Position==pos,]["P_melt"] <- max(dfn[dfn$Position==pos,]["P_melt"], dfo[dfo$Position==pos,]["P_melt"])
  # print("success1")
  dfo[dfo$Position==pos,]["P_Z"] <- max(dfn[dfn$Position==pos,]["P_Z"], dfo[dfo$Position==pos,]["P_Z"])
  # print("success2")
  dfo[dfo$Position==pos,]["P_cruciform"] <- max(dfn[dfn$Position==pos,]["P_cruciform"], dfo[dfo$Position==pos,]["P_cruciform"])
  print(pos)
  return(dfo)
}

#Mini-function to properly sort in increasing order the overlapping windows from SIST results (named with coordinate numbers ######.tsv)
list_order <- function(dir=dir_out){
  list_sist <- list.files(dir)
  reordering1 <- sapply(".tsv", gsub, "",  x=list_sist) #removing the ".tsv" from the filenames
  reordering2 <- reordering1[order(nchar(reordering1), reordering1)] 
  return(paste0(reordering2, ".tsv"))
} 

all_sist <- function(dir, saveTable=TRUE){
  sist_output <- list_order(dir)
  dfo <- read.csv(paste0(dir,sist_output[1]), sep="\t")       #initialize dfo with first window
  coord0 <- as.numeric(gsub(".tsv", "", gsub("U048-E_w", "", sist_output[1]))) #change GK18ext... character  to number coordinate
  print(paste0("coord0 is ", coord0))
  dfo$Position <- as.numeric(dfo$Position) + coord0 - 1
  for (window in sist_output[2:(length(sist_output)-1)]) {    #except last window, see further below
    dfn <- read.csv(paste0(dir, window),  sep="\t")
    coord <- as.numeric(gsub(".tsv", "", gsub("U048-E_w", "", window))) #change GK18ext... character  to number coordinate
    print(paste0("coord is ", coord))
    dfn$Position <- as.numeric(dfn$Position) + coord - 1
    #len <- length(dfn$Position)
    for (pos in dfn$Position[1:2500]){
      dfo <- higher(pos, dfn, dfo)
    }
    dfo <- rbind(dfo, dfn[2501:5000,])
  } 
  last_window <- sist_output[length(sist_output)]                  #adding last window
  dfn <- read.csv(paste0(dir, last_window),  sep="\t")
  coord <- as.numeric(gsub(".tsv", "", gsub("U048-E_w", "", last_window))) #convert to character name to coordinate number
  dfn$Position <- as.numeric(dfn$Position) + coord - 1
  if (length(dfn$Position) <= 2500) {
    print(paste("last window has less than 2500, is", length(dfn$Position)))
    for (pos in dfn$Position[1:length(dfn$Position)]){
      dfo <- higher(pos, dfn, dfo)
    }
  } else {
    print(paste("last window has MORE than 2500, is", length(dfn$Position)))
    for (pos in dfn$Position[1:2500]){
      dfo <- higher(pos, dfn, dfo)
    }
    dfo <- rbind(dfo, dfn[2501:length(dfn$Position),])
  }
  if (saveTable==TRUE){
    write.csv(dfo, paste0(docs, "total_sist2.csv"), row.names=FALSE)
  }
  return(dfo)
}

dir <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\output\\"
sist2 <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\output_GK18\\" #not including TR
gk18.sist2.df <- all_sist(dir)
min_reg.df <- all_sist(sist2)
