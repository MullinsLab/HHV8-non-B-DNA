library(ggplot2)
library(tidyr)
library(dplyr)
library(Biostrings)
library(ggthemes)

#This was used to make Figure 4
# SIST (Stress-Induced Structural Transitions) (Zhabinskaya, et al, 2015) for calculation of cruciform, denaturation and Z-DNA formation probabilities is run in Bash. R packages triplex (Hon, et al, 2013) and pqsfinder(Hon, et al, 2017) is used for calculating triplex and G-quadruplex scores, respectively.
# For functions used but not initialized here, refer to common_functions.R

# initialize variables for directories to be used
dir_in <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\input\\"
dir_out <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\output\\"
dir <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\"
docs <- "C:/Users/donju/OneDrive/Documents/"

bkpts <- read.csv(paste0(docs,"bkpts5.csv")) #table file containing list of breakpoints
#bkpts <- bkpts[18,]
HIPPOS_B <- readDNAStringSet(paste0(docs, "HIPPOSB.fasta"))
#source("C:/Users/donju/Dropbox/ggplot/B-DNA.R")
bin_size <- 500
g4_list <- vector(mode="list") #initialize list of G-quadruplex scores
tlist <- vector(mode="list") #initialize list of triplex scores

#Save FASTA files of 5 kb windows at given breakpoints in the table bkpts, to be input into SIST using Bash 
for (x in 1:length(bkpts$Coordinate)){
  s <- bkpts$Sequence[x]
  b <- bkpts$Coordinate[x]
  w <- DNAStringSet(HIPPOS_B[[s]][(b-(bin_size+1999)):(b+bin_size+2000)])   #to slice 5-kb windows from genome at given breakpoint
  n <- bkpts$Track[x]
  names(w) <- paste0(n, "_", b)
  writeXStringSet(w, paste0(docs, "windows/", n, "_", b, ".fasta"))

  #as well as calculate g4 scores at the windows  
  g4 <- g4_scores(b, bin_size, HIPPOS_B[[s]])#[101:900,] #can modify the window ssize here but record only inner 1000bp
  #g4 <- g4_scores(2500, 500, bkpt5kb.500shuff[[x]])#[2001:3000,]) #use this for permuted control
  g4$G4_scores <- g4$G4_scores/300 #to normalize to maximum of 300
  g4_list <- append(g4_list, list(g4))
  names(g4_list)[x] <- paste0(n, "_", b)
  # names(g4_list)[x] <- names(bkpt5kb.500shuff[x]) #for permuted control

  t3 <- triplex.search(w[[1]][2001:3000], max_len=100, max_loop=20, p_value=0.1, min_score = 15)
  # t3 <- triplex.search(bkpt5kb.500shuff[[x]][2001:3000], max_len=100, max_loop=20, p_value=0.1, min_score = 15) #for permuted control, looking at middle 1-kb window
  tlist <- append(tlist, list(t3))
  names(tlist)[x] <- paste0(n, "_", b)
  # names(tlist)[x] <- names(bkpt5kb.500shuff[x]) #for permuted control
}

tlist2 <- t3form(tlist, bin_size) #transform into list of length bin_size

#PLOT & SAVE EACH BREAKPOINT nonB-DNA PROBABILITIES USING ORIGINAL SEQUENCES
plotNBDNA <- function(sisto, g4_list, tlist2, seq, cutoff=FALSE, save_plot=FALSE, dir_out){
  df_list <- vector(mode="list", length=length(bkpts$Coordinate))
  for (t in sisto){ 
    if (seq=="obs"){
      df <- read.csv(paste0(dir_out, t), header=TRUE, sep="\t")[2001:3000,]
      dir_save <- paste0(docs,"obs_bkpts")
    }else if(seq=="perm") {
      df <- read.csv(paste0(dir_out, t), header=TRUE, sep="\t")[2001:3000,]
      dir_save <- paste0(docs,"bkpts_perm\\")
    }
    if (cutoff==TRUE){          #remove all probabilities less than 0.01
      df <- apply(df, 2, function(y) sapply(y, function(x) ifelse(x < 0.01, 0, x)))
    }
    
    # combine in G-quadruplex & triplex data, and name columns
    pic_lab <- gsub(".tsv", "", gsub("shuffled", "", t))
    print(paste("pic_lab assigned: ", pic_lab))
    df <- cbind(df, g4_list[pic_lab][[1]][,2])
    colnames(df)[5] <- "S_G4"
    df <- cbind(df, tlist2[pic_lab])
    colnames(df)[6] <- "S_triplex"
    
    #to change coordinate numbers, centering 0 at breakpoint
    numbering <- length(df$Position)/2        
    df$Position <- c((-numbering+1):numbering)
    
    if (save_plot==TRUE){
      p <- ggplot(df, aes(x=Position)) + geom_line(aes(y=P_melt, col="P_melt")) + 
        geom_line(aes(y=P_Z, col="P_Z")) + theme_light() + scale_color_discrete(name="Probabilities (P) \n or Scores (S)") +
        # geom_line(aes(y=P_cruciform, col="P_cruciform")) +
        theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text.x=element_text(size=18), axis.text.y=element_text(size=18), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=18), legend.direction="vertical", plot.title=element_text(size=20, face="bold")) + 
        guides(color=guide_legend(override.aes=list(size=3))) +
        geom_line(aes(y=S_G4, col="S_G4")) + geom_line(aes(y=S_triplex, col="S_triplex")) + geom_vline(xintercept=0) + ggtitle(pic_lab) + coord_cartesian(ylim=c(0,1)) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
      ggsave(paste0(pic_lab, ".png"), p, "png", dir_save)
    }
    
    x <- which(sisto==t)
    print(paste("t is", t, "; x is", x))
    df_list[x] <- list(df)
    names(df_list)[x] <- pic_lab
  }
  return(df_list)
}

sist_folder <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\output_bkpt5\\"

# Non-B DNA probabilities on the minimum region (mr) of over-coverage only
# Fig 4A
sist_mr <- list.files(paste0(dir, "\\output_mr")) #dir_out)
df_list_mr <- plotNBDNA(sist_mr, g4_list, tlist2, seq="obs", cutoff=FALSE, save_plot=FALSE, sist_folder)

# Non-B DNA probabilities on slices of windows centering on the 31 observed breakpoints
sist_bkpt5 <- list.files(paste0(dir, "\\output_bkpt5"))

df_list_all <- plotNBDNA(sist_bkpt5, g4_list, tlist2, seq="obs", cutoff=FALSE, save_plot=TRUE, sist_folder)
# df_list <- plotNBDNA(sist_bkpt7, g4_list, tlist2, seq="obs", cutoff=FALSE, save_plot=TRUE, sist_folder)

# Initialize data frames of SIST-calculated non-B DNA probablities (for permuted control dataset) 
sist_perm <- list.files(paste0(dir, "\\output_perm"))
# sist_perm <- list.files("C:\\Users\\donju\\Dropbox\\LAB\\SIST\\output_perm")
sist_folderp <- "\\\\themis.mullins.microbiol.washington.edu\\jcasjan\\SIST\\output_perm\\"
# sist_folderp <- "C:\\Users\\donju\\Dropbox\\LAB\\SIST\\output_perm\\"
df_listp <- plotNBDNA(sist_perm, g4_list, tlist2, seq="perm", cutoff=FALSE, save_plot=FALSE, sist_folderp)

#Combine the calculated probablities or scores of 5 non-B DNA types into one data frame (need to be of equal lengths)
comb_all <- function(df_list){
  melt <- comb_nbdna(df_list, "P_melt")
  zdna <- comb_nbdna(df_list, "P_Z")
  cruci <- comb_nbdna(df_list, "P_cruciform")
  g4 <- comb_nbdna(df_list, "S_G4")
  t3 <- comb_nbdna(df_list, "S_triplex")
  return(list(melt, zdna, cruci, g4, t3))
}


#Function for plotting the non-BDNA probailities or scores altogther on a (1 kb) window
# Fig 4B
plot_avg <- function(b_list, title="none"){
  mean_blist <- lapply(b_list, avg, mean) #max or mean or sum
  p_bdna <- cbind(mean_blist[[1]][1])
  for (a in mean_blist){
    p_bdna <- cbind(p_bdna, a[,2])
  }
  colnames(p_bdna) <- colnames(df_list[[1]])
  numbering <- length(p_bdna$Position)/2        
  p_bdna$Position <- c((-numbering+1):numbering)
  g <- ggplot(p_bdna, aes(x=Position)) + 
    geom_line(aes(y=P_melt, col="(P) melt region"), size=1) +
    geom_line(aes(y=P_Z, col="(P) Z-DNA"), size=1) +
    # geom_line(aes(y=P_cruciform, col="(P) cruciform")) + #no cruciforms found, so skipped
    geom_line(aes(y=S_G4, col="(S) G-quadruplex"), size=1) +
    geom_line(aes(y=S_triplex, col="(S) triplex DNA"), size=1)+
    geom_vline(xintercept=0) + 
    geom_hline(yintercept=gk18_avg_G4, linetype="dotdash", col="#00BFC4", size=1.2) + #for plotting the genome average of G4
    labs(y="probability/scores") + ggtitle(title) +
    #scale_color_discrete() +
    theme_classic() +
    theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=16), legend.title = element_blank(), legend.text = element_text(size=15), legend.direction="vertical", plot.title=element_text(size=20)) + guides(color=guide_legend(override.aes=list(size=4))) + coord_cartesian(xlim=c(-500,500), ylim=c(0,0.1)) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
  return(g)
}

b_list <- comb_all(df_list_all) #df_listp
bkplot_all <- plot_avg(b_list, "Observed")
# bkplot_all
# bkplot_mr
# bkplot_minreg
# bkplot_rand1
# bkplot_rand2
# bkplot_rand3
# bkplot_perm1

#Function to average all non-B DNA scores or probabilities in the window
avg_nbdna <- function(b_list, title="none"){
  mean_blist <- lapply(b_list, avg, mean) #b_list; max or mean or sum
  p_bdna <- cbind(mean_blist[[1]][1])
  for (a in mean_blist){
    p_bdna <- cbind(p_bdna, a[,2])
  }
  colnames(p_bdna) <- colnames(df_list[[1]])
  numbering <- length(p_bdna$Position)/2        
  p_bdna$Position <- c((-numbering+1):numbering)
  return(p_bdna)
}

# Plot comparison of the average probablities or scores between breakpoints in observed breakpoints, random control and permuted control
# For random and permuted controls, see controls_$_sims.R

avg_listr <- avg_nbdna(comb_all(df_listr)) #random control
avg_listp <- avg_nbdna(comb_all(df_listp)) #permuted control
avg_list_all <- avg_nbdna(comb_all(df_list_all)) #observed data

# For G4
# Genome-wide mean of G4 scores
gk18_avg_G4 <- mean(g4_gk18)/300 #mean of G4 scores along entire genome of strain GK18
avg_G4_list <- replicate(1000, gk18_avg_G4) #to create an array of length 1000 the genome-wide mean of G4 scores, to append to dataframe below
avg_g4 <- cbind(avg_list_all["Position"], avg_list_all["S_G4"], avg_listr["S_G4"], avg_listp["S_G4"], avg_G4_list)
names(avg_g4) <- c("Position", "obs", "rand", "perm", "avg")

# Fig 4D
g4_plot <- ggplot(avg_g4, aes(x=Position)) + 
  geom_line(aes(y=obs, col="Breakpoints observed"), size=1.25) +
  geom_line(aes(y=rand, col="Control random"), size=1) +
  geom_line(aes(y=perm, col="Control permuted"), size=1) +
  # geom_line(aes(y=avg, col="Genome average"), size=1, linetype="dashed") +
  geom_vline(xintercept=0) +
  scale_color_manual(values=c("#00bfc4", "lightslateblue", "indianred4", "#00bfc4")) +
  geom_hline(yintercept=gk18_avg_G4, linetype="dashed", col="#00bfc4", size=1.2) +
  labs(y="scores") + #ggtitle(title) +
  #scale_color_discrete() +
  theme_classic() +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=16), legend.title = element_blank(), legend.text = element_text(size=15), legend.direction="vertical", plot.title=element_text(size=20)) + guides(color=guide_legend(override.aes=list(size=4))) + coord_cartesian(xlim=c(-500,500), ylim=c(0,0.1)) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
g4_plot

#For Z-DNA
gk18_avg_zdna <- mean(gk18.sist2.df$P_Z) #mean of Z DNA probabilities along entire genome of strain GK18
avg_zdna_list <- replicate(1000, gk18_avg_zdna) #to create an array of length 1000 the genome-wide mean of Z DNA probabilities, to append to dataframe below
avg_zdna <- cbind(avg_list_all["Position"], avg_list_all["P_Z"], avg_listr["P_Z"], avg_listp["P_Z"], avg_zdna_list)
names(avg_zdna) <- c("Position", "obs", "rand", "perm", "avg")

# Fig 4C
zdna_plot <- ggplot(avg_zdna, aes(x=Position)) + 
  geom_line(aes(y=obs, col="Breakpoints observed"), size=1.5) +
  geom_line(aes(y=rand, col="Control random"), size=0.75) +
  geom_line(aes(y=perm, col="Control permuted"), size=0.75) +
  # geom_line(aes(y=avg, col="Genome average"), size=1, linetype="dashed") +
  scale_color_manual(values=c("#7CAE00", "lightslateblue", "indianred4")) +
  geom_vline(xintercept=0) + 
  geom_hline(yintercept=gk18_avg_zdna, linetype="dashed", col="#7CAE00", size=1.2) +
  labs(y="scores") + #ggtitle(title) +
  #scale_color_discrete() +
  theme_classic() +
  theme(axis.title.x=element_text(size=15), axis.title.y=element_blank(), axis.text.x=element_text(size=20), axis.text.y=element_text(size=16), legend.title = element_blank(), legend.text = element_text(size=15), legend.direction="vertical", plot.title=element_text(size=20)) + guides(color=guide_legend(override.aes=list(size=4))) + coord_cartesian(xlim=c(-500,500), ylim=c(0,0.1)) + scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))
zdna_plot
