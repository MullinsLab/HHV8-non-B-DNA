b_window <- 400 #size of window
docs <- "C:/Users/donju/OneDrive/Documents/"
fasta <- readDNAStringSet(paste0(docs, "GK18.txt")) #removed the TR sequences
GK18 <- DNAString(fasta[[1]])#[1:125000]
U048E <- DNAString(readDNAStringSet(paste0(docs, "U048-E.fasta"))[[1]])
# source("C:/Users/donju/Dropbox/ggplot/B-DNA.R")

#write out sequences of random slices of GK18 
for (x in bkptr1[,1]){
  w <- DNAStringSet(GK18[(x-999):x+1000])
  names(w) <- paste("GK18", x)
  writeXStringSet(w, paste0(docs, "windows_r/", "GK18_", x, ".fasta"))
}

#triplexes for whole of GK18
t <- triplex.search(GK18, max_len=100, max_loop=20, min_score = 10)
tv_gk18 <- rep(0, length(GK18))
for (y in (1:length(t))){             #for each triplex in t
    c <- start(t)[y]                  #c is start of triplex y
    wid <- width(t)[y]
    tv_gk18[c:(c+wid)] <- score(t)[y]#/50
}

#MELT, Z-DNA, CRUCI, G4 OF ENTIRE GK18 GENOME
g4_gk18 <- maxScores(pqsfinder(GK18[1:137168], deep=TRUE))

# g4_dr1dr2_007 <- maxScores(pqsfinder(dr1dr2[[1]], deep=TRUE))

# dfo <- read.csv(paste0(docs,"total_sist.csv"))
dfo <- gk18noTR.sist.df
dfo$Position <- dfo$Position-2403 #do adjust dfo coordinates (from no TR) to GK18 coordinates
dfog <- cbind(dfo, g4_gk18/300)
dfog <- cbind(dfog, tv_gk18/50)
names(dfog)[5] <- "S_g4"
names(dfog)[6] <- "S_triplex"
dfog_bin <- bin_df(dfog, 100)
ggplot(dfog_bin, aes(x=Position)) + scale_x_continuous(breaks = seq(0,140000,10000), expand = c(0,0)) +
  ggtitle(("40 random breakpoints (|) on GK18")) +
  geom_linerange(data=bkptr, aes(x=bkptr[,1], ymax=0.95, ymin=0.85)) + 
  # data=bkpts, aes(x=Coordinate ; data=bkptr, aes(x=bkptr[0,1]
  labs(y="probability (P) or score (S)") + 
  #coord_cartesian(xlim=c(20000,50000)) #+
  geom_line(aes(y=S_g4, col="G4_scores")) +
  geom_line(aes(y=P_melt, col="P_melt")) + 
  geom_line(aes(y=P_cruciform, col="P_cruciform")) +  
  geom_line(aes(y=P_Z, col="P_Z")) + 
  geom_line(aes(y=S_t3, col="S_triplex")) #+
  #guide_legend(keywidth=unit(2,"mm"))

#PLOT & SAVE EACH RANDOM BREAKPOINT nonB-DNA probabilities FROM GK18
plot_bkpt <- function(bkptr, dfo=dfo, b_window=400, write_plot=FALSE, template=GK18){
  df_listr <- vector(mode="list", length=length(bkptr))
  for (x in (1:nrow(bkptr))){
    breakpoint <- bkptr[x,1]
    dfw <- subset(dfo, Position>=(breakpoint-b_window) & Position<=(breakpoint+b_window-1))
    #add 2403 to adjust that dfo was made on GK18 with 3X 801-bp TR on both ends
    dfw$Position <- c((breakpoint-b_window):(breakpoint+b_window-1))

    g4r <- g4_scores(breakpoint, b_window, template)#[101:900,] #calc "b_window"-sized window, #but record only inner 800-bp
    g4r$G4_scores <- g4r$G4_scores/300
    dfw <- cbind(dfw, g4r$G4_scores)
    colnames(dfw)[5] <- "S_G4"

    window <- template[(breakpoint-b_window+1):(breakpoint+b_window)]
    t3r <- triplex.search(window, max_len=100, max_loop=20, p_value=0.1, min_score = 15)
    tlist2 <- t3form(list(t3r), b_window)
    dfw <- cbind(dfw, tlist2[1])
    colnames(dfw)[6] <- "S_triplex"

    pic_lab <- paste("GK18 position", breakpoint)
    if(write_plot==TRUE){
      p <- ggplot(dfw, aes(x=Position)) + 
        geom_line(aes(y=P_melt, col="(P) melt region"), size=1) + geom_line(aes(y=P_Z, col="(P) Z-DNA"), size=1) + geom_line(aes(y=S_G4, col="(S) G-quadruplex"), size=1) + geom_line(aes(y=S_triplex, col="(S) triplex DNA"), size=1) + 
        scale_color_discrete(name="Probabilities (P) or Relative Scores (S)") +
        ggtitle((pic_lab)) + theme_light() + 
        # geom_line(aes(y=P_cruciform, col="P_cruciform")) +
        coord_cartesian(xlim=c((breakpoint-b_window+1),(breakpoint+b_window)), ylim=c(0,1)) + 
        theme(axis.title.x=element_text(size=18), axis.title.y=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=14), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=18), legend.direction="horizontal", plot.title=element_text(size=18), legend.position = "top") +
        #+ geom_vline(xintercept=0)
        geom_hline(yintercept=gk18_avg_G4, linetype="dotdash", col="#00BFC4", size=1.2) +
        geom_linerange(data=bkpt_mr_all, aes(x=bkpt_mr_all[,1], ymax=-0.15, ymin=0), size=1) +
        scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0.25,0))
      # ggsave(paste0(pic_lab, ".png"), p, "png", paste0(docs,"min_agreement"))
    }
    df_listr[x] <- list(dfw)
    names(df_listr)[x] <- pic_lab
    print(pic_lab)
  }
  return(df_listr)
}

bkpt_mr <- data.frame(c(19426, 22800, 24500, 25318, 30219, 33741))
bkpt_mr_all <- data.frame(c(22442, 22991, 24550, 25099, 25153, 25632, 27832, 30062, 30252, 30376, 33139, 34106, 34109, 34343))
MR <- data.frame(c(28000))

df_list_48E <- plot_bkpt(bkpt_mr, dfo=min_reg.df, b_window=1000, write_plot=TRUE, template=U048E)
df_MR <- plot_bkpt(MR, dfo=min_reg.df, b_window=6767, write_plot=FALSE, template=U048E)

#Plot minimum region in GK18 sequence
x1 <- 22200
x2 <- 34500
ggplot(df_MR[[1]], aes(x=Position)) + 
  geom_line(aes(y=P_melt, col="(P) melt region "), size=1) + geom_line(aes(y=P_Z, col="(P) Z-DNA "), size=1) + geom_line(aes(y=S_G4, col="(S) G-quadruplex "), size=1) + geom_line(aes(y=S_triplex, col="(S) triplex DNA"), size=1) + 
  scale_color_discrete(name="Probabilities (P) or Relative Scores (S):") +
  # ggtitle(paste("U048E positions", x1, "to", x2)) + 
  theme_light() + xlab("Position in U048-E") +
  # geom_line(aes(y=P_cruciform, col="P_cruciform")) +
  coord_cartesian(xlim=c((x1),(x2)), ylim=c(-0.2,1)) + 
  theme(axis.title.x=element_text(size=18), axis.title.y=element_blank(), axis.text.x=element_text(size=16), axis.text.y=element_text(size=14), legend.title = element_text(size=16, face="bold"), legend.text = element_text(size=18), legend.direction="horizontal", plot.title=element_text(size=18), legend.position = "top") + guides(colour = guide_legend(override.aes = list(size = 5))) +
  #+ geom_vline(xintercept=0) 
  geom_linerange(data=bkpt_mr_all, aes(x=bkpt_mr_all[,1], ymax=-0.1, ymin=0), size=1) +
  scale_x_continuous(expand = c(0,0)) + scale_y_continuous(expand = c(0,0))