library(dada2);packageVersion("dada2")
library(Biostrings); packageVersion("Biostrings")
library(ShortRead); packageVersion("ShortRead")
library(ggplot2); packageVersion("ggplot2")
library(reshape2); packageVersion("reshape2")
library(gridExtra); packageVersion("gridExtra")
library(phyloseq); packageVersion("phyloseq")
library(dplyr)
library(tidyr)

#update path to folder containaing the .fastqgz files
path2 <- "/Users/jrabasc/Downloads/glubs"
path.out <- "Figures/"
path.rds <- "RDS/"
fns2 <- list.files(path2, pattern="fastq.gz", full.names=TRUE)

#primer removal: change primers F27 and R1492 to whichever primers were used
F27 <- "AGRGTTYGATYMTGGCTCAG"
R1492 <- "RGYTACCTTGTTACGACTT"
rc <- dada2:::rc
nops2 <- file.path(path2, "noprimers", basename(fns2))
prim2 <- removePrimers(fns2, nops2, primer.fwd=F27, primer.rev=dada2:::rc(R1492), orient=TRUE)

#visualize read length distribution expect peak at ~1400bp for full length 16s
lens.fn <- lapply(nops2, function(fn) nchar(getSequences(fn)))
lens <- do.call(c, lens.fn)
hist(lens, 100)

#after viz set read cut offs
filts2 <- file.path(path2, "noprimers", "filtered", basename(fns2))
track2 <- filterAndTrim(nops2, filts2, minQ=3, minLen=1000, maxLen=1600, maxN=0, rm.phix=FALSE, maxEE=2)
track2

#dada2 derep and learn errors, plot errors: make sure it looks like a long check mark
drp2 <- derepFastq(filts2, verbose=TRUE)
err2 <- learnErrors(drp2, errorEstimationFunction=PacBioErrfun, BAND_SIZE=32, multithread=TRUE)
#saveRDS(err2, file.path("/Users/jrabasc/Desktop/github/DMC_pacbio_read_fate_analysis/assets/error_model.rds"))
plotErrors(err2)

#dada2 denoise and bimera removal
dd2 <- dada(drp2, err=err2, BAND_SIZE=32, multithread=TRUE)
#saveRDS(dd2, file.path("/Users/jrabasc/Desktop/github/DMC_pacbio_read_fate_analysis/assets/dada_results.rds"))

read_ret_table<-as.data.frame(cbind(ccs=prim2[,1], primers=prim2[,2], filtered=track2[,2], denoised=sapply(dd2, function(x) sum(x$denoised))))
read_ret_table$read_ret_per<- read_ret_table$denoised / read_ret_table$ccs
#write.csv(read_ret_table, "/Users/jrabasc/Desktop/github/DMC_pacbio_read_fate_analysis/assets/read_retention_table.csv", row.names = FALSE)


graph_df<-read_ret_table
graph_df$rem_primer_removal<-read_ret_table$ccs - read_ret_table$primers
graph_df$rem_filtering<-read_ret_table$primers - read_ret_table$filtered
graph_df$rem_denoised<-read_ret_table$filtered - read_ret_table$denoised
graph_df$reads_retained<-read_ret_table$denoised
graph_df <- graph_df[, !colnames(graph_df) %in% c("ccs", "primers", "filtered", "denoised","read_ret_per")]
df_long <- graph_df %>%
  pivot_longer(
    cols = colnames(graph_df),
    names_to = "Read_Fate",
    values_to = "Reads"
  )
sample_names<-rownames(graph_df)
df_long$sample_names <- unlist(lapply(sample_names, function(x) rep(x, 4)))
ggplot(df_long, aes(x = sample_names, y = Reads, fill = Read_Fate)) +
  geom_bar(stat = "identity") +
  labs(title = "Read Fate", x = "Samples", y = "Read Numbers") +
  theme(axis.text.x = element_blank())
st2 <- makeSequenceTable(dd2); dim(st2)
bim2 <- isBimeraDenovo(st2, minFoldParentOverAbundance=3.5, multithread=TRUE)
sum(st2[,bim2])/sum(st2)






