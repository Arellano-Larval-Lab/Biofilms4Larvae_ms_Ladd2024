#### Microbial (16S) DADA2 Pipeline ####
### Final analysis for RR2102 Biofilms4larvae MS

## Ladd, T.M., Selci, M., Davis, D.J., Cannon, O., Plowman, C.Q., 
## Schlegel,I.,Inaba, A., Mills, S.W., Vetriani, C., Mullineaux,L.S.,
## Arellano, S.M. (2024).Faunal colonists, including mussel settlers, 
## respond to microbial biofilms at deep-sea hydrothermal vents. Deep-Sea
## Research Part I

## This code is used to take the demultiplexed paired end sequencing reads 
## from the sequencing facility and run them through the DADA2 pipeline.
## https://benjjneb.github.io/dada2/tutorial.html
## This includes filtering and trimming of primers and low quality regions,
## learning errors, sample inference, merging paired end reads, creating the
## sequence table, removing chimeras, and assigning taxonomy
## Note samples were sequenced in 2 separate runs - they are processed 
## separately until the chimera removal step

##### Setup #####

#clear workspace and set working directory
rm(list = ls())
#setwd() #CHANGE ME to desired directory

#load packages
library(dada2)
library(magrittr)
library(dplyr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Seq run 1 ####
#Define paths to raw (F and R) fasta files
path1 <- "~/run1" # CHANGE ME to the directory containing the fastq files
list.files(path1) # Verify the file list

fnFs1 <- sort(list.files(path1, pattern="_R1_001.fastq", full.names = TRUE))
fnRs1 <- sort(list.files(path1, pattern="_R2_001.fastq", full.names = TRUE))

sample.names1 <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)
(sample.names1)

filtFs1 <- file.path(path1, "filtered_run1", paste0(sample.names1, "_F_filt.fastq.gz"))
filtRs1 <- file.path(path1, "filtered_run1", paste0(sample.names1, "_R_filt.fastq.gz"))

##### Assess quality #####
png("run1_QplotF.png", width=1600, height=1200)
plotQualityProfile(fnFs1) # Forward sequences
dev.off()

png("run1_QplotR.png", width=1600, height=1200)
plotQualityProfile(fnRs1) # reverse sequences
dev.off()

##### Filter and trim #####
out1 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(210,240), maxEE=c(3,3), maxN=0, truncQ=2, rm.phix=TRUE, trimLeft=c(19,20), compress=TRUE, multithread=25)
head(out1)
mean(out1[,2]/out1[,1])

##### Learn errors #####

#Alter loess arguments (weights and span) & enforce monotonicity 
#(from https://github.com/benjjneb/dada2/issues/1307 JacobRPrice)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF1 <- learnErrors(
  filtFs1,
  multithread = 25,
  errorEstimationFunction = loessErrfun_mod, MAX_CONSIST = 20,
  verbose = TRUE
)  


errR1 <- learnErrors(
  filtRs1,
  multithread = 25,
  errorEstimationFunction = loessErrfun_mod,MAX_CONSIST = 20,
  verbose = TRUE
)

#save plots
png("run1_errF.png", width=1600, height=1200)
plotErrors(errF1, nominalQ=TRUE)
dev.off()

png("run1_errR.png", width=1600, height=1200)
plotErrors(errR1, nominalQ=TRUE)
dev.off()


##### Sample inference #####
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)

derepRs1 <- derepFastq(filtRs1, verbose=TRUE)

names(derepFs1) <- sample.names1
names(derepRs1) <- sample.names1

dadaFs1 <- dada(derepFs1, err=errF1, multithread=20)
dadaRs1 <- dada(derepRs1, err=errR1, multithread=20)

##### Merge the paired reads #####
mergers1 <- mergePairs(dadaFs1, derepFs1, dadaRs1, derepRs1, 
                       trimOverhang = TRUE, verbose=TRUE)
head(mergers1[[1]])

##### Create sequence table #####
seqtab1 <- makeSequenceTable(mergers1)
dim(seqtab1) # Gives you info on the number of Samples and ASVs identified in run_1

#check seq lengths
table(nchar(getSequences(seqtab1)))

#subset based on seq length
seqtab1_f <- seqtab1[,nchar(colnames(seqtab1)) %in% seq(252,255)]

#confirm seq lengths
table(nchar(getSequences(seqtab1_f)))

dim(seqtab1_F) # Gives you info on the number of Samples and ASVs after tail sequence dropping

getN <- function(x) sum(getUniques(x))
track1 <- cbind(out1, sapply(dadaFs1, getN), sapply(dadaRs1, getN), sapply(mergers1, getN))
colnames(track1) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track1) <- sample.names1
head(track1) 

#save tracking file
write.csv(track1, "run1_asv_stats.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Seq run 2 ####
#Define paths to raw (F and R) fasta files
path2 <- "~/run2" # CHANGE ME to the directory containing the fastq files
list.files(path2) # Verify the file list

FnFs2 <- sort(list.files(path2, pattern="_R1_001.fastq", full.names = TRUE))
FnRs2 <- sort(list.files(path2, pattern="_R2_001.fastq", full.names = TRUE))

sample.names2 <- sapply(strsplit(basename(FnFs2), "_"), `[`, 1)
(sample.names2)

filtFs2 <- file.path(path2, "filtered_run2", paste0(sample.names2, "_F_filt.fastq.gz"))
filtRs2 <- file.path(path2, "filtered_run2", paste0(sample.names2, "_R_filt.fastq.gz"))

##### Assess quality #####
png("run2_QplotF.png", width=1600, height=1200)
plotQualityProfile(FnFs2) # Forward sequences
dev.off()

png("run1_QplotR.png", width=1600, height=1200)
plotQualityProfile(FnRs2) # reverse sequences
dev.off()

##### Filter and trim #####
out2 <- filterAndTrim(FnFs2, filtFs2, FnRs2, filtRs2, truncLen=c(210,240), maxEE=c(3,3), maxN=0, truncQ=2, rm.phix=TRUE, trimLeft=c(19,20), compress=TRUE, multithread=25)
head(out2)
mean(out2[,2]/out2[,1])

##### Learn errors #####

#Alter loess arguments (weights and span) & enforce monotonicity 
#(from https://github.com/benjjneb/dada2/issues/1307 JacobRPrice)

loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https://github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https://github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

# check what this looks like
errF2 <- learnErrors(
  filtFs2,
  multithread = 25,
  errorEstimationFunction = loessErrfun_mod, MAX_CONSIST = 20,
  verbose = TRUE
)  


errR2 <- learnErrors(
  filtRs2,
  multithread = 25,
  errorEstimationFunction = loessErrfun_mod,MAX_CONSIST = 20,
  verbose = TRUE
)

#save plots
png("run2_errF.png", width=1600, height=1200)
plotErrors(errF2, nominalQ=TRUE)
dev.off()

png("run2_errR.png", width=1600, height=1200)
plotErrors(errR2, nominalQ=TRUE)
dev.off()


##### Sample inference #####
derepFs2 <- derepFastq(filtFs2, verbose=TRUE)

derepRs2 <- derepFastq(filtRs2, verbose=TRUE)

names(derepFs2) <- sample.names2
names(derepRs2) <- sample.names2

dadaFs2 <- dada(derepFs2, err=errF2, multithread=20)
dadaRs2 <- dada(derepRs2, err=errR2, multithread=20)

##### Merge the paired reads #####
mergers2 <- mergePairs(dadaFs2, derepFs2, dadaRs2, derepRs2, trimOverhang = TRUE, verbose=TRUE)
head(mergers2[[1]])

##### Create sequence table #####
seqtab2 <- makeSequenceTable(mergers2)
dim(seqtab2) # Gives you info on the number of Samples and ASVs identified in run_1

#check seq lengths
table(nchar(getSequences(seqtab2)))

#subset based on seq length
seqtab2_f <- seqtab2[,nchar(colnames(seqtab2)) %in% seq(252,255)]

#confirm seq lengths
table(nchar(getSequences(seqtab2_f)))

dim(seqtab2_f) # Gives you info on the number of Samples and ASVs after tail sequence dropping

getN <- function(x) sum(getUniques(x))
track2 <- cbind(out2, getN(dadaFs2), getN(dadaRs2), getN(mergers2))
colnames(track2) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
rownames(track2) <- sample.names2
track2

#save tracking file
write.csv(track2, "run2_asv_stats.csv")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Chimera Removal ####
#combine seqtabs from different runs
seqtab <- mergeSequenceTables(seqtab1_f, seqtab2_f) # Merging of different runs

## From here on we refer to the seqtab.nochim object obtained from the 
## merging done between seqtab1_f (#run1) and seqtab2_f (#run2) 

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, minFoldParentOverAbundance = 10, 
                                    method="consensus", multithread=20)
dim(seqtab.nochim)
dim(seqtab)
sum(seqtab.nochim)/sum(seqtab) 
# Gives you the percentage of sequences recovered after chimera removal
## This is a good check points. Even if a lot of ASVs have been removed, the majority of reads should
## be retained. Usually >0.80 (aka 80%) are retained

# 77.6%

saveRDS(seqtab.nochim,"seqtab.nochim_merged.rds")

#add chimera removal to tracking files
track <- rbind(track1, track2)
track$nonchim <- rowSums(seqtab.nochim)

write.csv(track, "track_reads.csv")

#### Assign Taxonomy ####

#Point to where the silva database actually is
set.seed(10)
taxa <- assignTaxonomy(seqtab.nochim,"silva_nr99_v138.1_train_set.fa.gz",
                       multithread=20)
unname(taxa)

taxa.plus <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz", verbose=TRUE)

taxa.print <- taxa.plus # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
write.csv(taxa.print, "taxa_print_taxonomy.csv") # For inspection in bash or excel

saveRDS(taxa.plus,"taxa_species_merged.rds")
