#### Microbial (16S) Post-DADA2 Create Phyloseq and Clean ####
### Final analysis for RR2102 Biofilms4larvae MS

## Ladd, T.M., Selci, M., Davis, D.J., Cannon, O., Plowman, C.Q., 
## Schlegel,I.,Inaba, A., Mills, S.W., Vetriani, C., Mullineaux,L.S.,
## Arellano, S.M. (2024).Faunal colonists, including mussel settlers, 
## respond to microbial biofilms at deep-sea hydrothermal vents. Deep-Sea
## Research Part I

## This code is used to take the raw sequencing table from the DADA2 pipeline
## output, create a phyloseq object, and remove Eukaryotic, chloroplast, 
## and mitochondrial sequences. ASV and read loss is tracked throughout
## The output phyloseq object will be used for all downstream analysis.

##### Setup #####

#load packages
library(phyloseq)
library(Biostrings)

#clear workspace and set working directory
rm(list=ls())
#setwd() #CHANGE ME to desired directory

##### create phyloseq obj #####
#create phyloseq obj from merged seqtab (DADA2 separately run on each library)
#open seqtab_nochim and taxa .rds file
seqtab.nochim<-readRDS("seqtab.nochim_merged.rds")
taxa<-readRDS("taxa_species_merged.rds")

#open sample metadata
sampdat <- read.csv("RR2102_sandwich_metadata.csv", header = TRUE)
rownames(sampdat) <- sampdat$ID

#create phyloseq obj
prok_data_raw<-phyloseq(sample_data(sampdat),
                        otu_table(seqtab_nochim, taxa_are_rows = F),
                        tax_table(taxa))
prok_data_raw
sample_sums(prok_data_raw)

##### Clean up #####
# Clean up unwanted sequences from Eukarya, mitochrondria and chloroplast
prok_data_euk <- subset_taxa(prok_data_raw,  (Kingdom != "Eukaryota") | is.na(Kingdom))
prok_data_euk
sample_sums(prok_data_euk)
sample_sums(prok_data_euk)/sample_sums(prok_data_raw) # check percentage of reads after Eukaryota removal

euk_cleaning<-(1-(sample_sums(prok_data_euk)/sample_sums(prok_data_raw)))*100 # check percentage of reads after Eukaryota removal
euk_cleaning

prok_data_chl <- subset_taxa(prok_data_euk, (Order!="Chloroplast") | is.na(Order))
prok_data_chl
sample_sums(prok_data_chl)
(sample_sums(prok_data_chl)/sample_sums(prok_data_raw)) # check percentage of reads after Chloroplast removal

chlorp_cleaning<-(1-(sample_sums(prok_data_chl)/sample_sums(prok_data_raw)))*100
chlorp_cleaning

prok_data_mit <- subset_taxa(prok_data_chl, (Family!="Mitochondria") | is.na(Family))
prok_data_mit
sample_sums(prok_data_mit)
(sample_sums(prok_data_mit)/sample_sums(prok_data_raw)) # check percentage of reads after Chloroplast removal

mitoc_cleaning<-(1-(sample_sums(prok_data_mit)/sample_sums(prok_data_raw)))*100 # check percentage of reads after Chloroplast removal
mitoc_cleaning

##### track read and ASV loss #####
#Cleaning summary table
trackreads <- cbind(euk_cleaning, chlorp_cleaning,mitoc_cleaning)
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(trackreads) <- c("Eukaryotes", "Chloroplasts", "Mitochondria")
head(trackreads)

write.csv(trackreads, "pscleaning_trackreadsbysamp.csv")

ASVsLost <- c(ntaxa(prok_data_raw), ntaxa(prok_data_euk), 
              ntaxa(prok_data_chl),ntaxa(prok_data_mit))

ReadsLost <- c(sum(sample_sums(prok_data_raw)),sum(sample_sums(prok_data_euk)),
               sum(sample_sums(prok_data_chl)), sum(sample_sums(prok_data_mit)))

trackall <- data.frame("TrackASVs" = ASVsLost, "TrackReads" = ReadsLost)

rownames(trackall) <- c("Raw", "Eukaryotes", "Chloroplasts", "Mitochondria")

write.csv(trackall, "pscleaning_trackreadsASVstotal.csv")

##### Rename ASVs and save FASTA file #####

#final "cleaned" ps is prok_data_mit
ps <- prok_data_mit

#simplify seq names from raw seqs - add as Biostrings set in ps obj
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)

ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

names(dna) <- taxa_names(ps)

#save ASV#s w/ sequences as a fasta file - can be used to BLAST
writeXStringSet(dna, filepath = "allprokASVs4blast.fasta", format = "fasta")

##### edit sample data frame for use #####
#add library size to sample data frame
sample_data(ps)$libsize <- sample_sums(ps)

#clean up other sample data - make some factors w/ specific orders
sample_data(ps)$Zone <- factor(sample_data(ps)$Zone, levels = c("Alvinellid", "Riftia", "Mussel"))
sample_data(ps)$Pursed <- factor(sample_data(ps)$Pursed, levels = c("Y", "N"))
levels(sample_data(ps)$Pursed) <- c("Established", "Fresh")
sample_data(ps)$ID.pairs <- factor(sample_data(ps)$ID.pairs)
sample_data(ps)$ID <- factor(sample_data(ps)$ID, levels = c("PW", "PP10","PX", "PP1", "PM", "PP2",
                                                            "PP3", "PR", "PP5", "PO", "PP11", "PQ", "PP12"))

#save phyloseq object for all downstream analysis
saveRDS(ps, "prokps_final.rds")
