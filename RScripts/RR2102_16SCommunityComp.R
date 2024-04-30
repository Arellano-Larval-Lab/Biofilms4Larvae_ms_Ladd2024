#### Microbial (16S) Analysis ####
### Final analysis for RR2102 Biofilms4larvae MS

## Ladd, T.M., Selci, M., Davis, D.J., Cannon, O., Plowman, C.Q., 
## Schlegel,I.,Inaba, A., Mills, S.W., Vetriani, C., Mullineaux,L.S.,
## Arellano, S.M. (2024).Faunal colonists, including mussel settlers, 
## respond to microbial biofilms at deep-sea hydrothermal vents. Deep 
## Sea Research Part I

## This code is used to analyze microbial 16S data on established and 
## fresh biofilm sandwiches across biogenic zones. It also includes 
## proportionality analysis for microbial-faunal interactions. It includes
## statistical analyses and figure creation.

##### Setup #####

#load packages
library(phyloseq)
library(ggplot2)
library(vegan)
library(car)
library(pairwiseAdonis)
library(dplyr)
library(ANCOMBC)
library(tibble)
library(Biostrings)
library(DECIPHER)
library(propr)
library(zCompositions)
library(reshape2)
library(circlize)
library(seecolor)

#clear workspace and set working directory
rm(list=ls())
#setwd() #CHANGE ME to desired directory

#open phyloseq obj
ps <- readRDS("prokps_final.rds")
#create relative abundance ps object
ps.ra <- transform_sample_counts(ps, function(otu) otu/sum(otu))

#extract sample data frame from ps obj (for easier manipulation)
sampdat <- data.frame(sample_data(ps))


#### Rarefaction curves ####
#plot rarefaction curves, color by biofilm age, line type by zone
Cmin = min(sample_sums(ps)) #minimum library size
otutab <- otu_table(ps)
class(otutab) <- "matrix"

#save rarefaction curve
tiff("Rarecurve_final.tiff", width = 2031, height = 1562, res = 300)

rarecurve(otutab, step=100, 
          col=c("#7e03a8", "#7e03a8", "#f89540", "#f89540",
                "#f89540","#f89540", "#f89540", "#f89540", "#7e03a8",
                "#7e03a8","#7e03a8","#7e03a8", "#f89540"), 
          lty=c(2,1,2,3,1,1,1,1,1,1,3,2,2),
          lwd=2, ylab="ASV Richness", label=F)

abline(v = Cmin)
legend("topright", legend=c("Established", "Fresh"),
       col=c("#7e03a8", "#f89540"), lty=1, cex=0.8, bty = "n")
legend("bottomright", legend=c("Alvinellid", "Riftia", "Mussel"),
       col="black", lty=c(3,2,1), cex=0.8, bty = "n")

dev.off()

#### Alpha diversity ####
#plot alpha diversity measures (Richness/Shannon) by zone and biofilm age
tiff("Alphadiv_byZoneandAge_final.tiff", width = 2187, height = 1719, res = 300)

p1 <- plot_richness(ps, x = "Zone", measures = c("Observed", "Shannon"), color = "Pursed") + 
  geom_boxplot(position = position_dodge(width = 0.75)) +
  scale_color_manual(values=c("Established" = "#7e03a8", "Fresh" = "#f89540"), name = "Biofilm") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA),
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black", angle = 60, hjust = 1),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))

p1$layers <- p1$layers[-1]
p1 + geom_point(position = position_dodge(width = 0.75))

dev.off()

#compare diversity measures statistically
divps <- estimate_richness(ps, measures=c("Observed", "Shannon"))

#add data to divps obj
divps$Pursed <- sampdat$Pursed
divps$ID.pairs <- sampdat$ID.pairs
divps$Zone <- sampdat$Zone

#choose div measure to test: Observed, Shannon
res.aov <- aov(Shannon ~ Zone*Pursed, data = divps)
plot(res.aov)
shapiro.test(res.aov$residuals) #normality
leveneTest(Shannon ~ Zone*Pursed, data = divps) #heterogeneity of variance
summary(res.aov)
#no interaction term
res.aov <- aov(Shannon~ Zone+Pursed, data = divps)
plot(res.aov)
shapiro.test(res.aov$residuals) #normality
leveneTest(Shannon ~ Zone, data = divps) #heterogeneity of variance
leveneTest(Shannon ~ Pursed, data = divps) #heterogeneity of variance
summary(res.aov)
#remove Alvinellid zone due to only 1 pair
res.aov <- aov(Shannon ~ Zone*Pursed, data = divps[divps$Zone %in% c("Riftia", "Mussel"),])
plot(res.aov)
shapiro.test(res.aov$residuals) #normality
leveneTest(Shannon ~ Zone*Pursed, data = divps[divps$Zone %in% c("Riftia", "Mussel"),]) #heterogeneity of variance
summary(res.aov)
#no interaction
res.aov <- aov(Shannon ~ Zone+Pursed, data = divps[divps$Zone %in% c("Riftia", "Mussel"),])
plot(res.aov)
shapiro.test(res.aov$residuals)
leveneTest(Shannon ~ Zone, data = divps[divps$Zone %in% c("Riftia", "Mussel"),])
leveneTest(Shannon ~ Pursed, data = divps[divps$Zone %in% c("Riftia", "Mussel"),])
summary(res.aov)

#Try paired test
#separate TicaSands into two dfs
divEstab <- subset(divps, divps$Pursed == "Established")
divEstab <- divEstab[order(divEstab$ID.pairs),]

divFresh <- subset(divps, divps$Pursed == "Fresh")
divFresh <- divFresh[order(divFresh$ID.pairs),]

#not all the same ID pairs...
sharepairs <- intersect(divFresh$ID.pairs, divEstab$ID.pairs)

divEstab <- divEstab[divEstab$ID.pairs %in% sharepairs,]
divFresh <- divFresh[divFresh$ID.pairs %in% sharepairs,]

PairDiff <- data.frame("ID.pairs"=divEstab[,4])
PairDiff$Obsdiff <- divEstab$Observed - divFresh$Observed
PairDiff$Shandiff <- divEstab$Shannon - divFresh$Shannon

(shapiro <- shapiro.test(PairDiff$Obsdiff)) #normal
(t<- t.test(PairDiff$Obsdiff, mu = 0, alternative = "two.sided"))#p=0.07303

(shapiro <- shapiro.test(PairDiff$Shandiff)) #not normal
#(t<- t.test(PairDiff$Shandiff, mu = 0, alternative = "two.sided"))
(wilcox <- exactRankTests::wilcox.exact(PairDiff$Shandiff, mu = 0, alternative = "two.sided"))#p=0.4375

#### Beta diversity - PCoA Bray ####
ord.PCoA.bray <- ordinate(ps.ra, method = "PCoA", distance = "bray")
#look at plot- get axes % variation explained: Axis.1: 28.2%, Axis.2: 19.9%
plot_ordination(ps.ra, ord.PCoA.bray)

ord.PCoA.bray$values$Eigenvalues

#PCoA get "scores" first 2 axes
PCoAaxes <- ord.PCoA.bray$vectors[,1:2]

#df with PCoA axes and ancil data
scoresdf <- data.frame(PCoAaxes)

scoresdf$ID <-sampdat$ID #sandwich IDs
scoresdf$Zone <- sampdat$Zone #biogenic zone
scoresdf$Pursed <- sampdat$Pursed #biofilm age
scoresdf$ID.pairs <-sampdat$ID.pairs # pairs of sands


#hulls by zone
Riftia <- scoresdf[scoresdf$Zone == "Riftia", ][chull(scoresdf[scoresdf$Zone == "Riftia", c("Axis.1", "Axis.2")]), ]  # hull values for Riftia Zone
Mussel <- scoresdf[scoresdf$Zone == "Mussel", ][chull(scoresdf[scoresdf$Zone == "Mussel", c("Axis.1", "Axis.2")]), ]  # hull values for Mussel Zone
Alvinellid <- scoresdf[scoresdf$Zone == "Alvinellid", ][chull(scoresdf[scoresdf$Zone == "Alvinellid", c("Axis.1", "Axis.2")]), ]  # hull values for Alvinellid Zone

hull.data <- rbind(Riftia, Mussel, Alvinellid)

p2 <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=Axis.1,y=Axis.2,fill=Zone,group=Zone),alpha=0.30) + # add the convex hulls
  geom_point(data=scoresdf,aes(x=Axis.1,y=Axis.2,shape=Pursed,fill=Zone),colour = "black", size=4) + # add the point markers
  geom_text(data=scoresdf,aes(x=Axis.1,y=Axis.2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_manual(values=c("Alvinellid" = "#fde725", "Riftia" = "#21918c", "Mussel" = "#440154"), name = "Zone") +
  scale_shape_manual(values=c("Established" = 21, "Fresh" = 25), name = "Biofilm", labels = c("Established", "Fresh")) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2), lim = c(-0.41, 0.25)) +
  #scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4), lim = c(-0.4, 0.5)) +
  labs(x = "Axis 1 [28.2%]", y = "Axis 2 [19.9%]") +
  coord_equal() +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1)), shape = guide_legend(order = 2)) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p2

tiff("PCoA_Ticasands_16SRA_byZone_final.tiff", width = 1565, height = 1370, res = 300)
p2
dev.off()

#hulls by filmed/unfilmed
Filmed <- scoresdf[scoresdf$Pursed == "Established", ][chull(scoresdf[scoresdf$Pursed == "Established", c("Axis.1", "Axis.2")]), ]  # hull values for filmed samples
Unfilmed <- scoresdf[scoresdf$Pursed == "Fresh", ][chull(scoresdf[scoresdf$Pursed == "Fresh", c("Axis.1", "Axis.2")]), ]  # hull values for unfilmed samples
hull.data <- rbind(Filmed, Unfilmed)

p3 <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=Axis.1,y=Axis.2,fill=Pursed,group=Pursed),alpha=0.30) + # add the convex hulls
  geom_point(data=scoresdf,aes(x=Axis.1,y=Axis.2,shape=Zone,fill=Pursed),colour = "black", size=4) + # add the point markers
  geom_text(data=scoresdf,aes(x=Axis.1,y=Axis.2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_manual(values=c("Established" = "#7e03a8", "Fresh" = "#f89540"), name = "Biofilm", labels = c("Established", "Fresh")) +
  scale_shape_manual(values=c("Alvinellid" = 22,"Riftia" = 23, "Mussel" = 24), name = "Zone") +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2), lim = c(-0.41, 0.25)) +
  #scale_x_continuous(breaks = c(-0.4, -0.2, 0, 0.2, 0.4), lim = c(-0.4, 0.5)) +
  labs(x = "Axis 1 [28.2%]", y = "Axis 2 [19.9%]") +
  coord_equal() +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1)), shape = guide_legend(order = 2)) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p3

tiff("PCoA_Ticasands_16SRA_bybiofilm_final.tiff", width = 1565, height = 1370, res = 300)
p3
dev.off()

##check if axes correlate with temp?
scoresdf$Exp.start.temp <- sampdat$Exp_start_temp.C
scoresdf$Recovery.temp <- sampdat$Recovery_temp.C


PCoA.fit <- envfit(PCoAaxes, scoresdf[,7:8], permutations = 999, na.rm = TRUE) 
PCoA.fit

en.data.scores = as.data.frame(scores(PCoA.fit, "vectors"))*0.5


p4 <- ggplot() + 
  geom_point(data=scoresdf,aes(x=Axis.1,y=Axis.2,shape=Zone,fill=Recovery.temp),colour = "black", size=4) + # add the point markers
  geom_text(data=scoresdf,aes(x=Axis.1,y=Axis.2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_viridis_c(option = "plasma", name = expression(paste("Temperature (",degree*C,")"))) +
  scale_shape_manual(values=c("Alvinellid" = 22,"Riftia" = 23, "Mussel" = 24), name = "Zone") +
  geom_segment(aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2), 
               data = en.data.scores, linewidth = 1, colour = "grey30", 
               arrow = arrow(length = unit(0.5,"cm"))) +
  geom_text(data = en.data.scores, aes(x = Axis.1 + 0.02*sign(Axis.1), y = Axis.2 + 0.02*sign(Axis.2)), 
            colour = "black", size = 3, hjust = 0.9,
            label = c("Start Temp", "Recovery Temp")) +
  scale_y_continuous(breaks = c(-0.4, -0.2, 0, 0.2), lim = c(-0.41, 0.25)) +
  labs(x = "Axis 1 [28.2%]", y = "Axis 2 [19.9%]") +
  coord_equal() +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p4

tiff("PCoA_Ticasands_16SRA_byZonewTemp_final.tiff", width = 1565, height = 1370, res = 300)
p4
dev.off()

##### Adonis - differences in community composition #####

#create dist obj
bray_dist <- phyloseq::distance(ps.ra, method="bray")

#just check on sig effects of libsize
(res_adonis <- adonis2(bray_dist ~ libsize, data = sampdat)) #non-sig

#full model
(res_adonis <- adonis2(bray_dist ~ Zone*Pursed, data = sampdat)) #marginally sig effects of zone and pursed, no interaction effects

#save stats table
data.frame(name = rownames(res_adonis), Df = res_adonis$Df,
           SumsOfSqs = res_adonis$SumOfSqs, R2 = res_adonis$R2, 
           F.model = res_adonis$F, 
           P = res_adonis$`Pr(>F)`)  %>% write.csv(file = "RR2102_PERMANOVA_braydistRA_full.csv")

#remove interaction
(res_adonis <- adonis2(bray_dist ~ Zone + Pursed, data = sampdat)) #sig effects of zone and pursed

data.frame(name = rownames(res_adonis), Df = res_adonis$Df,
           SumsOfSqs = res_adonis$SumOfSqs, R2 = res_adonis$R2, 
           F.model = res_adonis$F, 
           P = res_adonis$`Pr(>F)`) %>% write.csv(file = "RR2102_PERMANOVA_braydistRA_noint.csv")

#pairwise test to look at differences between zones
pairwise.adonis2(bray_dist ~ Zone, data = sampdat)

#check for differences in multivariate dispersion based on biofilm or zone (like homogeneity of variances)
beta <- betadisper(bray_dist, sampdat$Pursed)
permutest(beta) #no sig heterogeneity of dispersion

beta <- betadisper(bray_dist, sampdat$Zone)
permutest(beta) #no sig heterogeneity of dispersion

##### Top 10 families barchart #####

#make a plot that is ordered by relative abundance
#create agglom ps RA obj at Family and Phylum level
ps.family <- tax_glom(ps.ra, taxrank = "Family", NArm = F)
ps.phylum <- tax_glom(ps.ra, taxrank = "Phylum", NArm = F)

#see the top 10 phylums
topnphy <- names(sort(taxa_sums(ps.phylum), TRUE)[1:10])
tax_table(ps.phylum)[topnphy,]

#top 10 families
topnfam <- names(sort(taxa_sums(ps.family), TRUE)[1:10])
tax_table(ps.family)[topnfam,] #top 10 families are all within the top 4 phylums

#create new dataframe that stacks abundance by families
df <- psmelt(ps.family)

#create new groupings for phylum and families
phylums <- c("Campylobacterota","Proteobacteria","Bacteroidota", "Verrucomicrobiota")
#top four phylums %RA
sample_sums(subset_taxa(ps.phylum, Phylum %in% phylums))


df$Phylum[!df$Phylum %in% phylums] <- "Other"
df$Family[!df$Phylum %in% phylums] <- "Other"

df$Family[df$Phylum=="Campylobacterota" & 
            !df$Family %in% c("Sulfurovaceae","Sulfurimonadaceae",
                              "Arcobacteraceae")] <- "Other Campylobacterota"

df$Family[df$Phylum=="Proteobacteria" & df$Class== "Gammaproteobacteria" &
            !df$Family %in% c("Thiotrichaceae","Colwelliaceae")] <- "Other Gammaproteobacteria"

df$Family[df$Phylum=="Proteobacteria" & 
            !df$Class== "Gammaproteobacteria"] <- "Other Proteobacteria"

df$Family[df$Phylum=="Bacteroidota" &
            !df$Family %in% c("Crocinitomicaceae","Flavobacteriaceae","Marinifilaceae",
                              "Cryomorphaceae")] <- "Other Bacteroidota"

df$Family[df$Phylum=="Verrucomicrobiota" & 
            !df$Family %in% "Rubritaleaceae"] <- "Other Verrucomicrobiota"

#edit df to choose specific variables and make phylum and families factors
df2 <- dplyr::select(df, Sample, ID.pairs, Abundance, Zone, Pursed, Phylum, Family) %>%
  mutate(Phylum=factor(Phylum, levels=c(phylums, "Other")))

#make sample order correct
df2$Sample <- factor(df2$Sample , levels = levels(sample_data(ps)$ID))
df2$ID.pairs <- factor(df2$ID.pairs , levels = levels(sample_data(ps)$ID.pairs))
#order Family levels
df2$Family <- factor(df2$Family,
                     levels = c("Sulfurovaceae","Sulfurimonadaceae","Arcobacteraceae",
                                "Other Campylobacterota", "Thiotrichaceae",
                                "Colwelliaceae","Other Gammaproteobacteria", "Other Proteobacteria","Crocinitomicaceae",
                                "Flavobacteriaceae","Marinifilaceae","Cryomorphaceae", "Other Bacteroidota",
                                "Rubritaleaceae","Other Verrucomicrobiota", "Other"))

#summary of family RA
df2 %>% 
  group_by(Family) %>% 
  summarize(mean = mean(Abundance),
            min = min(Abundance),
            max = max(Abundance))



#color palette function to create a color palette based on nested groups
#https://copyprogramming.com/howto/stacked-bar-plot-r-with-age-groups
ColourPalleteMulti <- function(df, group, subgroup){
  
  # Find how many colour categories to create and the number of colours in each
  categories <- aggregate(as.formula(paste(subgroup, group, sep="~" )), df, function(x) length(unique(x)))
  category.start <- (scales::hue_pal(l = 80)(nrow(categories))) # Set the top of the colour pallete
  category.end  <- (scales::hue_pal(l = 30)(nrow(categories))) # set the bottom
  
  # Build Colour pallette
  colours <- unlist(lapply(1:nrow(categories),
                           function(i){
                             colorRampPalette(colors = c(category.start[i], category.end[i]))(categories[i,2])}))
  return(colours)
}

#get color palette for the top 10 families within 4 phylums
colours <-ColourPalleteMulti(df2, "Phylum", "Family")
names(colours) <- levels(df2$Family)
#change "Other" to grey
colours[16] <- "grey90"

p5 <- ggplot(df2, aes(fill=Family, y=Abundance*100, x=Sample)) + 
  geom_bar(position=position_stack(reverse = TRUE), stat="identity") +
  scale_fill_manual(name = "Family", values = colours) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 100), name = "Relative abundance (%)") +
  facet_grid(cols = vars(Zone), scales="free_x", space = "free_x") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_blank(),
        axis.text.x = element_text(size=12, colour="black", angle = 45,  hjust = 1),  
        axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14), 
        axis.title.x = element_text(size=14, colour="black"), 
        axis.title.y = element_text(size=14, colour="black"))
p5

tiff("Barchart_Ticasands_16S_top10Families_final.tiff", width = 2813, height = 1875, res = 300)
p5
dev.off()

##### Differential abundance analysis #####

#filter based on prevalence and max relative abundance
#Prevalence of a given ASV
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
#maximum relative abundance of each ASV
relabun = apply(otu_table(ps.ra), 
                MARGIN = ifelse(taxa_are_rows(ps.ra), yes = 1, no = 2),
                max)

#create prev and relabun data table w/ taxonomy
prevdf = data.frame(Prevalence = prevdf,
                    MaxRelAbundance = relabun,
                    tax_table(ps))


# Execute rel abun and prev filter, using `prune_taxa()` function
# Define abundance and prevelence threshold for taxa 
# must be >= 0.1% in at least 1 sample, and be in 3 or more samples
relAbundanceThreshold <- 0.001
prevalenceThreshold <- 3

keepTaxa <- rownames(prevdf)[(prevdf$MaxRelAbundance >= relAbundanceThreshold &
                                prevdf$Prevalence >= prevalenceThreshold)]

ps.f1 <- prune_taxa(keepTaxa, ps) #ps filt obj


ancom.out <- ancombc(phyloseq = ps.f1,formula = "Zone + Pursed",
                     p_adj_method = "BH", zero_cut = 1,
                     lib_cut = 0, group = NULL,
                     struc_zero = FALSE, neg_lb = FALSE,
                     tol = 1e-05, max_iter = 100,
                     conserve = FALSE, alpha = 0.05, global = FALSE)

# Extract results and add taxa designations to table
ancom.res <- data.frame(Representative_ASV=rownames(ancom.out$res$beta), Beta=ancom.out$res$beta[,3], 
                        SE=ancom.out$res$se[,3], P=ancom.out$res$p_val[,3], FDR_BH=ancom.out$res$q_val[,3],
                        DA=ancom.out$res$diff_abn[,3])
results <- inner_join(ancom.res, rownames_to_column(data.frame(tax_table(ps.f1)), "Representative_ASV"), by = "Representative_ASV")
results <- subset(results, results$DA == "TRUE")
results$Beta <- -results$Beta

write.csv(results, file = "RR2102_16S_ANCOMBC.csv")

#add new "phylum/class" column that can be changed to a factor
results$phyclass <- results$Phylum
unique(results$phyclass)

#split up Proteobacteria to show Gamma vs Alpha vs Zeta
for (i in 1:nrow(results)){
  if (results$Phylum[i] %in% "Proteobacteria" & results$Class[i] %in% "Gammaproteobacteria") {
    results$phyclass[i] <- results$Class[i]
  } else if (results$Phylum[i] %in% "Proteobacteria" & results$Class[i] %in% "Alphaproteobacteria") {
    results$phyclass[i] <- results$Class[i]
  } else if (results$Phylum[i] %in% "Proteobacteria" & results$Class[i] %in% "Zetaproteobacteria") {
    results$phyclass[i] <- results$Class[i]
  } else if (results$Phylum[i] %in% "Proteobacteria") {
    results$phyclass[i] <- "Other Proteobacteria"
  } else {
    results$phyclass[i] <- results$Phylum[i]
  } 
} 

#make phyclass a factor
unique(results$phyclass)
results$phyclass <- factor(results$phyclass, levels = c("Campylobacterota","Gammaproteobacteria",
                                                        "Alphaproteobacteria","Zetaproteobacteria",
                                                        "Bacteroidota", "Verrucomicrobiota",
                                                        "Desulfobacterota","Patescibacteria",
                                                        "Planctomycetota","Actinobacteriota"))


taxcolors = c("Campylobacterota" = "#C8564B","Gammaproteobacteria" = "#CCCE00",
              "Alphaproteobacteria"="#8B8D00","Zetaproteobacteria" = "#4B4D00",
              "Bacteroidota" ="#00BC79", "Verrucomicrobiota"= "#005FA7",
              "Desulfobacterota"= "violetred2", "Planctomycetota" = "purple4",
              "Patescibacteria" = "forestgreen","Actinobacteriota" = "darkred")

#plot log-linear coefficients (DA)
p6 <- ggplot(results, aes(x=Beta, y=phyclass, fill = phyclass)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(shape = 21, position = position_jitterdodge(jitter.width = 0.8, jitter.height = 0, seed = 2), size = 7,
             stroke = 1, alpha = 0.5) +
  scale_fill_manual(values = taxcolors, guide = "none") +
  scale_y_discrete(limits = rev(levels(results$phyclass))) +
  labs(y = "Taxonomic Phylum or Class", x = expression("Coefficient"), title = " ") +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        axis.title.x = element_text(size = 14), axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(size = 14, colour = "black")) 
p6

tiff("ANCOMBC_Ticasands_ASVprev3_final.tiff", width = 1400, height = 1600, res = 300)
p6
dev.off()

#DA trends by phylum/class
results[results$Class %in% "Gammaproteobacteria",]
results[results$Class %in% "Alphaproteobacteria",]
results[results$Phylum %in% "Actinobacteriota",]

##### Test differences in RA of top families across biofilm age #####

#top 10 families
topnfam <- names(sort(taxa_sums(ps.family), TRUE)[1:10])
fams <- as.vector(tax_table(ps.family)[topnfam, "Family"]) #top 10 families are all within the top 4 phylums


#subset df2 to only include top 10 fams, not others
famdf <- subset(df2, df2$Family %in% fams)

famdf$Family <- factor(famdf$Family, levels = fams)

#paired test
#separate TicaSands into two dfs
famEstab <- subset(famdf, famdf$Pursed == "Established")
famEstab <- famEstab[order(famEstab$Family, famEstab$ID.pairs),]

famFresh <- subset(famdf, famdf$Pursed == "Fresh")
famFresh <- famFresh[order(famFresh$Family, famFresh$ID.pairs),]

#not all the same ID pairs...
sharepairs <- intersect(famFresh$ID.pairs, famEstab$ID.pairs)

famEstab <- famEstab[famEstab$ID.pairs %in% sharepairs,]
famFresh <- famFresh[famFresh$ID.pairs %in% sharepairs,]

PairDiff <- famEstab[,c(2,4,6,7)]
PairDiff$RAdiff <- famEstab$Abundance - famFresh$Abundance

PairDiff$Family <- factor(PairDiff$Family,
                          levels = c("Sulfurovaceae","Sulfurimonadaceae","Arcobacteraceae",
                                     "Thiotrichaceae","Colwelliaceae","Crocinitomicaceae",
                                     "Flavobacteriaceae","Marinifilaceae","Cryomorphaceae", 
                                     "Rubritaleaceae"))


p7 <- ggplot(PairDiff, 
            aes(x=Family, y=RAdiff, fill = Family)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.shape = 21, width = 0.5) +
  stat_summary(fun=mean, geom="point", shape=4, size=1.5, color="black", stroke = 1) +
  scale_fill_manual(values = colours, guide = "none") +
  scale_x_discrete(limits = rev) +
  labs(x = " ", y = expression("RA difference (Established - Fresh)")) +
  coord_flip() +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))
p7

tiff("Boxplot_Ticasands_PairedFamilyRADiff_final.tiff", width = 1400, height = 1600, res = 300)
p7
dev.off()

PairDiffwide <- reshape(PairDiff, idvar = c("ID.pairs", "Zone"),
                        timevar = "Family", direction = "wide", drop = "Phylum")

TicaPairStats <- data.frame("Family" = fams, "shapiro_p" = NA, 
                            "t_p" = NA, "wilcox_p" = NA, "shapiro_p2" = NA,
                            "t_p2" = NA, "wilcox_p2" = NA)

for (i in 1:nrow(TicaPairStats)) {
  shapiro <- shapiro.test(PairDiffwide[,i+2])
  TicaPairStats$shapiro_p[i] <- shapiro$p.value
  t<- t.test(PairDiffwide[,i+2], mu = 0, alternative = "two.sided")
  TicaPairStats$t_p[i] <- t$p.value
  wilcox <- exactRankTests::wilcox.exact(PairDiffwide[,i+2], mu = 0, alternative = "two.sided")
  TicaPairStats$wilcox_p[i] <- wilcox$p.value
}

#excluding alvinellid zone?
for (i in 1:nrow(TicaPairStats)) {
  shapiro <- shapiro.test(PairDiffwide[-1,i+2])
  TicaPairStats$shapiro_p2[i] <- shapiro$p.value
  t<- t.test(PairDiffwide[-1,i+2], mu = 0, alternative = "two.sided")
  TicaPairStats$t_p2[i] <- t$p.value
  wilcox <- exactRankTests::wilcox.exact(PairDiffwide[-1,i+2], mu = 0, alternative = "two.sided")
  TicaPairStats$wilcox_p2[i] <- wilcox$p.value
}

write.csv(TicaPairStats, "RR2102_TicaPair16S_top10FamilyRA_final.csv")

##### compare ASVs to symbiont seqs #####
S <- readDNAStringSet("allprokASVs4blast.fasta") #fasta for our dataset

#Bathymodiolus thermophilus
#https://www.ncbi.nlm.nih.gov/nuccore/M99445.1
#Distel et al. 1988: https://doi.org/10.1128%2Fjb.170.6.2506-2510.1988
R <- readDNAStringSet("B_thermophilus_symbiont_16Ssequence.fasta")#dowloaded symbiont seq

pa <-  pairwiseAlignment(S, R, type = "overlap")
summary(pa)
pa #overlap starts at position 521
matchper <- data.frame("ASV" = names(S), "perID" = nmatch(pa)/width(S), "score" = score(pa))

#ASV w/ best match percent
matchper[matchper$perID == max(matchper$perID),]

#top 5
head(matchper[order(matchper$perID, decreasing = TRUE), ], n=5)
#   ASV     perID    score
# ASV284  0.9960474 493.5022
# ASV2832 0.9881423 477.7401
# ASV3986 0.9841897 469.8591
# ASV3279 0.9683794 438.3350
# ASV1223 0.9565217 414.6920

#highest score
matchper[matchper$score == max(matchper$score),]

head(matchper[order(matchper$score, decreasing = TRUE), ], n=5)
#   ASV     perID    score
# ASV284  0.9960474 493.5022
# ASV2832 0.9881423 477.7401
# ASV3986 0.9841897 469.8591
# ASV3279 0.9683794 438.3350
# ASV1223 0.9565217 414.6920

tiff("Boxplot_Ticasands_BthermsymRA_final.tiff", width = 1875, height = 1562, res = 300)
boxplot(otu_table(ps.ra)[,"ASV284"]*100 ~ sample_data(ps.ra)$Pursed, 
        ylab =" ", xlab = "Biofilm Age", las = 1, 
        col = c("#7e03a8", "#f89540"), par(mar = c(3, 6, 2, 2) + 0.1))
title(ylab = "B. thermophilus symbiont RA (%)", line = 3); 
dev.off()

#check out RAs for each sandwich
otu_table(ps.ra)[,"ASV284"]

#check out taxonomy
tax_table(ps.ra)["ASV284",]
tax_table(subset_taxa(ps.ra, Genus == "SUP05 cluster"))
#ASV284, ASV564, ASV2759, ASV3279, ASV3986

#look at sequences that match to B thermophilus symb
S2 <- S[c("ASV284", "ASV2832", "ASV3986", "ASV3279", "ASV1223")]
#SUP05 cluster
S3 <- S[c("ASV84", "ASV564", "ASV2759", "ASV3279", "ASV3986")]
#subset R to just region that aligns with seqs
Rsub <- subseq(R, 521, 773)

#top 1 match from pairwise alignment
BrowseSeqs(c(S["ASV284"], Rsub), "B.therm.bestmatch.html", openURL = TRUE)#save sequences to file- color coded to see trends
#top 5 matches from pairwise alignment
BrowseSeqs(c(S2, Rsub), "B.therm.matches.html", openURL = TRUE) #save sequences to file- color coded to see trends
#all "SUP05 cluster" seqs
BrowseSeqs(c(S3, Rsub), "B.therm.SUP05cluster.html", openURL = TRUE)#save sequences to file- color coded to see trends

#Riftia pachyptila
#https://www.ncbi.nlm.nih.gov/nuccore/175863
#Distel et al. 1988: https://doi.org/10.1128%2Fjb.170.6.2506-2510.1988
R <- readDNAStringSet("/Riftia_pachyptila_trophosome_symbiont_16S.fasta")#dowloaded symbiont seq

pa <-  pairwiseAlignment(S, R, type = "overlap")
summary(pa)
pa #overlap starts at position 424
matchper <- data.frame("ASV" = names(S), "perID" = nmatch(pa)/width(S), "score" = score(pa))

#ASV w/ best match percent
matchper[matchper$perID == max(matchper$perID),]#ASV3120

#top 5
head(matchper[order(matchper$perID, decreasing = TRUE), ], n=5)
#    ASV     perID    score
#   ASV3120 0.9723320 446.2161
#   ASV2040 0.9683794 438.3350
#   ASV1345 0.9644269 430.4539
#   ASV3704 0.9644269 430.4539
#   ASV2738 0.9604743 422.5729

#highest score
matchper[matchper$score == max(matchper$score),]

head(matchper[order(matchper$score, decreasing = TRUE), ], n=5)
#    ASV     perID    score
#   ASV3120 0.9723320 446.2161
#   ASV2040 0.9683794 438.3350
#   ASV3704 0.9644269 430.4539
#   ASV1345 0.9644269 430.4539
#   ASV2738 0.9604743 422.5729

tiff("figs/Boxplot_Ticasands_RiftiasymRA_final.tiff", width = 1875, height = 1562, res = 300)
boxplot(otu_table(ps.ra)[,"ASV3120"]*100 ~ sample_data(ps.ra)$Pursed, 
        ylab =" ", xlab = "Biofilm Age", las = 1, 
        col = c("#7e03a8", "#f89540"), par(mar = c(3, 6, 2, 2) + 0.1))
title(ylab = "R. pachyptila symbiont RA (%)", line = 4); 
dev.off()

#check out RAs for each sandwich
otu_table(ps.ra)[,"ASV3120"] #only in one sample, ~0.007%

#check out taxonomy
tax_table(ps.ra)["ASV3120",]#Sedimenticola sp. (Gammaproteobacteria)


#look at sequence that matches to R pachyptila symb
#subset R to just region that aligns with seqs
Rsub <- subseq(R, 424, 676)

#top 1 match from pairwise alignment
BrowseSeqs(c(S["ASV3120"], Rsub), "Riftia.bestmatch.html", openURL = TRUE)#save sequences to file- color coded to see trends

#Alvinella pompejana
#https://www.ncbi.nlm.nih.gov/nuccore/L35523
#Haddad et al. 1995: https://doi.org/10.1128/aem.61.5.1679-1687.1995.
R <- readDNAStringSet("../../Alvinella_pompejana_symbiont_16S_sequence_APG5.fasta")#dowloaded symbiont seq

pa <-  pairwiseAlignment(S, R, type = "overlap")
summary(pa)
pa #overlap starts at position 499
matchper <- data.frame("ASV" = names(S), "perID" = nmatch(pa)/width(S), "score" = score(pa))

#ASV w/ best match percent
matchper[matchper$perID == max(matchper$perID),]#ASV3236

#top 5
head(matchper[order(matchper$perID, decreasing = TRUE), ], n=5)
#    ASV     perID    score
#   ASV3236 0.9683794 438.3350
#   ASV3091 0.9525692 406.8108
#   ASV843 0.9486166 398.9298
#   ASV2918 0.9486166 398.9298
#   ASV401 0.9446640 391.0487

#highest score
matchper[matchper$score == max(matchper$score),]

head(matchper[order(matchper$score, decreasing = TRUE), ], n=5)
#    ASV     perID    score
#   ASV3236 0.9683794 438.3350
#   ASV3091 0.9525692 406.8108
#   ASV2918 0.9486166 398.9298
#   ASV843 0.9486166 398.9298
#   ASV401 0.9446640 391.0487

tiff("figs/Boxplot_Ticasands_AlvinellasymRA_final.tiff", width = 1875, height = 1562, res = 300)
boxplot(otu_table(ps.ra)[,"ASV3236"]*100 ~ sample_data(ps.ra)$Pursed, 
        ylab =" ", xlab = "Biofilm Age", las = 1, 
        col = c("#7e03a8", "#f89540"), par(mar = c(3, 6, 2, 2) + 0.1))
title(ylab = "A. pompejana symbiont RA (%)", line = 3); 
dev.off()

#check out RAs for each sandwich
otu_table(ps.ra)[,"ASV3236"] #only present in 2 samples (both fresh, in mussel zone)

#check out taxonomy
tax_table(ps.ra)["ASV3236",] #sulfurovum

#look at sequences that match to A pompejana symb
#subset R to just region that aligns with seqs
Rsub <- subseq(R, 499, 751)

#top 1 match from pairwise alignment
BrowseSeqs(c(S["ASV3236"], Rsub), "Alvinella.bestmatch.html", openURL = TRUE)#save sequences to file- color coded to see trends


##### correlations between animals and microbes? #####
#propr - after Murdock et al 2021: https://www.nature.com/articles/s43705-021-00031-1
#clr transform both 16S and faunal count data separately, then do propr

#### colonists
#open file with animal count data
colonists <- read.csv("RR2102_colonist_counts_SubKnownMorphs.csv")
rownames(colonists) <- colonists$ID

#subset useful information and only samples that also have 16S data
colsub <- colonists[rownames(colonists) %in% rownames(sampdat), c(1,6,11,22:76)]

#colonist counts raw
colcounts <- colsub[,c(5:58)]

#only look at morphs with greater than 5 counts and in at least 3 sandwiches
colcountkeep <- names(which(colSums(colcounts) >= 5))
colprevkeep <- names(which(colSums(colcounts > 0) >= 3))
colkeep <- intersect(colcountkeep, colprevkeep)

#only higher abundance and prevalence colonists
colcounts <- colcounts[,colkeep]

#get morphotype info for colonists
colmorphinfo <- read.csv("Morphcode_info.csv")

# Imputation of zeros
cols_cmultRepl<-cmultRepl(colcounts, method="CZM", output = "p-counts")
# Centered log-ratio transformation
clr_cols<-apply(cols_cmultRepl, 2, function(x){log(x)-mean(log(x))})

##### microbes, ASV level
#filter ASV level more than ps.f1
#Prevalence of a given ASV
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
#Max relative abundance of a given ASV
relabun = apply(otu_table(ps.ra), 
                MARGIN = ifelse(taxa_are_rows(ps.ra), yes = 1, no = 2),
                max)

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    MaxRelAbundance = relabun,
                    tax_table(ps))

# Execute rel abun filter, using `prune_taxa()` function
# Define abundance and prevelence threshold for taxa - must be at least 0.5% in at least one sample and be present in at least 3 samples
relAbundanceThreshold <- 0.005
prevalenceThreshold <- 3

keepTaxa <- rownames(prevdf)[(prevdf$MaxRelAbundance >= relAbundanceThreshold &
                                prevdf$Prevalence >= prevalenceThreshold)]

ps.f1.v2 <- prune_taxa(keepTaxa, ps) #ps filt obj

sample_sums(ps.f1.v2)/sample_sums(ps)#covers 53-76% of community with 0.005 RA, 64-94% with 0.001 RA

filteredtaxa <- data.frame(tax_table(ps.f1.v2))

# Imputation of zeros
ASVf1_cmultRepl<-cmultRepl(as.data.frame(otu_table(ps.f1.v2)), method="CZM", output = "p-counts")
# Centered log-ratio transformation
clr_ASVf1<-apply(ASVf1_cmultRepl, 2, function(x){log(x)-mean(log(x))})

#combine colonist and ASV transformed count tables
colASV_cmultRepl <- cbind(cols_cmultRepl, ASVf1_cmultRepl)


#### Proportionality analysis
colASV.propr<-propr(counts=colASV_cmultRepl, metric = "rho", p=1000)
colASVCorrALL<-getResults(colASV.propr)

write.csv(colASVCorrALL, "RR2102_ProprAnalysis.csv")

colASVCorrsub <- subset(colASVCorrALL, colASVCorrALL$propr >= 0.60 | colASVCorrALL$propr <= -0.60)


#circos plots to look at proportionality between important organisms
colASVcorrMatrix <- getMatrix(colASV.propr)

#change all values <abs(0.6) to NA??
Matrixsig <- colASVcorrMatrix
Matrixsig[abs(Matrixsig)<0.6] <- NA

Matrixsig["BIV1",]
Matrixsig["NEC1",]
Matrixsig["POL6",]
Matrixsig["POL16",]

colgroup = structure(as.vector(colmorphinfo[colmorphinfo$Code %in% colkeep,"Simp_morph"]), 
                     names = colnames(colcounts))

ASVgroup = structure(as.vector(filteredtaxa[keepTaxa,"Family"]), 
                     names = rownames(Matrixsig[keepTaxa,]))

#make matrix that is only colonists vs ASVs (non-symmetrical)
MatrixcolvASVsig <- Matrixsig[colnames(colcounts),keepTaxa]
#to save simplified matrix with only colonists vs ASVs with sig correlations
Matrix2save <- data.frame(MatrixcolvASVsig)
Matrix2save <- Matrix2save[colSums(!is.na(Matrix2save)) > 0]
Matrix2save <- Matrix2save[rowSums(is.na(Matrix2save)) != ncol(Matrix2save),  ]

write.csv(Matrix2save, "ProprMatrix_rho0.6_colvASVs_ASVRA0.005.csv")

#numb of interactions?
sum(!is.na(Matrix2save))#51 interactions
sum(!is.na(Matrix2save) & Matrix2save > 0) # 17 positive
sum(!is.na(Matrix2save) & Matrix2save < 0) # 34 negative

#fnc to get all unique combinations of groups
expand.grid.unique <- function(x, y, include.equals=FALSE)
{
  x <- unique(x)
  
  y <- unique(y)
  
  g <- function(i)
  {
    z <- setdiff(y, x[seq_len(i-include.equals)])
    
    if(length(z)) cbind(x[i], z, deparse.level=0)
  }
  
  do.call(rbind, lapply(seq_along(x), g))
}


combos <- expand.grid.unique(unique(colgroup), unique(ASVgroup[keepTaxa]), include.equals = TRUE)
combos[is.na(combos)] <- "Other"

df.int <- data.frame("group1" = combos[,1], "group2" = combos[,2])

for (i in 1:nrow(df.int)){
  g1 = colgroup %in% df.int$group1[i]
  g2 = ASVgroup %in% df.int$group2[i]
  num.pos <- sum(MatrixcolvASVsig[g1 ,g2] > 0, na.rm = TRUE)
  num.neg <- sum(MatrixcolvASVsig[g1 ,g2] < 0, na.rm = TRUE)
  df.int$pos[i] <- num.pos
  df.int$neg[i] <- num.neg*-1 #make this value negative
}

df.int <- df.int[rowSums(abs(df.int[,3:4])) > 0,]

stack <- melt(df.int, id.vars=1:2)

unique(stack$group1)
stack$group1 <- factor(stack$group1, 
                       levels = c("GAS","BV","NEC","POL","ISO","COP","FOR"))

unique(stack$group2)
survfams <- filteredtaxa[filteredtaxa$Family %in% unique(stack$group2),]
survfams <- survfams[,1:5]
survfams <- survfams[!duplicated(survfams), ]

stack$group2 <- factor(stack$group2, 
                       levels = c("Sulfurovaceae", "Sulfurimonadaceae", "Sulfurospirillaceae", #Campylo
                                  "Thiotrichaceae", "Colwelliaceae","Ectothiorhodospiraceae", #Gammaproteo
                                  "Crocinitomicaceae", "Flavobacteriaceae", "Cryomorphaceae", "Saprospiraceae", #Bacteroidota
                                  "Rubritaleaceae"))

#use family color pallet from RA figure
print_color(colours, type = "r")
colours
colours2 <- c(colours, "Sulfurospirillaceae" = "#910C00",
              "Ectothiorhodospiraceae" = "#4B4D00", "Saprospiraceae" = "#00641F") 

#copy colonist colors from barchart
#combine with microbe colors
grid.col = c('GAS' = "green3",'BV' = "gold2",'NEC' = "lightskyblue",
             'POL' = "dodgerblue3",'ISO' = "mediumpurple2", 'COP' = "navyblue", 
             'FOR' = 'lightsalmon1',colours2)

colASVorder = c(levels(stack$group1), levels(stack$group2))

#plot, both w/ and without labels - will fix up labels in another program
tiff("Propr_ASVvCol_NumPosNegInts_final.tiff", width = 2187, height = 2187, res = 300)
par(cex = 0.7, mar = c(0,0,0,0))
circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1.4, 1.4))

chordDiagram(stack[,c(1:2,4)], small.gap = 5, grid.col = grid.col,
             col = ifelse(stack$value > 0, "darkblue", "firebrick"), annotationTrack = "grid",
             order  = colASVorder)


for (si in get.all.sector.index()) {
  circos.axis(h = "top", major.at = seq(0,30, by=5), labels.cex = 1, major.tick.length = 0.5, 
              sector.index = si, track.index = 1)
}

# circos.track(track.index = 1, panel.fun = function(x, y) {
#   xlim = get.cell.meta.data("xlim")
#   xplot = get.cell.meta.data("xplot")
#   ylim = get.cell.meta.data("ylim")
#   sector.name = get.cell.meta.data("sector.index")
#   
#   circos.text(mean(xlim), 1, sector.name, facing = "clockwise", niceFacing = FALSE, adj = c(-0.2, 0))
# }, bg.border = NA)

#to check values for xplot
#print(get.cell.meta.data("xplot", sector.index = "BV"))

#circos.info()
circos.clear()
dev.off()

#any correlation between clr transformed pseudo-counts BV vs ASV284?
#recreate clr_ASVf1 (using maxRA = 0.0001 for ASVs - from ps.f1)

# Imputation of zeros
ASVf1_cmultRepl<-cmultRepl(as.data.frame(otu_table(ps.f1)), method="CZM", output = "p-counts")
# Centered log-ratio transformation
clr_ASVf1<-apply(ASVf1_cmultRepl, 2, function(x){log(x)-mean(log(x))})

colASV_cmultRepl <- cbind(cols_cmultRepl, ASVf1_cmultRepl)

#add to dataframe with metadata
colsub$clrBV <- as.numeric(clr_cols[,"BIV1"])
colsub$clrASV284 <- as.numeric(clr_ASVf1[,"ASV284"])

#correlation test
cor.test(colsub$clrASV284,colsub$clrBV, method = "spearman")

#plot clr transformed counts
p <- ggplot(data=colsub, aes(x=clrBV, y=clrASV284,shape=Zone,color=Pursed)) +
  geom_point(size=4) +
  scale_shape_manual(values=c("Alvinellid" = 15,"Riftia" = 18, "Mussel" = 17), name = "Zone") +
  scale_color_manual(values=c("Y" = "#7e03a8", "N" = "#f89540"), name = "Biofilm Age", labels = c( "Fresh", "Established")) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"))
p
