#### Faunal Colonist Analysis ####
### Final analysis for RR2102 Biofilms4larvae MS

## Ladd, T.M., Selci, M., Davis, D.J., Cannon, O., Plowman, C.Q., 
## Schlegel,I.,Inaba, A., Mills, S.W., Vetriani, C., Mullineaux,L.S.,
## Arellano, S.M. (2024).Faunal colonists, including mussel settlers, 
## respond to microbial biofilms at deep-sea hydrothermal vents. Deep-Sea 
## Research Part I

## This code is used to analyze colonist count data on established and 
## fresh biofilm sandwiches across biogenic zones. It includes
## statistical analyses and figure creation.

##### Setup #####

#load packages
library(vegan)
library(ggplot2)
library(pairwiseAdonis)
library(mvabund)
library(reshape2)
library(dplyr)

#clear workspace and set working directory
rm(list=ls())
#setwd() #CHANGE ME to desired directory

#open file - sandwich faunal colonist counts with metadata (not including unknown morphotypes)
TicaSandsSub <- read.csv("RR2102_colonist_counts_SubKnownMorphs.csv")

#make zone a factor with the correct order of levels
TicaSandsSub$Zone <- factor(TicaSandsSub$Zone, levels = c("Alvinellid", "Riftia", "Mussel"))
TicaSandsSub$Pursed <- factor(TicaSandsSub$Pursed, levels = c("Y", "N"))

#square-root transform
#function for any power transformation
pwr_trans<-function(x, trans){ 
  x<- ifelse(x>0,x^(1/trans),0)
  return(x)
}

##### PCoA #####
#take count table (none of the metadata) and sqrt transform all counts
TicaSandsSubSqRt <- apply(TicaSandsSub[,23:76],c(1,2), function(x) pwr_trans(x, 2))

#create a Bray-Curtis distance matrix from transformed count table
BraySqRt <- vegdist(TicaSandsSubSqRt, method = "bray")

#PCoA
PCoA.bray <- cmdscale(d = BraySqRt, eig = TRUE)

#look at plot- get axes % variation explained
ordiplot(PCoA.bray)

PCoA.bray$eig[1]/sum(PCoA.bray$eig)#0.327

PCoA.bray$eig[2]/sum(PCoA.bray$eig)#0.186

#extract PCoA axes
data.scores <- as.data.frame(PCoA.bray$points)
data.scores$ID <-TicaSandsSub$ID #sandwich IDs
data.scores$type <- TicaSandsSub$Pursed  #  biofilmed vs no
data.scores$zone <-TicaSandsSub$Zone # megafaunal zone
data.scores$ID.pairs <-TicaSandsSub$ID.pairs # pairs of sands
data.scores$Recovery.temp <- TicaSandsSub$Recovery_temp.C #temp recovery
data.scores$Start.temp <- TicaSandsSub$Exp_start_temp.C #temp experiment start

##### PCoA plots #####

#hulls by zone
Riftia <- data.scores[data.scores$zone == "Riftia", ][chull(data.scores[data.scores$zone == "Riftia", c("V1", "V2")]), ]  # hull values for Riftia zone
Mussel <- data.scores[data.scores$zone == "Mussel", ][chull(data.scores[data.scores$zone == "Mussel", c("V1", "V2")]), ]  # hull values for Mussel zone
Alvinellid <- data.scores[data.scores$zone == "Alvinellid", ][chull(data.scores[data.scores$zone == "Alvinellid", c("V1", "V2")]), ]  # hull values for Alvinellid zone

hull.data <- rbind(Riftia, Mussel, Alvinellid)

p1 <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=V1,y=V2,fill=zone,group=zone),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=V1,y=V2,shape=type,fill=zone),colour = "black", size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=V1,y=V2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_manual(values=c("Alvinellid" = "#fde725", "Riftia" = "#21918c", "Mussel" = "#440154"), name = "Zone") +
  scale_shape_manual(values=c("Y" = 21, "N" = 25), name = "Biofilm", labels = c("Established", "Fresh")) +
  scale_y_continuous(lim = c(-0.5, 0.25), breaks = c(-0.4, -0.2, 0.0, 0.2)) +
  #scale_x_continuous(lim = c(-0.8, 1.1), breaks = c(-0.5, 0, 0.5, 1.0)) +
  labs(x = "Axis 1 [32.7%]", y = "Axis 2 [18.6%]") +
  coord_equal() +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1)), shape = guide_legend(order = 2)) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p1

#save as tiff - res = 300 for pubs
tiff("Manuscript/Figs final/PCoA_Ticasands_knownmorphs_byzone_SqRtTRANSFORM.tiff", width = 1550, height = 1370, res = 300)
p1
dev.off()

#hulls by biofilm age
Filmed <- data.scores[data.scores$type == "Y", ][chull(data.scores[data.scores$type == "Y", c("V1", "V2")]), ]  # hull values for filmed samples
Unfilmed <- data.scores[data.scores$type == "N", ][chull(data.scores[data.scores$type == "N", c("V1", "V2")]), ]  # hull values for unfilmed samples
hull.data <- rbind(Filmed, Unfilmed)

p2 <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=V1,y=V2,fill=type,group=type),alpha=0.30) + # add the convex hulls
  geom_point(data=data.scores,aes(x=V1,y=V2,shape=zone,fill=type),colour = "black", size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=V1,y=V2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_manual(values=c("Y" = "#7e03a8", "N" = "#f89540"), name = "Biofilm", labels = c("Established", "Fresh")) +
  scale_shape_manual(values=c("Alvinellid" = 22,"Riftia" = 23, "Mussel" = 24), name = "Zone") +
  scale_y_continuous(lim = c(-0.5, 0.25), breaks = c(-0.4, -0.2, 0.0, 0.2)) +
  #scale_x_continuous(lim = c(-0.8, 1.1), breaks = c(-0.5, 0, 0.5, 1.0)) +
  labs(x = "Axis 1 [32.7%]", y = "Axis 2 [18.6%]") +
  coord_equal() +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1)), shape = guide_legend(order = 2)) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p2

#save as tiff - res = 300 for pubs
tiff("Manuscript/Figs final/PCoA_Ticasands_knownmorphs_bybiofilm_SqRtTRANSFORM.tiff", width = 1550, height = 1370, res = 300)
p2
dev.off()


#fit temp to PCoA axes
PCoA.fit <- envfit(data.scores[,1:2], data.scores[,7:8], permutations = 999, na.rm = TRUE) 
PCoA.fit

en.data.scores = as.data.frame(scores(PCoA.fit, "vectors"))*0.5


p3 <- ggplot() + 
  geom_point(data=data.scores,aes(x=V1,y=V2,shape=zone,fill=Recovery.temp),colour = "black", size=4) + # add the point markers
  geom_text(data=data.scores,aes(x=V1,y=V2,label=ID.pairs),size=2.5,vjust=-1,hjust=0.5) +  # add the site labels
  scale_fill_viridis_c(option = "plasma", name = expression(paste("Temperature (",degree*C,")"))) +
  scale_shape_manual(values=c("Alvinellid" = 22,"Riftia" = 23, "Mussel" = 24), name = "Zone") +
  geom_segment(aes(x = 0, y = 0, xend = V1, yend = V2), 
               data = en.data.scores, linewidth = 0.5, colour = "grey30", 
               arrow = arrow(length = unit(0.5,"cm"))) +
  geom_text(data = en.data.scores[1,], aes(x = V1 + 0.04*sign(V1), y = V2 + 0.04*sign(V2)), 
            colour = "black", size = 3, hjust = 0.1,
            label = "Start") +
  geom_text(data = en.data.scores[2,], aes(x = V1 + 0.04*sign(V1), y = V2 - 0.04*sign(V2)), 
            colour = "black", size = 3, hjust = 0.1,
            label = "Recovery") +
  scale_y_continuous(lim = c(-0.5, 0.25), breaks = c(-0.4, -0.2, 0.0, 0.2)) +
  #scale_x_continuous(lim = c(-0.8, 1.1), breaks = c(-0.5, 0, 0.5, 1.0)) +
  labs(x = "Axis 1 [32.7%]", y = "Axis 2 [18.6%]") +
  coord_equal() +
  #guides(fill = guide_legend(order = 1, override.aes = list(shape = NA, alpha = 1)), shape = guide_legend(order = 2)) +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_text(size=12, colour="black"),  axis.text.y = element_text(size=12, colour="black"), 
        axis.title.x = element_text(size=14, colour="black"), axis.title.y = element_text(size=14, colour="black"))
p3

tiff("Manuscript/Figs final/PCoA_Ticasands_knownmorphs_byzonewtemp_SqRtTRANSFORM.tiff", width = 1550, height = 1370, res = 300)
p3
dev.off()


##### Differences in community composition - Adonis #####

#Tica sands
#full model
(res_adonis <- adonis2(BraySqRt ~ Zone*Pursed, data = TicaSandsSub)) #sig effects of zone and pursed, no interaction effects
#remove interaction
(res_adonis <- adonis2(BraySqRt ~ Zone + Pursed, data = TicaSandsSub)) #sig effects of zone and pursed

#pairwise test to look at differences between zones
pairwise.adonis2(BraySqRt ~ Zone, data = TicaSandsSub)

#check for differences in multivariate dispersion based on biofilm or zone (like homogeneity of variances)
beta <- betadisper(BraySqRt, TicaSandsSub$Pursed)

permutest(beta) #no sig heterogeneity of dispersion

beta <- betadisper(BraySqRt, TicaSandsSub$Zone)

permutest(beta) #no sig heterogeneity of dispersion


##### Multivariate GLM methods #####
TicaSandsabund <- mvabund(TicaSandsSub[,23:76]) 

#explore
meanvar.plot(TicaSandsabund)
plot(TicaSandsabund ~ as.factor(TicaSandsSub$Pursed), cex.axis = 0.8, cex = 0.8)
plot(TicaSandsabund ~ as.factor(TicaSandsSub$Zone), cex.axis = 0.8, cex = 0.8)


#fit glm
TicaSands.nb <- manyglm(TicaSandsabund~TicaSandsSub$Zone*TicaSandsSub$Pursed,
                        family="negative.binomial")
#plots
plot(TicaSands.nb)

res <- residuals(TicaSands.nb)
qqnorm(res)
qqline(res,col="red")

#test effects
anova(TicaSands.nb, test = "score", cor.type = "shrink") #non significant interaction term
#univariate "species" tests
t <- anova(TicaSands.nb, p.uni = "adjusted", test = "score", cor.type = "shrink")
which(t$uni.p[2,] < 0.1) #zone
which(t$uni.p[3,] < 0.1) #pursed

#AIC
TicaSands.nb$AICsum

#compare to model without interaction?
TicaSandsXint.nb <- manyglm(TicaSandsabund~TicaSandsSub$Zone + TicaSandsSub$Pursed, 
                            family="negative.binomial")
#AIC
TicaSandsXint.nb$AICsum #lower AIC
TicaSands.nb$AICsum

#is this the best model??
anova(TicaSandsXint.nb, test = "score", cor.type = "shrink")
#univariate "species" tests
tXint <- anova(TicaSandsXint.nb, p.uni = "adjusted", test = "score", cor.type = "shrink")

which(tXint$uni.p[2,] < 0.1) #zone
which(tXint$uni.p[3,] < 0.1) #pursed

plot(TicaSandsXint.nb)

res <- residuals(TicaSandsXint.nb)
qqnorm(res)
qqline(res,col="red")

#test whether "full" model is different than after dropping terms?
(mcomp <- anova(TicaSands.nb, TicaSandsXint.nb)) #not sig different- can drop the interaction term

#does including pursed explain additional variation in the model?
TicaSandsZone.nb <- manyglm(TicaSandsabund~TicaSandsSub$Zone, 
                            family="negative.binomial")
anova(TicaSandsZone.nb, TicaSands.nb)
#AIC
TicaSandsZone.nb$AICsum

#test whether "full" model is different than after dropping terms?
mcomp <- anova(TicaSandsXint.nb, TicaSandsZone.nb) #significant

##### plot barcharts of community comp #####
#open morphotype file
morphinfo <- read.csv("Morphcode_info.csv")
#open count file that does not remove unknown groups
TicaSands <- read.csv("RR2102_colonist_counts.csv")

#create stacked df
stackmorph <- melt(TicaSands[,c(1,22,6,11,23:80)], id.vars=c("ID", "ID.pairs", "Zone", "Pursed"), 
                   value.name = "count", variable.name = "morphotype")

#add tax group to stacked df
stackmorph$group <- as.factor(morphinfo$Taxonomic_group[match(stackmorph$morphotype, morphinfo$Code)])
unique(stackmorph$group)
stackmorph$group <- factor(stackmorph$group, levels = c("Gastropoda", "Bivalvia","Polychaeta", "Crustacea",        
                                                        "Rhizaria", "Echinodermata","Nematoda", "Arachnida"))

stackmorph$simpmorph <- as.factor(morphinfo$Simp_morph_name[match(stackmorph$morphotype, morphinfo$Code)])
stackmorph$simpmorph <- factor(stackmorph$simpmorph, 
                               levels = c("Gastropod", "Gastropod veliger", "Juvenile mussel",
                                          "Bivalve veliger", "Polychaete", "Nectochaete",
                                          "Amphipod","Isopod", "Copepod" ,"Copepod nauplius","Ostracod",
                                          "Foraminifera","Ophiuroid","Nematode" ,"Marine mite"))



stackmorph$Zone <- factor(stackmorph$Zone, levels = c("Alvinellid", "Riftia", "Mussel"))
stackmorph$Pursed <- factor(stackmorph$Pursed, levels = c("Y", "N"))
levels(stackmorph$Pursed) <- c("Established", "Fresh")

stackmorph$ID.pairs <- factor(stackmorph$ID.pairs)
stackmorph$ID <- factor(stackmorph$ID, levels = c("PW", "PP10","PX", "PP1", "PM", "PP2", "PS", 
                                                  "PP3", "PR", "PP5", "PO", "PP11", "PQ", "PP12"))


#plots
morphcolors = c('Gastropod' = "green3",'Gastropod veliger' = "darkgreen",
                "Juvenile mussel" = "darkorange2",'Bivalve veliger' = "gold2", 'Polychaete'  = "dodgerblue3",
                 'Nectochaete'= "lightskyblue",'Amphipod' = "darkorchid4", 'Isopod' = "mediumpurple2", 'Copepod' = "navyblue", 
                'Copepod nauplius' ="turquoise" , 'Ostracod' = "maroon", 'Foraminifera' = 'lightsalmon1',
                'Nematode' = "hotpink",'Ophiuroid'= 'firebrick4', 'Marine mite' = 'thistle2')

#all colonists by zone and simp morph
p4 <- ggplot(stackmorph[!is.na(stackmorph$group),], aes(fill=simpmorph, y=count, x=ID)) + 
  geom_bar(position=position_stack(reverse = TRUE), stat="identity") +
  scale_fill_manual(name = "Morphotype group", values = morphcolors) +
  scale_y_continuous(expand = c(0,0), lim = c(0, 950), name = "Count") +
  scale_x_discrete(name = " ") +
  facet_grid(cols = vars(Zone), scales="free_x", space = "free_x") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=12, colour="black", angle = 45,  hjust = 1),  
        axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14), 
        axis.title.x = element_text(size=14, colour="black"), 
        axis.title.y = element_text(size=14, colour="black"))
p4

tiff("Barchart_Ticasands_Morphgroupcounts_knownMorphs.tiff", width = 2813, height = 1875, res = 300)
p4
dev.off()

#all colonists RA by zone and simp morph
p5 <- ggplot(stackmorph[!is.na(stackmorph$group),], aes(fill=simpmorph, y=count, x=ID)) + 
  geom_bar(position = position_fill(reverse = TRUE), stat="identity") +
  scale_fill_manual(name = "Morphotype group", values = morphcolors) +
  scale_y_continuous(expand = c(0,0), labels = c("0", "25", "50", "75", "100"), name = "Relative abundance (%)") +
  scale_x_discrete(name = " ") +
  facet_grid(cols = vars(Zone), scales="free_x", space = "free_x") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.text = element_text(size = 12), legend.title = element_text(size = 14),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size=12, colour="black", angle = 45,  hjust = 1),  
        axis.text.y = element_text(size=12, colour="black"), strip.text = element_text(size = 14), 
        axis.title.x = element_text(size=14, colour="black"), 
        axis.title.y = element_text(size=14, colour="black"))
p5

tiff("Barchart_Ticasands_MorphgroupRA_knownMorphs.tiff", width = 2813, height = 1875, res = 300)
p5
dev.off()

##### Paired tests for aggregated and indiv morphotypes #####
## create new df that has new columns with aggregated count data
TicaSandsAgg <- TicaSandsSub
TicaSandsAgg$totcounts <- rowSums(TicaSands[c(23:80)]) #all colonists (adults/juvenile immigrants and settlers)
#define "colonists" and "settlers" - confusing terminology, all are colonists but here C = immigrant colonists
#and S = potential settlers (mainly larvae and nectochaetes)
colonists <- morphinfo$Code[morphinfo$Colonist_Settler %in% "C"]
settlers <- morphinfo$Code[morphinfo$Colonist_Settler %in% "S"]

TicaSandsAgg$colcounts <- rowSums(TicaSands[colnames(TicaSands) %in% colonists])#adult colonists
TicaSandsAgg$larcounts <- rowSums(TicaSands[colnames(TicaSands) %in% settlers])#potential settlers

#aggregated by groups
TicaSandsAgg$GAScounts <- rowSums(TicaSands[colnames(TicaSands) %in% 
                                              morphinfo$Code[morphinfo$Simp_morph %in% "GAS"]])
TicaSandsAgg$GVcounts <- rowSums(TicaSands[colnames(TicaSands) %in% 
                                             morphinfo$Code[morphinfo$Simp_morph %in% "GV"]])
TicaSandsAgg$POLcounts <- rowSums(TicaSands[colnames(TicaSands) %in% 
                                              morphinfo$Code[morphinfo$Simp_morph %in% "POL"]])
TicaSandsAgg$NEcounts <- rowSums(TicaSands[colnames(TicaSands) %in% 
                                             morphinfo$Code[morphinfo$Simp_morph %in% "NEC"]])
TicaSandsAgg$perlar <- TicaSandsAgg$larcounts/TicaSandsAgg$totcounts

TicaSandsAgg$shannondiv <- diversity(TicaSandsSub[,23:76], index = "shannon")

#what is the maximum percent of counts in a sand that are "potential settlers"
max(TicaSandsAgg$perlar)*100

#separate TicaSands into two dfs to substract paired counts
TicaEst <- subset(TicaSandsAgg, TicaSandsAgg$Pursed == "Y")
TicaEst <- TicaEst[order(TicaEst$ID.pairs),]

TicaFres <- subset(TicaSandsAgg, TicaSandsAgg$Pursed == "N")
TicaFres <- TicaFres[order(TicaFres$ID.pairs),]

TicaPairDiff <- TicaEst[,23:83] - TicaFres[,23:83]
TicaPairDiff$ID.pairs <- TicaEst$ID.pairs
TicaPairDiff$IDEst <- TicaEst$ID
TicaPairDiff$IDFres <- TicaFres$ID
TicaPairDiff$Zone <- TicaEst$Zone
row.names(TicaPairDiff) <- TicaPairDiff$ID.pairs

#loop through columns and run tests
#only look at morphs with greater than 5 counts and in at least 3 sandwiches
morphcountkeep <- names(which(colSums(TicaSandsAgg[23:83]) >= 5))
morphprevkeep <- names(which(colSums(TicaSandsAgg[23:83] > 0) > 2))
morphkeep <- intersect(morphcountkeep, morphprevkeep)

TicaPairDiff2 <- TicaPairDiff[,morphkeep]

TicaPairStats <- data.frame("morph" = colnames(TicaPairDiff2), "shapiro_p" = NA, 
                            "t_p" = NA, "wilcox_p" = NA, "shapiro_p2" = NA,
                            "t_p2" = NA, "wilcox_p2" = NA)

for (i in 1:ncol(TicaPairDiff2)) {
  shapiro <- shapiro.test(TicaPairDiff2[,i])
  TicaPairStats$shapiro_p[i] <- shapiro$p.value
  t<- t.test(TicaPairDiff2[,i], mu = 0, alternative = "two.sided")
  TicaPairStats$t_p[i] <- t$p.value
  wilcox <- exactRankTests::wilcox.exact(TicaPairDiff2[,i], mu = 0, alternative = "two.sided")
  TicaPairStats$wilcox_p[i] <- wilcox$p.value
}

#excluding alvinellid zone
for (i in 1:ncol(TicaPairDiff2)) {
  shapiro <- shapiro.test(TicaPairDiff2[-1,i])
  TicaPairStats$shapiro_p2[i] <- shapiro$p.value
  t<- t.test(TicaPairDiff2[-1,i], mu = 0, alternative = "two.sided")
  TicaPairStats$t_p2[i] <- t$p.value
  wilcox <- exactRankTests::wilcox.exact(TicaPairDiff2[-1,i], mu = 0, alternative = "two.sided")
  TicaPairStats$wilcox_p2[i] <- wilcox$p.value
}

write.csv(TicaPairStats, "TicaPairCountStats_atleast5prev3.csv")

##### heatmap of morphotype counts by biofilm age #####

#make all 0 values NA
stackmorph$countx <- stackmorph$count
stackmorph$countx[stackmorph$countx == 0] <- NA

yorder <- unique(stackmorph$morphotype[order(stackmorph$simpmorph)])
stackmorph$morphotype <- factor(stackmorph$morphotype, levels = yorder)

#just morphotypes with >5 total counts in at least 3 sands
stackmorphsub <-stackmorph[stackmorph$morphotype %in% morphkeep,]
stackmorphsub <- droplevels(stackmorphsub)

#change order of pursed levels to change order on plot
stackmorphsub$Pursed <- factor(stackmorphsub$Pursed, levels = c("Fresh", "Established"))

#plot
p6 <- ggplot(stackmorphsub, aes(x=ID, y=morphotype, fill=countx)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "navy") +
  scale_fill_gradient(low = "#00127A", high = "#66CCFF", na.value = "black", trans = "log10",
                      breaks = c(1, 10, 100), name = "Count") +
  scale_x_discrete(expand = c(0, 0), labels = levels(stackmorph$ID.pairs)) +
  scale_y_discrete(limits = rev, expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 7.5), ylim = c(0.5,25.5), clip = "off") +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 0.5, yend = 1.5), colour = "hotpink", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 1.5, yend = 2.5), colour = "lightsalmon1", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 2.5, yend = 6.5), colour = "navyblue", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 6.5, yend = 7.5), colour = "mediumpurple2", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 7.5, yend = 8.5), colour = "darkorchid4", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 8.5, yend = 11.5), colour = "lightskyblue", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 11.5, yend = 18.5), colour = "dodgerblue3", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 18.5, yend = 19.5), colour = "gold2", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Fresh"),
               aes(x = -2.3, xend = -2.3,y = 19.5, yend = 25.5), colour = "green3", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"), #otherside
               aes(x = 7.9, xend = 7.9,y = 0.5, yend = 1.5), colour = "hotpink", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 1.5, yend = 2.5), colour = "lightsalmon1", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 2.5, yend = 6.5), colour = "navyblue", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 6.5, yend = 7.5), colour = "mediumpurple2", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 7.5, yend = 8.5), colour = "darkorchid4", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 8.5, yend = 11.5), colour = "lightskyblue", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 11.5, yend = 18.5), colour = "dodgerblue3", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 18.5, yend = 19.5), colour = "gold2", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub %>% filter(Pursed == "Established"),
               aes(x = 7.9, xend = 7.9,y = 19.5, yend = 25.5), colour = "green3", 
               size = 5, show.legend = FALSE, inherit.aes = FALSE) +
  geom_rect(data = stackmorphsub %>% filter(Pursed == "Established"),
            aes(xmin = -9.9, xmax = 8.2,ymin = 18.5, ymax = 19.5),
            fill = NA, color = "red", linewidth = 1) +
  geom_rect(data = stackmorphsub %>% filter(Pursed == "Established"),
            aes(xmin = -9.9, xmax = 8.2,ymin = 14.5, ymax = 15.5),
            fill = NA, color = "red", linewidth = 1) +
  geom_rect(data = stackmorphsub %>% filter(Pursed == "Established"),
            aes(xmin = -9.9, xmax = 8.2,ymin = 11.5, ymax = 12.5),
            fill = NA, color = "red", linewidth = 1) +
  geom_rect(data = stackmorphsub %>% filter(Pursed == "Established"),
            aes(xmin = -9.9, xmax = 8.2,ymin = 10.5, ymax = 11.5),
            fill = NA, color = "red", linewidth = 1) +
  geom_segment(data = stackmorphsub,
               aes(x = 0.5, xend = 1.5,y = 0.3, yend = 0.3), colour = "#fde725", 
               size = 2, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub,
               aes(x = 1.5, xend = 3.5,y = 0.3, yend = 0.3), colour = "#21918c", 
               size = 2, show.legend = FALSE, inherit.aes = FALSE) +
  geom_segment(data = stackmorphsub,
               aes(x = 3.5, xend = 7.5,y = 0.3, yend = 0.3), colour = "#440154", 
               size = 2, show.legend = FALSE, inherit.aes = FALSE) +
  theme_classic() +
  facet_grid(~Pursed, scales = "free_x", space = "free") +
  theme(panel.border = element_rect(colour = "black", fill = NA), 
        legend.title = element_text(size = 14), legend.text = element_text(size = 12),
        legend.key = element_blank(), legend.background=element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        strip.text = element_text(size = 14), axis.ticks.x=element_blank(),
        axis.text.x = element_text(size = 12,colour = "black", angle = 60, hjust = 1), 
        plot.margin = margin(20,5.5,5.5,30, "pt"),
        axis.text.y = element_text(size = 12, colour = "black"))
p6


tiff("Heatmap_Ticasands_morphsatleast5prev3_knownMorphs.tiff", width = 1719, height = 2188, res = 300)
p6
dev.off()
