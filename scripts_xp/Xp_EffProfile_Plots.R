####################################
### SUMARY ANALYSIS OF EFFECTOR PROFILE ###
####################################
# Description
# Obejctive
# Pseudocode
# Inputs
# Outputs

####################################

#########    IMPORT FILES AND LIBRARIES    #########  

####################################

### LOAD LIBRARIES
library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)



### LOAD FILES

### Set working directory
setwd("~/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/")



### EFFECTOR LENGTH table 
# Table with effector length, used to calculate query coverage
# effector_length.txt table
eff_length <- read.csv(file= "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/01152019_Xp_SCRI_269_vs_125_effectors/125_effector_length_2.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)  
colnames(eff_length)[which(names(eff_length) == "Effector")] <- "queryID"     # Change column name "Effector" for "queryID", to use it as Key
head(eff_length)

### STRAIN list and IMG key 
str_id <- read.csv(file= "strains_key2", sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
colnames(str_id) <- c("strain", "goldID")
head(str_id)
nrow(str_id)


### effector_allele-type_singletons.txt
alleles_all <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/effector_allele-type_singletons_clean2_editXopP.txt", 
                          sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
#alleles_all$haplotype <- as.factor(alleles_all$haplotype)
str(alleles_all)
nrow(alleles_all)
length(unique(alleles_all$strain))

### METADATA table
metadata <- read.csv(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/metadata_nmds2_Sujancolors.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(metadata)
#header <- c("strain", "country", "continent", "country_col", "continent_col")
#colnames(metadata) <- header
#rownames(metadata) <- metadata[,1]
head(metadata)

### ORDER of strains in PHYLOGENETIC TREE
phyl_order <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/Phylogenetic_tree_order.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(phyl_order)


### Freq Xp nucleotide 
### eff_freq_singletons.txt
file <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/eff_freq_singletons3.txt", #version 3 has the correction for XopP and the XopPcopy
                   sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)


### Numeric matrix (if available)
nmres <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/haplotypes_summary3.txt", 
                          sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)
str(nmres)
#row.names(nmres) <- metadata$strain # Make sure it has the correct order (same as phylogeneitc tree)
head(nmres)
nmres <- data.matrix(nmres)
str(nmres)
class(nmres)
ncol(nmres)
nrow(nmres)


####################################

### DATA FORMATING ###

####################################

##### REMOVE DUPLICATES if still present (Optional, OR upload "_clean" file)

### effector_allele-type_singletons.txt
res <- alleles_all
nrow(res)
# Clean results table from XopP, it has hit with two genes, I label one of them as Xop-2
res_xop <- res[which(res$effector=="XopP"),]    # Subset and edit duplicate hits from rows of effector XopP (2 hits)
nrow(res_xop)
xop1 <- res_xop[which(duplicated(res_xop$strain)==TRUE),]     # Duplicated rows, more than one hit per strain
nrow(xop1)  # 217
xop2 <- res_xop[which(duplicated(res_xop$strain)==FALSE),]
nrow(xop2)  # 269
head(res_xop)
head(xop1)
# Change the minoritary hit to XopP-2
res$effector[which(res$locus_tag%in%xop1$locus_tag)] <- "XopP-2"

# Clean results table from XopAQ
res_aq <- res[which(res$effector=="XopAQ"),]    # Subset and edit duplicate hits from rows of effector XopP (2 hits)
nrow(res_aq)
res_aq
aq2 <- res_aq[which(duplicated(res_aq$goldID)==TRUE),]     # Duplicated rows, more than one hit per strain
nrow(aq2)  # 217
res <- res[-(which(res$locus_tag %in% aq2$locus_tag)),]

# Clean results table from XopE2
res_e2 <- res[which(res$effector=="XopE2"),]    # Subset and edit duplicate hits from rows of effector XopP (2 hits)
nrow(res_e2)
res_e2[res_e2$strain=="Bzl13",]
e2 <- res_e2[which(duplicated(res_e2$goldID)==TRUE),]     # Duplicated rows, more than one hit per strain
nrow(e2)  # 2
res <- res[-(which(res$locus_tag %in% e2$locus_tag)),]
nrow(res)

#Remove row from GEV1063_2 (wrong assembly)
which(res$strain=="GEV1063_2")
nrow(res)
alleles_all <- res[-(which(res$strain=="GEV1063_2")),]
nrow(alleles_all) # Compare with nrow(res), should be less rows, if GEV1063_2 was still present
head(alleles_all)

# Back to object alleles_all
alleles_all <- res

# Save to file
write.table(alleles_all, file = "effector_allele-type_singletons_clean2.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE) 


####################################

### ORGANIZE RESULTS AS NUMERIC MATRIX

# Format table into numeric matric
# Object 'res' will be needed for plots, so run again the previous part if needed, or convert from object 'alleles_all'
# Or input haplotypes_summary table

res <- alleles_all

head(res)
nrow(res)
length(unique(res$locus_tag))
res_clean <- res[,c(1,2,5)]
head(res_clean)
str(res_clean)
results <- spread(data = res_clean, key = effector, value = haplotype, fill = "0")

{
  #res <- results    #alleles_tab object comes from alleles_summary_table.R, some strains have different name
  class(results)
  str(results)
  
  #change NAs for Os
  results[is.na(results)] <- 0
  
  #add strains as rownames and strains eliminate column from table
  row.names(results) <- results$strain
  head(results)
  results <- results[ , -which(names(results) %in% c("strain"))]
  
  #remove empty columns (no hit) and save to "empty" object
  #empty <- results[which(colSums(results)==0)]
  #res <- results[,-which(colSums(results)==0)]
  
  nmres <- data.matrix(nmres_clean)
  str(nmres)
  class(nmres)
  ncol(nmres)
  nrow(nmres)
}
head(results)
head(nmres)
nrow(nmres)


write.table(nmres, file = "haplotypes_summary3.txt", col.names = TRUE, row.names = TRUE, sep= "\t", quote = FALSE)


### ORDER STRAINS BY PHYLOGENETIC TREE

# Order numeric matric (object "nmres"), according to phylogenetic tree
head(nmres) 
nmres[which(!(rownames(nmres) %in% phyl_order$strain)),] # Confirm all strains have same name among datasets, Expect integer 0 if all match
rownames(nmres) <-  str_remove(rownames(nmres), "_2")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV839", "GEV839D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV893", "GEV893D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV904", "GEV904D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV909", "GEV909D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV915", "GEV915D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV917", "GEV917D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV936", "GEV936D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV940", "GEV940D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV968", "GEV968D")
rownames(nmres) <-  str_replace(rownames(nmres), "GEV993", "GEV993D")
nmres_clean <- nmres[match(phyl_order$strain,rownames(nmres)),]
nmres <- nmres_clean
head(nmres)


# Order metadata table to match the order and number of strains with respect with row names from matrix
head(metadata)
metadata[which(!(metadata$strain %in% phyl_order$strain)),] # Confirm all strains have same name among datasets, Expect integer 0 if all match
metadata_clean <- metadata[match(phyl_order$strain,metadata$strain),]
metadata <- metadata_clean
head(metadata)


### OPTIONAL ###
# Order numeric matrix and metadata table by phylogenetic group (not the phylogenetic tree)
# By BAPS levels
nmres_clean <- nmres[order(metadata$baps_levels),]
metadata_clean <- metadata[order(metadata$baps_levels),]

# By Effector profile
nmres_eff <- nmres_clean[order(metadata_clean$eff_cluster),]
metadata_eff <- metadata_clean[order(metadata_clean$eff_cluster),]
### ### ###


### CLEAN DATASET - Only XOP effectors (And a few Avr genes)

## Input "nmres"
head(nmres)
nrow(nmres)
ncol(nmres)
u_eff <- colnames(nmres)

# Subset effectors for analysis
#clean_eff <- as.character(u_eff[c(1,20:50)]) # Use for full dataset, not for files with "_xops"
clean_eff <- u_eff
length(clean_eff)

# Subset numeric matrix for effectors 
nmres_clean <- nmres[,which(colnames(nmres) %in% clean_eff)]
ncol(nmres[,which(colnames(nmres) %in% clean_eff)])     # 6512 sequences
ncol(nmres_clean)
nrow(nmres_clean)

#Save to file
write.table(nmres_clean, file = "haplotypes_summary_xops3.txt", col.names = TRUE, row.names = TRUE, sep= "\t", quote = FALSE)

## Input "alleles_all"
#head(alleles_all)
#nrow(alleles_all)
alleles_clean <- alleles_all
head(alleles_clean)

# Subset alleles table for effectors 
#alleles_clean <- alleles_all[(alleles_all$effector %in% clean_eff),]
nrow(alleles_all[(alleles_all$effector %in% clean_eff),])     # 6515 sequences
ncol(alleles_clean)
nrow(alleles_clean)


#Save to file
write.table(alleles_clean, file = "alleles_summary_xops3.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)




####################################

### DATA ANALYSIS AND PLOTS ###

####################################

##################
### PLOT FREQUENCY AND ALLELE DIVERSITY (allele number/type), 
# Xp aa sequence analysis
# From effector_frequency.R

# Input table
str(alleles_clean)
#alleles_clean <- alleles_all

# Use manual palette of colors, useful to include many colors
# (When palletes have fewer colors that needed, figures/plots return white where colors are not enough)
mycols <- c("mediumpurple4", "cornflowerblue", "mediumseagreen", "olivedrab3", "darkolivegreen3", "darkolivegreen1", "darkseagreen1", 
            "seashell", "lightgoldenrodyellow", "khaki2", "yellow2",  "goldenrod2", "darksalmon",  "coral2", "tomato2", "chocolate1",
            "palevioletred", "hotpink2", "coral3", "darkorchid1", "deeppink3", "firebrick3", "darkmagenta", "violetred4", "firebrick4", "red", "black")

# Histogram Plot 
pdf("Plot_Xp_effFrequency_All_Feb21.pdf", width= 10, height = 15)
ggplot(alleles_clean, aes(x = effector)) + 
  geom_histogram(aes(fill = as.factor(haplotype)), stat = "count") +
  coord_flip() +
  scale_fill_manual(values=mycols) +
  scale_x_discrete(name = "Effector") +
  scale_y_continuous(name = "Frequency", breaks = seq(0, 270, 50), limits = c(0, 270)) +
  theme_classic()
dev.off()



##################

### HEATMAP PLOT with EFFECTOR PROFILES
# From Heatmap2_052019_alleles_ordered-strains_continent.R


### Load libraries and files
# heatmap.2() from package 'gplots'
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

library(lattice)


### INPUTS

# Numeric matrix
#nmres <- nmres_clean
head(nmres)
class(nmres) # Make sure it is a numeric matrix

# Optional
#nmres <- nmres_clean    # Order by phylogenetic groups
nmres <- nmres_eff      # Order by effector profiles
head(nmres)
###

# Metadata information
head(metadata)

#Optional 
#metadata <- metadata_clean # Order by phylogenetic groups
metadata <- metadata_eff   # Order by effector profiles
head(metadata)
###

### LEGEND / COLORS
### Uncomment the corresponding objects with colors that will be needed on the legend
### Also modify in the code for legend when making the plot

## For country, continent, county
# Colors USA vs Others: "blue" and "darkgray"
#mycol1 <- c("blue", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray", "darkgray")

# Colors for all  13 countries
mycol2 <- c("blue", "coral", "goldenrod1", "red", "yellow", "tan", "turquoise4", "steelblue2", "sienna1", "pink3", "orchid3", "lightpink", "limegreen")

# Colors for all 7 continents
mycol3 <- c("blue", "yellow", "lightpink", "turquoise2", "orchid3", "steelblue2")

# Colors for BAPS groups
#mycolBaps <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#6A3D9A", "#00441B")
mycolBaps_Sujan <- c("#E4211C", "#4DAF4A", "#984EA3", "#FF7400", "#FFFB33", "#A65629", "#F781BF", "#6A3D9A", "#004B1C", "#377EB8") %>%
  setNames(1:10)
# Colors for effector clusters
ltr_ordr <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
mycoleff <- c("#FF7F50", "#0000FF", "#5CACEE", "#FFC125", "#B8DC8C", "#CD919E", "#FFFF00", "#FF0000", "#D2B48C")%>%
  setNames(ltr_ordr)

## Colors for different alleles

#coul =  colorRampPalette(c("white", "navajowhite", "goldenrod2", "orange", "darkorange", "orangered1", "darkred"))

coul = colorRampPalette(c("white", "mediumpurple4", "cornflowerblue", "mediumseagreen", "olivedrab3", "darkolivegreen3", "darkolivegreen1", "darkseagreen1", 
                          "seashell", "lightgoldenrodyellow", "khaki2", "yellow2",  "goldenrod2", "darksalmon",  "coral2", "tomato2", "chocolate1",
                          "palevioletred", "hotpink2", "coral3", "darkorchid1", "deeppink3", "firebrick3", "darkmagenta", "violetred4", "firebrick4"))



### Clustering and distance matrix functions

# Needed to calculate a distance matrix and the hierarchical clustering for the dendograms
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x, method="manhattan") 
# ?dist or ?Distance Matrix Computation 
# Euclidean: Usual distance between the two vectors (2 norm aka L_2), sqrt(sum((x_i - y_i)^2))
# Other methods: "manhattan", "maximum", "canberra", "binary", "minkowski" 

## Perform clustering on rows, to obtain the order of strains after clustering, this is useful for RowSideColors
distance <- distfunc(nmres)     #This distance matrix can be used to calcluate PCA, in script PCA_analysis.R
hcluster <- hclustfunc(distance)
dend1 <- as.dendrogram(hcluster)
#plot(dend1)
#
## To see the dendogram
#str(dend1)
## strains in the order of the dendogram
#labels(dend1)

## Colors for RowSideColors
# Phylogenetic groups
rowcolors <- c()
for(c in 1:length(metadata$baps_levels)){
  temp <- metadata$baps_levels[c]
  rowcolors[c] <- mycolBaps_Sujan[temp]
  
}
rowcolors
metadata$rowcolors <- rowcolors

## Effector groups
ltr_ordr <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")
rowcolors2 <- c()
for(d in 1:length(metadata$eff_cluster)){
  temp <- metadata$eff_cluster[d]
  numtemp <- which(ltr_ordr==temp)
  rowcolors2[d] <- mycoleff[numtemp]
  
}
rowcolors2
metadata$rowcolors2 <- rowcolors2

nmres <- data.matrix(nmres_clean)

### PLOT THE HEATMAP
### heatmap.2() function
{
  pdf(file = "Heatmap_Xp_SCRI_BAPS_All_Feb21.pdf") #Name of the file where the plot will be saved
  
  par(mar=c(8,10,7,7)) 
  h <-heatmap.2(nmres,  cexRow = 0.14, cexCol = 0.5, col = coul, trace = "none",  scale="none",
                key = TRUE, keysize = 1,  density.info = "none", key.ylab = TRUE, key.xlab= "Allele type" , key.title = NA, distfun=distfunc, hclustfun=hclustfunc, 
                RowSideColors = metadata$rowcolors, #For BAPS/phylogenetic groups
                #RowSideColors = metadata$rowcolors2, # For effector profile groups
                breaks = 32, Rowv = FALSE)
  # Legend for phylogenetic groups
  legend(x=-0.55, y=.95, legend=c(1:10), 
       fill = mycolBaps_Sujan, xpd = TRUE, cex=.5 ) #Note: Legend needs to be better fit on the plot!
  # Legend for effector groups
  #legend(x=-0.45, y=.95, legend=c(ltr_ordr), 
  #      fill = mycoleff, xpd = TRUE, cex=.5 ) #Note: Legend needs to be better fit on the plot!
  
  dev.off()
  
}



##################

### NMDS Analysis and Plots

### Load libraries and files
library(tidyr)
library(dplyr)
library(vegan)
library(ggplot2)
library(ggrepel)

### INPUTS
# Numeric matrix and metadata tables should have the strains inthe same order (in this case, according to phylogenetic tree)

# alleles or clean blast results: effector_allele-type_singletons.txt
alleles_all <- alleles_clean
head(alleles_all)
nrow(alleles_all)

# Metadata table
head(metadata)
nrow(metadata)
ncol(metadata)

# Numeric matrix
head(nmres)


### ANALIZE DATA
# Nonmetric Multidimensional Scaling (NMDS)

head(nmres)
str(nmres)

#Remove row for strain GEV1063 and empty effector columns caused by this strain
# Remove strain row
#which(rownames(nmres)=="GEV1063")
#nmres2 <- nmres[-(which(rownames(nmres)=="GEV1063")),]
# Remove effector column
#which(colSums(nmres2)=="0")
#nmres3 <- nmres2[,-(which(colSums(nmres2)==0))]
#ncol(nmres3)
#nrow(nmres3)

haplotypes <- t(nmres)
str(haplotypes)
head(haplotypes)
ncol(haplotypes)

# NMSD analysis, based on a dissimilarity index 
hap_dist <- vegdist(haplotypes, na.rm = T)
hap_nmds <- metaMDS(haplotypes, k=2, distance='bray', try = 30, trymax = 500, engine = "monoMDS") #trymax = 50, 150, 500 did not converged

# Checkpoint to confirm that my data has a linear fit, R2 is close to 1 (?)
stressplot(hap_nmds, hap_dist)

# Cleanup metadata table if needed
head(metadata)
#Sphap_nmds <- metadata   # If upload metadata_nmds2.txt table
Sphap_nmds <- cbind(hap_nmds$species, metadata[which(metadata$strain %in% rownames(hap_nmds$species)),])
head(Sphap_nmds)
nrow(Sphap_nmds)


### PLOT NMDS

shape_values <- seq(1,11) #assigning basic R shapes  

#Sphap_nmds <- metadata #Do this ONLY if input complete table for metadata

mycolyears <- c("brown2", "indianred2", "#D01C8B", "#7B3294", "mediumorchid3", "darkorange1", "#CA0020", "#FBB4AE",
                "cornflowerblue","#A6CEE3", "dodgerblue4", "#B3E2CD", "palegreen3", "darkcyan", "#8DD3C7",
                "#DFC27D", "rosybrown2", "sienna", "sandybrown", "tan", "slategray")

pdf(file = "NMDS_Xp_year_Feb21_2.pdf")
ggplot(Sphap_nmds, aes(shape= as.factor(year))) + 
  geom_text(data= Sphap_nmds, aes(x= MDS1, y= MDS2, label= as.factor(year)), alpha= 0.8) +  # add the species labels
  geom_point(data= Sphap_nmds, 
             aes(x= MDS1, y= MDS2, colour= as.factor(year)), 
             size= 6, 
             alpha = 0.8, 
             shape=20 ) + 
  ggtitle("NMDS of X. perforans")+
  guides( color =  guide_legend(c("BAPS cluster")),
          shape =  guide_legend( "phylogroup"))+
  #scale_colour_manual(values= c(mycolBaps_Sujan, unique(Sphap_nmds$rowcolors))) +  #For BAPS clusters
  #scale_colour_manual(values= c(mycoleff)) +  #For effector profile clusters
  scale_colour_manual(values=c(mycolyears))+ #For years
  coord_equal()+
  theme(panel.background = element_rect(fill = "white", 
                                        colour = "black", 
                                        size = 0.5, 
                                        linetype = "solid"),
        panel.grid.major = element_line(size = 0.5, 
                                        linetype = 'solid',
                                        colour = "gray94"), 
        axis.text=element_text(size=14),
        axis.title=element_text(size=12), 
        plot.title = element_text(hjust = 0.5))
dev.off()

#Save table to file
write.table(Sphap_nmds, file = "metadata_nmds2.txt", col.names = TRUE, row.names = TRUE, sep= "\t", quote = FALSE)




##################

### GET REPRESENTATIVES FROM EFFECTOR PROFILE CLUSTERS 
# From object mat:
head(nmres)

#Defining Effector clusters, from A to J
#metadata[which(metadata$eff_cluster == "J"),]

clust_A <- nmres[which(rownames(nmres) == "GEV2127"),]
clust_B <- nmres[which(rownames(nmres) == "Mexico3"),]
clust_C <- nmres[which(rownames(nmres) == "GEV936D"),]
clust_D <- nmres[which(rownames(nmres) == "ETH33"),]
clust_E <- nmres[which(rownames(nmres) == "THA-127"),]
clust_F <- nmres[which(rownames(nmres) == "SM-1812"),]
clust_G <- nmres[which(rownames(nmres) == "16-1205A"),]
clust_H <- nmres[which(rownames(nmres) == "GEV2125"),]
clust_I <- nmres[which(rownames(nmres) == "TB6"),]
clust_J <- nmres[which(rownames(nmres) == "GEV1063"),]

Effclusters <- rbind(clust_A, clust_B, clust_C, clust_D,clust_E, clust_F,clust_G,clust_H, clust_I, clust_J)
Effcl_col <- c(clust_A, clust_B, clust_C, clust_D, clust_E, clust_F, clust_G, clust_H, clust_I, clust_J)
colnames(Effclusters) <- Effcl_col
Effclusters

#sort(clust_A,decreasing =F)
order_clustA <- order(clust_A, decreasing =F)

order_Effclusters <- t(Effclusters)[order_clustA,]

t(order_Effclusters)
write.table(order_Effclusters, file = "EffProfile_strainclusters.txt", col.names = TRUE, row.names = T, sep= "\t")
getwd()




##################

###  MAP WITH EFFECTOR CLUSTERS
 

