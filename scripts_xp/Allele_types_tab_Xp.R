####################################
# ALLELE TYPES TABLE
# Fernanda Iruegas Bocardo
# April, 2019

# Inputs:  haplotypes_parsed.txt, table output from parse_haplotypes.R
#         alleles_type.txt, table output from parse_haplotypes.R
#         125effectors.p2o.cvs, table with Number assigned to effectors, from get_homologs (inputs_homologues/tmp/all.p2o.cvs). Copied the first 125 lines into a new file
#         clusters_summary.txt, Locus_tags organized by strain and effector, from scripts/Analysis/Get_CDS/locus_tags_tab.R
#         out_<eff_number>.faa, aligned multifasta file with sequences from effectors, outfut from MAFFT
#
# Outputs: effector_allele-type.txt, table with with effector name, allele/haplotype type, locus_tag and strain for each hit
#          redundant_hits.txt, coordinates (row, column) for each locus_tag/hit in the file cluster_summary.txt, notice redundant hits (more than one effector)
#          effector_allele-type_sequence.txt, Final table with effector name, allele/haplotype type, locus_tag, strain and sequence from MAFFT alignment
#


# Objective:
# Create table with effector name, strain, locus_tag, allele type number and aa sequence

# Pseudocode:
# 1. Get locus tags from object: all_hap
# 2. Assign allele type and effector number to each locus_tag from object: hap_type 
# 3. Get strain name for each locus_tag from object: loc_tag
# 4. Get effector name and replace effector number from object: eff_num
# 5. Get aligned sequence for each locus_tag and include in the table

####################################

###    IMPORT FILES AND LIBRARIES    


library(plyr)          # To concatenate dataframes from list
library(dplyr)         # **Double check if I need this library
library(tidyr)
library(Biostrings)    # To import multifasta files and manipulate as dataframe

## haplotypes_parsed.txt
# This .txt file is the output of parse_haplotypes.R
# Table contains the haplotypes for 70 effectors, organized by columns
# Columns contain locus tags for one haplotype and are named after the effector number


all_hap <- loc_tags # From script parse_haplotypes.R
# OR
all_hap <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes/haplotypes_parsed_singletons.txt", 
                   sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE, strip.white = T)
hap_hdr <- all_hap[1,]
colnames(all_hap) <- hap_hdr
all_hap <- all_hap[-1,]
nrow(all_hap)
head(all_hap)

# Optional!To remove (ref) 
all_hap_clean <-  all_hap
for(h in 1:ncol(all_hap_clean)){
  clean_col <- all_hap_clean[,h]
  clean_col <- gsub("\\s*\\([^\\)]+\\)", "", clean_col)
  all_hap_clean[,h] <- clean_col
}
all_hap <- all_hap_clean
head(all_hap)


## alleles_type.txt
# This .txt file is the output of parse_haplotypes.R
# Columns:
# e: effector
# p: position of column for the haplotype in file haplotypes_parsed.txt
# f: frequency or number of strains (nrow) with a type haplotype
# a: allele type, after sorting the alleles by frequency. 1 is the most abundant.

hap_type <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes/alleles_type_singletons.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE, strip.white = T)
# or
#hap_type <- alleles
head(hap_type)
  
## 125effectors.p2o.cvs
# Number assigned to effectors
# Created by get_homologs, as inputs_homologues/tmp/all.p2o.cvs
# The file is too large, I copy the first 125 lines into a new file, uploaded here
eff_num <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xe_All_analysis_070519/eff_number.txt", 
                       sep=",", stringsAsFactors=FALSE, fill = TRUE, header = FALSE, strip.white = T)
eff_hdr <- c("num", "gh_input", "effector")
colnames(eff_num) <- eff_hdr
eff_num

## OR
### STRAIN list and IMG key 
str_id <- read.csv(file= "strains_key", sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
colnames(str_id) <- c("strain", "goldID")
head(str_id)
nrow(str_id)


## clusters_summary.txt
# Locus_tags organized by strain and effector
# Created by scripts/Analysis/Get_CDS/locus_tags_tab.R

loc_tag <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/clusters_summary2.txt", 
                      sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE, strip.white = TRUE)

col_hdr <- loc_tag[1,]
row_hdr <- loc_tag[,1]
colnames(loc_tag) <- col_hdr
rownames(loc_tag) <- row_hdr
loc_tag <- loc_tag[-1,-1]
head(loc_tag)

## Alignned sequences from each locus tag
# Aligned multifasta files, in fasta format (.faa)
# Need to move working directory to where I have the output form MAFFT
# Do after organizing the rest of the information in the table




####################################

###    ANALYZE DATA THROUGH FILES   



### 1. and 2. Make a dataframe with effector name and haplotype type, for each locus tag

ncol(all_hap)
nrow(hap_type)
#i=1
#x=1
list <- list()
for(i in 1:nrow(hap_type)){
  eff.temp <- hap_type[i,"e"]
  hap.temp <- hap_type[i,"a"]
  pos.temp <- hap_type[i,"p"]
  col.temp <- as.character(all_hap[,pos.temp])
  col.temp[col.temp==""] <- NA
  col.temp <- na.omit(col.temp)
  df.temp <- data.frame(eff.temp, hap.temp, col.temp)
  list[i] <- list(df.temp)
}
#list
head(list)
str(list)

# Concatenate dataframes from list into a single dataframe
# Dataframe with effector, allele, locus_tag 
eal_df <- ldply(list, rbind) 
head(eal_df)
class(eal_df)
nrow(eal_df)
str(eal_df)
eal_df$col.temp <- as.character(eal_df$col.temp)
head(eal_df)


### 3. Get strain name for each locus_tag from the clusters_summary table
# Create dataframe with all the coordinates for each locus_tags from (eal_df), from the loc_tag object

# str(eal_df)
# x=1
# nrow(loc_tag)
str(loc_tag)
strains <- c()
redundant_hits <- c()
redundant_list <- list()
for(x in 1:nrow(eal_df)){
  #loc_temp <-  na.omit(as.character(eal_df$col.temp[x]))
  loc_temp <- gsub("\\s*\\([^\\)]+\\)", "", na.omit(as.character(eal_df$col.temp[x])))
  loc_coord <- data.frame(which(as.character(loc_tag) == as.character(loc_temp), arr.ind = TRUE))
  str.tmp <- rownames(which(loc_tag == loc_temp, arr.ind = TRUE))
  #redundant_list[x] <- list(loc_coord)
  if(length(str.tmp) == 1){
    strains[x] <- str.tmp
  }else{
    strains[x] <- str.tmp[1]
    #redundant_hits <- c(redundant_hits, loc_temp, str.tmp[2:length(str.tmp)])
  }
}
strains
length(strains)
redundant_hits
length(redundant_hits)
eal_df$strain <- strains

### OR 
### Use strains_key file
head(eal_df)
{
  tmp_loctag <- c(strsplit(eal_df$col.temp, split = "_"))
  out <- tmp_loctag[[1]]
  for(i in 2:length(tmp_loctag)) {
    out <- rbind(out, tmp_loctag[[i]])
  }
  eal_df$goldID <- out[,1]
  head(eal_df) 
}

# Add strain column to eal_df table 
eal_df$strainID <- str_id$strain[match(eal_df$goldID, str_id$goldID)]    
head(eal_df)
nrow(eal_df)
#eal_df$strainID[which(is.na(eal_df$strainID))] <- "GEV1063"

# Dataframe with effector, allele, locus_tag, strain 
eals_df <- eal_df

# List of coordinates for each locus_tag, or blast hit, to check for redundant hits
# redundant_list to redundant_dataframe
redundant_df <- ldply(redundant_list, rbind) #dataframe with effector, allele, locus_tag 
head(redundant_list)
nrow(redundant_df)

# - - - - - - - - - - - - - - -
### Alternative 3. Other method to include strain name into eal_df object
library(tidyr)
str(eal_df)
colnames(eal_df) <- c("eff.temp", "hap.temp", "locus_tag")
dat <- cbind(strains=rownames(loc_tag), loc_tag)
dat2 <- gather(dat, 
               key = "eff.temp", 
               value = "locus_tag", 
               -strains)
dat2$strains <- as.character(dat2$strains)
str(dat2)
dat3 <- dat2[!is.na(dat2$locus_tag),]
aver <- dat3[order(dat3$locus_tag) %in% order(eal_df$locus_tag),]
# aver$eff.temp <- as.numeric(aver$eff.temp)
head(aver)
head(eal_df)

eal_df2 <- cbind(eal_df[order(eal_df$locus_tag),], aver$strains[order(aver$locus_tag)])
eal_df2$id <- as.numeric(rownames(eal_df2))
eal_df2 <- eal_df2[order(eal_df2$id),]
eal_df2 <- eal_df2[-5]
write.table(eal_df2, "eal_df.txt", col.names = TRUE, row.names = FALSE, sep= "\t")
head(eal_df2, 40)

# Dataframe with effector, allele, locus_tag, strain 
eals_df <- eal_df2

# - - - - - - - - - - - - - - -

### 4. Change name of effectors from numbers to the effector name, from eff_num object

str(eals_df)

#y=1
eff.nm <- c()
for(y in 1:nrow(eals_df)){
  eff.tmp.nm <- eals_df$eff.temp[y]
  eff.tmp.nm <- eff_num[which(eff.tmp.nm==eff_num$num),3]
  eff.nm <- c(eff.nm, eff.tmp.nm)
}

#Add eff.nm column with the names of the effectors to the table eals_df
length(eff.nm)
eff.nm
eals_df$eff.temp <- eff.nm
eals_df
final_hdr <- c("effector", "haplotype", "locus_tag", "goldID", "strain")
colnames(eals_df) <- final_hdr
str(eals_df)
tail(eals_df)

### Save table to file
#Table with effector name, allele/haplotype type, locus_tag and strain for each hit
write.table(eals_df[,c(1:5)], file = "effector_allele-type_singletons2.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)

#Coordinates from cluster_summary.txt for all the hits, notice redundant hits
write.table(redundant_df, file = "redundant_hits_singletons.txt", col.names = TRUE, row.names = FALSE, sep= "\t")



### 5. Get aligned sequence for each locus_tag and add to the table in a new column
# Import multifasta files, with sequence and locus tag and incorporate in table
# Change working directory to where the alignment files are (out_<effector_number>.faa)
#library(Biostrings) #Loaded at the begining of the script

# Read file names from aligned multifasta files
path = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/alignments/"
out.file <- ""
file.names <- dir(path, pattern = "out_")

# Create list of dataframes from multifasta files, including locus_tag and sequence in different columns
lst.fasta <- list()
for(z in 1:length(file.names)){
  print(paste("Opening", file.names[z], sep = " "))
  
  #if((file.size(file.names[z]) != 0)==TRUE){
    fastaFile <- readDNAStringSet(file.names[z])
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)
    lst.fasta[z] <- list(df)
  #}else{
  #  next()
  #}
}
lst.fasta
head(lst.fasta)

# Eliminate null objects from list
lst.fasta <- lst.fasta[lengths(lst.fasta)>0L]
str(lst.fasta)
length(lst.fasta)

# Create single dataframe from list of dataframes with locus_tags and sequences
seq_df <- ldply(lst.fasta, rbind)
head(seq_df)
#write.csv(seq_df, "table_sequences.csv")
#seq_df <- read.csv("table_sequences.csv", as.is = T, stringsAsFactors = F)


# Match locus_tags from seq_df and eals_Df, and add sequence corresponding locus tag in a new column in table eals_df
# If using object eal_df or eals_df, *notice that they have different column names, modify if necessary

head(eals_df)
head(seq_df)
str(eals_df) #Notice column names and modify code if necessary
str(seq_df)
seq_df$seq_name <- as.character(seq_df$seq_name)
seq_df$sequence <- as.character(seq_df$sequence)
nrow(eals_df)
nrow(seq_df)

#s=1
seq <- c() 
eals_df$strain <- as.character(eals_df$strain)
for(s in 1:nrow(eals_df)){
  loc_temp <- gsub("\\s*\\([^\\)]+\\)", "", as.character(eals_df$locus_tag[s]))
  seq.row <- which(seq_df$seq_name %in% loc_temp)
  seq.tmp <- as.vector(seq_df[seq.row,2])
  # Some locus_tags have hits to more than one effector 
  # that information is already organized in table reduntant_df
  # here I ignore this and only use the first element
  # if(length(seq.tmp) == 1){
  seq <- c(seq, seq.tmp)
  # }else{
    # seq <- c(seq, seq.tmp[1])
  # }
}

length(seq)
eals_df$sequence <- seq
tail(eals_df)


### 
# Alternative 5.!!! Second part (after creating seq_df object): Add sequences to eals_df table     # <<== FINAL TABLE
eals_df$locus_tag <- sub("\\s*\\([^\\)]+\\)", "", eals_df$locus_tag) # remove extra space from locus tag
seq_df$seq_name %in% eals_df$locus_tag # check both names have equal names
colnames(seq_df) <- c("locus_tag","sequences") #match locus tag in both vectors
new_tabla <- merge(eals_df, seq_df, by = ("locus_tag") ) # comando merge mergea ambas columnas by "name"
dim(new_tabla)

head(new_tabla)
eals_df_final <- new_tabla



### WRITE FINAL TABLE
# Final table with effector name, allele/haplotype type, locus_tag, strain and sequence
# **Notice the working diretory where the file will be saved
write.table(eals_df, file = "effector_allele-type_sequence_singletons.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)












####################################

### TESTING

### 3. Alternative!!! to using the above loop (in case it does not work!)
locustags_alt <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/All_analysis_062719/haplotypes/strains_key.txt", 
                            sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE, strip.white = T)
str(eal_df)

strains_alt <- rownames(loc_tag)
gsub("_[0-9)])", "", locustags_alt)
grep()
grep(locustags_alt[y,], loc_temp)
substr(loc_temp,0, 9)
loc_coord <- as.character(which(locustags_alt == loc_temp, arr.ind = TRUE))
x=2
y=1

strains <- c()
for(x in 1:nrow(eal_df)){
  if(length(str.tmp) != 0){
    loc_temp <- as.character(na.omit(eal_df$col.temp[x]))
    loc_temp <-  substr(loc_temp, 1, 9)
    loc_coord <- as.character(which(locustags_alt$V2 == as.character(loc_temp)))
    str.tmp <- locustags_alt[loc_coord,1]
    strains[x] <- str.tmp
  }else{
    loc_temp <-  substr(loc_temp, 1, 10)
    loc_coord <- as.character(which(locustags_alt$V2 == loc_temp))
    str.tmp <- locustags_alt[loc_coord,1]
    strains[x] <- str.tmp[1]
  }
}
strains
length(strains)
redundant_hits
length(redundant_hits)
eal_df$strain <- strains


### Alternative 5. Add sequence to eals_df 
#Not working:
# Add seq vector to new column named sequence in eals_df table
eals_df_order <- eals_df[order(eals_df$locus_tag),] 
eals_df_order$sequence <- seq
seq
seq_df_order <- seq_df[order(seq_df$seq_name),]
head(seq_df_order)
head(eals_df_order)
cbind(as.character(seq_df_order$seq_name), as.character(eals_df_order$locus_tag))
which(eals_df_order$locus_tag %in% seq_df_order$seq_name)
seq # <<== FINAL TABLE
head(eals_df)
length(unique(eals_df$effector))



