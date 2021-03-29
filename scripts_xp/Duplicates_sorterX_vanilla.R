####################################
# # # DUPLICATES_SORTERX_VANILLA.R
####################################
# Sort and remove duplicated blast hits from Blast files, returning clean list of locus tags to extract sequences
# Fernanda Iruegas Bocardo
# August, 2020
# For SCRI-270 vs 126-effectors analysis
# Derived from script, Parse_blast-tab_Annotated-files_Xp.R



####################################

#########    IMPORT FILES AND LIBRARIES    

####################################

# Input: Blast output files in tabular format (-outfmt 6)
# Blast analysis was performed with 126 effector sequences as query, against a database of 270 genomes (worldwide areas) as subject
# The genomes were assembled with FIB pipeline Asembly_Ipe.sh ("_contigs_spades-qual-pol-rn.fasta") and renamed (file name and contigs) with scripts: update_contig_names.sh/.py by Sujan Timilsina
# Assemblies were annotated in IMG, nucl and prot annotations have same locus_tags
# *This previous information is relevant because of the specific file and contig name format adressed in the script
# Other FILES needed to execute this script:
#  Effector_length.txt
#  Strain list and IMG key 


library(plyr)
library(dplyr)
library(tidyverse)

### ARGUMENTS calling for R vanilla
FILE <- commandArgs(TRUE)
#FILE = c("RipS3,XopAD")
FILE <- unlist(str_split(FILE, pattern=","))

print(FILE)
file=list()
#Blast file name: Xp_270_vs_
for (o in seq_along(FILE)){
  file[[o]] <-  read.table(paste("Xp_270_vs_",FILE[[o]], sep = ""), sep="\t", header=F, stringsAsFactors = F) 
}
#file[[o]]

### EFFECTOR LENGTH table 
# Effector length from each effector will be added into a new column and then used to calculate query coverage, stored into a new column
# Load effector_length.txt table
eff_length <- read.csv(file= "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/01152019_Xp_SCRI_269_vs_125_effectors/125_effector_length_2.txt", 
                       sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)  
# Change the name of the column "Effector" for "queryID", to be able to use it as Key
colnames(eff_length)[which(names(eff_length) == "Effector")] <- "queryID"     
head(eff_length)

### STRAIN list and IMG key 
str_id <- read.csv(file= "strains_key", sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
colnames(str_id) <- c("strain", "goldID")
head(str_id)
nrow(str_id)


####################################

#########    ANALYZE DATA THROUGH FILES IN A LOOP    

####################################

# Read path for input files, and loop through each file

#path = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/200804_Xp_270_vs_125_effectors/"
#out.file <- ""
#file.names <- dir(path, pattern = "Xp_270_vs_")
#length(file.names)

strains <- c(str_id$strain)
hdr <- c("queryID", "subjectID", "percIdentity", "AlnLength", "mismatchCount", "gapOpenCount", "queryStart", 
         "queryEnd", "SubjectStart", "SubjectEnd", "eVal", "bitScore")
length(strains)

#####################################

### Open BLAST files with duplicated locustags (hits from morethan one effector to the same sequence)
blast <- data.frame()
for(n in seq_along(file)){
  temp <- file[[n]]
  blast <- rbind(blast, temp)
  #nrow(blast)
}
#unique(blast$V1)
#str(temp)
colnames(blast) <- hdr

# Split goldID number (IMG) from locus_Tag name and add as a column at the end of the table
tmp_loctag <- c(strsplit(blast$subjectID, split = "_"))
#out <- tmp_loctag[[1]]
out <- c()
for(i in 1:length(tmp_loctag)) {
    out <- rbind(out, tmp_loctag[[i]])
}
blast$goldID <- out[,1]
  
# Add column for effector length, coverage, strain name and weight calculation (to remove duplicates)
blast$queryLength <- eff_length$length[match(blast$queryID, eff_length$queryID)]
blast$queryCov <- (blast$AlnLength / blast$queryLength)*100
blast$strainId <- str_id$strain[match(blast$goldID, str_id$goldID)]
blast$Weight <- (100/as.numeric(blast$queryCov))*as.numeric(blast$percIdentity)   # Add column with weighted value considering Coverage and ID
#head(blast)
#str(blast)
#unique(blast[,1])

# Keep only columns needed and format as character
blast_clean <- data.frame(1:nrow(blast), blast$queryID, blast$strainId, blast$goldID, blast$subjectID, blast$percIdentity, blast$queryCov, blast$Weight)
colnames(blast_clean) <- c("IDabs", "queryId", "strainId", "goldID", "subjectId", "percIdentity", "queryCov", "weight")
#unique(blast_clean[2])
blast_clean$queryId <- as.character(blast_clean$queryId)
blast_clean$strainId <- as.character(blast_clean$strainId)
blast_clean$goldID <- as.character(blast_clean$goldID)
blast_clean$subjectId <- as.character(blast_clean$subjectId)
#head(blast_clean)
#unique(blast_clean[,2])
#str(blast_clean)

# Filter table by % ID, coverage and evalue
blast_clean <- filter(blast_clean, percIdentity >= 70 & queryCov >=50)
#blast_clean <- filter(blast_clean, percIdentity >= 70)
#unique(blast_clean[,1])
#str(blast_clean)

# Get effector name and strains with more than one hit to the same effector (after filtering by id)
#eff_name <- blast[1,1]
eff_name <- unique(blast$queryID)
dups <- blast_clean[duplicated(blast_clean$subjectId)|duplicated(blast_clean$subjectId, fromLast = TRUE),] # Duplicated strains, more than one hit per strain
#unique(dups$queryId)

if(nrow(dups)>0){
  #print(paste(eff_name, length(unique(dups$strainId))))
  #write.table(dups, file = paste("dups_", eff_name, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)
  # Remove duplicates hits from blast results
  uqsi <- unique(dups$subjectId)
  mins <- data.frame()
  #y=3
  for(y in 1:length(uqsi)) {
    aaa <- dups[dups$subjectId==uqsi[y],]
    tmp <- aaa[which(aaa$weight==min(aaa$weight)),]
    mins <- rbind(mins,tmp)
  }  
  #mins
  
  blast_clean <- subset(blast_clean, !(blast_clean$IDabs %in% mins$IDabs))
  #nrow(blast_clean)
  
} else{
  print(paste(eff_name, ", no duplicated hits", sep = ""))}

# select locus_tag and strain name and save to text file, to use for extracting sequences 
#filter(blast_clean, eff_name[[1]][1])
ls <- list()
for(z in seq_along(eff_name)){
  ls[[z]] <- blast_clean[which(blast_clean$queryId==eff_name[z]),] %>% 
    select(subjectId, strainId)
  # locustags <- select(blast_clean, subjectId, strainId)
 
  # Save list of locus_tags as text file 
  write.table(ls[[z]], file = paste(eff_name[z],"_loctag3.txt", sep = ""), col.names = FALSE, row.names = FALSE, sep= "\t", quote = FALSE)
}
  


#### 
# Duplicates found between files:
# RipS3,XopAD
# HopK1,XopAK
# RipAL,XopAP
# CT610,RipAL,XopAP
# AvrXccE1,XopE1,XopE2
# AvrXccE1,XopE2
# XopE2
# XopE1,XopE2
# XopZ1,XopZ

### TEST
#rips <- read.table("RipS3_loctag3.txt", sep="\t", stringsAsFactors = F)
#xopad <- read.table("XopAD_loctag3.txt", sep = "\t",stringsAsFactors = F)
#both <- rbind(rips, xopad)
#str(both)
#nrow(both)
#nrow(rips)
#nrow(xopad)
#table(both$V1)
#both$V1
#length(unique(both$V1))


