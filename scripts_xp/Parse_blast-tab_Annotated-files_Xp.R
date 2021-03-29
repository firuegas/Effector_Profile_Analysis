####################################
# Blast parsing by ID percent for tabular format
# Fernanda Iruegas Bocardo
# June, 2018
# Edited January 2019, for SCRI-269 vs 125-effectors analysis
# Edited July 2020, for SCRI-270 vs 125-effectors analysis


# Blast output files in tabular format (-outfmt 6)
# Blast analysis was performed with 125 effector sequences as query, against a database of 270 genomes (worldwide areas) as subject
# The genomes were assembled with FIB pipeline Asembly_Ipe.sh ("_contigs_spades-qual-pol-rn.fasta") and renamed (file name and contigs) with scripts: update_contig_names.sh/.py by Sujan Timilsina
# Assemblies were annotated in IMG, nucl and prot annotations have same locus_tags
# *This previous information is relevant because of the specific file and contig name format adressed in the script
# Other FILES needed to execute this script:
#  Effector_length.txt
#  Genome_db_metadata.txt     ##In progress, for color labels

####################################

#########    IMPORT FILES AND LIBRARIES    #########  

####################################

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)


### EFFECTOR LENGTH table 
# Effector length from each effector will be added into a new column and then used to calculate query coverage, stored into a new column

# Load effector_length.txt table
eff_length <- read.csv(file= "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/01152019_Xp_SCRI_269_vs_125_effectors/125_effector_length_2.txt", sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = TRUE)  

# Change the name of the column "Effector" for "queryID", to be able to use it as Key
colnames(eff_length)[which(names(eff_length) == "Effector")] <- "queryID"     
head(eff_length)


### STRAIN list and IMG key 
str_id <- read.csv(file= "strains_key", sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
colnames(str_id) <- c("strain", "goldID")
head(str_id)
nrow(str_id)


### Input single Blast file for Test
blast <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/200804_Xp_270_vs_125_effectors/GEV1063/gev1063_vs_125effectors", 
                    sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
hdr <- c("subjectID", "queryID", "percIdentity", "AlnLength", "mismatchCount", "gapOpenCount", "queryStart", 
         "queryEnd", "SubjectStart", "SubjectEnd", "eVal", "bitScore")
colnames(blast) <- hdr
head(blast)    


####################################

#########    Analysis of single Blast file

####################################

{
  # Split goldID number from locus_Tag name and add as a column at the end of the table
  {
    tmp_loctag <- c(strsplit(blast$subjectID, split = "_"))
    out <- tmp_loctag[[1]]
    for(i in 2:length(tmp_loctag)) {
      out <- rbind(out, tmp_loctag[[i]])
    }
    blast$goldID <- out[,1]
    #head(blast) 
  }
  
  #Add a column with effector length by indexing the information from the eff_length object
  #match() returns a vector of position, then I can use the position to  index[] the length column, outside the match function
  blast$queryLength <- eff_length$length[match(blast$queryID, eff_length$queryID)]
  
  # Calculate and add column with coverage
  blast$queryCov <- (blast$AlnLength / blast$queryLength)*100   
  
  # Add strain column to blast table 
  blast$strainID <- str_id$strain[match(blast$goldID, str_id$goldID)]    
  
  # Add column with weighted value considering Coverage and ID
  blast$Weight <- (100/as.numeric(blast$queryCov))*as.numeric(blast$percIdentity)
  head(blast)

    # Clean table
  blast_clean <- data.frame(blast$queryID, blast$strainID, blast$goldID, blast$subjectID, blast$percIdentity, blast$queryCov, blast$Weight)
  colnames(blast_clean) <- c("query", "strainId", "goldID", "subjectId", "percIdentity", "queryCov", "weight")
  head(blast_clean)
  
  ### Filter table by % ID, coverage and evalue
  # Comment the filters in order to identify partial hits (duplicated)
  blast_clean <- filter(blast_clean, percIdentity >= 70 & queryCov >=50)
  #blast_clean <- filter(blast_clean, percIdentity >= 70)
  
}

head(blast_clean)
blast_clean[order(blast_clean$strainId),]                 #order by strain name
nrow(blast_clean)                                         #number of hits for x effector
length(unique(blast_clean$strainId))                      #number of strains with hits
blast_clean$strainId[duplicated(blast_clean$strainId)]    #strains with duplicated hits

####################################

#########    ANALYZE DATA THROUGH FILES IN A LOOP    

####################################

# Read path for input files, and loop through each file

path = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Blast_results/200804_Xp_270_vs_125_effectors/"
out.file <- ""
file.names <- dir(path, pattern = "Xp_270_vs_")
length(file.names)

strains <- c(str_id$strain)
hdr <- c("queryID", "subjectID", "percIdentity", "AlnLength", "mismatchCount", "gapOpenCount", "queryStart", 
         "queryEnd", "SubjectStart", "SubjectEnd", "eVal", "bitScore")
length(strains)

####################################

# Gets 2 vector objects with lists of hits and no hits 
# The FOR loop uses an if() conditional statement that assess the size of the file 
# If the size of the file is 0 bits, then the file is empty meaning that there was no hit
# Then, replaces the values from the "file.names" vector with the content of "hits" vector

hits <- c()
no_hits <- c()

for(i in 1:length(file.names)){
  print(paste("Opening", file.names[i], sep = " "))
  
  if(file.size(file.names[i]) == 0){
    no_hits <- rbind(no_hits, file.names[i])
    hits <- hits
  }else{
    no_hits <- no_hits
    hits <- rbind(hits, file.names[i])
  }
}

# Updated list of files to analyze
file.names <- as.vector(hits)
length(file.names)

#i=13
list_blasts <- list()
for(i in 1:length(file.names)){
  # Open single file
  blast <- read.csv(file= file.names[i], sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
  colnames(blast) <- hdr 
  #nrow(blast)
  
  # Split goldID number (IMG) from locus_Tag name and add as a column at the end of the table
  tmp_loctag <- c(strsplit(blast$subjectID, split = "_"))
  #out <- tmp_loctag[[1]]
  out <- c()
  for(x in 1:length(tmp_loctag)) {
    out <- rbind(out, tmp_loctag[[x]])
  }
  blast$goldID <- out[,1]
  
  # Add column for effector length, coverage, strain name and weight calculation (to remove duplicates)
  blast$queryLength <- eff_length$length[match(blast$queryID, eff_length$queryID)]
  blast$queryCov <- (blast$AlnLength / blast$queryLength)*100
  blast$strainId <- str_id$strain[match(blast$goldID, str_id$goldID)]
  blast$Weight <- (100/as.numeric(blast$queryCov))*as.numeric(blast$percIdentity)   # Add column with weighted value considering Coverage and ID
  #head(blast)
  
  # Keep only columns needed
  blast_clean <- data.frame(blast$queryID, "GEV1063", blast$goldID, blast$subjectID, blast$percIdentity, blast$queryCov, blast$Weight)
  colnames(blast_clean) <- c("queryId", "strainId", "goldID", "subjectId", "percIdentity", "queryCov", "weight")
  head(blast_clean)
  blast_clean$queryId <- as.character(blast_clean$queryId)
  blast_clean$strainId <- as.character(blast_clean$strainId)
  blast_clean$goldID <- as.character(blast_clean$goldID)
  blast_clean$subjectId <- as.character(blast_clean$subjectId)
  
  # Filter table by % ID, coverage and evalue
  blast_clean <- filter(blast_clean, percIdentity >= 70 & queryCov >=50)
  #blast_clean <- filter(blast_clean, percIdentity >= 70)
  
  # Get effector name and strains with more than one hit to the same effector (after filtering by id)
  eff_name <- blast[1,1]
  dups <- blast_clean[duplicated(blast_clean$strainId)|duplicated(blast_clean$strainId, fromLast = TRUE),] # Duplicated strains, more than one hit per strain
  
  if(nrow(dups)>0){
    print(paste(eff_name, length(unique(dups$strainId))))
    write.table(dups, file = paste("dups_", eff_name, ".txt", sep = ""), col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)
  
    # Remove duplicates hits from blast results
    uqst <- unique(dups$strainId)
    mins <- data.frame()
    for(y in 1:length(uqst)) {
      aaa <- dups[dups$strainId==uqst[y],]
      tmp <- aaa[which(aaa$weight==min(aaa$weight)),]
      mins <- rbind(mins,tmp)
    }  
  #mins
  
  blast_clean <- subset(blast_clean, !(blast_clean$subjectId %in% mins$subjectId))
  #nrow(blast_clean)
  
  } else{
    print(paste(eff_name, ", no duplicated hits", sep = ""))}
  
  # select locus_tag and strain name and save to text file, to use for extracting sequences 
  locustags <- select(blast_clean, subjectId, strainId)
  #head(locustags)
  #nrow(locustags)
  
  #To save results in a list of lists
  #Comment the following code within the first loop, to save lists withouth writing the files
  list_blasts[i] <- list(blast_clean)
  
  # Save list of locus_tags as text file
  write.table(locustags, file = paste(eff_name,"_loctag2.txt", sep = ""), col.names = FALSE, row.names = FALSE, sep= "\t", quote = FALSE)
  
}

#### 

# - - - - - - - - -
# Make summary tables with clean results from BLAST in long and wide format (only locustags)

# Format lists and dataframe with all results, if it applies
list_blasts

# Remove empty objects from list
list_blasts <- list_blasts[lengths(list_blasts)>0L]
length(list_blasts)
#list_blasts

# Create dataframe from all lists
df_lists <- ldply(list_blasts, rbind)
head(df_lists)
nrow(df_lists)
length(unique(df_lists$subjectId))
df_lists_clean <- df_lists[,c(1,2,4)] #Clean dataframe with effector, strain and locustags only
head(df_lists_clean)

# Make WIDE format table, from LONG table with clean blast results
loctags_summary <- spread(df_lists_clean, queryId, subjectId)     #! this table still have duplicated locus_tags

# Save both tables to file
write.table(blast, file = "all_blast_clean.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)
write.table(loctags_summary, file = "clusters_summary.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE) 

