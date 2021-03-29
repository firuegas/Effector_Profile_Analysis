####################################
# Parse the allele or haplotype types from the log_pipe output from haplotype.pl
# Fernanda Iruegas Bocardo
# April, 2019

# Input:  tab delimited file with the log_pipe (ErrOut) information from the script haplotype.pl
#         This files contain the locus_tags from the effector sequences organized in columns:
#         The first locus_tag is recorded as (ref) or reference, the rest of the locus tags on the colum are identical to the reference (0) or with the no. of different aa (-/+ \d+)
# Outputs: haplotypes_parsed.txt, table with the locus tags organized by effectors (numeric name)
#          alleles_type.txt, table with the effector name(numeric), frequency, column position and allele type named numericly in descending order according to frequency


####################################

###    IMPORT FILES AND LIBRARIES    


library(plyr)          # To concatenate dataframes from list
library(dplyr)         # **Double check if I need this library



# This .txt file is the log_pipe_xxx file automatically outputed from running by sbatch the scripthaplotypes.pl
# In textwrangler:  " = " was substituted by "\t", and saved as log_pipe_xxx.ed
# In excel: table was transposed, and saved as log_pipe_xxx.ed.transpose.txt
# ** This previous steps can be done in R, maybe in future versions if necessary **

file <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes/haplotypes_with_gev1063/log_pipe_58023499_ed_transposed.txt", 
                      sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
head(file)


####################################

###    PARSE TABLE


### Parse effector names and locus tags into table 
# Find the positions of effector name and columns with locus_tags from alleles
opening <- c()
allele <- c()

for(i in 1:ncol(file)){
  if(!is.na(pmatch("opening", file[1,i])) == TRUE){
    opening[i] <- i
  }else {
    print("doesn't start with opening")
    #opening[i] <- NA
  }
  if(!is.na(pmatch("G", file[1,i])) == TRUE){
    allele[i] <-  i
  }else {
    print("doesn't start with G...")
    #allele[i] <-  NA
  }
  if(!is.na(pmatch("K", file[1,i])) == TRUE){
    allele[i] <-  i
  }else {
    print("doesn't start with K ...")
    #allele[i] <-  NA
  }
}

opening
allele

# Eliminating NAs from vectors
opening_clean <- opening[!is.na(opening)]
opening_clean
allele_clean <- allele[!is.na(allele)]
allele_clean
 

#########


### Organize locus_tags columns with the effector 
#i=4
#x=2

temp.col <- matrix()
loc_tags <- c()
effectors <- c()
for(i in 1:length(opening_clean)){
  eff.name <- file[1,opening_clean[i]]
  eff <- sub("opening ../alignments/out_*", "", eff.name)
  eff <- sub(".faa", "", eff)
  
  for(x in 1:length(allele_clean)){
      if( is.na(opening_clean[i+1]) == FALSE){
        if((isTRUE(opening_clean[i] < allele_clean[x]) &&  isTRUE(allele_clean[x] < opening_clean[i+1])) == TRUE){
          temp.col <- as.matrix(file[,allele_clean[x]])
        }  else{
          next()
        }
      }
      if(is.na(opening_clean[i+1]) == TRUE){
          if(isTRUE(opening_clean[i] < allele_clean[x]) == TRUE){
          temp.col <- as.matrix(file[,allele_clean[x]])
        }  else{
          next()
        }
      }
    loc_tags <- cbind(loc_tags,temp.col)
    effectors <- c(effectors,eff)
    }
}

#temp.col <- temp.col[grep("Ga", temp.col)]
temp.col
loc_tags
effectors

# Change name if needed
#effectors <- sub("103-51*", "103", effectors)
#effectors <- sub("124-67*", "67", effectors)
#effectors <- sub("99-36*", "99", effectors)
#effectors <- sub("93-54*", "93", effectors)
#effectors <- sub("80-8*", "8", effectors)
#effectors

str(loc_tags)

# Name colums of locus tags with the effector name 
#effectors <- as.vector(unlist(file[1,]))   #If importing directly from file haplotypes_parsed_singletons
#loc_tags <- file[-1,]    #If importing directly from file haplotypes_parsed_singletons
colnames(loc_tags) <- effectors
loc_tags


length(effectors)
ncol(loc_tags)
loc_tags


#Save table
write.table(loc_tags, file = "haplotypes_parsed.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)

# OR
#Open from file
loc_tags <- read.table(file = "/Users/FerIruegas/Documents/Documents_MacBook_Pro/UF_PhD_Xanthomonas/Sequences/Analysis_Sequence_MAFFT_RAxML/Xp_effProfile_2020/haplotypes_with_gev1063/haplotypes_parsed_singletons2.txt", 
                      sep="\t", stringsAsFactors=FALSE, fill = TRUE, header = FALSE)
nrow(loc_tags)
ncol(loc_tags)
head(loc_tags)
colnames(loc_tags) <- loc_tags[1,]
loc_tags <- loc_tags[-1,]
head(loc_tags)
#effectors <- colnames(loc_tags)
#u.effectors <- unique(effectors)

####################################

###    PARSE HAPLOTYPES

u.effectors <- unique(effectors)
length(u.effectors)
#loc_tags[which(colnames(loc_tags) == 103)]
#which(colnames(loc_tags) == 93)

#Create a list with the positions [2] for every effector [1]
#lst.lst <- list()
pos <- c()
for(i in 1:length(u.effectors)){
  df.eff <- u.effectors[i]
  col.eff <- (which(colnames(loc_tags) == df.eff))
  #lst.eff <- list(effector=df.eff, position= col.eff)   #Probably not needed
  #lst.lst <- list(lst.lst, lst.eff)                     #Probably not needed
  pos <- c(pos, col.eff)
}
#lst.eff 
#lst.lst
pos

#count frequency of allele, or row lenght per column
freq.all <- c()
for(i in 1:ncol(loc_tags)){
  freq <- sum((loc_tags[,i] != ""))
  freq.all <- c(freq.all,freq)
}
freq.all

sum(freq.all)

# Data frame with the effector name, its frequency and column position in the dataframe loc_tags

ef.fq.ps <- data.frame(effectors, freq.all, pos)

list.efas <- list()
name <- c()
for(i in 1:length(u.effectors)){
  
  which(ef.fq.ps$effectors==u.effectors[i])
  e <- ef.fq.ps$effectors[which(ef.fq.ps$effectors==u.effectors[i])]    #e: effector
  p <- which(ef.fq.ps$effectors==u.effectors[i])      #p:  position / column
  f <- ef.fq.ps$freq.all[p]                           #f: frequency
  a <- (1:length(p))                                  #a: allele number 
  ef <- data.frame(e,p,f)
  ef.s <- ef[order(-ef$f),]
  efa.s <- data.frame(ef.s,a)                       
  nm <- as.vector(u.effectors[i])                     #name
  name <- c(name, nm)
  list.efas[i] <- list(efa.s)

}

#unlist(list.efas)
e
f
a
ef
ef.s
efa.s
list.efas


# **Optional** Write out single data frames (Creates an object for each dataframe)
for (i in 1:length(list.efas)){
  assign(paste0("alleles_", name[i]), list.efas[[i]])
}

# Concatenate dataframes from list into a single dataframe <<<== FINAL TABLE!!!
alleles <- ldply(list.efas, rbind)
head(alleles)
sum(alleles$f)
tail(alleles)

#Save table
write.table(alleles, file = "alleles_type_singletons2.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)


### Effector total frequencies table
length(u.effectors)
u.effectors
alleles

# To get frequencies
freqs <- c()
for(h in 1:length(u.effectors)){
  each_eff <- alleles[which(alleles$e == u.effectors[h]),]
  frq <- sum(each_eff$f)
  freqs[h] <- frq
}
freqs

# To get no. alleles
no_alleles <- c()
for(a in 1:length(u.effectors)){
  each_eff <- alleles[which(alleles$e == u.effectors[a]),]
  alle <- max(each_eff$a)
  no_alleles[a] <- alle
}
no_alleles

# Create table
eff_freqs <- data.frame(u.effectors, freqs, no_alleles)

# Change name of effectors from numbers to the effector name, from eff_num object

str(eff_freqs)
#eff_freqs$u.effectors <- sub("-[0-9]*", "", eff_freqs$u.effectors)
#y=1
eff.nm <- c()
for(y in 1:nrow(eff_freqs)){
  eff.tmp.nm <- eff_freqs$u.effectors[y]
  eff.tmp.nm <- eff_num[which(eff.tmp.nm==eff_num$num),3]
  eff.nm <- c(eff.nm, eff.tmp.nm)
}

#Add eff.nm column with the names of the effectors to the table eals_df
length(eff.nm)
eff.nm
eff_freqs$e.name <- eff.nm
eff_freqs
final_hdr <- c("e", "f", "a", "e.name")
colnames(eff_freqs) <- final_hdr
str(eff_freqs)


### Save table to file
#Table with effector name, allele/haplotype type, locus_tag and strain for each hit
write.table(eff_freqs, file = "eff_freq_singletons2.txt", col.names = TRUE, row.names = FALSE, sep= "\t")



### - - - - - - - - - - Find missing sequences, and repeat first part of the code, to analyze with singletons included, change name of files to save

### Finding missing sequences: 9582 as initial input to the haplotype.pl script, and outputted 9469
# To include singletons from hap_clean_out...faa in the analysis 

# First, Identify which locus tags from the hap_clean<...>.faa files are not included in the log_pipe_<...> file (the singletons were not added to the haplotype analysis)

# Take the first line of loc_tags object, or haplotypes_parsed.txt file
# (loc_tags[1,] == hap_reference) #haplotypes_parsed.txt
hap_reference <- as.vector(loc_tags[1,])
hap_reference <- gsub("\\s*\\([^\\)]+\\)", "", hap_reference)
#which(hap_reference=='Ga0307256_104726')
length(hap_reference)

# locus tags from the sequences on the hap_out_xx.faa files, output of haplotypes.pl
hap_sequences <- names(all_clust)
hap_sequences <-  as.vector(c("Ga0373883_043_40091_40711", "Ga0307235_104410", "Ga0307220_100392", "Ga0307227_101093", "Ga0248036_11379", 
                     "Ga0248175_10280", "Ga0307206_100397", "Ga0307198_10281", "Ga0307256_104726", "Ga0307235_110323", 
                     "Ga0248044_1444", "Ga0307235_101139", "Ga0248175_11266", "Ga0136664_100312", "Ga0248088_10414", "Ga0307197_101254", 
                     "Ga0307185_103109", "Ga0077199_126393", "Ga0307235_105110", "Ga0307229_100472", "Ga0307155_103195", "Ga0307196_13442", 
                     "Ga0307235_104517", "Ga0248044_10963", "Ga0248072_108114", "Ga0248073_105415", "Ga0307180_105143", 
                     "Ga0136666_111516", "Ga0136678_113712", "Ga0307154_106143", "Ga0307196_103137", "Ga0248044_1445", 
                     "Ga0307235_100298", "Ga0307227_100413", "Ga0077200_123237", "Ga0307235_110213", "Ga0138109_111102", 
                     "Ga0248044_105176", "Ga0248095_103176", "Ga0307153_100387", "Ga0307235_10356", "Ga0307220_100316", 
                     "Ga0307229_12591", "Ga0307224_10121", "Ga0248124_1001157", "Ga0248080_101398", "Ga0136467_11561", 
                     "Ga0136663_11264", "Ga0307211_10172", "Ga0307214_11226", "Ga0248088_100210", "Ga0307153_100522", 
                     "Ga0307143_101157", "Ga0307196_10295", "Ga0307197_102217", "Ga0373888_04_168276_170603", "Ga0248072_108112", 
                     "Ga0248084_110618", "Ga0307235_10112", "Ga0373888_02_104944_107088", "Ga0373883_006_102612_104312", "Ga0248176_111103", "Ga0373881_008_129393_131537", 
                     "Ga0248124_100490", "Ga0307184_11299", "Ga0136465_111210", "Ga0373886_20_1323_3023", "Ga0307201_101167", "Ga0077649_10591", "Ga0077199_126355", "Ga0077201_119618", 
                     "Ga0373888_25_30516_31985", "Ga0307229_10584", "Ga0248122_101545", "Ga0307204_11464", "Ga0248088_10374", "Ga0248091_103945", 
                     "Ga0307235_101336", "Ga0307220_102168", "Ga0248072_102282", "Ga0248122_101843", "Ga0307235_101933", "Ga0307220_101332", "Ga0248044_10531", "Ga0307177_10232", 
                     "Ga0248079_102532", "Ga0248084_111313", "Ga0373887_22_31118_32104", "Ga0307242_10432", "Ga0307235_104334", "Ga0248044_10441", "Ga0136678_11151", "Ga0248088_104535", "Ga0307201_11276", 
                     "Ga0307202_11776", "Ga0307235_102243", "Ga0307226_105929", "Ga0307227_102915", "Ga0307232_10246", "Ga0248073_101343", 
                     "Ga0307184_12344", "Ga0136682_12641", "Ga0339969_102144", "Ga0307256_103515", "Ga0307235_103028", "Ga0256589_10681", "Ga0307235_101533", 
                     "Ga0307235_106210", "Ga0307220_101039", "Ga0307229_103718", "Ga0307226_100946", "Ga0248036_10713", "Ga0136667_114114", "Ga0136674_102818", "Ga0248084_104231", "Ga0248088_103139", "Ga0248085_13338", "Ga0248087_13239", 
                     "Ga0307152_13118", "Ga0307188_101216", "Ga0256605_10335", "Ga0307235_10352", "Ga0307220_100312", "Ga0307227_101012", "Ga0248073_105121", 
                     "Ga0248175_102163", "Ga0248124_1001161", "Ga0307179_10293", "Ga0248079_1001133", "Ga0126444_12068", "Ga0125296_121120", "Ga0136476_102814", 
                     "Ga0136663_11268", "Ga0136672_104612", "Ga0136680_10641", "Ga0136681_105931", "Ga0307205_103728", "Ga0307209_11514", "Ga0307219_10731", 
                     "Ga0248088_10026", "Ga0339969_101193", "Ga0307151_10672", "Ga0307239_106122", "Ga0307166_10772", "Ga0307187_101380", "Ga0256598_104160", 
                     "Ga0307235_10446", "Ga0307229_100424", "Ga0307227_101088", 
                     "Ga0248036_11384", "Ga0136674_102419", "Ga0307235_10559", "Ga0307223_10247", "Ga0307180_101190", "Ga0307235_100815", 
                     "Ga0248123_11815", "Ga0248124_101350", "Ga0307183_12115", "Ga0307185_12236", "Ga0307235_101069", 
                     "Ga0373882_22_63269_63499", "Ga0248036_12357", "Ga0248124_100623", "Ga0248172_1439", "Ga0307235_101054", "Ga0248091_10156", 
                     "Ga0373887_46_4741_7161", "Ga0307196_12752", "Ga0307235_102819", "Ga0307223_101618", "Ga0248094_10291", "Ga0307224_100491", 
                     "Ga0248088_101086", "Ga0248091_101317", "Ga0307196_1501", "Ga0307201_11592", "Ga0307235_103029", "Ga0256589_10680", 
                     "Ga0307235_10624", "Ga0248036_10718", "Ga0307183_11563", "Ga0307202_102287", "Ga0307235_101134", "Ga0248122_10752", 
                     "Ga0248088_10419", "Ga0307235_103027", "Ga0248124_10211", "Ga0248173_100972", "Ga0307206_10146", "Ga0248088_10176", 
                     "Ga0373886_06_188210_189769", "Ga0307193_1166", "Ga0307197_10836", "Ga0256589_10682", "Ga0248094_11631", "Ga0307220_101449", 
                     "Ga0248036_11731", "Ga0248073_108415", "Ga0256598_11031", "Ga0248088_100767", "Ga0248094_114113", "Ga0307226_10191", 
                     "Ga0248073_10291", "Ga0136465_112051", "Ga0136682_12831", "Ga0307210_105423", "Ga0307219_107117", "Ga0248085_1221", 
                     "Ga0307235_11823", "Ga0307229_11072", "Ga0248049_15416", "Ga0307235_101365", "Ga0248072_102253", "Ga0248088_101239", 
                     "Ga0307196_11149", "Ga0373881_025_69180_70064", "Ga0248036_12362", "Ga0307235_100770", "Ga0248086_12715", "Ga0307241_10622", 
                     "Ga0307220_10553", "Ga0307236_10693", "Ga0307228_10783", "Ga0307230_102923", "Ga0248036_10970", "Ga0248072_1374", 
                     "Ga0248074_14217", "Ga0248122_100472", "Ga0248173_10663", "Ga0307180_1343", "Ga0248079_10613", "Ga0136665_11241", 
                     "Ga0136684_11674", "Ga0248088_105714", "Ga0248086_1476", "Ga0307154_11680", "Ga0373887_14_12022_20022", "Ga0307242_1363", 
                     "Ga0307196_14914", "Ga0307201_11670", "Ga0307202_1403", "Ga0307235_101053", "Ga0248122_101243", "Ga0307196_12751", 
                     "Ga0307235_10195", "Ga0307229_12625", "Ga0077198_1121", "Ga0077199_118172", "Ga0248088_103715", "Ga0307235_101248", "Ga0307188_12529", 
                     "Ga0307235_101056", "Ga0248036_12356", "Ga0248124_100624", "Ga0373887_46_8486_10504", "Ga0307239_14723", "Ga0307196_12754", "Ga0307191_1338"))

hap_sequences[which(hap_sequences=="Ga0307256_104726")]

# Compare both objects
which(hap_reference %in% hap_sequences) #145 sequences, All hap_reference are included in hap_sequences
length(which(hap_reference %in% hap_sequences)) #145 sequences
which(!(hap_sequences %in% hap_reference)) # 113 sequences are not included in hap_references, or in the haplotypes_parsed file
length(which(!(hap_sequences %in% hap_reference))) #113 sequences
length(grep("GEV", hap_sequences))
length(grep("Xp", hap_sequences))
length(grep("DRAFT", hap_sequences))


# If no need to eliminate sequences (DRAFT, GEV or Xp)
# Create object with locus tags I need to add to the haplotypes parsed file
length(hap_sequences) #258 sequences
length(hap_reference) #145 sequences

length(which(!(hap_sequences %in% hap_reference)))
toADD <- hap_sequences[which(!(hap_sequences %in% hap_reference))]

#loc_tag[which(loc_tag[,]==toADD)]
a=1
loc_tag[which(loc_tag==toADD[a],arr.ind = T)]
rownames(loc_tag)
colnames(loc_tag)

strADD <- c()
effADD <- c()
for(a in 1:length(toADD)){
  coords_add <- which(loc_tag==toADD[a],arr.ind = T)[1,]
  effect_add <- rownames(loc_tag[coords_add[1],])
  effADD[a] <- effect_add
  strain_add <- colnames(loc_tag[coords_add[2]])
  strADD[a] <- strain_add
}

strADD
effADD

toADD_str_eff <- cbind(toADD, strADD, effADD)
toADD_eff <- matrix(toADD)
rownames(toADD_eff) <- effADD
toADD_eff

write.table(toADD_str_eff, file = "toAdd_loctags.txt", col.names = TRUE, row.names = FALSE, sep= "\t", quote = FALSE)
#ADD manually this locus tags into the log..._transpose file

# - - - - - - - - - - - - - - - 
# OR, if
# Eliminating DRAFT sequences to see how many I need to include in the file "manually"
length(hap_sequences) #258 sequences
hap_draft <- hap_sequences[grep("*DRAFT", hap_sequences)]
length(hap_draft) #164 DRAFT sequences
hap_seq_minDRAFT <- hap_sequences[which(!(hap_sequences %in% hap_draft))]
length(hap_seq_minDRAFT) #262 sequences that are not DRAFT

length(hap_reference) #145 sequences
hapRef_draft <- hap_reference[grep("*DRAFT", hap_reference)]
length(hapRef_draft) #37 sequences
hapRef_minDRAFT <- hap_reference[which(!(hap_reference %in% hapRef_draft))]
length(hapRef_minDRAFT) #147 sequences from the reference alleles that are not "DRAFT"

length(which((hap_seq_minDRAFT %in% hapRef_minDRAFT))) # 147 sequences from the hap_files are included in the references (all the sequences from refs without DRAFT)
length(which(!(hap_seq_minDRAFT %in% hapRef_minDRAFT))) #115 are not included in the references <-- to add!!

# Objet with locus_tags I need to add to the haplotypes parsed file
toADD <- hap_seq_minDRAFT[which(!(hap_seq_minDRAFT %in% hapRef_minDRAFT))]

loc_tag[which(loc_tag[,]==toADD)]
a=1
loc_tag[which(loc_tag==toADD[a],arr.ind = T)]

strADD <- c()
effADD <- c()
for(a in 1:length(toADD)){
  coords_add <- which(loc_tag==toADD[a],arr.ind = T)[1,]
  strain_add <- rownames(loc_tag[coords_add[1],])
  strADD[a] <- strain_add
  effect_add <- colnames(loc_tag[coords_add[2]])
  effADD[a] <- effect_add
}

strADD
effADD

toADD_str_eff <- cbind(toADD, strADD, effADD)
toADD_str_eff

write.table(toADD_str_eff, file = "toAdd_loctags.txt", col.names = TRUE, row.names = FALSE, sep= "\t")



###   TESTING

##Did not worked!!!:
#Add the singletons to the loc_tags matrix
#toADD_str_eff
#toADDmx <- as.matrix(t(toADD))
#colnames(toADDmx) <- effADD
#ncol(toADDmx)
#class(toADDmx)
#newcols <- c(colnames(loc_tags), effADD)
#colnames(loc_tags)
#test <- cbind(loc_tags,toADDmx)
#test <- merge(loc_tags,toADDmx, all=T)
#test <- merge(data.frame(loc_tags, row.names = NULL), data.frame(toADDmx, row.names = NULL), by = 0, all = TRUE)[-1]
#test <- merge(loc_tags, toADDmx, all=T)

### 

# This works but I don't need it right now

# Replace "NA" with NA
for(i in 1:length(allele)){
  if(allele[i] == "NA")
    allele[i] <- NA
}



