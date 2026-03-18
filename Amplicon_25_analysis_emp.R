#install packages ----
install.packages("writexl")
library(writexl)

install.packages("readxl")
library(readxl)

install.packages ("dplyr")
library("dplyr")

install.packages("reshape2")
library(reshape2)

install.packages("tidyr")
library("tidyr")

install.packages ("tidyverse")
library(tidyverse)

install.packages("ggplot2")
library(ggplot2)

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)

library(ShortRead)


BiocManager::install("DECIPHER")
library(DECIPHER)
packageVersion("DECIPHER")

library(Biostrings)

install.packages("seqinr")
library(seqinr)

install.packages("ape")
library(ape)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("msa")
library(msa)

BiocManager::install("remotes")
BiocManager::install("C://Users//erinv//Desktop//cutadapt (2).exe")

install.packages("phyloseq")
"C:\Users\erinv\Desktop\dada2-1.16"
"C:\Users\erinv\Desktop\dada2-master"

library(phyloseq)

BiocManager::install("phyloseq")

install.packages("phangorn")
library("phangorn")

install.packages("C:\\Users\\erinv\\Desktop\\dada2-master",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))
install.packages("C:\\Users\\erinv\\Desktop\\dada2-1.16",
                 repos = NULL,
                 type = "source",
                 dependencies = c("Depends", "Suggests","Imports"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead", version = "3.18") 

library("dada2")
packageVersion("dada2")

#check what version of R we're running
R.version

#save the whole environment
save.image("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//Amplicon_25_analysis_emp.RData")
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//Amplicon_25_analysis_emp.RData")

load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//emp_blast.RData")

#clear the environment
rm(list=ls())

#save individual objects
saveRDS(errF_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errF_para.rds")
saveRDS(errR_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errR_para.rds")


rm(qualityPlots, qualityPlots_paraF_filtered, qualityPlots_paraF_trimmed, qualityPlots_paraR_filtered, qualityPlots_paraR_trimmed)

#ctd----

#tracking reads through workflow ----
#that made it through the pipeline (EDIT this before running to make sure object names are accruate)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#I think this is the same code but just pasting it to be sure
getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(fnFs, getN), 
               sapply(fnRs, getN), 
               sapply(fnFs.noprim, getN), 
               sapply(fnRs.noprim, getN), 
               sapply(filtFs, getN), 
               sapply(filtFs, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("out1", "out2", "inputF", "inputR", "noprimF", "nprimR", "filteredF", "filteredR")
rownames(track) <- sample.names
head(track)
#ctd----

#filepath to raw unzipped fwd and reverse fastas----

#import files
#"C:\Users\erinv\Box\Cryptocercus_research\Amplicon_Data_25\amplicon_data\fastq"
sample_path <- "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//amplicon_data//fastq"
list.files(sample_path)

#when downloading files from file zilla be careful - they do not necessarily download in order
#so if I on;y want to use some files in a folder or sequencing run
#better to download them all then filter them in R rather than manually copy and pasting the chunk of files I need
#because some will get lost
#ctd----

#sort emp sequences into fwd and reverse reads----
#(assuming R1 is fwd and R2 is reverse)
fnFs_EMP <- sort(list.files(sample_path, pattern="-EMP_.*_R1_001.fastq", full.names = TRUE))
fnRs_EMP <- sort(list.files(sample_path, pattern="-EMP_.*_R2_001.fastq", full.names = TRUE))
print(fnFs_EMP)
print(fnRs_EMP)

#create short sample names for each fasta
EMP.namesF <- substr(basename(fnFs_EMP),1,nchar(basename(fnFs_EMP))-16)
print(EMP.namesF)

EMP.namesR <- substr(basename(fnRs_EMP),1,nchar(basename(fnRs_EMP))-16)
print(EMP.namesR)

#ctd----

#remove primers from EMP samples----
#assign fwd and reverse primers
#longer than traditional primers because they have the nex adapter on them for library prep
#nex-EMP-515F
FWD_EMP <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGYCAGCMGCCGCGGTAA" 
#nex-EMP-806R
REV_EMP <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACNVGGGTWTCTAAT" 

#verify presence and orientation of primers in data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD_EMP.orients <- allOrients(FWD_EMP)
REV_EMP.orients <- allOrients(REV_EMP)
FWD_EMP.orients #all orientations of the forward primer (fwd, comp, reverse, reverse comp)
REV_EMP.orients #all orientations of the reverse primer (fwd, comp, reverse, reverse comp)

#count number of primers in all orientations
#doing this for just one sample is fine because they were prepped using the same library
#files are zipped after the last step, do not unzip them or this chunk of code won't work

#function to count primers
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#if sread doesn't work reload the ShortRead package
library(ShortRead)

#count primers
rbind(FWD.ForwardReads_EMP= sapply(FWD_EMP.orients, 
                                   primerHits, 
                                   fn = fnFs_EMP[[6]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_EMP = sapply(FWD_EMP.orients,
                                    primerHits, 
                                    fn = fnRs_EMP[[6]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_EMP = sapply(REV_EMP.orients, 
                                    primerHits,
                                    fn = fnFs_EMP[[6]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_EMP = sapply(REV_EMP.orients, 
                                    primerHits, 
                                    fn = fnRs_EMP[[6]])) #count all orientations of reverse primer in reverse reads
#output is table of counts
#rows are primers (FWD, REV) and read direction (forward, or reverse)
#read direction as in there is a forward and reverse fasta file
#columns are orientations of primers (forward, complement, reverse, reverse comp)
#why are there y's in some sequences what does that mean? I mean seriously wtf

#there are reverse compliments of some forward primers in the reverse reads
#there are reverse compliments of some reverse primers in the forward reads
#consistent across multiple samples 
#(make sure to test this! change the number in the double brackets in the count primers code)

#                        Forward Complement Reverse RevComp
#FWD.ForwardReads_EMP       0          0       0       0
#FWD.ReverseReads_EMP       0          0       0     182
#REV.ForwardReads_EMP       0          0       0     159
#REV.ReverseReads_EMP       0          0       0       0

#remove primers using cutadapt
#https://benjjneb.github.io/dada2/ITS_workflow.html
#downloaded cutadapt single file executable for windows
#https://github.com/marcelm/cutadapt/releases
#stuck it in Box

#checking that cutadapt is installed and accessible in R
system('"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe" --version')
cutadapt <- "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe"

#create output directory for cutadapted files
path.cut_EMPF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnFs_noprim_EMP")
path.cut_EMPR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnRs_noprim_EMP")
if(!dir.exists(path.cut_EMPF)) dir.create(path.cut_EMPF)
if(!dir.exists(path.cut_EMPR)) dir.create(path.cut_EMPR)
fnFs.cut_EMP <- file.path(path.cut_EMPF, EMP.namesF)
fnRs.cut_EMP <- file.path(path.cut_EMPR, EMP.namesR)

FWD_EMP.RC <- dada2:::rc(FWD_EMP)
REV_EMP.RC <- dada2:::rc(REV_EMP)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags_EMP <- paste("-g", FWD_EMP, "-a", REV_EMP.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags_EMP <- paste("-G", REV_EMP, "-A", FWD_EMP.RC) 
# Run Cutadapt
for(i in seq_along(fnFs_EMP)) {
  system2(cutadapt, args = c(R1.flags_EMP, R2.flags_EMP, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_EMP[i], "-p", fnRs.cut_EMP[i], # output files
                             fnFs_EMP[i], fnRs_EMP[i])) # input files
}

#count primers after cutadapt to confirm they have been removed
#repeat for several samples just in case
rbind(FWD.ForwardReads_EMP = sapply(FWD_EMP.orients, 
                                    primerHits, 
                                    fn = fnFs.cut_EMP[[5]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_EMP = sapply(FWD_EMP.orients,
                                    primerHits, 
                                    fn = fnRs.cut_EMP[[5]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_EMP = sapply(REV_EMP.orients, 
                                    primerHits,
                                    fn = fnFs.cut_EMP[[5]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_EMP = sapply(REV_EMP.orients, 
                                    primerHits, 
                                    fn = fnRs.cut_EMP[[5]])) #count all orientations of reverse primer in reverse reads
#It worked. YAY!
#ctd----

#inspect quality plots for EMP samples ----

#store quality plots for pre-trimmed emp F samples
qualityPlots_empF <- plotQualityProfile(fnFs_EMP[1:145])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//qualityPlots_empF.png", plot = qualityPlots_empF, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#store quality plots for pre-trimmed oxy R samples
qualityPlots_empR <- plotQualityProfile(fnRs_EMP[1:145])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//qualityPlots_empR.png", plot = qualityPlots_empR, width = 25, height = 15, dpi = 300, limitsize = FALSE)
#these look good, quality scores are consistently around 35-40 
#and there is some quality decline in the middle or end of reads for a few samples but not most
#there is no need to trim or truncate reads for these samples


#visualize quality scores for cutadapted samples
#zero length reads were introduced by cutadapt (trouble shooting that I used to figure this out is in the para R file)
#these will be removed in the filtering step

#remove zero length reads for visualization purposes
#I did this for para and oxy but I'm going to skip it for emp to save time and go straight to the filter and trim step
#create output directory for zero filtered files
path.nozero_empF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp", "fnFs_nozero_emp")
path.nozero_empR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp", "fnRs_nozero_emp")
if(!dir.exists(path.nozero_empF)) dir.create(path.nozero_empF)
if(!dir.exists(path.nozero_empR)) dir.create(path.nozero_empR)
fnFs.nozero_emp <- file.path(path.nozero_empF, EMP.namesF)
fnRs.nozero_emp <- file.path(path.nozero_empR, EMP.namesR)

#filter out zero length reads for fwd reads
input_filesF <- fnFs.cut_EMP
output_filesF <- fnFs.nozero_emp

for (i in seq_along(input_filesF)) {
  # Read the FASTQ file
  fq <- readFastq(input_filesF[i])
  # Filter out zero-length reads
  fq_filtered <- fq[width(sread(fq)) > 0]
  # Write the filtered FASTQ file
  writeFastq(fq_filtered, output_filesF[i], compress = FALSE)
}

#filter out zero length reads for reverse reads
input_filesR <- fnRs.cut_oxy
output_filesR <- fnRs.nozero_oxy

for (i in seq_along(input_filesR)) {
  # Read the FASTQ file
  fq <- readFastq(input_filesR[i])
  # Filter out zero-length reads
  fq_filtered <- fq[width(sread(fq)) > 0]
  # Write the filtered FASTQ file
  writeFastq(fq_filtered, output_filesR[i], compress = FALSE)
}


qualityPlots <- plotQualityProfile(fnRs.nozero_para[1])
qualityPlots

#store quality plots for pre-trimmed para F samples
qualityPlots_oxyF_trimmed <- plotQualityProfile(fnFs.nozero_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyF_trimmed.png", plot = qualityPlots_oxyF_trimmed, width = 25, height = 15, dpi = 300, limitsize = FALSE)


#store quality plots for pre-trimmed para R samples
qualityPlots_oxyR_trimmed <- plotQualityProfile(fnRs.nozero_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyR_trimmed.png", plot = qualityPlots_oxyR_trimmed, width = 25, height = 15, dpi = 300, limitsize = FALSE)
#these plots don't look bad
#most samples have quality decline at around ~200-225 bp in the reverse reads
#but the decline is not severe enough  to warrant specific trimming
#same as the cutadapted para samples, now the bottom red line (indicating the lower 25th percentile) is much lower than in the raw reads
#this could mean that some reads were short and dropped off in quality when primers were trimmed
#the next trimming step should fix this
#ctd----

#filter emp samples----

# create directory for filtered files
path.filt_empF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp", "fnFs_filt_emp")
path.filt_empR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp", "fnRs_filt_emp")
if(!dir.exists(path.filt_empF)) dir.create(path.filt_empF)
if(!dir.exists(path.filt_empR)) dir.create(path.filt_empR)
fnFs.filt_emp <- file.path(path.filt_empF, EMP.namesF)
fnRs.filt_emp <- file.path(path.filt_empR, EMP.namesR)

#filter reads
#ok. so. removing zero length reads in the previous step was actually NOT the move
#it allowed me to plot quality scores and visualize teh post cutadapt data
#but it results in a different number of reads between the fwd and reverse files
#which messes up anything downstream including filter and trim
#so I'm going to do filter and trim on the noprim files, not the nozero files
#these are files that have had primers removed but still contain zero length reads
#running this step in open on demand r studio version 4.4.1 because it takes forever on my laptop

#create compressed versions of the input and output directories in Box
#upload those directories to workpace in OOD
#the files unzip automatically, but they lose their. fastq extension
#put the .fastq extension back on each file in the noprim directories
#(don't need to do anything for the output directories because tehre are no files in them yet - coudl aalaso just create output directories in OOD)


#rename fwd noprim files (OOD)
files <- list.files("/sfs/gpfs/tardis/home/avg6kb/fnFs_noprim_emp", full.names = TRUE)
# Create new names with the .fastq extension
new_names <- paste0(files, ".fastq")
# Rename all files
file.rename(files, new_names)

#rename rev noprim files (OOD)
files <- list.files("/sfs/gpfs/tardis/home/avg6kb/fnRs_noprim_emp", full.names = TRUE)
# Create new names with the .fastq extension
new_names <- paste0(files, ".fastq")
# Rename all files
file.rename(files, new_names)
list.files("/sfs/gpfs/tardis/home/avg6kb/fnRs_noprim_para", pattern = "\\.fastq$")

#filter and trim sequences (OOD: takes 15-20 mins)
filterAndTrim("fnFs_noprim_para", "fnFs_filt_para", "fnRs_noprim_para", "fnRs_filt_para", #specify input and output directories
              truncLen = c(0,0),  # don't truncate reads because qulaity plots showed no decline in quality at beginning or end and primers ahve already been removed
              maxN = 0,  # Discard any reads with ambiguous bases (Ns)
              maxEE = c(2, 2),  # Maximum expected errors allowed (standard parameter)
              truncQ = 2,  # Trim reads at the first instance of quality score ≤ 2 (fairly relaxed parameter, only kicks in if quality is very poor)
              rm.phix = TRUE,  # Remove reads matching the phiX genome, control genome commonly used in illumina sequencing (default setting)
              compress = TRUE,  # Compress output files
              multithread = FALSE)  # on windows set multithread = false

#the first time I ran this I got an error saying I had mismatched read lengths
#so I had to check to see which files had uneven numbers of fwd adn reverse reads

library(ShortRead)

# Count the number of reads in forward files
fwd_read_counts <- sapply(fnFs_noprim_emp, function(file) {
  fq <- readFastq(file)
  length(fq)  # Count the number of reads
})

# Count the number of reads in reverse files
rev_read_counts <- sapply(fnRs_noprim_emp, function(file) {
  fq <- readFastq(file)
  length(fq)  # Count the number of reads
})

# Print the read counts
print(fwd_read_counts)
print(rev_read_counts)

# Check if the forward and reverse files have matching read counts
matching_reads <- fwd_read_counts == rev_read_counts
print(matching_reads)

#check which read counts don't match between fwd adn reverse reads
mismatched_index <- which(fwd_read_counts != rev_read_counts)
print(mismatched_index)

#this is probably a result of something that happened during sequencing
#remove file 71 and try again
#sample 57H1
#sad
fnFs_noprim_emp <- fnFs_noprim_emp[-71]
fnRs_noprim_emp <- fnRs_noprim_emp[-71]
print(fnFs_noprim_emp)

fnFs.filt_emp <- fnFs.filt_emp[-71]
fnRs.filt_emp <- fnRs.filt_emp[-71]
#now I ran filter and trim again and it worked fine!

#export output files back to Box directory to continue working

#list files within directories
fnFs.filt_emp <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//fnFs.filt_emp", full.names = TRUE)
fnRs.filt_emp <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//fnRsfilt_emp", full.names = TRUE)
print(length(fnFs.filt_emp))

library(dada2)

#store quality plots for filtered para F samples
#print progress reports each time a plot is generated
qualityPlots_empF_filtered <- vector("list", length(fnFs.filt_emp))

for (i in seq_along(fnFs.filt_emp)) {
  qualityPlots_empF_filtered[[i]] <- plotQualityProfile(fnFs.filt_emp[i])
  message("Finished plot ", i, " of ", length(fnFs.filt_emp), ": ", basename(fnFs.filt_emp[i]))
}
library(ggplot2)
install.packages("patchwork")
library(patchwork)

combined_plot_filtF <- wrap_plots(qualityPlots_empF_filtered, ncol = 12)

ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//qualityPlots_empF_filtered.png",
       plot = combined_plot_filtF,
       width = 48, height = 48, dpi = 300, limitsize = FALSE)



#store quality plots for filtered para R samples
qualityPlots_empR_filtered <- vector("list", length(fnRs.filt_emp))

for (i in seq_along(fnRs.filt_emp)) {
  qualityPlots_empR_filtered[[i]] <- plotQualityProfile(fnRs.filt_emp[i])
  message("Finished plot ", i, " of ", length(fnRs.filt_emp), ": ", basename(fnRs.filt_emp[i]))
}

 combined_plot_filtR <- wrap_plots(qualityPlots_empR_filtered, ncol = 12)

ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//qualityPlots_empR_filtered.png",
       plot = combined_plot_filtR,
       width = 48, height = 48, dpi = 300, limitsize = FALSE)

#this is a different method for saving plots than I used for oxy and para,
#but I think it's better because I can zoom in on individual plots
#whereas with the other method I can't. If I remember correctly


#although the mean quality score across all reads is still high
#I think we're good to move on down the pipeline. YAY!
#ctd----

#merge sequences----

#before merging sequences and removing chimeras,
#track reads through workflow so far to figure out how many reads have been removed/ retained at each step
#takes about 20 minutes on my laptop for 140 samples
library(dada2)

#make sure these file paths are correctly listed
fnFs.cut_emp <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//fnFs_noprim_EMP", full.names = TRUE)
fnRs.cut_emp <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//fnRs_noprim_EMP", full.names = TRUE)

getN <- function(x) sum(getUniques(x))   #create function to count reads in a file
track_emp <- cbind(sapply(fnFs_EMP, getN),   #count reads in raw fwd files
                   sapply(fnRs_EMP, getN),   #count reads in raw reverse files
                   sapply(fnFs.cut_emp, getN),   #count reads in fwd cutadapted files
                   sapply(fnRs.cut_emp, getN),   #count reads in reverse cutadapted files
                   sapply(fnFs.filt_emp, getN),   #count reads in fwd filtered files
                   sapply(fnRs.filt_emp, getN))   #count reads in reverse filtered files

print(length(fnFs_EMP))
print(length(fnRs_EMP))
print(length(fnFs.cut_EMP))
print(length(fnRs.cut_EMP))
#this last two are shorter because 
#for sample 57H1 reads didn't match between fwd and reverse
#after primer removal so I removed the whole sample from the data set
print(length(fnFs.filt_emp))
print(length(fnRs.filt_emp))

#so instead we're going to use this code to account for differences in numbers of samples between the files
# Helper function to count reads and return a named vector
getN_named <- function(file_list) {
  counts <- sapply(file_list, getN)
  sample_names <- sub("_.*", "", basename(names(counts)))  # adjust to extract just the sample name
  names(counts) <- sample_names
  return(counts)
}

# Get named vectors (each takes 5 mins to run)
rawF  <- getN_named(fnFs_EMP)
rawR  <- getN_named(fnRs_EMP)
cutF  <- getN_named(fnFs.cut_emp) #warning: zero length sequences detected during dereplication
cutR  <- getN_named(fnRs.cut_emp) #warning: zero length sequences detected during dereplication
filtF <- getN_named(fnFs.filt_emp)
filtR <- getN_named(fnRs.filt_emp)

# Convert each named vector to a data frame with sample names as a column
df_rawF  <- data.frame(Sample = names(rawF), rawF = rawF)
df_rawR  <- data.frame(Sample = names(rawR), rawR = rawR)
df_cutF  <- data.frame(Sample = names(cutF), cutF = cutF)
df_cutR  <- data.frame(Sample = names(cutR), cutR = cutR)
df_filtF <- data.frame(Sample = names(filtF), filtF = filtF)
df_filtR <- data.frame(Sample = names(filtR), filtR = filtR)

# Use dplyr::full_join to merge them all by Sample name
library(dplyr)
library(purrr)

track_emp <- list(df_rawF, df_rawR, df_cutF, df_cutR, df_filtF, df_filtR) %>%
  reduce(full_join, by = "Sample")

View(track_emp)
#57H1 got duplicated because the numbers don't match for fwd and reverese
# so I'm just going to removed 57H1 entirely from the tracking frame
#ideally I wanted to track reads for that sample but it's being a pain
track_emp <- track_emp[-c(71, 72), ]

# Optional: Set Sample as rownames and drop the column
rownames(track_emp) <- track_emp$Sample
track_emp$Sample <- NULL
#renames rows
rownames(track_emp) <- gsub("-EMP$", "", rownames(track_emp))

print(nrow(track_emp))
#yes! finally the right number of rows

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_emp) <- c("rawF_emp", "rawR_emp", "noprimF_emp", "noprimR_emp", "filteredF_emp", "filteredR_emp")

#ok, this looks good! most samples decrease in read as we move through the pipeline,
#a few have a very small number of reads <1000 but most look good
#fwd and reverse reads  have the same numbers of reads for each sample after the filtering step

#estimate error rates
#not a problem that it's only using a subset of samples it just needs x number of bases
#because error rate is a feature of the sequencing run, not of individual samples or reads. I think?
#each error estimation take about 20 minutes on my laptop (40 minutes total)

errF_emp <- learnErrors(fnFs.filt_emp, multithread=FALSE)
#120890799 total bases in 527005 reads from 8 samples will be used for learning the error rates..

errR_emp <- learnErrors(fnRs.filt_emp, multithread=FALSE)
#109815834 total bases in 527005 reads from 8 samples will be used for learning the error rates.

#visualize error rates
plotErrors(errF_emp, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok

plotErrors(errR_emp, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok
#the error plots aren't identical for fwd and reverse but similar 


#takes 2 seconds on laptop
derepFs_emp <- derepFastq(fnFs.filt_emp, verbose=TRUE)
print(names(derepFs_emp))
names(derepFs_emp) <- sub("-EMP.*", "", rownames(derepFs_emp))
names(derepFs_emp) <- emp.namesF
print(EMP.namesF)
print(names(derepFs_emp))

#takes 2 seconds on laptop
derepRs_emp <- derepFastq(fnRs.filt_emp, verbose=TRUE)
names(derepRs_oxy) <- oxy.namesR
#Encountered 75454 unique sequences from 1019071 total sequences read

#estimate variants - takes 25 minutes on my laptop for each, 50 minutes total
dadaFs_emp<- dada(fnFs.filt_emp, err=errF_emp, multithread=FALSE)
dadaRs_emp<- dada(fnRs.filt_emp, err=errR_emp, multithread=FALSE)

#merge fwd and reverse sequences - takes 15 minutes on my laptop
mergers_emp <- mergePairs(dadaFs_emp, derepFs_emp, dadaRs_emp, derepRs_emp)
#got 4 of these printed when I ran that code
#Duplicate sequences in merged output.
#it wasn't labeled as a warning or an error though
#this means that sometimes different combinations of fwd adn reverse seqs gave the same seqeunces after merging
#it's not a problem
head(mergers_emp[[1]])

#construct an asv table (chimeras not yet removed)
seqtabAll_emp <- makeSequenceTable(mergers_emp[!grepl("Mock", names(mergers_emp))])

#shows how many instances of each sequence there are I think
table(nchar(getSequences(seqtabAll_emp)))

#remove chimeras - takes 5 minutes on my laptop
seqtabNoC_emp <- removeBimeraDenovo(seqtabAll_emp)

#calculate the percentage of merged sequences that were chimeric
sum(seqtabNoC_emp)/sum(seqtabAll_emp)
#0.8781728
#chimeras accounted for 13% of merged sequence reads

#get read counts for merged data
merged_emp <- as.data.frame(rowSums(seqtabAll_emp))
print(merged_emp[1])
nrow(merged_emp)
track_empdf <- track_emp
nrow(track_empdf)
track_empdf[7] <- merged_emp[1]

View(track_emp)

#get read counts for after chimera removal and add to log of read counts
nochim_emp <- as.data.frame(rowSums(seqtabNoC_emp))
print(nochim_emp[1])
track_empdf[8] <- nochim_emp[1]

#ctd----


#assign taxonomy silva ----
#assign taxonomy based on reference sequences 
#do taxonomy assignments using open on demand because they crashed my laptop

#export an object with the sequences so I can upload it to open on demand
saveRDS(seqtabNoC_emp, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//seqtabNoC_emp.rds")

seqtabNoC_emp <- readRDS("seqtabNoC_emp.rds")

library(dada2)

# Define your sequence variable (this is usually from your DADA2 output)
seqs <- getSequences(seqtabNoC_emp)  # or whatever your sequence table is called

# Assign taxonomy (up to genus + species)
#takes 20-30 mins OOD
taxa <- assignTaxonomy(seqs, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread = TRUE)
View(taxa)
#Warning message:
#In assignTaxonomy(seqs, "silva_nr99_v138.1_wSpecies_train_set.fa.gz",  :
#Some sequences were shorter than 50 nts and will not receive a taxonomic classification.

# Add species-level information
#takes 20-30 mins OOD
taxa2 <- addSpecies(taxa, "silva_species_assignment_v138.1.fa.gz")
View(taxa2)

#save taxtab file
saveRDS(taxa2, file = "taxTab_silva_emp.rds")

#ctd----


#construct phyloseq object----
#have to do phyloseq stuff on desktop because there are package dependency issues with biostrings in Rivanna
#unable to open biostrings package in rivanna

#load tax tab
taxTab_silva_emp <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//taxTab_silva_emp.rds")

#make data frame of samples names
samdf_emp <- data.frame(Sample=rownames(track_empdf))
head(samdf_emp)

#switch rows and columns of no chimera asv table to get it into phyloseq format
nochim_switch_emp <- data.frame(t(seqtabNoC_emp), check.names = FALSE)

#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_switch_emp) == rownames(taxTab_silva_emp))

#convert to matrices
nochim_mat_emp <- as.matrix(nochim_switch_emp) #convert asv table to matrix
#pr2_mat <- as.matrix(taxTab_pr2)
sam_mat_emp <- as.matrix(samdf_emp) #convert sample data to matrix

View(taxTab_silva_emp)
#convert to df
taxTab_test <- as.data.frame(taxTab_silva_emp)


#convert to matrix
silva_mat_emp <- as.matrix(taxTab_test)
View(silva_mat_emp)


#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_mat_emp) == rownames(silva_mat_emp))
#should print TRUE
setdiff(rownames(nochim_mat_emp), rownames(silva_mat_emp))
#should print character (0)

#check that samples names are the same between metadata and asvtable
all(colnames(nochim_mat_emp) == sam_mat_emp[,1])
#should print TRUE
#it print FALSE but I think that's because there are .fastq extensions still in the nochim_mat_para
#remove fastq extensions
colnames(nochim_mat_emp) <- gsub(".fastq", "", colnames(nochim_mat_emp))

#now check if all sample names are the sam between metadata and asv table
all(colnames(nochim_mat_emp) == sam_mat_emp[,1])
print(colnames(nochim_mat_emp))
print(sam_mat_emp[,1])
#the nochim_mat still has extra stuff in the colnames so let's remove it

colnames(nochim_mat_emp) <- sub("-EMP.*", "", colnames(nochim_mat_emp))
#now check if all sample names are the sam between metadata and asv table
all(colnames(nochim_mat_emp) == sam_mat_emp[,1])
#prints TRUE
#yay!

library(readxl)
#read in the metadata from the experiment
habitat_swap_meta_data <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//habitat_swap_meta_data.xlsx")
(View(habitat_swap_meta_data))

#merge experimental metadata with samdf
colnames(habitat_swap_meta_data)[4] <- "Sample"    #rename "roach_ID column to "Sample"
samdf_emp$Sample <- sub("-.*", "", samdf_emp$Sample)   #remove all but the roach id from the sample names in samdf
samdf_merge_emp <- merge(samdf_emp, habitat_swap_meta_data, by = "Sample")    #merge sam df with metadata


#create a rownames column
samdf_merge_emp$Extra_Column <- c(1:nrow(samdf_merge_emp))

#compiling phyloseq object

#remove all but roach ID from columns in nochim matrix
colnames(nochim_mat_emp) <- sub("-.*", "", colnames(nochim_mat_emp))

View(samdf_merge_emp)

#convert the first column of meta data into rownames (essential for phyloseq to run)
rownames(samdf_merge_emp) <- samdf_merge_emp[,1]
samdf_merge_emp[,1] <- NULL

#meta data has to be a data frame and the rownames have to be the samples for this to work. Jeez
#taxa are rows = TRUE if taxa are rownames in the asv matrix
library(phyloseq)
ps_emp <- phyloseq(
  otu_table(nochim_mat_emp, taxa_are_rows = TRUE), 
  sample_data(samdf_merge_emp), 
  tax_table(silva_mat_emp))
#phy_tree(fitJC_para$tree)) no tree for now
#ctd----

#filter phyloseq object and wrangle sample data-----
#filtering and sample wrangling used to be two separate sections 
#with filtering first, then calculate relative abundances, then filter out problem samples
#now I'm going to filter problem samples before calculating relative abundances
#I'm also going to track read counts per sample at every stage of the filtering process
#then after all of that I'm goingto filter by the 0.1 % relative abundance parameter

#ok, it's this point where I'm realizing I need to explore some different visualization and analysis approaches
#this rescource looks helpful for initial visualization of normalized data
#https://www.nicholas-ollberding.com/post/introduction-to-the-statistical-analysis-of-microbiome-data-in-r/

#but I need options for best way to normalize the data
#it seems pretty standard to just use relative abundances
#but can the relative abundances also be used for statistical analysis?
#and is there any better normalization than relative abundance?


#this rescoucre discusses different options for data normalization but doesn't discuss pros adn cons of each
#honestly I think I'm safe to use relative abundance though. idk why I was in my head about it
#https://mcic-osu.github.io/2020-12-microbiomics-workshop/08-postASV-analysis.html

#this rescource compares different normalization methods. 
#it looks like the deseq vst stuff is more commonly used for rna seq
#While the microbial ecologists prefer proportions to get a equivalent library size
#so we're gonna go with relative abundances
#https://scienceparkstudygroup.github.io/microbiome-lesson/05-data-filtering-and-normalisation/index.html


#after the taxonomy assignment step, I got a warning saying some taxonomy couldn't be assigned
#because of short sequences (though the majority of sequences are not short)
#Initially:
#to be conservative, I'm going to filter out any sequences less than 251 bp (the expected length)
#assuming anything shorter is due to sequencing error and not biologically relevant
#Edit 05_07_25: 
#I spoke with Gillian and she said we expect a range or lengths for this region
#and we should filter out anything less than 200 bp but leave the rest

#remove ASV's with sequences shorter than 200 bp
asv_seqs <- taxa_names(ps_emp)
print(range(nchar(asv_seqs)))
#21 461
valid_asvs <- asv_seqs[which(nchar(asv_seqs) >= 200)] # Keep only ASVs that are at least 200 bp long
print(length(valid_asvs))
print(length(asv_seqs))
ps_emp_filt <- prune_taxa(valid_asvs, ps_emp) # Prune phyloseq object to retain only valid ASVs

#inspect
print(nrow(tax_table(ps_emp_filt)))
print(nrow(tax_table(ps_emp)))

#track reads per sample 
read_counts <- sample_sums(ps_emp_filt)
track_empdf$filt_200bp <- read_counts[rownames(track_empdf)]

# create table to track ASV counts through filtering process (track_paradf already tracks reads per sample)
track_emp_asv <- data.frame(
  step = character(),
  n_ASVs = integer(),
  stringsAsFactors = FALSE
)


# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "merged", n_ASVs = ntaxa(ps_emp)))
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_200bp", n_ASVs = ntaxa(ps_emp_filt)))


#wrangle meta data

#inspect
View(sample_data(ps_emp_filt))

#change 45C3 to instar 3 (it was initially called instar 2 but later reevaluated and included in the experiment)
meta_filt <- sample_data(ps_emp_filt)
#meta_filt["45C3", "Instar"] <- 3
#View(meta_filt)
#not doing this becaseu 45C3 isn't in this data
View(meta_filt)

#update get rid of extra column
meta_filt$Extra_Column <- NULL

#add weight at death
meta_filt$death_weight <- NA

#convert to data frame so I can mutate
meta_filtdf <- data.frame(sample_data(meta_filt))
View(meta_filtdf)
print(class(meta_filtdf))

library(dplyr)
meta_filtdf<- meta_filtdf %>%
  mutate(death_weight = case_when(
    Dissected == "07_12_24" ~ weight_07_05_24,
    Dissected == "07_19_24"~ weight_07_19_24,
    Dissected == "07_26_24"~ weight_07_26_24,
    Dissected == "07_27_24"~ weight_07_27_24
  ))

meta_filtdf$site <- NA #make site column
colnames(meta_filtdf)[1] <- "log" #fix the name of the log column, it didn't read in right for some reason

#read in log site key so I can put a site variable in the data set
library(readxl)
sites <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Log_site_key.xlsx")
View(sites)
colnames(sites) <- c("log", "site") #rename sites columns so they start with lowercase  letters and match sample data

meta_filtdf <- meta_filtdf %>%
  left_join(sites %>% select(log, site), by = "log") %>%
  mutate(site = ifelse(is.na(site.x), site.y, site.x)) %>%
  select(-site.x, -site.y)

#put the edited sample data back in the phyloseq object
rownames(meta_filtdf) <- rownames(meta_filt)
sample_data(ps_emp_filt) <- sample_data(meta_filtdf)
#inspect
View(sample_data(ps_emp_filt))


#remove individuals accidentally swapped onto hardwood
#samples_to_remove <- c("41C1", "41F1", "41F2")  
#don't need to do this, these samples also aren't here

# Remove the selected samples from the phyloseq object
#ps_oxy_filt2 <- prune_samples(!(sample_names(ps_oxy_filt) %in% samples_to_remove), ps_oxy_filt)
# Remove ASV's associated with only those samples
#ps_oxy_filt2 <- prune_taxa(taxa_sums(ps_oxy_filt2) > 0, ps_oxy_filt2)


#inspect
#print(nrow(tax_table(ps_oxy_filt)))
#print(nrow(tax_table(ps_oxy_filt2)))

#track reads per sample
#read_counts <- sample_sums(ps_oxy_filt2)
#track_oxydf$filt_badswap <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
#track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_badswap", n_ASVs = ntaxa(ps_oxy_filt2)))



#filter out those families of instar 3 nymphs that violate independence

# Get sample names
all_samples <- sample_names(ps_emp_filt)
print(all_samples)
print(length(all_samples))

# Identify samples starting with "56A" or "57C"
target_samples <- grep("^56A|^57C", all_samples, value = TRUE)
print(target_samples)
print(length(target_samples))
# Subset phyloseq object to include only those samples
ps_violators_emp <- prune_samples(target_samples, ps_emp_filt)
# Optional: prune unused taxa
ps_violators_emp <- prune_taxa(taxa_sums(ps_violators_emp) > 0, ps_violators_emp)

#create column for treatment so that I can evaluate which groups were in what treatment
# Extract the sample data as a data frame
sample_data_df <- as.data.frame(sample_data(ps_violators_emp))
# Create new "treatment" column by pasting two existing columns
sample_data_df$treatment <- paste0(sample_data_df$type_original_log, ":", sample_data_df$type_log_new)
# Assign updated sample data back to the phyloseq object
sample_data(ps_violators_emp) <- sample_data_df
View(sample_data(ps_violators_emp))

#find the samples for each treatment group with the highest read counts
# Extract read counts
sample_data_df$reads <- sample_sums(ps_violators_emp)
View(sample_data_df)
#add back the version with read counts
sample_data(ps_violators_emp) <- sample_data_df

# Convert sample data to regular data frame
sample_df_raw <- as.data.frame(as(sample_data(ps_violators_emp), "data.frame"))
View(sample_df_raw)
sample_df_raw$sample_id <- rownames(sample_df_raw)

# Get top 1 sample (highest reads) per treatment group
top_sample_ids <- sample_df_raw %>%
  group_by(treatment) %>%
  slice_max(reads, n = 1, with_ties = FALSE) %>%
  pull(sample_id)

print(top_sample_ids)

print(sample_df_raw$reads)

#now remove all the samples from the phyloseq object that are from 56A or 57C 
#EXCEPT the ones with the highest read counts in each treatment group

# Get all sample names
all_samples <- sample_names(ps_emp_filt)
# Identify samples that start with 56A or 57C
to_remove <- grep("^56A|^57C", all_samples, value = TRUE)
# Subtract the ones you want to keep (top_sample_ids)
final_remove <- setdiff(to_remove, top_sample_ids)
print(final_remove) #inspect
print(length(final_remove))
# Prune the samples you want to remove
ps_emp_filt2 <- prune_samples(!(sample_names(ps_emp_filt) %in% final_remove), ps_emp_filt)
#prune ASV's that were assocaited with those samples
ps_emp_filt2 <- prune_taxa(taxa_sums(ps_emp_filt2) > 0, ps_emp_filt2)

#inspect
print(nrow(tax_table(ps_emp_filt)))
print(nrow(tax_table(ps_emp_filt2)))

#track reads per sample
read_counts <- sample_sums(ps_emp_filt2)
track_empdf$filt_violators <- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_violators", n_ASVs = ntaxa(ps_emp_filt2)))



#remove samples with low read counts
class(track_empdf)
str(track_empdf)
#make a list of samples that have 1000 or more reads
valid_samples <- rownames(track_empdf[!is.na(track_empdf$filt_violators) & track_empdf$filt_violators >= 1000, ]) 
print(length(valid_samples))
# Subset the phyloseq object to retain only valid samples
ps_emp_filt3 <- prune_samples(sample_names(ps_emp_filt2) %in% valid_samples, ps_emp_filt2) 
#prune ASV's not present in any remaining samples
ps_emp_filt3 <- prune_taxa(taxa_sums(ps_emp_filt3) > 0, ps_emp_filt3)

#inspect
print(nrow(tax_table(ps_emp_filt2)))
print(nrow(tax_table(ps_emp_filt3)))

#track reads per sample
read_counts <- sample_sums(ps_emp_filt3)
track_empdf$filt_lowreads <- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_lowreads", n_ASVs = ntaxa(ps_emp_filt3)))




#this should filter out OTU's that have fewer than 10 reads present across all samples
#not sure this is necessary because the relative abundance filtering step would probably catch all of these
#but hey I mean you never know just in case
#To increase the reliability of microbial composition, we advise removing OTUs with 
#<10 copies in individual samples, 
#particularly in studies where only one subsample per specimen is available for analysis. (like this one)
# doi: 10.3389/fcimb.2023.1165295.
#this paper says filtering based on variance decreases reliability compared to the recommended method
filter <- phyloseq::genefilter_sample(ps_emp_filt3, filterfun_sample(function(x) x >= 10))
ps_emp_filt4 <- prune_taxa(filter, ps_emp_filt3)

#inspect
print(nrow(tax_table(ps_emp_filt3)))
print(nrow(tax_table(ps_emp_filt4)))

#track reads per sample
read_counts <- sample_sums(ps_emp_filt4)
track_empdf$filt_10reads <- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_10reads", n_ASVs = ntaxa(ps_emp_filt4)))



#filtering according to Gillian's advice (remove potential contaminants and errors)
View(tax_table(ps_emp_filt4))

#first remove anything that only appears in one sample
#Looking at the taxonomy assignments, there isn’t anything unique or potentially 
#interesting that you don’t already have represented in your higher frequency ASVs.
#Plus, there’s a high likelihood that these are just inaccurate derivatives of real sequences. 
#And anyway I don’t think they’ll be informative in your analyses if they’re only in 1 sample. 
#This seems to be more useful than a read depth cutoff. That said, I only see one ASV remaining 
#in the 2+ samples set that has <100 reads, so you could raise your threshold to that if you wanted to,
#maybe it’s an easier filter. You won’t be losing anything significant. I don’t know if 
#this would apply to the oxy yet because the oxy counts sheet is still in proportional abundances.

# Compute how many samples each ASV appears in
asv_sample_counts <- apply(otu_table(ps_emp_filt4), 1, function(x) sum(x > 0))
head(asv_sample_counts)

# Identify ASVs that appear in more than one sample
asvs_to_keep <- names(asv_sample_counts[asv_sample_counts > 1])

# Prune the phyloseq object to keep only those ASVs
ps_emp_filt5 <- prune_taxa(asvs_to_keep, ps_emp_filt4)

#inspect
print(nrow(tax_table(ps_emp_filt4)))
print(nrow(tax_table(ps_emp_filt5)))
#a good number of these asvs only appeared in one sample

#track reads per sample
read_counts <- sample_sums(ps_emp_filt5)
track_empdf$filt_1_sample <- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_1_sample", n_ASVs = ntaxa(ps_emp_filt5)))


View(tax_table(ps_emp_filt5))


#ok now let's work on trying to sort out contamination for this dataset.
# I got chat to generate a list of bacterial families  that are commonly found as contaminants in emp datasets
# i checked sources and everything looks pretty legit



library(phyloseq)
library(dplyr)
library(tibble)

# Define list of contaminant families
contaminant_families <- c(
  "Staphylococcaceae", "Streptococcaceae", "Corynebacteriaceae", "Micrococcaceae",
  "Cutibacteriaceae", "Neisseriaceae", "Enterobacteriaceae", "Comamonadaceae",
  "Sphingomonadaceae", "Burkholderiaceae", "Bradyrhizobiaceae", "Methylobacteriaceae",
  "Delftiaceae", "Pseudomonadaceae", "Bacillaceae", "Paenibacillaceae"
)

# Set thresholds
min_prop_abundance <- 0.01  # e.g., 1% of total reads
min_sample_prevalence <- 0.1  # e.g., present in >10% of samples

# Convert to data frame
otu_table_df <- as.data.frame(otu_table(ps_emp_filt5))
if (taxa_are_rows(ps_emp_filt5)) otu_table_df <- t(otu_table_df)

# Relative abundance per taxon
otu_total <- colSums(otu_table_df)
otu_rel_abundance <- otu_total / sum(otu_total)

# Tax table
tax <- as.data.frame(tax_table(ps_emp_filt5))

# Filter for contaminant families
contam_taxa <- rownames(tax)[tax$Family %in% contaminant_families]

# Initialize list of taxa to remove
taxa_to_remove <- c()

library(tidyr)

# Create an empty data frame to store warning info
warning_table <- tibble(
  Taxon = character(),
  Family = character(),
  Prevalence = numeric(),
  RelativeAbundance = numeric()
)

# Loop over each potential contaminant taxon
for (taxon in contam_taxa) {
  abund <- otu_rel_abundance[taxon]
  prevalence <- sum(otu_table_df[, taxon] > 0) / nsamples(ps_emp_filt5)
  
  if (abund > min_prop_abundance | prevalence > min_sample_prevalence) {
    warning_table <- tibble::add_row(warning_table,
                             Taxon = taxon,
                             Family = tax[taxon, "Family"],
                             Prevalence = prevalence,
                             RelativeAbundance = abund
    )
  } else {
    taxa_to_remove <- c(taxa_to_remove, taxon)
  }
}

library(dplyr)

# Print warning summary
if (nrow(warning_table) > 0) {
  cat("\n⚠️  The following potential contaminants were retained due to high abundance or prevalence:\n\n")
  print(warning_table %>%
          mutate(
            Prevalence = sprintf("%.1f%%", 100 * Prevalence),
            RelativeAbundance = sprintf("%.2f%%", 100 * RelativeAbundance)
          ))
} else {
  cat("\n✅ No high-abundance contaminants were detected among flagged families.\n")
}

contam_warnings_df <- data.frame(warning_table)
View(contam_warnings_df)

# Filter the phyloseq object
ps_emp_filt6 <- prune_taxa(setdiff(taxa_names(ps_emp_filt5), taxa_to_remove), ps_emp_filt5)

#inspect
print(nrow(tax_table(ps_emp_filt5)))
print(nrow(tax_table(ps_emp_filt6)))
#not too many of these asvs were common lab or human-based contaminants

#track reads per sample
read_counts <- sample_sums(ps_emp_filt6)
track_empdf$filt_labcontam <- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_labcomtam", n_ASVs = ntaxa(ps_emp_filt6)))



#oops forgot to filter out eukaryotes and there are 20 left! let's get rid of them now
ps_emp_filt7 <- subset_taxa(ps_emp_filt6, Kingdom != "Eukaryota")

#inspect
print(nrow(tax_table(ps_emp_filt6)))
print(nrow(tax_table(ps_emp_filt7)))

#track reads per sample
read_counts <- sample_sums(ps_emp_filt7)
track_empdf$filt_euks<- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_euks", n_ASVs = ntaxa(ps_emp_filt7)))

library(phangorn)
library(phyloseq)


#from here on, we don't know how to tell what is a termite and what is from crypto
#I suspect that there is termite contamination, given the para and oxy data had contamination
#but the literture doesn't really have info on which taxa belong to termites vs. crypto
#so I have to blast them

library(phyloseq)
library(Biostrings)
library(dplyr)

#first export seqeunces, tax info, and unique ASV IDs to a fasta file

# Extract sequences (rownames of taxa)
seqs <- taxa_names(ps_emp_filt7)

# Create unique ASV IDs
asv_ids <- paste0("ASV_", seq_along(seqs))

# Get tax table as a data.frame
tax_df <- as.data.frame(tax_table(ps_emp_filt7))
tax_df$Sequence <- seqs
tax_df$ASV_ID <- asv_ids

# Collapse taxonomy into one string (optional: change separator if needed)
tax_df$Taxonomy <- apply(tax_df[, 1:(ncol(tax_df) - 2)], 1, function(x) paste(na.omit(x), collapse = ";"))

# Create DNAStringSet with FASTA headers
fasta_seqs <- DNAStringSet(seqs)
names(fasta_seqs) <- paste0(tax_df$ASV_ID, " | ", tax_df$Taxonomy)

# Write to FASTA
writeXStringSet(fasta_seqs, filepath = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//emp_contam_check.fasta")

View(tax_df)



#ok! Vera looked at the data and suspects there is very little contamination
#this matches with Gillian's expectations as well
#she said it is less likely there is contamination with the bacteria than the para or oxy
#because they weren't cloning that and haven't amplified emp sequences as often or recently

#she mentioned that some species level taxa assignments look funky, but I guess
#that's because the assignment database has more termite info so it probably assigned them as termites

#vera flagged one potential contaminant but also said it might not be
#then she said to be conservative I should filter out low abundance taxa
#anything that has less than 50 or 100 reads total
# so I think I'll do that


#so here's the plan
#filter out archea
#filter out low abundance (less than 100 reads)
#filter toss any samples with < 1000 reads


#remove archea
View(tax_table(ps_emp_filt7))
ps_emp_filt8 <- subset_taxa(ps_emp_filt7, Kingdom != "Archaea")

#inspect
print(nrow(tax_table(ps_emp_filt7)))
print(nrow(tax_table(ps_emp_filt8)))
#7 archea removed

#track reads per sample
read_counts <- sample_sums(ps_emp_filt8)
track_empdf$filt_archs<- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_archs", n_ASVs = ntaxa(ps_emp_filt8)))



#cool! now filter by abundance. remove anything with less than 100 total reads
ps_emp_filt9 <- filter_taxa(ps_emp_filt8, function(x) sum(x) >= 50, prune = TRUE)
#inspect
print(nrow(tax_table(ps_emp_filt8)))
print(nrow(tax_table(ps_emp_filt9)))
#332 taxa removed with 100 reads cutoff
#0 taxa removed with 10 reads cutoff (probably I already did this lol)
#163 taxa removed with 50 reads cutoff
#let's do 50
#so we're left with 882 taxa


#track reads per sample
read_counts <- sample_sums(ps_emp_filt9)
track_empdf$filt_50reads<- read_counts[rownames(track_empdf)]

# Add entries to ASV tracking data frame
track_emp_asv <- rbind(track_emp_asv, data.frame(step = "filt_50reads", n_ASVs = ntaxa(ps_emp_filt9)))


#now check if any samples have less than 1000 reads
#and they don't! so we're good



#transform sample counts into relative abundances
ps_relative_emp  = transform_sample_counts(ps_emp_filt9, function(x) x / sum(x) )
View(otu_table(ps_relative_emp))
View(tax_table(ps_relative_emp))

#add asv_ID to ps_relative_emp tax table
# Extract the tax_table and convert it to a data.frame
tax_tab <- as.data.frame(tax_table(ps_relative_emp))

# Make sure rownames of tax_tab match something in tax_df$Sequence_ID
# Merge using Sequence_ID as the key
tax_tab$Sequence <- rownames(tax_tab)
merged_tax <- merge(tax_tab, tax_df[, c("Sequence", "ASV_ID")],
                    by = "Sequence", all.x = TRUE)

View(tax_df)
View(merged_tax)

# Restore the rownames (you can choose Sequence_ID or ASV_ID depending on your use case)
rownames(merged_tax) <- merged_tax$Sequence  # or use Sequence_ID if you prefer

# 4. Remove the extra ID columns if not needed
merged_tax$Sequence <- NULL

# 5. Replace the tax_table in your phyloseq object
tax_table(ps_relative_emp) <- tax_table(as.matrix(merged_tax))

p <- plot_bar(ps_relative_emp, fill = "ASV_ID")
p

ps_small_emp = filter_taxa(ps_relative_emp, function(x) sum(x) > .1, TRUE)
print(nrow(tax_table(ps_small_emp)))
p <- plot_bar(ps_small_emp, fill = "ASV_ID")
p

saveRDS(ps_relative_emp, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//ps_relative_emp.rds")


library(ggplot2)
# Save as PNG
#this code isn't working because I can't get the size right so just export it manually as a png
#ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//asv_barplot_para.png", plot = p, width = 854, height = 419, dpi = 300, limitsize = FALSE)


#save tracking data as excel files 

# Load the library
library(tibble)

# Move row names to a new column called "Sample"
track_emp <- rownames_to_column(track_empdf, var = "Sample")
View(track_emp)

# Load the package
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "track_emp")
writeData(wb, sheet = "track_emp", x = track_emp)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//track_emp.xlsx", overwrite = TRUE)



wb <- createWorkbook()
addWorksheet(wb, "track_emp_asv")
writeData(wb, sheet = "track_emp_asv", x = track_emp_asv)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//track_emp_asv.xlsx", overwrite = TRUE)

#ctd----

#make phylogeny----

ps_relative_emp <- readRDS("ps_relative_emp.rds")
rm(ps_relative_emp)

#only include seqeunces in the phylogeny that are still here after filtering

#get ASV sequences from phyloseq object
seqs_rel <- taxa_names(ps_relative_emp)

#propogates labels to tips of the tree
names(seqs_rel) <- seqs_rel

#make alignment
library(DECIPHER)
aligned_rel <- AlignSeqs(DNAStringSet(seqs_rel), processors = 4, anchor=NA,verbose=TRUE)
names(aligned_rel) #check that ASV seqs are names in the alignment

library(phangorn)
# Read aligned sequences into phyDat format
phangAlign_rel <- phyDat(as(aligned_rel, "matrix"), type="DNA")

# Compute pairwise distances using maximum likelihood
dm_rel <- dist.ml(phangAlign_rel, processors = 1, verbose = TRUE)

# Build a neighbor-joining tree
tree_NJ_rel <- NJ(dm_rel)

class(tree_NJ_rel)
# Convert to pml object
fit_rel <- pml(tree_NJ_rel, data = phangAlign_rel)


# Fit the JC69 model
#fit_JC_rel <- optim.pml(update(fit_rel, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#took about an hour for the ~1000 emp seqs on OOD with 768GB mem
#skipping this because the previous version was significantly worse than GTR and we want to save time

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel <- optim.pml(update(fit_rel, model = "GTR", k = 4), 
                         optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#took 1 ish hours

#previous version before euk removal, for this version just went right away with GTR
AIC(fit_JC_rel, fit_GTR_rel)
#> AIC(fit_JC_rel, fit_GTR_rel)
#df       AIC
#fit_JC_rel  2141 112684.11
#fit_GTR_rel 2150  95760.45
#GTR is best, going with that

#check if the trees I plan on using are rooted
is.rooted(fit_GTR_rel$tree)

#the tree is not rooted, so we need to root it now
#because we do not have an outgroup, we will estimate the root using midpoint rooting
#this method places the root in the middle between the two longest branches

# If your tree is stored in a pml object (e.g., fit_GTR_rel$tree)
tree_unrooted <- fit_GTR_rel$tree

# Root using midpoint method
tree_rooted <- midpoint(tree_unrooted)

# Plot to check
plot(tree_rooted, show.tip.label = FALSE)
#this looks better but I'm still seeing a long branch and a small clade


#Identify the root node and its two immediate child clades
library(ape)
# Assuming tree_rooted is rooted and of class "phylo"
root_node <- length(tree_rooted$tip.label) + 1

# Get the two edges coming from the root
root_edges <- tree_rooted$edge[tree_rooted$edge[,1] == root_node, 2]

#Extract tip labels from each daughter clade
# Extract both clades
clade1_tips <- extract.clade(tree_rooted, root_edges[1])$tip.label
clade2_tips <- extract.clade(tree_rooted, root_edges[2])$tip.label

#Pick the smaller clade
small_branch_tips <- if (length(clade1_tips) < length(clade2_tips)) clade1_tips else clade2_tips
print(small_branch_tips)
print(length(small_branch_tips))

#Extract sequences for those tips (from phyloseq)
ps_branch_emp <- prune_taxa(small_branch_tips, ps_relative_emp)
View(tax_table(ps_branch_emp))

#oh ok, all of these are archea
#cool! they seem to be well-classified so I'm goingto leave them in for now
#they're probably serving as a pretty good outgroup for the tree anyway
#next we can pick back up with analysis from here


library(phangorn)
# Extract tax table
tax <- as.data.frame(tax_table(ps_relative_para))

# Match the order of tax_table to the tree tips
tax <- tax[tree_rooted$tip.label, ]

# Create tip labels: for example, use Genus or Genus + Species
# You can adjust this to another rank like Family or Order
tip_labels <- ifelse(
  !is.na(tax$Genus),
  as.character(tax$Genus),
  "Unclassified"
)

tree_label <- tree_rooted

# Assign new labels to tree
tree_label$tip.label <- tip_labels

# Plot with smaller text
plot(tree_label, show.tip.label = TRUE, cex = 0.5)

#plot a version of this that is more legible
# Get tip labels
labels <- tree_label$tip.label

# Replace labels that are exactly "Saccinobaculus" with empty string
labels[labels == "Saccinobaculus"] <- ""


# Assign back to the tree object
tree_label$tip.label <- labels

# Plot the tree with only non-Saccinobaculus labels showing
plot(tree_label, show.tip.label = TRUE, cex = 1)


library(ape)

# Create readable ASV IDs (ASV_1, ASV_2, ...) for each rowname
asv_ids <- paste0("ASV_", seq_len(nrow(tax)))
names(asv_ids) <- rownames(tax)  # rownames are long ASV sequences

# Copy tree
tree_label <- tree_rooted

# Extract current tip labels (long sequences)
tip_seqs <- tree_label$tip.label

# Get genus for each tip
genus <- tax[tip_seqs, "Genus"]

# Build custom labels: ASV_x | Genus for non-Saccinobaculus, blank otherwise
custom_labels <-  paste0(asv_ids[tip_seqs], " | ", genus)


print(custom_labels)

# Assign new labels
tree_label$tip.label <- custom_labels

# Plot
plot(tree_label, show.tip.label = TRUE, cex = 0.6)


#update phyloseq object with asv ids
label_df <- data.frame(
  ASV_Sequence = tree_rooted$tip.label,
  Custom_Label = tree_label$tip.label,
  stringsAsFactors = FALSE
)
View(label_df)


#Extract and convert taxonomy table
tax <- as.data.frame(tax_table(ps_relative_para))

# Add ASV_Sequence as a column for matching (taxa_names are ASV sequences)
tax$ASV_Sequence <- rownames(tax)

# Merge with custom labels
tax_labeled <- merge(tax, label_df, by = "ASV_Sequence", all.x = TRUE)

# Reorder rows to match original phyloseq object
tax_labeled <- tax_labeled[match(taxa_names(ps_relative_para), tax_labeled$ASV_Sequence), ]

# Set rownames and convert back to matrix
rownames(tax_labeled) <- tax_labeled$ASV_Sequence
tax_labeled_mat <- as.matrix(tax_labeled[, -which(names(tax_labeled) == "ASV_Sequence")])

# Replace tax_table in phyloseq
tax_table(ps_relative_para) <- tax_table(tax_labeled_mat)

View(tax_table(ps_relative_para))


ps_relative_emp_tree <- phyloseq(otu_table(ps_relative_emp), 
                                 tax_table(ps_relative_emp), 
                                 tree = tree_rooted)
saveRDS(ps_relative_emp_tree, file = "ps_relative_emp_tree.rds")





#make new tree in OOD after final tax filtering
ps_relative_emp <- readRDS("ps_relative_emp.rds")

#make a new alignment that reflects all the filtering we did
#going to try this with large and filtered datasets
#only include seqeunces in the phylogeny that are still here after filtering

#get ASV sequences from phyloseq object
seqs_rel <- taxa_names(ps_relative_emp)

#propogates labels to tips of the tree
names(seqs_rel) <- seqs_rel

#make alignment
library(DECIPHER)
aligned_rel <- AlignSeqs(DNAStringSet(seqs_rel), processors = 4, anchor=NA,verbose=TRUE)
names(aligned_rel) #check that ASV seqs are names in the alignment

library(phangorn)
# Read aligned sequences into phyDat format
phangAlign_rel <- phyDat(as(aligned_rel, "matrix"), type="DNA")

# Compute pairwise distances using maximum likelihood
dm_rel <- dist.ml(phangAlign_rel, processors = 1, verbose = TRUE)

# Build a neighbor-joining tree
tree_NJ_rel <- NJ(dm_rel)

class(tree_NJ_rel)
# Convert to pml object
fit_rel <- pml(tree_NJ_rel, data = phangAlign_rel)


# Fit the JC69 model
#fit_JC_rel <- optim.pml(update(fit_rel, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#took about an hour for the ~1000 emp seqs on OOD with 768GB mem
#skipping this because the previous version was significantly worse than GTR and we want to save time

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel <- optim.pml(update(fit_rel, model = "GTR", k = 4), 
                         optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#took 1 ish hours

#previous version before euk removal, for this version just went right away with GTR
AIC(fit_JC_rel, fit_GTR_rel)
#> AIC(fit_JC_rel, fit_GTR_rel)
#df       AIC
#fit_JC_rel  2141 112684.11
#fit_GTR_rel 2150  95760.45
#GTR is best, going with that

#check if the trees I plan on using are rooted
is.rooted(fit_GTR_rel$tree)

#the tree is not rooted, so we need to root it now
#because we do not have an outgroup, we will estimate the root using midpoint rooting
#this method places the root in the middle between the two longest branches

# If your tree is stored in a pml object (e.g., fit_GTR_rel$tree)
tree_unrooted <- fit_GTR_rel$tree

# Root using midpoint method
tree_rooted <- midpoint(tree_unrooted)

# Plot to check
plot(tree_rooted, show.tip.label = FALSE)
#this looks better but I'm still seeing a long branch and a small clade

View(tax_table(ps_relative_emp))


# Extract taxonomy table
tax <- as.data.frame(tax_table(ps_relative_emp))

# Handle missing Genus (optional, for cleaner labels)
tax$Genus[is.na(tax$Genus)] <- "Unclassified"

# Create Custom_Label
tax$Custom_Label <- paste(tax$ASV_ID, "|", tax$Genus)

View(tax)

# Assign back to tax_table
tax_table(ps_relative_emp) <- as.matrix(tax)
View(tax_table(ps_relative_emp))



library(ggtree)
library(phyloseq)
library(ggplot2)
library(dplyr)

# Extract tip labels from tree
tree <- tree_rooted

# Get the Custom_Label mapping from phyloseq object
tip_labels_df <- as.data.frame(tax_table(ps_relative_emp))
tip_labels_df$ASV_Sequence <- rownames(tip_labels_df)

# Prepare label lookup
label_lookup <- tip_labels_df %>%
  select(ASV_Sequence, Custom_Label)

# Plot tree and relabel tips
p <- ggtree(tree_rooted) %<+% label_lookup +
  geom_tiplab(aes(label = Custom_Label), size = 3)

p

#yikes ok this is a freaking mess. the tip labels aren't ledgible with
#so many taxa at all.
#but there aren't any small solo clades! 
#so we're gonna go with this tree

#add tree to phyloseq object
# Rebuild phyloseq object manually
ps_relative_emp_tree <- phyloseq(
  otu_table(ps_relative_emp),
  tax_table(ps_relative_emp),
  sample_data(ps_relative_emp),
  phy_tree(tree_rooted)
)

saveRDS(ps_relative_emp_tree, file = "ps_relative_emp_tree.rds")

save.image(file = "OOD_phylogeny_emp_08_07_25")


#load ps with tree back into environment
ps_relative_emp_tree <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//ps_relative_emp_tree.rds")

#load data objects from most recent OOD phylogeny
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//OOD_phylogeny_emp_08_07_25")

saveRDS(tax, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//tax_df_emp.rds")
#ctd----

#PERMANOVA (bray curtis)----
#extract asv table

bray_emp <- phyloseq::distance(ps_relative_emp_tree, method = "bray")
meta_emp <- data.frame(sample_data(ps_relative_emp_tree))
View(meta_emp)

??phyloseq::distance

(View(meta_emp))
#meta data missing site info for two samples
#add it back manually
#don't do this for para because these samples were already removed
meta_emp["49F1","site"] <- "pond_drain"
meta_emp["49G1","site"] <- "pond_drain"


sample_data(ps_relative_emp_tree) <- sample_data(meta_emp)

library(vegan)
#full models
permanova_emp_bray <- adonis2(bray_emp ~ Dissected + site/log + Instar + swapped*type_original_log, 
                               data = meta_emp, 
                               permutations = 999, 
                               method = "bray",
                               by = "terms")
print(permanova_emp_bray)


#let's check out the dispersion stuff
#check homogeneity of variances for site
disp1.1 <- betadisper(bray_emp, group=meta_emp$Instar)
permutest(disp1.1)
# p value > 0.05,  homogeneous variances 
#location but not variance change across levels of swapped*type_log_original


#visualization
library(vegan)

nmds_bray <- metaMDS(bray_emp, k = 3, trymax = 100)
nmds_bray$stress
# 0.1406287 acceptable!

#extract scores
nmds_bray_points <- as.data.frame(nmds_bray$points)
View(nmds_bray_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_bray_points$SampleID <- rownames(nmds_bray_points)
#remove scores previously added to meta_rel
View(meta_emp)
#add sample_ID column
meta_emp$Sample_ID <- rownames(meta_emp)
meta_emp <- meta_emp %>%
  rename(SampleID = Sample_ID)
sample_data(ps_relative_emp_tree) <- sample_data(meta_emp)

View(sample_data(ps_relative_emp_tree))


# Now merge by SampleID
nmds_bray_meta <- merge(nmds_bray_points, meta_emp, by = "SampleID")
View(nmds_bray_meta)

head(nmds_bray_meta)

library(scatterplot3d)
install.packages("rgl")
library(rgl)

library(scatterplot3d)
install.packages("rgl")
library(rgl)


library(plotly)

# Ensure group_factor and colors are set
group_factor <- as.factor(nmds_bray_meta$treatment)
colors <- c("dodgerblue", "orangered", "gold", "mediumspringgreen")

# Build 3D scatter plot
fig <- plot_ly(
  data = nmds_bray_meta,
  x = nmds_bray_meta[, 2],  # MDS1
  y = nmds_bray_meta[, 3],  # MDS2
  z = nmds_bray_meta[, 4],  # MDS3
  type = "scatter3d",
  mode = "markers",
  color = group_factor,
  colors = colors,
  marker = list(size = 4),
  text = ~rownames(nmds_bray_meta),
  hoverinfo = "text"
)

# Add title with stress value
stress_val <- round(nmds_bray$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Bray NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped*type_original_log"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//bray_nmds_treat_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for type_log original without hulls around groups
fig <- plot_ly(
  data = nmds_bray_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~type_original_log,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_bray_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_bray$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Bray NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "type_log_original"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//bray_nmds_log_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for swapped without hulls around groups
fig <- plot_ly(
  data = nmds_bray_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~swapped,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_bray_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_bray$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Bray NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped"))
)
# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//bray_nmds_swapped_emp.html", selfcontained = TRUE)


#ctd----

#PERMANOVA (jaccard)----

# Compute jaccard
jaccard_emp <- phyloseq::distance(ps_relative_emp_tree, method = "jaccard", binary = TRUE)


#we need meta data along with distance matrix to do pERmANOVA 
#we already have meta small and meta rel from our previous analysis
#actually meta small and meta rela should be the same data frame so we can use the interchangeably

#full models
permanova_emp_jaccard <- adonis2(jaccard_emp ~ Dissected + site/log + Instar + swapped*type_original_log, 
                                  data = meta_emp, 
                                  permutations = 999, 
                                  method = "jaccard",
                                  by = "terms")
print(permanova_emp_jaccard)
#eep no significant treatment effects


#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp2.2 <- betadisper(jaccard_emp, group=meta_emp$Instar)
permutest(disp2.2)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment



#visualization
library(vegan)
nmds_jaccard <- metaMDS(jaccard_emp, k = 3, trymax = 100)
nmds_jaccard$stress
#0.1766943 acceptable

#extract scores
nmds_jaccard_points <- as.data.frame(nmds_jaccard$points)
View(nmds_jaccard_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_jaccard_points$SampleID <- rownames(nmds_jaccard_points)
#remove scores previously added to meta_rel


# Now merge by SampleID
nmds_jaccard_meta <- merge(nmds_jaccard_points, meta_emp, by = "SampleID")
View(nmds_jaccard_meta)


library(scatterplot3d)
install.packages("rgl")
library(rgl)


library(plotly)

# Ensure group_factor and colors are set
group_factor <- as.factor(nmds_jaccard_meta$treatment)
colors <- c("dodgerblue", "orangered", "gold", "mediumspringgreen")

# Build 3D scatter plot
fig <- plot_ly(
  data = nmds_jaccard_meta,
  x = nmds_jaccard_meta[, 2],  # MDS1
  y = nmds_jaccard_meta[, 3],  # MDS2
  z = nmds_jaccard_meta[, 4],  # MDS3
  type = "scatter3d",
  mode = "markers",
  color = group_factor,
  colors = colors,
  marker = list(size = 4),
  text = ~rownames(nmds_jaccard_meta),
  hoverinfo = "text"
)

# Add title with stress value
stress_val <- round(nmds_jaccard$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Jaccard NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped*type_original_log"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//jaccard_nmds_treat_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for type_log original without hulls around groups
fig <- plot_ly(
  data = nmds_jaccard_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~type_original_log,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_jaccard_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_jaccard$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Jaccard NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "type_log_original"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//jaccard_nmds_log_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for swapped without hulls around groups
fig <- plot_ly(
  data = nmds_jaccard_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~swapped,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_jaccard_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_jaccard$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Jaccard NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped"))
)
# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//jaccard_nmds_swapped_emp.html", selfcontained = TRUE)

#ctd ----


#PERMANOVA (unifrac unweighted)----

# Compute Unweighted UniFrac
#this considers only the presence or absence of taxa.
#it measures the fraction of unique branch length in the phylogenetic tree that is exclusive to each community.
#Sensitive to rare taxa because all taxa are treated equally.
#Useful when you're interested in whether communities share the same taxa, regardless of how abundant they are.

unifrac_unweighted_rel <- phyloseq::distance(ps_relative_emp_tree, method = "unifrac", weighted = FALSE)


library(vegan)
permanova_rel_UU <- adonis2(unifrac_unweighted_rel ~ Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_emp, 
                            permutations = 999,
                            by = "terms")
print(permanova_rel_UU)
#slightly significant interaction
#not significant when I just do treatment tho

#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp4.4 <- betadisper(unifrac_unweighted_rel, group=meta_emp$Instar)
permutest(disp4.4)
# p value  =  0.028, non-homogeneous variances 
#variance changes across treatment could be contributing to significant PERMANOVA results



#visualization
library(vegan)
nmds_UU <- metaMDS(unifrac_unweighted_rel, k = 3, trymax = 100)
nmds_UU$stress
#stress = 0.1515661


#extract scores
nmds_UU_points <- as.data.frame(nmds_UU$points)
View(nmds_UU_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_UU_points$SampleID <- rownames(nmds_UU_points)
#remove scores previously added to meta_rel
#View(meta_oxy)
#meta_rel <- meta_rel[, -c(35:37)]
#meta_oxy$SampleID <- rownames(meta_oxy)

# Now merge by SampleID
nmds_UU_meta <- merge(nmds_UU_points, meta_emp, by = "SampleID")
View(nmds_UU_meta)


library(scatterplot3d)
install.packages("rgl")
library(rgl)

library(plotly)

# Ensure group_factor and colors are set
group_factor <- as.factor(nmds_UU_meta$treatment)
colors <- c("dodgerblue", "orangered", "gold", "mediumspringgreen")

# Build 3D scatter plot
fig <- plot_ly(
  data = nmds_UU_meta,
  x = nmds_UU_meta[, 2],  # MDS1
  y = nmds_UU_meta[, 3],  # MDS2
  z = nmds_UU_meta[, 4],  # MDS3
  type = "scatter3d",
  mode = "markers",
  color = group_factor,
  colors = colors,
  marker = list(size = 4),
  text = ~rownames(nmds_UU_meta),
  hoverinfo = "text"
)

# Add title with stress value
stress_val <- round(nmds_UU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Unweighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped*type_original_log"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//UU_nmds_treat_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for type_log original without hulls around groups
fig <- plot_ly(
  data = nmds_UU_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~type_original_log,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_UU_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_UU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Unweighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "type_log_original"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//UU_nmds_log_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for swapped without hulls around groups
fig <- plot_ly(
  data = nmds_UU_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~swapped,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_UU_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_UU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Unweighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped"))
)
# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//UU_nmds_swapped_emp.html", selfcontained = TRUE)


#ctd ----

#PERMANOVA (unifrac weighted)----
# Compute Weighted UniFrac
#considers both presence/absence and relative abundance of taxa.
#measures the fraction of branch length that is weighted by the difference in relative abundance between communities.
#Sensitive to abundant taxa.
#Gives more weight to dominant taxa and downplays rare taxa.

unifrac_weighted_rel <- phyloseq::distance(ps_relative_emp_tree, method = "unifrac", weighted = TRUE)


library(vegan)
#let's look at the control effects first

#full models
permanova_rel_WU <- adonis2(unifrac_weighted_rel ~  Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_emp, 
                            permutations = 999,
                            by = "terms")

print(permanova_rel_WU)

#this is fascinating. the relative abundance of lineages aren't really changing
#even tho the pressence of lineages
#relative abundance of taxa
#and pressence of taxa is changing
#this is telling me that even tho there is change of taxa and lineages present,
#major clade level restructuring of communities is not going on



#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp3.3 <- betadisper(unifrac_weighted_rel, group=meta_emp$type_original_log)
permutest(disp3.3)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment


#visualization
library(vegan)
nmds_WU <- metaMDS(unifrac_weighted_rel, k = 3, trymax = 100, verbose = FALSE)
nmds_WU$stress
#0.07380359 #really good!


#extract scores
nmds_WU_points <- as.data.frame(nmds_WU$points)

#merge with metadata for plotting
# First, add the rownames as a column in each data frame

nmds_WU_points$SampleID <- rownames(nmds_WU_points)

# Now merge by SampleID
nmds_WU_meta <- merge(nmds_WU_points, meta_emp, by = "SampleID")
View(nmds_WU_meta)


#visualization
library(plotly)
# Ensure group_factor and colors are set
group_factor <- as.factor(nmds_WU_meta$treatment)
colors <- c("dodgerblue", "orangered", "gold", "mediumspringgreen")

# Build 3D scatter plot
fig <- plot_ly(
  data = nmds_WU_meta,
  x = nmds_WU_meta[, 2],  # MDS1
  y = nmds_WU_meta[, 3],  # MDS2
  z = nmds_WU_meta[, 4],  # MDS3
  type = "scatter3d",
  mode = "markers",
  color = group_factor,
  colors = colors,
  marker = list(size = 4),
  text = ~rownames(nmds_WU_meta),
  hoverinfo = "text"
)

# Add title with stress value
stress_val <- round(nmds_WU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Weighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped*type_original_log"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//WU_nmds_treat_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for type_log original without hulls around groups
fig <- plot_ly(
  data = nmds_WU_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~type_original_log,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_WU_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_WU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Weighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "type_log_original"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//WU_nmds_log_emp.html", selfcontained = TRUE)



library(plotly)
#plot 3D NMDS ordination for swapped without hulls around groups
fig <- plot_ly(
  data = nmds_WU_meta,
  x = ~MDS1,  # replace with actual column names if different
  y = ~MDS2,
  z = ~MDS3,
  type = "scatter3d",
  mode = "markers",
  color = ~swapped,  # use column directly as formula
  colors = c("dodgerblue", "orangered"),
  marker = list(size = 4),
  text = ~rownames(nmds_WU_meta),
  hoverinfo = "text"
)
#nmds_bray_meta[,2], nmds_bray_meta[,3], nmds_bray_meta[,4],

# Optional: Add custom legend titles
# Add title with stress value
stress_val <- round(nmds_WU$stress, 4)
fig <- fig %>% layout(
  title = list(
    text = paste("Weighted Unifrac NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "swapped"))
)
# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//WU_nmds_swapped_emp.html", selfcontained = TRUE)


#ctd----

#multiple tests corrections####
#first we need to adjust global p-values across distance metrics
# I'm testing the same hypothesis (treatment effect) under 4 lenses
#bray-curtis, jaccard, weighted unifrac, unweighted unifrac,


# Named list of your models
models <- list(
  bray = permanova_emp_bray,
  jaccard = permanova_emp_jaccard,
  w_unifrac = permanova_rel_WU,
  u_unifrac = permanova_rel_UU
)

# Extract results from each and add term + metric labels
results_list <- lapply(names(models), function(metric) {
  df <- as.data.frame(models[[metric]])
  df$term <- rownames(df)
  rownames(df) <- NULL
  df$metric <- metric
  df
})

# Combine into one data frame
results_df <- do.call(rbind, results_list)

# Reorder columns for clarity
results_df <- results_df[, c("metric", "term", "Df", "SumOfSqs", "R2", "F", "Pr(>F)")]
View(results_df)


#apply fdr adjustment
results_df$p_adj_bonf <- ave(results_df$`Pr(>F)`, results_df$term, FUN = p.adjust, method = "bonferroni")
View(results_df)


#First, extract the 6 pairwise treatment comparisons:
meta_emp$treatment_paste <- paste(meta_emp$swapped, meta_emp$type_original_log)
View(meta_emp)
library(vegan)

treatment_levels <- unique(meta_emp$treatment_paste)
print(treatment_levels)
pairwise_comparisons <- combn(treatment_levels, 2, simplify = FALSE)
print(pairwise_comparisons)


#Define a Function to Run PERMANOVA per Pair and Metric
#Subset meta_oxy to just the two treatment levels.
#Subsets the corresponding distance matrix.
#Runs PERMANOVA with adonis2, including your control variables.
#Returns the p-value for treatment.

library(vegan)

run_pairwise_permanova <- function(dist_matrix, meta, treatment_var, control_vars, pair) {
  # -------------------------------
  # Function to run a PERMANOVA for one pairwise comparison of treatment levels
  # Arguments:
  #   dist_matrix   = full distance matrix (e.g., Bray, UniFrac, etc.)
  #   meta          = metadata dataframe with sample info
  #   treatment_var = name of the treatment column (as a string)
  #   control_vars  = character vector of control variable names (e.g., c("Dissected", "site/log", "Instar"))
  #   pair          = vector of 2 treatment levels to compare (e.g., c("A", "B"))
  # Returns:
  #   Raw p-value for treatment effect in PERMANOVA
  
  # Subset the metadata to only include samples from the two treatment levels in the current pair
  subset_meta <- meta[meta[[treatment_var]] %in% pair, ]
  
  # Get sample IDs that correspond to those rows
  sample_ids <- rownames(subset_meta)
  
  # Subset the distance matrix to only include those samples
  # (This assumes the row and column names of the distance matrix match the rownames of meta)
  sub_dist <- as.matrix(dist_matrix)[sample_ids, sample_ids]
  
  # Convert back to a 'dist' object (required by adonis2)
  sub_dist <- as.dist(sub_dist)
  
  # Construct the formula for adonis2
  # Combine control variables into a formula string: "Dissected + site/log + Instar"
  control_formula <- paste(control_vars, collapse = " + ")
  
  # Create the full formula string, e.g., "sub_dist ~ Dissected + site/log + Instar + treatment"
  full_formula_str <- paste("sub_dist ~", control_formula, "+", treatment_var)
  
  # Convert to a formula object
  formula <- as.formula(full_formula_str)
  
  # Run PERMANOVA using adonis2 with 999 permutations
  # We use `by = "terms"` to get individual p-values for each term
  mod <- adonis2(formula, data = subset_meta, permutations = 999, by = "terms")
  
  # Extract the p-value for the treatment term only
  # This assumes the row name in the result exactly matches the treatment_var name
  p_val <- mod[rownames(mod) == treatment_var, "Pr(>F)"]
  
  # Return the raw p-value
  return(p_val)
}



#apply this function across all community metrics 
# Run pairwise comparisons

# Prepare output holders
results_list <- list()

# Controls and distance matrices
control_vars <- c("Dissected", "site/log", "Instar")
distance_matrices <- list( #list only the distance matrices with significant interaction effects because these are the ones we need pairwise comparisons for
  bray = bray_emp,
  jaccard = jaccard_emp,
  u_unifrac = unifrac_unweighted_rel
)




# Run pairwise comparisons
for (metric in names(distance_matrices)) {
  dist_matrix <- distance_matrices[[metric]]
  
  res <- lapply(pairwise_comparisons, function(pair) {
    p <- run_pairwise_permanova(dist_matrix, meta_emp, "treatment_paste", control_vars, pair)
    
    # Inspect p in more detail
    cat("\n--- Metric:", metric, " | Pair:", pair[1], "vs", pair[2], "---\n")
    cat("Class of p:", class(p), "\n")
    cat("Length of p:", length(p), "\n")
    print(p)  # Will show full object
    
    # If p isn't a scalar numeric, skip
    if (!is.numeric(p) || length(p) != 1 || is.na(p)) {
      cat("⚠️ Skipping due to invalid p\n")
      return(NULL)
    }
    
    # Return clean data frame
    data.frame(
      metric = metric,
      group1 = pair[1],
      group2 = pair[2],
      p_value = p
    )
  })
  
  results_list[[metric]] <- do.call(rbind, res)
}




#Now pull all results together and apply bonferroni correction within each metric:

View(results_list[1])
View(results_list[2])
View(results_list[3])

# Combine all into one data frame
pairwise_results <- do.call(rbind, results_list)
View(pairwise_results)

# bonferroni adjustment within each metric
pairwise_results$p_adj_bonf <- ave(pairwise_results$p_value, pairwise_results$metric, FUN = p.adjust, method = "bonferroni")
View(pairwise_results)

#replace treatments with more intuitive labels
# Create a lookup named vector: names are the treatment_paste values, values are the actual treatments
lookup <- setNames(meta_emp$treatment, meta_emp$treatment_paste)

# Replace group1 and group2 using the lookup
pairwise_results$group1 <- lookup[pairwise_results$group1]
pairwise_results$group2 <- lookup[pairwise_results$group2]
View(pairwise_results)

# Function to round all numeric columns to 4 significant figures
round_sigfig <- function(df, digits = 4) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) signif(x, digits) else x
  })
  return(df)
}

# Apply rounding to both data frames
results_df <- round_sigfig(results_df)
pairwise_results <- round_sigfig(pairwise_results)


# Function to format p-values with asterisks (no separate column)
add_sig_stars <- function(p, digits = 4) {
  p_rounded <- signif(p, digits)
  stars <- ifelse(p < 0.001, "***",
                  ifelse(p < 0.01, "**",
                         ifelse(p < 0.05, "*", "")))
  return(paste0(p_rounded, stars))
}

# Apply to PERMANOVA summary sheet
if ("p_value" %in% names(results_df)) {
  results_df$p_value <- add_sig_stars(results_df$p_value)
}
if ("p_adj_fdr" %in% names(results_df)) {
  results_df$p_adj_fdr <- add_sig_stars(results_df$p_adj_fdr)
}
if ("p_adj_bonf" %in% names(results_df)) {
  results_df$p_adj_bonf <- add_sig_stars(results_df$p_adj_bonf)
}

# Apply to Pairwise Results sheet
if ("p_value" %in% names(pairwise_results)) {
  pairwise_results$p_value <- add_sig_stars(pairwise_results$p_value)
}
if ("p_adj_fdr" %in% names(pairwise_results)) {
  pairwise_results$p_adj_fdr <- add_sig_stars(pairwise_results$p_adj_fdr)
}
if ("p_adj_bonf" %in% names(pairwise_results)) {
  pairwise_results$p_adj_bonf <- add_sig_stars(pairwise_results$p_adj_bonf)
}

View(results_df)
# Blank out repeated metric values after the first occurrence
results_df$metric <- as.character(results_df$metric)
results_df$metric[duplicated(results_df$metric)] <- ""


pairwise_results$metric <- as.character(pairwise_results$metric)
pairwise_results$metric[duplicated(pairwise_results$metric)] <- ""

# Replace "NANA" strings in p-value columns with blank
cols_to_clean <- c("p_value", "p_adj_bonf")

for (col in cols_to_clean) {
  if (col %in% names(results_df)) {
    results_df[[col]][results_df[[col]] == "NANA"] <- ""
  }
}

#concaonate group columns to make a contrast column for pairwise comparisons
pairwise_results$Contrast_treatment <- paste(pairwise_results$group1, pairwise_results$group2, sep = "-")
View(pairwise_results)
pairwise_results <- pairwise_results %>% select(-group1, -group2)
pairwise_results <- pairwise_results[, c(1, ncol(pairwise_results), 2:(ncol(pairwise_results)-1))]
View(pairwise_results)

print(results_df)
print(class(results_df))

print(pairwise_results)
print(class(pairwise_results))

print(summary_table)
print(class(summary_table))


#create summary of results_df for results section
#create summary of results_df for results section
# Fill in the metric column
results_df_filled <- results_df %>%
  mutate(
    metric = trimws(as.character(metric)),
    term = trimws(as.character(term)),
    metric = ifelse(metric == "", NA, metric)
  ) %>%
  tidyr::fill(metric)

# Select variables and metrics of interest
vars <- c("type_original_log", "swapped", "swapped:type_original_log")
metrics <- c("bray", "jaccard", "w_unifrac", "u_unifrac")

# Filter and format combined string of R2 and p-value
filtered_df <- results_df_filled %>%
  filter(term %in% vars, metric %in% metrics) %>%
  mutate(
    R2 = round(R2, 3),
    pval_str = paste0("R² = ", R2, ", p = ", p_adj_bonf)
  ) %>%
  select(metric, term, pval_str)

# Pivot to wide format
summary_table <- filtered_df %>%
  pivot_wider(names_from = metric, values_from = pval_str) %>%
  column_to_rownames("term")

# Print result
print(summary_table)
View(summary_table)

print(summary_table)
print(class(summary_table))

print(results_df)
print(class(results_df))

print(pairwise_results)
print(class(pairwise_results))

summary_table <- summary_table %>%
  tibble::rownames_to_column(var = "treatment")

# 1. Round numeric columns in results_df to 2 significant figures
results_df <- results_df %>%
  mutate(across(where(is.numeric), ~ signif(.x, 2)))

# 2. Rename dissimilarity metrics in all 3 dataframes
metric_names <- c(
  "bray" = "Bray Curtis",
  "jaccard" = "Jaccard",
  "w_unifrac" = "weighted Unifrac",
  "u_unifrac" = "unweighted Unifrac"
)

results_df <- results_df %>%
  mutate(metric = dplyr::recode(metric, !!!as.list(metric_names)))

pairwise_results <- pairwise_results %>%
  mutate(metric = dplyr::recode(metric, !!!as.list(metric_names)))

colnames(summary_table) <- dplyr::recode(colnames(summary_table), !!!as.list(metric_names))

# 3. Rename terms in results_df
results_df <- results_df %>%
  mutate(term = case_when(
    term == "swapped" ~ "diet change",
    term == "type_original_log" ~ "original diet",
    term == "swapped:type_original_log" ~ "diet change x original diet",
    term == "Dissected" ~ "block",
    term == "site:log" ~ "site:original diet",
    TRUE ~ term
  ))


# 4. Replace treatment column names in summary_table
summary_table$treatment <- c(
  "diet change",
  "original diet",
  "diet change x original diet"
)

# 5. Rename contrasts in pairwise_results (change : to > within treatments)
pairwise_results <- pairwise_results %>%
  mutate(Contrast_treatment = gsub(":", ">", Contrast_treatment))

#reorder rows of summary table
# Desired order vector
desired_order <- c("original diet", "diet change", "diet change x original diet")

# Reorder rows by matching the treatment column to desired_order
summary_table <- summary_table[match(desired_order, summary_table$treatment), ]

print(summary_table)

# Print to verify
print(results_df)
print(pairwise_results)
print(summary_table)

#export both tables
# Install if needed
install.packages("openxlsx")

# Load package
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Add sheets
addWorksheet(wb, "PERMANOVA Summary")
writeData(wb, "PERMANOVA Summary", results_df)

addWorksheet(wb, "Pairwise Results")
writeData(wb, "Pairwise Results", pairwise_results)

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_table)

# Save workbook to working directory
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//permanova_all_results_emp.xlsx", overwrite = TRUE)
#ctd----

#diversity analysis----
library(phyloseq)

# Calculate shannon, and simpsons diversity for each samples
diversity_emp <- estimate_richness(ps_relative_emp, measures = c("Shannon", "Simpson"))
View(diversity_emp)


#this code cannot calculate richness - 
#it get's tripped up by relative abundances rather than raw counts
#throws this warning if I try to include richness
#Warning message:
#The data you have provided does not have any singletons. This is highly suspicious. Results of richness
#estimates (for example) are probably unreliable, or wrong, if you have already
#trimmed low-abundance taxa from the data. We recommended that you find the un-trimmed data and retry.

#so let's calculate richness manually by counting asvs with abundance > 0 for each sample

# Extract the ASV/OTU table (samples as rows)
otu_emp <- as(otu_table(ps_relative_emp_tree), "matrix")

# Make sure samples are rows and transpose if needed
if (!taxa_are_rows(ps_relative_emp_tree)) {
  otu_emp <- t(otu_emp)
}

View(otu_emp)

# Count number of ASVs with abundance > 0 per sample
richness_emp <- colSums(otu_emp > 0)

# Look at the result
head(richness_emp)



# Add diversity metrics to sample_data
library(phyloseq)

sample_data(ps_relative_emp_tree) <- cbind(sample_data(ps_relative_emp_tree), diversity_emp)
sample_data(ps_relative_emp_tree) <- cbind(sample_data(ps_relative_emp_tree), data.frame(richness_emp))
View(sample_data(ps_relative_emp_tree))

#edits treatment labels
# Extract the sample_data as a data.frame
sample_data_df <- as.data.frame(sample_data(ps_relative_emp_tree))
View(sample_data_df)
# Modify the 'treatment' column
sample_data_df$treatment <- gsub("wood", "", sample_data_df$treatment)
sample_data_df$type_original_log <- gsub("wood", "", sample_data_df$type_original_log)
sample_data_df$treatment <- gsub("_", ":", sample_data_df$treatment)
# Reassign the modified sample_data back into the phyloseq object
sample_data(ps_relative_emp_tree) <- sample_data(sample_data_df)
View(sample_data(ps_relative_emp_tree))

#plot raw diversity metrics by treatment
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract sample metadata from phyloseq object
meta_emp <- data.frame(sample_data(ps_relative_emp_tree))
View(meta_emp)

#rename richness column
meta_emp <- meta_emp %>%
  rename(Richness = richness_emp)

#add SampleID column
meta_emp$SampleID <- rownames(meta_emp)

library(dplyr)
library(tidyr)
#pivot data to long format
meta_long <- meta_emp %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness),
               names_to = "Metric",
               values_to = "Value")
View(meta_long)

# Make three separate boxplots, one for each grouping variable
p1 <- ggplot(meta_long, aes(x = treatment, y = Value, fill = treatment)) +
  geom_boxplot(color = "black") +
  facet_wrap(~ Metric, nrow = 3, scales = "free_y") +
  scale_fill_grey(start = 0.3, end = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Treatment", y = "Diversity", title = "By Treatment")

p2 <- ggplot(meta_long, aes(x = type_original_log, y = Value, fill = type_original_log)) +
  geom_boxplot(color = "black") +
  facet_wrap(~ Metric, nrow = 3, scales = "free_y") +
  scale_fill_grey(start = 0.3, end = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Original Log Type", y = "Diversity", title = "By Original Log Type")

p3 <- ggplot(meta_long, aes(x = swapped, y = Value, fill = swapped)) +
  geom_boxplot(color = "black") +
  facet_wrap(~ Metric, nrow = 3, scales = "free_y") +
  scale_fill_grey(start = 0.3, end = 0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Swapped", y = "Diversity", title = "By Swapped")

# Arrange into a 3x3 grid
library(patchwork)

(p1 | p2 | p3) + plot_layout(guides = "collect") & theme(legend.position = "right")

#there seem to be some differences in means by treatment for all diversity metrics

#now let's do the glms to test for significant differences between treatments

#establish assumptions test function
library(car)
library(DHARMa)
library(glmmTMB)
dharma.test<-function(x,n=1000){
  q <- simulateResiduals(x, n)
  par(mfrow = c(2, 2), # 2 x 2 pictures on one plot
      pty = "s")
  plotQQunif(q)
  UR <- recordPlot()
  plotResiduals(q)
  UL <- recordPlot()
  print(testDispersion(q))
  LL <- recordPlot()
  print(testZeroInflation(q))
  LR <- recordPlot()
  par(mfrow=c(1,1))
}


#test simpsons

hist(meta_emp$Simpson)
#data is a little right skewed
meta_emp$log_Simpsons <- log(meta_emp$Simpson)
hist(meta_emp$log_Simpson)

#Gamma(link = "log")

#trying a gamma distribution because chat gpt said it was good for positive right skewed data
#Gaussian (normal) distribution didn't pass the dHARMA tests

glm_simp_emp <- glmmTMB(as.numeric(Simpson) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_emp)
summary(glm_simp_emp) 
aov_simp_emp <- Anova(glm_simp_emp) 
aov_simp_emp

length(unique(meta_emp$site))
length(unique(meta_emp$log))
length(unique(meta_emp$Dissected))
length(unique(meta_emp$Instar))
length(unique(meta_emp$swapped))
length(unique(meta_emp$type_original_log))

nrow(meta_emp)

dharma.test(glm_simp_emp,n=1000)
#assumptions test looks good!
# significant effect of swapping. interesting

#pairwise comparisons
pairwise_simp_emp_instar <- emmeans(glm_simp_emp, pairwise ~ Instar, adjust = "tukey")
pairwise_simp_emp_instar$contrasts

#4 and A are barely different but the others aren't

pairwise_simp_emp_swap <- emmeans(glm_simp_emp, pairwise ~ swapped, adjust = "tukey")
pairwise_simp_emp_swap$contrasts
#"no" treatment has lower Simpson values → higher diversity
#"yes" treatment has higher Simpson values → lower diversity

#test shannon

hist(meta_emp$Shannon)
#looks super normal

glm_shan_emp <- glmmTMB(as.numeric(Shannon) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_emp)
summary(glm_shan_emp) 
aov_shan_emp <- Anova(glm_shan_emp) 
aov_shan_emp

dharma.test(glm_shan_emp,n=1000)
#assumptions test looks good!
# effects of swapping and log


#pairwise comparisons
pairwise_shan_emp_instar <- emmeans(glm_shan_emp, pairwise ~ Instar, adjust = "tukey")
pairwise_shan_emp_instar$contrasts

#again 4 and A are barely different but the others aren't

pairwise_shan_emp_swap <- emmeans(glm_shan_emp, pairwise ~ swapped, adjust = "tukey")
pairwise_shan_emp_swap$contrasts
#Simpson index says: "no" treatment has higher diversity.
#Shannon index says: "yes" treatment has higher diversity.

#Shannon is more sensitive to richness and rare species.
#Simpson puts more weight on dominance or common species.

#This suggests:
#The "no" treatment community likely has fewer dominant species and fewer rare species — 
#maybe a more even community but with less richness overall.
#The "yes" treatment community probably has more rare species, 
#increasing Shannon diversity, but also some dominant species that lower Simpson diversity.

#"No" treatment: community is more even but less rich.
#"Yes" treatment: community is richer with more rare species but also has dominant species.

#test richness

hist(meta_emp$Richness)
#looks normal

glm_rich_emp <- glmmTMB(as.numeric(Richness) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_emp)
summary(glm_rich_emp) 
aov_rich_emp <- Anova(glm_rich_emp) 
aov_rich_emp

#hmm interesting
#slight effect of original log
#huuuuuge effect of interaction

dharma.test(glm_rich_emp,n=1000)
#assumptions test looks good!
# so no significant effects 


#pairwise comparisons
pairwise_rich_emp_original <- emmeans(glm_rich_emp, pairwise ~ type_original_log, adjust = "tukey")
pairwise_rich_emp_original$contrasts

#on average, softwood has 27 more taxa than hardwood. interesting!

#not different in the contrasts comparison

pairwise_rich_emp_treat <- emmeans(glm_rich_emp, pairwise ~ swapped*type_original_log, adjust = "tukey")
pairwise_rich_emp_treat$contrasts
# Extract the contrasts table
contrast_df <- as.data.frame(pairwise_rich_emp_treat$contrasts)

# Define replacement mapping
rename_treatment <- function(label) {
  label <- str_replace_all(label, "no hard", "hard:hard")
  label <- str_replace_all(label, "yes hard", "hard:soft")
  label <- str_replace_all(label, "no soft", "soft:soft")
  label <- str_replace_all(label, "yes soft", "soft:hard")
  return(label)
}

library(stringr)
# Apply the renaming
contrast_df$contrast <- rename_treatment(contrast_df$contrast)

# View updated table
pairwise_rich_emp_treat$contrasts <- contrast_df
pairwise_rich_emp_treat$contrasts

#it seems there are more softwood associated taxa than hardwood
#and that they can be regained? I need to wrap my head around  this more


library(glmmTMB)
library(car) # for Anova()
library(openxlsx)
library(dplyr)
library(broom) # for tidying model output
library(tibble)

# Helper function to format ANOVA tables with a model heading
prepare_anova_df <- function(anova_obj, sigfigs = 2) {
  df <- as.data.frame(anova_obj) %>%
    tibble::rownames_to_column("rowname")
  
  # Apply signif() to all numeric columns
  df <- df %>%
    mutate(across(where(is.numeric), ~ signif(.x, sigfigs)))
  
  # Add significance stars to p-values
  if ("Pr(>Chisq)" %in% names(df)) {
    df <- df %>%
      mutate(
        sig_stars = case_when(
          `Pr(>Chisq)` < 0.001 ~ "***",
          `Pr(>Chisq)` < 0.01 ~ "**",
          `Pr(>Chisq)` < 0.05 ~ "*",
          `Pr(>Chisq)` < 0.1 ~ ".",
          TRUE ~ ""
        ),
        `Pr(>Chisq)` = paste0(`Pr(>Chisq)`, sig_stars)
      ) %>%
      select(-sig_stars)
  }
  
  return(df)
}

# Helper function to rename variables in formula string
rename_formula_vars <- function(formula_str) {
  formula_str %>%
    gsub("Dissected", "block", .) %>%
    gsub("swapped", "diet change", .) %>%
    gsub("type_original_log", "original diet", .) %>%
    gsub("\\*", " x ", .)  # Optional: make interaction operator more readable
}

# Define models, ANOVA objects, and formulas in a list
anova_models <- list(
  Simpson = list(
    formula = "Simpson ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_simp_emp
  ),
  Shannon = list(
    formula = "Shannon ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_shan_emp
  ),
  Richness = list(
    formula = "Richness ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_rich_emp
  )
)

# Create Excel workbook
wb <- createWorkbook()

# Loop over each model to add sheets
for (name in names(anova_models)) {
  addWorksheet(wb, name)
  
  # Rename variables in formula and write at the top
  renamed_formula <- rename_formula_vars(anova_models[[name]]$formula)
  writeData(wb, name, renamed_formula, startRow = 1, startCol = 1)
  
  # Prepare and write the ANOVA table starting at row 3
  anova_df <- prepare_anova_df(anova_models[[name]]$table)
  
  # Rename variable labels for readability in table
  anova_df$rowname <- anova_df$rowname %>%
    dplyr::recode(
      "swapped" = "diet change",
      "type_original_log" = "original diet",
      "swapped:type_original_log" = "diet change x original diet"
    )
  
  writeData(wb, name, anova_df, startRow = 3, startCol = 1)
}

# Save the workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//diversity_stats_tables_emp.xlsx", overwrite = TRUE)

#export contrasts
#contrasts:
#pairwise_simp_emp_swap$contrasts
#pairwise_shan_emp_swap$contrasts
#pairwise_rich_emp_treat$contrasts
#pairwise_rich_emp_original$contrasts

contrast_tables <- list(
  Simpson_diet_change = pairwise_simp_emp_swap$contrasts,
  Shannon_diet_change = pairwise_shan_emp_swap$contrasts,
  Richness_change_x_original = pairwise_rich_emp_treat$contrasts,
  Richness_original_diet = pairwise_rich_emp_original$contrasts
)

# Create a new workbook
wb <- createWorkbook()

for (name in names(contrast_tables)) {
  addWorksheet(wb, name)
  
  # Process the data frame
  contrast_df <- as.data.frame(contrast_tables[[name]]) %>%
    # Round all numeric columns except p.value
    mutate(across(where(is.numeric) & !matches("p\\.value"), ~ signif(.x, 2))) %>%
    
    # Format p.value separately and append asterisks
    mutate(p.value = case_when(
      p.value < 0.001 ~ paste0(signif(p.value, 2), " ***"),
      p.value < 0.01  ~ paste0(signif(p.value, 2), " **"),
      p.value < 0.05  ~ paste0(signif(p.value, 2), " *"),
      TRUE            ~ as.character(signif(p.value, 2))
    ))
  
  # For the Richness interaction contrasts, replace : with > in contrast column
  if (name == "Richness_change_x_original" && "contrast" %in% names(contrast_df)) {
    contrast_df$contrast <- gsub(":", ">", contrast_df$contrast)
  }
  
  writeData(wb, name, contrast_df)
}

# Save the workbook
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//diversity_contrasts_tables_emp.xlsx", overwrite = TRUE)

#plot emmeans from all models for visualizing pairwise contrasts
library(emmeans)
library(patchwork)
library(glmmTMB) 
library(dplyr)

#Get EMMs for each metric and predictor combination 

get_emm_df <- function(model, variable, metric_name) {
  var_formula <- if (variable == "treatment") {
    ~ swapped * type_original_log
  } else {
    reformulate(variable)
  }
  
  # Get EMMs
  emm <- emmeans(model, var_formula, infer = c(TRUE, TRUE)) %>%
    as.data.frame()
  
  # Fix confidence limits if needed
  if (!"lower.CL" %in% names(emm)) emm$lower.CL <- emm$asymp.LCL
  if (!"upper.CL" %in% names(emm)) emm$upper.CL <- emm$asymp.UCL
  
  # Add labels for plotting
  emm <- emm %>%
    mutate(
      Metric = metric_name,
      Variable = variable
    )
  
  # Add Level column safely
  if (variable == "treatment") {
    emm <- emm %>%
      mutate(
        Level = case_when(
          type_original_log == "hard" & swapped == "no"  ~ "hard>hard",
          type_original_log == "hard" & swapped == "yes" ~ "hard>soft",
          type_original_log == "soft" & swapped == "no"  ~ "soft>soft",
          type_original_log == "soft" & swapped == "yes" ~ "soft>hard",
          TRUE ~ NA_character_
        ),
        Level = factor(Level, levels = c("hard>hard", "soft>soft", "hard>soft", "soft>hard"))
      )
  } else {
    emm <- emm %>%
      mutate(
        Level = as.factor(.[[1]])
      )
  }
  
  return(emm)
}

# Get EMMs for all combos
emm_list <- list(
  get_emm_df(glm_simp_emp, "treatment", "Simpson"),
  get_emm_df(glm_simp_emp, "swapped", "Simpson"),
  get_emm_df(glm_simp_emp, "type_original_log", "Simpson"),
  
  get_emm_df(glm_shan_emp, "treatment", "Shannon"),
  get_emm_df(glm_shan_emp, "swapped", "Shannon"),
  get_emm_df(glm_shan_emp, "type_original_log", "Shannon"),
  
  get_emm_df(glm_rich_emp, "treatment", "Richness"),
  get_emm_df(glm_rich_emp, "swapped", "Richness"),
  get_emm_df(glm_rich_emp, "type_original_log", "Richness")
)

# Combine into one dataframe
emm_all <- bind_rows(emm_list)

# Rename variables for plotting
emm_all <- emm_all %>%
  mutate(Variable = case_when(
    Variable == "treatment" ~ "diet change x original diet",
    Variable == "type_original_log" ~ "original diet",
    Variable == "swapped" ~ "diet change",
    TRUE ~ Variable
  ))

# Plotting function for one predictor variable
library(ggplot2)

library(stringr)

# Function to plot one predictor
plot_emm <- function(df, metric, variable) {
  
  # Titles without diversity metric, rename "diet change" to "diet swap"
  title_text <- switch(variable,
                       "diet change x original diet" = "diet swap x original diet",  # <-- updated
                       "diet change" = "diet swap",                                   # <-- updated
                       "original diet" = "original diet",
                       variable
  )
  
  df %>%
    filter(Metric == metric, Variable == variable) %>%
    ggplot(aes(x = Level, y = emmean)) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.2,
      color = "black",
      na.rm = TRUE
    ) +
    theme_bw() +
    labs(
      title = title_text,
      y = metric,
      x = NULL
    ) +
    theme(
      plot.title.position = "panel",           # <-- title stays inside the panel
      axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14, margin = margin(r = 10)),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      # Title styling and position
      plot.title = element_text(
        size = 14,
        hjust = 0.25,                         # <-- nudges title slightly right from panel left edge
        margin = margin(t = 5, b = 5)
      ),
      
      # Optional: panel margin inside each plot (vertical spacing)
      plot.margin = margin(t = 5, r = 5, b = 5, l = 1)
    )
}

# Relative buffer function (percentage of the data range)
get_ylim_with_relative_buffer <- function(df, metric, variable, buffer_pct = 0.1) {
  sub_df <- df %>% filter(Metric == metric, Variable == variable)
  max_upper <- max(sub_df$upper.CL, na.rm = TRUE)
  min_lower <- min(sub_df$lower.CL, na.rm = TRUE)
  
  buffer <- (max_upper - min_lower) * buffer_pct
  c(min_lower, max_upper + buffer)
}

# Plot order: original diet, diet change, interaction
p1 <- plot_emm(emm_all, "Simpson", "original diet")
p2 <- plot_emm(emm_all, "Simpson", "diet change")
p3 <- plot_emm(emm_all, "Simpson", "diet change x original diet")


p4 <- plot_emm(emm_all, "Shannon", "original diet")
p5 <- plot_emm(emm_all, "Shannon", "diet change")

p6 <- plot_emm(emm_all, "Shannon", "diet change x original diet")

p7 <- plot_emm(emm_all, "Richness", "original diet")
p8 <- plot_emm(emm_all, "Richness", "diet change")
p9 <- plot_emm(emm_all, "Richness", "diet change x original diet")

# Apply dynamic y-limits with relative buffer
plots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9)
metrics <- rep(c("Simpson", "Shannon", "Richness"), each = 3)
variables <- c("original diet", "diet change", "diet change x original diet")

for (i in seq_along(plots)) {
  plots[[i]] <- plots[[i]] + coord_cartesian(
    ylim = get_ylim_with_relative_buffer(
      emm_all,
      metrics[i],
      variables[(i-1) %% 3 + 1],
      buffer_pct = 0.1  # 10% of data range above tallest error bar
    )
  )
}

# Reassign plots
p1 <- plots[[1]]; p2 <- plots[[2]]; p3 <- plots[[3]]
p4 <- plots[[4]]; p5 <- plots[[5]]; p6 <- plots[[6]]
p7 <- plots[[7]]; p8 <- plots[[8]]; p9 <- plots[[9]]

# Combine plots with patchwork
final_plot <- (
  ((p1 | p2 | p3) + plot_layout(widths = c(1, 1, 3))) /
    ((p4 | p5 | p6) + plot_layout(widths = c(1, 1, 3))) /
    ((p7 | p8 | p9) + plot_layout(widths = c(1, 1, 3)))
) +
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "",
    tag_suffix = ")"
  ) & 
  theme(
    plot.tag = element_text(face = "bold", size = 14, hjust = 0),
    plot.tag.position = c(0.075, 0.99),  # ← moved right
    plot.margin = margin(t = 5, r = 15, b = 5, l = 2)
  )
# Save final plot
ggsave(
  "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//diversity_plot.png",
  final_plot,
  width = 7,
  height = 10,
  dpi = 600
)
#so these don't look different visually even tho some of them are
#according to chat gpt, it may be better to plot pairwise contrasts with their associated confidence intervals
#rather than plotting emmeans

#so I may do that eventually for clearer visualization
#but I'll need to do it for my para and oxy plots as well
#and I don't have time for that now so I'll need to come back to it.

#ctd----


#processing raw metadata ----
#!!!! not emp specific, should probably live in a separate file
habitat_swap_meta_data <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//habitat_swap_meta_data.xlsx")
View(habitat_swap_meta_data)

#change 45C3 to instar 3 (it was initially called instar 2 but later reevaluated and included in the experiment)
habitat_swap_methods <- data.frame(habitat_swap_meta_data)
habitat_swap_methods <- habitat_swap_methods %>%
  mutate(Instar = ifelse(Sample == "45C3", 3, Instar))
View(habitat_swap_methods)

#remove anything with instar = 2
habitat_swap_methods <- habitat_swap_methods %>%
  filter(Instar != 2)
View(habitat_swap_methods)

#add weight at death
habitat_swap_methods$death_weight <- NA

library(dplyr)
habitat_swap_methods<- habitat_swap_methods %>%
  mutate(death_weight = case_when(
    Dissected == "07_12_24" ~ weight_07_05_24,
    Dissected == "07_19_24"~ weight_07_19_24,
    Dissected == "07_26_24"~ weight_07_26_24,
    Dissected == "07_27_24"~ weight_07_27_24
  ))

habitat_swap_methods$site <- NA #make site column
colnames(habitat_swap_methods)[1] <- "log" #fix the name of the log column, it didn't read in right for some reason

#read in log site key so I can put a site variable in the data set
library(readxl)
sites <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Log_site_key.xlsx")
View(sites)
colnames(sites) <- c("log", "site") #rename sites columns so they start with lowercase  letters and match sample data

habitat_swap_methods <- habitat_swap_methods %>%
  left_join(sites %>% select(log, site), by = "log") %>%
  mutate(site = ifelse(is.na(site.x), site.y, site.x)) %>%
  select(-site.x, -site.y)
View(habitat_swap_methods)

#remove anything with death weight = NA (means it was never dissected and died during the experiment)
habitat_swap_methods <- habitat_swap_methods %>%
  filter(!is.na(death_weight))
View(habitat_swap_methods)

#remove anything with instar = 2
habitat_swap_methods <- habitat_swap_methods %>%
  filter(Instar != 2)
View(habitat_swap_methods)

#remove badswap
#remove individuals accidentally swapped onto hardwood
samples_to_remove <- c("41C1", "41F1", "41F2")  
habitat_swap_methods <- habitat_swap_methods %>%
  filter(!Sample %in% samples_to_remove)


#remove all but 4 of the weird grouped instar 3
#families 56A and 57C
#keep the first one in the family where type_log_new = "softwood" 
#and the first one where type_log_new = "hardwood"

#do 56A first
# Identify the rows to keep (first softwood and first hardwood for family_ID == "56A")
rows_to_keep <- habitat_swap_methods %>%
  filter(Family_ID == "56A") %>%
  group_by(type_log_new) %>%
  slice_head(n = 1) %>%
  ungroup()

# Filter out all other 56A rows and then add the two you want back in
habitat_swap_methods <- habitat_swap_methods %>%
  filter(Family_ID != "56A") %>%       # remove all 56A
  bind_rows(rows_to_keep)              # add back the two to keep

View(habitat_swap_methods)
#note, it added them back to the bottom of the df so they are in a differnt order now

#now do the same for 57C
# Identify the rows to keep (first softwood and first hardwood for family_ID == "57C")
rows_to_keep <- habitat_swap_methods %>%
  filter(Family_ID == "57C") %>%
  group_by(type_log_new) %>%
  slice_head(n = 1) %>%
  ungroup()

# Filter out all other 56A rows and then add the two you want back in
habitat_swap_methods <- habitat_swap_methods %>%
  filter(Family_ID != "57C") %>%       # remove all 56A
  bind_rows(rows_to_keep)              # add back the two to keep

View(habitat_swap_methods)

#meta data missing site info for two samples
#add it back manually
habitat_swap_methods$site[grepl("49G1|49F1", habitat_swap_methods$Sample)] <- "pond_drain"

unique(habitat_swap_methods$site)

table(habitat_swap_methods$Instar)
table(habitat_swap_methods$type_original_log)

View(habitat_swap_methods)
print(nrow(habitat_swap_methods))

#remove last two rows they are empty for some reason
habitat_swap_methods <- habitat_swap_methods[1:(nrow(habitat_swap_methods) - 2), ]
print(nrow(habitat_swap_methods))

#remove decimals from log labels
habitat_swap_methods$log <- as.integer(habitat_swap_methods$log)

#count how many hardwood and softwood logs

counts <- table(
  unique(habitat_swap_methods[, c("log", "type_original_log")])$type_original_log
)
counts

#count the average, max, and min number of samples per log
library(dplyr)

habitat_swap_methods %>%
  count(log, name = "n_samples") %>%
  summarise(
    avg_samples = mean(n_samples),
    min_samples = min(n_samples),
    max_samples = max(n_samples)
  )

head(habitat_swap_methods)

unique(habitat_swap_methods$log[
  habitat_swap_methods$site %in% c("white_rock_campground", "pond_drain")
])

unique(habitat_swap_methods$log[
  habitat_swap_methods$site %in% c("white_rock_campground")
])

unique(habitat_swap_methods$log[
  habitat_swap_methods$site %in% c("pond_drain")
])
#ctd----


