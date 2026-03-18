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
"C:\Users\erinv\Box\Cryptocercus_research\Amplicon_Data_25\para\Amplicon_25_analysis_para.RData"
save.image("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//Amplicon_25_analysis_para.RData")
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//Amplicon_25_analysis_para.RData")
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//nmds.RData")

#clear the environment
rm(list=ls())

#save individual objects
saveRDS(errF_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errF_para.rds")
saveRDS(errR_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errR_para.rds")


#remove individual objects from environment
rm(readability_results)

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

#sort para sequences into fwd and reverse reads---- 
#(assuming R1 is fwd and R2 is reverse)
fnFs_para <- sort(list.files(sample_path, pattern="-Para_.*_R1_001.fastq", full.names = TRUE))
fnRs_para <- sort(list.files(sample_path, pattern="-Para_.*_R2_001.fastq", full.names = TRUE))
print(fnFs_para)
print(fnRs_para)

#keep only my samples (file contained a bunch of Gillian's samples as well)
fnFs_para <- fnFs_para[1:140]
fnRs_para <- fnRs_para[1:140]

#create short sample names for each fasta
para.namesF <- substr(basename(fnFs_para),1,nchar(basename(fnFs_para))-16)
print(para.namesF)

para.namesR <- substr(basename(fnRs_para),1,nchar(basename(fnRs_para))-16)
print(para.namesR)


#ctd----

#remove primers from para samples----
#assign fwd and reverse primers
#longer than traditional primers because they have the nex adapter on them for library prep
#nexParaV45F
FWD_para <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGCYGCGGTAATWCCAGCTCT" 
#nexParaV45R
REV_para <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTGCNCTTCCGTCAATTYCTT" 

#verify presence and orientation of primers in data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD_para.orients <- allOrients(FWD_para)
REV_para.orients <- allOrients(REV_para)
FWD_para.orients #all orientations of the forward primer (fwd, comp, reverse, reverse comp)
REV_para.orients #all orientations of the reverse primer (fwd, comp, reverse, reverse comp)

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
rbind(FWD.ForwardReads_para = sapply(FWD_para.orients, 
                                primerHits, 
                                fn = fnFs_para[[10]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_para = sapply(FWD_para.orients,
                                primerHits, 
                                fn = fnRs_para[[10]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_para = sapply(REV_para.orients, 
                                primerHits,
                                fn = fnFs_para[[10]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_para = sapply(REV_para.orients, 
                                primerHits, 
                                fn = fnRs_para[[10]])) #count all orientations of reverse primer in reverse reads
#output is table of counts
#rows are primers (FWD, REV) and read direction (forward, or reverse)
#read direction as in there is a forward and reverse fasta file
#columns are orientations of primers (forward, complement, reverse, reverse comp)
#why are there y's in some sequences what does that mean? I mean seriously wtf

#there are reverse compliments of some forward primers in the reverse reads
#there are reverse compliments of some reverse primers in the forward reads
#consistent across multiple samples 
#(make sure to test this! change the number in the double brackets in the count primers code)

#                   Forward Complement Reverse RevComp
#FWD.ForwardReads       0          0       0       0
#FWD.ReverseReads       0          0       0   21754
#REV.ForwardReads       0          0       0   19062
#REV.ReverseReads       0          0       0       0

#remove primers using cutadapt
#https://benjjneb.github.io/dada2/ITS_workflow.html
#downloaded cutadapt single file executable for windows
#https://github.com/marcelm/cutadapt/releases
#stuck it in Box

#checking that cutadapt is installed and accessible in R
system('"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe" --version')
cutadapt <- "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe"

#create output directory for cutadapted files
path.cut_paraF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnFs_noprim_para")
path.cut_paraR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnRs_noprim_para")
if(!dir.exists(path.cut_paraF)) dir.create(path.cut_paraF)
if(!dir.exists(path.cut_paraR)) dir.create(path.cut_paraR)
fnFs.cut_para <- file.path(path.cut_paraF, para.namesF)
fnRs.cut_para <- file.path(path.cut_paraR, para.namesR)

FWD_para.RC <- dada2:::rc(FWD_para)
REV_para.RC <- dada2:::rc(REV_para)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags_para <- paste("-g", FWD_para, "-a", REV_para.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags_para <- paste("-G", REV_para, "-A", FWD_para.RC) 
# Run Cutadapt
for(i in seq_along(fnFs_para)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_para[i], "-p", fnRs.cut_para[i], # output files
                             fnFs_para[i], fnRs_para[i])) # input files
}

#checking that all the file paths for input and output are correct
for(i in seq_along(fnFs_para)) {
  # Print input and output file paths for debugging
  message("Input file for forward read: ", fnFs_para[i])
  message("Input file for reverse read: ", fnRs_para[i])
  message("Output file for forward read: ", fnFs.cut_para[i])
  message("Output file for reverse read: ", fnRs.cut_para[i])
}

#count primers after cutadapt to confirm they have been removed
#repeat for several samples just in case
rbind(FWD.ForwardReads_para = sapply(FWD_para.orients, 
                                primerHits, 
                                fn = fnFs.cut_para[[10]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_para = sapply(FWD_para.orients,
                                primerHits, 
                                fn = fnRs.cut_para[[10]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_para = sapply(REV_para.orients, 
                                primerHits,
                                fn = fnFs.cut_para[[10]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_para = sapply(REV_para.orients, 
                                primerHits, 
                                fn = fnRs.cut_para[[10]])) #count all orientations of reverse primer in reverse reads
#It worked. YAY!

#ctd----

#inspect quality plots for para samples ----

#store quality plots for pre-trimmed para F samples
qualityPlots_paraF <- plotQualityProfile(fnFs.cut_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraF.png", plot = qualityPlots_paraF, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#store quality plots for pre-trimmed para R samples
qualityPlots_paraFR <- plotQualityProfile(fnRs_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraR.png", plot = qualityPlots_paraFR, width = 25, height = 15, dpi = 300, limitsize = FALSE)
#these look amazing, quality scores are consistently around 35-40
#and there is only rarely quality decline at the beginning or end of the reads
#there is some quality decline in the middle and at the end of a few reads but the scores don't drop below 25-30
#according to chat gpt the filter  step should catch any problematic reads and remove them
#there is no need to trim or truncate reads for these samples

#quality plots are not generating for the cutadapted files, but are generating for the pre-trim files
#something went wrong during teh cutadapt step
#troubleshooting that:

# Extract quality scores
qscores <- as(quality(fnFs.cut_para[1]), "matrix")
#check if files exist
file.exists(fnFs.cut_para[1])
fq <- readFastq(fnFs.cut_para[1])
# Check if the object is correctly loaded
fq
class(fq)
# Check for missing values
sum(is.na(qscores))
head_lines <- readLines(fnFs.cut_para[1], n = 8)
# Print the output
print(head_lines)
qscores <- as(fnFs.cut_para[9], "matrix")
# Check the structure of the quality score matrix
str(qscores)
# Count missing values
sum(is.na(qscores))
read_lengths <- width(sread(readFastq(fnFs.cut_para[23])))
# Check if there are zero-length reads
sum(read_lengths == 0)
#ok! So some files have zero length reads and it looks like those files are the ones 
#that are throwing errors when I try to make quality plots
#these could be introduced if reads were very short to begin with or reads contained only adapters
#to proceed, I need to remove zero length reads
#filtering out zero length reads is quicker than rerunning cutadapt

#create output directory for zero filtered files
path.nozero_paraF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnFs_nozero_para")
path.nozero_paraR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnRs_nozero_para")
if(!dir.exists(path.nozero_paraF)) dir.create(path.nozero_paraF)
if(!dir.exists(path.nozero_paraR)) dir.create(path.nozero_paraR)
fnFs.nozero_para <- file.path(path.nozero_paraF, para.namesF)
fnRs.nozero_para <- file.path(path.nozero_paraR, para.namesR)

#filter out zero length reads for fwd reads
input_filesF <- fnFs.cut_para
output_filesF <- fnFs.nozero_para

for (i in seq_along(input_filesF)) {
  # Read the FASTQ file
  fq <- readFastq(input_filesF[i])
  # Filter out zero-length reads
  fq_filtered <- fq[width(sread(fq)) > 0]
  # Write the filtered FASTQ file
  writeFastq(fq_filtered, output_filesF[i], compress = FALSE)
}

#filter out zero length reads for reverse reads
input_filesR <- fnRs.cut_para
output_filesR <- fnRs.nozero_para

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
qualityPlots_paraF_trimmed <- plotQualityProfile(fnFs.nozero_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraF_trimmed.png", plot = qualityPlots_paraF_trimmed, width = 25, height = 15, dpi = 300, limitsize = FALSE)
 

#store quality plots for pre-trimmed para R samples
qualityPlots_paraR_trimmed <- plotQualityProfile(fnRs.nozero_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraR.png", plot = qualityPlots_paraR_trimmed, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#the green line on these plots indicates the average quality score at each position across all reads
#the red lines indicate 75th and 25th percentiles for that 
#after primer trimming, the 25th percentile line is really low for all samples
#this could mean that some reads were short and dropped off in quality when primers were trimmed
#the next trimming step should fix this
#ctd----

#filter para samples----

# create directory for filtered files
path.filt_paraF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnFs_filt_para")
path.filt_paraR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25", "fnRs_filt_para")
if(!dir.exists(path.filt_paraF)) dir.create(path.filt_paraF)
if(!dir.exists(path.filt_paraR)) dir.create(path.filt_paraR)
fnFs.filt_para <- file.path(path.filt_paraF, para.namesF)
fnRs.filt_para <- file.path(path.filt_paraR, para.namesR)

#filter reads
#ok. so. removig zero length reads in the previous step was actually NOT the move
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
files <- list.files("/sfs/gpfs/tardis/home/avg6kb/fnFs_noprim_para", full.names = TRUE)
# Create new names with the .fastq extension
new_names <- paste0(files, ".fastq")
# Rename all files
file.rename(files, new_names)

#rename rev noprim files (OOD)
files <- list.files("/sfs/gpfs/tardis/home/avg6kb/fnRs_noprim_para", full.names = TRUE)
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

#export output files back to Box directory to continue working

#list files within directories
fnFs.filt_para <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25/fnFs_filt_para", full.names = TRUE)
fnRs.filt_para <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25/fnRs_filt_para", full.names = TRUE)

#store quality plots for filtered para F samples
qualityPlots_paraF_filtered <- plotQualityProfile(fnFs.filt_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraF_filtered.png", plot = qualityPlots_paraF_filtered, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#store quality plots for filtered para R samples
qualityPlots_paraR_filtered <- plotQualityProfile(fnRs.filt_para[1:140])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//qualityPlots_paraR_filtered.png", plot = qualityPlots_paraR_filtered, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#these quality plots don't look like the changed that much in the trimming step 
#not totally sure that step was necessary because there are still reads with low quality scores, 
#although the mean quality score across all reads is still high
#but I think we're good to move on down the pipeline. YAY!
#ctd----

saveRDS(ps_relative_para_tree, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//ps_relative_para_tree.rds")
saveRDS(tax_labeled, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//tax_df_para.rds")

View(tax_labeled)
#merge sequences----

#before merging seqeunces and removing chimeras,
#track reads through workflow so far to figure out how many reads have been removed/ retained at each step
#takes about 20 minutes on my laptop for 140 samples
getN <- function(x) sum(getUniques(x))   #create function to count reads in a file
track_para <- cbind(sapply(fnFs_para, getN),   #count reads in raw fwd files
               sapply(fnRs_para, getN),   #count reads in raw reverse files
               sapply(fnFs.cut_para, getN),   #count reads in fwd cutadapted files
               sapply(fnRs.cut_para, getN),   #count reads in reverse cutadapted files
               sapply(fnFs.filt_para, getN),   #count reads in fwd filtered files
               sapply(fnRs.filt_para, getN))   #count reads in reverse filtered files
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_para) <- c("rawF_para", "rawR_para", "noprimF_para", "noprimR_para", "filteredF_para", "filteredR_para")
rownames(track_para) <- para.namesF
head(track_para)
tail(track_para)

#ok, this looks good! most samples decrease in read as we move through the pipeline,
#a few have a very small number of reads but most look good
#fwd and reverse reads  have the same numbers of reads for each sample after the filtering step

#estimate error rates
#not a problem that it's only using a subset of samples it just needs x number of bases
#because error rate is a feature of the sequencing run, not of individual samples or reads. I think?
#each error estimation take about 20 minutes on my laptop (40 minutes total)

errF_para <- learnErrors(fnFs.filt_para, multithread=FALSE)
#101579803 total bases in 536026 reads from 8 samples will be used for learning the error rates.

errR_para <- learnErrors(fnRs.filt_para, multithread=FALSE)
#102566875 total bases in 588782 reads from 9 samples will be used for learning the error rates.

#visualize error rates
plotErrors(errF_para, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok

plotErrors(errR_para, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok
#the error plots aren't identical for fwd and reverse but similar 

#takes 2 seconds on laptop
derepFs_para <- derepFastq(fnFs.filt_para, verbose=TRUE)
names(derepFs) <- sample.namesF

#takes 2 seconds on laptop
derepRs_para <- derepFastq(fnRs.filt_para, verbose=TRUE)
names(derepRs) <- sample.namesR
#Encountered 75454 unique sequences from 1019071 total sequences read

#estimate variants - takes 25 minutes on my laptop for each, 50 minutes total
dadaFs_para<- dada(fnFs.filt_para, err=errF_para, multithread=FALSE)
dadaRs_para<- dada(fnRs.filt_para, err=errR_para, multithread=FALSE)

#merge fwd and reverse sequences - takes 15 minutes on my laptop
mergers_para <- mergePairs(dadaFs_para, derepFs_para, dadaRs_para, derepRs_para)
#got 8 of these printed when I ran that code
#Duplicate sequences in merged output.
#it wasn't labeled as a warning or an error though
head(mergers_para[[1]])

#construct an asv table (chimeras not yet removed)
seqtabAll_para <- makeSequenceTable(mergers_para[!grepl("Mock", names(mergers_para))])

#shows how many instances of each sequence there are I think
table(nchar(getSequences(seqtabAll_para)))

#remove chimeras - takes 5 minutes on my laptop
seqtabNoC_para <- removeBimeraDenovo(seqtabAll_para)

#calculate the percentage of merged sequences that were chimeric
sum(seqtabNoC_para)/sum(seqtabAll_para)
#0.9519227
#chimeras accounted for 5% of merged sequence reads

#get read counts for merged data
track_paradf <- as.data.frame(track_para)
merged_para <- as.data.frame(rowSums(seqtabAll_para))
print(merged_para[1])
track_paradf[7] <- merged_para[1]

library(dada2)
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//seqtabNoC_para.rds")

#get read counts for after chimera removal and add to log of read counts
nochim <- as.data.frame(rowSums(seqtabNoC_para))
print(nochim[1])
track_df[5] <- nochim[1]


#realized that I never added nochim data to the current version of my tracking data frame so doing that now
#track reads per sample
nochim <- as.data.frame(nochim)
print(head(nochim))
colnames(nochim[1]) <- "nochim"
track_paradf$nochim <- nochim[1]


#Rename the last column to "nochim"
colnames(track_paradf)[ncol(track_paradf)] <- "nochim"

# 2. Reorder the columns: move "nochim" between "merged" and "filt_200bp"
track_paradf <- track_paradf[, c(
  setdiff(names(track_paradf), "nochim") |> 
    (\(x) append(x, "nochim", after = which(x == "merged")))()
)]

#renaming and reordering the columns isn't working so I'm gonna deal with it later

#there are very few reads for 49F1-55I1 (8 samples) and 63J4. I'll probably remove these later
#might arbitrarily remove anything with fewer than 1000 reads
#ctd----

#assign taxonomy PR2----
#assign taxonomy based on reference sequences 
#do taxonomy assignments using open on demand because they crashed my laptop

#export an object with the sequences so I can upload it to open on demand
saveRDS(seqtabNoC_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//seqtabNoC_para.rds")

#not working with data frame object but does work with matrix
#just upload unzipped fasta file (.fasta)
#.tgz file doesn't work

#run this code in open on demand to get taxonomy assignments
#download pr2 database from online source
#Gillian says this is the most up to date resource for protist sequence taxonomic assignments (all of  her data is in there)
#takes about an hour and a half to run in OOD
fastaRef_pr2 <- "pr2_version_5.0.0_SSU_dada2.fasta.gz"
taxTab_pr2_para <- assignTaxonomy(seqtabNoC_para, 
                             refFasta = fastaRef_pr2,  
                             taxLevels = c("Domain", "SubDomain","Phylum", "Class", "Order", "SuperFamily", "Family", "Genus", "Species"), 
                            multithread=TRUE) #set multithread = TRUE on OOD, because it will automatically allow the code to access multiple cores
class(seqtabNoC_para)
#Warning message:
#In assignTaxonomy(seqtabNoC_para, refFasta = fastaRef_pr2, taxLevels = c("Domain",  :
#Some sequences were shorter than 50 nts and will not receive a taxonomic classification.

#check the number of sequences shorter than 250 and 50 (sequences should be 251 bp)
#sum(nchar(colnames(seqtabNoC_para)) > 250)

#to be conservative, I will probably remove these short sequences. Most of them are unclassified,
#suggesting they are not biologically relevant, but rather an artifact of sequencing errors and trimming
#any anyway there are still so many ASV's of appropriate length (5375) so it's not hurting the data
#they were probably sequences that got trimmed from some reads that had quality score dips in the middle of the reads
#ctd----

#construct phyloseq object----
#have to do phyloseq stuff on desktop because there are package dependency issues with biostrings in Rivanna
#unable to open biostrings package in rivanna

#make data frame of samples names
samdf_para <- data.frame(Sample=para.namesF)

#switch rows and columns of no chimera asv table to get it into phyloseq format
seqtabNoC_edit_para <- cbind(rownames(seqtabNoC_para), data.frame(seqtabNoC_para, row.names=NULL))
colnames(seqtabNoC_edit_para)[1] <- "Sample"
nochim_switch_para <- data.frame(t(seqtabNoC_edit_para[-1]))
colnames(nochim_switch_para) <- seqtabNoC_edit_para[, 1]

#load tax tab
taxTab_pr2_para <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//taxTab_pr2_para.rds")

#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_switch_para) == rownames(taxTab_pr2_para))

#convert asv and taxonomy table to matrices
nochim_mat_para <- as.matrix(nochim_switch_para)
#pr2_mat <- as.matrix(taxTab_pr2)
sam_mat_para <- as.matrix(samdf_para)

View(taxTab_pr2_para)
#edit column names in taxonomy table and add a column to indicate uniqe asv sequences
taxTab_test <- as.data.frame(taxTab_pr2_para)
taxTab_test$asv <- c(1:nrow(taxTab_test))

#convert to matrix
pr2_mat_para <- as.matrix(taxTab_test)
View(pr2_mat_para)


#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_mat_para) == rownames(pr2_mat_para))
#should print TRUE
setdiff(rownames(nochim_mat_para), rownames(pr2_mat_para))
#should print character (0)

#check that samples names are the same between metadata and asvtable
all(colnames(nochim_mat_para) == sam_mat_para[,1])
#should print TRUE
#it print FALSE but I think that's because there are .fastq extensions still in the nochim_mat_para
#remove fastq extensions
colnames(nochim_mat_para) <- gsub(".fastq", "", colnames(nochim_mat_para))

#now check if all sample names are the sam between metadata and asv table
all(colnames(nochim_mat_para) == sam_mat_para[,1])
#yay they are

#read in the metadata from the experiment
habitat_swap_meta_data <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//habitat_swap_meta_data.xlsx")

#merge experimental metadata with samdf
colnames(habitat_swap_meta_data)[4] <- "Sample"    #rename "roach_ID column to "Sample"
samdf_para$Sample <- sub("-.*", "", samdf_para$Sample)   #remove all but the roach id from the sample names in samdf
samdf_merge_para <- merge(samdf_para, habitat_swap_meta_data, by = "Sample")    #merge sam df with metadata


#create a rownames column
samdf_merge_para$Extra_Column <- c(1:nrow(samdf_merge_para))

#compiling phyloseq object

#remove all but roach ID from columns in nochim matrix
colnames(nochim_mat_para) <- sub("-.*", "", colnames(nochim_mat_para))

#convert the first column of meta data into rownames (essential for phyloseq to run)
rownames(samdf_merge_para) <- samdf_merge_para[,1]
samdf_merge_para[,1] <- NULL

#meta data has to be a data frame and the rownames have to be the samples for this to work. Jeez
#taxa are rows = TRUE if taxa are rownames in the asv matrix
ps_para <- phyloseq(
               otu_table(nochim_mat_para, taxa_are_rows = TRUE), 
               sample_data(samdf_merge_para), 
               tax_table(pr2_mat_para),
               phy_tree(fitJC_para$tree))
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


#remove column from track para df
#track_paradf$filt_10_reads <- NULL
#df <- df[-row_number, ]

#colnames(track_paradf)[7] <- "merged" #rename the last column of track_paradf
#rownames(track_paradf) <- sub("-.*", "", rownames(track_paradf)) #rename samples to just the roach ID


# create table to track ASV counts through filtering process (track_paradf already tracks reads per sample)
track_para_asv <- data.frame(
  step = character(),
  n_ASVs = integer(),
  stringsAsFactors = FALSE
)



#after the taxonomy assignment step, I got a warning saying some taxonomy couldn't be assigned
#because of short sequences (though the majority of sequences are not short)
#Initially:
#to be conservative, I'm going to filter out any sequences less than 251 bp (the expected length)
#assuming anything shorter is due to sequencing error and not biologically relevant
#Edit 05_07_25: 
#I spoke with Gillian and she said we expect a range or lengths for this region
#and we should filter out anything less than 200 bp but leave the rest

#remove ASV's with sequences shorter than 200 bp
asv_seqs <- taxa_names(ps_para)
print(range(nchar(asv_seqs)))
#20 465
valid_asvs <- asv_seqs[which(nchar(asv_seqs) >= 200)] # Keep only ASVs that are at least 200 bp long
print(length(valid_asvs))
ps_para_filt <- prune_taxa(valid_asvs, ps_para) # Prune phyloseq object to retain only valid ASVs

#inspect
print(nrow(tax_table(ps_para_filt)))
print(nrow(tax_table(ps_para)))

#track reads per sample
read_counts <- sample_sums(ps_para_filt)
track_paradf$filt_200bp <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "merged", n_ASVs = ntaxa(ps_para)))
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_200bp", n_ASVs = ntaxa(ps_para_filt)))


#wrangle meta data

#inspect
View(sample_data(ps_para_filt))

#change 45C3 to instar 3 (it was initially called instar 2 but later reevaluated and included in the experiment)
meta_filt <- sample_data(ps_para_filt)
meta_filt["45C3", "Instar"] <- 3
View(meta_filt)

#update get rid of extra column
meta_filt$Extra_Column <- NULL

#add weight at death
meta_filt$death_weight <- NA

#convert to data frame so I can mutate
meta_filtdf <- data.frame(sample_data(meta_filt))
View(meta_filtdf)
print(class(meta_smalldf))

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
sample_data(ps_para_filt) <- sample_data(meta_filtdf)
#inspect
View(sample_data(ps_para_filt))


#remove individuals accidentally swapped onto hardwood
samples_to_remove <- c("41C1", "41F1", "41F2")  

# Remove the selected samples from the phyloseq object
ps_para_filt2 <- prune_samples(!(sample_names(ps_para_filt) %in% samples_to_remove), ps_para_filt)
# Remove ASV's associated with only those samples
ps_para_filt2 <- prune_taxa(taxa_sums(ps_para_filt2) > 0, ps_para_filt2)


#inspect
print(nrow(tax_table(ps_para_filt)))
print(nrow(tax_table(ps_para_filt2)))

#track reads per sample
read_counts <- sample_sums(ps_para_filt2)
track_paradf$filt_badswap <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_badswap", n_ASVs = ntaxa(ps_para_filt2)))



#filter out those families of instar 3 nymphs that violate independence

# Get sample names
all_samples <- sample_names(ps_para_filt2)
print(all_samples)
print(length(all_samples))

# Identify samples starting with "56A" or "57C"
target_samples <- grep("^56A|^57C", all_samples, value = TRUE)
print(target_samples)
print(length(target_samples))
# Subset phyloseq object to include only those samples
ps_violators_para <- prune_samples(target_samples, ps_para_filt2)
# Optional: prune unused taxa
ps_violators_para <- prune_taxa(taxa_sums(ps_violators_para) > 0, ps_violators_para)

#create column for treatment so that I can evaluate which groups were in what treatment
# Extract the sample data as a data frame
sample_data_df <- as.data.frame(sample_data(ps_violators_para))
# Create new "treatment" column by pasting two existing columns
sample_data_df$treatment <- paste0(sample_data_df$type_original_log, ":", sample_data_df$type_log_new)
# Assign updated sample data back to the phyloseq object
sample_data(ps_violators_para) <- sample_data_df
View(sample_data(ps_violators_para))

#find the samples for each treatment group with the highest read counts
# Extract read counts
sample_data_df$reads <- sample_sums(ps_violators_para)
View(sample_data_df)
#add back the version with read counts
sample_data(ps_violators_para) <- sample_data_df

# Convert sample data to regular data frame
sample_df_raw <- as.data.frame(as(sample_data(ps_violators_para), "data.frame"))
sample_df_raw$sample_id <- rownames(sample_df_raw)

# Get top 1 sample (highest reads) per treatment group
top_sample_ids <- sample_df_raw %>%
  group_by(treatment) %>%
  slice_max(reads, n = 1, with_ties = FALSE) %>%
  pull(sample_id)

print(top_sample_ids)

#now remove all the samples from the phyloseq object that are from 56A or 57C 
#EXCEPT the ones with the highest read counts in each treatment group

# Get all sample names
all_samples <- sample_names(ps_para_filt2)
# Identify samples that start with 56A or 57C
to_remove <- grep("^56A|^57C", all_samples, value = TRUE)
# Subtract the ones you want to keep (top_sample_ids)
final_remove <- setdiff(to_remove, top_sample_ids)
print(final_remove) #inspect
print(length(final_remove))
# Prune the samples you want to remove
ps_para_filt3 <- prune_samples(!(sample_names(ps_para_filt2) %in% final_remove), ps_para_filt2)
#prune ASV's that were assocaited with those samples
ps_para_filt3 <- prune_taxa(taxa_sums(ps_para_filt3) > 0, ps_para_filt3)

#inspect
print(nrow(tax_table(ps_para_filt2)))
print(nrow(tax_table(ps_para_filt3)))

#track reads per sample
read_counts <- sample_sums(ps_para_filt3)
track_paradf$filt_violators <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_violators", n_ASVs = ntaxa(ps_para_filt3)))



#remove samples with low read counts
class(track_paradf)
str(track_paradf)
#make a list of samples that have 1000 or more reads
valid_samples <- rownames(track_paradf[!is.na(track_paradf$filt_violators) & track_paradf$filt_violators >= 1000, ]) 
print(length(valid_samples))
# Subset the phyloseq object to retain only valid samples
ps_para_filt4 <- prune_samples(sample_names(ps_para_filt3) %in% valid_samples, ps_para_filt3) 
#prune ASV's not present in any remaining samples
ps_para_filt4 <- prune_taxa(taxa_sums(ps_para_filt4) > 0, ps_para_filt4)

#inspect
print(nrow(tax_table(ps_para_filt3)))
print(nrow(tax_table(ps_para_filt4)))

#track reads per sample
read_counts <- sample_sums(ps_para_filt4)
track_paradf$filt_lowreads <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_lowreads", n_ASVs = ntaxa(ps_para_filt4)))




#this should filter out OTU's that have fewer than 10 reads present across all samples
#not sure this is necessary because the relative abundance filtering step would probably catch all of these
#but hey I mean you never know just in case
#To increase the reliability of microbial composition, we advise removing OTUs with 
#<10 copies in individual samples, 
#particularly in studies where only one subsample per specimen is available for analysis. (like this one)
# doi: 10.3389/fcimb.2023.1165295.
#this paper says filtering based on variance decreases reliability compared to the recommended method
filter <- phyloseq::genefilter_sample(ps_para_filt4, filterfun_sample(function(x) x >= 10))
ps_para_filt5 <- prune_taxa(filter, ps_para_filt4)

print(nrow(tax_table(ps_para_filt)))
print(nrow(tax_table(ps_para)))

#inspect
print(nrow(tax_table(ps_para_filt4)))
print(nrow(tax_table(ps_para_filt5)))

#track reads per sample
read_counts <- sample_sums(ps_para_filt5)
track_paradf$filt_10reads <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_10reads", n_ASVs = ntaxa(ps_para_filt5)))



#filtering according to Gillian's advice (remove potential contaminants and errors)
View(tax_table(ps_para_filt5))

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
asv_sample_counts <- apply(otu_table(ps_para_filt5), 1, function(x) sum(x > 0))
head(asv_sample_counts)

# Identify ASVs that appear in more than one sample
asvs_to_keep <- names(asv_sample_counts[asv_sample_counts > 1])

# Prune the phyloseq object to keep only those ASVs
ps_para_filt6 <- prune_taxa(asvs_to_keep, ps_para_filt5)

#inspect
print(nrow(tax_table(ps_para_filt5)))
print(nrow(tax_table(ps_para_filt6)))
#a huuuge number of these asv's only appeared in one sample

#track reads per sample
read_counts <- sample_sums(ps_para_filt6)
track_paradf$filt_1_sample <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_1_sample", n_ASVs = ntaxa(ps_para_filt6)))



#remove anything that didn't get assigned to at least parabasalia
#There are a few ASVs in the oxy sheet that only come out to Excavata and when I 
#blast them the closest hit is Saccinobaculus, but the % identity is low, and it’s 
#only 3 ASVs once we get rid of the ASVs that are only in 1 sample each. For 
#simplicity might as well just ignore (toss) those. Also note I spot checked a bunch 
#of the ASVs that came out as just Eukaryota (in both oxy and para), and they’re just junk, 
#like bits of bacterial genomes and stuff. Not even 18S or 16S. So I think it’s safe 
#to toss anything that doesn’t come out all the way to phylum. You have quite a few 
#random things like bacteria, archaea, fungi, etc., especially in the oxy set.

# Subset taxa where the Class is "Parabasalia"
ps_para_filt7<- subset_taxa(ps_para_filt6, Class == "Parabasalia")

#inspect
print(nrow(tax_table(ps_para_filt6)))
print(nrow(tax_table(ps_para_filt7)))
#a huuuge number of these asv's weren't parabasalia

#track reads per sample
read_counts <- sample_sums(ps_para_filt7)
track_paradf$filt_para <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_para", n_ASVs = ntaxa(ps_para_filt7)))




#remove any known contaminants
#For contamination in the para set, get rid of anything Spirotrichonymphea and Cristamonadea.
#This way you will have taken care of 99% of the known contaminants. I only saw a 
#few other ASVs in there that are for sure from termites but not Cristamonadea or Spirotrichonymphea 
#(they were Trichonympha agilis and Pseudotrichonympha). To keep your methods section simple, 
#and to keep this repeatable, you could just stick to the no Crista no Spiro rule 
#and not worry about cutting anything else. Or, if you want, I could take another look at 
#your filtered set and flag the few remaining known contaminants.

View(tax_table(ps_para_filt7))

#ok so the firs thing I'm gonna do is look at the relative abundances of known contaminants
#I'm going to save a merged table with read counts per sample and total read counts

# Extract OTU and taxonomy tables for your object
otus <- as.data.frame(otu_table(ps_para_filt7))
taxa <- as.data.frame(tax_table(ps_para_filt7))

# Add ASV sequence
otus$ASV <- rownames(otus)
taxa$ASV <- rownames(taxa)

# Merge tax and otu
ps_para_merged <- merge(taxa, otus, by = "ASV")
View(ps_para_merged)

#remove random irrelevant indexing column
ps_para_merged$asv <- NULL

#compute total reads per ASV
ps_para_merged$total_abundance <- rowSums(ps_para_merged[, sapply(ps_para_merged, is.numeric)])

#sort by ASV abundance
ps_para_sorted <- ps_para_merged[order(-ps_para_merged$total_abundance), ]
View(ps_para_sorted)




#save as two separate sheets in an excel file
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

addWorksheet(wb, "para_contam")
writeData(wb, "para_contam", ps_para_sorted, na.string = "NA")


# Save workbook
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//Taxonomy_table_para_contam.xlsx", overwrite = TRUE)
#went in an manually his the columns with the actual sequences because they were getting in the way
#yay


ps_para_filt8 <- subset_taxa(ps_para_filt7, !(Order %in% c("Cristamonadea", "Spirotrichonymphea")))

#inspect
print(nrow(tax_table(ps_para_filt7)))
print(nrow(tax_table(ps_para_filt8)))
#a large number of these were contaminants yikes


#track reads per sample
read_counts <- sample_sums(ps_para_filt8)
track_paradf$filt_contam <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_contam", n_ASVs = ntaxa(ps_para_filt8)))



#all this additional filtering brought several samples read counts lower than 1000
#so gonna remove those
ps_para_filt9 <- prune_samples(sample_sums(ps_para_filt8) >= 1000, ps_para_filt8)
ps_para_filt9 <- prune_taxa(taxa_sums(ps_para_filt9) > 0, ps_para_filt9)

#inspect
print(nrow(tax_table(ps_para_filt8)))
print(nrow(tax_table(ps_para_filt9)))


#track reads per sample
read_counts <- sample_sums(ps_para_filt9)
track_paradf$filt_lowreads2 <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_lowreads2", n_ASVs = ntaxa(ps_para_filt9)))




#remove other contaminants
#remove a few additional unidentified ASV that Gillian flagged as a contaminant
#i think she got to that by blasting the sequence and it hit a match with a contaminant


#sequence of ASV to remove (specified by Gillian)
bad_seqs <- c("ACAAGTTTGTTTCTGTTGTTACGCTTAAAACGCGGGTAGTCTGGATTTTGAGCGTCTTAATGGCGC
  GTTCACCGTAATTAAATGAGAGCGCTTAAAGCAGTTCATAACTGAATCGCTTAGTGCGGAATGAAGGCTTTACCCCAT
  TAAGGGTCAGCCCAAAGGAAGCCATCGGGGGCAGGGGTACTTGAGGGCGAGCGGCGGAATGTGTTGACCCCTCAGGGA
  CCAACGAATGCGAAGGCACCTGCCTAGAAGGTTTTCGTCGATCAACCGCGAGAGCAGCAGGATCCAACCGGATTAGAGA
  CCCGGGTAGTTGCTGCCTTAAACGATGCCGACAGGGCCACGTCCCTTAACGGGGCGGGGTCTCAGGAGAAATCATAGT
  CTCTGGGCTGCGGGGGAACTACGACCGCAAGGCTGAAACTTG", "GCGAGTTTGCTCCCATATTGTTGCAGTTAAAAC
  GCCCGTAGTCAGAAGTCCTCGGATCAGGGAACACCCGTCCACCGCGAACAAAACAAAACGCTTAGAGTATTTTCGAAGAA
  TGATTAAGCGCGGTCTGGAGTATTTACCAAGTCAGGTAAGGTCAAGGAGAGCGAGTGGTGCTGGTAGTATTCCTGAGCGA
  GCGGTGGAATGCATTGACCTTAGGGAGACTAACGAAAGCGAAAGCAGCCAGCAAGTAGGTTTCTATCGATCAAGGGCGA
  GAGTAGGGGTAGCGAACCGGATTAGAGACCCGGGTAGTCCCTACAGTAAACGATGTCGACAAGGGATTGGTTAAGAAAGA
  AGATCAGTACCTTATTGAGAAATCATAGTTCATGGACTCTGGGGGAACTACGATCGCAAGACAGAAACTTG", "GCGAG
  TTTGCTCCCATATTGTTGCAGTTAAAACGCCCGTAGTCAGAACTCAGTAGAAAAGGTACATCCGTCCACCGCGAATAA
  AACAGGACGCTTAGGGTAAACACAGGAATGATAGAGCGCGGTCTGGAGTATTTACTCACTCGAGTAAGGTCAAAGAGAG
  CGAGTGGTGCTGGTAGTATTCCTAAGCGAGCGGTGGAATGCATTGACCTTAGGGAGACTAACGAAAGCGAAGGCAGCCA
  GCAAGTAGGTTTCTATCGATCAAGGGCGAGAGTAGGGGTAGCGAACCGGATTAGAGACCCGGGTAGTCCCTACAGTAAA
  CGATGTCGACAAGGGAATGGTTTTTTTGAATCAGTACCTTTGGAGAAATCATAGTTCATGGACTCTGGGGGAACTACGA
  TCGCAAGACAGAGACTTG", "GTAAGTGTGCTCCCATATTGTTGCAGTTAAAAAGCCCGTAGTCGGATTTCAATTGA
  CTAATTATTCACTGTGAATAAATTAGGACGCTTAAAGTATGGTTGCATGAATACTGTAGCGCAGTATGAAAGATTTTGC
  TCTTGTGGCAAGATCAAAGAGAGCCATTGGGGGTATTTCTATTTCATGGCGAGCGGTGAAATGCGTTGACCCATGGGAGAG
  AAACGAAAGCGAAGGCAAATACCCAAAGGGTTTCTGTCGATCAAGGGCGAGAGTAGGGGGAGCGAACCGGATCAGAGACCC
  GGGTAGTCCCTACTGTAAACGATGTCGACAGGGGATTGTCTACTAATGTAGGCAGAACCTTAGCAAAAATGATAGTTCATGG
  ACTCTGGGGGAACTACGACTGCAAGGCTGAAACTTG", "ACAAGTTTGTTTCTGTTGTTACGCTTAAAACGCGGGTAGTC
  TGGATTTTGAGCGTCTTAATGGCGCGTTCACCGTAATTAAATGAGAGCGCTTAAAGCAGTTCATAACTGAATCGCTTAGT
  GCGGAATGAAGGCTTTACCCCATTAGGGGTCAGCCCAAAGGAAGCCATCGGGGGCAGGGGTACTTGAGGGCGAGCGG
  CGGAATGTGTTGACCCCTCAGGGACCAACGAATGCGAAGGCACCTGCCTAGAAGGTTTTCGTCGATCAACCGCGAGAG
  CAGCAGGATCCAACCGGATTAGAGACCCGGGTAGTTGCTGCCTTAAACGATGCCGACAGGGCCACGTCCCTTAACGGGG
  CGGGGTCTCAGGAGAAATCATAGTCTCTGGGCTGCGGGGGAACTACGACCGCAAGGCTGAAACTTG", "ACAAGTTTG
  TTTCTGTTGTTACGCTTAAAACGCGGGTAGTCTGGATTTTGAGCGTCTTAATGGCGCGTTCACCGTAATTAAATGAGAG
  CGCTTAAAGCAGTTCATAACTGAATCGCTTAGTGCGGAATGAAGGCTTTACCCCATTAAGGGTCAGCCCAAAGGAAGCC
  ATCGGGGGCAGGGGTACTTGAGGGCGAGCGGCGGAATGTGTTGACCCCTCAGGGACCAACGAATGCGAAGGCACCTGCCT
  AGAAGGTTTTCGTCGATCAACCGCGAGAGCAGCAGGATCCAACCGGATTAGAGACCCGGGTAGTTGCTGCCTTAAACGAT
  GCCGACAGGGCCACGTCCCTTAACGGGGCGGGGTCTCAGGAGAAATCATAGTCTCTGGGCTGCAGGGGAACTACGACCGC
  AAGGCTGAAACTTG", "GCAAGTTTGCTCCCATATTGTTGCAGTTAAAACGCCCGTAGTCTGAATTTGTGGCTTGCT
  ACACCGTCTTTTAGACGTTCACTGTGAACAAATCAGGACGCTTAGAGTATGCAATTGAATGACTTAGCGCAGTATGATG
  TTTTTACCCTAGGTAAAATCAAAGAGAGCCACCGGGGGTAGATCTATTTCATGGCGTAGCGGTGGAATGTTCTGACCCA
  TGAGAGAGAAACGAAGGCGAAGGCATCTACCTAGAGGGTTTCTGTCGATCAAGGGCGAGAGTAGAAGTATCCAACCGGAT
  TAGAGACCCGGGTAGTTTCTACCTTAAACAATGCCGACAGGGGCTTGTCTAACAAGGCAGGACCTTAGGAGAAATCAT
  AGTTCTTGGGCTCTGGGGGAACTACGACCGCAAGGCTGAAACTTG", "ACAAGTTTGTTTCTGTTGTTACGCTTAAAA
  CGCGGGTAGTCTGGATTTTGAGCGTCTTAATGGCGCGTTCGCCGTAATTAAATGAGAGCGCTTAAAGCAGTTCATAACT
  GAATCGCTTAGTGCGGAATGAAGGCTTTACCCCATTAAGGGTCAGCCCAAAGGAAGCCATCGGGGGCAGGGGTACTTGA
  GGGCGAGCGGCGGAATGTGTTGACCCCTCAGGGACCAACGAATGCGAAGGCACCTGCCTAGAAGGTTTTCGTCGATCAAC
  CGCGAGAGCAGCAGGATCCAACCGGATTAGAGACCCGGGTAGTTGCTGCCTTAAACGATGCCGACAGGGCCACGTCCCTT
  AACGGGGCGGGGTCTCAGGAGAAATCATAGTCTCTGGGCTGCGGGGGAACTACGACCGCAAGGCTGAAACTTG")



# Remove the ASV from your phyloseq object
ps_para_filt10 <- subset_taxa(ps_para_filt9, !taxa_names(ps_para_filt9) %in% bad_seqs)
rm(ps_oxy_filt10)

#inspect
print(nrow(tax_table(ps_para_filt9)))
print(nrow(tax_table(ps_para_filt10)))

#it looks like all of these sequences have also already been removed via other filtering
#so we're good!

#we missed some spirotrichonympha when filtering previously
#these are likely contaminants so we gonna take them out

View(tax_table(ps_para_filt9))

#make a fasta of seqs that contain these taxa for blasting
#Tritrichomonadida
#Tritrichomonadea

#we need to blast anything that has order = Tritrichomonadea
library(phyloseq)
library(Biostrings)


# Extract taxonomy table
tax <- tax_table(ps_para_filt9)

# Ensure we remove NA values before subsetting
asv_seqs <- rownames(tax)[!is.na(tax[, "Order"]) & tax[, "Order"] == "Tritrichomonadea"]

# Now convert to DNAStringSet
asv_fasta <- DNAStringSet(asv_seqs)
names(asv_fasta) <- paste0("ASV_", seq_along(asv_fasta))  # or use asv_seqs if you want the sequence as name

print(length(asv_fasta))

# Write to FASTA
writeXStringSet(asv_fasta, filepath = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//Tritrichomonadea_ASVs.fasta")

#Here's what Gillian said about these 5 sequences:
#Interesting, they’re not that similar to anything in genbank. I ran a quick tree
#with some of our unpublished refs, and they’re not very similar to anything we’ve 
#cloned and sequenced in our lab either. ASV 5 is pretty close to sequences from 
#Kalotermitidae, but not identical, and the rest are distant enough that I could 
#believe they're from Cryptocercus. So, I would say either keep them all, 
#or keep 1-4 and toss 5. Sound good?

#so I'll be keeping the first 4 and tossing the 5th just to be conservative

#Check the names of the sequences in the fasta file
names(asv_fasta)  # Should show ASV_1, ASV_2, ..., ASV_5
# Get the sequence to remove (the actual sequence, not the label)
seq_to_remove <- as.character(asv_fasta[[5]])
# Get all ASV sequences from ps_para_filt9
asv_seqs <- taxa_names(ps_para_filt9)  # should be DNA strings if ASVs are sequences
# Filter the phyloseq object
ps_para_filt10 <- prune_taxa(asv_seqs != seq_to_remove, ps_para_filt9)

#inspect
print(nrow(tax_table(ps_para_filt9)))
print(nrow(tax_table(ps_para_filt10)))
#1 asv removed as it should


#track reads per sample
read_counts <- sample_sums(ps_para_filt10)
track_paradf$filt_trit <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_trit", n_ASVs = ntaxa(ps_para_filt10)))

#filter pseudotrichonympha
View(tax_table(ps_para_filt10))
ps_para_filt11 <- subset_taxa(
  ps_para_filt10,
  is.na(Genus) | Genus != "Pseudotrichonympha"
)

#inspect
print(nrow(tax_table(ps_para_filt10)))
print(nrow(tax_table(ps_para_filt11)))
#2 asvs removed as it should


#track reads per sample
read_counts <- sample_sums(ps_para_filt11)
track_paradf$filt_pseudo <- read_counts[rownames(track_paradf)]

# Add entries to ASV tracking data frame
track_para_asv <- rbind(track_para_asv, data.frame(step = "filt_pseudo", n_ASVs = ntaxa(ps_para_filt11)))




#transform sample counts into relative abundances
ps_relative_para  = transform_sample_counts(ps_para_filt11, function(x) x / sum(x) )
View(otu_table(ps_relative_para))

#filter to only the most abundant taxa for small version of data for visualization
#don't try and make plots with more ASV's than this or R freaks out
ps_small_para = filter_taxa(ps_relative_para, function(x) sum(x) > .1, TRUE)
ps_lesssmall_para = filter_taxa(ps_relative_para, function(x) sum(x) > .05, TRUE)

print(nrow(tax_table(ps_relative_para)))
#now we have 112 ASV's. I'm going to use this for initial visualization but I will
#probably include more ASV's in my NMDS later

#visualize small filtered
#don't plot the bigger unfiltered data set because laptop can't deal and it doesn't show anyway
p <- plot_bar(ps_relative_para, fill = "asv")
p

library(ggplot2)
# Save as PNG
#this code isn't working because I can't get the size right so just export it manually as a png
#ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//asv_barplot_para.png", plot = p, width = 854, height = 419, dpi = 300, limitsize = FALSE)

#save tracking data as excel files 

# Load the library
library(tibble)

# Move row names to a new column called "Sample"
track_para <- rownames_to_column(track_paradf, var = "Sample")
View(track_para)
class(track_para)

str(track_para)

#there's a woerd column that's actually a list length 1 for some reason
#let's fix that 
# Replace the nested data.frame column with its vector
track_para$nochim <- track_para$nochim[[1]]

# Load the package
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "track_para")
writeData(wb, sheet = "track_para", x = track_para)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//track_para.xlsx", overwrite = TRUE)



wb <- createWorkbook()
addWorksheet(wb, "track_para_asv")
writeData(wb, sheet = "track_para_asv", x = track_para_asv)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//track_para_asv.xlsx", overwrite = TRUE)


#ctd----


#make phylogeny----

ps_relative_para <- readRDS("ps_relative_para.rds")

#save updated ps_relative_oxy phyloseq object
#this object has gillian's contamination edits and updated filtering pipelin
#delted previous versions with incomplete filteirng pipeline and contaminants
saveRDS(ps_relative_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//ps_relative_para.rds")
rm(list=ls())

load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//OOD_phylogeny_para.RData")
save.image("OOD_phylogeny_para.RData")
ps_relative_oxy <- readRDS("ps_relative_oxy.rds")
"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//OOD_phylogeny_para.RData"
ps_oxy_filt9 <- readRDS("ps_oxy_filt9.rds")
track_oxydf <- readRDS("track_oxydf.rds")
track_oxy_asv <- readRDS("track_oxy_asv.rds")

file.access(getwd(), mode = 2)

file.remove("OOD_phylogeny_para.RData")


#only include seqeunces in the phylogeny that are still here after filtering

#get ASV sequences from phyloseq object
seqs_rel <- taxa_names(ps_relative_para)

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
fit_JC_rel <- optim.pml(update(fit_rel, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#ran wayyy quicker with only 131 seqs basically instantaneous on OOD
#also basically instant on laptop yay! 131 seqs is so muh better than 5000
#15 minutes if doing stochastic rearrangement

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel <- optim.pml(update(fit_rel, model = "GTR", k = 4), 
                         optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#also basically instantaneous with 768 GB memory
#15 minutes if doing stochastic rearrangement

AIC(fit_JC_rel, fit_GTR_rel)
#> AIC(fit_JC_rel, fit_GTR_rel)
#df      AIC
#fit_JC_rel  321 20226.74
#fit_GTR_rel 330 18815.57

#GTR model has a lower AIC so we'l go with that one

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

View(tax_labeled)


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


ps_relative_para_tree <- phyloseq(otu_table(ps_relative_para), 
                                 tax_table(ps_relative_para), 
                                 tree = tree_rooted)

#make a new alignment that reflects all the filtering we did
#going to try this with large and filtered datasets
#only include seqeunces in the phylogeny that are still here after filtering








#after removing that one trit potential contaminant
#get ASV sequences from phyloseq object
seqs_rel2 <- taxa_names(ps_relative_para)

#propogates labels to tips of the tree
names(seqs_rel2) <- seqs_rel2

#make alignment
library(DECIPHER)
aligned_rel2 <- AlignSeqs(DNAStringSet(seqs_rel2), processors = 4, anchor=NA,verbose=TRUE)
names(aligned_rel2) #check that ASV seqs are names in the alignment

library(phangorn)
# Read aligned sequences into phyDat format
phangAlign_rel2 <- phyDat(as(aligned_rel2, "matrix"), type="DNA")

# Compute pairwise distances using maximum likelihood
dm_rel2 <- dist.ml(phangAlign_rel2, processors = 1, verbose = TRUE)

# Build a neighbor-joining tree
tree_NJ_rel2 <- NJ(dm_rel2)

class(tree_NJ_rel2)
# Convert to pml object
fit_rel2 <- pml(tree_NJ_rel2, data = phangAlign_rel2)


# Fit the JC69 model
#fit_JC_rel2 <- optim.pml(update(fit_rel2, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#not doing this because the GTR has been always fitting better

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel2 <- optim.pml(update(fit_rel2, model = "GTR", k = 4), 
                         optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#also basically instantaneous with 768 GB memory
#15 minutes if doing stochastic rearrangement

#AIC(fit_JC_rel, fit_GTR_rel)


#check if the trees I plan on using are rooted
is.rooted(fit_GTR_rel2$tree)

#the tree is not rooted, so we need to root it now
#because we do not have an outgroup, we will estimate the root using midpoint rooting
#this method places the root in the middle between the two longest branches

# If your tree is stored in a pml object (e.g., fit_GTR_rel$tree)
tree_unrooted2 <- fit_GTR_rel2$tree

# Root using midpoint method
tree_rooted2 <- midpoint(tree_unrooted2)

# Plot to check
plot(tree_rooted2, show.tip.label = FALSE)
#ok that looks good, it doesn't seem like removing those few seqeuencess changed things too much

#change tip labels to custom label
library(ape)

# Check current tip labels
head(tree_rooted2$tip.label)
#seqeunces

# Create a named vector mapping old → new labels
label_map <- setNames(tax_labeled$Custom_Label, tax_labeled$ASV_Sequence)

tree_label2 <- tree_rooted2

# Replace tip labels using the mapping
tree_label2$tip.label <- label_map[tree_label2$tip.label]

# Plot the relabeled tree
plot(tree_label2, cex = 0.6, no.margin = TRUE)


#add rooted tree to phyloseq object
#make sure to add the version with seqeunce IDs as tip labels
ape::is.rooted(tree_rooted2)
tree_rooted2$rooted <- TRUE

ps_relative_para_tree <- merge_phyloseq(ps_relative_para, tree_rooted2)

#ctd----


#PERMANOVA (bray curtis)----
#extract asv table

#edit treatment labels
meta_para$treatment <- paste(meta_para$type_original_log, meta_para$type_log_new, sep = ":")
meta_para$treatment <- gsub("wood", "", meta_para$treatment)
meta_para$type_original_log <- gsub("wood", "", meta_para$type_original_log)
sample_data(ps_relative_para_tree) <- sample_data(meta_para)
View(sample_data(ps_relative_para_tree))

bray_para <- phyloseq::distance(ps_relative_para_tree, method = "bray")

View(meta_para)

library(vegan)
#full models
permanova_para_bray <- adonis2(bray_para ~ Dissected + site/log + Instar + swapped*type_original_log, 
                              data = meta_para, 
                              permutations = 999, 
                              method = "bray",
                              by = "terms")
print(permanova_para_bray)


#let's check out the dispersion stuff
#check homogeneity of variances for site
disp1.1 <- betadisper(bray_para, group=meta_para$treatment)
permutest(disp1.1)
# p value > 0.05,  homogeneous variances 
#location but not variance change across levels of swapped*type_log_original


#visualization
library(vegan)
nmds_bray <- metaMDS(bray_para, k = 3, trymax = 100)
nmds_bray$stress
# 0.1978647 acceptable!

#extract scores
nmds_bray_points <- as.data.frame(nmds_bray$points)
View(nmds_bray_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_bray_points$SampleID <- rownames(nmds_bray_points)
#remove scores previously added to meta_rel
View(meta_para)
#meta_rel <- meta_rel[, -c(35:37)]
meta_para$SampleID <- rownames(meta_para)


# Now merge by SampleID
nmds_bray_meta <- merge(nmds_bray_points, meta_para, by = "SampleID")
View(nmds_bray_meta)

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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//bray_nmds_treat_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//bray_nmds_log_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//bray_nmds_swapped_para.html", selfcontained = TRUE)




#ctd----

#PERMANOVA (jaccard)----

# Compute jaccard
jaccard_para <- phyloseq::distance(ps_relative_para_tree, method = "jaccard", binary = TRUE)


#we need meta data along with distance matrix to do pERmANOVA 
#we already have meta small and meta rel from our previous analysis
#actually meta small and meta rela should be the same data frame so we can use the interchangeably

#full models
permanova_para_jaccard <- adonis2(jaccard_para ~ Dissected + site/log + Instar + swapped*type_original_log, 
                                 data = meta_para, 
                                 permutations = 999, 
                                 method = "jaccard",
                                 by = "terms")
print(permanova_para_jaccard)
#eep no significant treatment effects


#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp2.2 <- betadisper(jaccard_rel, group=meta_rel$treatment)
permutest(disp2.2)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment



#visualization
library(vegan)
nmds_jaccard <- metaMDS(jaccard_para, k = 3, trymax = 100)
nmds_jaccard$stress
#0.1907383 acceptable

#extract scores
nmds_jaccard_points <- as.data.frame(nmds_jaccard$points)
View(nmds_jaccard_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_jaccard_points$SampleID <- rownames(nmds_jaccard_points)


# Now merge by SampleID
nmds_jaccard_meta <- merge(nmds_jaccard_points, meta_para, by = "SampleID")
View(nmds_jaccard_meta)


library(scatterplot3d)
install.packages("rgl")
library(rgl)


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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//jaccard_nmds_treat_para.html", selfcontained = TRUE)



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
    text = paste("Bray NMDS Ordination (Stress =", stress_val, ")"),
    x = 0.5,
    xanchor = "center"
  ),
  legend = list(title = list(text = "type_log_original"))
)

# Save to HTML
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//jaccard_nmds_log_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//jaccard_nmds_swapped_para.html", selfcontained = TRUE)


#ctd ----


#PERMANOVA (unifrac unweighted)----

# Compute Unweighted UniFrac
#this considers only the presence or absence of taxa.
#it measures the fraction of unique branch length in the phylogenetic tree that is exclusive to each community.
#Sensitive to rare taxa because all taxa are treated equally.
#Useful when you're interested in whether communities share the same taxa, regardless of how abundant they are.


library(vegan)
permanova_rel_UU <- adonis2(unifrac_unweighted_rel ~ Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_para, 
                            permutations = 999,
                            by = "terms")
print(permanova_rel_UU)
#slightly significant interaction
#not significant when I just do treatment tho

#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp4.4 <- betadisper(unifrac_unweighted_rel, group=meta_oxy$treatment)
permutest(disp4.4)
# p value  =  0.028, non-homogeneous variances 
#variance changes across treatment could be contributing to significant PERMANOVA results



#visualization
library(vegan)
nmds_UU <- metaMDS(unifrac_unweighted_rel, k = 3, trymax = 100)
nmds_UU$stress
#stress = 0.1447448


#extract scores
nmds_UU_points <- as.data.frame(nmds_UU$points)
View(nmds_UU_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_UU_points$SampleID <- rownames(nmds_UU_points)
#remove scores previously added to meta_rel
#View(meta_oxy)
#meta_rel <- meta_rel[, -c(35:37)]
#add sample_ID column
meta_para$SampleID <- rownames(meta_para)

# Now merge by SampleID
nmds_UU_meta <- merge(nmds_UU_points, meta_para, by = "SampleID")
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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//UU_nmds_treat_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//UU_nmds_log_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//UU_nmds_swapped_para.html", selfcontained = TRUE)




#ctd ----


#PERMANOVA (unifrac weighted)----
# Compute Weighted UniFrac
#considers both presence/absence and relative abundance of taxa.
#measures the fraction of branch length that is weighted by the difference in relative abundance between communities.
#Sensitive to abundant taxa.
#Gives more weight to dominant taxa and downplays rare taxa.

unifrac_weighted_rel <- phyloseq::distance(ps_relative_para_tree, method = "unifrac", weighted = TRUE)


library(vegan)
#full models
permanova_rel_WU <- adonis2(unifrac_weighted_rel ~  Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_para, 
                            permutations = 999,
                            by = "terms")

print(permanova_rel_WU)
#significant effects of type_original log and interaction
#weird that some of the blocking variables got close rto significant in the full model




#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp3.3 <- betadisper(unifrac_weighted_rel, group=meta_para$treatment)
permutest(disp3.3)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment

disp3.2 <- betadisper(unifrac_weighted_rel, group=meta_oxy$type_original_log)
permutest(disp3.2)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment


??pairwise.adonis2



#visualization
library(vegan)
nmds_WU <- metaMDS(unifrac_weighted_rel, k = 3, trymax = 100, verbose = FALSE)
nmds_WU$stress
#0.05368489 #really good!


#extract scores
nmds_WU_points <- as.data.frame(nmds_WU$points)

#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_WU_points$SampleID <- rownames(nmds_WU_points)

# Now merge by SampleID
nmds_WU_meta <- merge(nmds_WU_points, meta_para, by = "SampleID")
View(nmds_WU_meta)


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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//WU_nmds_treat_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//WU_nmds_log_para.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//WU_nmds_swapped_para.html", selfcontained = TRUE)


#ctd----

#multiple tests corrections####
#first we need to adjust global p-values across distance metrics
# I'm testing the same hypothesis (treatment effect) under 4 lenses
#bray-curtis, jaccard, weighted unifrac, unweighted unifrac,

# Named list of your models
models <- list(
  bray = permanova_para_bray,
  jaccard = permanova_para_jaccard,
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

#apply bonferroni adjustment
results_df$p_adj_bonf <- ave(results_df$`Pr(>F)`, results_df$term, FUN = p.adjust, method = "bonferroni")
View(results_df)


#not calculating pairwise comparisons because there we no significant interactions for para data

# Function to round all numeric columns to 4 significant figures
round_sigfig <- function(df, digits = 4) {
  df[] <- lapply(df, function(x) {
    if (is.numeric(x)) signif(x, digits) else x
  })
  return(df)
}

# Apply rounding to both data frames
results_df <- round_sigfig(results_df)
View(results_df)


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

if ("p_adj_bonf" %in% names(results_df)) {
  results_df$p_adj_bonf <- add_sig_stars(results_df$p_adj_bonf)
}

View(results_df)
# Blank out repeated metric values after the first occurrence
results_df$metric <- as.character(results_df$metric)
results_df$metric[duplicated(results_df$metric)] <- ""

# Replace "NANA" strings in p-value columns with blank
cols_to_clean <- c("p_value", "p_adj_bonf")

for (col in cols_to_clean) {
  if (col %in% names(results_df)) {
    results_df[[col]][results_df[[col]] == "NANA"] <- ""
  }
}

print(results_df)


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

summary_table <- summary_table %>%
  tibble::rownames_to_column(var = "treatment")



print(summary_table)
print(class(summary_table))

print(results_df)
print(class(results_df))


# 1. Round numeric columns in results_df to 2 significant figures
results_df <- results_df %>%
  mutate(across(where(is.numeric), ~ signif(.x, 2)))

summary_table <- summary_table %>%
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

# Print to verify
print(results_df)
print(summary_table)

#it didn't really round anything. I think because some of the number columns are actually strings
#here are code snippets to fix that
# I'll need to go back and redo the emp as well, I don't think those got fully rounded

library(stringr)
#results_df
# Function to add significance stars
add_stars <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# Convert Pr(>F) and p_adj_bonf to numeric (remove any existing stars first)
results_df <- results_df %>%
  mutate(
    `Pr(>F)` = as.numeric(str_remove(`Pr(>F)`, "\\*+$")),
    p_adj_bonf = as.numeric(str_remove(p_adj_bonf, "\\*+$"))
  ) %>%
  # Round numeric columns to 2 significant figures
  mutate(across(where(is.numeric), ~ signif(.x, 2))) %>%
  # Re-add significance stars
  mutate(
    `Pr(>F)` = paste0(`Pr(>F)`, add_stars(`Pr(>F)`)),
    p_adj_bonf = paste0(p_adj_bonf, add_stars(p_adj_bonf))
  )


#summary table
# For summary_table, extract numeric values from strings like "R² = 0.006, p = 1"
summary_table <- summary_table %>%
  mutate(
    across(-treatment, ~ {
      # Extract R² and p
      matches <- str_match(.x, "R² = ([0-9.e-]+), p = ([0-9.e-]+)")
      R2 <- as.numeric(matches[,2])
      p <- as.numeric(matches[,3])
      # Round R² and p to 2 significant figures
      R2 <- signif(R2, 2)
      p <- signif(p, 2)
      # Add stars based on p-value
      stars <- add_stars(p)
      # Reconstruct string
      paste0("R² = ", R2, ", p = ", p, stars)
    })
  )

# Print to verify
print(results_df)
print(summary_table)

#reorder rows of summary table
# Desired order vector
desired_order <- c("original diet", "diet change", "diet change x original diet")

# Reorder rows by matching the treatment column to desired_order
summary_table <- summary_table[match(desired_order, summary_table$treatment), ]

print(summary_table)
View(summary_table)

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

addWorksheet(wb, "Summary")
writeData(wb, "Summary", summary_table)

# Save workbook to working directory
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//permanova_all_results_para.xlsx", overwrite = TRUE)

#ctd----

#diversity analysis----
library(phyloseq)

# Calculate shannon, and simpsons diversity for each samples
diversity_para <- estimate_richness(ps_relative_para, measures = c("Shannon", "Simpson"))
View(diversity_para)


#this code cannot calculate richness - 
#it get's tripped up by relative abundances rather than raw counts
#throws this warning if I try to include richness
#Warning message:
#The data you have provided does not have any singletons. This is highly suspicious. Results of richness
#estimates (for example) are probably unreliable, or wrong, if you have already
#trimmed low-abundance taxa from the data. We recommended that you find the un-trimmed data and retry.

#so let's calculate richness manually by counting asvs with abundance > 0 for each sample

# Extract the ASV/OTU table (samples as rows)
otu_para <- as(otu_table(ps_relative_para), "matrix")

# Make sure samples are rows and transpose if needed
if (!taxa_are_rows(ps_relative_para)) {
  otu_para<- t(otu_para)
}

View(otu_para)

# Count number of ASVs with abundance > 0 per sample
richness_para <- colSums(otu_para > 0)

# Look at the result
head(richness_para)



# Add diversity metrics to sample_data

sample_data(ps_relative_para) <- cbind(sample_data(ps_relative_para), diversity_para)
sample_data(ps_relative_para) <- cbind(sample_data(ps_relative_para), data.frame(richness_para))
View(sample_data(ps_relative_para))

#plot raw diversity metrics by treatment
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract sample metadata from phyloseq object
meta_para <- data.frame(sample_data(ps_relative_para))
View(meta_para)

#rename richness column
meta_para <- meta_para %>%
  rename(Richness = richness_para)

#add SampleID column
meta_para$SampleID <- rownames(meta_para)

#add treatment column to metadata
meta_para$treatment <- paste(meta_para$type_original_log, meta_para$type_log_new, sep = "_")

library(dplyr)
library(tidyr)
#pivot data to long format
meta_long <- meta_para %>%
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

#there don't seem to be any obvious differences in mean diversity metrics between treatments
#however, there does seem to be differences in variance in some places

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

dharma.test(glm3,n=1000)


#test simpsons

View(meta_para)

hist(meta_para$Simpson)
#data is heavily right skewed
meta_para$log_Simpsons <- log(meta_para$Simpson)
hist(meta_para$log_Simpson)
#even more right skewed

#trying a gamma distribution because chat gpt said it was good for positive right skewed data
#Gaussian (normal) distribution didn't pass the dHARMA tests

meta_para$site <- factor(meta_para$site)
meta_para <- droplevels(meta_para)

glm_simp_para <- glmmTMB(as.numeric(Simpson) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  Gamma(link = "log"), 
                        data = meta_para)

#dropping columns from rank-deficient conditional model: sitewhite_rock_campground:log
#there was probably only one log represented at this site. what this means:
#The main effect of site stays.
#The main effect of log stays.
#The other interaction terms (e.g., swapped * type_original_log) stay.
#Any non-problematic levels of the site:log interaction also stay.

#so few logs (i think only 1) in white rock campground is causing rank-defficieny issues for dharma testing
#since it;s not a main hypothesis driven effect and since it is not a significant effect, 
#I'm dropping log from these models completely

# I asked chat gpt why the permanovas were able to handle this data just fine
#Why PERMANOVA handled the small group just fine
#PERMANOVA (e.g., adonis() in vegan):
#Is a non-parametric, permutation-based method.
#It doesn't estimate coefficients for individual factor levels.
#Instead, it partitions total variance in a distance matrix across groups using permutations.
#This means:
#PERMANOVA doesn’t require estimating effects like "sitewhite_rock_campground:log" — so it doesn't care if a level has only 5 rows.
##It can still shuffle group labels and calculate sums of squares as long as each group has at least one observation.
#It is robust to imbalance, though power may be lower for small groups.
#Why GLMs (like glmmTMB) struggled
#GLMs (especially with interactions or hierarchical models):
#Estimate one parameter per level or interaction (e.g., sitewhite_rock_campground:log).
#If there's not enough data for a particular combination (e.g., only 5 observations of one site * log pair), 
#R can't estimate a unique effect → the model matrix becomes rank deficient.
#glmmTMB internally drops those unidentifiable columns, but downstream tools like DHARMa break because of this mismatch.

# so I think we're safe to drop the nested log term

install.packages("statmod")
library(statmod)
meta_para$log_Simpson <- log(meta_para$Simpson + 1e-6)  # add tiny offset if zeros exist
meta_para$one_minus_Simpson <- 1 - meta_para$Simpson #changing simpson to 1-simpson because zero inflation is easier for glmmTMB to handle apparently

library(glmmTMB)
library(car)
glm_simp_para <- glmmTMB(
  Simpson ~ Dissected + site/log + Instar + swapped * type_original_log,
  family = beta_family(link = "cloglog"),
  data = meta_para
)

summary(glm_simp_para) 
aov_simp_para <- Anova(glm_simp_para) 
aov_simp_para
#still no significant terms

table(meta_para$site)       # how many rows for each site
levels(meta_para$site) 

dharma.test(glm_simp_para,n=1000)



#this looks whack with the gamma distribution. maybe let's try a different distribution
#tweedie model didn't work with glmmTMB
#guassian had a bad qqplot
#beta_family(link = logit) slightly better but still have significant outliers and deviation

sim_res <- simulateResiduals(glm_simp_para)
plot(sim_res)
testOutliers(sim_res)

#so there are 4 outliers among the residuals
#let's visualize them

outlier_obs <- which(sim_res$scaledResiduals > 0.89)

meta_para$outlier <- FALSE
meta_para$outlier[outlier_obs] <- TRUE
table(meta_para$outlier)
View(meta_para)

# Basic histogram of your response variable (Simpson)
p <- ggplot(meta_para, aes(x = Simpson)) + 
  geom_histogram(binwidth = 0.05, fill = "lightblue", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Simpson with Outliers Highlighted",
       x = "Simpson Index",
       y = "Count")

# Overlay the outliers as red points along the x-axis at y=0 (just below histogram)
p + geom_point(data = subset(meta_para, outlier == TRUE), 
               aes(x = Simpson, y = 0), 
               color = "red", size = 3) +
  geom_text(data = subset(meta_para, outlier == TRUE),
            aes(x = Simpson, y = 0, label = rownames(subset(meta_para, outlier == TRUE))),
            color = "red", vjust = 1.5, size = 3)





#test shannon

hist(meta_para$Shannon)
#less heavily right skewed than simpsons
#using a normal dist for this data

glm_shan_para <- glmmTMB(as.numeric(Shannon) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_para)
summary(glm_shan_para) 
aov_shan_para <- Anova(glm_shan_para) 
aov_shan_para

dharma.test(glm_shan_para,n=1000)
#assumptions test looks good!
# so no significant effects 



#test richness

hist(meta_para$Richness)
#this one is a scooch left skewed
#which I would expect for richness data

glm_rich_para <- glmmTMB(as.numeric(Richness) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_para)
summary(glm_rich_para) 
aov_rich_para <- Anova(glm_rich_para) 
aov_rich_para

dharma.test(glm_rich_para,n=1000)
#assumptions test looks good!
# so no significant effects 

pairwise_rich_para_instar <- emmeans(glm_rich_para, pairwise ~ Instar, adjust = "tukey")
pairwise_rich_para_instar$contrasts
#3 and A are very different




#export anova tables as an excel file
#no contrasts because there were no significant swap, original log type, or treatment effects


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
    table = aov_simp_para
  ),
  Shannon = list(
    formula = "Shannon ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_shan_para
  ),
  Richness = list(
    formula = "Richness ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_rich_para
  )
)
print(aov_simp_para)
print(aov_shan_para)
print(aov_rich_para)

library(openxlsx)
library(dplyr)

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
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//diversity_stats_tables_para.xlsx", overwrite = TRUE)




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
  get_emm_df(glm_simp_para, "treatment", "Simpson"),
  get_emm_df(glm_simp_para, "swapped", "Simpson"),
  get_emm_df(glm_simp_para, "type_original_log", "Simpson"),
  
  get_emm_df(glm_shan_para, "treatment", "Shannon"),
  get_emm_df(glm_shan_para, "swapped", "Shannon"),
  get_emm_df(glm_shan_para, "type_original_log", "Shannon"),
  
  get_emm_df(glm_rich_para, "treatment", "Richness"),
  get_emm_df(glm_rich_para, "swapped", "Richness"),
  get_emm_df(glm_rich_para, "type_original_log", "Richness")
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

# Function to plot one predictor
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
  "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//diversity_plot.png",
  final_plot,
  width = 7,
  height = 10,
  dpi = 600
)
#ctd----


#assign taxonomy LabInfo----
#assign taxonomy based on reference sequences 
#do taxonomy assignments using open on demand because they crashed my laptop
#this is an updated reference from gillian that is supposed to be better/ more complete than pr2
#I've done a lot of filtering to my phyloseq object since the original taxonomy assignment
#so instead of assigning taxonomy to my original setab_NoC objetc, I'm going to extract 
#sequences from the existing phyloseq object, do taxonomy assignment, and then add it back to the phyloseq obejct
#this should hopefully mak ethe runtime quicker and avoid the need to refilter
rm(list=ls())

#export current phyloseq object so I can upload it to OOD
#current phyloseq object doesn't have an associated tree, I need to remake that
#because I changed the filtering pararamters so now the taxa present are probably different somewhat
saveRDS(ps_relative_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//ps_relative_para.rds")

#import phyloseq object into R
ps_relative_para <- readRDS("ps_relative_para.rds")

library(Biostrings)

#this reference contains two files:
#a fasta file with sequences (.fasta)
#a text file with taxonomy (.tsv)
#we need to merge them

# Specify the path to the FASTA file
fasta_file <- "LabInfo_sequences.fasta"
file.exists("LabInfo_sequences.fasta")

# Read the FASTA file
sequences <- readDNAStringSet(fasta_file)

# Check the first few sequences
head(sequences)
##oooo pretty colors cool
#looks good and normal

# Read the tax file
#even tho ths is a tsv file it should read in the same as a txt
taxonomy <- read.table("LabInfo_taxonomy.tsv",
                       sep = "\t",
                       header = FALSE,
                       stringsAsFactors = FALSE,
                       quote = "",
                       fill = TRUE)
#rename columns
colnames(taxonomy) <- c("id", "taxonomy")

#inspect
head(taxonomy)
head(taxonomy$id)
head(taxonomy$taxonomy)

#inspect sequence IDS
names(sequences)

#count IDS for txt and fasta file and make sure they match
nrow(taxonomy)
length(sequences)
#they are the same length, now make sure they all match

# Match IDs between FASTA and taxonomy table
matched_tax <- taxonomy[match(names(sequences), taxonomy$id), ]

# Confirm matches are good
if (any(is.na(matched_tax$taxonomy))) {
  warning("Some sequences do not have matching taxonomy.")
}
#no warning so they all match. yay!


# Construct new headers
new_headers <- paste(matched_tax$id, matched_tax$taxonomy)

# Set new headers
names(sequences) <- new_headers

head(sequences)

#get seqeunces from phyloseq object
seqs <- taxa_names(ps_relative_para)

#write reference to fasta file
writeXStringSet(sequences, filepath = "LabInfo_full_reference.fasta")

#assign taxonomy
library(dada2)
taxTab <- assignTaxonomy(seqs, refFasta = "LabInfo_full_reference.fasta", multithread = TRUE)
#takes about 15 minutes

#save taxonomy as the tax tab for the ps object
ps_relative_para_LabInfo <- ps_relative_para
View(tax_table(ps_relative_para))
tax_table(ps_relative_para_LabInfo) <- tax_table(as.matrix(taxTab))
ps_relative_para_LabInfo <- ps_relative_para
View(tax_table(ps_relative_para_LabInfo))
saveRDS(ps_relative_para_LabInfo, file = "ps_relative_para_LabInfo.rds")

#ok I think the existing version had taxonomy assignment from the original PR2 
#without lab data so we should be good

#download phyloseq object then
#import phyloseq object into R environment
#now whene I save my environment, the new object will be there
ps_relative_para_LabInfo <- readRDS("ps_relative_para_LabInfo.rds")

#ctd----

