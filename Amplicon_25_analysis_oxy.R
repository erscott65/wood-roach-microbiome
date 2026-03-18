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
save.image("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//Amplicon_25_analysis_oxy.RData")
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//Amplicon_25_analysis_oxy.RData")
gc()


load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//OOD_phylogeny_oxy.RData")

#clear the environment
rm(list=ls())
rm(aligned_rel, fit_GTR_rel, fit_JC_rel, fit_rel, phangAlign_rel, ps_relative_oxy, tree_NJ_rel, tree_unrooted, dm_rel, seqs_rel)

#save individual objects
saveRDS(errF_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errF_para.rds")
saveRDS(errR_para, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//errR_para.rds")

rm(qualityPlots_oxyFR, qualityPlots_oxyR)

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

#sort oxy sequences into fwd and reverse reads----
#(assuming R1 is fwd and R2 is reverse)
fnFs_oxy <- sort(list.files(sample_path, pattern="-Oxy_.*_R1_001.fastq", full.names = TRUE))
fnRs_oxy <- sort(list.files(sample_path, pattern="-Oxy_.*_R2_001.fastq", full.names = TRUE))
print(fnFs_oxy)
print(fnRs_oxy)

#keep only my samples (file contained a bunch of Gillian's samples as well)
fnFs_oxy <- fnFs_oxy[1:135]
fnRs_oxy <- fnRs_oxy[1:135]

#create short sample names for each fasta
oxy.namesF <- substr(basename(fnFs_oxy),1,nchar(basename(fnFs_oxy))-16)
print(oxy.namesF)

warnings()
oxy.namesR <- substr(basename(fnRs_oxy),1,nchar(basename(fnRs_oxy))-16)
print(oxy.namesR)

#print without quotes
cat(oxy.namesF, sep="\n")
#ctd----

#remove primers from oxy samples----
#assign fwd and reverse primers
#longer than traditional primers because they have the nex adapter on them for library prep
#nexOxyV4F
FWD_oxy <- "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGAAGTCTGGTGCCAGCAG" 
#nexOxyV4R
REV_oxy <- "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTTTATTATTCCATGCTAATGTGTTC" 

#verify presence and orientation of primers in data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = Biostrings::reverse(dna),
               RevComp = Biostrings::reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD_oxy.orients <- allOrients(FWD_oxy)
REV_oxy.orients <- allOrients(REV_oxy)
FWD_oxy.orients #all orientations of the forward primer (fwd, comp, reverse, reverse comp)
REV_oxy.orients #all orientations of the reverse primer (fwd, comp, reverse, reverse comp)

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
rbind(FWD.ForwardReads_oxy= sapply(FWD_oxy.orients, 
                                   primerHits, 
                                   fn = fnFs_oxy[[1]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_oxy = sapply(FWD_oxy.orients,
                                    primerHits, 
                                    fn = fnRs_oxy[[1]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_oxy = sapply(REV_oxy.orients, 
                                    primerHits,
                                    fn = fnFs_oxy[[1]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_oxy = sapply(REV_oxy.orients, 
                                    primerHits, 
                                    fn = fnRs_oxy[[1]])) #count all orientations of reverse primer in reverse reads
#output is table of counts
#rows are primers (FWD, REV) and read direction (forward, or reverse)
#read direction as in there is a forward and reverse fasta file
#columns are orientations of primers (forward, complement, reverse, reverse comp)
#why are there y's in some sequences what does that mean? I mean seriously wtf

#there are reverse compliments of some forward primers in the reverse reads
#there are reverse compliments of some reverse primers in the forward reads
#consistent across multiple samples 
#(make sure to test this! change the number in the double brackets in the count primers code)

#                         Forward Complement Reverse RevComp
#FWD.ForwardReads_oxy       0          0       0       0
#FWD.ReverseReads_oxy       0          0       0   17114
#REV.ForwardReads_oxy       0          0       0   14101
#REV.ReverseReads_oxy       0          0       0       0

#remove primers using cutadapt
#https://benjjneb.github.io/dada2/ITS_workflow.html
#downloaded cutadapt single file executable for windows
#https://github.com/marcelm/cutadapt/releases
#stuck it in Box

#checking that cutadapt is installed and accessible in R
system('"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe" --version')
cutadapt <- "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//cutadapt.exe.exe"

#create output directory for cutadapted files
path.cut_oxyF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnFs_noprim_oxy")
path.cut_oxyR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnRs_noprim_oxy")
if(!dir.exists(path.cut_oxyF)) dir.create(path.cut_oxyF)
if(!dir.exists(path.cut_oxyR)) dir.create(path.cut_oxyR)
fnFs.cut_oxy <- file.path(path.cut_oxyF, oxy.namesF)
fnRs.cut_oxy <- file.path(path.cut_oxyR, oxy.namesR)

FWD_oxy.RC <- dada2:::rc(FWD_oxy)
REV_oxy.RC <- dada2:::rc(REV_oxy)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags_oxy <- paste("-g", FWD_oxy, "-a", REV_oxy.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags_oxy <- paste("-G", REV_oxy, "-A", FWD_oxy.RC) 
# Run Cutadapt
for(i in seq_along(fnFs_oxy)) {
  system2(cutadapt, args = c(R1.flags_oxy, R2.flags_oxy, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut_oxy[i], "-p", fnRs.cut_oxy[i], # output files
                             fnFs_oxy[i], fnRs_oxy[i])) # input files
}

#count primers after cutadapt to confirm they have been removed
#repeat for several samples just in case
rbind(FWD.ForwardReads_oxy = sapply(FWD_oxy.orients, 
                                    primerHits, 
                                    fn = fnFs.cut_oxy[[5]]), #count all orientations of fwd primer in fwd reads
      FWD.ReverseReads_oxy = sapply(FWD_oxy.orients,
                                    primerHits, 
                                    fn = fnRs.cut_oxy[[5]]), #count all orientations of fwd primer in reverse reads
      REV.ForwardReads_oxy = sapply(REV_oxy.orients, 
                                    primerHits,
                                    fn = fnFs.cut_oxy[[5]]), #count all orientations of reverse primer in fwd reads
      REV.ReverseReads_oxy = sapply(REV_oxy.orients, 
                                    primerHits, 
                                    fn = fnRs.cut_oxy[[5]])) #count all orientations of reverse primer in reverse reads

#It worked. YAY!


#rename with fwd and reverse identifiers for ncbi sra upload
#check current filenames
fnF_path <- "C:/Users/erinv/Box/Cryptocercus_research/Amplicon_Data_25/oxy/fnFs_noprim_oxy"
fnR_path <- "C:/Users/erinv/Box/Cryptocercus_research/Amplicon_Data_25/oxy/fnRs_noprim_oxy"

fnFs_noprim <- list.files(fnF_path)
fnRs_noprim <- list.files(fnR_path)

fnFs_noprim
fnRs_noprim

length(fnFs_noprim)
length(fnRs_noprim)

#rename files with fwd adn reverse identifiers as well as .fastq.gz extensions
# forward files
oldF <- file.path(fnF_path, fnFs_noprim)
newF <- file.path(fnF_path, paste0(fnFs_noprim, "_R1.fastq.gz"))

file.rename(oldF, newF)

# reverse files
oldR <- file.path(fnR_path, fnRs_noprim)
newR <- file.path(fnR_path, paste0(fnRs_noprim, "_R2.fastq.gz"))

file.rename(oldR, newR)

#verfiy new names 
list.files(fnF_path)[1:10]
list.files(fnR_path)[1:10]
#ctd----

#inspect quality plots for oxy samples ----

#store quality plots for pre-trimmed oxy F samples
qualityPlots_oxyF <- plotQualityProfile(fnFs_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyF.png", plot = qualityPlots_oxyF, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#store quality plots for pre-trimmed oxy R samples
qualityPlots_oxyR <- plotQualityProfile(fnRs_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyR.png", plot = qualityPlots_oxyR, width = 25, height = 15, dpi = 300, limitsize = FALSE)
#these look good, quality scores are consistently around 35-40 
#and there is some quality decline in the middle or end of reads for a few samples but not most
#there is no need to trim or truncate reads for these samples


#visualize quality scores for cutadapted samples
#zero length reads were introduced by cutadapt (trouble shooting that I used to figure this out is in the para R file)
#these will be removed in the filtering step

#remove zero length reads for visualization purposes
#create output directory for zero filtered files
path.nozero_oxyF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnFs_nozero_oxy")
path.nozero_oxyR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnRs_nozero_oxy")
if(!dir.exists(path.nozero_oxyF)) dir.create(path.nozero_oxyF)
if(!dir.exists(path.nozero_oxyR)) dir.create(path.nozero_oxyR)
fnFs.nozero_oxy <- file.path(path.nozero_oxyF, oxy.namesF)
fnRs.nozero_oxy <- file.path(path.nozero_oxyR, oxy.namesR)

#filter out zero length reads for fwd reads
input_filesF <- fnFs.cut_oxy
output_filesF <- fnFs.nozero_oxy

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

#filter oxy samples----

# create directory for filtered files
path.filt_oxyF <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnFs_filt_oxy")
path.filt_oxyR <- file.path("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy", "fnRs_filt_oxy")
if(!dir.exists(path.filt_oxyF)) dir.create(path.filt_oxyF)
if(!dir.exists(path.filt_oxyR)) dir.create(path.filt_oxyR)
fnFs.filt_oxy <- file.path(path.filt_oxyF, oxy.namesF)
fnRs.filt_oxy <- file.path(path.filt_oxyR, oxy.namesR)

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
fnFs.filt_oxy <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//fnFs_filt_oxy", full.names = TRUE)
fnRs.filt_oxy <- list.files("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//fnRs_filt_oxy", full.names = TRUE)

#store quality plots for filtered para F samples
qualityPlots_oxyF_filtered <- plotQualityProfile(fnFs.filt_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyF_filtered.png", plot = qualityPlots_oxyF_filtered, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#store quality plots for filtered para R samples
qualityPlots_oxyR_filtered <- plotQualityProfile(fnRs.filt_oxy[1:135])  # Store plot in variable
ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//qualityPlots_oxyR_filtered.png", plot = qualityPlots_oxyR_filtered, width = 25, height = 15, dpi = 300, limitsize = FALSE)

#the trimming step seems to have removed the quality drop off at the end of the reverse reads
#but otherwise not changed too much
#not totally sure that step was necessary because there are still reads with low quality scores, 
#although the mean quality score across all reads is still high
#but I think we're good to move on down the pipeline. YAY!
#ctd----

#merge sequences----

#before merging sequences and removing chimeras,
#track reads through workflow so far to figure out how many reads have been removed/ retained at each step
#takes about 20 minutes on my laptop for 140 samples
library(dada2)
getN <- function(x) sum(getUniques(x))   #create function to count reads in a file
track_oxy <- cbind(sapply(fnFs_oxy, getN),   #count reads in raw fwd files
                    sapply(fnRs_oxy, getN),   #count reads in raw reverse files
                    sapply(fnFs.cut_oxy, getN),   #count reads in fwd cutadapted files
                    sapply(fnRs.cut_oxy, getN),   #count reads in reverse cutadapted files
                    sapply(fnFs.filt_oxy, getN),   #count reads in fwd filtered files
                    sapply(fnRs.filt_oxy, getN))   #count reads in reverse filtered files
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track_oxy) <- c("rawF_oxy", "rawR_oxy", "noprimF_oxy", "noprimR_oxy", "filteredF_oxy", "filteredR_oxy")
rownames(track_oxy) <- oxy.namesF
head(track_oxy)
tail(track_oxy)
#ok, this looks good! most samples decrease in read as we move through the pipeline,
#a few have a very small number of reads <1000 but most look good
#fwd and reverse reads  have the same numbers of reads for each sample after the filtering step

#estimate error rates
#not a problem that it's only using a subset of samples it just needs x number of bases
#because error rate is a feature of the sequencing run, not of individual samples or reads. I think?
#each error estimation take about 20 minutes on my laptop (40 minutes total)

errF_oxy <- learnErrors(fnFs.filt_oxy, multithread=FALSE)
#116747470 total bases in 531748 reads from 6 samples will be used for learning the error rates.

errR_oxy <- learnErrors(fnRs.filt_oxy, multithread=FALSE)
#103638569 total bases in 531748 reads from 6 samples will be used for learning the error rates.

#visualize error rates
plotErrors(errF_oxy, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok

plotErrors(errR_oxy, nominalQ=TRUE)
#not sure what I'm supposed to see but the black lines seem to fit the data ok
#the error plots aren't identical for fwd and reverse but similar 


#takes 2 seconds on laptop
derepFs_oxy <- derepFastq(fnFs.filt_oxy, verbose=TRUE)
names(derepFs_oxy) <- oxy.namesF

#takes 2 seconds on laptop
derepRs_oxy <- derepFastq(fnRs.filt_oxy, verbose=TRUE)
names(derepRs_oxy) <- oxy.namesR
#Encountered 75454 unique sequences from 1019071 total sequences read

#estimate variants - takes 25 minutes on my laptop for each, 50 minutes total
dadaFs_oxy<- dada(fnFs.filt_oxy, err=errF_oxy, multithread=FALSE)
dadaRs_oxy<- dada(fnRs.filt_oxy, err=errR_oxy, multithread=FALSE)

#merge fwd and reverse sequences - takes 15 minutes on my laptop
mergers_oxy <- mergePairs(dadaFs_oxy, derepFs_oxy, dadaRs_oxy, derepRs_oxy)
#got 4 of these printed when I ran that code
#Duplicate sequences in merged output.
#it wasn't labeled as a warning or an error though
head(mergers_oxy[[1]])

#construct an asv table (chimeras not yet removed)
seqtabAll_oxy <- makeSequenceTable(mergers_oxy[!grepl("Mock", names(mergers_oxy))])

#shows how many instances of each sequence there are I think
table(nchar(getSequences(seqtabAll_oxy)))

#remove chimeras - takes 5 minutes on my laptop
seqtabNoC_oxy <- removeBimeraDenovo(seqtabAll_oxy)

#calculate the percentage of merged sequences that were chimeric
sum(seqtabNoC_oxy)/sum(seqtabAll_oxy)
#0.9785458
#chimeras accounted for 2% of merged sequence reads

#get read counts for merged data
track_oxydf <- as.data.frame(track_oxy)
merged_oxy <- as.data.frame(rowSums(seqtabAll_oxy))
print(merged_oxy[1])
track_oxydf[7] <- merged_oxy[1]

#get read counts for after chimera removal and add to log of read counts
nochim_oxy <- as.data.frame(rowSums(seqtabNoC_oxy))
print(nochim_oxy[1])
track_oxydf[8] <- nochim_oxy[1]
#there are very few reads for 49F1-55I1 (8 samples) and 63J4. I'll probably remove these later
#might arbitrarily remove anything with fewer than 1000 reads
#ctd----

#assign taxonomy PR2 ----
#assign taxonomy based on reference sequences 
#do taxonomy assignments using open on demand because they crashed my laptop

#export an object with the sequences so I can upload it to open on demand
saveRDS(seqtabNoC_oxy, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//seqtabNoC_oxy.rds")
saveRDS(seqtabNoC_oxy, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//seqtabNoC_oxy.rds")

#not working with data frame object but does work with matrix
#just upload unzipped fasta file (.fasta)
#.tgz file doesn't work

#run this code in open on demand to get taxonomy assignments
#download pr2 database from online source
#Gillian says this is the most up to date resource for protist sequence taxonomic assignments (all of  her data is in there)
#takes about an hour and a half to run in OOD
fastaRef_pr2 <- "pr2_version_5.0.0_SSU_dada2.fasta.gz"
taxTab_pr2_oxy <- assignTaxonomy(seqtabNoC_oxy, 
                                  refFasta = fastaRef_pr2,  
                                  taxLevels = c("Domain", "SubDomain","Phylum", "Class", "Order", "SuperFamily", "Family", "Genus", "Species"), 
                                  multithread=TRUE) #set multithread = TRUE on OOD, because it will automatically allow the code to access multiple cores
??assignTaxonomy
#Warning message:
#In assignTaxonomy(seqtabNoC_para, refFasta = fastaRef_pr2, taxLevels = c("Domain",  :
#Some sequences were shorter than 50 nts and will not receive a taxonomic classification.

#check the number of sequences shorter than 250 and 50 (sequences should be 251 bp)
#sum(nchar(colnames(seqtabNoC_oxy)) > 250)

#to be conservative, I will probably remove these short sequences. Most of them are unclassified,
#suggesting they are not biologically relevant, but rather an artifact of sequencing errors and trimming
#any anyway there are still so many ASV's of appropriate length (5375) so it's not hurting the data
#they were probably sequences that got trimmed from some reads that had quality score dips in the middle of the reads
#ctd----

#construct phyloseq object----
#have to do phyloseq stuff on desktop because there are package dependency issues with biostrings in Rivanna
#unable to open biostrings package in rivanna

#read in no chimera ASV file -- not sure why this was a separate file oh well
seqtabNoC_oxy <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//seqtabNoC_oxy.rds")

#make data frame of samples names
samdf_oxy <- data.frame(Sample=oxy.namesF)
head(samdf_oxy)

#switch rows and columns of no chimera asv table to get it into phyloseq format
nochim_switch_oxy <- data.frame(t(seqtabNoC_oxy), check.names = FALSE)

#load tax tab (if not already)
taxTab_pr2_oxy <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//taxTab_pr2_oxy.rds")

#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_switch_oxy) == rownames(taxTab_pr2_oxy))

#convert to matrices
nochim_mat_oxy <- as.matrix(nochim_switch_oxy) #convert asv table to matrix
#pr2_mat <- as.matrix(taxTab_pr2)
sam_mat_oxy <- as.matrix(samdf_oxy) #convert sample data to matrix

View(taxTab_pr2_oxy)
#edit column names in taxonomy table and add a column to indicate uniqe asv sequences
taxTab_test <- as.data.frame(taxTab_pr2_oxy)
taxTab_test$asv <- c(1:nrow(taxTab_test))

#convert to matrix
pr2_mat_oxy <- as.matrix(taxTab_test)
View(pr2_mat_oxy)


#check that sequences are the same between the asv table and the taxonomy table
all(rownames(nochim_mat_oxy) == rownames(pr2_mat_oxy))
#should print TRUE
setdiff(rownames(nochim_mat_oxy), rownames(pr2_mat_oxy))
#should print character (0)

#check that samples names are the same between metadata and asvtable
all(colnames(nochim_mat_oxy) == sam_mat_oxy[,1])
#should print TRUE
#it print FALSE but I think that's because there are .fastq extensions still in the nochim_mat_para
#remove fastq extensions
colnames(nochim_mat_oxy) <- gsub(".fastq", "", colnames(nochim_mat_oxy))

#now check if all sample names are the sam between metadata and asv table
all(colnames(nochim_mat_oxy) == sam_mat_oxy[,1])
#prints TRUE
#yay!

library(readxl)
#read in the metadata from the experiment
habitat_swap_meta_data <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//habitat_swap_meta_data.xlsx")

#merge experimental metadata with samdf
colnames(habitat_swap_meta_data)[4] <- "Sample"    #rename "roach_ID column to "Sample"
samdf_oxy$Sample <- sub("-.*", "", samdf_oxy$Sample)   #remove all but the roach id from the sample names in samdf
samdf_merge_oxy <- merge(samdf_oxy, habitat_swap_meta_data, by = "Sample")    #merge sam df with metadata


#create a rownames column
samdf_merge_oxy$Extra_Column <- c(1:nrow(samdf_merge_oxy))

#compiling phyloseq object

#remove all but roach ID from columns in nochim matrix
colnames(nochim_mat_oxy) <- sub("-.*", "", colnames(nochim_mat_oxy))

#convert the first column of meta data into rownames (essential for phyloseq to run)
rownames(samdf_merge_oxy) <- samdf_merge_oxy[,1]
samdf_merge_oxy[,1] <- NULL

#meta data has to be a data frame and the rownames have to be the samples for this to work. Jeez
#taxa are rows = TRUE if taxa are rownames in the asv matrix
library(phyloseq)
ps_oxy <- phyloseq(
  otu_table(nochim_mat_oxy, taxa_are_rows = TRUE), 
  sample_data(samdf_merge_oxy), 
  tax_table(pr2_mat_oxy))
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


#remove column from track para df
#track_paradf$filt_10_reads <- NULL
#df <- df[-row_number, ]

#colnames(track_paradf)[7] <- "merged" #rename the last column of track_paradf
#rownames(track_paradf) <- sub("-.*", "", rownames(track_paradf)) #rename samples to just the roach ID


# create table to track ASV counts through filtering process (track_paradf already tracks reads per sample)
track_oxy_asv <- data.frame(
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
asv_seqs <- taxa_names(ps_oxy)
print(range(nchar(asv_seqs)))
#20 486
valid_asvs <- asv_seqs[which(nchar(asv_seqs) >= 200)] # Keep only ASVs that are at least 200 bp long
print(length(valid_asvs))
ps_oxy_filt <- prune_taxa(valid_asvs, ps_oxy) # Prune phyloseq object to retain only valid ASVs

#inspect
print(nrow(tax_table(ps_oxy_filt)))
print(nrow(tax_table(ps_oxy)))

#track reads per sample 
read_counts <- sample_sums(ps_oxy_filt)
track_oxydf$filt_200bp <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "merged", n_ASVs = ntaxa(ps_oxy)))
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_200bp", n_ASVs = ntaxa(ps_oxy_filt)))


#wrangle meta data

#inspect
View(sample_data(ps_oxy_filt))

#change 45C3 to instar 3 (it was initially called instar 2 but later reevaluated and included in the experiment)
meta_filt <- sample_data(ps_oxy_filt)
meta_filt["45C3", "Instar"] <- 3
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
sample_data(ps_oxy_filt) <- sample_data(meta_filtdf)
#inspect
View(sample_data(ps_oxy_filt))


#remove individuals accidentally swapped onto hardwood
samples_to_remove <- c("41C1", "41F1", "41F2")  

# Remove the selected samples from the phyloseq object
ps_oxy_filt2 <- prune_samples(!(sample_names(ps_oxy_filt) %in% samples_to_remove), ps_oxy_filt)
# Remove ASV's associated with only those samples
ps_oxy_filt2 <- prune_taxa(taxa_sums(ps_oxy_filt2) > 0, ps_oxy_filt2)


#inspect
print(nrow(tax_table(ps_oxy_filt)))
print(nrow(tax_table(ps_oxy_filt2)))

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt2)
track_oxydf$filt_badswap <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_badswap", n_ASVs = ntaxa(ps_oxy_filt2)))



#filter out those families of instar 3 nymphs that violate independence

# Get sample names
all_samples <- sample_names(ps_oxy_filt2)
print(all_samples)
print(length(all_samples))

# Identify samples starting with "56A" or "57C"
target_samples <- grep("^56A|^57C", all_samples, value = TRUE)
print(target_samples)
print(length(target_samples))
# Subset phyloseq object to include only those samples
ps_violators_oxy <- prune_samples(target_samples, ps_oxy_filt2)
# Optional: prune unused taxa
ps_violators_oxy <- prune_taxa(taxa_sums(ps_violators_oxy) > 0, ps_violators_oxy)

#create column for treatment so that I can evaluate which groups were in what treatment
# Extract the sample data as a data frame
sample_data_df <- as.data.frame(sample_data(ps_violators_oxy))
# Create new "treatment" column by pasting two existing columns
sample_data_df$treatment <- paste0(sample_data_df$type_original_log, ":", sample_data_df$type_log_new)
# Assign updated sample data back to the phyloseq object
sample_data(ps_violators_oxy) <- sample_data_df
View(sample_data(ps_violators_oxy))

#find the samples for each treatment group with the highest read counts
# Extract read counts
sample_data_df$reads <- sample_sums(ps_violators_oxy)
View(sample_data_df)
#add back the version with read counts
sample_data(ps_violators_oxy) <- sample_data_df

# Convert sample data to regular data frame
sample_df_raw <- as.data.frame(as(sample_data(ps_violators_oxy), "data.frame"))
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
all_samples <- sample_names(ps_oxy_filt2)
# Identify samples that start with 56A or 57C
to_remove <- grep("^56A|^57C", all_samples, value = TRUE)
# Subtract the ones you want to keep (top_sample_ids)
final_remove <- setdiff(to_remove, top_sample_ids)
print(final_remove) #inspect
print(length(final_remove))
# Prune the samples you want to remove
ps_oxy_filt3 <- prune_samples(!(sample_names(ps_oxy_filt2) %in% final_remove), ps_oxy_filt2)
#prune ASV's that were assocaited with those samples
ps_oxy_filt3 <- prune_taxa(taxa_sums(ps_oxy_filt3) > 0, ps_oxy_filt3)

#inspect
print(nrow(tax_table(ps_oxy_filt2)))
print(nrow(tax_table(ps_oxy_filt3)))

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt3)
track_oxydf$filt_violators <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_violators", n_ASVs = ntaxa(ps_oxy_filt3)))



#remove samples with low read counts
class(track_oxydf)
str(track_oxydf)
#make a list of samples that have 1000 or more reads
valid_samples <- rownames(track_oxydf[!is.na(track_oxydf$filt_violators) & track_oxydf$filt_violators >= 1000, ]) 
print(length(valid_samples))
# Subset the phyloseq object to retain only valid samples
ps_oxy_filt4 <- prune_samples(sample_names(ps_oxy_filt3) %in% valid_samples, ps_oxy_filt3) 
#prune ASV's not present in any remaining samples
ps_oxy_filt4 <- prune_taxa(taxa_sums(ps_oxy_filt4) > 0, ps_oxy_filt4)

#inspect
print(nrow(tax_table(ps_oxy_filt3)))
print(nrow(tax_table(ps_oxy_filt4)))

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt4)
track_oxydf$filt_lowreads <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_lowreads", n_ASVs = ntaxa(ps_oxy_filt4)))




#this should filter out OTU's that have fewer than 10 reads present across all samples
#not sure this is necessary because the relative abundance filtering step would probably catch all of these
#but hey I mean you never know just in case
#To increase the reliability of microbial composition, we advise removing OTUs with 
#<10 copies in individual samples, 
#particularly in studies where only one subsample per specimen is available for analysis. (like this one)
# doi: 10.3389/fcimb.2023.1165295.
#this paper says filtering based on variance decreases reliability compared to the recommended method
filter <- phyloseq::genefilter_sample(ps_oxy_filt4, filterfun_sample(function(x) x >= 10))
ps_oxy_filt5 <- prune_taxa(filter, ps_oxy_filt4)

#inspect
print(nrow(tax_table(ps_oxy_filt4)))
print(nrow(tax_table(ps_oxy_filt5)))

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt5)
track_oxydf$filt_10reads <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_10reads", n_ASVs = ntaxa(ps_oxy_filt5)))



#filtering according to Gillian's advice (remove potential contaminants and errors)
View(tax_table(ps_oxy_filt5))

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
asv_sample_counts <- apply(otu_table(ps_oxy_filt5), 1, function(x) sum(x > 0))
head(asv_sample_counts)

# Identify ASVs that appear in more than one sample
asvs_to_keep <- names(asv_sample_counts[asv_sample_counts > 1])

# Prune the phyloseq object to keep only those ASVs
ps_oxy_filt6 <- prune_taxa(asvs_to_keep, ps_oxy_filt5)

#inspect
print(nrow(tax_table(ps_oxy_filt5)))
print(nrow(tax_table(ps_oxy_filt6)))
#a huuuge number of these asv's only appeared in one sample

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt6)
track_oxydf$filt_1_sample <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_1_sample", n_ASVs = ntaxa(ps_oxy_filt6)))



#remove anything that didn't get assigned to at least preaxostyla
#There are a few ASVs in the oxy sheet that only come out to Excavata and when I 
#blast them the closest hit is Saccinobaculus, but the % identity is low, and it’s 
#only 3 ASVs once we get rid of the ASVs that are only in 1 sample each. For 
#simplicity might as well just ignore (toss) those. Also note I spot checked a bunch 
#of the ASVs that came out as just Eukaryota (in both oxy and para), and they’re just junk, 
#like bits of bacterial genomes and stuff. Not even 18S or 16S. So I think it’s safe 
#to toss anything that doesn’t come out all the way to phylum. You have quite a few 
#random things like bacteria, archaea, fungi, etc., especially in the oxy set.

# Subset taxa where the Class is "Parabasalia"
ps_oxy_filt7<- subset_taxa(ps_oxy_filt6, Class == "Preaxostyla")

#inspect
print(nrow(tax_table(ps_oxy_filt6)))
print(nrow(tax_table(ps_oxy_filt7)))
#a huuuge number of these asv's weren't oxys

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt7)
track_oxydf$filt_preax <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_preax", n_ASVs = ntaxa(ps_oxy_filt7)))




#remove any known contaminants
#For contamination in the oxy set, get rid of anything Pyrsonymphidae and Oxymonadidae.
#All of the known termite contamination falls into one of those two families so 
#that’s pretty straightforward. All the other families are things that might plausibly 
#be in Cryptocercus (and anything Saccinobaculus is 100% from Cryptocercus).


View(tax_table(ps_oxy_filt7))

#ok so the firs thing I'm gonna do is look at the relative abundances of known contaminants
#I'm going to save a merged table with read counts per sample and total read counts

# Extract OTU and taxonomy tables for your object
otus <- as.data.frame(otu_table(ps_oxy_filt7))
taxa <- as.data.frame(tax_table(ps_oxy_filt7))

# Add ASV sequence
otus$ASV <- rownames(otus)
taxa$ASV <- rownames(taxa)

# Merge tax and otu
ps_oxy_merged <- merge(taxa, otus, by = "ASV")
View(ps_oxy_merged)

#remove random irrelevant indexing column
ps_oxy_merged$asv <- NULL

#compute total reads per ASV
ps_oxy_merged$total_abundance <- rowSums(ps_oxy_merged[, sapply(ps_oxy_merged, is.numeric)])

#sort by ASV abundance
ps_oxy_sorted <- ps_oxy_merged[order(-ps_oxy_merged$total_abundance), ]
View(ps_oxy_sorted)




#save as two separate sheets in an excel file
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

addWorksheet(wb, "oxy_contam")
writeData(wb, "oxy_contam", ps_oxy_sorted, na.string = "NA")


# Save workbook
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//Taxonomy_table_oxy_contam.xlsx", overwrite = TRUE)
#went in an manually his the columns with the actual sequences because they were getting in the way
#yay



ps_oxy_filt8 <- subset_taxa(ps_oxy_filt7, !(Family %in% c("Pyrsonymphidae", "Oxymonadidae")))

#inspect
print(nrow(tax_table(ps_oxy_filt7)))
print(nrow(tax_table(ps_oxy_filt8)))
#only a small number were known contaminants


#track reads per sample
read_counts <- sample_sums(ps_oxy_filt8)
track_oxydf$filt_contam <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_contam", n_ASVs = ntaxa(ps_oxy_filt8)))



#all this additional filtering brought several samples read counts lower than 1000
#so gonna remove those
ps_oxy_filt9 <- prune_samples(sample_sums(ps_oxy_filt8) >= 1000, ps_oxy_filt8)
ps_oxy_filt9 <- prune_taxa(taxa_sums(ps_oxy_filt9) > 0, ps_oxy_filt9)

#inspect
print(nrow(tax_table(ps_oxy_filt8)))
print(nrow(tax_table(ps_oxy_filt9)))


#track reads per sample
read_counts <- sample_sums(ps_oxy_filt9)
track_oxydf$filt_lowreads2 <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_lowreads2", n_ASVs = ntaxa(ps_oxy_filt9)))



#remove other contaminants
#remove one additional unidentified ASV that Gillian flagged as a contaminant
#i think she got to that by blasting the sequence and it hit a match with a contaminant

#sequence of ASV to remove (specified by Gillian)
bad_seq <- "CCGCGGTAATTCCAGCTCTAATAGTGTATACTAACATTGCTGCAGTTAAAAAGCTCGTAGTTGAATTTT
    GGGGGGCTATGTTCATTTTGTTGCGGCCTTTGTGGCCGGGCAGAGTATATGGGCTTGTCCTCCTCTCTCTTCCCAAT
    CTTTTGGGAGGAGATAGTTTTCGGGGCATGGTCATCTCCTTCAGTGGTGGATGGGTGTGCGCCGGAGTGTTTACTTT
    GAAGAAAATGGAGTGTTCAAAGCAGGCGTTTCGCAGCT"

bad_seq %in% taxa_names(ps_oxy_filt9)

# Remove the ASV from your phyloseq object
ps_oxy_filt10 <- subset_taxa(ps_oxy_filt9, taxa_names(ps_oxy_filt9) != bad_seq)

#inspect
print(nrow(tax_table(ps_oxy_filt9)))
print(nrow(tax_table(ps_oxy_filt10)))

#this sequence is no longer in the dataset so it must have been filtered out already during one of the other filtering steps


#remove oxymonas and opisthomitus seqeunces
#these were sequences masquerading as valid asvs's based on tax assignment 
#that got contaminant hits when blasted

#look in the better phylogeny section for the code where I found  the bad_asvs

print(bad_asvs)
print(length(bad_asvs))

View(tax)

# Filter tax for sequences that match bad_seqs
bad_seqs <- tax$ASV_sequence[tax$ASV_ID %in% bad_asvs]

# Print the ASV IDs
print(bad_seqs)
print(length(bad_seqs))

# Remove bad sequences from ps_para_filt9
library(phyloseq)
library(phangorn)
ps_oxy_filt10 <- prune_taxa(!(taxa_names(ps_oxy_filt9) %in% bad_seqs), ps_oxy_filt9)

#inspect
print(nrow(tax_table(ps_oxy_filt9)))
print(nrow(tax_table(ps_oxy_filt10)))

#track reads per sample
read_counts <- sample_sums(ps_oxy_filt10)
track_oxydf$filt_sneakycontam <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_sneakycontam", n_ASVs = ntaxa(ps_oxy_filt10)))


#this last filtering step took some samples down below 1000 total reads threshold (i think only 1)
#lowreads removal
#it looks like this additional filtering brought one sample to less than 1000 reads so let's remove that 
ps_oxy_filt11 <- prune_samples(sample_sums(ps_oxy_filt10) >= 1000, ps_oxy_filt10)
#remove any asvs that were only present in that sample
ps_oxy_filt11 <- prune_taxa(taxa_sums(ps_oxy_filt11) > 0, ps_oxy_filt11)

#inspect
print(nrow(tax_table(ps_oxy_filt10)))
print(nrow(tax_table(ps_oxy_filt11)))


#track reads per sample
read_counts <- sample_sums(ps_oxy_filt11)
track_oxydf$filt_lowreads3 <- read_counts[rownames(track_oxydf)]

# Add entries to ASV tracking data frame
track_oxy_asv <- rbind(track_oxy_asv, data.frame(step = "filt_lowreads3", n_ASVs = ntaxa(ps_oxy_filt11)))


#transform sample counts into relative abundances
ps_relative_oxy  = transform_sample_counts(ps_oxy_filt11, function(x) x / sum(x) )
View(otu_table(ps_relative_oxy))

saveRDS(ps_oxy_filt11, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//ps_oxy_filt11.rds")
saveRDS(tax_full, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//tax_full.rds")

View(track_oxydf)
View(track_oxy_asv)


#filter to only the most abundant taxa for small version of data for visualization
#don't try and make plots with more ASV's than this or R freaks out
ps_small_oxy = filter_taxa(ps_relative_oxy, function(x) sum(x) > .1, TRUE)
ps_lesssmall_oxy = filter_taxa(ps_relative_oxy, function(x) sum(x) > .05, TRUE)

print(nrow(tax_table(ps_lesssmall_oxy)))
#195
print(nrow(tax_table(ps_small_oxy)))
#112
#now we have 112 ASV's. I'm going to use this for initial visualization but I will
#probably include more ASV's in my NMDS later

#visualize small filtered
#don't plot the bigger unfiltered data set because laptop can't deal and it doesn't show anyway
p <- plot_bar(ps_small_oxy, fill = "asv")
p

p <- plot_bar(ps_lesssmall_oxy, fill = "asv")
p

p <- plot_bar(ps_relative_oxy, fill = "asv")
p

library(ggplot2)
# Save as PNG
#this code isn't working because I can't get the size right so just export it manually as a png
#ggsave("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//asv_barplot_para.png", plot = p, width = 854, height = 419, dpi = 300, limitsize = FALSE)



#save tracking data as excel files 

# Load the library
library(tibble)

# Move row names to a new column called "Sample"
track_oxy <- rownames_to_column(track_oxydf, var = "Sample")
View(track_oxy)

# Load the package
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "track_oxy")
writeData(wb, sheet = "track_oxy", x = track_oxy)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//track_oxy.xlsx", overwrite = TRUE)



wb <- createWorkbook()
addWorksheet(wb, "track_oxy_asv")
writeData(wb, sheet = "track_oxy_asv", x = track_oxy_asv)

# Save workbook
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//track_oxy_asv.xlsx", overwrite = TRUE)


#ctd----


#make phylogeny----
load("OOD_phylogeny_oxy.RData")
save.image("OOD_phylogeny_oxy.RData")
ps_relative_oxy <- readRDS("ps_relative_oxy.rds")

#save updated ps_relative_oxy phyloseq object
#this object has gillian's contamination edits and updated filtering pipeline
#deleted previous versions with incomplete filtering pipeline and contaminants
saveRDS(ps_relative_oxy, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25////oxy//ps_relative_oxy.rds")
saveRDS(ps_oxy_filt9, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25////oxy//ps_oxy_filt9.rds")
saveRDS(track_oxydf, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25////oxy//track_oxydf.rds")
saveRDS(track_oxy_asv, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25////oxy//track_oxy_asv.rds")

load("OOD_phylogeny_oxy.RData")
save.image("OOD_phylogeny_oxy.RData")
ps_relative_oxy <- readRDS("ps_relative_oxy.rds")
ps_oxy_filt9 <- readRDS("ps_oxy_filt9.rds")
track_oxydf <- readRDS("track_oxydf.rds")
track_oxy_asv <- readRDS("track_oxy_asv.rds")

#save updated ps_relative_oxy phyloseq object
#this object has gillian's contamination edits and updated filtering pipelin
#delted previous versions with incomplete filteirng pipeline and contaminants
saveRDS(ps_relative_oxy, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25////oxy//ps_relative_oxy.rds")


#only include sequences in the phylogeny that are still here after filtering

#get ASV sequences from phyloseq object
seqs_rel <- taxa_names(ps_relative_oxy)

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
#> AIC(fit_JC_rel2, fit_GTR_rel)
#df      AIC
#fit_JC_rel  259 14648.44
#fit_GTR_rel 268 13373.36

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

# Extract tax table
tax <- as.data.frame(tax_table(ps_relative_oxy))

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
custom_labels <- ifelse(
  is.na(genus) | genus != "Saccinobaculus",
  paste0(asv_ids[tip_seqs], " | ", ifelse(is.na(genus), "Unclassified", genus)),
  ""
)

# Assign new labels
tree_label$tip.label <- custom_labels

# Plot
plot(tree_label, show.tip.label = TRUE, cex = 0.6)



#create consistent ASV mapping across all trees and 
#create table with ASV ids, tax info, and seqs

# STEP 1: Prepare taxonomy info and ASV ID mapping
tax <- as.data.frame(tax_table(ps_relative_oxy))
asv_seqs <- taxa_names(ps_relative_oxy)
asv_ids <- paste0("ASV_", seq_along(asv_seqs))
names(asv_ids) <- asv_seqs  # map ASV sequence to ID

# Add ASV metadata
tax$ASV_sequence <- rownames(tax)
tax$ASV_ID <- asv_ids[rownames(tax)]
View(tax)

# Replace NAs
tax$Genus[is.na(tax$Genus)] <- "Unclassified"
tax$Species[is.na(tax$Species)] <- "Unclassified"
tax$Family[is.na(tax$Family)] <- "Unclassified"

tax_full <- tax

#Create labels for tree plots (all ASVs)
tax_full$Tree_Label <- paste0(tax_full$ASV_ID, " | ", tax_full$Genus)
View(tax_full)

# Create labels for tree plots (Genus only)
tax$Tree_Label <- ifelse(
  tax$Genus != "Saccinobaculus",
  paste0(tax$ASV_ID, " | ", tax$Genus),
  ""  # blank for Saccinobaculus
)

View(tax)

View(tax_full)

# STEP 2: write a function to Relabel tips for each tree
label_tree <- function(tree, tax_df, label_col = "Tree_Label") {
  tree_copy <- tree
  tip_seqs <- tree_copy$tip.label
  tree_copy$tip.label <- tax_df[tip_seqs, label_col]
  return(tree_copy)
}

tree_label_rooted <- label_tree(tree_rooted, tax)
tree_label_JC     <- label_tree(tree_rooted_JC, tax)
tree_label_rel    <- label_tree(tree_rooted_rel, tax)

# STEP 3: Plot trees with genus-only labels
par(mfrow = c(1, 3))  # side-by-side layout
par(mfrow = c(1, 1), mar = c(1, 1, 4, 1)) #reset layout to one plot per window
plot(tree_label_rooted, show.tip.label = TRUE, cex = 0.6, main = "GTR")
plot(tree_label_JC,     show.tip.label = TRUE, cex = 0.6, main = "JC")
plot(tree_label_rel,    show.tip.label = TRUE, cex = 0.6, main = "NJ")

# STEP 4: Create final table including species info (for non-Saccinobaculus only)
tree_table <- tax[tax$Genus != "Saccinobaculus", 
                   c("ASV_ID", "ASV_sequence", "Genus", "Species", "Tree_Label")]
View(tax)

#ok so these trees are all dealing with monocercomonoides, blattamonas, and unclassified ASV's very differently
# I find it weird that the monocercomonoides genus is getting consistently assigned to two different clades


#first make a tree of just the moncercomonoides, blattamonas, and other unclassified seqeunces
# Assuming your table is named `final_table`
asv_names <- tree_table$Tree_Label
asv_seqs <- tree_table$ASV_sequence

# Create DNAStringSet
asv_dna <- DNAStringSet(asv_seqs)
names(asv_dna) <- asv_names

# Write to FASTA
writeXStringSet(asv_dna, filepath = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//monocercomonoides_subset_oxy.fasta")

#make alignment
aligned_subset <- AlignSeqs(asv_dna, anchor = NA)
names(aligned_subset)

library(phangorn)
# Read aligned sequences into phyDat format
phangAlign_sub <- phyDat(as(aligned_subset, "matrix"), type="DNA")

# Compute pairwise distances using maximum likelihood
dm_sub <- dist.ml(phangAlign_sub)

# Build a neighbor-joining tree
tree_NJ_sub <- NJ(dm_sub)
class(tree_NJ_sub)

tree_NJ_sub_rooted <-midpoint(tree_NJ_sub)

head(tree_NJ_sub_rooted$tip.label)

# Fit JC model
fit_JC_sub <- pml(tree_NJ_sub, phangAlign_sub)
fit_JC_sub <- optim.pml(fit_JC_sub, model = "JC", rearrangement = "stochastic") 
#this is the first time I've used rearrangement = stochastic and chat gpt says its slower but make better fit trees

fit_GTR_sub <- pml(tree_NJ_sub, phangAlign_sub)
fit_GTR_sub <- optim.pml(fit_GTR_sub, model = "GTR", rearrangement = "stochastic")


AIC(fit_JC_sub, fit_GTR_sub)
#fit_JC_sub  61 7995.449
#fit_GTR_sub 70 7359.054

#GTR model has a lower AIC so we'l go with that one

#check if the trees I plan on using are rooted
is.rooted(fit_GTR_sub$tree)


#the tree is not rooted, so we need to root it now
#because we do not have an outgroup, we will estimate the root using midpoint rooting
#this method places the root in the middle between the two longest branches

# If your tree is stored in a pml object (e.g., fit_GTR_rel$tree)
tree_unrooted_sub <- fit_GTR_sub$tree
head(tree_unrooted_sub$tip.label) 

# Root using midpoint method
tree_rooted_sub <- midpoint(tree_unrooted_sub)
head(tree_rooted_sub$tip.label) 

#now gnerate bootstrap values 
bs_sub <- bootstrap.pml(fit_GTR_sub, bs = 100, optNni = TRUE, multicore = TRUE, mc.cores = 4)
head(bs_sub$tip.label) 

#ok. everything has tip labels
  
#prepare named vector for tip labels
tip_labels <- tree_table$Tree_Label
names(tip_labels) <- tree_table$ASV_ID

#check that all tip labels exist in names tip labels
#should return TRUE
all(tree_rooted_sub$tip.label %in% names(tip_labels))  # Should be TRUE

#rename main tree tips with tax adn ASV ID, only do this once
tree_rooted_sub$tip.label <- tip_labels[tree_rooted_sub$tip.label]
head(tree_rooted_sub$tip.label) 

print(colnames(tree_table))
# Rename bootstrap tree tips first, BEFORE midpoint rooting
bs_sub_renamed <- lapply(bs_sub, function(tree) {
  tree$tip.label <- tip_labels[tree$tip.label]  # replaces "ASV_5" with "ASV_5 | Genus"
  tree
})
class(bs_sub_renamed) <- "multiPhylo"
head(bs_sub_renamed[[1]]$tip.label) # this is a list of phylo objects so index the first one

#root the bs trees
bs_sub_rooted <- lapply(bs_sub_renamed, midpoint)
class(bs_sub_rooted) <- "multiPhylo"

# Then plot
plotBS(tree_rooted_sub, bs_sub_rooted, p = 0, type = "phylogram", cex = 0.6)

#Here's a quick interpretation guide:
#≥ 95%: Very strong support
#85–94%: Strong support
#70–84%: Moderate support
#< 70%: Weak support (often considered unreliable)
#So, if a clade has a bootstrap value of 98, 
#it's very likely that grouping is biologically meaningful. 
#If it's 55, you should be cautious about interpreting that split.

#so it's looking like there is a lot of confidence in the more external nodes but less confidence in the more basal nodes
#which makes sense. I think. Given the weirdness with monocercomonoides.


#now make a version of the original GTR tree with bootstrapping
#generate bootstrap values 
bs <- bootstrap.pml(fit_GTR_rel, bs = 100, optNni = TRUE, multicore = TRUE, mc.cores = 4)
head(bs[[1]]$tip.label) # this is a list of phylo objects so index the first one 
#ok so tip labels are seqeunces for the bootstrapped trees
head(fit_GTR_rel$tip.label)
#and this has no tip labels.
head(tree_rooted$tip.label)

library(phangorn)

#root the GTR tree
tree_unrooted <- fit_GTR_rel$tree
tree_rooted <- midpoint(tree_unrooted)


#Create the tip label mapping from your tree_table
View(tree_table)
tip_labels <- tree_table$Tree_Label
names(tip_labels) <- tree_table$ASV_sequence  # assuming this matches what's in the tree

#Rename tips in the bootstrap trees before rooting:
#only rename tip labels for seqquences represented in tree_table
#the non-sacc seqs
bs_renamed <- lapply(bs, function(tree) {
  tree$tip.label <- ifelse(
    tree$tip.label %in% names(tip_labels),
    tip_labels[tree$tip.label],
    ""  # blank label
  )
  tree
})
class(bs_renamed) <- "multiPhylo"


#Midpoint-root the renamed bootstrap trees:
#bs_rooted <- lapply(bs_renamed, midpoint)
#class(bs_rooted) <- "multiPhylo"

#relabel the main tree:
#and leave labels blank if not in tree_table
tree_rooted$tip.label <- ifelse(
  tree_rooted$tip.label %in% names(tip_labels),
  tip_labels[tree_rooted$tip.label],
  ""  # blank for ASVs not in tree_table
)

#plot without bootstraps to confirm its working
plot(tree_rooted, type = "phylogram", cex = 0.6)
# Convert to ggtree object for readability
library(ggtree)
p <- ggtree(tree_rooted) +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.3)

p
#try this for readability also 
# Convert to ggtree object
plot(tree_rooted, type = "phylogram", cex = 0.5, label.offset = 0.2)

# Generate a majority-rule consensus tree
cons_tree <- consensus(bs_rooted, p = 0.5)

# Optional: map support values to original tree
tree_with_support <- plotBS(tree_rooted, bs_rooted, type = "none", p = 50)

# Plot the consensus tree with support values
pdf("tree_with_consensus_support.pdf", width = 12, height = 18)
plot(tree_with_support, show.node.label = TRUE, cex = 0.6)
dev.off()

#ok it looks like this isn't working because R doesn't have enough memory
#it's a whole thing not gonna get into it right now
#gonna go blast those sequences instead

pdf("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//tree_bootstrap_oxy.pdf", width = 12, height = 18)

plotBS(tree_rooted, bs_rooted, p = 50, type = "phylogram", cex = 0.6)
dev.off()

#blasted the sequences
#some of them may have been miss classified in the tax assignment? need to ask Gillian
#for example some things that were unclassfied in the tax assignment are getting assigned to Saccinobaculus by blast
#but after a cursory look through it doesn't look like they are mistakes
#so I'm going to keep the GTR tree as is and use that in my downstream analysis

#add GTR rooted tree to ps_relative_oxy
#make sure to regenerate tree_rooted before doing this because we messed with tip labels
ps_relative_oxy_tree <- phyloseq(otu_table(ps_relative_oxy), 
                                 tax_table(ps_relative_oxy), 
                                 tree = tree_rooted)

#check that the tree was successfully added
phy_tree(ps_relative_oxy_tree)
#yes! looks good

View(tree_table)


#plot tree with outlier sequences that contributed the most heavily to
# between treatment differences in weighted unifrac PERMANOVA analysis

#check that tip labels are ASV seqs
tree_rooted$tip.label
#they are not, we relabeled them. so let's go back to tree unrooted

tree_unrooted$tip.label
#yes! these are ASV seqs. so let's go from here

#midponit root the tree
tree_rooted_outlier <- midpoint(tree_unrooted)
tree_rooted_outliers$tip.label


# Get the current tip labels (these are ASV sequences)
original_tip_labels <- tree_rooted_outlier$tip.label

# Create a named vector to map sequences to informative labels (e.g., "ASV_001_Bacillus")
# This assumes you have a data frame `tax_full` with columns: ASV_sequence and Tree_Label
label_map <- setNames(tax_full$Tree_Label, tax_full$ASV_sequence)

# Get the set of ASV sequences that are significant
sig_seqs <- rownames(signif_taxa_sig)

# Initialize blank labels for all tips
tip_labels_outliers <- rep("", length(original_tip_labels))
names(tip_labels_outliers) <- original_tip_labels  # names = ASV sequences

# Identify the subset of sequences that are both significant and in the tree
matching_seqs <- intersect(sig_seqs, names(tip_labels_outliers))

# Assign new labels for matching tips using label_map
tip_labels_outliers[matching_seqs] <- label_map[matching_seqs]

# Apply updated labels to the tree
tree_rooted_outlier$tip.label <- tip_labels_outliers

# Optional: verify the result
head(tree_rooted_outlier$tip.label[tree_rooted_outlier$tip.label != ""])

#check
head(tree_rooted_outlier$tip.label)
#ok this worked



#now plot
# Load required libraries
library(ggtree)
library(treeio)
library(ggplot2)

# Build ggtree object
p <- ggtree(tree_rooted_outlier) + 
  geom_tree() +
  theme_tree2() +
  ggtitle("Phylogenetic Tree with Significant ASV Labels")

# Access ggtree's internal data frame
p_data <- p$data

# Modify only significant labels (i.e., those not empty)
p_data$label[p_data$isTip & p_data$label == ""] <- NA  # hide non-significant labels

# Re-plot using modified data
p %+% p_data + geom_tiplab(align = TRUE, linetype = "dashed", size = 3)


#simpler plot
plot(tree_rooted_outlier, show.tip.label = TRUE, cex = 0.6)


#after blasting the moncercomonoides, unclassified, and blattamonas taxa in this dataset
#it looks like many are getting hits for opithomitus
#which gillian said is probably from a termite
#additionally, several of the most significant contributors to the funky outliers in the weighted unifrac data
#are from that subset
#so I want to remove anything getting opithomitus hits on blast
#but when I download blast hits from ncbi the tax name of the hits doesn't download
#so I'm going to download them and try to access tax names using the accession number

#read in the blast results
# Set the path to your CSV file
monocercomonoides_subset_blast <- "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//monocercomonoides_subset_blast.csv"  # Change this if your file has a different name
"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//monocercomonoides_subset_blast.csv"

# add colnames (colnames didn't come with the download but chat gpt generated these based on standard output)
blast_colnames <- c(
  "qseqid",    # Query sequence ID
  "sseqid",    # Subject sequence ID
  "pident",    # % identity
  "length",    # Alignment length
  "mismatch",  # Number of mismatches
  "gapopen",   # Number of gap openings
  "qstart",    # Start of alignment in query
  "qend",      # End of alignment in query
  "sstart",    # Start of alignment in subject
  "send",      # End of alignment in subject
  "evalue",    # E-value
  "bitscore"   # Bit score
)

# Read the CSV file (no header)
mono_blast_df <- read.csv(monocercomonoides_subset_blast, header = FALSE, stringsAsFactors = FALSE)

# Assign column names
colnames(mono_blast_df) <- blast_colnames
View(mono_blast_df)

#now I need to use accession numbers to find the tax names for each hit
#but I'm going to do this in open on demand, so I need to export the file first
library(writexl)

write_xlsx(mono_blast_df, path = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//monocercomonoides_subset_blast.xlsx")

library(readxl)

# Read the Excel file to open on deman (assuming it's in your working directory)
mono_blast_df <- read_excel("monocercomonoides_subset_blast.xlsx")


library(rentrez)

# Get unique accession numbers (to avoid redundant queries)
accessions <- unique(mono_blast_df$sseqid)

# Look up organism names
get_organism_name <- function(acc) {
  # Fetch GenBank record summary
  tryCatch({
    summary <- entrez_summary(db = "nucleotide", id = acc)
    summary$organism  # Extract the organism name
  }, error = function(e) {
    NA  # Return NA if the accession fails
  })
}

# Run lookup (may take a few minutes if many accessions)
organism_names <- sapply(accessions, get_organism_name)

# Create a lookup table
organism_df <- data.frame(sseqid = accessions, organism = organism_names, stringsAsFactors = FALSE)

#merge organims names into table
mono_blast_wnames <- merge(mono_blast_df, organism_df, by = "sseqid", all.x = TRUE)
View(mono_blast_wnames)

library(writexl)
library(readxl)

write_xlsx(mono_blast_wnames, path = "mono_blast_wnames.xlsx")

mono_blast_wnames <- read_excel("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//mono_blast_wnames.xlsx")

View(mono_blast_wnames)

#this version has the top 10 hits for each asv 
#now I want to find all the ASV IDs that have opisthomitus or oxymonas in the organism column

bad_asvs <- mono_blast_wnames %>%
  filter(grepl("Oxymonas|Opisthomitus", organism, ignore.case = TRUE)) %>%
  pull(qseqid)
print(bad_asvs)

library(dplyr)

mono_blast_wnames <- mono_blast_wnames %>%
  arrange(qseqid)

bad_asvs <- unique(bad_asvs)
print(bad_asvs)



#remove oxymonas and poisthomitus
#everything in bad_asvs
#the remake tree and redo analysis and see how shee looks


#ok now that those problem seqeunces have been removed, I need to make a new phylogeny to reflect the filtered data

#get ASV sequences from phyloseq object
seqs_rel2 <- taxa_names(ps_relative_oxy)

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
fit_JC_rel2 <- optim.pml(update(fit_rel2, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#ran wayyy quicker with only 131 seqs basically instantaneous on OOD
#also basically instant on laptop yay! 131 seqs is so muh better than 5000
#15 minutes if doing stochastic rearrangement

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel2 <- optim.pml(update(fit_rel2, model = "GTR", k = 4), 
                         optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#also basically instantaneous with 768 GB memory
#15 minutes if doing stochastic rearrangement

AIC(fit_JC_rel2, fit_GTR_rel2)
#> AIC(fit_JC_rel2, fit_GTR_rel2)
#df      AIC
#fit_JC_rel2  237 11049.74
#fit_GTR_rel2 246 10056.82

#GTR model has a lower AIC so we'l go with that one

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

#now plot with tip labels
# Create custom tip labels like: ASV_10 | Oxymonas
tax$ASV_labels <- paste(tax$ASV_ID, tax$Genus, sep = " | ")

# Create a named vector: names = sequences, values = new labels
label_map <- setNames(tax$ASV_labels, tax$ASV_sequence)

tree_rooted2$tip.label <- ifelse(tree_rooted2$tip.label %in% names(label_map),
                                 label_map[tree_rooted2$tip.label],
                                 tree_rooted2$tip.label)

plot(tree_rooted2, cex = 0.5)

#ok, it looks like there are still two ASV's here with that long branch
# i guess they weren't visually accessible on the previous version of the plot
# so I missed them but we still need to remove them

#11 and 128 so going back to add those to the list of bad seqs in the filter sand wrangle sample data section

#replotted without those and now there isn't one log weird branch but the other branches resolved a lot differently
#there are two clades for saccinolbacculas but I'm just going to leave it like that
#I don't trust tree models
# also need to go back and calculate branch lengths and remove that one sample with very few asvs and remake everything from there
#if removing it removed any other asvs


#ugh ok, so I removed that one lowread sample. There were no unique asvs for that sample,
#so the tree shouldn't have changed. But I remade it anyway just to be safe and it looks wayyy different now
#which it shouldn't!
#i think we have a problem of it just can't resolve major clades based on the data
# so we need to include long read sequences
#i'm not even going to bother reanalyzing that until Gillian sends be the long read stuff 
#because I don't really trust the tree at all.

#given that I'm also not totally sure that 11 and 128 need to be removed either so 
#reserving judgement on that
#they don't meet any other elimination criteria, just appeared on a long branch in one version of a tree
#if they stay then that one lowreads sample can also stay which would be clutch. we'll see



rm(list=ls())

tax_full <- readRDS("tax_full.rds")
ps_oxy_filt10 <- readRDS("ps_oxy_filt10.rds")
#ok remakeing tree now that sneaky contam has been removed

#get ASV sequences from phyloseq object
seqs_rel3 <- taxa_names(ps_oxy_filt11)


#propogates labels to tips of the tree
names(seqs_rel3) <- seqs_rel3

#make alignment
library(DECIPHER)
aligned_rel3 <- AlignSeqs(DNAStringSet(seqs_rel3), processors = 4, anchor=NA,verbose=TRUE)
names(aligned_rel3) #check that ASV seqs are names in the alignment

library(phangorn)
# Read aligned sequences into phyDat format
phangAlign_rel3 <- phyDat(as(aligned_rel3, "matrix"), type="DNA")

# Compute pairwise distances using maximum likelihood
dm_rel3 <- dist.ml(phangAlign_rel3, processors = 1, verbose = TRUE)

# Build a neighbor-joining tree
tree_NJ_rel3 <- NJ(dm_rel3)

class(tree_NJ_rel3)
# Convert to pml object
fit_rel3 <- pml(tree_NJ_rel3, data = phangAlign_rel3)


# Fit the JC69 model
#fit_JC_rel2 <- optim.pml(update(fit_rel2, model = "JC"), rearrangement = "stochastic", optNni = TRUE) 
#only doing GTR model now cause its always better

# Fit the GTR model with Gamma rate heterogeneity
fit_GTR_rel3 <- optim.pml(update(fit_rel3, model = "GTR", k = 4), 
                          optNni = TRUE, optBf = TRUE, optQ = TRUE, optGamma = TRUE, rearrangement = "stochastic")
#also basically instantaneous with 768 GB memory
#15 minutes if doing stochastic rearrangement

#check if the trees I plan on using are rooted
is.rooted(fit_GTR_rel3$tree)

#the tree is not rooted, so we need to root it now
#because we do not have an outgroup, we will estimate the root using midpoint rooting
#this method places the root in the middle between the two longest branches

# If your tree is stored in a pml object (e.g., fit_GTR_rel$tree)
tree_unrooted3 <- fit_GTR_rel3$tree

# Root using midpoint method
tree_rooted3 <- midpoint(tree_unrooted3)

# Plot to check
plot(tree_rooted3, show.tip.label = FALSE)

#now plot with tip labels

# Create a named vector: names = sequences, values = new labels
label_map <- setNames(tax_full$Tree_Label, tax_full$ASV_sequence)

tree_rooted_label <- tree_rooted3

tree_rooted_label$tip.label <- ifelse(tree_rooted_label$tip.label %in% names(label_map),
                                      label_map[tree_rooted_label$tip.label],
                                      tree_rooted_label$tip.label)

plot(tree_rooted_label, cex = 0.5)

#ok you know what on second thoughts, I think this tree looks just fine
#there's a long branch but that just looks like where all the monocercomonoides and blattamonas got placed
#they are different from everything else so it makes sense they should be in their own claade
#they're fine right there and no more filtering required I think.
# we should be all good to go with analysis from here
#going back to make a relative abundance oxy ps from ps_oxy_filt10
#and that 's the one we'll be adding the tree to

load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//OOD_phylogeny_oxy3.RData")

#add tree to phyloseq object
ps_relative_oxy_tree <- ps_relative_oxy
ps_relative_oxy_tree <- phy_tree(ps_relative_oxy_tree) <- tree_rooted3

#check
phy_tree(ps_relative_oxy_tree)

saveRDS(ps_relative_oxy_tree, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//ps_relative_oxy_tree.rds")
saveRDS(tax_full, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//tax_df_oxy.rds")



#ctd----


#PERMANOVA (bray curtis)----
#extract asv table

meta_oxy <- data.frame(sample_data(ps_relative_oxy))

#add treatment column to metadata
meta_oxy$treatment <- paste(meta_oxy$type_original_log, meta_oxy$type_log_new, sep = ":")
meta_oxy$treatment <- gsub("wood", "", meta_oxy$treatment)
meta_oxy$type_original_log <- gsub("wood", "", meta_oxy$type_original_log)
(View(meta_oxy))

print(unique(meta_oxy$treatment))
#meta data missing site info for two samples
#add it back manually
meta_oxy["49F1","site"] <- "pond_drain"
meta_oxy["49G1","site"] <- "pond_drain"

sample_data(ps_relative_oxy) <- sample_data(meta_oxy)

bray_oxy <- phyloseq::distance(ps_relative_oxy, method = "bray")


library(vegan)
#full models
permanova_oxy_bray <- adonis2(bray_oxy ~ Dissected + site/log + Instar + swapped*type_original_log, 
                                data = meta_oxy, 
                                permutations = 999, 
                                method = "bray",
                                by = "terms")
print(permanova_oxy_bray)
str(permanova_oxy_bray$aov.tab)
class(permanova_oxy_bray)

#neither swapped nor original log type significant on their own but together there is an interaction effect


#let's check out the dispersion stuff
#check homogeneity of variances for site
disp1.1 <- betadisper(bray_oxy, group=meta_oxy$treatment)
permutest(disp1.1)
# p value > 0.05,  homogeneous variances 
#location but not variance change across levels of swapped*type_log_original


#visualization
library(vegan)
nmds_bray <- metaMDS(bray_oxy, k = 3, trymax = 100)
nmds_bray$stress
#0.0.1660138 acceptable!

#extract scores
nmds_bray_points <- as.data.frame(nmds_bray$points)
View(nmds_bray_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_bray_points$SampleID <- rownames(nmds_bray_points)
#remove scores previously added to meta_rel
View(meta_oxy)
#meta_rel <- meta_rel[, -c(35:37)]
meta_oxy$SampleID <- rownames(meta_oxy)


# Now merge by SampleID
nmds_bray_meta <- merge(nmds_bray_points, meta_oxy, by = "SampleID")
View(nmds_bray_meta)

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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//bray_nmds_treat_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//bray_nmds_log_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//bray_nmds_swapped_oxy.html", selfcontained = TRUE)



#ctd----

#PERMANOVA (jaccard)----

# Compute jaccard
jaccard_oxy <- phyloseq::distance(ps_relative_oxy, method = "jaccard", binary = TRUE)


#we need meta data along with distance matrix to do pERmANOVA 
#we already have meta small and meta rel from our previous analysis
#actually meta small and meta rela should be the same data frame so we can use the interchangeably


#full models
permanova_oxy_jaccard <- adonis2(jaccard_oxy ~ Dissected + site/log + Instar + swapped*type_original_log, 
                                   data = meta_oxy, 
                                   permutations = 999, 
                                   method = "jaccard",
                                   by = "terms")
print(permanova_oxy_jaccard)
#eep no significant treatment effects


#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp2.2 <- betadisper(jaccard_rel, group=meta_rel$treatment)
permutest(disp2.2)
# p value > 0.05, homogeneous variances 
#location but not variances change across treatment



#visualization
library(vegan)
nmds_jaccard <- metaMDS(jaccard_oxy, k = 3, trymax = 100)
nmds_jaccard$stress
#0.2044216 kind of high take visualization with a grain of salt

#extract scores
nmds_jaccard_points <- as.data.frame(nmds_jaccard$points)
View(nmds_jaccard_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_jaccard_points$SampleID <- rownames(nmds_jaccard_points)
#remove scores previously added to meta_rel
View(meta_oxy)
#meta_rel <- meta_rel[, -c(35:37)]
meta_oxy$SampleID <- rownames(meta_oxy)

# Now merge by SampleID
nmds_jaccard_meta <- merge(nmds_jaccard_points, meta_oxy, by = "SampleID")
View(nmds_jaccard_meta)


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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//jaccard_nmds_treat_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//jaccard_nmds_log_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//jaccard_nmds_swapped_oxy.html", selfcontained = TRUE)


#ctd ----


#PERMANOVA (unifrac unweighted)----

# Compute Unweighted UniFrac
#this considers only the presence or absence of taxa.
#it measures the fraction of unique branch length in the phylogenetic tree that is exclusive to each community.
#Sensitive to rare taxa because all taxa are treated equally.
#Useful when you're interested in whether communities share the same taxa, regardless of how abundant they are.
ps_relative_oxy_tree <- merge_phyloseq(ps_relative_oxy, phy_tree(tree_rooted3))
sample_data(ps_relative_oxy_tree) <- sample_data(meta_oxy)

unifrac_unweighted_rel <- phyloseq::distance(ps_relative_oxy_tree, method = "unifrac", weighted = FALSE)



library(vegan)
permanova_rel_UU <- adonis2(unifrac_unweighted_rel ~ Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_oxy, 
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

#extract scores
nmds_UU_points <- as.data.frame(nmds_UU$points)
View(nmds_UU_points)
#merge with metadata for plotting
# First, add the rownames as a column in each data frame
nmds_UU_points$SampleID <- rownames(nmds_UU_points)
#remove scores previously added to meta_rel
View(meta_oxy)
#meta_rel <- meta_rel[, -c(35:37)]
meta_oxy$SampleID <- rownames(meta_oxy)

# Now merge by SampleID
nmds_UU_meta <- merge(nmds_UU_points, meta_oxy, by = "SampleID")
View(nmds_UU_meta)


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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//UU_nmds_treat_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//UU_nmds_log_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//UU_nmds_swapped_oxy.html", selfcontained = TRUE)




#ctd ----

#PERMANOVA (unifrac weighted)----
# Compute Weighted UniFrac
#considers both presence/absence and relative abundance of taxa.
#measures the fraction of branch length that is weighted by the difference in relative abundance between communities.
#Sensitive to abundant taxa.
#Gives more weight to dominant taxa and downplays rare taxa.

unifrac_weighted_rel <- phyloseq::distance(ps_relative_oxy_tree, method = "unifrac", weighted = TRUE)


library(vegan)
#let's look at the control effects first

#full models
permanova_rel_WU <- adonis2(unifrac_weighted_rel ~  Dissected + site/log + Instar + swapped*type_original_log, 
                            data = meta_oxy, 
                            permutations = 999,
                            by = "terms")

print(permanova_rel_WU)
#significant effects of type_original log and interaction
#weird that some of the blocking variables got close rto significant in the full model



#let's check out the dispersion stuff
#check homogeneity of variances for treatment
disp3.3 <- betadisper(unifrac_weighted_rel, group=meta_oxy$treatment)
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


#extract scores
nmds_WU_points <- as.data.frame(nmds_WU$points)
nmds_WU_points$SampleID <- rownames(nmds_WU_points)

View(meta_oxy)

# Now merge by SampleID
nmds_WU_meta <- merge(nmds_WU_points, meta_oxy, by = "SampleID")


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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//WU_nmds_treat_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//WU_nmds_log_oxy.html", selfcontained = TRUE)



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
htmlwidgets::saveWidget(fig, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//WU_nmds_swapped_oxy.html", selfcontained = TRUE)


#ctd----


#multiple tests corrections####
#first we need to adjust global p-values across distance metrics
# I'm testing the same hypothesis (treatment effect) under 4 lenses
#bray-curtis, jaccard, weighted unifrac, unweighted unifrac,

# Named list of your models
models <- list(
  bray = permanova_oxy_bray,
  jaccard = permanova_oxy_jaccard,
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


#extract the 6 pairwise treatment comparisons:

treatment_levels <- unique(meta_oxy$treatment)
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
distance_matrices <- list(
  bray = bray_oxy,
  w_unifrac = unifrac_weighted_rel
)

# Run pairwise comparisons
for (metric in names(distance_matrices)) {
  dist_matrix <- distance_matrices[[metric]]
  
  res <- lapply(pairwise_comparisons, function(pair) {
    p <- run_pairwise_permanova(dist_matrix, meta_oxy, "treatment", control_vars, pair)
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

# Combine all into one data frame
pairwise_results <- do.call(rbind, results_list)

# bonferroni adjustment within each metric
pairwise_results$p_adj_bonf <- ave(pairwise_results$p_value, pairwise_results$metric, FUN = p.adjust, method = "bonf")
View(pairwise_results)

#Ok cool. Now that I've done that I need to do pairwise comparisons for my treatment 
#variable. For each treatment, I have 4 levels, and I want to compare between all 
#6 possible level combinations. The treatment and sample info is in a table called meta_oxy. 
#With the fdr correction, I got a significant effect of treatment for bray, weighted unifrac, 
#nd unweighted unifrac so I want to do pairwise comparisons for those community metrics. 
#The ordinations are called bray_oxy, unifrac_unweighted_rel, and unifrac_weighted_rel.
#For each pairwise test, I also want to include control variables Dissected, site/log, 
#and Instar in the model in that order. I want a way to organize all the results of 
#the pairwise comparisons together into a table, or at least 3 tables, 
#one for each community metric. In the results tables I want both p values and fdr 
#adjusted p values. 


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

View(results_df)
View(pairwise_results)

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

# Apply to Pairwise Results sheet
if ("p_value" %in% names(pairwise_results)) {
  pairwise_results$p_value <- add_sig_stars(pairwise_results$p_value)
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
View(pairwise_results)

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


# Print result
print(summary_table)
View(summary_table)

print(summary_table)
print(class(summary_table))

print(results_df)
print(class(results_df))

print(pairwise_results)
print(class(pairwise_results))


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

# Print to verify
print(results_df)
print(pairwise_results)
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


# Clean, round, and re-add stars for pairwise_results
pairwise_results <- pairwise_results %>%
  mutate(
    # Remove existing stars and convert to numeric
    p_value = as.numeric(str_remove(p_value, "\\*+$")),
    p_adj_bonf = as.numeric(str_remove(p_adj_bonf, "\\*+$"))
  ) %>%
  # Round numeric columns to 2 significant figures
  mutate(across(c(p_value, p_adj_bonf), ~ signif(.x, 2))) %>%
  # Re-add significance stars
  mutate(
    p_value = paste0(p_value, add_stars(p_value)),
    p_adj_bonf = paste0(p_adj_bonf, add_stars(p_adj_bonf))
  )

#reorder rows of summary table
# Desired order vector
desired_order <- c("original diet", "diet change", "diet change x original diet")

# Reorder rows by matching the treatment column to desired_order
summary_table <- summary_table[match(desired_order, summary_table$treatment), ]

print(summary_table)
View(summary_table)

# Print to verify
print(results_df)
print(summary_table)
print(pairwise_results)

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
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//permanova_all_results_oxy.xlsx", overwrite = TRUE)

#ctd----


#diversity analysis----
library(phyloseq)

# Calculate shannon, and simpsons diversity for each samples
diversity_oxy <- estimate_richness(ps_relative_oxy, measures = c("Shannon", "Simpson"))
View(diversity_oxy)

View(meta_oxy)
#this code cannot calculate richness - 
#it get's tripped up by relative abundances rather than raw counts
#throws this warning if I try to include richness
#Warning message:
  #The data you have provided does not have any singletons. This is highly suspicious. Results of richness
  #estimates (for example) are probably unreliable, or wrong, if you have already
  #trimmed low-abundance taxa from the data. We recommended that you find the un-trimmed data and retry.

#so let's calculate richness manually by counting asvs with abundance > 0 for each sample

# Extract the ASV/OTU table (samples as rows)
otu_oxy <- as(otu_table(ps_relative_oxy), "matrix")

# Make sure samples are rows and transpose if needed
if (!taxa_are_rows(ps_relative_oxy)) {
  otu_oxy <- t(otu_oxy)
}

View(otu_oxy)

# Count number of ASVs with abundance > 0 per sample
richness_oxy <- colSums(otu_oxy > 0)

# Look at the result
head(richness_oxy)



# Add diversity metrics to sample_data

sample_data(ps_relative_oxy) <- cbind(sample_data(ps_relative_oxy), diversity_oxy)
sample_data(ps_relative_oxy) <- cbind(sample_data(ps_relative_oxy), data.frame(richness_oxy))
View(sample_data(ps_relative_oxy))

#plot raw diversity metrics by treatment
library(phyloseq)
library(ggplot2)
library(tidyr)
library(dplyr)

# Extract sample metadata from phyloseq object
meta_oxy <- data.frame(sample_data(ps_relative_oxy))
View(meta_oxy)

#rename richness column
meta_oxy <- meta_oxy %>%
  rename(Richness = richness_oxy)

#add SampleID column
meta_oxy$SampleID <- rownames(meta_oxy)

library(ggplot2)
library(dplyr)
library(tidyr)
#pivot data to long format
meta_long <- meta_oxy %>%
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

remove.packages("glmmTMB")
install.packages("glmmTMB", dependencies = TRUE)
install.packages("glmmTMB", dependencies = TRUE)

install.packages("glmmTMB")
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

hist(meta_oxy$Simpson)
#data is heavily right skewed
meta_oxy$log_Simpsons <- log(meta_oxy$Simpson)
hist(meta_oxy$log_Simpson)

#trying a gamma distribution because chat gpt said it was good for positive right skewed data
#Gaussian (normal) distribution didn't pass the dHARMA tests

glm_simp_oxy <- glmmTMB(as.numeric(Simpson) ~ Dissected + site/log + Instar + swapped*type_original_log,
                    family = Gamma(link = "log"),
                    data = meta_oxy)
summary(glm_simp_oxy) 
aov_simp_oxy <- Anova(glm_simp_oxy) 
aov_simp_oxy

dharma.test(glm_simp_oxy,n=1000)
#assumptions test looks good!
# so no significant effects 

packageVersion("glmmTMB")

#test shannon

hist(meta_oxy$Shannon)
#less heavily right skewed than simpsons
#using a normal dist for this data

glm_shan_oxy <- glmmTMB(as.numeric(Shannon) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_oxy)
summary(glm_shan_oxy) 
aov_shan_oxy <- Anova(glm_shan_oxy) 
aov_shan_oxy

dharma.test(glm_shan_oxy,n=1000)
#assumptions test looks good!
# so no significant effects 

pairwise_shan_oxy_instar <- emmeans(glm_shan_oxy, pairwise ~ Instar, adjust = "tukey")
pairwise_shan_oxy_instar$contrasts
#3 and 4 are different from A but not each other


#test richness

hist(meta_oxy$Richness)
#this one is a scooch left skewed
#which I would expect for richness data

glm_rich_oxy <- glmmTMB(as.numeric(Richness) ~ Dissected + site/log + Instar + swapped*type_original_log,
                        family =  gaussian, 
                        data = meta_oxy)
summary(glm_rich_oxy) 
aov_rich_oxy <- Anova(glm_rich_oxy) 
aov_rich_oxy

dharma.test(glm_rich_oxy,n=1000)
#assumptions test looks good!
# so no significant effects 

pairwise_rich_oxy_instar <- emmeans(glm_rich_oxy, pairwise ~ Instar, adjust = "tukey")
pairwise_rich_oxy_instar$contrasts
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
    table = aov_simp_oxy
  ),
  Shannon = list(
    formula = "Shannon ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_shan_oxy
  ),
  Richness = list(
    formula = "Richness ~ Dissected + site/log + Instar + swapped*type_original_log",
    table = aov_rich_oxy
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
saveWorkbook(wb, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//diversity_stats_tables_oxy.xlsx", overwrite = TRUE)




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
  get_emm_df(glm_simp_oxy, "treatment", "Simpson"),
  get_emm_df(glm_simp_oxy, "swapped", "Simpson"),
  get_emm_df(glm_simp_oxy, "type_original_log", "Simpson"),
  
  get_emm_df(glm_shan_oxy, "treatment", "Shannon"),
  get_emm_df(glm_shan_oxy, "swapped", "Shannon"),
  get_emm_df(glm_shan_oxy, "type_original_log", "Shannon"),
  
  get_emm_df(glm_rich_oxy, "treatment", "Richness"),
  get_emm_df(glm_rich_oxy, "swapped", "Richness"),
  get_emm_df(glm_rich_oxy, "type_original_log", "Richness")
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
  "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//diversity_plot.png",
  final_plot,
  width = 7,
  height = 10,
  dpi = 600
)
#ctd----









