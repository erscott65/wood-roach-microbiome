# wood-roach-microbiome

READ ME

Introduction to the data
This document describes analyses and data associated with the article “Microbiome communities with different transmission modes respond differently to host diet changes in a wood-feeding insect”.

Three different microbial groups from the same set of wood roach guts were analyzed for this study: Parabasalia, Oxymonadida, and bacteria. The R files Amplicon_25_analysis_para.R, Amplicon_25_analysis_oxy.R, and Amplicon_25_analysis_emp.R respectively contain analyses associated with each microbial group. These three microbial groups were analyzed separately using identical methods, and chunks of code correspond between the R files associated with each microbial group.

In the article, we also analyze covariance between ecological features of each pair of microbial groups. The R file Amplicon_25_analysis_all. R contains analyses relating microbial groups to each other.

Raw meta data describing each wood roach sample in this study is in habitat_swap_meta_data.xlsx. Log_site_key.xslx links each log in the study to the site at which it was found.

For each microbial group, we created a phyloseq object (ps_relative_emp_tree, ps_relative_oxy_tree, ps_relative_para_tree). These objects contain four components: metadata for each roach sample, which can be accessed using sample_data(ps); an ASV table describing the relative abundance of each unique amplicon sequence variant in the dataset, which can be accessed using otu_table(ps); a taxonomy table describing the taxonomic identification of each unique amplicon sequence variant in the dataset, which can be accessed using tax_table(ps); and a phylogeny describing the evolutionary relationships between each unique amplicon sequence variant in the dataset, as predicted by a GTR model, which can be accessed using phy_tree(ps). These phyloseq objects can be used to recreate the data analysis associated with this study.

For each microbial group, we also created supplementary taxonomy tables (tax_df_oxy, tax_df_para, and tax_df_emp). These contain the taxonomic information from tax_table(ps) as well as unique tip labels containing numeric identifiers and genus names for each unique amplicon sequence variant.



Description of columns in habitat_swap_meta_data.xlsx
From left to right in this table, we first note the log that each roach was found on (Log). We also note the “family” the roach belongs to (Family) – roaches found in the same gallery were given the same family designation and presumed to be related. The Family_ID combines the log and family designations, and the Roach_ID is a unique identifier given to each roach in the experiment with the first number and letter originating from the Family_ID. We note the sex (if the roach is instar 4 or adult) and the Instar, or life stage of each roach. Life stages are noted as instar 2, 3, 4 or A for adult. We note the weight (g) of the roach at the time it was collected (Intake_weight), and we describe any unique features of the roach in the Notes column. We note the date the roach was collected (Collected) the type of wood it was found on (type_original_log), and whether it was swapped on to a new diet during the experiment (swapped).  We note which log the roach’s new diet was taken from (log_new), which is the smae as its original log if swapped = no). We also note the type of wood the roach was swapped onto (type_log_new), which is the same as type_original_log if swapped = no. We note the mass of the roach at the time the experiment began (starting_mass), as well as the mass of the woody substrate the roach was provided (wood starting_mass) and the date the experiment started (swap_date). Columns P- R, and U-Z note the mass (g) of each roach at different time points during the experiment. Not every roach has a mass noted at every time point, as the experiment start dates for different groups of roaches were staggered. We also not the date roaches molted (Molted) as well as their new lifestage designation post-molt (New Instar). We note the date any roaches were found dead in the experiment (found_dead) and the date roaches were sacrificed and dissected (Dissected). Roaches dissected on the same day had experiment start dates within a week of each other, and were considered part of the same block in the experiment. The remaining columns pertain to storage of the gut and roach samples. 
The meta data contained within the phyloseq objects ps_relative_para_tree, ps_relative_oxy_tree, and ps_relative_emp_tree are subsets of habitat_swap_meta_data.xlsx which contain data only for samples retained throughout the entire experiment, sequencing, and analysis pipeline for each microbial group. These meta data also contain the site data from log_site_key.


Code chunks within Amplicon_25_analysis_para.R, Amplicon_25_analysis_oxy.R, and Amplicon_25_analysis_emp.R

#install packages
In this section we install packages. This section also contains code to save the environment to an Rdata file and load the environment from that Rdata file. *Note: this code saves the environment to a specific filepath, which you will need to modify if replicating the analyses. There is also code in this section to save individual data objects as .rds files, and code to clear the environment.

#tracking reads through workflow
In this section we generate a function to track the number of reads in the dataset at different points in the analysis pipeline

#filepath to raw unzipped fwd and reverse fastas
In this section we initialize a filepath that will be used in future code chunks to read groups of fasta files in and out of the environment. Raw, infiltered demultiplexed fastq files for EMP, oxy, and para amplicon sequences can be accessed on the NCBI SRA  under the accession PRJNA1438676.

#sort sequences into fwd and reverse reads
The original data came back from sequencing with forwards and reverse reads all together in one folder along with reads from another analysis which had been sequencing on the same sequencing run. In this section we create lists of forward and reverse amplicon reads by identifying patterns in file names. 

#remove primers from samples
In this section we remove primers from each amplicon sequence.

#inspect quality plots
In this section we creates plots to visualize average quality scores at each position in the sequence. This allows us to decide how to trim and filter sequences in subsequent steps in the pipeline.

#filter samples
Here we trims reads at the first instance of a quality score below a specified threshold, and discards reads with ambiguous bases, more than the maximum expected number of errors, and any matches to the phix genome.

#merge sequences
Here we merge forward and reverse reads. Then, we merge unique sequence variants so that we have counts of each unique sequence in the dataset (seqtabAll) After that, we remove chimeric sequences (seqtabNoC). In this section also we create a table tracking the number of reads retained at each step in the pipeline (track_oxydf, track_paradf, track_empdf), which we continually add to in subsequent steps.

#assign taxonomy
Here we assign taxonomy to each sequence in the dataset. For Parabasalia and Oxymonadida, the PR2 database was used for taxonomy assignments. For bacteria, the silva database was used for taxonomy assignments.

#construct phyloseq object
Here we read in sample data (habitat_swap_meta_data.xlsx). We then merge meta data describing wood roach samples and treatments, taxonomic information, and sequence data to create a phyloseq object. This object consists of 3 components at this stage: sample data, a tax table, and an out table. We use this data object for all downstream analyses. 

#filter phyloseq object and wrangle sample data
Here we filter the data to remove contamination and low abundance features from the dataset. First, we create a table (track_oxy_asv, track_para_asv, track_emp_asv) to track the number of unique amplicon sequence variants (ASV’s) at each stage of the filtering process. We then remove sequences shorter than 200 bp, sequences associated with samples that were placed in the wrong treatment group, sequences associated with all but one sample for samples that stayed in a family group during the experiment (as we expected family groups would have similar microbial responses to treatments, which could have biases our data), samples with low read counts (<1000 reads), ASV;s represented by fewer than 10 total reads throughout the dataset, and potentially contaminating sequences based on taxonomic identification and blast data. We then transformed ASV counts for filtered data into relative abundance, creating  new filtered, relative-abundance based phyloseq objects (ps_relative_oxy, ps_relative_para, ps_relative_emp). Also in this section we generate meta data sheets which contain the same information as the sample data in the phyloseq objects (meta_oxy, meta_para, meta_emp) which we use for downstream analyses.

#make phylogeny
Here we create a phylogenetic tree describing predicted evolutionary relationships between all sequences retained after filtering. First, we align all sequences, then we compute pairwise distances between them using a maximum likelihood method and use that information to create an initial neighbor-joined tree. We then build two trees using two GTR and JC69 models. We compare AIC values for the GTR and JC69 trees, and we retain the tree with the lower AIC (it was the GTR tree for all 3 microbial groups). We root the tree using midpoint rooting. Finally, we integrate the midpoint-rooted GTR tree with taxonomy, ASV, and metadata to create new phyloseq objects: ps_relative_emp_tree, ps_relative_oxy_tree, ps_relative_para_tree.
In this section we also generate tax tables which contain custom tip labels with numeric identifiers and genus names for each ASV (tax_df_oxy, tax_df_para, tax_df_emp). These are used for visualizing the phylogeny. 
	 


#PERMANOVA (bray curtis)
Here we use PERMANOVA to compare microbial community structure across treatment groups. First, we generate a distance matrix containing bray-curtis similarities for each pair of samples in the dataset. Then we use PERMANOVA to test whether average bray-curtis distances between treatment groups were greater than average distances among treatment groups. We include covariates (Dissected, site/log, and Instar) in the model. We used the by = “terms” parameter to ensure that terms are tested in sequential order (variance explained by terms earlier in the mode is removed before calculating the contribution of later terms). We constructed the model this way to ensure that variance explained by non-treatment covariates be accounted for when evaluating the effects of our treatments. We then tested whether dispersion (i.e., average distance from centroid of an ordination of samples) differed across treatment groups. Finally, we created a 3-dimensional NMDS ordination and plotted that ordination with colors indicating treatment groups. A 2-dimensional ordination was insufficient to visualize our data as it had a stress above 0.2.

#PERMANOVA (jaccard)
In this section we repeat the same steps as in PERMANOVA (bray curtis), except we generate the original distance matrix using jaccard distances.

#PERMANOVA (unifrac unweighted)
In this section we repeat the same steps as in PERMANOVA (bray curtis), except we generate the original distance matrix using unweighted unifrac distances.


#PERMANOVA (unifrac weighted)
In this section we repeat the same steps as in PERMANOVA (bray curtis), except we generate the original distance matrix using weighted unifrac distances.

#multiple tests corrections
Here we adjusted p-values from the PERMANOVA models across community structure metrics using a Bonferroni correction. We used this correction to account for testing the same hypothesis across multiple community structure metrics (bray-curtis, jaccard, weighted unifrac, unweighted unifrac). For metrics where the interaction between diet swapping and original diet was significant after correction, pairwise comparisons among all interaction levels were conducted to determine which levels of diet swapping x original diet different significantly from each other. We generate a new column in the meta data for each microbial group called “treatment” which contains unique levels of diet swap x original diet interaction and describes the diet regime experimented by each roach during the experiment. We conducted these pairwise comparisons by subsetting the distance matrix and metadata to the two treatment levels being compared while retaining the same covariates in the model (Dissected, site/log, and Instar). We used that data subset to perform PERMANOVA tests comparing specific treatment levels. P-values from pairwise comparisons were adjusted for multiple testing within each distance metric using Bonferroni correction.

#diversity analysis
Here we calculate the richness, Shannon, and Simspon’s microbial diversity for each sample. We then use GLM’s to test for the effect of diet treatments on these alpha diversity metrics. In this section we also create boxplots for visualizing differences in raw diversity across treatment groups, as well as dot and whisker plots for visualizing different in emmeans for diversity across treatment groups after GLM’s have been used to account for covariates.

Code chunks within Amplicon_25_analysis_all.R
At the beginning of the file we import RDS objects: phyloseq objects for each microbial group (ps_relative_oxy_tree, ps_relative_emp_tree, ps_relative_para_tree) containing taxonomy, ASV, meta data, and phylogenies for each microbial group. We also import tax_df_emp, tax_df_para, and tax_df_oxy, which contain taxonomic information and custom tree labels for sequences in each microbial group.

#spearman ASV correlation test
Here we use a spearman test to identify covariance in individual ASV abundances between microbial groups. We create a summary table showing the percentages of taxa which covary between each pair of microbial groups. The table also shows the percentage of total samples in each microbial group shared with the paired group, since due to filtering not all samples in the original experiment were retained throughout the analysis pipeline for each microbial group.

#processing raw meta data
Here we process the meta data by editing some data typos, removing instar 2 nymphs, and renaming/ rearranging certain columns. We also add a site column to the metadata.

#table to summarize sample counts by treatment and dataset
Here we count the number of samples in each treatment for each microbial group

#procrustes covariance analysis
Here we use the Procrustes package to overlay ordinations of samples from each pair of microbial groups based on each community structure distance metric (bray-curtis, jaccard, weighted unifrac, and unweighted unifrac). We then rotate and scale the ordinations to get the best possible alignment between them, then test whether each microbial group pair significantly aligns with each other for each distance metric. 

#does community alignment differ across diet
Here we test whether Procrustes alignment differs across levels of treatment. We do this  by evaluating whether distances between paired samples in the aligned ordinations differ across treatment groups.

#plot top 10 most abundant
Here we generate barplots showing the top 10 most abundant ASV’s in each dite treatment for each microbial group.
