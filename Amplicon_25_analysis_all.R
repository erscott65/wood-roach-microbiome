#import phyloseq objects
ps_relative_emp_tree <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//ps_relative_emp_tree.rds")
ps_relative_oxy_tree <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//ps_relative_oxy_tree.rds")
ps_relative_para_tree <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//ps_relative_para_tree.rds")


#import tax info
tax_df_emp <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//tax_df_emp.rds")
tax_df_oxy <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//oxy//tax_df_oxy.rds")
tax_df_para <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//tax_df_para.rds")
#View(tax_df_emp) Sequence Label
#View(tax_df_oxy) ASV_sequence Tree_Label
#View(tax_df_para) ASV_sequence Custom_Label

save.image(file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Amplicon_25_analysis.RData")
load("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Amplicon_25_analysis.RData")

file.info("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Amplicon_25_analysis.RData")$size

readBin("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//Amplicon_25_analysis.RData", what = "raw", n = 20)

#processing raw metadata ----
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

####table to summarize sample counts by treatment and dataset----

library(phyloseq)
library(dplyr)
library(tidyr)

# List of phyloseq objects
phyloseq_list <- list(
  emp = ps_relative_emp_tree,
  oxy = ps_relative_oxy_tree,
  para = ps_relative_para_tree
)

# Function to count samples per treatment
count_treatment_samples <- function(ps_obj, dataset_name) {
  sample_data(ps_obj) %>%
    data.frame() %>%
    group_by(treatment) %>%
    summarise(Sample_Count = n(), .groups = "drop") %>%
    mutate(Dataset = dataset_name)
}

# Combine counts for all datasets
treatment_summary <- bind_rows(
  count_treatment_samples(phyloseq_list$emp, "emp"),
  count_treatment_samples(phyloseq_list$oxy, "oxy"),
  count_treatment_samples(phyloseq_list$para, "para")
)

# Pivot the table so treatments are rows and datasets are columns
treatment_pivot <- treatment_summary %>%
  pivot_wider(names_from = Dataset, values_from = Sample_Count, values_fill = 0)

# View result
print(treatment_pivot)

library(stringr)

treatment_pivot <- treatment_pivot %>%
  mutate(treatment = str_replace_all(treatment, "hardwood", "hard"),
         treatment = str_replace_all(treatment, "softwood", "soft"))


# Load the package
library(writexl)

# Export to Excel
write_xlsx(treatment_pivot, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//treatment_summary.xlsx")


View(sample_data(ps_relative_oxy_tree))
####ctd----

####spearman ASV correlation test----

#ok! the goal here is to determine if there's any covariance between para, oxy, and emp asvs
#let's get into it.

#You have three phyloseq objects (one per amplicon), with overlapping but not identical sample sets, and you want to:
  
  #Identify ASVs whose abundance patterns co-vary across datasets (e.g., high abundance in emp when another ASV is abundant in para).

#Perform this pairwise across all three dataset combinations.

#Statistically assess whether these co-abundances are meaningful.



library(phyloseq)
library(dplyr)

# Test function for ps_relative_oxy_tree vs ps_relative_emp_tree
correlate_phyloseq_verbose <- function(ps1, ps2, 
                                       method = "spearman", 
                                       fdr_cutoff = 0.05, 
                                       cor_cutoff = 0.6, 
                                       label1 = "ASV_1", 
                                       label2 = "ASV_2") {
  # 1. Get shared samples
  shared_samples <- intersect(sample_names(ps1), sample_names(ps2))
  ps1_shared <- prune_samples(shared_samples, ps1)
  ps2_shared <- prune_samples(shared_samples, ps2)
  
  # 2. Extract OTU tables
  otu1 <- as(otu_table(ps1_shared), "matrix")
  otu2 <- as(otu_table(ps2_shared), "matrix")
  
  # 3. Transpose if taxa are rows
  if (taxa_are_rows(ps1)) otu1 <- t(otu1)
  if (taxa_are_rows(ps2)) otu2 <- t(otu2)
  
  # 4. Match sample order
  otu1 <- otu1[shared_samples, ]
  otu2 <- otu2[shared_samples, ]
  
  # 5. Initialize results list
  results <- list()
  total_1 <- length(colnames(otu1))
  total_2 <- length(colnames(otu2))
  
  message("Starting correlation analysis:")
  message("  ", label1, " ASVs: ", total_1)
  message("  ", label2, " ASVs: ", total_2)
  message("  Total comparisons: ", total_1 * total_2)
  
  # 6. Loop through ASVs and compute correlations
  count <- 0
  for (asv1 in colnames(otu1)) {
    for (asv2 in colnames(otu2)) {
      test <- suppressWarnings(cor.test(otu1[, asv1], otu2[, asv2], method = method))
      results[[paste0(asv1, "_vs_", asv2)]] <- data.frame(
        ASV_1 = asv1,
        ASV_2 = asv2,
        correlation = test$estimate,
        p_value = test$p.value
      )
      count <- count + 1
    }
    message("Finished ", asv1, " (", count, "/", total_1 * total_2, " tests completed)")
  }
  
  # 7. Combine and adjust p-values
  cor_df <- bind_rows(results)
  cor_df$p_adj <- p.adjust(cor_df$p_value, method = "fdr")
  
  # 8. Filter significant correlations
  cor_sig <- cor_df %>%
    filter(p_adj < fdr_cutoff, abs(correlation) > cor_cutoff) %>%
    arrange(desc(abs(correlation)))
  
  message("✅ Correlation analysis complete.")
  message("  Significant correlations (|r| > ", cor_cutoff, ", FDR < ", fdr_cutoff, "): ", nrow(cor_sig))
  
  return(list(
    all_correlations = cor_df,
    significant = cor_sig
  ))
}

??cor.test
# cor.test# Run the test correlation function for oxy-emp comparisons and process results
test_results <- correlate_verbose(ps_relative_oxy_tree, ps_relative_emp_tree)

# View significant results
head(test_results$significant)

nrow(test_results$significant)
View(test_results$significant)
# Count unique oxy ASVs
length(unique(test_results$significant$ASV_oxy))

# Count unique emp ASVs
length(unique(test_results$significant$ASV_emp))
#ok so basically there are a lot of correlations, way more correlations than significantly differential taxa. interesting

View(tax_df_emp)

annotated_results <- test_results$significant %>%
  # Join Tree_Label from tax_df_oxy, matching ASV_oxy to ASV_sequence
  left_join(
    tax_df_oxy[, c("ASV_sequence", "Tree_Label")],
    by = c("ASV_oxy" = "ASV_sequence")
  ) %>%
  # Join Custom_Label from tax_df_emp, matching ASV_emp to Seqeunce
  left_join(
    tax_df_emp[, c("Sequence", "Label")],
    by = c("ASV_emp" = "Sequence")
  )
colnames(tax_df_emp)

View(annotated_results)
annotated_results <- annotated_results %>%
  rename(
    Oxy_Label = Tree_Label,
    Emp_Label = Label
  )

annotated_results_oxy_emp <- annotated_results
View(annotated_results_oxy_emp)







# Run the test correlation function for para- emp comparisons and process results
test_results <- correlate_phyloseq_verbose(ps_relative_para_tree, ps_relative_emp_tree)

# View significant results
head(test_results$significant)

nrow(test_results$significant)
View(test_results$significant)

#rename to match dataset names
test_results$significant <- test_results$significant %>%
  rename(
    ASV_para = ASV_1,
    ASV_emp = ASV_2
  )

View(test_results$significant)

# Count unique oxy ASVs
length(unique(test_results$significant$ASV_para))

# Count unique emp ASVs
length(unique(test_results$significant$ASV_emp))
#ok so basically there are a lot of correlations, way more correlations than significantly differential taxa. interesting

#annotate with custom lables containing asv ID and tax info
View(tax_df_emp)
View(tax_df_para)

annotated_results_para_emp <- test_results$significant %>%
  # Join Tree_Label from tax_df_oxy, matching ASV_oxy to ASV_sequence
  left_join(
    tax_df_para[, c("ASV_Sequence", "Custom_Label")],
    by = c("ASV_para" = "ASV_Sequence")
  ) %>%
  # Join Custom_Label from tax_df_emp, matching ASV_emp to Seqeunce
  left_join(
    tax_df_emp[, c("Sequence", "Label")],
    by = c("ASV_emp" = "Sequence")
  )
View(annotated_results_para_emp)


annotated_results_para_emp <- annotated_results_para_emp %>%
  rename(
    Para_Label = Custom_Label,
    Emp_Label = Label
  )

View(annotated_results_para_emp)







# Run the test correlation function for para- oxy comparisons and process results
test_results <- correlate_phyloseq_verbose(ps_relative_para_tree, ps_relative_oxy_tree)

# View significant results
head(test_results$significant)

nrow(test_results$significant)
View(test_results$significant)

#rename to match dataset names
test_results$significant <- test_results$significant %>%
  rename(
    ASV_para = ASV_1,
    ASV_oxy = ASV_2
  )

View(test_results$significant)

# Count unique oxy ASVs
length(unique(test_results$significant$ASV_para))

# Count unique emp ASVs
length(unique(test_results$significant$ASV_oxy))
#ok so basically there are a lot of correlations, way more correlations than significantly differential taxa. interesting

#annotate with custom lables containing asv ID and tax info
View(tax_df_para)
View(tax_df_oxy)

annotated_results_para_oxy <- test_results$significant %>%
  # Join Tree_Label from tax_df_oxy, matching ASV_oxy to ASV_sequence
  left_join(
    tax_df_para[, c("ASV_Sequence", "Custom_Label")],
    by = c("ASV_para" = "ASV_Sequence")
  ) %>%
  # Join Custom_Label from tax_df_emp, matching ASV_emp to Seqeunce
  left_join(
    tax_df_oxy[, c("ASV_sequence", "Tree_Label")],
    by = c("ASV_oxy" = "ASV_sequence")
  )
View(annotated_results_para_oxy)


annotated_results_para_oxy <- annotated_results_para_oxy %>%
  rename(
    Para_Label = Custom_Label,
    Oxy_Label = Tree_Label
  )

View(annotated_results_para_oxy)






#generate summary table
# Sample counts
n_samples_para <- nsamples(ps_relative_para_tree)
n_samples_oxy  <- nsamples(ps_relative_oxy_tree)
n_samples_emp  <- nsamples(ps_relative_emp_tree)

# Taxa counts
n_taxa_para <- ntaxa(ps_relative_para_tree)
n_taxa_oxy  <- ntaxa(ps_relative_oxy_tree)
n_taxa_emp  <- ntaxa(ps_relative_emp_tree)

# Overlapping samples
n_samples_para_emp <- length(intersect(sample_names(ps_relative_para_tree), sample_names(ps_relative_emp_tree)))
n_samples_para_oxy <- length(intersect(sample_names(ps_relative_para_tree), sample_names(ps_relative_oxy_tree)))
n_samples_oxy_emp  <- length(intersect(sample_names(ps_relative_oxy_tree),  sample_names(ps_relative_emp_tree)))

# Unique significant ASVs per comparison
n_sig_para_pe <- length(unique(annotated_results_para_emp$Para_Label))
n_sig_emp_pe  <- length(unique(annotated_results_para_emp$Emp_Label))

n_sig_para_po <- length(unique(annotated_results_para_oxy$Para_Label))
n_sig_oxy_po  <- length(unique(annotated_results_para_oxy$Oxy_Label))

n_sig_oxy_oe  <- length(unique(annotated_results_oxy_emp$Oxy_Label))
n_sig_emp_oe  <- length(unique(annotated_results_oxy_emp$Emp_Label))


summary_table <- data.frame(
  Comparison = c("para–emp", "para–oxy", "oxy–emp"),
  
  Common_Samples = c(n_samples_para_emp, n_samples_para_oxy, n_samples_oxy_emp),
  
  Total_Para_Samples = c(n_samples_para, n_samples_para, NA),
  Total_Para_Taxa    = c(n_taxa_para, n_taxa_para, NA),
  
  Total_Oxy_Samples  = c(NA, n_samples_oxy, n_samples_oxy),
  Total_Oxy_Taxa     = c(NA, n_taxa_oxy, n_taxa_oxy),
  
  Total_Emp_Samples  = c(n_samples_emp, NA, n_samples_emp),
  Total_Emp_Taxa     = c(n_taxa_emp, NA, n_taxa_emp),
  
  Covarying_Para_Taxa = c(n_sig_para_pe, n_sig_para_po, NA),
  Covarying_Oxy_Taxa  = c(NA, n_sig_oxy_po, n_sig_oxy_oe),
  Covarying_Emp_Taxa  = c(n_sig_emp_pe, NA, n_sig_emp_oe)
)

View(summary_table)

#reorganized column names
summary_table <- summary_table[, c("Comparison",
                                   "Covarying_Emp_Taxa",
                                   "Covarying_Para_Taxa", 
                                   "Covarying_Oxy_Taxa",
                                   "Total_Emp_Taxa",
                                   "Total_Para_Taxa",
                                   "Total_Oxy_Taxa",
                                   "Common_Samples",
                                   "Total_Emp_Samples",
                                   "Total_Para_Samples",
                                   "Total_Oxy_Samples")]
View(summary_table)



print(summary_table)



#configure percentage summary table
library(dplyr)
library(tidyr)

# Create a tidy data frame with one row per dataset per comparison
percentage_table <- summary_table %>%
  transmute(
    Comparison,
    emp_covarying     = Covarying_Emp_Taxa,
    emp_total_taxa    = Total_Emp_Taxa,
    emp_total_samples = Total_Emp_Samples,
    
    para_covarying     = Covarying_Para_Taxa,
    para_total_taxa    = Total_Para_Taxa,
    para_total_samples = Total_Para_Samples,
    
    oxy_covarying     = Covarying_Oxy_Taxa,
    oxy_total_taxa    = Total_Oxy_Taxa,
    oxy_total_samples = Total_Oxy_Samples,
    
    Common_Samples
  ) %>%
  pivot_longer(
    cols = -c(Comparison, Common_Samples),
    names_to = c("Dataset", ".value"),
    names_pattern = "(.+?)_(.+)"
  ) %>%
  # Keep only rows where Total Taxa is NOT NA (dataset actually in comparison)
  filter(!is.na(total_taxa)) %>%
  mutate(
    `% Covarying Taxa` = 100 * covarying / total_taxa,
    `% Common Samples` = 100 * Common_Samples / pmin(total_samples, na.rm = TRUE)
  ) %>%
  rename(
    `Covarying Taxa` = covarying,
    `Total Taxa` = total_taxa,
    `Total Samples` = total_samples,
    `Common Samples` = Common_Samples
  ) %>%
  select(Comparison, Dataset, `Covarying Taxa`, `Total Taxa`,
         `% Covarying Taxa`, `Total Samples`, `Common Samples`, `% Common Samples`) %>%
  arrange(Comparison, Dataset)

View(percentage_table)
print(percentage_table)

#round percentages
percentage_table <- percentage_table %>%
  mutate(
    `% Covarying Taxa` = round(`% Covarying Taxa`, 1),
    `% Common Samples` = round(`% Common Samples`, 1)
  )

print(percentage_table)


#create barplot
library(ggplot2)

ggplot(percentage_table, aes(x = Comparison, y = `% Covarying Taxa`, fill = Dataset)) +
  geom_col(position = position_dodge()) +  # side-by-side bars
  geom_text(aes(label = sprintf("%.1f%%", `% Covarying Taxa`)),
            position = position_dodge(width = 0.9), vjust = -0.5, size = 3) + # labels on top
  labs(
    title = "% Covarying Taxa by Dataset and Comparison",
    x = "Comparison",
    y = "% Covarying Taxa",
    fill = "Dataset"
  ) +
  theme_minimal() +
  ylim(0, max(percentage_table$`% Covarying Taxa`) * 1.15)


#add percentage signs to values in percentage column

percentage_table <- percentage_table %>%
  mutate(
    `% Covarying Taxa` = paste0(round(`% Covarying Taxa`, 1), "%"),
    `% Common Samples` = paste0(round(`% Common Samples`, 1), "%")
  )

View(percentage_table)




# Create numeric version of the percent column
percentage_table <- percentage_table %>%
  mutate(percent_numeric = as.numeric(sub("%", "", `% Covarying Taxa`)))

# Plot
ggplot(percentage_table,
       aes(x = Comparison, y = percent_numeric, fill = Dataset)) +
  geom_col(position = position_dodge()) +  # side-by-side bars
  geom_text(aes(label = `% Covarying Taxa`), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  labs(
    title = "% Covarying Taxa by Dataset and Comparison",
    x = "Comparison",
    y = "% Covarying Taxa",
    fill = "Dataset"
  ) +
  theme_minimal() +
  ylim(0, max(percentage_table$percent_numeric) * 1.15)

#export as excel file
install.packages("openxlsx")  # if not already installed
library(openxlsx)

# Create a new Excel workbook
wb <- createWorkbook()

# Add sheets with each data frame
addWorksheet(wb, "para_vs_oxy")
writeData(wb, "para_vs_oxy", annotated_results_para_oxy)

addWorksheet(wb, "para_vs_emp")
writeData(wb, "para_vs_emp", annotated_results_para_emp)

addWorksheet(wb, "oxy_vs_emp")
writeData(wb, "oxy_vs_emp", annotated_results_oxy_emp)

addWorksheet(wb, "summary_table")
writeData(wb, "summary_table", summary_table)

# Add new worksheet for percentage_table
addWorksheet(wb, "percentage_table")
writeData(wb, "percentage_table", percentage_table)

# Save the workbook, overwriting the old file
saveWorkbook(wb, file = "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//asv_correlation_summary.xlsx", overwrite = TRUE)






#check overlap between correlating taxa and DA taxa
ps_relative_para_tree <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//para//ps_relative_para_tree.rds")
sig_treatment_taxa_passed_emp <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//sig_treatment_taxa_passed.rds")
sig_Instar_taxa_passed_emp <- readRDS("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//emp//sig_Instar_taxa_passed.rds")

View(sig_Instar_taxa_passed_emp)

View(annotated_results_oxy_emp)
View(annotated_results_para_emp)


# make a list of significant taxa in the oxy-emp analysis that match significant Instar taxa emp
oxy_emp_labels_matched <- annotated_results_oxy_emp$Emp_Label[
  annotated_results_oxy_emp$ASV_emp %in% sig_Instar_taxa_passed_emp$taxon
]
print(oxy_emp_labels_matched)

# make a list of significant taxa in the para-emp analysis that match significant Instar taxa emp
para_emp_labels_matched <- annotated_results_para_emp$Emp_Label[
  annotated_results_para_emp$ASV_emp %in% sig_Instar_taxa_passed_emp$taxon
]

print(para_emp_labels_matched)

####ctd----

####procrustes covariance analysis----

#first, create a list of dataset pairs and then subset the datasets in each pair to contain only common samples and associated taxa
library(phyloseq)

# Your phyloseq objects
datasets <- list(
  emp = ps_relative_emp_tree,
  oxy = ps_relative_oxy_tree,
  para = ps_relative_para_tree
)

View(sample_data(ps_relative_emp_tree))
View(sample_data(ps_relative_oxy_tree))
View(sample_data(ps_relative_para_tree))


# Generate all unique pairs of dataset names
pairs <- combn(names(datasets), 2, simplify = FALSE)
print(pairs)

# Function to subset to common samples AND remove zero-count taxa
subset_common_samples <- function(ps1, ps2) {
  common_samples <- intersect(sample_names(ps1), sample_names(ps2))
  
  # Print before subsetting
  cat("Before subsetting:\n")
  cat("ps1:", nsamples(ps1), "samples,", ntaxa(ps1), "taxa\n")
  cat("ps2:", nsamples(ps2), "samples,", ntaxa(ps2), "taxa\n")
  cat("Common samples:", length(common_samples), "\n")
  
  # Subset to common samples and remove unobserved taxa
  ps1_sub <- prune_samples(common_samples, ps1)
  ps1_sub <- prune_taxa(taxa_sums(ps1_sub) > 0, ps1_sub)
  
  ps2_sub <- prune_samples(common_samples, ps2)
  ps2_sub <- prune_taxa(taxa_sums(ps2_sub) > 0, ps2_sub)
  
  # Print after subsetting
  cat("After subsetting:\n")
  cat("ps1:", nsamples(ps1_sub), "samples,", ntaxa(ps1_sub), "taxa\n")
  cat("ps2:", nsamples(ps2_sub), "samples,", ntaxa(ps2_sub), "taxa\n")
  cat("-----\n")
  
  return(list(ps1 = ps1_sub, ps2 = ps2_sub))
}

# Create a list to store results for each pair
paired_datasets <- list()

# Loop over pairs and subset to common samples
for (pair in pairs) {
  ps1_name <- pair[1]
  ps2_name <- pair[2]
  
  cat("Processing pair:", ps1_name, "vs", ps2_name, "\n")
  
  ps1 <- datasets[[ps1_name]]
  ps2 <- datasets[[ps2_name]]
  
  subsetted <- subset_common_samples(ps1, ps2)
  
  paired_datasets[[paste(ps1_name, ps2_name, sep = "_vs_")]] <- subsetted
}


library(phyloseq)
library(vegan)        # for distances like Jaccard, Bray
library(ape) 

#compute distance matrices
# Create a list to store distance matrices
distances <- list()

for (pair_name in names(paired_datasets)) {
  cat("Computing distances for pair:", pair_name, "\n")
  
  ps1 <- paired_datasets[[pair_name]]$ps1
  ps2 <- paired_datasets[[pair_name]]$ps2
  
  # Distance functions
  dist_list <- function(ps) {
    list(
      bray      = phyloseq::distance(ps, method = "bray"),
      jaccard   = phyloseq::distance(ps, method = "jaccard"),
      wunifrac  = phyloseq::distance(ps, method = "wunifrac"),
      unifrac   = phyloseq::distance(ps, method = "unifrac")
    )
  }
  
  distances[[pair_name]] <- list(
    ps1 = dist_list(ps1),
    ps2 = dist_list(ps2)
  )
}


distance_types <- c("bray", "jaccard", "unifrac", "wunifrac")

#inspect distance matrices
# Loop through all pairs
for (pair_name in names(distances)) {
  cat("\n========== Pair:", pair_name, "==========\n")
  
  for (ps_label in c("ps1", "ps2")) {
    cat("\n-- Dataset:", ps_label, "--\n")
    
    for (dist_type in distance_types) {
      cat("\nDistance type:", dist_type, "\n")
      
      # Try to access the matrix and convert to full matrix form
      dist_obj <- distances[[pair_name]][[ps_label]][[dist_type]]
      
      # Check if it exists (in case something failed earlier)
      if (!is.null(dist_obj)) {
        dist_mat <- as.matrix(dist_obj)
        print(dist_mat[1:min(6, nrow(dist_mat)), 1:min(6, ncol(dist_mat))])
      } else {
        cat("Distance matrix not found.\n")
      }
    }
  }
}



#plot overlaid distance matrix pairs
library(vegan)
library(ggplot2)
library(gridExtra)

# Run NMDS on each distance matrix in your distances list
nmds_results <- list()

for (pair_name in names(distances)) {
  nmds_results[[pair_name]] <- list()
  
  for (ps_label in c("ps1", "ps2")) {
    nmds_results[[pair_name]][[ps_label]] <- list()
    
    for (dist_type in distance_types) {
      dist_mat <- distances[[pair_name]][[ps_label]][[dist_type]]
      
      if (!is.null(dist_mat)) {
        # Run NMDS; try k=2 dimensions and a reasonable number of tries
        nmds <- metaMDS(dist_mat, k=3, trymax=100, autotransform=FALSE)
        nmds_results[[pair_name]][[ps_label]][[dist_type]] <- nmds
      }
    }
  }
}



#plot raw NMDS overlays
save_raw_nmds_plots_3d <- function(nmds_results, output_dir = ".", method_filter = NULL) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  for (pair_name in names(nmds_results)) {
    for (method in names(nmds_results[[pair_name]]$ps1)) {
      
      if (!is.null(method_filter) && method != method_filter) next
      
      nmds1 <- nmds_results[[pair_name]]$ps1[[method]]
      nmds2 <- nmds_results[[pair_name]]$ps2[[method]]
      
      if (is.null(nmds1) || is.null(nmds2)) next
      
      # Extract 3D coordinates directly
      coords1 <- as.data.frame(scores(nmds1, display = "sites"))
      coords2 <- as.data.frame(scores(nmds2, display = "sites"))
      
      # Ensure they both have 3 columns
      if (ncol(coords1) != 3 || ncol(coords2) != 3) {
        warning("Skipping ", pair_name, " - ", method, " (NMDS not 3D)")
        next
      }
      
      colnames(coords1) <- paste0("NMDS", 1:3)
      colnames(coords2) <- paste0("NMDS", 1:3)
      coords1$Sample <- rownames(coords1)
      coords2$Sample <- rownames(coords2)
      
      # Merge on Sample
      joined <- merge(coords1, coords2, by = "Sample", suffixes = c("_1", "_2"))
      
      # Extract dataset names
      pair_split <- strsplit(pair_name, "_vs_")[[1]]
      dataset1 <- pair_split[1]
      dataset2 <- pair_split[2]
      
      # Create base plot
      p <- plot_ly(
        data = joined,
        x = ~NMDS1_1, y = ~NMDS2_1, z = ~NMDS3_1,
        type = "scatter3d",
        mode = "markers",
        marker = list(size = 5, color = "#1f77b4"),
        name = dataset1,
        text = ~Sample,
        hoverinfo = "text"
      ) %>%
        add_trace(
          x = ~NMDS1_2, y = ~NMDS2_2, z = ~NMDS3_2,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 5, color = "#ff7f0e"),
          name = dataset2,
          text = ~Sample,
          hoverinfo = "text"
        )
      
      # Add lines between matched samples
      lines_x <- unlist(mapply(function(x1, x2) c(x1, x2, NA), joined$NMDS1_1, joined$NMDS1_2))
      lines_y <- unlist(mapply(function(y1, y2) c(y1, y2, NA), joined$NMDS2_1, joined$NMDS2_2))
      lines_z <- unlist(mapply(function(z1, z2) c(z1, z2, NA), joined$NMDS3_1, joined$NMDS3_2))
      
      p <- add_trace(
        p,
        x = lines_x, y = lines_y, z = lines_z,
        type = "scatter3d",
        mode = "lines",
        line = list(color = 'gray', width = 1),
        showlegend = FALSE,
        hoverinfo = "none"
      )
      
      # Final layout
      p <- layout(
        p,
        title = list(text = paste("Raw NMDS:", pair_name, "-", method), x = 0.5),
        scene = list(
          xaxis = list(title = "NMDS1"),
          yaxis = list(title = "NMDS2"),
          zaxis = list(title = "NMDS3"),
          aspectmode = "data"
        ),
        showlegend = TRUE
      )
      
      # Save to file
      filename <- file.path(
        output_dir,
        paste0(gsub("[^a-zA-Z0-9_]", "_", paste(pair_name, method, "rawNMDS", sep = "_")), ".html")
      )
      
      saveWidget(p, filename, selfcontained = TRUE)
      message("✅ Saved raw NMDS plot: ", filename)
    }
  }
}

save_raw_nmds_plots_3d(
  nmds_results,
  output_dir = "C:/Users/erinv/Box/Cryptocercus_research/Amplicon_Data_25/raw_nmds_plots",
  method_filter = NULL
)

names(nmds_results[["emp_vs_oxy"]]$ps1)

#procrustes analysis
# Store results here
distance_methods <- c("bray", "jaccard", "unifrac", "wunifrac")
procrustes_results <- list()
summary_table <- list()

for (pair_name in names(distances)) {
  for (method in distance_methods) {
    cat("\n🔍 Running Procrustes for pair:", pair_name, "| method:", method, "\n")
    
    dist1 <- distances[[pair_name]]$ps1[[method]]
    dist2 <- distances[[pair_name]]$ps2[[method]]
    
    # Run NMDS
    nmds1 <- try(metaMDS(dist1, k = 3, trymax = 100), silent = TRUE)
    nmds2 <- try(metaMDS(dist2, k = 3, trymax = 100), silent = TRUE)
    
    if (inherits(nmds1, "try-error") || inherits(nmds2, "try-error")) {
      cat("❌ NMDS failed for", pair_name, "method:", method, "\n")
      summary_table[[paste(pair_name, method, sep = "_")]] <- list(
        pair = pair_name,
        method = method,
        R = NA,
        p_value = NA,
        mean_disp1 = NA,
        mean_disp2 = NA,
        dispersion_ratio = NA,
        note = "NMDS failed"
      )
      next
    }
    
    # Procrustes alignment
    proc <- procrustes(nmds1, nmds2, symmetric = TRUE)
    ??procrustes
    
    # Run protest
    protest_result <- try(protest(nmds1, nmds2, permutations = 999), silent = TRUE)
    
    if (inherits(protest_result, "try-error")) {
      cat("❌ protest() failed for", pair_name, "method:", method, "\n")
      summary_table[[paste(pair_name, method, sep = "_")]] <- list(
        pair = pair_name,
        method = method,
        R = NA,
        p_value = NA,
        mean_disp1 = NA,
        mean_disp2 = NA,
        dispersion_ratio = NA,
        note = "protest() failed"
      )
      next
    }
    
    # Store successful results
    if (is.null(procrustes_results[[pair_name]])) {
      procrustes_results[[pair_name]] <- list()
    }
    
    procrustes_results[[pair_name]][[method]] <- list(
      procrustes = proc,
      protest = protest_result,
      nmds1 = nmds1,
      nmds2 = nmds2
    )
    
    # Calculate dispersion (mean distance to centroid) in NMDS space
    coords1 <- scores(nmds1)
    coords2 <- scores(nmds2)
    
    dist_coords1 <- vegdist(coords1)
    dist_coords2 <- vegdist(coords2)
    
    disp1 <- betadisper(dist_coords1, rep("group1", nrow(coords1)))
    disp2 <- betadisper(dist_coords2, rep("group2", nrow(coords2)))
    
    mean_disp1 <- mean(disp1$distances)
    mean_disp2 <- mean(disp2$distances)
    dispersion_ratio <- mean_disp1 / mean_disp2
    
    # Warn if strong imbalance
    note <- ifelse(dispersion_ratio > 2 | dispersion_ratio < 0.5,
                   "⚠️ Dispersion imbalance",
                   "")
    
    # Save summary
    R <- sqrt(1 - protest_result$ss)
    pval <- protest_result$signif
    
    cat("✅ R =", round(R, 3), "| p =", signif(pval, 3), 
        "| Disp1 =", round(mean_disp1, 3), 
        "| Disp2 =", round(mean_disp2, 3),
        "| Ratio =", round(dispersion_ratio, 2), note, "\n")
    
    summary_table[[paste(pair_name, method, sep = "_")]] <- list(
      pair = pair_name,
      method = method,
      R = R,
      p_value = pval,
      mean_disp1 = mean_disp1,
      mean_disp2 = mean_disp2,
      dispersion_ratio = dispersion_ratio,
      note = note
    )
  }
}


# Convert summary to data frame
summary_df <- do.call(rbind, lapply(summary_table, as.data.frame))

# Round for readability
summary_df$R <- signif(summary_df$R, 4)
summary_df$p_value <- signif(summary_df$p_value, 4)
summary_df$mean_disp1 <- round(summary_df$mean_disp1, 4)
summary_df$mean_disp2 <- round(summary_df$mean_disp2, 4)
summary_df$dispersion_ratio <- round(summary_df$dispersion_ratio, 3)

# Add asterisks to p-values
summary_df$p_value <- ifelse(as.numeric(summary_df$p_value) <= 0.001, 
                             paste0(summary_df$p_value, "***"),
                             ifelse(as.numeric(summary_df$p_value) <= 0.01, 
                                    paste0(summary_df$p_value, "**"),
                                    ifelse(as.numeric(summary_df$p_value) <= 0.05, 
                                           paste0(summary_df$p_value, "*"),
                                           as.character(summary_df$p_value))))

# Remove row names
rownames(summary_df) <- NULL

# View results
cat("\n📊 Procrustes + Dispersion Summary:\n")
print(summary_df)


#Load the package
library(writexl)

# Write to Excel
write_xlsx(summary_df, "C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//procrustes_summary.xlsx")
####ctd----

####does community alignment differ across diet#####
library(nlme)

model_results <- list()

for (pair_name in names(procrustes_results)) {
  for (method in distance_methods) {
    
    result <- procrustes_results[[pair_name]][[method]]
    if (is.null(result)) next
    
    # Extract Procrustes residuals
    proc <- result$procrustes
    resid_vals <- residuals(proc)
    
    # Extract metadata
    ps1 <- paired_datasets[[pair_name]]$ps1
    meta <- data.frame(sample_data(ps1))
    meta$SampleID <- rownames(meta)
    
    # Combine residuals with metadata
    df <- data.frame(
      SampleID = names(resid_vals),
      Residual = as.numeric(resid_vals),
      Method = method,
      Pair = pair_name
    )
    df <- merge(df, meta, by = "SampleID")
    
    # Check required columns
    required_cols <- c("type_original_log", "swapped", "site", "Instar", "Dissected")
    if (!all(required_cols %in% colnames(df))) {
      warning("⚠️ Missing required metadata columns for: ", pair_name)
      next
    }
    
    # Drop rows with missing values
    df <- droplevels(df[complete.cases(df[, required_cols]), ])
    
    # Ensure at least 2 levels per factor
    if (length(unique(df$type_original_log)) < 2 || length(unique(df$swapped)) < 2) next
    
    # Interaction factor for variance modeling
    df$treatment <- interaction(df$type_original_log, df$swapped, drop = TRUE)
    
    # Fit GLS with heteroscedasticity across treatment levels
    gls_fit <- try(
      gls(
        Residual ~ type_original_log * swapped + site + Instar + Dissected,
        data = df,
        weights = varIdent(form = ~1 | treatment),
        method = "REML"
      ),
      silent = TRUE
    )
    
    if (inherits(gls_fit, "try-error")) {
      warning("⚠️ GLS failed for: ", pair_name, " ", method)
      next
    }
    
    # Extract ANOVA table for fixed effects
    anova_table <- anova(gls_fit)
    anova_df <- data.frame(
      Pair = pair_name,
      Method = method,
      Term = rownames(anova_table),
      F_value = signif(anova_table$F, 4),
      P_value = signif(anova_table$p, 4),
      Model = "GLS"
    )
    
    model_results[[paste(pair_name, method, sep = "_")]] <- anova_df
  }
}

# Combine all results
results_df <- do.call(rbind, model_results)

#check for significnat diet regime effects
# Define main hypotheses terms
main_terms <- c("type_original_log", "swapped", "type_original_log:swapped")

# Filter for significant effects
significant_results <- subset(
  results_df,
  Term %in% main_terms & P_value < 0.05
)

# Optional: sort by p-value
significant_results <- significant_results[order(significant_results$P_value), ]

# Display
cat("\n📌 Significant Effects (p < 0.05) for Main Hypotheses:\n")
print(significant_results)
#no significant effects of diet regime!

#check for significant other effects
# Filter for any significant effects (p < 0.05)
all_significant <- subset(
  results_df,
  P_value < 0.05
)

# Optional: sort by p-value
all_significant <- all_significant[order(all_significant$P_value), ]

# Display
cat("\n📌 All Significant Effects (p < 0.05):\n")
print(all_significant)

#export summary table
# Load required package
library(openxlsx)

# 1. Subset relevant columns
summary_df <- results_df[, c("Pair", "Method", "Term", "F_value", "P_value")]

# 2. Relabel Pair names
summary_df$Pair <- gsub("emp",  "bacteria",     summary_df$Pair)
summary_df$Pair <- gsub("oxy",  "oxymonadida",  summary_df$Pair)
summary_df$Pair <- gsub("para", "parabasalia", summary_df$Pair)

# 3. Relabel model terms
summary_df$Term <- gsub("^type_original_log$", "original diet", summary_df$Term)
summary_df$Term <- gsub("^swapped$", "diet swap", summary_df$Term)
summary_df$Term <- gsub("^type_original_log:swapped$",
                        "original diet:diet swap",
                        summary_df$Term)
summary_df$Term <- gsub("^Dissected$", "block", summary_df$Term)

# 4. Optional: clean capitalization
summary_df$Term <- gsub("^Instar$", "instar", summary_df$Term)

# 5. Optional: sort table for readability
summary_df <- summary_df[order(summary_df$Pair,
                               summary_df$Method,
                               summary_df$Term), ]

summary_df$P_value <- paste0(
  formatC(summary_df$P_value, format = "f", digits = 4),
  ifelse(summary_df$P_value < 0.001, "***",
         ifelse(summary_df$P_value < 0.01, "**",
                ifelse(summary_df$P_value < 0.05, "*", "")))
)
View(summary_df)

# Remove intercept term from summary table
summary_df <- subset(summary_df, Term != "(Intercept)")

#"C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//alignment across diet regimes.xlsx"

#summary_df <- read.xlsx("C://Users//erinv//Box//Cryptocercus_research//Amplicon_Data_25//alignment across diet regimes.xlsx", sheet = 1)


# Helper for capitalization rules
format_group <- function(x) {
  if (tolower(x) == "bacteria") {
    "bacteria"
  } else {
    tools::toTitleCase(tolower(x))
  }
}

summary_df$Pair <- vapply(
  strsplit(summary_df$Pair, "_vs_"),
  function(groups) {
    paste(format_group(groups[1]),
          "vs.",
          format_group(groups[2]))
  },
  character(1)
)

head(summary_df)

# Relabel distance methods
summary_df$Method <- gsub("^bray$",     "Bray-Curtis",        summary_df$Method, ignore.case = TRUE)
summary_df$Method <- gsub("^jaccard$",  "Jaccard",            summary_df$Method, ignore.case = TRUE)
summary_df$Method <- gsub("^unifrac$",  "Unweighted UniFrac", summary_df$Method, ignore.case = TRUE)
summary_df$Method <- gsub("^wunifrac$", "Weighted UniFrac",   summary_df$Method, ignore.case = TRUE)



# 6. Export to Excel
output_file <- "alignment across diet regimes.xlsx"

write.xlsx(
  summary_df,
  file = output_file,
  sheetName = "alignment across diet regimes",
  rowNames = FALSE
)


####ctd----

####plot top 10 most abundant----

#I want to generate relative abundance plots showing the top 10 taxa by abundance in each treatment
#for each dataset

library(tidyr)
library(dplyr)

library(phyloseq)
library(tidyverse)


#filter the top 10 most abundant asv's for each treatment and then plot
plot_top_asvs_labeled <- function(ps_rel, tax_df, treatment_var, top_n = 10,
                                  seq_col, label_col, dataset_name = NULL) {
  library(stringr)  # for str_replace_all
  
  # Extract and format OTU table
  otu <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  otu <- as.data.frame(otu)
  
  otu$SampleID <- rownames(otu)
  
  # Extract metadata
  meta <- as.data.frame(sample_data(ps_rel))
  meta$SampleID <- rownames(meta)
  
  # Convert to long format and join metadata
  otu_long <- otu %>%
    pivot_longer(cols = -SampleID, names_to = "ASV", values_to = "Abundance") %>%
    left_join(meta[, c("SampleID", treatment_var)], by = "SampleID")
  
  # Calculate mean abundance per treatment × ASV
  grouped <- otu_long %>%
    group_by(!!sym(treatment_var), ASV) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop")
  
  # Filter to top N ASVs per treatment group
  plot_df <- grouped %>%
    group_by(!!sym(treatment_var)) %>%
    slice_max(mean_abund, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  # Map ASV sequences to short labels
  label_map <- tax_df %>%
    select(seq = !!sym(seq_col), label = !!sym(label_col)) %>%
    distinct()
  
  plot_df <- plot_df %>%
    left_join(label_map, by = c("ASV" = "seq")) %>%
    mutate(ASV_label = ifelse(is.na(label), ASV, label))
  
  # Clean up treatment labels for x-axis by removing "wood"
  plot_df <- plot_df %>%
    mutate(
      treatment_clean = str_replace_all(!!sym(treatment_var), "wood", "")
    )
  
  # Build plot title
  plot_title <- dataset_name
  
  # Plot
  ggplot(plot_df, aes(x = treatment_clean, y = mean_abund, fill = ASV_label)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      x = treatment_var,
      y = "Mean Relative Abundance",
      fill = "ASV",
      title = plot_title
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_top_asvs_labeled(
  ps_relative_emp_tree, tax_df_emp,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "Sequence",
  label_col = "Label",
  dataset_name = "Emp Dataset"
)

plot_top_asvs_labeled(
  ps_relative_oxy_tree, tax_df_oxy,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "ASV_sequence",
  label_col = "Tree_Label",
  dataset_name = "Oxy Dataset"
)

plot_top_asvs_labeled(
  ps_relative_para_tree, tax_df_para,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "ASV_Sequence",
  label_col = "Custom_Label",
  dataset_name = "Para Dataset"
)

library(tidyr)
library(dplyr)
library(ggplot2)

#filter the top 10 most abundant asvs in each dataset and then plot
plot_top_asvs_total <- function(ps_rel, tax_df, treatment_var, top_n = 10,
                                  seq_col, label_col, dataset_name = NULL) {
  library(stringr)  # for str_replace_all
  
  # Extract and format OTU table
  otu <- as(otu_table(ps_rel), "matrix")
  if (taxa_are_rows(ps_rel)) {
    otu <- t(otu)
  }
  otu <- as.data.frame(otu)
  otu$SampleID <- rownames(otu)
  
  # Extract metadata
  meta <- as.data.frame(sample_data(ps_rel))
  meta$SampleID <- rownames(meta)
  
  # Convert to long format and join metadata
  otu_long <- otu %>%
    pivot_longer(cols = -SampleID, names_to = "ASV", values_to = "Abundance") %>%
    left_join(meta[, c("SampleID", treatment_var)], by = "SampleID")
  
  # --- NEW STEP: Identify top N ASVs across the whole dataset ---
  top_asvs <- otu_long %>%
    group_by(ASV) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop") %>%
    slice_max(mean_abund, n = top_n, with_ties = FALSE) %>%
    pull(ASV)
  
  # Filter to only top N ASVs
  otu_top <- otu_long %>%
    filter(ASV %in% top_asvs)
  
  # Calculate mean abundance per treatment × ASV
  plot_df <- otu_top %>%
    group_by(!!sym(treatment_var), ASV) %>%
    summarise(mean_abund = mean(Abundance), .groups = "drop")
  
  # Map ASV sequences to short labels
  label_map <- tax_df %>%
    select(seq = !!sym(seq_col), label = !!sym(label_col)) %>%
    distinct()
  
  plot_df <- plot_df %>%
    left_join(label_map, by = c("ASV" = "seq")) %>%
    mutate(ASV_label = ifelse(is.na(label), ASV, label))
  
  # Clean up treatment labels for x-axis by removing "wood"
  plot_df <- plot_df %>%
    mutate(
      treatment_clean = str_replace_all(!!sym(treatment_var), "wood", "")
    )
  
  # Build plot title
  plot_title <- dataset_name
  
  # Plot
  ggplot(plot_df, aes(x = treatment_clean, y = mean_abund, fill = ASV_label)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(
      x = treatment_var,
      y = "Mean Relative Abundance",
      fill = "ASV",
      title = plot_title
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



plot_top_asvs_total(
  ps_relative_emp_tree, tax_df_emp,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "Sequence",
  label_col = "Label",
  dataset_name = "Emp Dataset"
)

plot_top_asvs_total(
  ps_relative_oxy_tree, tax_df_oxy,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "ASV_sequence",
  label_col = "Tree_Label",
  dataset_name = "Oxy Dataset"
)

plot_top_asvs_total(
  ps_relative_para_tree, tax_df_para,
  treatment_var = "treatment",
  top_n = 10,
  seq_col = "ASV_Sequence",
  label_col = "Custom_Label",
  dataset_name = "Para Dataset"
)
####ctd----


