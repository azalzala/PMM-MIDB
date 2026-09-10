making the dataframe
# first thing was to get nuccore_acc to match with the nuccore_UID - to do this we had to expand the initial biosample dataframe NUCCORE_UID column 
# first thing was to reduce as much as possible the number of rows to annotate 
# annotating on the ECOSYSTEM_QUERY column - has to be present 
# then merged with multiple data tables including the cluster (mob) and plasmidfinder (nucc_acc) 
# with nuccore_acc present you could begin to merge with other dataframes 
# filter out miscallaneous annotations from the integronfinder dataset to then annotate the nuccore_ACC  
# merge annotated dataframe with tree dataframe to know the phylogeny distribution of the plasmids in the dataframe 
#filter out unknowns 
# filter out n >= 100 

library(ggplot2)
library("tidyverse")
library("dplyr")
library("readr")
library("zoo")
library(readr)

# first biosamples should be merged 
biosample <- read.csv('~/shared-team/people/aya/thesis/files/plsdb/biosample.csv')
nuccore <- read.csv('~/shared-team/people/aya/thesis/files/nuccore.csv')
mob <- read.csv('~/shared-team/people/aya/thesis/files/plsdb/typing.csv')
integrons <- read_delim("~/shared-team/people/aya/thesis/files/integronfinder_results_integrons.tsv", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
View(integrons)
bakta_annotations <- read_delim("~/shared-team/people/aya/thesis/files/bakta_annotations.tsv", delim = "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
mmseqs2_clusters <- read_delim("~/shared-team/people/aya/thesis/files/mmseqs2_clusters.tsv", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>% rename("representative" = X1, "ID" = X2)


# need to explode on NUCCORE_ACC for exact primary cluster ID relationship 

biosample_exploded <- biosample %>%
  mutate(NUCCORE_UID = str_remove_all(NUCCORE_UID, "\\[|\\]")) %>%
  separate_rows(NUCCORE_UID, sep = ",\\s*") %>%
  mutate(NUCCORE_UID = str_trim(NUCCORE_UID))


biosample_exploded <- biosample_exploded %>% mutate(NUCCORE_UID = as.numeric(NUCCORE_UID))

nuccore_2 <- nuccore %>% select(NUCCORE_UID, NUCCORE_ACC)
nuccore_2 <- nuccore_2 %>%
  mutate(NUCCORE_UID = as.numeric(NUCCORE_UID))

merge_biosample_nuccore <-  biosample_exploded %>%
  left_join(nuccore_2, by = "NUCCORE_UID")


dataset_C_filtered <- mob %>%
  group_by(primary_cluster_id) %>%
  filter(n() >= 5) %>%
  ungroup()


merge_biosamples_clusters <- left_join(merge_biosample_nuccore, dataset_C_filtered, by="NUCCORE_ACC")

# filtering integrons dataframe to keep only important annotations for later mapping 

integrons_complete <- integrons %>% filter(type=="complete") %>%
  filter(type_elt %in% c("attC", "attI") | (type_elt == "protein" & annotation == "protein")) %>%
  mutate(midpoint = (pos_beg + pos_end) / 2) %>%
  group_by(ID_replicon, ID_integron) %>%
  arrange(midpoint, .by_group = TRUE) %>%
  mutate(
    is_boundary = type_elt %in% c("attC", "attI"),
    last_boundary = na.locf(ifelse(is_boundary, type_elt, NA), na.rm = FALSE),
    next_boundary = na.locf(ifelse(is_boundary, type_elt, NA), fromLast = TRUE, na.rm = FALSE)
  ) %>%
  ungroup() %>%
  filter(
    type_elt == "protein",
    annotation == "protein",
    ((last_boundary == "attC" & next_boundary %in% c("attC", "attI")) |
       (next_boundary == "attC" & last_boundary %in% c("attC", "attI")))
  ) %>%
  mutate(ID = paste(ID_replicon, ID_integron, element, sep="|"))

# error with merging by.x = "ID" - solved by importing and renaming

cds_to_keep <- merge(bakta_annotations, mmseqs2_clusters, by.x="ID") %>%
  group_by(representative) %>%
  filter(!any(Gene == "sul1", na.rm = TRUE)) %>%
  filter(!any(Gene == "qacE", na.rm = TRUE)) %>%
  filter(!any(Gene == "qacEdelta1", na.rm = TRUE)) %>%
  filter(!any(Product == "QACEdelta1", na.rm = TRUE)) %>%
  ungroup() %>%
  pull(ID)

accessions_to_keep <- integrons_complete %>%
  filter(ID %in% cds_to_keep)
dim(integrons_complete)

dim(accessions_to_keep)

# lose the representative column from mmseqs2_clusters so need to merge on ID column following cleanup 
accessions_to_keep <- left_join(accessions_to_keep, mmseqs2_clusters, by = "ID")

# assign niches 
merge_biosample_integrons <- inner_join(merge_biosamples_clusters, accessions_to_keep,by = c("NUCCORE_ACC" = "ID_replicon"))
dim(merge_biosample_integrons)

merge_safe_integrons <- left_join(accessions_to_keep, merge_biosamples_clusters, by = c("ID_replicon" = "NUCCORE_ACC"))


# try ECOSYSTEM_query and ECOSYSTEM_tags for unknowns whilst keeping existing classifications
niche <- merge_biosample_integrons  %>%
  mutate(niche = sapply(ECOSYSTEM_query, assign_niche))

niche <- niche %>%
  mutate(
    niche = ifelse(
      niche == "unknown",
      sapply(ECOSYSTEM_tags, assign_niche),  
      niche                              
    )
  )

# check the counts 

# clusters 
clusters_df <- known_niche_clusters_tip_df %>%
  count(representative, sort = TRUE, name = "Frequency") %>%
  filter(Frequency > 10) %>%
  arrange(desc(Frequency)) %>%
  mutate(representative = factor(representative, levels = representative))

# plots clusters frequency (>10) by their representative group (mixed unknown and known niches) 
ggplot(clusters_df, aes(x = representative, y = Frequency)) +
  geom_col() +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 50, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 100, color = "green", linetype = "dashed") +
  annotate("text", x = Inf, y = 10, label = "n = 10", color = "red", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 50, label = "n = 50", color = "blue", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 100, label = "n = 100", color = "green", hjust = 1.1, vjust = -0.5, size = 3) +
  scale_y_log10() +
  labs(x = "Representative", y = "Count (log scale)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
# after plotting,  I want a list of the same for each 
replicon_counts_unknown <- niche %>% filter(niche == "" | niche == "unknown") %>% count(primary_cluster_id, sort = TRUE)
replicon_counts_known <- niche %>% filter(niche != "unknown") %>% count(primary_cluster_id, sort = TRUE)

replicon_counts_known
barplot(
  replicon_counts_unknown$n,
  names.arg = replicon_counts_unknown$primary_cluster_id,
  las = 2,          # rotate x labels vertically
  col = "steelblue",
  ylab = "Count",
  xlab = "",
  main = "MASH Cluster Frequency of Unknown Niche Plasmids"
)

barplot(
  replicon_counts_known$n,
  names.arg = replicon_counts_known$primary_cluster_id,
  las = 2,          # rotate x labels vertically
  col = "steelblue",
  ylab = "Count",
  xlab = "",
  main = "MASH Cluster Frequency of Known Niche Plasmids"
)

# after plotting the top 10 clusters in known / unknown niches, i want to know for future reference, the number of entries for each niche 
known_niche_df %>% count(niche, sort = TRUE)

# remove unknowns - niche, and primary_cluster_id 
known_niche_df <- niche %>% filter(niche != "unknown")
known_niche_df <- known_niche_df %>% filter(niche != "")

known_niche_clusters_df <- known_niche_df %>% filter(!is.na(primary_cluster_id))
colnames(known_niche_clusters_df)

# merge with GTDB tree  on NUCCORE_ACC 
tip_df <- read.csv('~/shared-team/people/aya/thesis/files/PLSDB_GTDB_phylogeny.csv')
colnames(known_niche_clusters_df)
known_niche_clusters_tip_df <- inner_join(known_niche_clusters_df, tip_df, by = "NUCCORE_ACC")

# remove missing tips (if inner_join didnt work)
known_niche_clusters_tip_df <- known_niche_clusters_tip_df %>% filter(!is.na(GTDB_GENOME))


# Quality Control for Including Rep. Genes in the Dataframe: Determining frequency of rep. gene i.e. gene family for >10, >50, >100 occurences 
known_niche_clusters_tip_df %>% count(representative) %>% filter (n >= 10) # 78 rows
known_niche_clusters_tip_df %>% count(representative) %>% filter (n >= 50) # 31 rows
known_niche_clusters_tip_df %>% count(representative) %>% filter (n >= 100) # 23 rows 
colnames(known_niche_clusters_tip_df)

known_niche_clusters_tip_df_gene <- left_join(known_niche_clusters_tip_df, bakta_annotations, by = "ID")
subset_gene_clusters <- known_niche_clusters_tip_df_gene %>% select(ID, representative, niche, Gene)


# Graph plot to show rep gene frequency and abundance across dataframe
rep_qc_plot <- known_niche_clusters_tip_df_gene %>%
  group_by(representative) %>%
  summarise(
    Frequency = n(),
    Gene = first(Gene),
    .groups = "drop"
  ) %>%
  filter(Frequency > 10) %>%
  arrange(desc(Frequency)) %>%
  mutate(representative = factor(representative, levels = representative))

# Plot representative x frequency with gene labels 
ggplot(rep_qc_plot, aes(x = representative, y = Frequency)) + xlab("Gene") + 
  geom_col() +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 50, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 100, color = "green", linetype = "dashed") +
  annotate("text", x = Inf, y = 10, label = "n = 10", color = "red", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 50, label = "n = 50", color = "blue", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 100, label = "n = 100", color = "green", hjust = 1.1, vjust = -0.5, size = 3) +
  scale_y_log10() +
  scale_x_discrete(labels = rep_qc_plot$Gene) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# Plot gene x frequency 
ggplot(rep_qc_plot, aes(x = Gene, y = Frequency)) +
  geom_col() +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 50, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 100, color = "green", linetype = "dashed") +
  annotate("text", x = Inf, y = 10, label = "n = 10", color = "red", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 50, label = "n = 50", color = "blue", hjust = 1.1, vjust = -0.5, size = 3) +
  annotate("text", x = Inf, y = 100, label = "n = 100", color = "green", hjust = 1.1, vjust = -0.5, size = 3) +
  scale_y_log10() +
  labs(x = "Gene", y = "Count (log scale)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )
# Quick look at the distribution of genes across plasmids from different sampling niches 
known_niche_clusters_tip_df %>% group_by(representative) %>% count(niche)
my_data <- known_niche_clusters_tip_df %>% group_by(representative) %>% filter(n() >= 10) %>% count(niche)
niche_distribution <- model_data_100 %>% count(niche)

# reduce dataframe to select number of rep genes 
known_rep_counts_100_df <- known_niche_clusters_tip_df %>%  group_by(representative) %>% filter(n() >= 100) %>% ungroup()
# forgot a step: the hardest ONE ! I need to make the rep column , binary, somehow! 
model_data_100 <- known_rep_counts_100_df %>% select(representative, niche, primary_cluster_id, GTDB_GENOME)

model_data_100$gene_id <- sub(".*\\|", "", model_data_100$representative)

# increasing size of dataframe 23X Fold BECAUSE each entry will now be multiplied for all unique values of representative across the entire dataset. this is to make the response variable a binary value i.e. gene present 1, gene absent 0

model_data_100_binary <- tidyr::crossing(model_data_100, gene_id_all = unique(model_data_100$gene_id)) %>%
  mutate(present = as.integer(gene_id == gene_id_all))

model_data_100_binary$niche <- factor(model_data_100_binary$niche)
levels(model_data_100_binary$niche)
model_data_100_binary$niche <- relevel(factor(model_data_100_binary$niche), ref = "human")

write.csv(model_data_100_binary, "model_data_binary_100.csv")
