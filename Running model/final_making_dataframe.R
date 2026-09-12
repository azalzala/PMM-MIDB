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

# annotate 
# annotate
niche_lookup <- tribble(
  ~keyword,                ~broad_niche,
  # Human
  "homo sapiens",          "human",
  "Homo", "human", 
  "homo", "human", 
  "Homosapiens", "human", 
  "human",                 "human",
  "patient",               "human",
  "clinical",              "human",
  "anthropogenic",         "human", 
  "Children", "human", 
  "blood", "human", 
  "bile", "human", 
  "Hono sapiens", "human", 
  "ward", "human", 
  "Homo-sapiens", "human",
  "9 year old girl suffering from paratyphoid fever", "human",
  "Korean adult feces", "human",
  "Korean infant feces", "human",
  "dental root canal", "human", 
  "Dialysis fluid", "human", 
  "Swab from nose slime of 72 years old man", "human", 
  # Livestock
  "farm",                 "livestock", 
  "sus scrofa",            "livestock",
  "bos taurus",            "livestock",
  "gallus gallus",         "livestock",
  "chicken",               "livestock",
  "bovine",                "livestock",
  "pig",                   "livestock",
  "swine",                 "livestock",
  "cattle",                "livestock",
  "poultry",               "livestock",
  "duck",                 "livestock",
  "Gallus", "livestock", 
  "hen", "livestock", 
  "Gallas gallas", "livestock", 
  "broiler", "livestock", 
  "Livestock", "livestock", 
  "Calf", "livestock", 
  "calf", "livestock", 
  "Capra aegagrus hircus", "livestock", 
  "goat", "livestock", 
  "Goat", "livestock", 
  "Capra hircus", "livestock", 
  "Caprinae (caprine)", "livestock",
  "Sheep", "livestock", 
  "sheep", "livestock", 
  "cow", "livestock", 
  "Cow", "livestock", 
  "cows", "livestock", 
  "lung of a calf", "livestock", 
  "Apis mellifera", "livestock", 
  "Bos primigenius taurus", "livestock", 
  "porcine", "livestock",
  "Triticum turgidum subsp. durum", "livestock", 
  
  # Wildlife
  "Spermophilus beecheyi", "wildlife", 
  "Spodoptera frugiperda", "wildlife", 
  
  "Macaca assamensis", "wildlife",
  "rabbit", "wildlife", 
  "Rabbit", "wildlife", 
  "Bombyx mori", "wildlife", 
  "Caenorhabditis elegans",  "wildlife", 
  "caprine", "wildlife", 
  "Macaca mulatta",  "wildlife", 
  "Mammal", "wildlife", 
  "Manis javanica", "wildlife", 
  "larus",                 "wildlife",
  "ips typographus",       "wildlife",
  "crow", "wildlife", 
  "Crow", "wildlife", 
  "Chroicocephalus novaehollandiae", "wildlife", 
  "thaumetopoea",          "wildlife",
  "rodents",               "wildlife",
  "rodent",                "wildlife",
  "mouse", "wildlife", 
  "Mus musculus", "wildlife", 
  "Siamese crocodiles", "wildlife", 
  "giant panda", "wildlife", 
  "housefly", "wildlife", 
  "Junco", "wildlife", 
  "Giant Panda", "wildlife", 
  "dog",                    "wildlife",
  "insect",                "wildlife",
  "Marmot",                "wildlife",
  "marmot", "wildlife", 
  "Mastacembelus", "wildlife", 
  "Takifugu obscurus", "wildlife", 
  "Meles meles", "wildlife", 
  "Marmota", "wildlife", 
  "Marmota baibacina", "wildlife",
  "Aepyceros melampus", "wildlife",
  "Tenebrio Molitor",      "wildlife",
  "Tenebrio molitor", "wildlife",
  "Wild boar", "wildlife", 
  "Bat", "wildlife", 
  "weever", "wildlife", 
  "Vulpes vulpes", "wildlife", 
  "Camelus dromedarius", "wildlife", 
  "Ceratitis capitata", "wildlife", 
  "Chamaeleonidae", "wildlife", 
  "fleas", "wildlife", 
  "goose", "wildlife",
  "deer", "wildlife", 
  "Aves", "wildlife", 
  "avian", "wildlife",
  "Avian", "wildlife", 
  "Anabas testudineus", "wildlife", 
  "Salmo salar", "wildlife", 
  "honey bee", "wildlife", 
  "Bivalve mollusc", "wildlife", 
  "elephant", "wildlife", 
  "Coturnix coturnix", "wildlife", 
  "Blattella germanica", "wildlife",
  "Horse", "wildlife", 
  "Grus grus", "wildlife", 
  "Delichon urbicum", "wildlife", 
  "bird", "wildlife", 
  "Camel", "wildlife",
  "animal", "wildlife", 
  "Anser", "wildlife",
  "feline", "wildlife", 
  "Frog", "wildlife", 
  "Bos mutus", "wildlife", 
  "Aratinga solstitialis", "wildlife",
  "Branta leucopsis", "wildlife", 
  "Bubalus bubalis", "wildlife", 
  "Brook charr", "wildlife", 
  "Buffalo milk", "wildlife", 
  "Moschus berezovskii", "wildlife", 
  "mushroom", "wildlife", 
  "Nezara viridula" , "wildlife", 
  "bull manure", "wildlife", 
  "Canis lupus", "wildlife", 
  "Canis latrans", "wildlife", 
  "horse", "wildlife", 
  "Gopherus berlandieri", "wildlife", 
  "Deinagkistrodon acutus", "wildlife", 
  "Dendrolimus sibericus", "wildlife", 
  "Gull", "wildlife", 
  "silver gull", "wildlife", 
  "Hawk", "wildlife", 
  "Helicoverpa armigera", "wildlife", 
  "Hippopotamus amphibius", "wildlife", 
  "Larimichthys crocea", "wildlife", 
  "Lates calcarifer", "wildlife", 
  "Oropsylla silantiewi", "wildlife", 
  "Ochotona", "wildlife", 
  "Anguilla marmorata", "wildlife", 
  "Amazona aestiva", "wildlife", 
  "Bactrocera dorsalis", "wildlife", 
  "Bighead Carp", "wildlife", 
  "black bass", "wildlife", 
  "Black kite", "wildlife", 
  "Black-collared Starling", "wildlife",
  "Zophobas atratus", "wildlife", 
  "Cervus elaphus (wild-living)", "wildlife", 
  "Chaeturichthys stigmatias", "wildlife", 
  "Correlophus ciliatus", "wildlife", 
  "Costelytra zealandica", "wildlife", 
  "Culicoides impunctatus", "wildlife", 
  "Eclectus roratus", "wildlife", 
  "equine abscess isolate", "wildlife", 
  "Equus", "wildlife", 
  "Galleria mellonella", "wildlife", 
  "grass carp", "wildlife", 
  "Ictalurus", "wildlife", 
  "Salvelinus fontinalis", "wildlife", 
  "Sturnus vulgaris", "wildlife", 
  "Procyon lotor", "wildlife", 
  "Pseudocaranx dentex", "wildlife", 
  "Pteropus poliocephalus", "wildlife", 
  "pheasant", "wildlife", 
  "Pelodiscus sinensis", "wildlife", 
  
  # Food
  "apple", "food", 
  "agricultural plant", "food", 
  "Cucurbita pepo ssp. texana", "food", 
  "braised burdock", "food", 
  "food",                  "food",
  "carrot", "food", 
  "Morone chrysops x Morone saxatilis", "food", 
  "Diced lamb", "food", 
  "lamb", "food", 
  "meat",                  "food",
  "blue mussel", "food", 
  "dairy",                 "food",
  "eggshell rinse", "food", 
  "produce",               "food",
  "cucumber", "food", 
  "pork",                 "food", 
  "frozen peas", "food", 
  "Frozen yellow eel chunk", "food", 
  "turkey",                "food", 
  "(Sweet Corn)",            "food", 
  "Bean paste", "food", 
  "Triticum aestivum",     "food",
  "Clam", "food", 
  "coriander",  "food",
  "Raw mutton", "food", 
  "blueberry", "food", 
  "cabbage", "food", 
  "red spinach", "food",
  "betel leaf", "food", 
  "Lactuca sativa", "food", 
  "Spinacia oleracea", "food", 
  "White radish", "food", 
  "Oyster", "food", 
  "smoked salmon", "food", 
  "cooked salmon", "food",
  "fish", "food", 
  "Atlantic salmon", "food",
  "beef", "food",
  "brewing yeast sample", "food", 
  "Capsicum annuum", "food", 
  "sprouts", "food", 
  "lettuce", "food", 
  "shrimp", "food", 
  "shrimip", "food", 
  "apple", "food", 
  
  
  # Wastewater
  "drain",                 "wastewater",
  "sewage",                "wastewater",
  "wastewater",            "wastewater",
  "effluent",              "wastewater", 
  "sludge",                "wastewater",
  # Environmental
  "soil",                  "environmental",
  "sediment",              "environmental",
  "river",                 "environmental",
  "marine",                "environmental",
  "water",                 "environmental",
  "environment",           "environmental",
  "environmental",         "environmental",
  "150 km offshore",       "environmental", 
  "aquatic",               "environmental",
  "forest", "environmental", 
  "terrestrial", "environmental", 
  "Mud", "environmental",
  "mud", "environmental", 
  "faeces", "environmental", 
  "fecal", "environmental",
  "Fecal", "environmental", 
  "feces", "environmental",
  "faecal", "environmental",
  "Faecal", "environmental", 
  "Dust", "environmental", 
  "stool", "environmental",
  "koi carp bred in the Czech Republic", "environmental", 
  "litter", "environmental", 
  "Litter", "environmental", 
  
  
  # Pets 
  "Canis lupus familiaris", "domesticated animal",
  "Cat", "domesticated animal",
  "cat", "domesticated animal", 
  "Musca domestica", "domesticated animal",
  "Columba livia", "domesticated animal",
  "canine", "domesticated animal", 
  "Canine", "domesticated animal",
  "Equus caballus", "domesticated animal", 
  "Equus ferus caballus", "domesticated animal", 
  "Red-breasted Parakeet", "domesticated animal", 
  
  
  # Plants
  "mulberry", "plant",
  "Brassica napus", "plant", 
  "Plant", "plant",
  "Alhagi sparsifolia Shap.", "plant",
  "almond drupe", "plant", 
  "Bursaphelenchus xylophilus", "plant", 
  "Cornus officinalis", "plant",
  "Dracaena sanderiana", "plant",
  "Eucalyptus grandis", "plant", 
  "Lolium perenne", "plant", 
  "Rice seed", "plant", 
  'Solanum tuberosum', 'plant',
  "Cotinus coggygria", "plant", 
  "Crotalaria pallida", "plant", 
  "Euphorbia granulata", "plant",
  "Glycine max", "plant", 
  "Flower", "plant", 
  "Stachytarpheta glabra", "plant", 
  "Suaeda salsa", "plant", 
  "Phaseolus vulgaris", "plant", 
  
  
  # Cell culture 
  #"cell_culture", "cell culture",
  #"cell culture", "cell culture"
  
  "laboratory strain", "laboratory",
  "Laboratory strain", "laboratory",
  "Laboratory isolate", "laboratory", 
  "lab strain", "laboratory", 
  "Transcojugate labstrain", "laboratory", 
  "transconjugant", "laboratory", 
  "Clone isolated from the evolution experiment", "laboratory", 
  "Cell culture", "laboratory",
  "Derived from parental isolate EC0880B in laboratory.", "laboratory", 
  "Derived from parental isolate EC0026B in laboratory.", "laboratory", 
  "derived from strain s4454", "laboratory", 
  "Drosophila melanogaster", "laboratory", 
  "DSMZ Isolate", "laboratory", 
  "EQA test strain", "laboratory", 
  "NIST Mixed Microbial RM strain", "laboratory", 
  "Seth Lab strain", "laboratory", 
  # Bacteria
  "bacteria", "Bacteria", 
  "Bacteria", "Bacteria",
  "bateria", "Bacteria",
  
  "not isolated", "",
  "not available: to be reported later", ""
)

# final dataset : biosample_UID +nuccore_acc, broad_niche, host lineage, plasmid cluster, integron data
assign_niche <- function(source) {
  source <- str_trim(source)
  match <- niche_lookup %>% filter(str_detect(source, keyword))
  if (nrow(match) > 0) return(match$broad_niche[1])
  return("unknown")
}

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
