figure 1  - formula 2

# process results 
library(brms)
library(tidyverse)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)

# model findings - formula 1 

# converting 
gene_effects_2 <- as.data.frame(coef(fit_integron_2.1)$gene_id_all) %>%
  rownames_to_column("gene_id_all")

# outputs are the same as hits, except hits includes effect, OR, and probability (include in final table)
gene_effects_long_3 <- gene_effects_2 %>%
  pivot_longer(
    cols = -gene_id_all,
    names_to = c(".value", "niche"),
    names_pattern = "(Estimate|Est.Error|Q2.5|Q97.5)\\.(.*)"
  ) %>% filter(niche != "Intercept") %>% 
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5,
    se = Est.Error
  ) %>%   mutate(
    niche = str_remove(niche, "^niche")
  ) %>% 
  filter(lower > 0 | upper < 0)

# appendix version (includes non-significant hits)
gene_effects_long_appendix_3 <- gene_effects_2 %>%
  pivot_longer(
    cols = -gene_id_all,
    names_to = c(".value", "niche"),
    names_pattern = "(Estimate|Est.Error|Q2.5|Q97.5)\\.(.*)"
  ) %>% filter(niche != "Intercept") %>%
  rename(
    estimate = Estimate,
    lower = Q2.5,
    upper = Q97.5,
    se = Est.Error
  ) %>%
  mutate(
    niche = str_remove(niche, "^niche")) 

# calculating OR, effects, probability 
draws_2.1 <- as_draws_df(fit_integron_2.1)

genes_2.1 <- names(draws_2.1)[grepl("^r_gene_id_all\\[.*?,Intercept\\]$", names(draws_2.1))] %>%
  str_extract("(?<=\\[).*(?=,Intercept\\])")

niches <- c( 
  "nicheBacteria",
  "nichedomesticatedanimal",
  "nicheenvironmental",
  "nichefood", 
  "nichelaboratory", 
  "nichelivestock",  
  "nicheplant", 
  "nichewastewater",
  "nichewildlife"   
)

results_2.1 <- map_dfr(genes_2.1, function(gene){
  
  map_dfr(niches, function(niche){
    
    int_name <- paste0("r_gene_id_all[", gene, ",Intercept]")
    slope_name <- paste0("r_gene_id_all[", gene, ",", niche, "]")
    
    stopifnot(int_name %in% names(draws_2.1))
    stopifnot(slope_name %in% names(draws_2.1))
    
    effect <- draws_2.1[[paste0("b_", niche)]] +
      draws_2.1[[slope_name]]
    
    log_odds <- draws_2.1$b_Intercept +
      draws_2.1[[int_name]] +
      effect
    
    
    tibble(
      gene_id_all = gene,
      niche = niche,
      effect = effect,
      OR = exp(effect),
      prob = plogis(log_odds)
    )
    
  })
  
})

tibble(
  gene_id_all = gene,
  niche = niche,
  effect = effect,
  OR = exp(effect),
  prob = plogis(log_odds)
)

})

})

gene_effects_table_3 <- results_2.1 %>%
  group_by(gene_id_all, niche) %>%
  summarise(
    effect_median = median(effect),
    effect_lower  = quantile(effect, 0.025),
    effect_upper  = quantile(effect, 0.975),
    
    OR_median = median(OR),
    OR_lower  = quantile(OR, 0.025),
    OR_upper  = quantile(OR, 0.975),
    
    prob_median = median(prob),
    prob_lower  = quantile(prob, 0.025),
    prob_upper  = quantile(prob, 0.975),
    
    .groups = "drop"
  ) %>%
  arrange(desc(effect_median))


# hits- reported niche associations
hits_3 <- gene_effects_table_3 %>%
  filter(effect_lower > 0 | effect_upper < 0) %>% 
  mutate(niche = str_remove(niche, "^niche"))


# Gene 
hits_3 <- left_join(hits_3, bakta_annotations %>% select(Gene, Product, gene_id), by = c("gene_id_all" = "gene_id"))

hits_3 <- hits_3 %>% arrange(desc(effect_median))

# saving tables needed for report 
write.csv(gene_effects_table_3, "~/shared-team/people/aya/thesis/appendix_formula2.csv") # ordered by OR
write.csv(hits_3, "~/shared-team/people/aya/thesis/hits_formula2.csv") # ordered by OR 

# making the dataframe for visualisations

seed_genes_3 <- hits_3 %>%
  pull(gene_id_all) %>% unique()

# need to use first part of the representative "ID_replicon" to find all genes in acccessions to keep and then merge with bakta_annotations and bring those genes
# into the dataset with present = 0

all_seed_genes_3 <- model_data_100_binary %>% filter(gene_id %in% seed_genes_3)
all_seed_genes_3 <- all_seed_genes_3 %>% select(-gene_id_all, -present)
all_seed_genes_unique_3 <- unique(all_seed_genes_3)

mmseqs2_clusters_mod <- mmseqs2_clusters %>%
  transmute(
    representative,
    ID,
    representative_gene_id = sub(
      "^[^|]+\\|[^|]+\\|",
      "",
      representative
    ),
    member_accession = sub(
      "\\|.*",
      "",
      ID
    )
  ) %>%
  distinct()

# i want to make a change so that gene_id.x is populated with the gene_id as it merges the mmseqs and all_seed_genes_unique 

representatives_mmseqs_3 <- mmseqs2_clusters_mod %>%
  left_join(
    all_seed_genes_unique_3 %>%
      select(
        representative,
        gene_id,
        niche,
        primary_cluster_id, 
        GTDB_GENOME
      ) %>%
      distinct(),
    by = "representative",
    suffix = c("", "_model")
  ) %>%
  mutate(
    gene_id = gene_id
  ) %>%
  distinct()

# join with bakta_annotations and tip_df to get the remaining tree tips then add present column where all members of the same representative_gene_id should be present (genes) and those that have been identified as enriched/depleted are the columns /accessions which have a niche assigned to them . 

# checking which columns bakta_annotations$ID is in mmseqs2_clusters_mod$ID or mmseqs2_clusters_mod$representative
all(bakta_annotations$ID %in% mmseqs2_clusters_mod$ID) # TRUE
# same for tip_df and representatives_mmseqs_genes_tips
all(representatives_mmseqs_genes$member_accession %in% tip_df$NUCCORE_ACC) # FALSE

representatives_mmseqs_genes_3 <- representatives_mmseqs_3 %>% left_join(bakta_annotations %>% select(ID, Gene, Product),  by = "ID") 
representatives_mmseqs_genes_tips_3 <- left_join(representatives_mmseqs_genes_3, tip_df %>% select(NUCCORE_ACC, GTDB_GENOME, GTDB_TAXONOMY), by = c("member_accession" = "NUCCORE_ACC"))
representatives_mmseqs_genes_tips_clusters_3 <- left_join(representatives_mmseqs_genes_tips_3, mob %>% select(NUCCORE_ACC, primary_cluster_id), by = c("member_accession" = "NUCCORE_ACC"))
representatives_mmseqs_genes_tips_clusters_final_3 <- representatives_mmseqs_genes_tips_clusters_3 %>% select(representative, ID, representative_gene_id, member_accession, Gene, Product, niche, primary_cluster_id.y, GTDB_GENOME.y, GTDB_TAXONOMY)


# join with hits by gene_id_all + gene !
hits_3 %>%
  count(gene_id_all, Gene) %>%
  filter(n > 1)

representatives_figure_df_3 <- representatives_mmseqs_genes_tips_clusters_final_3 %>%
  left_join(hits_3,by = c(
    "representative_gene_id" = "gene_id_all",
    "Gene"))

representatives_figure_df_3 <- representatives_figure_df_3 %>% rename(
  primary_cluster = primary_cluster_id.y, 
  tree_tip = GTDB_GENOME.y,  
  taxonomy = GTDB_TAXONOMY
)

representatives_figure_df_3 <- representatives_figure_df_3 %>%
  mutate(
    effect_label = case_when(
      effect_median < 0 ~ paste0("Depleted in ", niche),
      effect_median > 0 ~ paste0("Enriched in ", niche),
      effect_median == 0 ~ "No effect",
      TRUE ~ "No effect"
    )
  )

effect_levels <- sort(unique(
  na.omit(heatmap_df$effect_label)
))

effect_levels

effect_colours <- c(
  "Depleted in environmental" = "#66C2A5",
  "Depleted in livestock" = "#8DA0CB",
  "Enriched in livestock" = "#E78AC2",
  "Depleted in wastewater" = "#A6D854",
  "Enriched in wastewater" = "#FFD92F", 
  "No effect" = "#FC8D62"
)

my_data_2 <- representatives_figure_df_3 %>% filter(representative_gene_id %in% seed_genes_3)
my_data_2 <- my_data_2 %>%
  left_join(hits,by = c(
    "representative_gene_id" = "gene_id_all",
    "Gene"))

# effect_label missing
write.csv(my_data_2, "~/shared-team/people/aya/thesis/representatives_df_2.csv")

# create the table with representatives mapped to independent accessions 

# list of gtdb_genome that map to all members of the hit representative_gene_id 

# final dataframe download for building the tree and heatmap 
