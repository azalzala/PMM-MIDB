library(brms)
library(tidyverse)
library(tibble)
library(dplyr)
library(stringr)
library(tidyr)

# model findings - formula 1 
fit_integron_2 <- readRDS()
# converting 
gene_effects <- as.data.frame(coef(fit_integron_2)$gene_id_all) %>%
  rownames_to_column("gene_id_all")

# outputs are the same as hits, except hits includes effect, OR, and probability (include in final table)
gene_effects_long_2 <- gene_effects %>%
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
gene_effects_long_appendix <- gene_effects %>%
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
draws <- as_draws_df(fit_integron_2)

genes <- names(draws)[grepl("^r_gene_id_all\\[.*?,Intercept\\]$", names(draws))] %>%
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

results <- map_dfr(genes, function(gene){
  
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

gene_effects_table <- results %>%
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
hits <- gene_niche_table %>%
  filter(effect_lower > 0 | effect_upper < 0)

# Gene 
bakta_annotations$gene_id <- sub(".*\\|", "", bakta_annotations$ID)
hits <- inner_join(hits, bakta_annotations, by = c("gene_id_all" = "gene_id"))
hits <- hits %>% select(gene_id_all, Gene, niche, effect_median, effect_lower, effect_upper, OR_median, OR_lower, OR_upper, prob_median, prob_lower, prob_upper)

hits <- hits %>% arrange(desc(effect_median))
hits$Gene[5] <- "HAD family hydrolase"

# saving tables needed for report 
write.csv(gene_niche_table, "~/shared-team/people/aya/thesis/appendix_formula1.csv") # ordered by OR
write.csv(hits, "~/shared-team/people/aya/thesis/hits_formula1.csv") # ordered by OR 


# making the dataframe for visualisations

seed_genes <- hits %>%
  pull(gene_id_all) %>% unique()

# need to use first part of the representative "ID_replicon" to find all genes in acccessions to keep and then merge with bakta_annotations and bring those genes
# into the dataset with present = 0

all_seed_genes <- model_data_100_binary %>% filter(gene_id %in% seed_genes)
all_seed_genes <- all_seed_genes %>% select(-gene_id_all, -present)
all_seed_genes_unique <- unique(all_seed_genes)

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

representatives_mmseqs <- mmseqs2_clusters_mod %>%
  left_join(
    all_seed_genes_unique %>%
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

representatives_mmseqs_genes <- representatives_mmseqs %>% left_join(bakta_annotations %>% select(ID, Gene, Product),  by = "ID") 
representatives_mmseqs_genes_tips <- left_join(representatives_mmseqs_genes, tip_df %>% select(NUCCORE_ACC, GTDB_GENOME, GTDB_TAXONOMY), by = c("member_accession" = "NUCCORE_ACC"))
representatives_mmseqs_genes_tips_clusters <- left_join(representatives_mmseqs_genes_tips, mob %>% select(NUCCORE_ACC, primary_cluster_id), by = c("member_accession" = "NUCCORE_ACC"))
representatives_mmseqs_genes_tips_clusters_final <- representatives_mmseqs_genes_tips_clusters %>% select(representative, ID, representative_gene_id, member_accession, Gene, Product, niche, primary_cluster_id.y, GTDB_GENOME.y, GTDB_TAXONOMY)


# join with hits by gene_id_all + gene !
hits %>%
  count(gene_id_all, Gene) %>%
  filter(n > 1)

hits <- hits %>%   mutate(niche = str_remove(niche, "^niche"))

representatives_figure_df <- representatives_mmseqs_genes_tips_clusters_final %>%
  left_join(hits,by = c(
    "representative_gene_id" = "gene_id_all",
    "Gene",
    "niche"))

representatives_figure_df <- representatives_figure_df %>% rename(
  primary_cluster = primary_cluster_id.y, 
  tree_tip = GTDB_GENOME.y,  
  taxonomy = GTDB_TAXONOMY
)

representatives_figure_df <- representatives_figure_df %>%
  mutate(
    effect_label = case_when(
      effect_median < 0 ~ paste0("Depleted in ", niche),
      effect_median > 0 ~ paste0("Enriched in ", niche),
      effect_median == 0 ~ "No effect",
      TRUE ~ "No effect"
    )
  )

effect_levels <- sort(unique(
  na.omit(representatives_figure_df$effect_label)
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


write.csv(representatives_figure_df, "~/shared-team/people/aya/thesis/representatives_df.csv")