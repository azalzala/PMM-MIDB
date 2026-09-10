library("ape")
BiocManager::install("brms")
library("brms")
library("dplyr")
library("tidyr")
library(purrr)

# for working with processed data
model_data_100_binary <- read.csv('~/shared-team/people/aya/thesis/model_data_100_binary.csv')
# load GTDB tree with taxonomy
tree <- read.tree("~/shared-team/people/aya/thesis/files/PLSDB_GTDB_phylogeny.tree")

# dropping tips using ape

# checking the tips that are different i.e. confirming labels
setdiff(tree$tip.label, model_data_100_binary$GTDB_GENOME)

tips_to_drop <- setdiff(tree$tip.label, model_data_100_binary$GTDB_GENOME)

tree_pruned <- drop.tip(tree, tips_to_drop)

ids <- sort(model_data_100_binary$GTDB_GENOME)
phylo_covar <- vcv.phylo(tree_pruned)
phylo_covar <- phylo_covar / max(phylo_covar)

# troubleshoot 
isSymmetric(phylo_covar) # should be true 
eigen(phylo_covar)$values # must be positive values 
any(eigen(phylo_covar)$values <= 0) # should be FALSE
any(is.na(model_data_binary$representative)) # should be FALSE 

#trying everything with a smaller dataframe >= 100 rep groups

# model 2 - gene specific niche effects
brms_formula_2 <- bf(present ~ niche + (1 + niche | gene_id_all),family = bernoulli(link = "logit"))
fit_integron_2 <-  brm(formula = brms_formula_2, data = model_data_100_binary, chains = 4,cores = 4,iter = 10000,warmup = 2500, backend="cmdstanr")
fit_integron_2
summary(fit_integron_2)

saveRDS(fit_integron_2, "fit_integron_2_100_rep_group.rds")
# interpret results 
# subset based on CI intervals, and others, calculate odds ratio

# levels(model_data_100_binary$niche)

# add in population structure controls 
brms_formula_2.1 <- bf(present ~ niche + (1 + niche | gene_id_all) + (1 | gr(GTDB_GENOME, cov = A)) + (1 | primary_cluster_id),family = bernoulli(link = "logit"))
fit_integron_2.1 <- brm(formula = brms_formula_2.1, data = model_data_100_binary, data2 = list(A = phylo_covar), chains = 4, cores = 4, iter = 10000, warmup = 2500, backend = "cmdstanr")

saveRDS(fit_integron_2.1, "fit_integron_2.1_100_rep.rds")
summary(fit_integron_2.1)

