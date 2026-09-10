# tree final 
# issue: some tips showing up in my_data that have not been included in the tip_df_subset

#loading packages 
library(ape)
library(dplyr)
library(stringr)
library(ggtree)
library(tidyr)
library(tibble)
library(ggplot2)



# loading data 
tree <- read.tree("/Users/ayazalzala/Downloads/PLSDB_GTDB_phylogeny.tree")
tip_df <- read.csv('/Users/ayazalzala/Downloads/PLSDB_GTDB_phylogeny.csv')
my_data  <- read.csv( '/Users/ayazalzala/THESIS_FINAL SCRIPTS/representatives_df.csv')
my_data_2 <- read.csv( '/Users/ayazalzala/THESIS_FINAL SCRIPTS/representatives_df_2.csv' )

# subsetting tree metadata to match uncorrected and corrected genomes 

# find unique tips 
tips_to_keep <- unique(c(
  my_data$tree_tip,
  my_data_2$tree_tip
))

# reduce dataframe to keep only unique entries
tip_df_subset <- tip_df %>%
  select(NUCCORE_ACC, GTDB_GENOME, GTDB_TAXONOMY) %>%
  filter(GTDB_GENOME %in% tips_to_keep) %>%
  mutate(in_set = TRUE)

write.csv(tip_df_subset, "/Users/ayazalzala/Downloads/tips.csv")

# adding tree label "Genus species" 
tip_df_subset <- tip_df_subset %>% mutate(tree_label = str_remove(as.character(GTDB_TAXONOMY),"^.*s__" ))

tip_df_subset <- tip_df_subset %>% mutate(family_name = str_remove(as.character(GTDB_TAXONOMy), "^.*g"))

# adding column to tree metadata subset of formula 1 and formula 2 tips that describes presence across these two dataframes
tip_df_subset <- tip_df_subset %>%
  mutate(
    in_2_0 = GTDB_GENOME %in% my_data$tree_tip, 
    in_2_1 = GTDB_GENOME %in% my_data_2$tree_tip,
    hit  = case_when(
      in_2_0 & !in_2_1  ~ "Uncorrected",
      !in_2_0 & in_2_1  ~ "Corrected",
      in_2_0 & in_2_1   ~ "Both"
    )
  )


# using tree metadata to build final tree 
keep_tips <- intersect(tree$tip.label, tip_df_subset$GTDB_GENOME)
tree_sub <- drop.tip(tree, setdiff(tree$tip.label, keep_tips))


# order for genomes that matches the tree - useful for matchin heatmap and bar charts later
tip_order <- fortify(tree_sub) %>%
  filter(isTip) %>%
  arrange(desc(y)) %>%
  pull(label)

tree_levels <- tibble(GTDB_GENOME = tip_order, tree_pos = seq_along(tip_order))


# relocating tree_tip to beginning of the dataframe
tip_df_subset <- tip_df_subset %>%
  relocate(
    GTDB_GENOME,
    hit,
    in_2_0,
    in_2_1,
    NUCCORE_ACC,
    GTDB_TAXONOMY
  )

head(tip_df_subset)
tree_sub$tip.label
tree$tip.label
names(tip_df_subset)

tip_df_for_plot <- tip_df_subset %>%
  mutate(tree_tip = GTDB_GENOME)

# sizing and spacing 
subtree_p <- ggtree(tree_sub) %<+% tip_df_for_plot +
  geom_tippoint(aes(colour = hit), size = 1) +
  geom_tiplab(
    aes(label = tree_label),
    align = TRUE,
    linetype = "dotted",
    linesize = 0.25,
    offset = 0.01,
    size = 2,
    colour = "black"
  ) +
  xlim_tree(xlim = c(0, 1.5)) +
  scale_colour_manual(
    values = c(
      "None"        = "lightblue4",
      "Uncorrected" = "#A01A9C",
      "Corrected"   = "#3FBC73",
      "Both"        = "orange2"
    ),
    name = "Model hit"
  ) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(10, 0, 0, 0)
  ) +
  geom_treescale(
    width = 0.1,
    offset = 0.2,
    x = 0,
    y = 50
  )




