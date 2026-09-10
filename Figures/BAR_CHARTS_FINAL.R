# loading packages 
library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(ggplot2)
library(viridis)
library(patchwork)

# loading data -- tree_final.r should be run before this for ordering the genomes along y axis 
tip_df <- read.csv('/Users/ayazalzala/Downloads/PLSDB_GTDB_phylogeny.csv')
my_data  <- read.csv('/Users/ayazalzala/THESIS_FINAL SCRIPTS/representatives_df.csv')
my_data_2 <- read.csv('/Users/ayazalzala/Downloads/representatives_figure_df_2.csv')

# BAR CHARTS 

# theme for bar charts 
theme_bar_viridis <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(color = "black"),
      axis.title.x = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# switch between my_data and my_data_2 
niche_df <- my_data %>%
  select(member_accession, tree_tip, niche.x) %>%
  filter(!is.na(tree_tip)) %>%
  distinct()

tip_df_subset <- tip_df_subset %>%
  distinct(GTDB_TAXONOMY, .keep_all = TRUE)

niche_df <- niche_df %>% left_join(tip_df_subset %>% select(NUCCORE_ACC, GTDB_GENOME, tree_label), by = c("tree_tip" = "GTDB_GENOME"), relationship = "many-to-one")
niche_counts <- niche_df %>% distinct(tree_tip, niche.x) %>% count(tree_tip, niche.x)
niche_counts <- niche_counts %>%
  mutate(
    tree_tip = factor(
      tree_tip,
      levels = rev(tip_order),
      ordered = TRUE
    )
  )

bar1 <- ggplot(niche_counts, aes(x = tree_tip, y = n, fill = niche.x)) +
  geom_col(width= 0.95) + coord_flip() +
  labs(x = "Genome", y = "Niche count", x = NULL, fill = NULL) +
  theme_bar_viridis()  + theme(axis.text.y = element_text(size = 8), legend.position = "none") + scale_fill_brewer(palette = "Set1") + scale_y_continuous(expand = c(0, 0))

bar1
# representative IDs barchart 


representative_df <- my_data %>%
  select(tree_tip, Gene) %>%
  filter(!is.na(tree_tip), !is.na(Gene)) %>%
  distinct() %>%
  left_join(tree_levels, by = c("tree_tip" = "GTDB_GENOME"))

representative_counts <- representative_df %>% count(tree_tip, Gene) 
representative_counts <- representative_counts %>% 
  mutate(
    tree_tip = factor(
      tree_tip,
      levels = rev(tip_order),   # reverse so top tip = leftmost bar
      ordered = TRUE
    )
  )

bar2 <- ggplot(
  representative_counts,
  aes(x = tree_tip, y = n, fill = Gene)
) +
  geom_col(width= 0.95) +
  scale_fill_viridis_d(option = "H", end = 0.9) +
  labs(
    y = "Gene count",
    x = NULL,
    fill = NULL
  ) +
  scale_y_continuous(
    breaks = scales::breaks_width(1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.05))
  ) +
  coord_flip() +
  theme_bar_viridis() + theme(legend.position= "none", axis.text.y = element_blank())

bar2

#cluster  
cluster_df <- my_data %>% 
  select(tree_tip, primary_cluster) %>%
  filter(!is.na(tree_tip), !is.na(primary_cluster)) %>%
  distinct() %>%
  left_join(tree_levels, by = c("tree_tip" = "GTDB_GENOME")) %>%
  group_by(tree_tip, tree_pos) %>%
  summarise(
    n = n_distinct(primary_cluster),
    cluster_seq = paste(unique(primary_cluster), collapse = ", "),
    .groups = "drop"
  )

cluster_df <- cluster_df %>% 
  mutate(
    tree_tip = factor(
      tree_tip,
      levels = rev(tip_order),   # reverse so top tip = leftmost bar
      ordered = TRUE
    )
  )
# dont like that NA is white and overall the colourscheme is not really working, thickness of bars

bar3 <- ggplot(
  cluster_df,
  aes(
    x = tree_tip,
    y = n
  )
) +
  geom_col(width = 0.95) + coord_flip() + 
  scale_y_continuous(
    limits = c(0, 90),
    breaks = seq(0, 90, by = 5),
    expand = expansion(mult = c(0, 0.02))) +
  scale_x_discrete(
    expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    x = NULL,
    y = "Cluster count",
    fill = NULL) +
  theme_bar_viridis() +
  theme(legend.position = "none", axis.text.y = element_blank())

# adding captions/ modifications

bar1 <- bar1 + labs(tag = "C")
bar2 <- bar2 + labs(tag = "D")
bar3 <- bar3 + labs(caption = "(Based on data from MIDB)", tag = "E")


# layout 
bottom_row <- (bar1 | bar2 | bar3) +
  plot_layout(widths = c(1, 1, 1))

bottom_row



