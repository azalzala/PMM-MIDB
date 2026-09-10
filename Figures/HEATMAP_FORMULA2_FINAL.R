library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)

# loading data , run after tree_final.R 
my_data_2 <- read.csv('/Users/ayazalzala/THESIS_FINAL SCRIPTS/representatives_df_2.csv')

# checking right data is present in the final dataframe
unique(my_data_2$representative)
genes_na <- unique(my_data) %>% filter(is.na(Gene.x))

#figuring out cols
names(my_data_2)
names(my_data)
my_data_2 <- my_data_2 %>% select(representative, 
                              ID, 
                              representative_gene_id, 
                              member_accession, 
                              Gene,
                              Product.x,
                              niche.x, 
                              niche.y, 
                              effect_label, 
                              primary_cluster, 
                              tree_tip, 
                              taxonomy,
                              effect_median.x, 
                              effect_lower.x, 
                              effect_upper.x,
                              OR_median.x, 
                              OR_lower.x,
                              OR_upper.x, 
                              prob_median.x, 
                              prob_lower.x, 
                              prob_upper.x)




# trying something else 

gene_effect_labels <- tibble(
  Gene = c(
    "dfrA12",
    "catB3",
    "blaVIM",
    "blaOXA",
    "blaIMP-1",
    "aadA25", 
    "qacL", 
    "lnu(F)", 
    "estX", 
    "HAD family hydrolase",
    "cmlA1", 
    "aadA5", 
    "aac(6')-IIc"
  ),
  combined_effect_label = c(
    "Enriched in wastewater and livestock",
    "Depleted in environmental and livestock",
    "Depleted in wastewater and livestock",
    "Depleted in environmental, wastewater and livestock", 
    "Depleted in environmental, wastewater and livestock", 
    "Enriched in livestock", 
    "Enriched in livestock", 
    "Enriched in livestock", 
    "Enriched in livestock",
    "Enriched in livestock", 
    "Enriched in livestock", 
    "Depleted in livestock", 
    "Depleted in livestock"
  )
)

unique(my_data_2$Gene)
na_gene <- my_data_2 %>% filter(is.na(Gene))
na_gene %>%
  count(representative_gene_id) %>%
  arrange(desc(n))

# confirmed that all NA (gene) are for the gene family AP026508.1_104
# replacing NA with had family hydrolase 

my_data_2 <- my_data_2 %>%
  mutate(
    Gene = if_else(is.na(Gene), "HAD family hydrolase", Gene)
  )

my_data_2 %>%
  count(Gene) %>%
  filter(Gene == "HAD family hydrolase")

# join to plot data 
my_data_2_plot <- my_data_2 %>% 
  filter(!is.na(tree_tip)) %>% 
  # Join your manual labels
  left_join(
    gene_effect_labels,
    by = c("Gene")   # or just "Gene" if label is per gene
  ) %>% 
  mutate(
    tree_tip = factor(
      tree_tip,
      levels = rev(tip_order),
      ordered = TRUE
    )
  )

unique(my_data_2_plot$combined_effect_label)

effect_colours <- c(
  "Depleted in environmental, wastewater and livestock" = "#66C2A5",
  "Depleted in environmental and livestock" = "#8DA0CB",
  "Enriched in wastewater and livestock" = "#E78AC2",
  "Depleted in wastewater and livestock" = "#A6D854",
  "Enriched in livestock" = "#FFD92F", 
  "Depleted in livestock" = "#FFB29A"
)

effect_levels <- names(effect_colours)

my_data_2_plot <- my_data_2_plot %>% 
  mutate(
    combined_effect_label = factor(
      combined_effect_label,
      levels = effect_levels
    )
  )


# plot 
heatmap <- ggplot(my_data_2_plot, aes(x = Gene, y = tree_tip)) +
  geom_tile(
    aes(fill = representative),
    colour = NA
  ) +
  geom_tile(
    aes(colour = combined_effect_label),
    fill = NA,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(14),
    guide = "none",
    na.value = "white"
  ) +
  scale_colour_manual(
    values = effect_colours,
    limits = effect_levels,
    breaks = effect_levels[effect_levels != "No effect"],
    na.value = NA, 
    na.translate = TRUE, 
    name = "Effect",
    drop = FALSE
  ) +
  labs(
    x = "Gene"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(), 
    legend.position = "bottom", 
    panel.grid.major.y = element_line(
      colour = "grey85",
      linewidth = 0.25
    ), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# after checks
library(dplyr)
my_data %>%
  filter(tree_tip == "RS_GCF_0000_18045.1") %>% count()

heatmap
# final tree outlay 
subtree_p | heatmap
top_half <- 

