# ============================================================
# Turtle OR Repertoire Ecology Analysis
# Author: Lindsay Roberts
# Course: Advanced Statistics for Genomics
# ============================================================

# ----------------------------
# 1. Load packages
# ----------------------------

library(tidyverse)
library(ggplot2)
library(vegan)
library(ape)
library(ggtree)
library(patchwork)

# ----------------------------
# 2. Set working directory
# ----------------------------

setwd("~/Desktop/YoheLab/turtle_OR_repertoire_ecology_analysis")

# Make sure output folders exist
dir.create("data", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# ----------------------------
# 3. Load data
# ----------------------------

or_data <- read.csv("data/turtle_traits_and_ORs.csv")

cat("Data loaded successfully\n")
cat("Number of species:", nrow(or_data), "\n")
print(table(or_data$Habitat))

# ----------------------------
# 4. Reshape OR data
# ----------------------------

coding_long <- or_data %>%
  select(Species, Habitat, Diet, ends_with(".CODING")) %>%
  pivot_longer(
    cols = ends_with(".CODING"),
    names_to = "OR_Family",
    values_to = "Coding_Count"
  ) %>%
  mutate(OR_Family = gsub("\\.CODING$", "", OR_Family))

pseudo_long <- or_data %>%
  select(Species, Habitat, Diet, ends_with(".PSEUDOGENE")) %>%
  pivot_longer(
    cols = ends_with(".PSEUDOGENE"),
    names_to = "OR_Family",
    values_to = "Pseudogene_Count"
  ) %>%
  mutate(OR_Family = gsub("\\.PSEUDOGENE$", "", OR_Family))

or_long <- coding_long %>%
  left_join(pseudo_long, by = c("Species", "Habitat", "Diet", "OR_Family")) %>%
  mutate(
    Total_Count = Coding_Count + Pseudogene_Count,
    Pseudogene_Proportion = ifelse(
      Total_Count > 0,
      Pseudogene_Count / Total_Count,
      NA
    )
  )

# ----------------------------
# 5. Species-level pseudogene summary
# ----------------------------

species_summary <- or_long %>%
  group_by(Species, Habitat, Diet) %>%
  summarise(
    Total_Coding = sum(Coding_Count, na.rm = TRUE),
    Total_Pseudogenes = sum(Pseudogene_Count, na.rm = TRUE),
    Total_OR_Genes = Total_Coding + Total_Pseudogenes,
    Overall_Pseudogene_Proportion = Total_Pseudogenes / Total_OR_Genes,
    .groups = "drop"
  )

species_summary$Habitat <- factor(
  species_summary$Habitat,
  levels = c("Freshwater", "Marine", "Terrestrial")
)

# ----------------------------
# 6. Aim 1: Pseudogene proportion by habitat
# ----------------------------

kw_pseudo <- kruskal.test(
  Overall_Pseudogene_Proportion ~ Habitat,
  data = species_summary
)

cat("\nKruskal-Wallis test for pseudogene proportion:\n")
print(kw_pseudo)

pseudo_by_habitat <- species_summary %>%
  group_by(Habitat) %>%
  summarise(
    n_species = n(),
    Total_Pseudogenes = sum(Total_Pseudogenes),
    Total_Coding = sum(Total_Coding),
    Total_Genes = Total_Pseudogenes + Total_Coding,
    Prop_Pseudo = Total_Pseudogenes / Total_Genes,
    CI_lower = qbeta(0.025, Total_Pseudogenes + 1, Total_Coding + 1),
    CI_upper = qbeta(0.975, Total_Pseudogenes + 1, Total_Coding + 1),
    .groups = "drop"
  )

overall_p0 <- sum(species_summary$Total_Pseudogenes) /
  sum(species_summary$Total_OR_Genes)

binom_results <- pseudo_by_habitat %>%
  rowwise() %>%
  mutate(
    binom_p = binom.test(
      Total_Pseudogenes,
      Total_Genes,
      p = overall_p0
    )$p.value
  ) %>%
  ungroup()

write.csv(
  binom_results,
  "data/Aim1_binomial_results.csv",
  row.names = FALSE
)

habitat_cols <- c(
  "Freshwater" = "#E76F51",
  "Marine" = "#2A9D8F",
  "Terrestrial" = "#3A86FF"
)

p1_box <- ggplot(
  species_summary,
  aes(x = Habitat, y = Overall_Pseudogene_Proportion, fill = Habitat)
) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.5) +
  scale_fill_manual(values = habitat_cols) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Pseudogene proportion by habitat",
    x = "Habitat",
    y = "Pseudogene proportion"
  )

p1_ci <- ggplot(
  pseudo_by_habitat,
  aes(x = Habitat, y = Prop_Pseudo, color = Habitat)
) +
  geom_point(size = 3.5) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.12) +
  scale_color_manual(values = habitat_cols) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Habitat pseudogene proportion with 95% CI",
    x = "Habitat",
    y = "Pseudogene proportion"
  )

fig1 <- p1_box | p1_ci

ggsave(
  "figures/Figure1_Pseudogene_Habitat.pdf",
  fig1,
  width = 11,
  height = 5
)

# ----------------------------
# 7. Aim 2: OR subfamily composition
# ----------------------------

major_families <- c("OR12", "OR14", "OR4", "OR5.8.9", "OR51", "OR55")

composition_species <- or_long %>%
  filter(OR_Family %in% major_families) %>%
  group_by(Species, Habitat, OR_Family) %>%
  summarise(
    Coding_Count = sum(Coding_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Species) %>%
  mutate(Proportion = Coding_Count / sum(Coding_Count)) %>%
  ungroup()

composition_species$Habitat <- factor(
  composition_species$Habitat,
  levels = c("Freshwater", "Marine", "Terrestrial")
)

composition_habitat <- composition_species %>%
  group_by(Habitat, OR_Family) %>%
  summarise(
    Total_Coding = sum(Coding_Count, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(Habitat) %>%
  mutate(Proportion = Total_Coding / sum(Total_Coding)) %>%
  ungroup()

subfamily_cols <- c(
  "OR12" = "orange",
  "OR14" = "forestgreen",
  "OR4" = "purple",
  "OR5.8.9" = "red3",
  "OR51" = "dodgerblue3",
  "OR55" = "saddlebrown"
)

p2_bar <- ggplot(
  composition_habitat,
  aes(x = Habitat, y = Proportion, fill = OR_Family)
) +
  geom_col(color = "black", width = 0.7) +
  scale_fill_manual(values = subfamily_cols) +
  theme_classic(base_size = 12) +
  labs(
    title = "OR subfamily composition by habitat",
    x = "Habitat",
    y = "Proportion of coding OR genes",
    fill = "OR subfamily"
  )

p2_heat <- ggplot(
  composition_habitat,
  aes(x = OR_Family, y = Habitat, fill = Proportion)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(Proportion, 2)), size = 3.5) +
  scale_fill_gradient(low = "white", high = "red3") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Mean OR subfamily proportions",
    x = "OR subfamily",
    y = "Habitat",
    fill = "Proportion"
  )

fig2 <- p2_bar | p2_heat

ggsave(
  "figures/Figure2_OR_Composition_Habitat.pdf",
  fig2,
  width = 12,
  height = 5
)

# ----------------------------
# 8. PCoA analysis
# ----------------------------

composition_matrix <- composition_species %>%
  select(Species, OR_Family, Proportion) %>%
  pivot_wider(
    names_from = OR_Family,
    values_from = Proportion,
    values_fill = 0
  )

metadata <- composition_species %>%
  select(Species, Habitat) %>%
  distinct()

mat <- composition_matrix %>%
  column_to_rownames("Species") %>%
  as.matrix()

dist_mat <- vegdist(mat, method = "bray")

pcoa <- cmdscale(dist_mat, eig = TRUE, k = 2)

pcoa_data <- data.frame(
  Species = rownames(pcoa$points),
  PCoA1 = pcoa$points[, 1],
  PCoA2 = pcoa$points[, 2]
) %>%
  left_join(metadata, by = "Species")

var_explained <- round(
  100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]),
  1
)

cat("\nPCoA variance explained:\n")
print(var_explained[1:2])

p_pcoa <- ggplot(
  pcoa_data,
  aes(x = PCoA1, y = PCoA2, color = Habitat)
) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(aes(group = Habitat), linewidth = 0.7, linetype = 2) +
  scale_color_manual(values = habitat_cols) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCoA of OR subfamily composition",
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)"),
    color = "Habitat"
  )

ggsave(
  "figures/Figure3_PCoA_OR_Composition.pdf",
  p_pcoa,
  width = 6.5,
  height = 5
)

# ----------------------------
# 9. Subfamily tests with BH correction
# ----------------------------

subfamily_tests <- composition_species %>%
  group_by(OR_Family) %>%
  summarise(
    p_value = kruskal.test(Proportion ~ Habitat)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj_BH)

cat("\nSubfamily tests with BH correction:\n")
print(subfamily_tests)

write.csv(
  subfamily_tests,
  "data/Subfamily_Kruskal_BH_results.csv",
  row.names = FALSE
)

p_tests <- ggplot(
  subfamily_tests,
  aes(x = reorder(OR_Family, p_adj_BH), y = -log10(p_adj_BH))
) +
  geom_col(fill = "gray50", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "Subfamily-by-subfamily habitat tests",
    x = "OR subfamily",
    y = "-log10(BH-adjusted p-value)"
  )

ggsave(
  "figures/Figure4_Subfamily_BH_Tests.pdf",
  p_tests,
  width = 6.5,
  height = 4.5
)

# ----------------------------
# 10. Pairwise Wilcoxon tests for significant families
# ----------------------------

sig_families <- subfamily_tests %>%
  filter(p_adj_BH < 0.05) %>%
  pull(OR_Family)

pairwise_results <- composition_species %>%
  filter(OR_Family %in% sig_families) %>%
  group_by(OR_Family) %>%
  group_modify(~{
    pw <- pairwise.wilcox.test(
      x = .x$Proportion,
      g = .x$Habitat,
      p.adjust.method = "BH"
    )
    
    as.data.frame(as.table(pw$p.value)) %>%
      filter(!is.na(Freq)) %>%
      rename(
        Habitat_1 = Var1,
        Habitat_2 = Var2,
        p_adj_BH = Freq
      )
  })

cat("\nPairwise Wilcoxon tests for significant families:\n")
print(pairwise_results)

write.csv(
  pairwise_results,
  "data/Pairwise_Subfamily_BH_results.csv",
  row.names = FALSE
)

# ----------------------------
# 11. Significant subfamily boxplots
# ----------------------------

p_sig_box <- composition_species %>%
  filter(OR_Family %in% sig_families) %>%
  ggplot(aes(x = Habitat, y = Proportion, fill = Habitat)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.4) +
  facet_wrap(~ OR_Family, scales = "free_y") +
  scale_fill_manual(values = habitat_cols) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none") +
  labs(
    title = "Significant OR subfamilies after BH correction",
    x = "Habitat",
    y = "Species-level OR subfamily proportion"
  )

ggsave(
  "figures/Figure7_Significant_Subfamilies_Boxplots.pdf",
  p_sig_box,
  width = 9,
  height = 6
)

# ----------------------------
# 12. Phylogenetic visualization
# ----------------------------

tree <- read.tree("data/june_tree_turtles.tre")

species_summary <- species_summary %>%
  mutate(Tree_Label = gsub(" ", "_", Species))

shared_species <- intersect(tree$tip.label, species_summary$Tree_Label)

cat("\nShared species between tree and OR dataset:\n")
print(length(shared_species))

cat("\nSpecies in OR data but not in tree:\n")
print(setdiff(species_summary$Tree_Label, tree$tip.label))

tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, shared_species))

tree_data <- species_summary %>%
  filter(Tree_Label %in% shared_species)

p_tree_final <- ggtree(tree_pruned, size = 0.4) %<+% tree_data +
  geom_tippoint(
    aes(color = Overall_Pseudogene_Proportion, shape = Habitat),
    size = 3
  ) +
  geom_tiplab(
    aes(label = gsub("_", " ", Tree_Label)),
    size = 1.7,
    align = FALSE
  ) +
  scale_color_gradient(low = "gray85", high = "red3") +
  scale_shape_manual(values = c(
    "Freshwater" = 16,
    "Marine" = 17,
    "Terrestrial" = 15
  )) +
  theme_tree2() +
  labs(
    title = "Phylogenetic distribution of OR pseudogene proportion",
    subtitle = paste("n =", length(shared_species), "species"),
    color = "Pseudogene\nproportion",
    shape = "Habitat"
  )

ggsave(
  "figures/Figure6_Final_Tree.pdf",
  p_tree_final,
  width = 16,
  height = 10
)

# ----------------------------
# 13. Save session info
# ----------------------------

sink("data/session_info.txt")
sessionInfo()
sink()

cat("\nAnalysis complete. Figures saved in figures/. Tables saved in data/.\n")
