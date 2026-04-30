# Turtle Olfactory Receptor Repertoire Ecology Analysis

This repository contains R code, figures, and summary tables for an Advanced Statistics project analyzing olfactory receptor (OR) repertoires across turtle habitats.

---

## Project Question

Do freshwater, marine, and terrestrial turtles differ in:
1. OR pseudogene proportion (functional integrity)
2. OR subfamily composition (repertoire structure)

---

## Data

The dataset includes OR coding and pseudogene counts across ~50 turtle species.  
Genomes were obtained from NCBI and processed through a custom olfactory receptor annotation pipeline.  

Ecological metadata (habitat and diet) were curated from literature sources, including Animal Diversity Web.
Key steps in the pipeline included:

- Sequence similarity searches (BLAST-based filtering of candidate OR genes)

- Translation-guided alignment and curation

- Phylogenetic clustering and orthogroup inference using IQ-TREE and UPhO

- Classification of sequences into OR subfamilies

- Identification of pseudogenes based on disruptions in coding sequence

The resulting dataset provides standardized OR gene counts across species, enabling direct comparison of both:

1. OR repertoire integrity (coding vs. pseudogene)

2. OR subfamily composition

Ecological metadata (habitat and diet) were curated from primary literature and online databases, including Animal Diversity Web (https://animaldiversity.org/). Habitat categories were standardized into freshwater, marine, and terrestrial groups for analysis.

Only species with:

- High-quality genome assemblies  

- Consistent OR annotations using the same pipeline  

- Available ecological metadata  

were included in the final dataset to ensure comparability across species.
---

## Methods Overview

- Pseudogene proportion modeled as a binomial process
- Habitat differences tested using Kruskal-Wallis tests
- OR subfamily composition normalized per species
- Multivariate structure analyzed using Bray-Curtis distance and PCoA
- Subfamily-specific differences tested with Kruskal-Wallis + Benjamini-Hochberg correction
- Pairwise Wilcoxon tests for significant families
- Phylogenetic visualization using a pruned turtle tree

---

## Results

### 1. Pseudogene Proportion Differs by Habitat

Kruskal-Wallis test showed a significant difference in pseudogene proportion across habitats  
(p ≈ 0.0007), indicating that OR repertoire integrity varies by ecological niche.

Interpretation:  
Habitat is associated with OR gene degradation patterns, suggesting ecological differences in olfactory reliance.

---

### 2. OR Subfamily Composition Shows Habitat Structure

PCoA analysis explained:
- Axis 1: ~54.8%
- Axis 2: ~24.2%

Species show partial clustering by habitat, indicating that OR repertoires are structured but not completely separated.

Interpretation:  
Habitat influences OR composition, but variation is continuous rather than discrete.

---

### 3. Specific OR Subfamilies Differ Across Habitats

After multiple testing correction (BH):

Significant families:
- OR4
- OR51
- OR12
- OR14

Interpretation:  
Ecological differences are driven by specific OR gene families rather than the entire repertoire.

---

### 4. Phylogenetic Signal is Present

Mapping traits onto the phylogeny shows clustering of similar OR profiles among closely related species.

Interpretation:  
Phylogeny and habitat are confounded — closely related species share similar OR repertoires.

---

## Figures

- **Figure 1** – Pseudogene proportion by habitat  
- **Figure 2** – OR subfamily composition by habitat  
- **Figure 3** – PCoA of OR subfamily composition  
- **Figure 4** – Subfamily-level statistical tests (BH corrected)  
- **Figure 5** – Phylogenetic distribution of pseudogene proportion  
- **Figure 6** – Significant OR subfamilies across habitats  

---

## Repository Structure

- `data/` – input datasets and statistical results  
- `figures/` – final figures used in analysis  
- `scripts/` – R script for full analysis pipeline  

---

## Key Takeaways

- OR pseudogene proportion differs significantly across habitats  
- OR composition shows partial ecological clustering  
- Specific gene families drive ecological differences  
- Phylogenetic relatedness limits statistical independence  

---


## Notes

AI was used for:

- Code troubleshooting

- Debugging

- Formatting and clarity improvements

All analyses and interpretations were performed by the author.
