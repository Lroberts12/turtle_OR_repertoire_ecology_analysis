# Turtle Olfactory Receptor Repertoire Ecology Analysis

This repository contains R code, figures, and summary tables for an Advanced Statistics final project analyzing olfactory receptor (OR) repertoires across turtle habitats.

## Project Question

Do freshwater, marine, and terrestrial turtles differ in OR pseudogene proportion and OR subfamily composition?

## Data

The dataset contains OR coding and pseudogene counts across approximately 50 turtle species. Genome assemblies were obtained from NCBI and processed through a custom OR annotation pipeline. Ecological metadata were compiled from literature and resources including Animal Diversity Web.

## Analyses

- Species-level pseudogene proportion by habitat
- Binomial confidence intervals
- Kruskal-Wallis tests
- OR subfamily composition plots
- PCoA of OR subfamily composition
- Subfamily-by-subfamily tests with Benjamini-Hochberg correction
- Pairwise Wilcoxon tests
- Phylogenetic visualization of habitat and pseudogene proportion

## Repository Structure

- `data/` input datasets and results tables
- `figures/` generated figures
- `scripts/` R analysis script

## Notes

AI was used for code troubleshooting, grammar editing, and project organization.