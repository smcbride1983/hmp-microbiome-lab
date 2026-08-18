# HMP Microbial Ecology Lab

An interactive teaching app built from processed Human Microbiome Project V3–V5 16S data. The classroom dataset contains 40 visit-one samples from each of five habitats: gut, nasal, oral, skin, and vaginal.

## Run locally

Open the project directory in R and install the two app dependencies:

```r
install.packages(c("shiny", "ggplot2"))
source("check_app.R")
shiny::runApp("app")
```

## Export for GitHub Pages

Install Shinylive and export the app:

```r
install.packages("shinylive")
shinylive::export("app", "site")
httpuv::runStaticServer("site")
```

The included GitHub Actions workflow performs the export and deploys the resulting static site automatically whenever the `main` branch is updated.

In the repository's GitHub settings, choose **Pages → Source → GitHub Actions**.

## Teaching scope

The app supports:

- a guided question → hypothesis → analysis → interpretation workflow;
- sequencing-depth exploration;
- observed richness, Shannon diversity, and Simpson diversity;
- Kruskal–Wallis or ANOVA with optional pairwise comparisons;
- Bray–Curtis or Jaccard PCoA;
- global and optional pairwise permutation-based PERMANOVA;
- a permutation test of multivariate dispersion;
- genus, family, and phylum composition;
- exploratory comparisons of selected taxa among habitats;
- habitat and participant-sex filtering where supported by the data;
- guided ecological interpretation prompts.

## Data provenance

The source data were accessed with the Bioconductor `HMP16SData` package and reduced to a balanced classroom subset. Counts are aggregated to taxonomic ranks and low-prevalence genera were removed for the browser-based activity.
