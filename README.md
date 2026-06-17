# RNA-seq Explorer

A browser-based tool for exploring bulk RNA-seq differential expression. It runs a complete DESeq2-style workflow — differential testing, PCA, heatmaps, pathway enrichment, and transcription-factor activity — entirely in your browser, with no installation and no server.

**▶ Use it here: https://zoeliu-duke.github.io/RNAseqExplorer/**

---

## What it does

Upload a raw count matrix, choose your control and treatment groups, and run the analysis. Everything happens client-side in JavaScript — your count data never leaves your machine.

- **Differential expression** — DESeq2-style Wald test with median-of-ratios normalization, dispersion-trend shrinkage, and BH-corrected FDR (padj).
- **Outlier handling** — optional Cook's-distance filtering with four selectable modes.
- **Volcano plot, expression scatter, and DEG table** — interactive, with per-gene NCBI/UniProt lookups and CSV export.
- **PCA & heatmap** — built on a blind variance-stabilizing transform (VST).
- **Pathway enrichment** — over-representation against Reactome, KEGG, and WikiPathways, plus Hallmark enrichment.
- **Transcription-factor activity** — inferred from CollecTRI regulons using the univariate linear model (ULM) from *decoupleR*.
- **Session save/load** — store and restore your samples, settings, and results as a small JSON file.

Mouse genes are mapped to MGI symbols (~77,600 Ensembl IDs resolved).

## How to use

1. Open the [live site](https://zoeliu-duke.github.io/RNAseqExplorer/).
2. Upload your count files (or drag a saved session onto the upload zone).
3. Select the **control** and **treatment** groups.
4. Adjust the log₂FC, FDR, and count-filter thresholds in the sidebar.
5. Click **Run Analysis** and explore the tabs.

> **Note:** pathway and transcription-factor analyses call external services (Reactome, Enrichr, OmniPath). An internet connection is required for those tabs; the core differential-expression workflow runs fully offline.

## Running locally

No build step is required. Download `index.html` and open it in any modern browser, or clone the repo:

```bash
git clone https://github.com/Zoeliu-duke/RNAseqExplorer.git
cd RNAseqExplorer
# then open index.html in your browser
```

## Methods

Counts are pre-filtered by a user-set `rowSums` threshold and normalized with the median-of-ratios method. Gene-wise dispersions are shrunk toward a fitted mean–dispersion trend, and differential expression is called with a DESeq2-style Wald test; p-values are adjusted across all tested genes by the Benjamini–Hochberg procedure. log₂ fold change is computed on size-factor-normalized counts (0.5 pseudocount), with positive values indicating up-regulation in treatment relative to control. PCA and the heatmap use a blind VST independent of group labels. Pathway over-representation is tested by a hypergeometric test against Reactome, KEGG, and WikiPathways; transcription-factor activity uses CollecTRI regulons with the *decoupleR* ULM (Müller-Dott et al., 2023). A full description is available in the **Building Method** panel inside the app.

## AI Disclosure

This tool was developed with the assistance of Claude (Anthropic), a large language model, which generated the application code. Yaxin Liu specified the analytical requirements and design, conducted testing and debugging, and takes responsibility for the tool's implementation and outputs.

## Contact

**Yaxin Liu** — zoe.liu@duke.edu
Edward Miao Lab, Duke University

---

*Version 35 · June 2026*
