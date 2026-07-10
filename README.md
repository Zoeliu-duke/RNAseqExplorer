# RNA-seq Explorer

A browser-based tool for exploring bulk RNA-seq differential expression. It runs a complete DESeq2-style workflow — differential testing, PCA, heatmaps, pathway enrichment, and transcription-factor activity — entirely in your browser, with no installation and no server.

**▶ Use it here: https://zoeliu-duke.github.io/RNAseqExplorer/**

---
## What's new in v37

- **Authoritative DESeq2 in R** — a new round-trip bridge lets you swap the in-browser approximation for reference DESeq2 numbers. **Run real DESeq in R** exports a bundle (contrast count matrix, all-sample matrix, and a ready-to-run `deseq2_<contrast>.R` script); run `Rscript deseq2_<contrast>.R` locally, then drop the results table and blind-VST matrix onto **Import DESeq2 output**. Every view (Volcano, DEG Table, Scatter, PCA, Heatmap, Pathways, TF) switches to the real results, with each file auto-detected from its columns. The whole exchange stays on your own machine.
- **Sample QC tab** — per-sample summaries computed straight from raw counts: library size, detected genes, median-of-ratios size factor, library complexity (top-50 share), mitochondrial fraction, and nearest same-group replicate correlation, with heuristic flags for samples worth inspecting.
- **More heatmap gene-selection modes** — beyond top contrast DEGs, the heatmap can now show most-variable genes, one-vs-rest marker genes per selected group, or a multi-group likelihood-ratio test (LRT) across arbitrary user-chosen groups.
- **Expanded Cook's outlier options** — the four modes are documented in more detail, with the **robust leave-one-out filter** as the recommended default (it only flips a gene's significance flag and never alters reported log₂FC, p, or padj).


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

## Input files

Accepted formats: `.txt`, `.tsv`, `.tab`, `.count`, `.csv` (and `.json` for saved sessions). You can upload multiple files at once. Columns may be tab-, comma-, or whitespace-separated — the delimiter is detected automatically, and lines beginning with `#` are ignored, so raw featureCounts output works as-is.

Two layouts are supported:

- **Per-sample (long):** two columns — gene ID and count — with one file per sample.
- **Count matrix (wide):** a gene-ID column followed by one count column per sample, with an optional header row of sample names.

Gene IDs should be Ensembl IDs; if a gene-symbol column (e.g. `mgi_symbol`) is present, it is used as well. Counts should be raw, un-normalized integers — normalization is performed inside the tool. Dropping a saved `.json` session onto the upload zone restores your samples, settings, and results.

## Running locally

No build step is required. Download `index.html` and open it in any modern browser, or clone the repo:

```
git clone https://github.com/Zoeliu-duke/RNAseqExplorer.git
cd RNAseqExplorer
# then open index.html in your browser
```

## Methods

Counts are pre-filtered by a user-set `rowSums` threshold and normalized with the median-of-ratios method. Gene-wise dispersions are shrunk toward a fitted mean–dispersion trend, and differential expression is called with a DESeq2-style Wald test; p-values are adjusted across all tested genes by the Benjamini–Hochberg procedure. log₂ fold change is computed on size-factor-normalized counts (0.5 pseudocount), with positive values indicating up-regulation in treatment relative to control. PCA and the heatmap use a blind VST independent of group labels. Pathway over-representation is tested by a hypergeometric test against Reactome, KEGG, and WikiPathways; transcription-factor activity uses CollecTRI regulons with the *decoupleR* ULM (Müller-Dott et al., 2023). A full description is available in the **Building Method** panel inside the app.

## AI Disclosure

This tool was developed with the assistance of Claude (Anthropic), a large language model, which generated the application code. Yaxin Liu specified the analytical requirements and design, conducted testing and debugging, and takes responsibility for the tool's implementation and outputs.

## Contact

**Yaxin (Zoe) Liu** — zoe.liu@duke.edu
Edward Miao Lab, Duke University

---

*Version 37 · July 2026*
