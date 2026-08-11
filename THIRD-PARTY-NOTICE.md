# Third-party notices

RNA-seq Explorer bundles third-party software and redistributes data derived
from third-party knowledge resources. Each is the property of its authors and is
used under its own license. This file lists those components with their licenses
and the citations to use in a methods section.

If you publish results produced with this tool, please cite the underlying
methods/resources below in addition to the tool itself (see `CITATION.cff`).

---

## Bundled software (vendored into the single HTML file)

### Chart.js
- Purpose: plotting (volcano, PCA, TF activity, etc.).
- License: MIT.
- Copyright © Chart.js Contributors. https://www.chartjs.org

### canvas2svg
- Purpose: exporting canvas figures as SVG.
- License: MIT.
- Copyright © Gliffy, Inc. and contributors. https://github.com/gliffy/canvas2svg

---

## Analyses run in R (via the "Run real DESeq2 in R" export)

The exported `deseq2_run.R` calls, on the user's own machine:

### DESeq2
- Differential expression, size factors, dispersion, VST.
- License: LGPL (>= 3) (Bioconductor).
- Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change
  and dispersion for RNA-seq data with DESeq2. *Genome Biology* 15:550.
  doi:10.1186/s13059-014-0550-8

### apeglm (optional log2 fold-change shrinkage)
- License: GPL-2.
- Zhu, A., Ibrahim, J.G., Love, M.I. (2019). Heavy-tailed prior distributions for
  sequence count data: removing the noise and preserving large differences.
  *Bioinformatics* 35(12):2084–2092. doi:10.1093/bioinformatics/bty895

### ashr (fallback log2 fold-change shrinkage)
- License: GPL (>= 3).
- Stephens, M. (2017). False discovery rates: a new deal. *Biostatistics*
  18(2):275–294. doi:10.1093/biostatistics/kxw041

---

## Transcription-factor regulons and orthology

### CollecTRI
- TF–target regulon network. The bundled **mouse** network was derived from the
  human CollecTRI, ortholog-translated to mouse (see MGI below) and validated
  against the MGI symbol set. Regulon content © the CollecTRI authors.
- Müller-Dott, S., et al. (2023). Expanding the coverage of regulons from
  high-confidence prior knowledge for accurate estimation of transcription factor
  activities. *Nucleic Acids Research* 51(20):10934–10949. doi:10.1093/nar/gkad841

### decoupleR / OmniPath
- Retrieval of CollecTRI and the ULM activity method.
- decoupleR license: GPL (>= 3). OmniPath data is subject to the licenses of its
  constituent resources; see https://omnipathdb.org/#about for redistribution
  terms and per-resource attribution.
- Badia-i-Mompel, P., et al. (2022). decoupleR: ensemble of computational methods
  to infer biological activities from omics data. *Bioinformatics Advances*
  2(1):vbac016. doi:10.1093/bioadv/vbac016
- Türei, D., et al. (2021). Integrated intra- and intercellular signaling
  knowledge for multicellular omics analysis. *Molecular Systems Biology*
  17:e9923. doi:10.15252/msb.20209923

### MGI — Mouse Genome Informatics (The Jackson Laboratory)
- Human↔mouse ortholog table (`HOM_MouseHumanSequence.rpt`) used to translate the
  network to mouse and to validate gene symbols.
- Data © The Jackson Laboratory; used per MGI's data-use terms
  (https://www.informatics.jax.org). Please cite MGI:
  Baldarelli, R.M., et al. Mouse Genome Informatics: an integrated knowledgebase
  system for the laboratory mouse. *Nucleic Acids Research* (current release).

### Ensembl
- Gene identifier / symbol mapping (`gene_map.csv`).
- Data © EMBL-EBI, released under a no-restriction license with attribution.
  https://www.ensembl.org . Please cite the current Ensembl NAR database paper.

---

## Pathway / enrichment web services (queried at runtime, if online)

These are called only when the Pathways / TF tabs are used, and receive derived
gene identifiers (symbols), never raw counts.

### Reactome
- Pathway over-representation (AnalysisService).
- Gillespie, M., et al. (2022). The Reactome pathway knowledgebase 2022.
  *Nucleic Acids Research* 50(D1):D687–D692. doi:10.1093/nar/gkab1028

### Enrichr (Ma'ayan Lab)
- Hallmark / gene-set over-representation.
- Chen, E.Y., et al. (2013). Enrichr: interactive and collaborative HTML5 gene
  list enrichment analysis tool. *BMC Bioinformatics* 14:128.
- Kuleshov, M.V., et al. (2016). Enrichr: a comprehensive gene set enrichment
  analysis web server 2016 update. *Nucleic Acids Research* 44(W1):W90–W97.
  doi:10.1093/nar/gkw377

### MSigDB Hallmark gene sets (queried via Enrichr)
- Liberzon, A., et al. (2015). The Molecular Signatures Database (MSigDB)
  hallmark gene set collection. *Cell Systems* 1(6):417–425.
  doi:10.1016/j.cels.2015.12.004

---

*Citations above are provided for convenience; verify volume/issue/year against
the primary source before formal citation. Redistribution of any bundled data
(e.g. `collectri_mouse.csv`) remains subject to the source resource's license.*
