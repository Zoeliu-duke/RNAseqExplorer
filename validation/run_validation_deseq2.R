# =============================================================================
# run_validation_deseq2.R
#   Reproducible validation backbone for RNA-seq Explorer.
#   Takes a real mouse dataset (GSE60450, basal vs luminal mammary epithelium),
#   maps Entrez -> MGI symbol, and runs GENUINE DESeq2. The outputs are the
#   reference ("expected") numbers that RNA-seq Explorer is validated against.
#
# Input : GSE60450_counts.txt.gz  (GEO supplementary gene-wise counts)
# Output: mouse_counts.csv             tool-ready raw count matrix (symbol + samples)
#         sample_metadata.csv          sample -> group
#         expected_deseq2_results.csv  genuine DESeq2 results (the reference)
#         sessionInfo.txt              exact package versions for provenance
# Run   : Rscript run_validation_deseq2.R
# =============================================================================
suppressMessages({library(DESeq2); library(org.Mm.eg.db)})

# Fetch the public GEO counts if not already present (self-contained reproduce)
gz  <- "GSE60450_counts.txt.gz"
url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60450/suppl/GSE60450_Lactation-GenewiseCounts.txt.gz"
if (!file.exists(gz)) download.file(url, gz, mode = "wb")

raw <- read.delim(gz, check.names = FALSE)
entrez <- as.character(raw$EntrezGeneID)
cnt <- as.matrix(raw[, -(1:2)])                      # drop EntrezGeneID, Length
rownames(cnt) <- entrez
colnames(cnt) <- sub("_.*", "", colnames(cnt))       # MCL1-DG_BC2...  -> MCL1-DG
cat("samples:", ncol(cnt), "| genes(entrez):", nrow(cnt), "\n")

# Entrez -> MGI symbol; drop unmapped; sum counts over duplicate symbols
sym <- mapIds(org.Mm.eg.db, keys = entrez, column = "SYMBOL",
              keytype = "ENTREZID", multiVals = "first")
ok  <- !is.na(sym)
cntsym <- rowsum(cnt[ok, ], group = sym[ok])         # rows = unique MGI symbols
cat("mapped to", nrow(cntsym), "unique MGI symbols\n")

# Cell type from sample tag: MCL1-D* = basal, MCL1-L* = luminal (GSE60450 design)
tag  <- sub("^MCL1-", "", colnames(cntsym))
grp  <- factor(ifelse(startsWith(tag, "D"), "basal", "luminal"),
               levels = c("basal", "luminal"))       # basal = reference level
print(table(colnames(cntsym), grp))

# Light pre-filter (both methods see the same matrix)
cm <- round(cntsym[rowSums(cntsym) >= 10, ])
cat("genes after rowSums>=10 filter:", nrow(cm), "\n")

# ---- tool-ready inputs ----
write.csv(data.frame(gene = rownames(cm), cm, check.names = FALSE),
          "mouse_counts.csv", row.names = FALSE, quote = FALSE)
write.csv(data.frame(sample = colnames(cm), group = as.character(grp)),
          "sample_metadata.csv", row.names = FALSE, quote = FALSE)

# ---- GENUINE DESeq2 (positive log2FC = up in luminal vs basal) ----
dds <- DESeqDataSetFromMatrix(cm, DataFrame(group = grp), ~group)
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "luminal", "basal"))
out <- data.frame(gene = rownames(res),
                  baseMean = res$baseMean,
                  log2FoldChange = res$log2FoldChange,
                  lfcSE = res$lfcSE, stat = res$stat,
                  pvalue = res$pvalue, padj = res$padj)
out <- out[order(out$padj), ]
write.csv(out, "expected_deseq2_results.csv", row.names = FALSE, quote = FALSE)

sig <- sum(out$padj < 0.05 & abs(out$log2FoldChange) > 1, na.rm = TRUE)
cat(sprintf("DESeq2 DEGs (padj<0.05 & |lfc|>1): %d\n", sig))
cat("--- marker sanity check (expect luminal-high +, basal-high -) ---\n")
for (g in c("Krt8","Krt18","Csn2","Elf5",  "Krt5","Krt14","Acta2","Oxtr"))
  if (g %in% out$gene) cat(sprintf("  %-7s lfc=% .2f  padj=%.1e\n",
        g, out$log2FoldChange[out$gene==g], out$padj[out$gene==g]))

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
cat("DONE\n")
