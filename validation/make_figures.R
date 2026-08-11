# Independent concordance check + figure from tool_results.csv (tool vs DESeq2)
d <- read.csv("tool_results.csv", stringsAsFactors = FALSE)
d <- d[is.finite(d$deseq_lfc) & is.finite(d$tool_lfc), ]
n <- nrow(d)

pear  <- cor(d$deseq_lfc, d$tool_lfc, method = "pearson")
spear <- cor(d$deseq_lfc, d$tool_lfc, method = "spearman")

deseq_sig <- with(d, is.finite(deseq_padj) & deseq_padj < 0.05 & abs(deseq_lfc) > 1)
tool_sig  <- with(d, is.finite(tool_padj)  & tool_padj  < 0.05 & abs(tool_lfc)  > 1)
inter <- sum(deseq_sig & tool_sig); uni <- sum(deseq_sig | tool_sig)
sens  <- inter / sum(deseq_sig); prec <- inter / sum(tool_sig); jacc <- inter / uni
TN <- sum(!deseq_sig & !tool_sig); spec <- TN / (TN + sum(tool_sig & !deseq_sig))
either <- deseq_sig | tool_sig
signc  <- mean((d$deseq_lfc[either] > 0) == (d$tool_lfc[either] > 0))

cat(sprintf("n=%d\nPearson=%.4f  Spearman=%.4f\nsign_concordance=%.4f\n", n, pear, spear, signc))
cat(sprintf("DESeq2 DEGs=%d  tool DEGs=%d  overlap=%d\nJaccard=%.3f sens=%.3f prec=%.3f spec=%.3f\n",
            sum(deseq_sig), sum(tool_sig), inter, jacc, sens, prec, spec))

# ---- figure: LFC concordance ----
png("lfc_concordance.png", width = 1500, height = 760, res = 150)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
lim <- c(-13, 13)
col_pt <- ifelse(deseq_sig & tool_sig, "#1f77b4",
          ifelse(deseq_sig | tool_sig, "#ff7f0e", "#999999"))
plot(d$deseq_lfc, d$tool_lfc, pch = 16, cex = .28,
     col = adjustcolor(col_pt, .35), xlim = lim, ylim = lim,
     xlab = "genuine DESeq2  log2FC", ylab = "RNA-seq Explorer  log2FC",
     main = "log2 fold-change concordance")
abline(0, 1, col = "black", lwd = 1.4, lty = 2)
legend("topleft", bty = "n", cex = .8,
       legend = c(sprintf("Pearson r = %.3f", pear),
                  sprintf("Spearman = %.3f", spear),
                  sprintf("n = %s genes", format(n, big.mark = ","))))
legend("bottomright", bty = "n", cex = .72, pch = 16,
       col = c("#1f77b4", "#ff7f0e", "#999999"),
       legend = c("DEG in both", "DEG in one", "not a DEG"))

# ---- figure: significance (-log10 padj) concordance ----
tp <- pmin(-log10(pmax(d$tool_padj,  1e-300)), 60)
dp <- pmin(-log10(pmax(d$deseq_padj, 1e-300)), 60)
ok <- is.finite(tp) & is.finite(dp)
plot(dp[ok], tp[ok], pch = 16, cex = .28, col = adjustcolor("#444444", .3),
     xlab = "genuine DESeq2  -log10(padj)", ylab = "RNA-seq Explorer  -log10(padj)",
     main = "significance concordance")
abline(0, 1, col = "black", lwd = 1.2, lty = 2)
abline(h = -log10(.05), v = -log10(.05), col = "#cc0000", lwd = .9, lty = 3)
dev.off()
cat("wrote lfc_concordance.png\n")

writeLines(c(
  sprintf("n_common_genes\t%d", n),
  sprintf("pearson_lfc\t%.4f", pear),
  sprintf("spearman_lfc\t%.4f", spear),
  sprintf("sign_concordance\t%.4f", signc),
  sprintf("deseq2_DEGs\t%d", sum(deseq_sig)),
  sprintf("tool_DEGs\t%d", sum(tool_sig)),
  sprintf("overlap_DEGs\t%d", inter),
  sprintf("jaccard\t%.4f", jacc),
  sprintf("sensitivity\t%.4f", sens),
  sprintf("precision\t%.4f", prec),
  sprintf("specificity\t%.4f", spec)
), "concordance_metrics.tsv")
