suppressMessages({library(decoupleR)})
GM  <- "/Users/zoeliu/Documents/Miaolab/Bulk RNAseq/rnaseq explorer/gene_map.csv"
HOM <- "hom.rpt"

# --- MGI symbol universe (what the user's DEG table uses) ---
mgi <- read.csv(GM, stringsAsFactors=FALSE)
mgi_up <- toupper(unique(mgi$symbol))
mgi_by_up <- setNames(unique(mgi$symbol), toupper(unique(mgi$symbol)))
mgiset <- new.env(hash=TRUE); for(s in mgi_up) assign(s, TRUE, mgiset)
inMGI <- function(u) exists(u, envir=mgiset, inherits=FALSE)

# --- MGI HOM human->mouse ortholog map (authoritative; e.g. TP53->Trp53) ---
hom <- read.delim(HOM, stringsAsFactors=FALSE, check.names=FALSE, quote="")
key <- hom[["DB Class Key"]]; org <- hom[["Common Organism Name"]]; sym <- hom[["Symbol"]]
isM <- grepl("mouse", org, ignore.case=TRUE); isH <- grepl("human", org, ignore.case=TRUE)
mouseByKey <- split(sym[isM], key[isM]); humanByKey <- split(sym[isH], key[isH])
orth <- new.env(hash=TRUE)
for(k in intersect(names(humanByKey), names(mouseByKey))){
  ms <- mouseByKey[[k]][1]                       # mouse ortholog symbol
  for(hs in humanByKey[[k]]) assign(toupper(hs), ms, orth)
}
cat("HOM ortholog pairs:", length(ls(orth)), "\n")

# --- Human CollecTRI (same provenance as user's existing human file) ---
h <- get_collectri(organism="human", split_complexes=FALSE)
wcol <- if("weight" %in% names(h)) "weight" else "mor"
h <- data.frame(source=h$source, target=h$target, w=ifelse(h[[wcol]]<0,-1,1), stringsAsFactors=FALSE)
cat("human edges:", nrow(h), " TFs:", length(unique(h$source)), "\n")

META <- c("AP1","NFKB")
src_orth <- 0; src_conv <- 0; src_drop <- 0
tr <- function(sym, tally=FALSE){
  u <- toupper(sym)
  if(u %in% META) return(sym)
  if(exists(u, envir=orth, inherits=FALSE)){                 # (1) MGI orthology
    m <- get(u, orth); if(inMGI(toupper(m))){ if(tally) src_orth<<-src_orth+1; return(mgi_by_up[[toupper(m)]]) }
  }
  if(inMGI(u)){ if(tally) src_conv<<-src_conv+1; return(mgi_by_up[[u]]) }  # (2) convention, MGI-validated
  if(tally) src_drop<<-src_drop+1; NA_character_
}
syms <- unique(c(h$source,h$target)); mp <- setNames(vapply(syms,tr,character(1)), syms)
h$ms <- ifelse(h$source %in% META, h$source, mp[h$source]); h$mt <- mp[h$target]
keep <- !is.na(h$ms) & !is.na(h$mt)
mnet <- unique(data.frame(source=h$ms[keep], target=h$mt[keep], weight=h$w[keep], stringsAsFactors=FALSE))
mnet <- mnet[!is.na(mnet$source)&!is.na(mnet$target),]
# provenance tally over distinct human symbols
for(s in syms) invisible(tr(s, tally=TRUE))
write.csv(mnet, "collectri_mouse.csv", row.names=FALSE, quote=FALSE)

cat("\n=== FINAL collectri_mouse.csv ===\n")
cat("edges:", nrow(mnet), " TFs:", length(unique(mnet$source)),
    sprintf(" (%.1f%% of human edges kept)\n",100*nrow(mnet)/nrow(h)))
cat(sprintf("symbol provenance (distinct): orthology=%d  convention=%d  dropped=%d\n", src_orth,src_conv,src_drop))
cat(sprintf("TF sources MGI-valid: %.1f%%  targets MGI-valid: %.1f%%\n",
    100*mean(toupper(setdiff(unique(mnet$source),META))%in%mgi_up),
    100*mean(toupper(unique(mnet$target))%in%mgi_up)))
for(g in c("Trp53","Stat1","Stat3","Nfkb1","Rela","Irf1","Irf3","Myc","Jun","Fos","Cebpb","Spi1","Nr3c1"))
  cat(sprintf("  %-7s targets: %d\n", g, sum(mnet$source==g)))
cat("Cxcl2 regulators:", paste(head(sort(unique(mnet$source[mnet$target=="Cxcl2"])),14),collapse=", "),"\n")
