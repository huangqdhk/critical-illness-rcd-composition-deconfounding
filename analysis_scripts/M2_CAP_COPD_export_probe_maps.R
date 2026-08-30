# ============================================================
# M2_CAP_COPD_export_probe_maps.R  (v3: pure ASCII, paths via env vars)
# Exports probe->symbol maps for GPL570 (HG-U133 Plus 2.0) and
# GPL23159 (Clariom S Human) as CSV for the Python pipeline.
# Run via _run_R_export.py (sets TMP/TEMP + M2_OUTDIR + M2_MANIFEST).
# ============================================================
options(timeout = 900)
suppressPackageStartupMessages({
  library(hgu133plus2.db)
  library(clariomshumantranscriptcluster.db)
  library(AnnotationDbi)
})

outdir <- Sys.getenv("M2_OUTDIR")
if (!nzchar(outdir)) stop("M2_OUTDIR env var not set")
dir.create(outdir, showWarnings = FALSE)

# ---------- 1. GPL570 (HG-U133 Plus 2.0) ----------
cat("GPL570 keytypes:", paste(keytypes(hgu133plus2.db), collapse = ","), "\n")
keys570 <- keys(hgu133plus2.db, keytype = "PROBEID")
cat("GPL570 n probes:", length(keys570), "\n")
map570 <- select(hgu133plus2.db, keys = keys570, columns = "SYMBOL", keytype = "PROBEID")
map570 <- map570[!is.na(map570$SYMBOL), ]
write.csv(map570, file.path(outdir, "GPL570_probe2symbol_hgu133plus2db.csv"),
          row.names = FALSE, quote = FALSE)
cat("GPL570 exported:", nrow(map570), "probe-gene pairs\n")

# ---------- 2. GPL23159 (Clariom S Human) ----------
keys23159 <- keys(clariomshumantranscriptcluster.db, keytype = "PROBEID")
cat("GPL23159 n probes:", length(keys23159), "\n")
map23159 <- select(clariomshumantranscriptcluster.db, keys = keys23159,
                   columns = "SYMBOL", keytype = "PROBEID")
map23159 <- map23159[!is.na(map23159$SYMBOL), ]
write.csv(map23159, file.path(outdir, "GPL23159_probe2symbol_clariomsdb.csv"),
          row.names = FALSE, quote = FALSE)
cat("GPL23159 exported:", nrow(map23159), "probe-gene pairs\n")

# ---------- 3. coverage self-check vs 80-gene manifest ----------
mf <- Sys.getenv("M2_MANIFEST")
if (nzchar(mf) && file.exists(mf)) {
  manifest <- read.csv(mf, stringsAsFactors = FALSE)
  genes80 <- unique(c(toupper(manifest$hgnc_symbol), toupper(manifest$gene_symbol)))
  cat("\n=== 80-gene coverage ===\n")
  cat("GPL570   covered:", sum(genes80 %in% toupper(map570$SYMBOL)), "/ 80\n")
  cat("GPL23159 covered:", sum(genes80 %in% toupper(map23159$SYMBOL)), "/ 80\n")
  cat("GPL570   missing:", paste(setdiff(genes80, toupper(map570$SYMBOL)), collapse = ","), "\n")
  cat("GPL23159 missing:", paste(setdiff(genes80, toupper(map23159$SYMBOL)), collapse = ","), "\n")
}
cat("\nDONE\n")
