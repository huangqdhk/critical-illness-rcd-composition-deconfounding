# M10B GSE216009 RDS extraction (ASCII-only script; run via Rscript with args)
# Usage: Rscript extract_gse216009.R <rds_path> <outdir>
# Extracts: meta.data subset, RNA data/counts submatrix for manifest genes + MT genes.
args <- commandArgs(trailingOnly = TRUE)
rds_path <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE)

suppressMessages(library(Seurat))
cat("reading RDS ...\n")
obj <- readRDS(rds_path)
cat("RDS loaded\n")

meta <- obj@meta.data
keep_cols <- intersect(c("sample_id", "source", "diagnosis", "age", "sex", "batch",
                         "fine_annot", "broad_annot", "active.ident"), colnames(meta))
write.csv(meta[, keep_cols, drop = FALSE], file.path(outdir, "meta.csv"), row.names = TRUE)
cat("meta saved\n")

man <- read.csv("E:/SCI/_rtmp/manifest_ascii.csv", stringsAsFactors = FALSE)
aliases <- unlist(strsplit(paste(man$aliases, collapse = ";"), "[;,]"))
targets <- toupper(unique(c(
  man$gene_symbol, man$hgnc_symbol, aliases,
  c("MT-ATP6", "MT-ATP8", "MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB",
    "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6")
)))
targets <- targets[!is.na(targets) & nchar(targets) > 0]
cat("target symbols:", length(targets), "\n")

rn <- rownames(obj[["RNA"]])
hit <- which(toupper(rn) %in% targets)
cat("genes matched:", length(hit), "\n")
write.csv(data.frame(gene = rn[hit], stringsAsFactors = FALSE),
          file.path(outdir, "genes_hit.csv"), row.names = FALSE)

sub_data <- GetAssayData(obj, assay = "RNA", layer = "data")[hit, , drop = FALSE]
write.csv(as.matrix(sub_data), file.path(outdir, "data_sub.csv"), row.names = TRUE)
cat("data_sub saved\n")

sub_counts <- GetAssayData(obj, assay = "RNA", layer = "counts")[hit, , drop = FALSE]
write.csv(as.matrix(sub_counts), file.path(outdir, "counts_sub.csv"), row.names = TRUE)
cat("counts_sub saved\n")

cat("DONE\n")
