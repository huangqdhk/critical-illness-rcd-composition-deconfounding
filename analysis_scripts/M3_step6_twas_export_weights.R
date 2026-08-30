## M3 腿3 步骤6：导出 FUSION 权重（80 基因 × 肺/全血）
## 输出：_intermediate/M3_twas_weights_{tissue}.csv
##  每行：gene, model_method, rsid, chr, pos_hg38, a1, a2, weight
##  含全部方法（top1/enet/susie/lasso）与 cv rsq（供 best-CV 选择）
library(data.table)
PROJ <- "E:/SCI/proj"
INT  <- file.path(PROJ, "_intermediate")
man <- fread(file.path(PROJ, "04_AUDIT_GOVERNANCE", "Mitoxyperilysis_Gene_Manifest_v1.0.csv"))
targets <- unique(man$ensembl_gene_id)
cat("target genes:", length(targets), "\n")

for (tiss in c("GTExv8.EUR.Lung", "GTExv8.EUR.Whole_Blood")) {
  d <- file.path(PROJ, "00_RAW_DATA/GWAS/gtex_v8_eQTL", tiss, tiss)
  allf <- list.files(d, pattern = "wgt\\.RDat$")
  out_rows <- list()
  n_gene <- 0
  for (g in targets) {
    hit <- grep(paste0("^", g, "\\."), allf, value = TRUE)
    if (length(hit) == 0) next
    n_gene <- n_gene + 1
    e <- new.env()
    load(file.path(d, hit[1]), envir = e)
    snps <- e$snps
    w <- e$wgt.matrix
    cv <- e$cv.performance
    if (!is.data.frame(snps) || nrow(snps) == 0) next
    for (meth in colnames(w)) {
      wt <- w[, meth]
      keep <- which(!is.na(wt) & wt != 0)
      if (length(keep) == 0) next
      rsq <- NA_real_
      if (meth %in% colnames(cv)) rsq <- cv["rsq", meth]
      out_rows[[length(out_rows) + 1]] <- data.frame(
        gene = g, model = meth, rsid = snps$V2[keep],
        chr = snps$V1[keep], pos_hg38 = snps$V4[keep],
        a1 = snps$V5[keep], a2 = snps$V6[keep],
        weight = as.numeric(wt[keep]), cv_rsq = rsq,
        stringsAsFactors = FALSE)
    }
  }
  cat(tiss, "genes with weights:", n_gene, "rows:", length(out_rows), "\n")
  if (length(out_rows)) {
    dt <- rbindlist(out_rows)
    fwrite(dt, file.path(INT, paste0("M3_twas_weights_", tiss, ".csv")))
  }
}
cat("DONE\n")
