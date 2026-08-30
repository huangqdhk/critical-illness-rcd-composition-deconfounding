## M3 腿1 步骤4b：coloc.abf 共定位（deCODE pQTL × 结局）
## 输入：_intermediate/M3_coloc_{gene}_{outcome}.csv
## 输出：Table_S15e_coloc_{Results,snp_pp}.csv
## p12=1e-5（主）、1e-6/5e-5（敏感性）；deCODE type=quant；结局 type=cc
## 口径与 Table S15c 对齐；coloc 5.2.3
suppressPackageStartupMessages({library(data.table); library(coloc)})
PROJ <- "E:/SCI/proj"
INT  <- file.path(PROJ, "_intermediate")
OUT  <- file.path(PROJ, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")

NCASE <- c(FinnGen_R10_ARDS = 357, IEU_sepsis_ieu_b_69 = 10154)
NCONTROL <- c(FinnGen_R10_ARDS = 406536, IEU_sepsis_ieu_b_69 = 454764)
CC_PRIOR <- c("FinnGen_R10_ARDS" = 357 / (357 + 406536),
              "IEU_sepsis_ieu-b-69" = 10154 / (10154 + 454764))

run_coloc <- function(df, p12, oc) {
  d1 <- list(snp = df$snp, beta = df$pqtl_beta, varbeta = df$pqtl_se^2,
             MAF = df$pqtl_maf, N = df$pqtl_n, sdY = 1, type = "quant")
  d2 <- list(snp = df$snp, beta = df$out_beta, varbeta = df$out_se^2,
             MAF = df$out_maf, s = CC_PRIOR[[oc]], type = "cc")
  res <- coloc.abf(d1, d2, p12 = p12)
  list(summary = res$summary, results = res$results)
}

res_rows <- list(); pp_rows <- list()
for (gene in c("ITPR1", "SOD2", "TFRC")) {
  for (oc in c("FinnGen_R10_ARDS", "IEU_sepsis_ieu-b-69")) {
    f <- file.path(INT, sprintf("M3_coloc_%s_%s.csv", gene, oc))
    if (!file.exists(f)) next
    df <- fread(f)
    df <- df[is.finite(pqtl_beta) & is.finite(out_beta) &
             is.finite(pqtl_se) & is.finite(out_se) & pqtl_se > 0 & out_se > 0]
    df <- df[!duplicated(snp)]
    if (nrow(df) < 50) next
    for (p12 in c(1e-5, 1e-6, 5e-5)) {
      rr <- run_coloc(df, p12, oc)
      s <- as.list(rr$summary)
      res_rows[[length(res_rows) + 1]] <- data.frame(
        gene = gene, outcome = oc, p12 = p12,
        n_snps_used = s$nsnps,
        PP.H0 = s$PP.H0.abf, PP.H1 = s$PP.H1.abf, PP.H2 = s$PP.H2.abf,
        PP.H3 = s$PP.H3.abf, PP.H4 = s$PP.H4.abf,
        conclusion = ifelse(s$PP.H4.abf >= 0.8, "colocalization supported (PP.H4>=0.8)",
                     ifelse(s$PP.H4.abf >= 0.5, "suggestive (0.5<=PP.H4<0.8)",
                     "no colocalization (PP.H4<0.5)")),
        gene_set_version = "Mitoxy-80_v1.0", score_version = "m3pqtl_v1.0",
        stringsAsFactors = FALSE)
      if (p12 == 1e-5) {
        pp <- rr$results
        pp_rows[[length(pp_rows) + 1]] <- data.frame(
          gene = gene, outcome = oc, p12 = p12,
          snp = pp$snp,
          snp_pp_h4 = pp$SNP.PP.H4,
          gene_set_version = "Mitoxy-80_v1.0", score_version = "m3pqtl_v1.0",
          stringsAsFactors = FALSE)
      }
    }
  }
}
resDT <- rbindlist(res_rows, fill = TRUE)
ppDT  <- rbindlist(pp_rows, fill = TRUE)
fwrite(resDT, file.path(OUT, "Table_S15e_coloc_Results.csv"))
fwrite(ppDT,  file.path(OUT, "Table_S15e_coloc_snp_pp.csv"))
cat("\n===== M3 pQTL coloc DONE =====\n")
print(resDT[order(outcome, gene, p12),
            .(gene, outcome, p12, n_snps_used, PP.H4 = round(PP.H4, 4),
              PP.H3 = round(PP.H3, 4), PP.H1 = round(PP.H1, 4), conclusion)])
