## M3 腿1 步骤3：pQTL-MR 估计（IVW / Wald / 加权中位数 / MR-Egger / Q / LOO）
## 输入：_intermediate/M3_pqtl_harmonized.csv（step2 已 harmonize）
## 输出：02_SUPPLEMENTARY_TABLES/SUPPLEMENTARY_Tables_CSV/Table_S15e_pQTL_MR_{Results,Sensitivity,Instruments}.csv
## 口径与 Table S15a/b（twmr_finngen.R）一致；score_version=m3pqtl_v1.0
library(data.table)
suppressPackageStartupMessages(library(data.table))
set.seed(42)
PROJ <- "E:/SCI/proj"
INT  <- file.path(PROJ, "_intermediate")
OUT  <- file.path(PROJ, "02_SUPPLEMENTARY_TABLES", "SUPPLEMENTARY_Tables_CSV")

d <- fread(file.path(INT, "M3_pqtl_harmonized.csv"))
d[, `:=`(gene = as.character(gene), source = as.character(source),
         rsid = as.character(rsid), outcome = as.character(outcome))]

## ---- MR 函数（与 twmr_finngen.R 同口径）----
mr_ivw <- function(b_e, b_o, se_o) {
  w <- 1 / se_o^2; b <- sum(w * b_e * b_o) / sum(w * b_e^2)
  se <- sqrt(1 / sum(w * b_e^2))
  Q <- sum(w * (b_o - b_e * b)^2)
  list(b = b, se = se, Q = Q)
}
mr_egger <- function(b_e, b_o, se_o) {
  w <- 1 / se_o^2
  fit <- lm(b_o ~ b_e, weights = w)
  s <- summary(fit)$coefficients
  list(b = s["b_e", "Estimate"], se = s["b_e", "Std. Error"],
       p = s["b_e", "Pr(>|t|)"],
       ic = s["(Intercept)", "Estimate"], ic_se = s["(Intercept)", "Std. Error"],
       ic_p = s["(Intercept)", "Pr(>|t|)"])
}
mr_wald <- function(b_e, b_o, se_o, se_e) {
  b <- b_o / b_e
  se <- abs(b) * sqrt((se_o / b_o)^2 + (se_e / b_e)^2)
  list(b = b, se = se)
}
wmedian <- function(r, ww) {
  o <- order(r); cw <- cumsum(ww[o] / sum(ww))
  r[o][min(which(cw >= 0.5))]
}
mr_wmed <- function(b_e, b_o, se_o) {
  r <- b_o / b_e; w <- 1 / (se_o / b_e)^2
  b <- wmedian(r, w)
  nn <- length(b_e); B <- 1000; boot <- numeric(B)
  for (j in seq_len(B)) {
    s <- sample(nn, replace = TRUE)
    boot[j] <- wmedian(b_o[s] / b_e[s], 1 / (se_o[s] / b_e[s])^2)
  }
  list(b = b, se = sd(boot))
}

mk <- function(exposure, outcome, method, nsnp, b, se, Q = NA, Qp = NA, Fm = NA,
               ic = NA, ic_se = NA, ic_p = NA, exp_id) {
  p <- 2 * pnorm(-abs(b / se))
  list(exposure = exposure, outcome = outcome, method = method, nsnp = nsnp,
       b = b, se = se, pval = p, OR = exp(b),
       OR_CI_low = exp(b - 1.96 * se), OR_CI_high = exp(b + 1.96 * se),
       exposure_id = exp_id, outcome_id = outcome,
       Q = Q, Q_pval = Qp, F_stat = Fm,
       intercept = ic, intercept_se = ic_se, intercept_pval = ic_p,
       gene_set_version = "Mitoxy-80_v1.0", score_version = "m3pqtl_v1.0")
}

res <- list(); sens <- list(); loo <- list()
for (g in sort(unique(d$gene))) {
  for (src in sort(unique(d$source))) {
    sub_all <- d[gene == g & source == src]
    if (nrow(sub_all) == 0) next
    for (oc in sort(unique(sub_all$outcome))) {
      sub <- sub_all[outcome == oc]
      k <- nrow(sub)
      if (k == 0) next
      exp_id <- paste0(src, ":", g)
      Fm <- mean(sub$F_e)
      sens[[length(sens) + 1]] <- list(
        exposure = paste0(src, ":", g), outcome = oc, source = src,
        n_instruments = k, mean_F = Fm, min_F = min(sub$F_e))
      if (k == 1) {
        w <- mr_wald(sub$beta_e, sub$beta_o, sub$se_o, sub$se_e)
        res[[length(res) + 1]] <- mk(paste0(src, ":", g), oc, "Wald ratio", 1,
                                     w$b, w$se, Fm = Fm, exp_id = exp_id)
      } else {
        iv <- mr_ivw(sub$beta_e, sub$beta_o, sub$se_o)
        Qp <- pchisq(iv$Q, k - 1, lower.tail = FALSE)
        res[[length(res) + 1]] <- mk(paste0(src, ":", g), oc,
                                     "Inverse variance weighted", k,
                                     iv$b, iv$se, Q = iv$Q, Qp = Qp, Fm = Fm,
                                     exp_id = exp_id)
        if (k >= 3) {
          eg <- mr_egger(sub$beta_e, sub$beta_o, sub$se_o)
          res[[length(res) + 1]] <- mk(paste0(src, ":", g), oc, "MR Egger", k,
                                       eg$b, eg$se, Fm = Fm, ic = eg$ic,
                                       ic_se = eg$ic_se, ic_p = eg$ic_p,
                                       exp_id = exp_id)
          wm <- mr_wmed(sub$beta_e, sub$beta_o, sub$se_o)
          res[[length(res) + 1]] <- mk(paste0(src, ":", g), oc,
                                       "Weighted median", k, wm$b, wm$se,
                                       Fm = Fm, exp_id = exp_id)
        }
        ## leave-one-out（IVW）
        for (i in seq_len(k)) {
          iv_i <- mr_ivw(sub$beta_e[-i], sub$beta_o[-i], sub$se_o[-i])
          p_i <- 2 * pnorm(-abs(iv_i$b / iv_i$se))
          loo[[length(loo) + 1]] <- list(
            exposure = paste0(src, ":", g), outcome = oc,
            removed_snp = sub$rsid[i], nsnp = k - 1,
            b = iv_i$b, se = iv_i$se, pval = p_i)
        }
      }
    }
  }
}
resDT <- rbindlist(res, fill = TRUE)
sensDT <- rbindlist(sens, fill = TRUE)
looDT <- rbindlist(loo, fill = TRUE)

## ---- BH 校正（预注册口径：按结局分层）----
main_method <- function(x) {
  x[method == "Inverse variance weighted" | method == "Wald ratio"]
}
bh_layer <- function(dt, label) {
  out <- list()
  for (oc in unique(dt$outcome)) {
    sub <- dt[outcome == oc]
    sub[, bh_p := p.adjust(pval, method = "BH")]
    out[[length(out) + 1]] <- sub
  }
  rbindlist(out)
}
res_main <- bh_layer(main_method(resDT), "deCODE")
## 合并 UKB-PPP 敏感性层做 pooled BH
pooled <- main_method(resDT)
pooled[, layer := ifelse(grepl("^UKB-PPP", exposure), "UKB-PPP", "deCODE")]
pooled_out <- list()
for (oc in unique(pooled$outcome)) {
  sub <- pooled[outcome == oc]
  sub[, bh_p_pooled := p.adjust(pval, method = "BH")]
  pooled_out[[length(pooled_out) + 1]] <- sub
}
pooledDT <- rbindlist(pooled_out)
resDT <- merge(resDT, pooledDT[, .(exposure, outcome, bh_p_pooled)],
               by = c("exposure", "outcome"), all.x = TRUE)
resDT <- merge(resDT, res_main[, .(exposure, outcome, bh_p)],
               by = c("exposure", "outcome"), all.x = TRUE)

setcolorder(resDT, c("exposure", "outcome", "method", "nsnp", "b", "se", "pval",
                     "OR", "OR_CI_low", "OR_CI_high", "bh_p", "bh_p_pooled",
                     "exposure_id", "outcome_id", "Q", "Q_pval", "F_stat",
                     "intercept", "intercept_se", "intercept_pval",
                     "gene_set_version", "score_version"))

fwrite(resDT, file.path(OUT, "Table_S15e_pQTL_MR_Results.csv"))
fwrite(sensDT, file.path(OUT, "Table_S15e_pQTL_MR_Sensitivity.csv"))
fwrite(looDT, file.path(OUT, "Table_S15e_pQTL_MR_LeaveOneOut.csv"))
fwrite(d, file.path(OUT, "Table_S15e_pQTL_MR_Instruments.csv"))

cat("\n===== M3 pQTL-MR DONE =====\n")
cat("Method rows:", nrow(resDT), "\n")
cat("\n--- 主方法（IVW/Wald）结果 ---\n")
print(resDT[method == "Inverse variance weighted" | method == "Wald ratio",
            .(exposure, outcome, method, nsnp, b = round(b, 4),
              OR = round(OR, 3), pval = signif(pval, 3),
              bh_p = signif(bh_p, 3), bh_p_pooled = signif(bh_p_pooled, 3),
              Qp = ifelse(is.na(Q_pval), NA, signif(Q_pval, 3)))])
