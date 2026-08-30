## TWMR rerun: 15 gene-expression exposures -> FinnGen R10 ARDS outcome
## Output: Table_S15b_MR_FinnGenARDS_{Results,Instruments,Sensitivity}.csv
## Exposure beta in SD units (eQTLGen normalized) -> b is log-odds ARDS per 1-SD expression increase.
library(data.table)
set.seed(42)
OUT  <- "E:/SCI/SCI论文1黄裕荣_Mitoxyperilysis_ARDS/01_RESULTS_TABLES"
WORK <- "E:/SCI/_mr_work"

## ---- 1. instruments (gene expression only; drop CRP/Ferritin biomarker) ----
instr <- fread(file.path(OUT, "Table_S15_MR_Instruments.csv"))
instr <- instr[exposure_type == "gene_expression"]
setnames(instr, c("effect_allele","other_allele","beta","se","pval","eaf"),
         c("EA_e","OA_e","beta_e","se_e","pval_e","eaf_e"))
instr[, SNP := as.character(SNP)]

## ---- 2. FinnGen outcome (local) ----
fg <- fread(file.path(WORK, "finngen_instruments.csv"))
fg[, SNP := as.character(rsids)]
fg <- fg[!duplicated(SNP)]   # drop multi-allelic dup rsIDs
fg <- fg[, .(SNP, REF_o=ref, ALT_o=alt, beta_o=beta, se_o=sebeta, eaf_o=af_alt, pval_o=pval)]

## ---- 3. merge ----
d <- merge(instr, fg, by="SNP")

## ---- 4. harmonize alleles (outcome aligned to exposure effect allele) ----
comp <- function(x){ m <- c("A"="T","C"="G","G"="C","T"="A"); unname(m[toupper(x)]) }
n <- nrow(d); keep <- logical(n); flip <- logical(n)
for(i in seq_len(n)){
  e1 <- toupper(d$EA_e[i]); e2 <- toupper(d$OA_e[i])
  o1 <- toupper(d$ALT_o[i]); o2 <- toupper(d$REF_o[i])   # outcome effect=ALT
  f <- FALSE; k <- TRUE
  if(!is.na(e1) && !is.na(o1)){
    if(e1==o1 && e2==o2){ }
    else if(e1==o2 && e2==o1){ f <- TRUE }
    else {
      o1c <- comp(o1); o2c <- comp(o2)
      if(e1==o1c && e2==o2c){ }
      else if(e1==o2c && e2==o1c){ f <- TRUE }
      else k <- FALSE
    }
  } else k <- FALSE
  pal <- (e1=="A"&&e2=="T")||(e1=="T"&&e2=="A")||(e1=="C"&&e2=="G")||(e1=="G"&&e2=="C")
  if(k && pal){
    eo <- if(f) 1-d$eaf_o[i] else d$eaf_o[i]
    if(is.na(d$eaf_e[i]) || is.na(eo)){ k <- FALSE }
    else if(abs(d$eaf_e[i]-eo) <= 0.2){ }
    else if(abs(d$eaf_e[i]-(1-eo)) <= 0.2){ f <- !f }
    else k <- FALSE
  }
  keep[i] <- k; flip[i] <- f
}
d[, c("keep","flip") := .(keep, flip)]
dh <- d[keep==TRUE]
dh[, beta_o_h := ifelse(flip, -beta_o, beta_o)]
dh[, eaf_o_h := ifelse(flip, 1-eaf_o, eaf_o)]
dropped <- d[keep==FALSE, .(SNP, exposure, EA_e, OA_e, ALT_o, REF_o, eaf_e, eaf_o)]

## ---- 5. MR methods ----
mr_ivw <- function(b_e,b_o,se_o){
  w <- 1/se_o^2; b <- sum(w*b_e*b_o)/sum(w*b_e^2); se <- sqrt(1/sum(w*b_e^2))
  Q <- sum(w*(b_o-b_e*b)^2); list(b=b,se=se,Q=Q)
}
mr_egger <- function(b_e,b_o,se_o){
  w <- 1/se_o^2; fit <- lm(b_o ~ b_e, weights=w); s <- summary(fit)$coefficients
  list(b=s["b_e","Estimate"], se=s["b_e","Std. Error"], p=s["b_e","Pr(>|t|)"],
       ic=s["(Intercept)","Estimate"], ic_se=s["(Intercept)","Std. Error"], ic_p=s["(Intercept)","Pr(>|t|)"])
}
mr_wald <- function(b_e,b_o,se_o,se_e){
  b <- b_o/b_e; se <- abs(b)*sqrt((se_o/b_o)^2+(se_e/b_e)^2); list(b=b,se=se)
}
wmedian <- function(r,ww){ o<-order(r); cw<-cumsum(ww[o]/sum(ww)); r[o][min(which(cw>=0.5))] }
mr_wmed <- function(b_e,b_o,se_o){
  r <- b_o/b_e; w <- 1/(se_o/b_e)^2; b <- wmedian(r,w)
  nn <- length(b_e); B <- 1000; boot <- numeric(B)
  for(j in seq_len(B)){ s<-sample(nn,replace=TRUE); boot[j]<-wmedian(b_o[s]/b_e[s], 1/(se_o[s]/b_e[s])^2) }
  list(b=b, se=sd(boot))
}
mk <- function(g,method,nsnp,b,se,Q=NA,Qp=NA,Fm,ic=NA,ic_se=NA,ic_p=NA,exp_id){
  p <- 2*pnorm(-abs(b/se))
  list(exposure=g, outcome="FinnGen_R10_ARDS", method=method, nsnp=nsnp,
       b=b, se=se, pval=p, OR=exp(b), OR_CI_low=exp(b-1.96*se), OR_CI_high=exp(b+1.96*se),
       exposure_id=exp_id, outcome_id="finngen_R10_J10_ARDS",
       Q=Q, Q_pval=Qp, F_stat=Fm, intercept=ic, intercept_se=ic_se, intercept_pval=ic_p)
}

res <- list(); sens <- list()
for(g in sort(unique(dh$exposure))){
  sub <- dh[exposure==g]; k <- nrow(sub)
  if(k==0) next
  exp_id <- sub$exposure_id[1]
  Fm <- mean((sub$beta_e/sub$se_e)^2)
  sens[[length(sens)+1]] <- list(exposure=g, n_instruments=nrow(instr[exposure==g]),
                                 n_after_harmonization=k, n_dropped=nrow(dropped[exposure==g]),
                                 mean_F=Fm)
  if(k==1){
    w <- mr_wald(sub$beta_e, sub$beta_o_h, sub$se_o, sub$se_e)
    res[[length(res)+1]] <- mk(g,"Wald ratio",1,w$b,w$se,Fm=Fm,exp_id=exp_id)
  } else {
    iv <- mr_ivw(sub$beta_e, sub$beta_o_h, sub$se_o)
    Qp <- pchisq(iv$Q, k-1, lower.tail=FALSE)
    res[[length(res)+1]] <- mk(g,"Inverse variance weighted",k,iv$b,iv$se,Q=iv$Q,Qp=Qp,Fm=Fm,exp_id=exp_id)
    if(k>=3){
      eg <- mr_egger(sub$beta_e, sub$beta_o_h, sub$se_o)
      res[[length(res)+1]] <- mk(g,"MR Egger",k,eg$b,eg$se,Fm=Fm,ic=eg$ic,ic_se=eg$ic_se,ic_p=eg$ic_p,exp_id=exp_id)
      wm <- mr_wmed(sub$beta_e, sub$beta_o_h, sub$se_o)
      res[[length(res)+1]] <- mk(g,"Weighted median",k,wm$b,wm$se,Fm=Fm,exp_id=exp_id)
    }
  }
}
resDT <- rbindlist(res, fill=TRUE)
setcolorder(resDT, c("exposure","outcome","method","nsnp","b","se","pval","OR","OR_CI_low","OR_CI_high",
                     "exposure_id","outcome_id","Q","Q_pval","F_stat","intercept","intercept_se","intercept_pval"))
sensDT <- rbindlist(sens, fill=TRUE)
instr_out <- dh[, .(SNP, exposure, exposure_id, EA_e, OA_e, beta_e, se_e, eaf_e, ALT_o, REF_o, beta_o_h, se_o, eaf_o_h, pval_o, harmonize=ifelse(flip,"flip","same"))]

fwrite(resDT,   file.path(OUT, "Table_S15b_MR_FinnGenARDS_Results.csv"))
fwrite(instr_out, file.path(OUT, "Table_S15b_MR_FinnGenARDS_Instruments.csv"))
fwrite(sensDT,  file.path(OUT, "Table_S15b_MR_FinnGenARDS_Sensitivity.csv"))
fwrite(dropped, file.path(WORK, "dropped_instruments.csv"))

cat("\n=== TWMR FinnGen ARDS: DONE ===\n")
cat("Genes analyzed:", length(unique(resDT$exposure)), "\n")
cat("Total method-rows:", nrow(resDT), "\n")
cat("\n--- IVW results (log-odds ARDS per 1-SD expression) ---\n")
print(resDT[method=="Inverse variance weighted" | method=="Wald ratio",
            .(exposure, nsnp, b=round(b,4), OR=round(OR,3), pval=signif(pval,3), F=round(F_stat,1))])
cat("\nDropped (incompatible/ambiguous) SNPs:", nrow(dropped), "\n")
cat("Sensitivity:\n"); print(sensDT)
