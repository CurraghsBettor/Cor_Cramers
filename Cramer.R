rm(list=ls()) ## clear workspace 

CramersV <- function(data) {
  estimate <- as.numeric(chisq.test(data, correct = F)$statistic)
  phi <- sqrt(estimate/sum(data))
  ## estimate biased Cramer's V; Cramer (1946)
  V <- sqrt(phi^2/min(ncol(data)-1, nrow(data)-1))
  fz <- atanh(V) ## transform Cramer's V into Fisher's z using inverse hyperbolic tangent transformation
  se <- 1/sqrt(sum(data)-3)
  Ll_fz <- sum(fz, qnorm(0.025, mean = 0, sd = 1)*se)
  Ul_fz <- sum(fz, qnorm(1-0.025, mean = 0, sd = 1)*se)
  Ll_V <- tanh(Ll_fz) # transform both lower and upper limits into both lower and upper limits for Cramers'V using hyperbolic tangent transformation
  Ul_V <- tanh(Ul_fz)
  ## estimate biased Tschuprow's T; Tschuprow, 1925, 1939
  T_eff <- sqrt(phi^2/sqrt(nrow(data)*ncol(data)))
  fzT <- atanh(T_eff) ## transform Tschuprow's T into Fisher's z using inverse hyperbolic tangent transformation
  seT <- 1/sqrt(sum(data)-3)
  Ll_fzT <- sum(fzT, qnorm(0.025, mean = 0, sd = 1)*seT)
  Ul_fzT <- sum(fzT, qnorm(1-0.025, mean = 0, sd = 1)*seT)
  Ll_T <- tanh(Ll_fzT) # transform both lower and upper limits into both lower and upper limits for Tschuprow's T using hyperbolic tangent transformation
  Ul_T <- tanh(Ul_fzT)
  # unbiased estimator phi squared and cie...
  phi_sq_c <- phi^2 - 1/(sum(data)-1)*(ncol(data)-1)*(nrow(data)-1) ## bias-corrected version of phi^2
  phi_sq_unbiased <- max(0, phi_sq_c)
  r <- nrow(data) - 1/(sum(data)-1)*((nrow(data)-1)^2)
  c <- ncol(data) - 1/(sum(data)-1)*((ncol(data)-1)^2)
  ## estimate unbiased Cramer's V Bergsma (2013)
  V_unbiased <- sqrt(phi_sq_unbiased/min(c-1, r-1))
  fz_u <- atanh(V_unbiased) ## transform Cramer's V (unbiased) into Fisher's z using inververse hyperbolic tangent transformation
  se_u <- 1/sqrt(sum(data)-3)
  Ll_fz_u <- sum(fz_u, qnorm(0.025, mean = 0, sd = 1)*se_u)
  Ul_fz_u <- sum(fz_u, qnorm(1-0.025, mean = 0, sd = 1)*se_u)
  Ll_Vu <- tanh(Ll_fz_u) # transform both lower and upper limits into both lower and upper limits for Cramers'V (unibiased) using hyperbolic tangent transformation
  Ul_Vu <- tanh(Ul_fz_u)
  ## estimate unbiased Tschuprow's T
  T_unbiased <- sqrt(phi_sq_unbiased/sqrt(nrow(data)*ncol(data)))
  fzTu <- atanh(T_unbiased) ## transform Tschuprow's T into Fisher's z using inververse hyperbolic tangent transformation
  seTu <- 1/sqrt(sum(data)-3)
  Ll_fzTu <- sum(fzTu, qnorm(0.025, mean = 0, sd = 1)*seTu)
  Ul_fzTu <- sum(fzTu, qnorm(1-0.025, mean = 0, sd = 1)*seTu)
  Ll_Tu <- tanh(Ll_fzTu) # transform both lower and upper limits into both lower and upper limits for Tschuprow's T using hyperbolic tangent transformation
  Ul_Tu <- tanh(Ul_fzTu)
  
  cat("Cramer's V =", V, "95%CI [", Ll_V, Ul_V, "]\n")
  cat("Tschuprow's T =", T_eff, "95%CI [", Ll_T, Ul_T, "]\n")
  cat("unbiased Cramer's V =", V_unbiased, "95%CI [", Ll_Vu, Ul_Vu, "]\n")
  cat("unbiased Tschuprow's T =", T_unbiased, "95%CI [", Ll_Tu, Ul_Tu, "]\n")
  
}

data <- matrix(c(29, (72-29), 43, (160-43), 13, (161-13), 4, (104), 3, (112-3)),
          ncol = 2, byrow = T); print(data)
fisher.test(data, simulate.p.value = T)
CramersV(data)
