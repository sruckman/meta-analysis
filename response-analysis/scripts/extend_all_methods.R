# Extend the four-method comparison to source sensitivity and both pub bias tests.
# Fills in three cells per remaining method (raw rho + MCMCglmm, Fisher Z + metafor,
# raw r + metafor). Fisher Z + MCMCglmm is already covered by the main script.
# Corrected dataset only. Applies escalc rbis fix.
#
# Runtime: MCMCglmm part takes 2-3 hours at full chain. metafor parts are near-instant.

# Set working directory to the meta-analysis repo root before running.
# Example (edit for your machine): setwd("~/meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(meta)      # trim and fill
library(metafor)   # escalc and rma.mv

# === Settings ===
CHAIN_LENGTH <- "full"   # "short" or "full"

if (CHAIN_LENGTH == "short") {
  N_ITT <- 500000; N_BURN <- 100000; N_THIN <- 200
} else {
  N_ITT <- 5000000; N_BURN <- 2500000; N_THIN <- 1000
}

CORR_CSV <- "response-analysis/data/meta_complete_data2_corrected.csv"
OUT_DIR  <- "response-analysis/results/"

# === Phylogeny ===
tree <- read.tree("Excel Sheets/list.nwk")
tree <- root(tree, "Phymactis_clematis")
tree <- force.ultrametric(tree)
for (i in seq_along(tree$edge.length)) {
  if (tree$edge.length[i] == 0) tree$edge.length[i] <- 0.000001
}
tree$tip.label <- gsub("_", " ", tree$tip.label)
phylo_cor <- vcv(tree, corr = TRUE)

# === Priors ===
prior.ex <- list(
  G = list(
    G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)
  ),
  R = list(V = 1, nu = 0.02)
)

# === rbis fix (same as main script) ===
apply_rbis_fix <- function(meta) {
  mask <- meta$Stat.Test == "mean" &
          !is.na(meta$mean1) & !is.na(meta$mean2) &
          !is.na(meta$sd1)   & !is.na(meta$sd2) &
          !is.na(meta$n1)    & !is.na(meta$n2)
  if (!any(mask)) return(meta)
  es <- escalc(measure = "RBIS",
               m1i = meta$mean1[mask], m2i = meta$mean2[mask],
               sd1i = meta$sd1[mask],  sd2i = meta$sd2[mask],
               n1i = meta$n1[mask],    n2i = meta$n2[mask])
  new_rho <- as.numeric(es$yi)
  for (j in seq_along(new_rho)) {
    idx <- which(mask)[j]
    if (!is.na(meta$rho[idx]) && sign(meta$rho[idx]) != 0) {
      new_rho[j] <- abs(new_rho[j]) * sign(meta$rho[idx])
    }
  }
  meta$rho[mask]  <- new_rho
  meta$SE_r[mask] <- sqrt(as.numeric(es$vi))
  rho_b <- pmin(pmax(new_rho, -0.99), 0.99)
  meta$Fisher_Z[mask] <- 0.5 * log((1 + rho_b) / (1 - rho_b))
  meta$SE_Z[mask]     <- 1 / sqrt(meta$Sample.Size[mask] - 3)
  meta
}

load_clean <- function(csv_path, fix_rbis = TRUE) {
  meta <- read.csv(csv_path, header = TRUE)
  meta$Eu_Pheomelanin <- ifelse(meta$Eu_Pheomelanin == "N/A",
                                meta$Classification, meta$Eu_Pheomelanin)
  names(meta)[5] <- "animal"
  meta <- meta[meta$animal %in% tree$tip.label, ]
  if (any(meta$Classification == "pteridine", na.rm = TRUE)) {
    meta <- meta[meta$Classification != "pteridine", ]
  }
  if (fix_rbis) meta <- apply_rbis_fix(meta)
  meta$Stat.Test <- factor(meta$Stat.Test)
  meta[, c(2, 4:7, 12:23)] <- lapply(meta[, c(2, 4:7, 12:23)], factor)
  meta$SE_rho    <- sqrt((1 - meta$rho^2)^2 / (meta$Sample.Size - 1))
  meta$es_id     <- seq_len(nrow(meta))
  meta$sqrt_inv_n <- sqrt(1 / meta$Sample.Size)
  meta
}

# ==========================================================================
# METHOD 2: Raw rho + MCMCglmm
# ==========================================================================
run_method2 <- function() {
  cat("\n\n===== METHOD 2: Raw rho + MCMCglmm =====\n")
  m <- load_clean(CORR_CSV, fix_rbis = TRUE)

  # Intercept-only (needed for residuals)
  cat("  Intercept-only...\n")
  m2_int <- MCMCglmm(rho ~ 1,
                     random = ~ animal + Authors + us(SE_rho):units,
                     data = m, pedigree = tree,
                     nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
                     prior = prior.ex, verbose = FALSE)

  # Source sensitivity
  cat("  Source sensitivity...\n")
  m2_src <- MCMCglmm(rho ~ Stat.Test - 1,
                     random = ~ animal + Authors + us(SE_rho):units,
                     data = m, pedigree = tree,
                     nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
                     prior = prior.ex, verbose = FALSE)

  # Pub bias new (Nakagawa 2022 sqrt(1/n) moderator)
  cat("  Pub bias new (sqrt(1/n) moderator)...\n")
  m2_pbn <- MCMCglmm(rho ~ sqrt_inv_n,
                     random = ~ animal + Authors + us(SE_rho):units,
                     data = m, pedigree = tree,
                     nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
                     prior = prior.ex, verbose = FALSE)

  # Pub bias old (Egger on residuals + trim/fill)
  pred <- predict(m2_int, interval = "confidence")
  pred_vec <- if (is.matrix(pred)) pred[, "fit"] else pred
  Precision <- 1 / m$SE_rho
  MR <- m$rho - pred_vec[seq_len(nrow(m))]
  zMR <- MR * Precision
  egger <- glm(zMR ~ Precision, family = "gaussian")
  e_int <- summary(egger)$coefficients["(Intercept)", ]
  e_slp <- summary(egger)$coefficients["Precision", ]
  tf <- meta::trimfill(MR, m$SE_rho)

  # Extract summaries
  fmt_mod <- function(mod, method_label) {
    out <- data.frame(Method = character(), Term = character(),
                      Estimate = character(), CI = character(),
                      pMCMC = character(), stringsAsFactors = FALSE)
    for (term in colnames(mod$Sol)) {
      est <- mean(mod$Sol[, term])
      ci  <- quantile(mod$Sol[, term], probs = c(0.025, 0.975))
      p   <- 2 * min(mean(mod$Sol[, term] > 0), mean(mod$Sol[, term] < 0))
      out <- rbind(out, data.frame(
        Method = method_label, Term = term,
        Estimate = sprintf("%.4f", est),
        CI = sprintf("[%.4f, %.4f]", ci[1], ci[2]),
        pMCMC = sprintf("%.4f", p),
        stringsAsFactors = FALSE
      ))
    }
    out
  }

  src_out <- fmt_mod(m2_src, "Raw rho + MCMCglmm")
  pbn_out <- fmt_mod(m2_pbn, "Raw rho + MCMCglmm")

  pbo_out <- data.frame(
    Method = "Raw rho + MCMCglmm",
    Test = c("Egger intercept", "Egger slope (Precision)", "Trim&Fill k_after"),
    Estimate = c(sprintf("%.4f", e_int["Estimate"]),
                 sprintf("%.4f", e_slp["Estimate"]),
                 sprintf("%d", tf$k)),
    SE_or_p = c(sprintf("SE=%.4f p=%.4f", e_int["Std. Error"], e_int["Pr(>|t|)"]),
                sprintf("SE=%.4f p=%.4f", e_slp["Std. Error"], e_slp["Pr(>|t|)"]),
                sprintf("added=%d", tf$k - sum(!is.na(MR)))),
    stringsAsFactors = FALSE
  )

  list(src = src_out, pb_old = pbo_out, pb_new = pbn_out)
}

# ==========================================================================
# METAFOR helper (methods 3 and 4)
# ==========================================================================
run_metafor <- function(yi_col, V_col, se_col, method_label) {
  cat("\n\n===== METHOD:", method_label, "=====\n")
  m <- load_clean(CORR_CSV, fix_rbis = TRUE)
  species_in <- intersect(rownames(phylo_cor), as.character(m$animal))
  P <- phylo_cor[species_in, species_in]
  m <- m[as.character(m$animal) %in% species_in, ]

  yi <- m[[yi_col]]
  V  <- m[[V_col]]^2

  # Intercept-only
  cat("  Intercept-only...\n")
  mod_int <- rma.mv(yi = yi, V = V,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m,
                    method = "REML", sparse = TRUE)

  # Source sensitivity
  cat("  Source sensitivity...\n")
  mod_src <- rma.mv(yi = yi, V = V,
                    mods = ~ Stat.Test - 1,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m,
                    method = "REML", sparse = TRUE)

  # Pub bias new (sqrt(1/n) moderator)
  cat("  Pub bias new (sqrt(1/n) moderator)...\n")
  mod_pbn <- rma.mv(yi = yi, V = V,
                    mods = ~ sqrt_inv_n,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m,
                    method = "REML", sparse = TRUE)

  # Pub bias old (Egger on residuals + trim/fill)
  resid_val <- residuals(mod_int)
  Precision <- 1 / m[[se_col]]
  zMR <- resid_val * Precision
  egger <- glm(zMR ~ Precision, family = "gaussian")
  e_int <- summary(egger)$coefficients["(Intercept)", ]
  e_slp <- summary(egger)$coefficients["Precision", ]
  tf <- meta::trimfill(resid_val, m[[se_col]])

  # Format
  fmt_rma <- function(mod, method_label) {
    out <- data.frame(Method = character(), Term = character(),
                      Estimate = character(), CI = character(),
                      pval = character(), stringsAsFactors = FALSE)
    coefs <- coef(summary(mod))
    for (i in seq_len(nrow(coefs))) {
      term <- rownames(coefs)[i]
      out <- rbind(out, data.frame(
        Method = method_label, Term = term,
        Estimate = sprintf("%.4f", coefs[i, "estimate"]),
        CI = sprintf("[%.4f, %.4f]", coefs[i, "ci.lb"], coefs[i, "ci.ub"]),
        pval = sprintf("%.4f", coefs[i, "pval"]),
        stringsAsFactors = FALSE
      ))
    }
    out
  }

  src_out <- fmt_rma(mod_src, method_label)
  pbn_out <- fmt_rma(mod_pbn, method_label)

  pbo_out <- data.frame(
    Method = method_label,
    Test = c("Egger intercept", "Egger slope (Precision)", "Trim&Fill k_after"),
    Estimate = c(sprintf("%.4f", e_int["Estimate"]),
                 sprintf("%.4f", e_slp["Estimate"]),
                 sprintf("%d", tf$k)),
    SE_or_p = c(sprintf("SE=%.4f p=%.4f", e_int["Std. Error"], e_int["Pr(>|t|)"]),
                sprintf("SE=%.4f p=%.4f", e_slp["Std. Error"], e_slp["Pr(>|t|)"]),
                sprintf("added=%d", tf$k - sum(!is.na(resid_val)))),
    stringsAsFactors = FALSE
  )

  list(src = src_out, pb_old = pbo_out, pb_new = pbn_out)
}

# ==========================================================================
# RUN ALL
# ==========================================================================
m2 <- run_method2()

# Redefined helper below with cleaner signature (yi_col, se_col, method_label)
run_metafor <- function(yi_col, se_col, method_label) {
  cat("\n\n===== METHOD:", method_label, "=====\n")
  m <- load_clean(CORR_CSV, fix_rbis = TRUE)
  species_in <- intersect(rownames(phylo_cor), as.character(m$animal))
  P <- phylo_cor[species_in, species_in]
  m <- m[as.character(m$animal) %in% species_in, ]
  yi <- m[[yi_col]]
  V  <- m[[se_col]]^2
  mod_int <- rma.mv(yi = yi, V = V,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m, method = "REML", sparse = TRUE)
  mod_src <- rma.mv(yi = yi, V = V, mods = ~ Stat.Test - 1,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m, method = "REML", sparse = TRUE)
  mod_pbn <- rma.mv(yi = yi, V = V, mods = ~ sqrt_inv_n,
                    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
                    R = list(animal = P), data = m, method = "REML", sparse = TRUE)
  resid_val <- residuals(mod_int)
  Precision <- 1 / m[[se_col]]
  zMR <- resid_val * Precision
  egger <- glm(zMR ~ Precision, family = "gaussian")
  e_int <- summary(egger)$coefficients["(Intercept)", ]
  e_slp <- summary(egger)$coefficients["Precision", ]
  tf <- meta::trimfill(resid_val, m[[se_col]])
  fmt_rma <- function(mod, method_label) {
    out <- data.frame(Method = character(), Term = character(),
                      Estimate = character(), CI = character(),
                      pval = character(), stringsAsFactors = FALSE)
    coefs <- coef(summary(mod))
    for (i in seq_len(nrow(coefs))) {
      term <- rownames(coefs)[i]
      out <- rbind(out, data.frame(
        Method = method_label, Term = term,
        Estimate = sprintf("%.4f", coefs[i, "estimate"]),
        CI = sprintf("[%.4f, %.4f]", coefs[i, "ci.lb"], coefs[i, "ci.ub"]),
        pval = sprintf("%.4f", coefs[i, "pval"]),
        stringsAsFactors = FALSE
      ))
    }
    out
  }
  list(
    src = fmt_rma(mod_src, method_label),
    pb_new = fmt_rma(mod_pbn, method_label),
    pb_old = data.frame(
      Method = method_label,
      Test = c("Egger intercept", "Egger slope (Precision)", "Trim&Fill k_after"),
      Estimate = c(sprintf("%.4f", e_int["Estimate"]),
                   sprintf("%.4f", e_slp["Estimate"]),
                   sprintf("%d", tf$k)),
      SE_or_p = c(sprintf("SE=%.4f p=%.4f", e_int["Std. Error"], e_int["Pr(>|t|)"]),
                  sprintf("SE=%.4f p=%.4f", e_slp["Std. Error"], e_slp["Pr(>|t|)"]),
                  sprintf("added=%d", tf$k - sum(!is.na(resid_val)))),
      stringsAsFactors = FALSE
    )
  )
}

m3 <- run_metafor("Fisher_Z", "SE_Z", "Fisher Z + metafor")
m4 <- run_metafor("rho",      "SE_r", "Raw r + metafor")

# ==========================================================================
# COMBINE AND SAVE
# ==========================================================================
# MCMCglmm uses pMCMC and metafor uses pval. Rename to a common column
# so rbind works. Remember when reading: p is Bayesian pMCMC for method 2
# and frequentist p for methods 3 and 4.
names(m2$src)[names(m2$src) == "pMCMC"]      <- "p"
names(m3$src)[names(m3$src) == "pval"]       <- "p"
names(m4$src)[names(m4$src) == "pval"]       <- "p"
names(m2$pb_new)[names(m2$pb_new) == "pMCMC"] <- "p"
names(m3$pb_new)[names(m3$pb_new) == "pval"]  <- "p"
names(m4$pb_new)[names(m4$pb_new) == "pval"]  <- "p"

src_all <- rbind(m2$src, m3$src, m4$src)
pbo_all <- rbind(m2$pb_old, m3$pb_old, m4$pb_old)
pbn_all <- rbind(m2$pb_new, m3$pb_new, m4$pb_new)

cat("\n\n=========== SOURCE SENSITIVITY (3 methods) ===========\n")
print(src_all, row.names = FALSE)

cat("\n=========== PUB BIAS - YOUR WAY (3 methods) ===========\n")
print(pbo_all, row.names = FALSE)

cat("\n=========== PUB BIAS - THEIR WAY (3 methods) ===========\n")
print(pbn_all, row.names = FALSE)

write.csv(src_all, paste0(OUT_DIR, "results_source_extended.csv"), row.names = FALSE)
write.csv(pbo_all, paste0(OUT_DIR, "results_pubbias_old_extended.csv"), row.names = FALSE)
write.csv(pbn_all, paste0(OUT_DIR, "results_pubbias_new_extended.csv"), row.names = FALSE)

cat("\nSaved to:", OUT_DIR, "\n")
cat("Combine with results_source.csv, results_pubbias_old.csv, results_pubbias_new.csv\n")
cat("from the main run to see all four methods side by side.\n")
