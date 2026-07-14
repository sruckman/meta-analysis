# Response analysis: original vs corrected, all six tests
# 1. Intercept-only meta-analysis (MCMCglmm)
# 2. Publication year moderator (over-time analysis)
# 3. Effect-size source sensitivity (Stat.Test moderator)
# 4. Publication bias - Egger's regression on meta-analytic residuals + trim and fill
# 5. Publication bias - Nakagawa et al. 2022 sqrt(1/n) moderator
# 6. metafor::rma.mv sanity check with phylogenetic correlation matrix (their framework)
#
# Note: the corrected dataset also gets its rbis values recomputed via metafor::escalc()
# to fully address Critique 2 (the parenthesis bug in the original custom rbis function).

# Set working directory to the meta-analysis repo root before running.
# Example (edit for your machine): setwd("~/meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(meta)         # for trim and fill
library(metafor)      # for escalc and rma.mv

# === Settings ===
CHAIN_LENGTH <- "full"   # "short" or "full"

if (CHAIN_LENGTH == "short") {
  N_ITT <- 500000; N_BURN <- 100000; N_THIN <- 200
} else {
  N_ITT <- 5000000; N_BURN <- 2500000; N_THIN <- 1000
}

ORIG_CSV <- "Excel Sheets/meta_complete_data2.csv"
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

# Phylogenetic correlation matrix for metafor
phylo_cor <- vcv(tree, corr = TRUE)

# === Priors (matches your original) ===
prior.ex <- list(
  G = list(
    G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)
  ),
  R = list(V = 1, nu = 0.02)
)

# === Helper: recompute rbis for "mean" rows using metafor::escalc() ===
# Fixes the parenthesis bug Sánchez-Tójar flagged in the original custom function.
apply_rbis_fix <- function(meta) {
  mask <- meta$Stat.Test == "mean" &
          !is.na(meta$mean1) & !is.na(meta$mean2) &
          !is.na(meta$sd1)   & !is.na(meta$sd2) &
          !is.na(meta$n1)    & !is.na(meta$n2)
  if (!any(mask)) return(meta)
  cat("  Applying escalc(measure='RBIS') fix to", sum(mask), "rows\n")
  es <- escalc(
    measure = "RBIS",
    m1i = meta$mean1[mask], m2i = meta$mean2[mask],
    sd1i = meta$sd1[mask],  sd2i = meta$sd2[mask],
    n1i = meta$n1[mask],    n2i = meta$n2[mask]
  )
  # Preserve original sign convention (positive = darker more aggressive)
  # by checking the existing sign in rho and matching it
  new_rho <- as.numeric(es$yi)
  for (j in seq_along(new_rho)) {
    idx <- which(mask)[j]
    if (!is.na(meta$rho[idx]) && sign(meta$rho[idx]) != 0) {
      new_rho[j] <- abs(new_rho[j]) * sign(meta$rho[idx])
    }
  }
  meta$rho[mask]      <- new_rho
  meta$SE_r[mask]     <- sqrt(as.numeric(es$vi))
  # Fisher Z transform with bounds
  rho_b <- pmin(pmax(new_rho, -0.99), 0.99)
  meta$Fisher_Z[mask] <- 0.5 * log((1 + rho_b) / (1 - rho_b))
  meta$SE_Z[mask]     <- 1 / sqrt(meta$Sample.Size[mask] - 3)
  meta
}

# === Helper to load and clean ===
load_clean <- function(csv_path, fix_rbis = FALSE) {
  meta <- read.csv(csv_path, header = TRUE)
  meta$Eu_Pheomelanin <- ifelse(meta$Eu_Pheomelanin == "N/A",
                                meta$Classification, meta$Eu_Pheomelanin)
  names(meta)[5] <- "animal"
  meta <- meta[meta$animal %in% tree$tip.label, ]
  if (any(meta$Classification == "pteridine", na.rm = TRUE)) {
    meta <- meta[meta$Classification != "pteridine", ]
  }
  if (fix_rbis) meta <- apply_rbis_fix(meta)
  meta$Publication.Year <- as.numeric(as.character(meta$Publication.Year))
  meta$Stat.Test <- factor(meta$Stat.Test)
  meta[, c(2, 4:7, 12:23)] <- lapply(meta[, c(2, 4:7, 12:23)], factor)
  meta
}

# === Helper: run one MCMCglmm model with a given formula ===
run_mod <- function(meta, formula_str, label) {
  cat("  Running:", label, "\n")
  fml <- as.formula(formula_str)
  MCMCglmm(
    fml,
    random = ~ animal + Authors + us(SE_Z):units,
    data = meta, pedigree = tree,
    nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
    prior = prior.ex, verbose = FALSE
  )
}

# === Helper: extract intercept summary ===
intercept_summary <- function(mod, label) {
  z_mean <- mean(mod$Sol[, "(Intercept)"])
  z_ci   <- quantile(mod$Sol[, "(Intercept)"], probs = c(0.025, 0.975))
  r_mean <- (exp(2 * z_mean) - 1) / (exp(2 * z_mean) + 1)
  r_lo   <- (exp(2 * z_ci[1]) - 1) / (exp(2 * z_ci[1]) + 1)
  r_hi   <- (exp(2 * z_ci[2]) - 1) / (exp(2 * z_ci[2]) + 1)
  rand <- mod$VCV / apply(mod$VCV, 1, sum)
  rand_means <- apply(rand, 2, mean)
  total_var <- mod$VCV[, "animal"] + mod$VCV[, "Authors"] +
               mod$VCV[, ncol(mod$VCV)]
  pi_lo <- z_mean - 1.96 * sqrt(mean(total_var))
  pi_hi <- z_mean + 1.96 * sqrt(mean(total_var))
  I2_total <- (rand_means[1] + rand_means[2] + rand_means[length(rand_means)]) / sum(rand_means)
  data.frame(
    Dataset = label,
    k       = nrow(mod$X),
    Z_mean  = sprintf("%.4f", z_mean),
    Z_CI    = sprintf("[%.4f, %.4f]", z_ci[1], z_ci[2]),
    r_mean  = sprintf("%.4f", r_mean),
    r_CI    = sprintf("[%.4f, %.4f]", r_lo, r_hi),
    PI_Z    = sprintf("[%.4f, %.4f]", pi_lo, pi_hi),
    I2      = sprintf("%.4f", I2_total),
    DIC     = sprintf("%.2f", mod$DIC),
    stringsAsFactors = FALSE
  )
}

# === Helper: extract moderator coefficients ===
moderator_summary <- function(mod, label) {
  out <- data.frame(Dataset = character(), Term = character(),
                    Estimate = character(), CI = character(),
                    pMCMC = character(), stringsAsFactors = FALSE)
  for (term in colnames(mod$Sol)) {
    est <- mean(mod$Sol[, term])
    ci  <- quantile(mod$Sol[, term], probs = c(0.025, 0.975))
    p   <- 2 * min(mean(mod$Sol[, term] > 0), mean(mod$Sol[, term] < 0))
    out <- rbind(out, data.frame(
      Dataset = label, Term = term,
      Estimate = sprintf("%.4f", est),
      CI = sprintf("[%.4f, %.4f]", ci[1], ci[2]),
      pMCMC = sprintf("%.4f", p),
      stringsAsFactors = FALSE
    ))
  }
  out
}

# === Helper: publication bias via Egger regression on residuals + trim and fill ===
pubbias_old <- function(mod, meta, label) {
  pred_matrix <- predict(mod, interval = "confidence")
  pred <- if (is.matrix(pred_matrix)) pred_matrix[, "fit"] else pred_matrix
  Precision <- 1 / meta$SE_Z
  MR <- meta$Fisher_Z - pred[seq_len(nrow(meta))]
  zMR <- MR * Precision
  egger <- glm(zMR ~ Precision, family = "gaussian")
  e_int <- summary(egger)$coefficients["(Intercept)", ]
  e_slp <- summary(egger)$coefficients["Precision", ]
  tf <- meta::trimfill(MR, meta$SE_Z)
  data.frame(
    Dataset = label,
    Test = c("Egger intercept", "Egger slope (Precision)", "Trim&Fill: k_after"),
    Estimate = c(sprintf("%.4f", e_int["Estimate"]),
                 sprintf("%.4f", e_slp["Estimate"]),
                 sprintf("%d", tf$k)),
    SE_or_p = c(sprintf("SE=%.4f p=%.4f", e_int["Std. Error"], e_int["Pr(>|t|)"]),
                sprintf("SE=%.4f p=%.4f", e_slp["Std. Error"], e_slp["Pr(>|t|)"]),
                sprintf("studies_added=%d", tf$k - sum(!is.na(MR)))),
    stringsAsFactors = FALSE
  )
}

# === Helper: publication bias via Nakagawa et al. 2022 sqrt(1/n) moderator ===
pubbias_new <- function(meta, label) {
  meta$sqrt_inv_n <- sqrt(1 / meta$Sample.Size)
  m <- MCMCglmm(
    Fisher_Z ~ sqrt_inv_n,
    random = ~ animal + Authors + us(SE_Z):units,
    data = meta, pedigree = tree,
    nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
    prior = prior.ex, verbose = FALSE
  )
  out <- moderator_summary(m, label)
  out$Test <- "Nakagawa 2022: Fisher_Z ~ sqrt(1/n)"
  out
}

# === Helper: metafor::rma.mv sanity check with phylogeny ===
# Matches the framework used by Sánchez-Tójar and D'Amelio (2026) so we can compare I² apples-to-apples.
metafor_check <- function(meta, label) {
  cat("  Running: metafor::rma.mv sanity check\n")
  meta$es_id <- seq_len(nrow(meta))
  # subset phylo matrix to species in this dataset
  species_in <- intersect(rownames(phylo_cor), as.character(meta$animal))
  P <- phylo_cor[species_in, species_in]
  meta_sub <- meta[as.character(meta$animal) %in% species_in, ]
  mod <- rma.mv(
    yi = Fisher_Z,
    V = SE_Z^2,
    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
    R = list(animal = P),
    data = meta_sub,
    method = "REML",
    sparse = TRUE
  )
  # Variance components
  s2 <- mod$sigma2
  # Higgins-style I² using sum of sampling variances
  W <- diag(1 / meta_sub$SE_Z^2)
  X <- model.matrix(mod)
  P_mat <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (mod$k - mod$p) / sum(diag(P_mat))
  I2_total <- 100 * sum(s2) / (sum(s2) + typical_v)
  PI <- predict(mod, level = 95)
  data.frame(
    Dataset = label,
    k = mod$k,
    Estimate = sprintf("%.4f", mod$beta[1]),
    CI = sprintf("[%.4f, %.4f]", mod$ci.lb, mod$ci.ub),
    PI = sprintf("[%.4f, %.4f]", PI$pi.lb, PI$pi.ub),
    sigma2_phylo = sprintf("%.4f", s2[1]),
    sigma2_study = sprintf("%.4f", s2[2]),
    sigma2_resid = sprintf("%.4f", s2[3]),
    I2_total_pct = sprintf("%.2f", I2_total),
    stringsAsFactors = FALSE
  )
}

# === Run everything for one dataset ===
run_all <- function(csv_path, label, fix_rbis = FALSE) {
  cat("\n============================\n")
  cat("DATASET:", label, "\n")
  cat("============================\n")
  meta <- load_clean(csv_path, fix_rbis = fix_rbis)

  m_int  <- run_mod(meta, "Fisher_Z ~ 1", "intercept")
  s_int  <- intercept_summary(m_int, label)

  m_year <- run_mod(meta, "Fisher_Z ~ Publication.Year", "year moderator")
  s_year <- moderator_summary(m_year, label)

  m_src  <- run_mod(meta, "Fisher_Z ~ Stat.Test - 1", "source sensitivity")
  s_src  <- moderator_summary(m_src, label)

  s_pb_old <- pubbias_old(m_int, meta, label)
  s_pb_new <- pubbias_new(meta, label)
  s_metafor <- metafor_check(meta, label)

  list(intercept = s_int, year = s_year, source = s_src,
       pb_old = s_pb_old, pb_new = s_pb_new, metafor = s_metafor)
}

# Original: as published (no rbis fix)
orig_res <- run_all(ORIG_CSV, "Original", fix_rbis = FALSE)
# Corrected: all data fixes + escalc rbis recomputation
corr_res <- run_all(CORR_CSV, "Corrected", fix_rbis = TRUE)

# === Print and save ===
cat("\n\n=============================================\n")
cat("           COMPARISON TABLES\n")
cat("=============================================\n")

cat("\n--- 1. Intercept-only ---\n")
print(rbind(orig_res$intercept, corr_res$intercept), row.names = FALSE)

cat("\n--- 2. Year moderator ---\n")
print(rbind(orig_res$year, corr_res$year), row.names = FALSE)

cat("\n--- 3. Effect-size source sensitivity ---\n")
print(rbind(orig_res$source, corr_res$source), row.names = FALSE)

cat("\n--- 4. Publication bias (Egger + Trim/Fill) ---\n")
print(rbind(orig_res$pb_old, corr_res$pb_old), row.names = FALSE)

cat("\n--- 5. Publication bias (Nakagawa 2022 sqrt(1/n) moderator) ---\n")
print(rbind(orig_res$pb_new, corr_res$pb_new), row.names = FALSE)

cat("\n--- 6. metafor::rma.mv sanity check (their framework) ---\n")
print(rbind(orig_res$metafor, corr_res$metafor), row.names = FALSE)

# Save all tables
write.csv(rbind(orig_res$intercept, corr_res$intercept), paste0(OUT_DIR, "results_intercept.csv"), row.names = FALSE)
write.csv(rbind(orig_res$year, corr_res$year),           paste0(OUT_DIR, "results_year.csv"), row.names = FALSE)
write.csv(rbind(orig_res$source, corr_res$source),       paste0(OUT_DIR, "results_source.csv"), row.names = FALSE)
write.csv(rbind(orig_res$pb_old, corr_res$pb_old),       paste0(OUT_DIR, "results_pubbias_old.csv"), row.names = FALSE)
write.csv(rbind(orig_res$pb_new, corr_res$pb_new),       paste0(OUT_DIR, "results_pubbias_new.csv"), row.names = FALSE)
write.csv(rbind(orig_res$metafor, corr_res$metafor),     paste0(OUT_DIR, "results_metafor_check.csv"), row.names = FALSE)
save(orig_res, corr_res, file = paste0(OUT_DIR, "response_models_full.RDATA"))

cat("\nAll tables saved to:", OUT_DIR, "\n")
