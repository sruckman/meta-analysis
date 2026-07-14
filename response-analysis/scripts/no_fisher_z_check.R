# Sensitivity check: MCMCglmm WITHOUT Fisher Z transformation
# Runs the same intercept-only model on both datasets using raw correlation
# (rho) instead of Fisher Z. Purpose: demonstrate the bounded-distribution
# problem in Bayesian inference and support the defense of Fisher Z.
#
# Expect: poor mixing, wide credible intervals, possibly biased point estimates
# compared to the Fisher Z version.

# Set working directory to the meta-analysis repo root before running.
# Example (edit for your machine): setwd("~/meta-analysis")

library(ape)
library(MCMCglmm)
library(phytools)
library(metafor)

# === Settings ===
# Full chain to match the published analysis
N_ITT  <- 5000000
N_BURN <- 2500000
N_THIN <- 1000

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

# === Priors ===
prior.ex <- list(
  G = list(
    G1 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G2 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5),
    G3 = list(V = 1, nu = 0.02, alpha.mu = 0.4, alpha.V = 0.5)
  ),
  R = list(V = 1, nu = 0.02)
)

# === rbis fix (same as in main script) ===
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
  meta
}

# === Loader ===
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
  meta$Stat.Test <- factor(meta$Stat.Test)
  meta[, c(2, 4:7, 12:23)] <- lapply(meta[, c(2, 4:7, 12:23)], factor)
  meta
}

# === Run model on raw rho instead of Fisher_Z ===
run_raw_rho <- function(csv_path, label, fix_rbis = FALSE) {
  cat("\n=== Running MCMCglmm with raw rho (no Fisher Z) on:", label, "===\n")
  meta <- load_clean(csv_path, fix_rbis = fix_rbis)
  # SE for raw rho: sqrt((1-rho^2)^2 / (n-1))
  meta$SE_rho <- sqrt((1 - meta$rho^2)^2 / (meta$Sample.Size - 1))
  mod <- MCMCglmm(
    rho ~ 1,
    random = ~ animal + Authors + us(SE_rho):units,
    data = meta, pedigree = tree,
    nitt = N_ITT, thin = N_THIN, burnin = N_BURN,
    prior = prior.ex, verbose = FALSE
  )
  rho_mean <- mean(mod$Sol[, "(Intercept)"])
  rho_ci   <- quantile(mod$Sol[, "(Intercept)"], probs = c(0.025, 0.975))
  rand     <- mod$VCV / apply(mod$VCV, 1, sum)
  rand_m   <- apply(rand, 2, mean)
  I2_total <- (rand_m[1] + rand_m[2] + rand_m[length(rand_m)]) / sum(rand_m)
  # Diagnostics: autocorrelation and effective sample size
  auto1    <- autocorr(mod$Sol)[2, 1, 1]  # lag-1 autocorrelation of intercept
  ess      <- effectiveSize(mod$Sol)["(Intercept)"]
  data.frame(
    Dataset = label,
    k       = nrow(mod$X),
    rho_mean = sprintf("%.4f", rho_mean),
    rho_CI  = sprintf("[%.4f, %.4f]", rho_ci[1], rho_ci[2]),
    I2      = sprintf("%.4f", I2_total),
    DIC     = sprintf("%.2f", mod$DIC),
    lag1_autocor = sprintf("%.4f", auto1),
    eff_sample   = sprintf("%.1f", ess),
    note    = ifelse(abs(auto1) > 0.1, "high autocorr (poor mixing)", "OK"),
    stringsAsFactors = FALSE
  )
}

orig_raw <- run_raw_rho(ORIG_CSV, "Original", fix_rbis = FALSE)
corr_raw <- run_raw_rho(CORR_CSV, "Corrected (rbis fixed)", fix_rbis = TRUE)

cat("\n=============== NO FISHER Z RESULTS ===============\n")
print(rbind(orig_raw, corr_raw), row.names = FALSE)

write.csv(rbind(orig_raw, corr_raw),
          paste0(OUT_DIR, "results_no_fisher_z.csv"),
          row.names = FALSE)
cat("\nSaved to:", paste0(OUT_DIR, "results_no_fisher_z.csv"), "\n")
