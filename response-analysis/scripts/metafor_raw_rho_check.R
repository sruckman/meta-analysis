# metafor::rma.mv with raw r (no Fisher Z) - exact replication of the framework used by Sánchez-Tójar and D'Amelio (2026)
# Standalone script. Fast. Just two rma.mv calls, one per dataset.

# Set working directory to the meta-analysis repo root before running.
# Example (edit for your machine): setwd("~/meta-analysis")

library(ape)
library(phytools)
library(metafor)

ORIG_CSV <- "Excel Sheets/meta_complete_data2.csv"
CORR_CSV <- "response-analysis/data/meta_complete_data2_corrected.csv"
OUT_DIR  <- "response-analysis/results/"

tree <- read.tree("Excel Sheets/list.nwk")
tree <- root(tree, "Phymactis_clematis")
tree <- force.ultrametric(tree)
for (i in seq_along(tree$edge.length)) {
  if (tree$edge.length[i] == 0) tree$edge.length[i] <- 0.000001
}
tree$tip.label <- gsub("_", " ", tree$tip.label)
phylo_cor <- vcv(tree, corr = TRUE)

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
  meta
}

run_metafor_raw <- function(csv_path, label, fix_rbis = FALSE) {
  cat("\n=== metafor::rma.mv with RAW r (no Fisher Z) on:", label, "===\n")
  meta <- load_clean(csv_path, fix_rbis = fix_rbis)
  meta$es_id <- seq_len(nrow(meta))
  species_in <- intersect(rownames(phylo_cor), as.character(meta$animal))
  P <- phylo_cor[species_in, species_in]
  meta_sub <- meta[as.character(meta$animal) %in% species_in, ]
  mod <- rma.mv(
    yi = rho,
    V = SE_r^2,
    random = list(~ 1 | animal, ~ 1 | Authors, ~ 1 | es_id),
    R = list(animal = P),
    data = meta_sub,
    method = "REML",
    sparse = TRUE
  )
  s2 <- mod$sigma2
  W <- diag(1 / meta_sub$SE_r^2)
  X <- model.matrix(mod)
  P_mat <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  typical_v <- (mod$k - mod$p) / sum(diag(P_mat))
  I2_total <- 100 * sum(s2) / (sum(s2) + typical_v)
  PI <- predict(mod, level = 95)
  data.frame(
    Dataset = label,
    k = mod$k,
    r_mean = sprintf("%.4f", mod$beta[1]),
    r_CI = sprintf("[%.4f, %.4f]", mod$ci.lb, mod$ci.ub),
    PI = sprintf("[%.4f, %.4f]", PI$pi.lb, PI$pi.ub),
    sigma2_phylo = sprintf("%.4f", s2[1]),
    sigma2_study = sprintf("%.4f", s2[2]),
    sigma2_resid = sprintf("%.4f", s2[3]),
    I2_total_pct = sprintf("%.2f", I2_total),
    stringsAsFactors = FALSE
  )
}

orig_raw_metafor <- run_metafor_raw(ORIG_CSV, "Original", fix_rbis = FALSE)
corr_raw_metafor <- run_metafor_raw(CORR_CSV, "Corrected (rbis fixed)", fix_rbis = TRUE)

cat("\n=========== metafor + raw r (Sánchez-Tójar and D'Amelio 2026 setup) ===========\n")
print(rbind(orig_raw_metafor, corr_raw_metafor), row.names = FALSE)

write.csv(rbind(orig_raw_metafor, corr_raw_metafor),
          paste0(OUT_DIR, "results_metafor_raw_r.csv"),
          row.names = FALSE)
cat("\nSaved to:", paste0(OUT_DIR, "results_metafor_raw_r.csv"), "\n")
