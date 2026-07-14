# Forest plot for the four-method comparison
# Builds two versions:
#   (A) Unified r-scale plot (Fisher Z values back-transformed)
#   (B) Two-panel plot showing native scales (Fisher Z vs rho/r)

library(ggplot2)
library(patchwork)  # install.packages("patchwork", repos = "https://cloud.r-project.org") if needed

OUT_DIR <- "response-analysis/figures/"

# Native-scale estimates from the corrected dataset runs
results <- data.frame(
  Method = c("Fisher Z + MCMCglmm",
             "Raw rho + MCMCglmm",
             "Fisher Z + metafor",
             "Raw r + metafor (Sánchez-Tójar and D'Amelio 2026 setup)"),
  estimate_native = c(0.2702, 0.2257, 0.2721, 0.2224),
  lower_native    = c(-0.0556, 0.0182, 0.0043, 0.0519),
  upper_native    = c(0.6066, 0.4503, 0.5398, 0.3928),
  scale = c("Fisher Z", "r", "Fisher Z", "r")
)

# Order so the bottom row is the "Sánchez-Tójar setup" highlight
results$Method <- factor(results$Method, levels = rev(results$Method))

# Convert all to r for unified plot
back_to_r <- function(x, scale) ifelse(scale == "Fisher Z", tanh(x), x)
results$estimate_r <- back_to_r(results$estimate_native, results$scale)
results$lower_r    <- back_to_r(results$lower_native,    results$scale)
results$upper_r    <- back_to_r(results$upper_native,    results$scale)

# === Plot A: Unified r scale ===
pA <- ggplot(results, aes(x = estimate_r, y = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower_r, xmax = upper_r), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(limits = c(-0.1, 0.6), breaks = seq(-0.1, 0.6, 0.1)) +
  labs(x = "Effect size r (95% CI)",
       y = NULL,
       title = "Four-method comparison, corrected dataset",
       subtitle = "All estimates on the r (correlation) scale; Fisher Z values back-transformed via tanh()") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(paste0(OUT_DIR, "forest_plot_unified_r.png"), pA, width = 8, height = 3.5, dpi = 300)
ggsave(paste0(OUT_DIR, "forest_plot_unified_r.pdf"), pA, width = 8, height = 3.5)

# === Plot B: Two-panel native scales ===
pB_top <- ggplot(subset(results, scale == "Fisher Z"),
                 aes(x = estimate_native, y = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower_native, xmax = upper_native), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(limits = c(-0.1, 0.7)) +
  labs(x = "Fisher Z (95% CI)", y = NULL, title = "Fisher Z scale") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

pB_bot <- ggplot(subset(results, scale == "r"),
                 aes(x = estimate_native, y = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower_native, xmax = upper_native), height = 0.18, linewidth = 0.6) +
  geom_point(size = 3.5, color = "black") +
  scale_x_continuous(limits = c(-0.1, 0.7)) +
  labs(x = "r (95% CI)", y = NULL, title = "r scale (raw)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), panel.grid.minor = element_blank())

pB <- pB_top / pB_bot +
  plot_annotation(title = "Four-method comparison, corrected dataset",
                  subtitle = "Native estimation scales preserved",
                  theme = theme(plot.title = element_text(face = "bold")))

ggsave(paste0(OUT_DIR, "forest_plot_native_scales.png"), pB, width = 8, height = 5, dpi = 300)
ggsave(paste0(OUT_DIR, "forest_plot_native_scales.pdf"), pB, width = 8, height = 5)

# Print summary to console so you can sanity-check the back-transformations
cat("\n=== Back-transformed r values used in the unified plot ===\n")
print(results[, c("Method", "scale", "estimate_native", "estimate_r", "lower_r", "upper_r")],
      row.names = FALSE)

cat("\nFiles saved:\n",
    OUT_DIR, "forest_plot_unified_r.png and .pdf\n",
    OUT_DIR, "forest_plot_native_scales.png and .pdf\n", sep = "")
