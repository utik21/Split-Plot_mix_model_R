# =============================================================================
# Getting the Split-Plot Right: Oat Variety x Nitrogen Field Trial
# -----------------------------------------------------------------------------
# Self-initiated portfolio project using public data.
# Author : Utik Tri Wulan Cahya, PhD
# Data   : Yates (1935) oats split-plot experiment, distributed as `oats`
#          in the R package MASS (Venables & Ripley).
#
# Purpose: show the correct analysis of a split-plot field trial with a
#          linear mixed model, and how the common mistake of analysing it as
#          an ordinary randomized complete block design (RCBD) changes the
#          conclusions.
#
# Run from the project root:  Rscript R/splitplot_analysis.R
# Outputs are written to output/figures, output/tables and output/model_output.txt
# =============================================================================

suppressPackageStartupMessages({
  library(MASS)        # data
  library(lme4)        # mixed models
  library(lmerTest)    # F tests with Kenward-Roger / Satterthwaite df
  library(emmeans)     # estimated marginal means and contrasts
  library(multcomp)    # compact letter display
  library(multcompView)
  library(car)         # Levene's test
  library(ggplot2)
  library(patchwork)
})

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("data",           showWarnings = FALSE)

log_file <- "output/model_output.txt"
sink(log_file, split = TRUE)
section <- function(title) cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")

# -----------------------------------------------------------------------------
# 1. Data and design
# -----------------------------------------------------------------------------
section("1. DATA AND DESIGN")

data(oats, package = "MASS")
oats <- data.frame(
  Block    = oats$B,
  Variety  = factor(oats$V, labels = c("Golden Rain", "Marvellous", "Victory")),
  Nitrogen = factor(oats$N, labels = c("0", "0.2", "0.4", "0.6")),
  N_rate   = as.numeric(sub("cwt", "", as.character(oats$N))),
  Yield    = oats$Y
)
# Each whole plot is one variety within one block; nitrogen is applied to
# the four subplots inside it.
oats$WholePlot <- interaction(oats$Block, oats$Variety, drop = TRUE)
write.csv(oats, "data/oats.csv", row.names = FALSE)

cat("Observations:", nrow(oats), "\n")
cat("Blocks:", nlevels(oats$Block), "| Varieties (whole-plot factor):", nlevels(oats$Variety),
    "| Nitrogen rates (subplot factor):", nlevels(oats$Nitrogen), "\n")
cat("Whole plots:", nlevels(oats$WholePlot), "(", nlevels(oats$Block), "blocks x",
    nlevels(oats$Variety), "varieties ), each split into 4 nitrogen subplots\n\n")
cat("Replication check (subplots per block x variety):\n")
print(table(oats$Block, oats$Variety))
cat("\nNitrogen is in cwt per acre; yield is in the original units of the dataset",
    "(quarter-pounds per subplot).\n")

# -----------------------------------------------------------------------------
# 2. The common mistake: analysing the split-plot as an RCBD
# -----------------------------------------------------------------------------
section("2. INCORRECT ANALYSIS - treated as an ordinary RCBD")
cat("All effects are tested against one pooled residual error.\n\n")
m_wrong <- aov(Yield ~ Block + Variety * Nitrogen, data = oats)
print(summary(m_wrong))

# -----------------------------------------------------------------------------
# 3. Correct analysis: split-plot ANOVA (two error strata)
# -----------------------------------------------------------------------------
section("3. CORRECT ANALYSIS - split-plot ANOVA with two error strata")
cat("Variety is tested against the whole-plot error (Block x Variety);\n",
    "Nitrogen and Variety x Nitrogen against the subplot error.\n\n", sep = "")
m_split <- aov(Yield ~ Variety * Nitrogen + Error(Block / Variety), data = oats)
print(summary(m_split))

# -----------------------------------------------------------------------------
# 4. Same model as a linear mixed model (the modern, general approach)
# -----------------------------------------------------------------------------
section("4. LINEAR MIXED MODEL - random block and whole-plot effects")
m_mixed <- lmer(Yield ~ Variety * Nitrogen + (1 | Block) + (1 | Block:Variety),
                data = oats, REML = TRUE)
cat("Model: Yield ~ Variety * Nitrogen + (1 | Block) + (1 | Block:Variety)\n\n")
cat("Type III F tests, Kenward-Roger degrees of freedom:\n")
an_mixed <- anova(m_mixed, ddf = "Kenward-Roger")
print(an_mixed)

vc <- as.data.frame(VarCorr(m_mixed))
vc <- data.frame(Component = c("Block:Variety (whole-plot error)", "Block", "Residual (subplot error)"),
                 Variance  = round(vc$vcov[match(c("Block:Variety", "Block", "Residual"), vc$grp)], 1))
vc$Percent <- round(100 * vc$Variance / sum(vc$Variance), 1)
cat("\nVariance components:\n"); print(vc, row.names = FALSE)
write.csv(vc, "output/tables/variance_components.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Side-by-side comparison of the two analyses
# -----------------------------------------------------------------------------
section("5. COMPARISON - what the wrong analysis gets wrong")

get_row <- function(tab, term) {
  rn <- trimws(rownames(tab))
  tab[rn == term, , drop = FALSE]
}
wrong_tab <- summary(m_wrong)[[1]]
split_wp  <- summary(m_split)[["Error: Block:Variety"]][[1]]
split_sp  <- summary(m_split)[["Error: Within"]][[1]]
wrong_res_df <- get_row(wrong_tab, "Residuals")$Df

terms <- c("Variety", "Nitrogen", "Variety:Nitrogen")
comparison <- do.call(rbind, lapply(terms, function(tm) {
  w <- get_row(wrong_tab, tm)
  s <- if (tm == "Variety") get_row(split_wp, tm) else get_row(split_sp, tm)
  s_err_df <- if (tm == "Variety") get_row(split_wp, "Residuals")$Df else get_row(split_sp, "Residuals")$Df
  m <- an_mixed[tm, ]
  data.frame(
    Effect          = sub(":", " x ", tm),
    RCBD_error_df   = wrong_res_df,
    RCBD_F          = round(w$`F value`, 2),
    RCBD_p          = signif(w$`Pr(>F)`, 3),
    SplitPlot_error_df = s_err_df,
    SplitPlot_F     = round(s$`F value`, 2),
    SplitPlot_p     = signif(s$`Pr(>F)`, 3),
    MixedModel_F    = round(m$`F value`, 2),
    MixedModel_p    = signif(m$`Pr(>F)`, 3)
  )
}))
print(comparison, row.names = FALSE)
write.csv(comparison, "output/tables/anova_comparison_rcbd_vs_splitplot.csv", row.names = FALSE)

cat("\nKey point: in the RCBD analysis Variety is tested against the small subplot\n",
    "error with ", wrong_res_df, " df. In the correct analysis it is tested against the\n",
    "much larger whole-plot error with only ", get_row(split_wp, "Residuals")$Df, " df.\n", sep = "")

# -----------------------------------------------------------------------------
# 6. Assumption checks (mixed model residuals)
# -----------------------------------------------------------------------------
section("6. ASSUMPTION CHECKS")
oats$fitted <- fitted(m_mixed)
oats$resid  <- residuals(m_mixed)
sw  <- shapiro.test(oats$resid)
lev <- leveneTest(resid ~ Nitrogen, data = oats)
lev_v <- leveneTest(resid ~ Variety, data = oats)
cat("Shapiro-Wilk test of residuals: W =", round(sw$statistic, 3), ", p =", signif(sw$p.value, 3), "\n")
cat("Levene's test, residuals by nitrogen rate: F =", round(lev$`F value`[1], 2),
    ", p =", signif(lev$`Pr(>F)`[1], 3), "\n")
cat("Levene's test, residuals by variety: F =", round(lev_v$`F value`[1], 2),
    ", p =", signif(lev_v$`Pr(>F)`[1], 3), "\n")
cat("Random-effect estimates (BLUPs) are also inspected in the diagnostic figure.\n")

# -----------------------------------------------------------------------------
# 7. Treatment comparisons
# -----------------------------------------------------------------------------
section("7. TREATMENT COMPARISONS (emmeans, Kenward-Roger df)")
emm_options(lmer.df = "kenward-roger")

# The interaction is not significant, so main-effect means are interpreted.
emm_N <- emmeans(m_mixed, ~ Nitrogen)
cld_N <- cld(emm_N, adjust = "tukey", Letters = letters, reversed = TRUE)
cat("Nitrogen means (averaged over varieties), Tukey-adjusted letters:\n")
print(cld_N)

emm_V <- emmeans(m_mixed, ~ Variety)
cld_V <- cld(emm_V, adjust = "tukey", Letters = letters, reversed = TRUE)
cat("\nVariety means (averaged over nitrogen rates), Tukey-adjusted letters:\n")
print(cld_V)

cat("\nNitrogen is quantitative, so the response shape is tested with polynomial contrasts:\n")
poly_N <- contrast(emm_N, "poly")
print(poly_N)

clean_cld <- function(x, factor_name) {
  x <- as.data.frame(x)
  data.frame(Level = x[[factor_name]], Mean = round(x$emmean, 1), SE = round(x$SE, 2),
             df = round(x$df, 1), Lower_95CI = round(x$lower.CL, 1),
             Upper_95CI = round(x$upper.CL, 1), Group = trimws(x$.group))
}
means_N <- clean_cld(cld_N, "Nitrogen"); names(means_N)[1] <- "Nitrogen_cwt_per_acre"
means_V <- clean_cld(cld_V, "Variety");  names(means_V)[1] <- "Variety"
means_N <- means_N[order(as.numeric(as.character(means_N[[1]]))), ]
write.csv(means_N, "output/tables/means_nitrogen_tukey.csv", row.names = FALSE)
write.csv(means_V, "output/tables/means_variety_tukey.csv", row.names = FALSE)
write.csv(as.data.frame(summary(poly_N)), "output/tables/nitrogen_polynomial_contrasts.csv", row.names = FALSE)

emm_VN <- as.data.frame(emmeans(m_mixed, ~ Variety * Nitrogen))
emm_VN <- data.frame(Variety = emm_VN$Variety, Nitrogen = emm_VN$Nitrogen,
                     emmean = emm_VN$emmean, SE = emm_VN$SE)
cell <- data.frame(Variety = emm_VN$Variety, Nitrogen_cwt_per_acre = emm_VN$Nitrogen,
                   Mean = round(emm_VN$emmean, 1), SE = round(emm_VN$SE, 2))
write.csv(cell, "output/tables/cell_means_variety_x_nitrogen.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 8. Figures
# -----------------------------------------------------------------------------
# Categorical colours: first three slots of a colour-vision-deficiency
# validated palette (all pairs separable). Text stays in neutral ink.
pal  <- c("Golden Rain" = "#2a78d6", "Marvellous" = "#eb6834", "Victory" = "#1baf7a")
ink  <- "#0b0b0b"; ink2 <- "#52514e"; grid_col <- "#e6e5e1"

theme_pub <- function(base = 11) {
  theme_minimal(base_size = base, base_family = "DejaVu Sans") +
    theme(
      text = element_text(colour = ink),
      axis.text = element_text(colour = ink2),
      axis.title = element_text(colour = ink2),
      plot.title = element_text(face = "bold", size = base + 3),
      plot.subtitle = element_text(colour = ink2),
      plot.caption = element_text(colour = ink2, size = base - 2, hjust = 0),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.4),
      legend.position = "top", legend.justification = "left",
      legend.title = element_text(colour = ink2),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title.position = "plot", plot.caption.position = "plot"
    )
}

# --- Figure 1: variety x nitrogen response with standard errors -------------
emm_VN$N_rate <- as.numeric(as.character(emm_VN$Nitrogen))
dodge <- position_dodge(width = 0.04)
ends  <- emm_VN[emm_VN$N_rate == max(emm_VN$N_rate), ]
ends$label_y <- ends$emmean + c("Golden Rain" = -0.5, "Marvellous" = 4.5, "Victory" = -4)[as.character(ends$Variety)]
fig1 <- ggplot(emm_VN, aes(N_rate, emmean, colour = Variety, group = Variety)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.025,
                linewidth = 0.5, position = dodge) +
  geom_line(linewidth = 0.9, position = dodge) +
  geom_point(size = 3, position = dodge, stroke = 1.2, shape = 21, fill = "white") +
  geom_text(data = ends, aes(y = label_y, label = Variety), hjust = 0, nudge_x = 0.035,
            size = 3.6, colour = ink, show.legend = FALSE) +
  annotate("text", x = means_N$Nitrogen_cwt_per_acre |> as.character() |> as.numeric(),
           y = min(emm_VN$emmean - emm_VN$SE) - 6, label = means_N$Group,
           colour = ink, fontface = "bold", size = 4) +
  scale_colour_manual(values = pal) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), expand = expansion(mult = c(0.05, 0.28))) +
  labs(title = "Oat yield responds to nitrogen; varieties do not differ significantly",
       subtitle = "Means +/- SE from the split-plot mixed model. Letters: nitrogen groups (Tukey, p < 0.05)",
       x = "Nitrogen rate (cwt per acre)", y = "Yield (1/4 lb per subplot)", colour = "Variety",
       caption = "Data: Yates (1935) oats split-plot trial, R package MASS. Self-initiated portfolio project using public data.") +
  theme_pub()
ggsave("output/figures/fig1_variety_nitrogen_response.png", fig1, width = 8, height = 5.2, dpi = 300, bg = "white")

# --- Figure 2: p-values, wrong vs correct analysis --------------------------
pv <- rbind(
  data.frame(Effect = comparison$Effect, Analysis = "Incorrect: analysed as RCBD", p = comparison$RCBD_p),
  data.frame(Effect = comparison$Effect, Analysis = "Correct: split-plot mixed model", p = comparison$MixedModel_p)
)
pv$Effect   <- factor(pv$Effect, levels = rev(comparison$Effect))
pv$Analysis <- factor(pv$Analysis, levels = c("Incorrect: analysed as RCBD", "Correct: split-plot mixed model"))
pv$p_plot   <- pmax(pv$p, 1e-12)
pv$lab      <- ifelse(pv$p < 0.001, formatC(pv$p, format = "e", digits = 1), formatC(pv$p, format = "f", digits = 3))
pal2 <- c("Incorrect: analysed as RCBD" = "#eb6834", "Correct: split-plot mixed model" = "#2a78d6")
seg <- reshape(pv[, c("Effect", "Analysis", "p_plot")], idvar = "Effect", timevar = "Analysis", direction = "wide")
names(seg) <- c("Effect", "p_wrong", "p_right")
fig2 <- ggplot(pv, aes(p_plot, Effect)) +
  scale_y_discrete() +
  annotate("rect", xmin = 0.05, xmax = 1.5, ymin = -Inf, ymax = Inf, fill = "#f3f2ef") +
  geom_vline(xintercept = 0.05, colour = ink2, linetype = "dashed", linewidth = 0.5) +
  geom_segment(data = seg, aes(x = p_wrong, xend = p_right, y = Effect, yend = Effect),
               colour = "#b9b8b2", linewidth = 0.8, inherit.aes = FALSE) +
  geom_point(aes(colour = Analysis), size = 4) +
  geom_text(aes(label = lab, vjust = ifelse(Analysis == "Incorrect: analysed as RCBD", -1.1, 2.1)), size = 3.3, colour = ink) +
  scale_colour_manual(values = pal2) +
  scale_x_log10(breaks = c(1e-10, 1e-6, 1e-3, 0.05, 1), labels = c("1e-10", "1e-6", "0.001", "0.05", "1"),
                limits = c(1e-12, 1.5)) +
  labs(title = "The wrong (RCBD) analysis makes the variety effect falsely significant",
       subtitle = "p-values for the same data under both analyses. Shaded area: p > 0.05 (not significant)",
       x = "p-value (log scale)", y = NULL, colour = NULL,
       caption = "The RCBD analysis tests Variety against one pooled error (55 df) instead of the whole-plot error (Block x Variety, 10 df).") +
  theme_pub() + theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.4))
ggsave("output/figures/fig2_rcbd_vs_splitplot_pvalues.png", fig2, width = 8.5, height = 4.4, dpi = 300, bg = "white")

# --- Figure 3: diagnostics ---------------------------------------------------
pt_col <- "#52514e"
d1 <- ggplot(oats, aes(fitted, resid)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.4) +
  geom_point(colour = pt_col, size = 2.2, alpha = 0.8) +
  labs(title = "Residuals vs fitted", x = "Fitted yield", y = "Residual") + theme_pub(10)
d2 <- ggplot(oats, aes(sample = resid)) +
  stat_qq_line(colour = ink2, linewidth = 0.5) +
  stat_qq(colour = pt_col, size = 2.2, alpha = 0.8) +
  labs(title = "Normal Q-Q of residuals", subtitle = sprintf("Shapiro-Wilk p = %.2f", sw$p.value),
       x = "Theoretical quantiles", y = "Sample quantiles") + theme_pub(10)
d3 <- ggplot(oats, aes(Nitrogen, resid)) +
  geom_hline(yintercept = 0, colour = ink2, linewidth = 0.4) +
  geom_boxplot(fill = "#f3f2ef", colour = ink2, width = 0.55, outlier.shape = NA) +
  geom_jitter(width = 0.12, height = 0, colour = pt_col, size = 1.8, alpha = 0.8) +
  labs(title = "Residual spread by nitrogen rate", subtitle = sprintf("Levene p = %.2f", lev$`Pr(>F)`[1]),
       x = "Nitrogen rate (cwt per acre)", y = "Residual") + theme_pub(10)
re_wp <- ranef(m_mixed)[["Block:Variety"]][, 1]
d4 <- ggplot(data.frame(b = re_wp), aes(sample = b)) +
  stat_qq_line(colour = ink2, linewidth = 0.5) +
  stat_qq(colour = pt_col, size = 2.2, alpha = 0.8) +
  labs(title = "Normal Q-Q of whole-plot effects", x = "Theoretical quantiles", y = "Estimated effect") + theme_pub(10)
fig3 <- (d1 + d2) / (d3 + d4) +
  plot_annotation(title = "Model diagnostics: no strong evidence against normality or equal variance",
                  theme = theme(plot.title = element_text(face = "bold", size = 13, colour = ink,
                                                          family = "DejaVu Sans")))
ggsave("output/figures/fig3_model_diagnostics.png", fig3, width = 9, height = 7, dpi = 300, bg = "white")

# --- Figure 4: nitrogen main effect with Tukey letters ----------------------
mN <- means_N; mN$N_rate <- as.numeric(as.character(mN$Nitrogen_cwt_per_acre))
fig4 <- ggplot(mN, aes(factor(N_rate), Mean)) +
  geom_col(fill = "#2a78d6", width = 0.6) +
  geom_errorbar(aes(ymin = Lower_95CI, ymax = Upper_95CI), width = 0.15, colour = ink, linewidth = 0.5) +
  geom_text(aes(y = Upper_95CI, label = Group), vjust = -0.7, fontface = "bold", size = 4.2, colour = ink) +
  geom_text(aes(y = 4, label = sprintf("%.1f", Mean)), vjust = 0, size = 3.5, colour = "white") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Yield rises with nitrogen; 0.4 and 0.6 cwt do not differ",
       subtitle = "Means over varieties with 95% CI. Different letters differ (Tukey, p < 0.05)",
       x = "Nitrogen rate (cwt per acre)", y = "Yield (1/4 lb per subplot)",
       caption = "CIs include block-to-block variation; Tukey tests use the within-block error, so overlapping CIs can still differ.\nData: Yates (1935) oats split-plot trial, R package MASS.") +
  theme_pub()
ggsave("output/figures/fig4_nitrogen_means_tukey.png", fig4, width = 7.5, height = 5, dpi = 300, bg = "white")

# --- Cover image for the Upwork portfolio (16:9) ----------------------------
cover <- fig1 +
  labs(title = "Getting the Split-Plot Right: Mixed Model Analysis in R",
       subtitle = "Oat variety x nitrogen field trial | correct error terms, Tukey letters, reproducible script",
       caption = "Self-initiated portfolio project using public data (Yates 1935, R package MASS)") +
  theme(plot.title = element_text(size = 17, face = "bold"),
        plot.subtitle = element_text(size = 12),
        plot.margin = margin(18, 24, 12, 24))
ggsave("output/figures/cover_upwork.png", cover, width = 10, height = 5.625, dpi = 200, bg = "white")

# -----------------------------------------------------------------------------
section("SESSION INFO")
print(sessionInfo())
sink()
message("Done. Outputs in output/")
