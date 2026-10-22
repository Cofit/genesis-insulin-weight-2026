# ==============================================================================
# ANALYSIS CODE — Submission Version
# ==============================================================================
# Title:   Insulin response and protein intake are associated with distinct
#          aspects of weight loss during obesity treatment
# Authors: Yi Jing Ngu, Tang-Yueh Li, Ya-Wen Chen, et al.
# Journal: Submitted to Nature Health (under review)
# Date:    2026-10-22
#
# Description:
#   This script reproduces all tables and figures in the manuscript.
#   It requires a single cleaned input CSV (sample_data.csv)
#   and writes all outputs (XLSX, PNG, PDF) to OUT_DIR.
#   Run the script sequentially from top to bottom; later sections depend on
#   objects created in earlier sections (see dependency notes in section headers).
#
# Input:
#   DATA_DIR/sample_data.csv
#     — Analysis-ready dataset: all clinical variables (W0 + W8),
#       3D circumference columns (W0/W8), derived delta variables,
#       and disease groupings. All cleaning steps have been applied.
#     — Column dictionary is provided as a supplementary data file.
#
# Software:
#   R version 4.5.1
#   Key packages: dplyr, ggplot2, patchwork, writexl, car, uwot,
#                 ggnewscale, MASS (see Section 0 for full list)
# ==============================================================================

# ==============================================================================
# 0. Packages
# ==============================================================================
library(pacman)
pacman::p_load(
  readxl, dplyr, writexl, lubridate, readr, tidyr, forcats, purrr,
  ggplot2, patchwork, cowplot, ggpubr, scales, showtext,
  table1, openxlsx, gt,
  rstatix, FSA, car,
  uwot, ggnewscale, MASS
)

# ==============================================================================
# I. Data Preparation
# ==============================================================================
# Edit DATA_DIR and OUT_DIR to match your local environment.
#   DATA_DIR : directory containing sample_data.csv
#   OUT_DIR  : directory where all output files will be written

DATA_DIR <- ""   # <-- edit this path
OUT_DIR  <- ""   # <-- edit this path

setwd(OUT_DIR)

df_complete_full <- read.csv(file.path(DATA_DIR, "sample_data.csv"))

# Convenience subsets (used by multiple downstream sections)
df_nondm <- df_complete_full %>% filter(diabetes == "No")
df_dm    <- df_complete_full %>% filter(diabetes == "Yes")

# Reconstruct body_3d_wide proxy (3D scanner data pre-merged into main dataset)
# Required by Figure 5 (UMAP) section for column-name identification.
body_3d_wide <- df_complete_full

# Pre-compute 3D circumference base variable names
# (Columns matching "circumference" and ending in "_w0", stripped of "_w0" suffix)
circ_base_vars <- names(df_complete_full)[
  grepl("circumference", names(df_complete_full), ignore.case = TRUE) &
    grepl("_w0$",        names(df_complete_full)) &
    grepl("\\.",         names(df_complete_full))   # 3D scanner cols use CamelCase.Dot naming
]
circ_base_vars <- sub("_w0$", "", circ_base_vars)
cat("3D circumference variables identified:", length(circ_base_vars), "\n")

# ==============================================================================
# II. Helper Functions
# ==============================================================================

# ── GLM extractor ────────────────────────────────────────────────────────────
# Adj R² computed directly from deviance
# For Gaussian GLM: R² = 1 - Deviance/Null_Deviance  (equivalent to OLS R²)
# Adj R² = 1 - (1 - R²) × (n-1) / (n-k-1)
extract_glm <- function(model, model_name) {
  coef_table <- summary(model)$coefficients
  ci_table   <- confint(model)
  n      <- length(model$fitted.values)
  k      <- length(model$coefficients) - 1   # predictors, exclude intercept
  r2     <- 1 - model$deviance / model$null.deviance
  adj_r2 <- round(1 - (1 - r2) * (n - 1) / (n - k - 1), 4)
  data.frame(
    Model    = model_name,
    Variable = rownames(coef_table),
    Beta     = round(coef_table[, "Estimate"],   3),
    SE       = round(coef_table[, "Std. Error"], 3),
    CI_95    = paste0(round(ci_table[, 1], 3), " to ", round(ci_table[, 2], 3)),
    P_value  = ifelse(coef_table[, "Pr(>|t|)"] < 0.001, "<.001",
                      as.character(round(coef_table[, "Pr(>|t|)"], 3))),
    Sig      = ifelse(coef_table[, "Pr(>|t|)"] < 0.001, "***",
                      ifelse(coef_table[, "Pr(>|t|)"] < 0.01,  "**",
                             ifelse(coef_table[, "Pr(>|t|)"] < 0.05,  "*",
                                    ifelse(coef_table[, "Pr(>|t|)"] < 0.1, ".", "ns")))),
    AIC    = round(AIC(model), 1),
    BIC    = round(BIC(model), 1),
    Adj_R2 = adj_r2,
    N      = n
  )
}

# ── Post-hoc standardized betas: β_std = β_raw × (SD_X / SD_Y) ───────────
# SD_X: from model.matrix columns (dummy variables use SD of 0/1 column).
# SD_Y: from model.frame first column (log-transformed Y uses log-scale SD).
# Intercept row carries NA for standardized columns.
standardize_beta_posthoc <- function(model) {
  mf    <- model.frame(model)
  mm    <- model.matrix(model)                  # dummy-expanded predictors
  sd_y  <- sd(mf[[1]], na.rm = TRUE)            # SD of outcome (log scale if logged)
  sd_x  <- apply(mm[, -1, drop = FALSE], 2, sd, na.rm = TRUE)
  sf    <- sd_x / sd_y                          # scale factors

  coef_raw <- coef(model)[-1]
  se_raw   <- sqrt(diag(vcov(model)))[-1]
  ci_raw   <- suppressMessages(confint(model))[-1, , drop = FALSE]

  data.frame(
    Variable  = names(coef_raw),
    Beta_std  = round(coef_raw * sf, 4),
    SE_std    = round(se_raw   * sf, 4),
    CI_std_95 = paste0(round(ci_raw[, 1] * sf, 4), " to ",
                       round(ci_raw[, 2] * sf, 4)),
    SD_X      = round(sd_x, 4),
    SD_Y      = round(sd_y, 4),
    stringsAsFactors = FALSE
  )
}

# Wrapper: extract_glm output + standardized beta columns merged
extract_glm_with_std <- function(model, model_name) {
  base <- extract_glm(model, model_name)
  std  <- standardize_beta_posthoc(model)
  sd_y <- round(sd(model.frame(model)[[1]], na.rm = TRUE), 4)
  # Intercept row: NAs for std columns (intercept is not standardized)
  intercept_std <- data.frame(
    Variable  = "(Intercept)",
    Beta_std  = NA_real_, SE_std = NA_real_,
    CI_std_95 = NA_character_,
    SD_X      = NA_real_,  SD_Y  = sd_y,
    stringsAsFactors = FALSE
  )
  std_full <- rbind(intercept_std, std)
  dplyr::left_join(base, std_full, by = "Variable")
}

# ── Forest plot data prep ─────────────────────────────────────────────────────
prepare_forest <- function(model, var_labels) {
  coef_table <- summary(model)$coefficients
  ci_table   <- confint(model)
  data.frame(
    Variable = rownames(coef_table),
    Beta     = coef_table[, "Estimate"],
    Lower    = ci_table[, 1],
    Upper    = ci_table[, 2],
    P_value  = coef_table[, "Pr(>|t|)"]
  ) %>%
    filter(Variable != "(Intercept)") %>%
    mutate(
      Variable = dplyr::recode(Variable, !!!var_labels),
      Variable = factor(Variable, levels = rev(var_labels))
    )
}

# ── Forest plot draw ──────────────────────────────────────────────────────────
# Significance stars (from P_value column) are shown to the right of each CI bar:
#   *** p < .001  |  ** p < .01  |  * p < .05  |  . p < .1  |  (blank) ns
plot_forest <- function(df, title) {
  ggplot(df, aes(x = Beta, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "solid",
               color = "gray40", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                   height = 0.25, linewidth = 0.8, color = "#C0392B") +
    geom_point(size = 3, color = "#C0392B") +
    geom_text(
      aes(x = Upper,
          label = dplyr::case_when(
            P_value < 0.001 ~ "***",
            P_value < 0.01  ~ "**",
            P_value < 0.05  ~ "*",
            P_value < 0.1   ~ ".",
            TRUE            ~ ""
          )),
      hjust = -0.25, vjust = 0.4,
      size = 4, color = "gray20", fontface = "bold"
    ) +
    labs(title = title, x = "Beta Coefficient (95% CI)", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      plot.title          = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.y         = element_text(size = 11, hjust = 1, color = "black"),
      axis.ticks.y        = element_line(linewidth = 0.4),
      axis.ticks.length.y = unit(0.3, "cm"),
      panel.grid          = element_blank(),
      legend.position     = "none"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.12)))  # right margin for stars
}

# ── Half violin helper (no extra packages — uses base density + geom_polygon) ─
# side = "r": density extends rightward (Baseline)
# side = "l": density extends leftward  (Week 8)
.half_dens_poly <- function(values, x_center, side = "r", width = 0.38) {
  vals <- na.omit(values)
  if (length(vals) < 5) return(NULL)
  d <- density(vals, bw = "nrd0")
  s <- d$y / max(d$y) * width                     # scale to [0, width]
  n <- length(d$x)
  if (side == "r")
    data.frame(px = c(x_center, x_center + s, x_center),
               py = c(d$x[1],   d$x,          d$x[n]))
  else
    data.frame(px = c(x_center, x_center - s, x_center),
               py = c(d$x[1],   d$x,          d$x[n]))
}

# Gender colours
# Non-T2DM: original shades; T2DM: one shade lighter
.gender_colors_nondm <- c("Male" = "#4E79A7", "Female" = "#E15759")
.gender_colors_dm    <- c("Male" = "#9BBFE0", "Female" = "#F0A6A7")

# Combined lookup for faceted plots (diabetes_label + Gender → colour)
# Keys must match paste(diabetes_label, Gender) exactly;
# diabetes_label uses "non-T2DM" (lowercase n) from the data pipeline (line 64)
.gender_colors_by_dm <- c(
  "non-T2DM Male"   = "#4E79A7",
  "non-T2DM Female" = "#E15759",
  "T2DM Male"       = "#9BBFE0",
  "T2DM Female"     = "#F0A6A7"
)

# ── Violin panel (single DM group, no facet, W0 vs W8) ───────────────────────
# Half violin split by Gender (Male = left, Female = right) at each timepoint
# dm_group: "T2DM" or "non-T2DM"
# ylim: optional c(lower, upper) — coord_cartesian clips display only
make_violin_single <- function(df, var_w0, var_w8, y_label, dm_group,
                               ylim = NULL) {
  # Select colour palette based on DM group
  gender_cols <- if (dm_group == "T2DM") .gender_colors_dm else .gender_colors_nondm

  long_data <- df %>%
    filter(!is.na(.data[[var_w0]]), !is.na(.data[[var_w8]]),
           diabetes_label == dm_group, !is.na(gender)) %>%
    dplyr::select(client_id, gender,
           Baseline = all_of(var_w0), `Week 8` = all_of(var_w8)) %>%
    pivot_longer(cols = c(Baseline, `Week 8`),
                 names_to = "Timepoint", values_to = "Value") %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("Baseline", "Week 8")),
      Gender    = factor(
        ifelse(tolower(as.character(gender)) == "male", "Male", "Female"),
        levels = c("Male", "Female")
      )
    )

  # Half violin polygons: Male = left, Female = right, per Timepoint
  tp_map <- c("Baseline" = 1, "Week 8" = 2)
  poly_df <- bind_rows(lapply(names(tp_map), function(tp) {
    xp <- tp_map[[tp]]
    bind_rows(
      .half_dens_poly(
        long_data$Value[long_data$Timepoint == tp & long_data$Gender == "Male"],
        x_center = xp, side = "l") %>%
        mutate(Gender = "Male",   Timepoint = tp),
      .half_dens_poly(
        long_data$Value[long_data$Timepoint == tp & long_data$Gender == "Female"],
        x_center = xp, side = "r") %>%
        mutate(Gender = "Female", Timepoint = tp)
    )
  })) %>%
    mutate(
      Timepoint = factor(Timepoint, levels = c("Baseline", "Week 8")),
      Gender    = factor(Gender, levels = c("Male", "Female"))
    )

  p <- ggplot() +
    geom_polygon(data = poly_df,
                 aes(x = px, y = py, fill = Gender,
                     group = interaction(Gender, Timepoint)),
                 alpha = 0.7) +
    # Male boxplot: nudged left
    geom_boxplot(data = long_data %>% filter(Gender == "Male"),
                 aes(x = as.numeric(Timepoint), y = Value, group = Timepoint),
                 width = 0.05, fill = gender_cols["Male"], color = "black",
                 outlier.shape = NA, alpha = 0.9,
                 position = position_nudge(x = -0.12)) +
    # Female boxplot: nudged right
    geom_boxplot(data = long_data %>% filter(Gender == "Female"),
                 aes(x = as.numeric(Timepoint), y = Value, group = Timepoint),
                 width = 0.05, fill = gender_cols["Female"], color = "black",
                 outlier.shape = NA, alpha = 0.9,
                 position = position_nudge(x = +0.12)) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Baseline", "Week 8"),
                        limits = c(0.4, 2.6)) +
    scale_fill_manual(values = gender_cols, name = "Gender") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    labs(x = NULL, y = y_label) +
    theme_minimal(base_size = 16) +
    theme(
      legend.position    = "bottom",
      legend.text        = element_text(size = 15),
      axis.text.x        = element_text(color = "black", size = 13),
      axis.text.y        = element_text(color = "black", size = 14),
      axis.title.y       = element_text(size = 15),
      panel.grid.major.x = element_blank()
    )

  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  return(p)
}

# ── Median (IQR) string ───────────────────────────────────────────────────────
fmt_med_iqr <- function(x) {
  x <- na.omit(x)
  if (length(x) == 0) return("—")
  sprintf("%.2f (%.2f–%.2f)", median(x), quantile(x, 0.25), quantile(x, 0.75))
}

# ── W0 vs W8 paired Wilcoxon table builder ────────────────────────────────────
# Delta = W8 - W0; uses pre-computed delta_* columns from df when available
# (e.g. weight_w0 -> delta_weight), falls back to inline W8-W0 otherwise.
make_paired_table <- function(df, var_pairs, var_labels) {
  purrr::map2_dfr(var_pairs, var_labels, function(vp, lbl) {
    var_w0 <- vp[1]; var_w8 <- vp[2]

    d <- df %>%
      filter(!is.na(.data[[var_w0]]), !is.na(.data[[var_w8]])) %>%
      dplyr::select(w0 = all_of(var_w0), w8 = all_of(var_w8))
    if (nrow(d) < 5) return(NULL)

    # Look up pre-computed delta_* column; fall back to inline computation
    base_name <- sub("_w0$", "", var_w0)
    delta_col <- paste0("delta_", base_name)
    if (delta_col %in% colnames(df)) {
      delta_vals <- df %>%
        filter(!is.na(.data[[var_w0]]), !is.na(.data[[var_w8]])) %>%
        pull(all_of(delta_col))
    } else {
      delta_vals <- d$w8 - d$w0
    }

    # % change per person: (W8 - W0) / |W0| * 100
    # W0 = 0 is set to NA to avoid division by zero
    pct_change_vals <- ifelse(
      d$w0 != 0 & !is.na(d$w0) & !is.na(d$w8),
      (d$w8 - d$w0) / abs(d$w0) * 100,
      NA_real_
    )

    wt <- wilcox.test(d$w0, d$w8, paired = TRUE, exact = FALSE)
    n  <- nrow(d)

    # Compute Z directly from Wilcoxon V statistic (sum of positive signed ranks)
    # Avoids Inf from back-calculating Z via qnorm(p.value/2) when p ≈ 0
    V     <- as.numeric(wt$statistic)
    diffs <- d$w8 - d$w0
    n_eff <- sum(diffs != 0, na.rm = TRUE)   # non-zero differences used by signed-rank test
    z_stat <- (V - n_eff * (n_eff + 1) / 4) /
              sqrt(n_eff * (n_eff + 1) * (2 * n_eff + 1) / 24)
    r <- round(min(abs(z_stat) / sqrt(n), 1), 3)  # r in [0, 1]
    data.frame(
      Variable                    = lbl,
      N                           = n,
      `Baseline Median (IQR)`     = fmt_med_iqr(d$w0),
      `Week 8 Median (IQR)`       = fmt_med_iqr(d$w8),
      `Delta Median (IQR)`        = fmt_med_iqr(delta_vals),
      `% Change Median (IQR)`     = fmt_med_iqr(pct_change_vals),
      `p-value`                   = ifelse(wt$p.value < 0.001, "<.001",
                                           sprintf("%.3f", wt$p.value)),
      `Effect size (r)`           = r,
      check.names = FALSE
    )
  })
}

# ── Standardized-beta forest plots ───────────────────────────────────────────
# prepare_forest_std: builds forest df from Beta_std / CI_std_95 columns
# plot_forest_std:    wraps plot_forest with x-axis re-labelled "Std. Beta"
# Used by: Supplementary Table 4, Figure 3
prepare_forest_std <- function(model, var_labels) {
  std    <- standardize_beta_posthoc(model)
  p_vals <- summary(model)$coefficients[-1, "Pr(>|t|)"]   # drop intercept

  # parse CI string "x to y" produced by standardize_beta_posthoc
  ci_parts <- strsplit(std$CI_std_95, " to ")
  lower    <- as.numeric(sapply(ci_parts, `[`, 1))
  upper    <- as.numeric(sapply(ci_parts, `[`, 2))

  data.frame(
    Variable = std$Variable,
    Beta     = std$Beta_std,
    Lower    = lower,
    Upper    = upper,
    P_value  = p_vals[std$Variable]
  ) %>%
    mutate(
      Variable = dplyr::recode(Variable, !!!var_labels),
      Variable = factor(Variable, levels = rev(var_labels))
    )
}

# plot_forest_std: same aesthetics as plot_forest but x-axis label differs
plot_forest_std <- function(df, title) {
  plot_forest(df, title) +
    labs(x = "Standardized Beta (95% CI)")
}


# ==============================================================================
# Table 1 — Baseline Demographics (Total / Non-T2DM / T2DM)
# Output: table1_baseline_clinical_characteristics.xlsx
# ==============================================================================

tbl1_cont_vars <- list(
  list(var = "age",                               label = "Age (years)"),
  list(var = "weight_w0",                         label = "Weight (kg)"),
  list(var = "bmi_w0",                            label = "BMI (kg/m²)"),
  list(var = "waist_circumference_w0",             label = "Waist Circumference (cm)"),
  list(var = "body_fat_mass_w0",                  label = "Fat Mass (kg)"),
  list(var = "body_fat_mass_percentage_w0",        label = "Fat Mass (%)"),
  list(var = "weight_without_fat_w0",              label = "Fat-Free Mass (kg)"),
  list(var = "muscle_mass_w0",                    label = "Muscle Mass (kg)"),
  list(var = "hba1c_w0",                          label = "HbA1c (%)"),
  list(var = "glucose_ac_w0",                     label = "Fasting Glucose (mg/dL)"),
  list(var = "insulin_w0",                        label = "Fasting Insulin (μU/mL)"),
  list(var = "homa_ir_w0",                        label = "HOMA-IR"),
  list(var = "triglyceride_w0",                   label = "Triglyceride (mg/dL)"),
  list(var = "total_cholesterol_w0",              label = "Total Cholesterol (mg/dL)"),
  list(var = "high_density_lipoprotein_w0",       label = "HDL-C (mg/dL)"),
  list(var = "low_density_lipoprotein_cholesterol_w0", label = "LDL-C (mg/dL)"),
  list(var = "systolic_blood_pressure_w0",        label = "Systolic BP (mmHg)"),
  list(var = "diastolic_blood_pressure_w0",       label = "Diastolic BP (mmHg)"),
  list(var = "glucose_pc_2hr_w0",                 label = "2-hr OGTT Glucose (mg/dL)"),
  list(var = "insulin_pc_1hr_w0",                 label = "1-hr OGTT Insulin (μU/mL)"),
  list(var = "insulin_pc_2hr_w0",                 label = "2-hr OGTT Insulin (μU/mL)"),
  list(var = "total_energy",                       label = "Total energy intake (kcal/day)"),
  list(var = "diet_w8_protein",                   label = "Protein Intake (g/day)"),
  list(var = "diet_w8_fat",                       label = "Fat Intake (g/day)"),
  list(var = "diet_w8_carbohydrate",              label = "Carbohydrate Intake (g/day)")
)

df_nondm <- df_complete_full %>% filter(diabetes == "No")
df_dm    <- df_complete_full %>% filter(diabetes == "Yes")

tbl1_rows_cont <- purrr::map_dfr(tbl1_cont_vars, function(item) {
  v <- item$var
  if (!v %in% colnames(df_complete_full)) return(NULL)

  x_all   <- na.omit(df_complete_full[[v]])
  x_nondm <- na.omit(df_nondm[[v]])
  x_dm    <- na.omit(df_dm[[v]])

  # Wilcoxon unpaired test (Non-T2DM vs T2DM)
  if (length(x_nondm) >= 3 && length(x_dm) >= 3) {
    wt      <- wilcox.test(x_nondm, x_dm, exact = FALSE)
    n1      <- length(x_nondm)
    n2      <- length(x_dm)
    n_cmp   <- n1 + n2
    # Direct Z from Wilcoxon W statistic (avoids Inf when p ≈ 0)
    # Unpaired Mann–Whitney: W ~ N(n1*n2/2, sqrt(n1*n2*(n1+n2+1)/12))
    W       <- as.numeric(wt$statistic)
    mu_W    <- n1 * n2 / 2
    sigma_W <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    z_val   <- (W - mu_W) / sigma_W
    r_val   <- round(min(abs(z_val) / sqrt(n_cmp), 1), 3)
    p_str   <- ifelse(wt$p.value < 0.001, "<.001", sprintf("%.3f", wt$p.value))
  } else {
    p_str <- "—"; r_val <- NA
  }

  data.frame(
    Variable             = item$label,
    Type                 = "Continuous",
    `Total N`            = length(x_all),
    `Total Median (IQR)` = fmt_med_iqr(x_all),
    `Non-T2DM N`         = length(x_nondm),
    `Non-T2DM Median (IQR)` = fmt_med_iqr(x_nondm),
    `T2DM N`             = length(x_dm),
    `T2DM Median (IQR)`  = fmt_med_iqr(x_dm),
    `p-value`            = p_str,
    `Effect size (r)`    = ifelse(is.na(r_val), "—", as.character(r_val)),
    check.names = FALSE
  )
})

# Gender row (categorical)
# Build explicit 2×2 matrix to avoid named-vector issues from table()
n_male_total   <- sum(df_complete_full$gender == "male",   na.rm = TRUE)
n_female_total <- sum(df_complete_full$gender == "female", na.rm = TRUE)
n_male_nondm   <- sum(df_nondm$gender == "male",           na.rm = TRUE)
n_female_nondm <- sum(df_nondm$gender == "female",         na.rm = TRUE)
n_male_dm      <- sum(df_dm$gender == "male",              na.rm = TRUE)
n_female_dm    <- sum(df_dm$gender == "female",            na.rm = TRUE)
n_total_gender <- n_male_total + n_female_total

gmat <- matrix(
  c(n_male_nondm, n_female_nondm,
    n_male_dm,    n_female_dm),
  nrow = 2, byrow = FALSE,
  dimnames = list(c("male", "female"), c("No", "Yes"))
)

chi_test <- chisq.test(gmat, correct = FALSE)
phi_coef <- round(sqrt(unname(chi_test$statistic) / sum(gmat)), 3)
p_gender <- ifelse(chi_test$p.value < 0.001, "<.001",
                   sprintf("%.3f", chi_test$p.value))

gender_rows <- data.frame(
  Variable             = c("Gender — Male", "Gender — Female"),
  Type                 = "Categorical",
  `Total N`            = c(n_male_total,   n_female_total),
  `Total Median (IQR)` = c(
    sprintf("%d (%.1f%%)", n_male_total,   n_male_total   / n_total_gender * 100),
    sprintf("%d (%.1f%%)", n_female_total, n_female_total / n_total_gender * 100)
  ),
  `Non-T2DM N`         = c(n_male_nondm,   n_female_nondm),
  `Non-T2DM Median (IQR)` = c(
    sprintf("%d (%.1f%%)", n_male_nondm,   n_male_nondm   / nrow(df_nondm) * 100),
    sprintf("%d (%.1f%%)", n_female_nondm, n_female_nondm / nrow(df_nondm) * 100)
  ),
  `T2DM N`             = c(n_male_dm,   n_female_dm),
  `T2DM Median (IQR)`  = c(
    sprintf("%d (%.1f%%)", n_male_dm,   n_male_dm   / nrow(df_dm) * 100),
    sprintf("%d (%.1f%%)", n_female_dm, n_female_dm / nrow(df_dm) * 100)
  ),
  `p-value`         = c(p_gender, ""),
  `Effect size (r)` = c(as.character(phi_coef), ""),
  check.names = FALSE
)

table_1 <- bind_rows(gender_rows, tbl1_rows_cont)

write_xlsx(
  list(
    "Table_1" = table_1,
    "Notes"   = data.frame(Note = c(
      "Continuous variables: Median (IQR); Wilcoxon rank-sum test (Non-T2DM vs T2DM).",
      "Effect size r = |Z| / sqrt(N); |r| < 0.1 negligible, 0.1 small, 0.3 medium, 0.5 large.",
      "Gender: n (%), Chi-square test; effect size = phi coefficient.",
      paste0("Total N = ", nrow(df_complete_full),
             " | Non-T2DM N = ", nrow(df_nondm),
             " | T2DM N = ", nrow(df_dm))
    ))
  ),
  path = "table1_baseline_clinical_characteristics.xlsx"
)


# ==============================================================================
# Supplementary Table 1 — Changes in clinical parameters of T2D participants 
#                         in week 0 (baseline) and week 8.
# Groups: Total / T2DM / Non-T2DM
# Variables: body composition + metabolic markers
# Statistics: Median (IQR), paired Wilcoxon p, effect size r
# Output: Ex_data_table1_paired_comparison.xlsx
# ==============================================================================
s1_var_pairs <- list(
  list(c("weight_w0",                          "weight_w8"),                         "Weight (kg)"),
  list(c("bmi_w0",                             "bmi_w8"),                            "BMI (kg/m²)"),
  list(c("waist_circumference_w0",             "waist_circumference_w8"),            "Waist Circumference (cm)"),
  list(c("body_fat_mass_w0",                   "body_fat_mass_w8"),                  "Fat Mass (kg)"),
  list(c("body_fat_mass_percentage_w0",         "body_fat_mass_percentage_w8"),       "Fat Mass (%)"),
  list(c("muscle_mass_w0",                     "muscle_mass_w8"),                    "Muscle Mass (kg)"),
  list(c("muscle_mass_percentage_w0",          "muscle_mass_percentage_w8"),         "Muscle Mass (%)"),
  list(c("hba1c_w0",                           "hba1c_w8"),                          "HbA1c (%)"),
  list(c("glucose_ac_w0",                      "glucose_ac_w8"),                     "Fasting Glucose (mg/dL)"),
  list(c("insulin_w0",                         "insulin_w8"),                        "Fasting Insulin (μU/mL)"),
  list(c("homa_ir_w0",                         "homa_ir_w8"),                        "HOMA-IR"),
  list(c("triglyceride_w0",                    "triglyceride_w8"),                   "Triglyceride (mg/dL)"),
  list(c("total_cholesterol_w0",               "total_cholesterol_w8"),              "Total Cholesterol (mg/dL)"),
  list(c("high_density_lipoprotein_w0",        "high_density_lipoprotein_w8"),       "HDL-C (mg/dL)"),
  list(c("low_density_lipoprotein_cholesterol_w0",
         "low_density_lipoprotein_cholesterol_w8"),                                  "LDL-C (mg/dL)"),
  list(c("systolic_blood_pressure_w0",         "systolic_blood_pressure_w8"),        "Systolic BP (mmHg)"),
  list(c("diastolic_blood_pressure_w0",        "diastolic_blood_pressure_w8"),       "Diastolic BP (mmHg)"),
  list(c("glucose_pc_2hr_w0",                  "glucose_pc_2hr_w8"),                 "2-hr OGTT Glucose (mg/dL)")
)

s1_pairs  <- map(s1_var_pairs, ~ .x[[1]])
s1_labels <- map_chr(s1_var_pairs, ~ .x[[2]])

tbl_s1a <- make_paired_table(df_complete_full, s1_pairs, s1_labels)
tbl_s1b <- make_paired_table(df_dm,            s1_pairs, s1_labels)
tbl_s1c <- make_paired_table(df_nondm,         s1_pairs, s1_labels)

write_xlsx(
  list(
    "Ex_Table1_Total"   = tbl_s1a,
    "Ex_Table1_T2DM"    = tbl_s1b,
    "Ex_Table1_NonT2DM" = tbl_s1c
    ),
  path = "Ex_data_table1_paired_comparison.xlsx"
)


# ==============================================================================
# Supplementary Table 2 — Gender-stratified changes in weight, BMI, and 
#                         body composition at weeks 0 and 8.
# Subgroups: T2DM Male | T2DM Female | Non-T2DM Male | Non-T2DM Female
# Statistics: Median (IQR), paired Wilcoxon p, effect size r
# Output: Ex_data_table2_paired_comparison_by_T2D_Gender.xlsx
# ==============================================================================
df_dm_male    <- df_dm    %>% filter(tolower(as.character(gender)) == "male")
df_dm_female  <- df_dm    %>% filter(tolower(as.character(gender)) == "female")
df_ndm_male   <- df_nondm %>% filter(tolower(as.character(gender)) == "male")
df_ndm_female <- df_nondm %>% filter(tolower(as.character(gender)) == "female")

tbl_s2a <- make_paired_table(df_dm_male,    s1_pairs, s1_labels)
tbl_s2b <- make_paired_table(df_dm_female,  s1_pairs, s1_labels)
tbl_s2c <- make_paired_table(df_ndm_male,   s1_pairs, s1_labels)
tbl_s2d <- make_paired_table(df_ndm_female, s1_pairs, s1_labels)

write_xlsx(
  list(
    "Sup_Table1_T2DM_Male"     = tbl_s2a,
    "Sup_Table1_T2DM_Female"   = tbl_s2b,
    "Sup_Table1_NonT2DM_Male"  = tbl_s2c,
    "Sup_Table1_NonT2DM_Female"= tbl_s2d
  ),
  path = "Sup_data_table2_paired_comparison_by_T2D_Gender.xlsx"
)


# ==============================================================================
# Supplementary Table 3 — Changes in metabolic syndrome, and hyperlipidemia
#                         status from baseline (week 0) to week 8 by
#                         baseline diabetes status.
# Conditions: Metabolic syndrome | Hyperlipidemia
# Groups: Non-T2D / T2D
# Output: Sup_data_table3_met_lipid_status.xlsx
# ==============================================================================
# Complete-case subsets for each disease condition 
df_met_complete <- df_complete_full %>% filter(!is.na(met_gp), !is.na(met_gp_w8))
df_hlp_complete <- df_complete_full %>% filter(!is.na(hlp_gp), !is.na(hlp_gp_w8))

make_disease_table <- function(df_met_c, df_hlp_c, group_label) {
  bind_rows(
    # MetS
    bind_rows(
      df_met_c %>% count(group = met_gp)    %>% mutate(timepoint = "baseline"),
      df_met_c %>% count(group = met_gp_w8) %>% mutate(timepoint = "W8")
    ) %>%
      group_by(timepoint) %>% mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
      mutate(disease = "Metabolic Syndrome", Group = group_label),
    # HLP
    bind_rows(
      df_hlp_c %>% count(group = hlp_gp)    %>% mutate(timepoint = "baseline"),
      df_hlp_c %>% count(group = hlp_gp_w8) %>% mutate(timepoint = "W8")
    ) %>%
      group_by(timepoint) %>% mutate(pct = n / sum(n) * 100) %>% ungroup() %>%
      mutate(disease = "Hyperlipidemia", Group = group_label)
  ) %>%
    mutate(`n (%)` = sprintf("%d (%.1f%%)", n, pct)) %>%
    dplyr::select(Group, disease, group, timepoint, `n (%)`) %>%
    pivot_wider(names_from = timepoint, values_from = `n (%)`,
                names_prefix = "") %>%
    rename(`Classification` = group, Disease = disease,
           `baseline n (%)` = baseline, `W8 n (%)` = W8)
}

# Split by diabetes_label (T2DM vs non-T2DM) for stratified breakdown
met_nondm  <- df_met_complete %>% filter(diabetes_label == "non-T2DM")
met_t2dm   <- df_met_complete %>% filter(diabetes_label == "T2DM")
hlp_nondm  <- df_hlp_complete %>% filter(diabetes_label == "non-T2DM")
hlp_t2dm   <- df_hlp_complete %>% filter(diabetes_label == "T2DM")

tbl_disease_status <- bind_rows(
  make_disease_table(met_nondm, hlp_nondm, "Non-T2DM"),
  make_disease_table(met_t2dm,  hlp_t2dm,  "T2DM")
)

write_xlsx(list("Table_Disease_Status" = tbl_disease_status),
           path = "Sup_data_table3_met_lipid_status.xlsx")


# ==============================================================================
# Supplementary Table 4 — Multivariable linear regression analysis
# Predictor of interest: % Change in Fasting Insulin (insulin_change_pct)
# Three models:
#   m_pct        : Y = weight_loss_pct    ~ age + gender + diabetes + total_energy + Δins%
#   m_muscle_pct : Y = muscle_change_pct  ~ age + gender + diabetes + total_energy + Δins%
#   m_log_fat_pct: Y = log(fat_w8)        ~ log(fat_w0) + age + gender + diabetes + total_energy + Δins%
# Output: Sup_table4_models.xlsx
# ==============================================================================

# ── Fit % Δinsulin models ──────────────────────────────────────────────────────
m_pct <- glm(weight_loss_pct ~ age + gender + diabetes + total_energy + insulin_change_pct,
              data = df_complete_full %>%
                filter(complete.cases(insulin_change_pct)),
              family = gaussian())

m_muscle_pct <- glm(muscle_change_pct ~ age + gender + diabetes + total_energy + insulin_change_pct,
                     data = df_complete_full %>%
                       filter(complete.cases(insulin_change_pct)),
                     family = gaussian())

m_log_fat_pct <- glm(log(body_fat_mass_w8) ~ log(body_fat_mass_w0) +
                        age + gender + diabetes + total_energy + insulin_change_pct,
                      data = df_complete_full %>%
                        filter(body_fat_mass_w8 > 0, body_fat_mass_w0 > 0,
                               complete.cases(total_energy, insulin_change_pct)),
                      family = gaussian())



# ── Variable labels for % Δinsulin models ─────────────────────────────────────
labels_m_XY_pct <- c(
  "age"                    = "Age",
  "gendermale"             = "Gender",
  "diabetesYes"            = "Diabetes",
  "total_energy"               = "Total energy intake (kcal/day)",
  "insulin_change_pct"= "Δ Fasting insulin (%)"
)
labels_m_log_fat_pct <- c(
  "log(body_fat_mass_w0)"  = "log(Baseline Fat Mass)",
  "age"                    = "Age",
  "gendermale"             = "Gender",
  "diabetesYes"            = "Diabetes",
  "total_energy"               = "Total energy intake (kcal/day)",
  "insulin_change_pct"= "Δ Fasting insulin (%)"
)

# ── Export Supplementary Table 4 ──────────────────────────────────────────────
# Columns: β, SE, 95% CI, p-value, AIC, BIC, Adj R², N,
#          standardized β, SE_std, 95% CI_std
write_xlsx(
  list(
    "m_pct_weight_loss"   = extract_glm_with_std(m_pct,
                               "m_pct: weight_loss_pct [% Δinsulin]"),
    "m_pct_muscle_change" = extract_glm_with_std(m_muscle_pct,
                               "m_pct_muscle: muscle_change_pct [% Δinsulin]"),
    "m_pct_log_fat_w8"    = extract_glm_with_std(m_log_fat_pct,
                               "m_pct_log_fat: log(fat_w8) [% Δinsulin]")
  ),
  path = "Sup_table4_models.xlsx"
)


# ==============================================================================
# Supplementary Table 5&6 — Multivariable linear regression analysis
#                           (dietary macronutrient composition)
# Models retained: M1 (protein g/kg/day + total energy intake) and
#                  M2 (carbohydrate % energy + fat % energy)
# Y variables:
#   (a) weight_loss_pct          (b) log(insulin_w8)  (c) log(homa_ir_w8)
#   (d) muscle_change_pct        (e) log(body_fat_mass_w8)
#   (f) muscle_change_pct  + insulin_change_pct
#   (g) log(fat_w8)        + insulin_change_pct
# Workflow: (1) pre-GLM skewness/SW on raw outcomes → (2) model fitting → (3) VIF
# Output: exploratory_GLM_results_20260512.xlsx, exploratory_outcome_SW.csv,
#         exploratory_sample_sizes.csv, exploratory_VIF.csv
# ==============================================================================

# ── Compute macronutrient energy fractions ──────────────────────────────────
df_complete_full <- df_complete_full %>%
  mutate(
    carb_pct_energy = carb_g * 4 / total_energy * 100,
    fat_pct_energy  = fat_g  * 9 / total_energy * 100
  )

# ── Complete-case analysis datasets (consistent N across M1 and M2) ──────────
df_W <- df_complete_full %>%
  filter(complete.cases(
    weight_loss_pct, age, gender, diabetes,
    weight_w0, weight_w8,               # weight_w8 explicit for log model
    carb_g, fat_g,
    protein_per_kg, total_energy,
    carb_pct_energy, fat_pct_energy
  ), weight_w0 > 0, weight_w8 > 0)

df_I <- df_complete_full %>%
  filter(complete.cases(
    age, gender, diabetes,
    bmi_w0, insulin_w0, insulin_w8,     # insulin_w8 explicit for log model
    carb_g, fat_g,
    protein_per_kg, total_energy,
    carb_pct_energy, fat_pct_energy
  ), insulin_w0 > 0, insulin_w8 > 0)

df_H <- df_complete_full %>%
  filter(complete.cases(
    age, gender, diabetes,
    bmi_w0, homa_ir_w0, homa_ir_w8,    # homa_ir_w8 explicit for log model
    carb_g, fat_g,
    protein_per_kg, total_energy,
    carb_pct_energy, fat_pct_energy
  ), homa_ir_w0 > 0, homa_ir_w8 > 0)

df_Mu <- df_complete_full %>%
  filter(complete.cases(
    age, gender, diabetes,
    bmi_w0, muscle_mass_w0, muscle_mass_w8,
    carb_g, fat_g,
    protein_per_kg, total_energy,
    carb_pct_energy, fat_pct_energy
  ), muscle_mass_w0 > 0, muscle_mass_w8 > 0)

df_Fa <- df_complete_full %>%
  filter(complete.cases(
    age, gender, diabetes,
    body_fat_mass_w0, body_fat_mass_w8, # bmi_w0 removed (collinear)
    carb_g, fat_g,
    protein_per_kg, total_energy,
    carb_pct_energy, fat_pct_energy
  ), body_fat_mass_w0 > 0, body_fat_mass_w8 > 0)

# ── Sample size summary ───────────────────────────────────────────────────────
sample_sizes <- data.frame(
  Y_original  = c("weight_loss_pct", "log(insulin_w8)",
                  "log(homa_ir_w8)", "log(muscle_mass_w8)", "log(body_fat_mass_w8)"),
  N           = c(nrow(df_W), nrow(df_I), nrow(df_H), nrow(df_Mu), nrow(df_Fa)),
  Baseline_covariate = c("log(weight_w0)", "log(insulin_w0)", "log(homa_ir_w0)",
                          "log(muscle_mass_w0)", "log(body_fat_mass_w0)"),
  BMI_in_model = c("No (weight_w0 used instead)",
                   "Yes (bmi_w0)", "Yes (bmi_w0)", "Yes (bmi_w0)",
                   "No (collinear with body_fat_mass_w0, removed)")
)
print(sample_sizes)
readr::write_csv(sample_sizes, "exploratory_sample_sizes.csv")

# ── Pre-GLM: Skewness & Shapiro-Wilk check on raw outcome variables ───────────
# Methods: "Outcomes with markedly skewed residuals (fasting insulin, HOMA-IR,
# and fat mass) were log-transformed prior to modeling."
# Checked on raw (untransformed) outcomes from complete-case datasets to inform
# log-transformation decisions before model fitting.
pre_glm_outcomes <- list(
  list(df_W$weight_loss_pct,         "weight_loss_pct        (→ no transform)"),
  list(df_I$insulin_w8,              "insulin_w8             (→ log-transformed)"),
  list(df_H$homa_ir_w8,              "homa_ir_w8             (→ log-transformed)"),
  list(df_Mu$muscle_change_pct,      "muscle_change_pct      (→ no transform)"),
  list(df_Fa$body_fat_mass_w8,       "body_fat_mass_w8       (→ log-transformed)")
)

pre_glm_sw <- purrr::map_dfr(pre_glm_outcomes, function(item) {
  x   <- item[[1]][!is.na(item[[1]])]
  n   <- length(x)
  set.seed(42)
  x_s  <- if (n > 5000) sample(x, 5000) else x
  sw   <- shapiro.test(x_s)
  skew <- mean(((x - mean(x)) / sd(x))^3)
  kurt <- mean(((x - mean(x)) / sd(x))^4) - 3
  data.frame(
    Outcome         = item[[2]],
    N               = n,
    SW_W            = round(sw$statistic, 4),
    SW_p            = ifelse(sw$p.value < 0.001, "<.001", sprintf("%.4f", sw$p.value)),
    Normality       = ifelse(sw$p.value > 0.05, "Pass", "Fail"),
    Skewness        = round(skew, 3),
    Excess_Kurtosis = round(kurt, 3),
    Skew_severity   = dplyr::case_when(
      abs(skew) < 0.5 ~ "Mild",
      abs(skew) < 1.0 ~ "Moderate",
      TRUE            ~ "Severe"
    )
  )
})

cat("\n=== Pre-GLM Outcome Distribution Check ===\n")
print(pre_glm_sw, row.names = FALSE)
readr::write_csv(pre_glm_sw, "exploratory_outcome_SW.csv")
cat("exploratory_outcome_SW.csv saved.\n")

# ── Fit M1 and M2 models ──────────────────────────────────────────────────────
# M1: y ~ covariates + protein_per_kg + total_energy
# M2: y ~ covariates + carb_pct_energy + fat_pct_energy
#
# Log-transformation rationale (Methods: "Outcomes with markedly skewed residuals
# [fasting insulin, HOMA-IR, and fat mass] were log-transformed prior to modeling,
# with corresponding baseline values included as log-transformed covariates"):
#   → log(insulin_w8)        | baseline: log(insulin_w0)
#   → log(homa_ir_w8)        | baseline: log(homa_ir_w0)
#   → log(body_fat_mass_w8)  | baseline: log(body_fat_mass_w0)
# Transformation decisions were based on residual skewness diagnostics run after
# an initial untransformed fit; see Section 6 (residual diagnostics) below.

# (a) Y = weight_loss_pct
mod_W_M1 <- glm(weight_loss_pct ~ age + gender + diabetes +
                  weight_w0 + protein_per_kg + total_energy,
                data = df_W, family = gaussian())
mod_W_M2 <- glm(weight_loss_pct ~ age + gender + diabetes +
                  weight_w0 + carb_pct_energy + fat_pct_energy,
                data = df_W, family = gaussian())

# (b) Y = log(fasting insulin W8)
mod_I_M1 <- glm(log(insulin_w8) ~ age + gender + diabetes +
                  bmi_w0 + log(insulin_w0) + weight_loss_pct + protein_per_kg + total_energy,
                data = df_I, family = gaussian())
mod_I_M2 <- glm(log(insulin_w8) ~ age + gender + diabetes +
                  bmi_w0 + log(insulin_w0) + weight_loss_pct + carb_pct_energy + fat_pct_energy,
                data = df_I, family = gaussian())

# (c) Y = log(HOMA-IR W8)
mod_H_M1 <- glm(log(homa_ir_w8) ~ age + gender + diabetes +
                  bmi_w0 + log(homa_ir_w0) + weight_loss_pct + protein_per_kg + total_energy,
                data = df_H, family = gaussian())
mod_H_M2 <- glm(log(homa_ir_w8) ~ age + gender + diabetes +
                  bmi_w0 + log(homa_ir_w0) + weight_loss_pct + carb_pct_energy + fat_pct_energy,
                data = df_H, family = gaussian())

# (d) Y = muscle_change_pct  [log(muscle_w8) excluded: |skew| = 1.2]
mod_Mu_orig_M1 <- glm(muscle_change_pct ~ age + gender + diabetes +
                        bmi_w0 + muscle_mass_w0 + weight_loss_pct + protein_per_kg + total_energy,
                      data = df_Mu, family = gaussian())
mod_Mu_orig_M2 <- glm(muscle_change_pct ~ age + gender + diabetes +
                        bmi_w0 + muscle_mass_w0 + weight_loss_pct + carb_pct_energy + fat_pct_energy,
                      data = df_Mu, family = gaussian())

# (e) Y = log(body_fat_mass_w8)
mod_Fa_M1 <- glm(log(body_fat_mass_w8) ~ age + gender + diabetes +
                   log(body_fat_mass_w0) + weight_loss_pct + protein_per_kg + total_energy,
                 data = df_Fa, family = gaussian())
mod_Fa_M2 <- glm(log(body_fat_mass_w8) ~ age + gender + diabetes +
                   log(body_fat_mass_w0) + weight_loss_pct + carb_pct_energy + fat_pct_energy,
                 data = df_Fa, family = gaussian())

cat("Core exploratory models (M1 + M2) fitted.\n")

# ── Insulin-predictor sub-datasets (% change only) ───────────────────────────
df_Mu_ins_pct <- df_Mu %>%
  filter(complete.cases(insulin_change_pct))
df_Fa_ins_pct <- df_Fa %>%
  filter(complete.cases(insulin_change_pct))

cat("N for insulin_change_pct models:",
    "Mu =", nrow(df_Mu_ins_pct), "| Fa =", nrow(df_Fa_ins_pct), "
")

# (h) Y = muscle_change_pct  + Δ fasting insulin (%)
mod_Mu_orig_InsPct_M1 <- glm(muscle_change_pct ~ age + gender + diabetes +
                                bmi_w0 + muscle_mass_w0 + weight_loss_pct +
                                insulin_change_pct + protein_per_kg + total_energy,
                              data = df_Mu_ins_pct, family = gaussian())
mod_Mu_orig_InsPct_M2 <- glm(muscle_change_pct ~ age + gender + diabetes +
                                bmi_w0 + muscle_mass_w0 + weight_loss_pct +
                                insulin_change_pct + carb_pct_energy + fat_pct_energy,
                              data = df_Mu_ins_pct, family = gaussian())

# (i) Y = log(body_fat_mass_w8)  + Δ fasting insulin (%)
mod_Fa_InsPct_M1 <- glm(log(body_fat_mass_w8) ~ age + gender + diabetes +
                           log(body_fat_mass_w0) + weight_loss_pct +
                           insulin_change_pct + protein_per_kg + total_energy,
                         data = df_Fa_ins_pct, family = gaussian())
mod_Fa_InsPct_M2 <- glm(log(body_fat_mass_w8) ~ age + gender + diabetes +
                           log(body_fat_mass_w0) + weight_loss_pct +
                           insulin_change_pct + carb_pct_energy + fat_pct_energy,
                         data = df_Fa_ins_pct, family = gaussian())

cat("InsPct models (h/i) fitted.\n")

# ── Collect all exploratory models for All_Models sheet ──────────────────────
expl_model_list <- list(
  list(mod_W_M1,              "W_M1:           weight_loss_pct ~ weight_w0 + protein_pkg + total energy intake"),
  list(mod_W_M2,              "W_M2:           weight_loss_pct ~ weight_w0 + carb% + fat%"),
  list(mod_I_M1,              "I_M1:           log(insulin_w8) ~ bmi + log(insulin_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_I_M2,              "I_M2:           log(insulin_w8) ~ bmi + log(insulin_w0) + wloss + carb% + fat%"),
  list(mod_H_M1,              "H_M1:           log(homa_ir_w8) ~ bmi + log(homa_ir_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_H_M2,              "H_M2:           log(homa_ir_w8) ~ bmi + log(homa_ir_w0) + wloss + carb% + fat%"),
  list(mod_Mu_orig_M1,        "Mu_orig_M1:     muscle_change_pct ~ bmi + muscle_w0 + wloss + protein_pkg + total energy intake"),
  list(mod_Mu_orig_M2,        "Mu_orig_M2:     muscle_change_pct ~ bmi + muscle_w0 + wloss + carb% + fat%"),
  list(mod_Fa_M1,             "Fa_M1:          log(fat_w8) ~ log(fat_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_Fa_M2,             "Fa_M2:          log(fat_w8) ~ log(fat_w0) + wloss + carb% + fat%"),
  list(mod_Mu_orig_InsPct_M1, "Mu_InsPct_M1:   muscle_change_pct + Δinsulin(%) ~ bmi + muscle_w0 + wloss + protein_pkg + total energy intake"),
  list(mod_Mu_orig_InsPct_M2, "Mu_InsPct_M2:   muscle_change_pct + Δinsulin(%) ~ bmi + muscle_w0 + wloss + carb% + fat%"),
  list(mod_Fa_InsPct_M1,      "Fa_InsPct_M1:   log(fat_w8) + Δinsulin(%) ~ log(fat_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_Fa_InsPct_M2,      "Fa_InsPct_M2:   log(fat_w8) + Δinsulin(%) ~ log(fat_w0) + wloss + carb% + fat%")
)

all_expl_results <- purrr::map_dfr(expl_model_list,
                                   ~ extract_glm(.x[[1]], .x[[2]]))

write_xlsx(
  list(
    "All_Models"           = all_expl_results,
    "W_weight_loss_pct"    = bind_rows(
      extract_glm(mod_W_M1,   "W_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_W_M2,   "W_M2: carb% + fat%")
    ),
    "I_log_insulin_w8"     = bind_rows(
      extract_glm(mod_I_M1,   "I_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_I_M2,   "I_M2: carb% + fat%")
    ),
    "H_log_homa_ir_w8"     = bind_rows(
      extract_glm(mod_H_M1,   "H_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_H_M2,   "H_M2: carb% + fat%")
    ),
    "Mu_muscle_change_pct" = bind_rows(
      extract_glm(mod_Mu_orig_M1, "Mu_orig_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_Mu_orig_M2, "Mu_orig_M2: carb% + fat%")
    ),
    "Fa_log_fat_w8"        = bind_rows(
      extract_glm(mod_Fa_M1,  "Fa_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_Fa_M2,  "Fa_M2: carb% + fat%")
    ),
    "Mu_orig_InsPct"       = bind_rows(
      extract_glm(mod_Mu_orig_InsPct_M1, "Mu_orig_InsPct_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_Mu_orig_InsPct_M2, "Mu_orig_InsPct_M2: carb% + fat%")
    ),
    "Fa_log_InsPct"        = bind_rows(
      extract_glm(mod_Fa_InsPct_M1, "Fa_InsPct_M1: protein_per_kg + total energy intake"),
      extract_glm(mod_Fa_InsPct_M2, "Fa_InsPct_M2: carb% + fat%")
    )
  ),
  path = "exploratory_GLM_results_20260512.xlsx"
)
cat("exploratory_GLM_results_20260512.xlsx saved.\n")

# ── VIF check (M1 + M2 models only) ─────────────────────────────────────────
# Methods: "Multicollinearity was assessed using variance inflation factors (VIFs)."
# Threshold: VIF < 5 acceptable; 5–10 moderate concern; > 10 high concern.
# Results reported in Methods; outputs saved to exploratory_VIF.csv.
extract_vif <- function(model, model_name) {
  v <- car::vif(model)
  if (is.matrix(v)) {
    vif_adj <- (v[, "GVIF^(1/(2*Df))"])^2
    vars     <- rownames(v)
  } else {
    vif_adj <- v
    vars     <- names(v)
  }
  data.frame(
    Model    = model_name,
    Variable = vars,
    VIF      = round(vif_adj, 3),
    Flag     = dplyr::case_when(
      vif_adj >= 10 ~ "HIGH (>10)",
      vif_adj >= 5  ~ "Moderate (5-10)",
      TRUE          ~ "OK (<5)"
    ),
    row.names = NULL
  )
}

all_expl_models <- list(
  # ── Primary models (Supplementary Table 4) ───────────────────────────────
  list(m_pct,              "Primary_m_pct:       weight_loss_pct ~ total_energy + insulin_change_pct"),
  list(m_muscle_pct,       "Primary_m_muscle:    muscle_change_pct ~ total_energy + insulin_change_pct"),
  list(m_log_fat_pct,      "Primary_m_log_fat:   log(fat_w8) ~ log(fat_w0) + total_energy + insulin_change_pct"),
  # ── Exploratory M1/M2 core models ────────────────────────────────────────
  list(mod_W_M1,           "W_M1:    weight_loss_pct ~ weight_w0 + protein_pkg + total energy intake"),
  list(mod_W_M2,           "W_M2:    weight_loss_pct ~ weight_w0 + carb% + fat%"),
  list(mod_I_M1,           "I_M1:    log(insulin_w8) ~ bmi + log(insulin_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_I_M2,           "I_M2:    log(insulin_w8) ~ bmi + log(insulin_w0) + wloss + carb% + fat%"),
  list(mod_H_M1,           "H_M1:    log(homa_ir_w8) ~ bmi + log(homa_ir_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_H_M2,           "H_M2:    log(homa_ir_w8) ~ bmi + log(homa_ir_w0) + wloss + carb% + fat%"),
  list(mod_Mu_orig_M1,     "Mu_orig_M1: muscle_change_pct ~ bmi + muscle_w0 + wloss + protein_pkg + total energy intake"),
  list(mod_Mu_orig_M2,     "Mu_orig_M2: muscle_change_pct ~ bmi + muscle_w0 + wloss + carb% + fat%"),
  list(mod_Fa_M1,          "Fa_M1:   log(fat_w8) ~ log(fat_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_Fa_M2,          "Fa_M2:   log(fat_w8) ~ log(fat_w0) + wloss + carb% + fat%"),
  # ── InsPct models (muscle / fat + Δinsulin %) ────────────────────────────
  list(mod_Mu_orig_InsPct_M1, "Mu_InsPct_M1: muscle_change_pct + Δinsulin(%) ~ bmi + muscle_w0 + wloss + protein_pkg + total energy intake"),
  list(mod_Mu_orig_InsPct_M2, "Mu_InsPct_M2: muscle_change_pct + Δinsulin(%) ~ bmi + muscle_w0 + wloss + carb% + fat%"),
  list(mod_Fa_InsPct_M1,      "Fa_InsPct_M1: log(fat_w8) + Δinsulin(%) ~ log(fat_w0) + wloss + protein_pkg + total energy intake"),
  list(mod_Fa_InsPct_M2,      "Fa_InsPct_M2: log(fat_w8) + Δinsulin(%) ~ log(fat_w0) + wloss + carb% + fat%")
)
vif_all <- purrr::map_dfr(all_expl_models,
                           ~ extract_vif(.x[[1]], .x[[2]]))

cat("\n=== VIF Summary (flagged variables only) ===\n")
vif_flagged <- vif_all %>% filter(Flag != "OK (<5)")
if (nrow(vif_flagged) > 0) {
  print(vif_flagged, row.names = FALSE)
} else {
  cat("All VIF values < 5. No collinearity concern.\n")
}

readr::write_csv(vif_all, "exploratory_VIF.csv")
cat("exploratory_VIF.csv saved.\n")

# ── Post-VIF alert summary ────────────────────────────────────────────────────
cat("\n=== Pre-GLM Skewness Alert ===\n")
cat("Severe   (|skew| >= 1.0):", sum(pre_glm_sw$Skew_severity == "Severe"),   "outcomes\n")
cat("Moderate (0.5–1.0)      :", sum(pre_glm_sw$Skew_severity == "Moderate"), "outcomes\n")
cat("Mild     (< 0.5)        :", sum(pre_glm_sw$Skew_severity == "Mild"),     "outcomes\n")

cat("\n=== VIF Alert ===\n")
cat("Variables with VIF >= 10:", sum(vif_all$Flag == "HIGH (>10)"), "\n")
cat("Variables with VIF 5-10 :", sum(vif_all$Flag == "Moderate (5-10)"), "\n")



# ==============================================================================
# Supplementary Table 7&8 — Changes in regional circumference of participants 
#                         in week 0 (baseline) and week 8
# 19 circumference sites × 4 subgroups
# (T2D×Male | T2D×Female | Non-T2D×Male | Non-T2D×Female)
# Note: Requires circ_base_vars (computed in Section I) and .circ_label_map /
#       .anatomical_levels (defined immediately below; also appear in Figure 5).
# Output: Ex_data_table7_3D_circumference_paired_comparison.xlsx
# ==============================================================================

# ── Circumference labels and anatomical ordering ──────────────────────────────
.circ_label_map <- c(
  "Chest.Circumference"           = "Chest",
  "Narrow.Waist.Circumference"    = "Narrow Waist",
  "Low.Hip.Circumference"         = "Low Hip",
  "High.Hip.Circumference"        = "High Hip",
  "Belly.Circumference"           = "Belly",
  "Left.Thigh.Circumference"      = "Left Thigh",
  "Right.Thigh.Circumference"     = "Right Thigh",
  "Left.Knee.Circumference"       = "Left Knee",
  "Right.Knee.Circumference"      = "Right Knee",
  "Left.Calf.Circumference"       = "Left Calf",
  "Right.Calf.Circumference"      = "Right Calf",
  "Left.Ankle.Circumference"      = "Left Ankle",
  "Right.Ankle.Circumference"     = "Right Ankle",
  "Neck.Base.Circumference"       = "Neck Base",
  "Collar.Circumference"          = "Collar",
  "Left.Upper.Arm.Circumference"  = "Left Upper Arm",
  "Right.Upper.Arm.Circumference" = "Right Upper Arm",
  "Left.Wrist.Circumference"      = "Left Wrist",
  "Right.Wrist.Circumference"     = "Right Wrist"
)

# Anatomical order: foot-to-head as factor levels
# (ggplot y-axis: level 1 = bottom → Collar ends up at top of chart)
.anatomical_levels <- c(
  "Right Ankle", "Left Ankle",
  "Right Calf",  "Left Calf",
  "Right Knee",  "Left Knee",
  "Right Thigh", "Left Thigh",
  "Low Hip",     "High Hip",
  "Belly",       "Narrow Waist",
  "Right Wrist", "Left Wrist",
  "Right Upper Arm", "Left Upper Arm",
  "Chest",
  "Neck Base",   "Collar"
)

# Build var_pairs and labels from circ_base_vars + .circ_label_map
# circ_base_vars uses dot-notation (e.g. "Chest.Circumference")
# Columns in df_complete_full: "<var>_w0" / "<var>_w8"
circ_s_pairs <- lapply(circ_base_vars, function(v) {
  c(paste0(v, "_w0"), paste0(v, "_w8"))
})
circ_s_labels <- dplyr::recode(circ_base_vars, !!!.circ_label_map, .default = circ_base_vars)

# Reorder pairs + labels by anatomical order (head-to-foot, matching lollipop4)
.circ_order_idx <- order(match(circ_s_labels, rev(.anatomical_levels)))
circ_s_pairs  <- circ_s_pairs[.circ_order_idx]
circ_s_labels <- circ_s_labels[.circ_order_idx]

tbl_circ_total      <- make_paired_table(df_complete_full, circ_s_pairs, circ_s_labels)
tbl_circ_t2d        <- make_paired_table(df_dm,            circ_s_pairs, circ_s_labels)
tbl_circ_nont2d     <- make_paired_table(df_nondm,         circ_s_pairs, circ_s_labels)

# DM × Gender subgroups
df_t2d_male    <- df_dm    %>% filter(tolower(as.character(gender)) == "male")
df_t2d_female  <- df_dm    %>% filter(tolower(as.character(gender)) == "female")
df_non_male    <- df_nondm %>% filter(tolower(as.character(gender)) == "male")
df_non_female  <- df_nondm %>% filter(tolower(as.character(gender)) == "female")

tbl_circ_t2d_male   <- make_paired_table(df_t2d_male,   circ_s_pairs, circ_s_labels)
tbl_circ_t2d_female <- make_paired_table(df_t2d_female, circ_s_pairs, circ_s_labels)
tbl_circ_non_male   <- make_paired_table(df_non_male,   circ_s_pairs, circ_s_labels)
tbl_circ_non_female <- make_paired_table(df_non_female, circ_s_pairs, circ_s_labels)

write_xlsx(
  list(
    "Total"          = tbl_circ_total,
    "T2D"            = tbl_circ_t2d,
    "Non-T2D"        = tbl_circ_nont2d,
    "T2D_Male"       = tbl_circ_t2d_male,
    "T2D_Female"     = tbl_circ_t2d_female,
    "Non-T2D_Male"   = tbl_circ_non_male,
    "Non-T2D_Female" = tbl_circ_non_female,
    "Notes"          = data.frame(Note = c(
      "Paired Wilcoxon signed-rank test (W0 vs W8).",
      "Delta = Week 8 − Baseline (negative = reduction).",
      "Effect size r = |Z| / sqrt(N); Z from Wilcoxon V statistic.",
      paste0("Total N = ",          nrow(df_complete_full %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("T2D N = ",            nrow(df_dm            %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("Non-T2D N = ",        nrow(df_nondm         %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("T2D Male N = ",       nrow(df_t2d_male      %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("T2D Female N = ",     nrow(df_t2d_female    %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("Non-T2D Male N = ",   nrow(df_non_male      %>% filter(!is.na(date_w0), !is.na(date_w8)))),
      paste0("Non-T2D Female N = ", nrow(df_non_female    %>% filter(!is.na(date_w0), !is.na(date_w8))))
    ))
  ),
  path = "Ex_data_table7_3D_circumference_paired_comparison.xlsx"
)


# ── 3D scan analysis dataset (shared by ET9 and Figure 5) ────────────────────
# Requires: date_w0, date_w8, insulin_change (non-NA); circ_base_vars (Section I)
# df_3d is created here so ET9 below can use it; Figure 5 references it directly.
df_3d <- df_complete_full %>%
  filter(!is.na(date_w0), !is.na(date_w8), !is.na(insulin_change)) %>%
  mutate(
    Gender = factor(
      ifelse(tolower(as.character(gender)) == "male", "Male", "Female"),
      levels = c("Male", "Female")
    ),
    diabetes_label = factor(as.character(diabetes_label),
                            levels = c("non-T2DM", "T2DM")),
    fill_group = factor(
      paste(as.character(diabetes_label), as.character(Gender)),
      levels = names(.gender_colors_by_dm)
    )
  )

cat("\ndf_3d (W0+W8+insulin_change complete):\n")
cat("  N clients  :", nrow(df_3d), "\n")
cat("  T2DM       :", sum(df_3d$diabetes_label == "T2DM",     na.rm = TRUE), "\n")
cat("  Non-T2DM   :", sum(df_3d$diabetes_label == "non-T2DM", na.rm = TRUE), "\n")
cat("  Male       :", sum(df_3d$Gender == "Male",   na.rm = TRUE), "\n")
cat("  Female     :", sum(df_3d$Gender == "Female", na.rm = TRUE), "\n")

# Δ circumferences (W8 − W0) — used by ET9 Spearman and Figure 5 lollipop/dumbbell
for (.v in circ_base_vars) {
  df_3d[[paste0("delta_", .v)]] <-
    df_3d[[paste0(.v, "_w8")]] - df_3d[[paste0(.v, "_w0")]]
}
delta_cols <- paste0("delta_", circ_base_vars)

# ==============================================================================
# Supplementary Table 9 — Spearman r: Δcircumference vs Δinsulin by DM group
# Circumference sites × DM status; Spearman r with BH-corrected p-value
# Output: lollipop1_data.csv
# ==============================================================================
# lollipop1_data is also used by Figure 5d/e; BH correction applied within each DM group.
.lollipop1_rows <- vector("list", 0L)
for (.dm in c("T2DM", "non-T2DM")) {
  .sub <- df_3d[as.character(df_3d$diabetes_label) == .dm, ]
  for (.v in circ_base_vars) {
    .dv  <- .sub[[paste0("delta_", .v)]]
    .ins <- .sub$insulin_change
    .ok  <- !is.na(.dv) & !is.na(.ins)
    if (sum(.ok) < 5L) next
    .ct <- cor.test(.dv[.ok], .ins[.ok], method = "spearman", exact = FALSE)
    .lollipop1_rows[[length(.lollipop1_rows) + 1L]] <- data.frame(
      diabetes_label = .dm,
      variable       = .v,
      r              = as.numeric(.ct$estimate),
      p_raw          = .ct$p.value,
      n              = sum(.ok),
      stringsAsFactors = FALSE
    )
  }
}
lollipop1_data <- do.call(rbind, .lollipop1_rows)

# BH correction within each DM group
lollipop1_data <- do.call(rbind, lapply(
  split(lollipop1_data, lollipop1_data$diabetes_label),
  function(x) {
    x$p_adj <- p.adjust(x$p_raw, method = "BH")
    x
  }
))

write_csv(lollipop1_data, "lollipop1_data.csv")


# ==============================================================================
# Figure 2 — Body Composition and Metabolic Marker Changes (W0 vs W8)
# Panels:
#   Fig2a: Body weight change — Non-T2DM
#   Fig2b–e: Body composition (Fat mass %, Muscle mass %, BMI, Waist) — Non-T2DM
#   Fig2f: Body weight change — T2DM
#   Fig2g–j: Body composition (Fat mass %, Muscle mass %, BMI, Waist) — T2DM
#   Fig2k: Metabolic markers (Glucose, Insulin, HOMA-IR) — Non-T2DM (3 files)
#   Fig2l: Metabolic markers (Glucose, Insulin, HOMA-IR) — T2DM (3 files)
#   Fig2m: Disease reclassification (glycaemic status, metabolic syndrome, HLP)
# Note: DM group split uses diabetes_label ("non-T2DM" / "T2DM").
# Output: Fig2a–l_*.pdf (16 panel files), Fig2m.pdf
# ==============================================================================

# ── Fig2a–l: Body weight, body composition + metabolic markers, by DM group ───
# Specs list ordered Fig2a (Non-T2DM weight), Fig2b → Fig2e (Non-T2DM body comp),
# Fig2f (T2DM weight), Fig2g → Fig2j (T2DM body comp),
# Fig2k (Non-T2DM metabolic ×3), Fig2l (T2DM metabolic ×3).
.fig2_main_specs <- list(
  # --- Fig2a: Non-T2DM Weight ---
  list(panel = "a", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "weight_w0", var_w8 = "weight_w8",
       ylabel = "Weight (kg)", suffix = "weight", ylim = c(35, 135)),
  # --- Fig2b: Non-T2DM Fat mass ---
  list(panel = "b", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "body_fat_mass_percentage_w0", var_w8 = "body_fat_mass_percentage_w8",
       ylabel = "Fat Mass (%)", suffix = "fat_pct", ylim = c(10, 65)),
  # --- Fig2c: Non-T2DM Muscle mass ---
  list(panel = "c", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "muscle_mass_percentage_w0", var_w8 = "muscle_mass_percentage_w8",
       ylabel = "Muscle Mass (%)", suffix = "muscle_pct", ylim = c(18, 52)),
  # --- Fig2d: Non-T2DM BMI ---
  list(panel = "d", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "bmi_w0", var_w8 = "bmi_w8",
       ylabel = "BMI (kg/m²)", suffix = "bmi", ylim = c(15, 52)),
  # --- Fig2e: Non-T2DM Waist ---
  list(panel = "e", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "waist_circumference_w0", var_w8 = "waist_circumference_w8",
       ylabel = "Waist Circumference (cm)", suffix = "waist", ylim = c(50, 155)),
  # --- Fig2f: T2DM Weight ---
  list(panel = "f", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "weight_w0", var_w8 = "weight_w8",
       ylabel = "Weight (kg)", suffix = "weight", ylim = c(35, 135)),
  # --- Fig2g: T2DM Fat mass ---
  list(panel = "g", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "body_fat_mass_percentage_w0", var_w8 = "body_fat_mass_percentage_w8",
       ylabel = "Fat Mass (%)", suffix = "fat_pct", ylim = c(10, 65)),
  # --- Fig2h: T2DM Muscle mass ---
  list(panel = "h", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "muscle_mass_percentage_w0", var_w8 = "muscle_mass_percentage_w8",
       ylabel = "Muscle Mass (%)", suffix = "muscle_pct", ylim = c(18, 52)),
  # --- Fig2i: T2DM BMI ---
  list(panel = "i", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "bmi_w0", var_w8 = "bmi_w8",
       ylabel = "BMI (kg/m²)", suffix = "bmi", ylim = c(15, 52)),
  # --- Fig2j: T2DM Waist ---
  list(panel = "j", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "waist_circumference_w0", var_w8 = "waist_circumference_w8",
       ylabel = "Waist Circumference (cm)", suffix = "waist", ylim = c(50, 155)),
  # --- Fig2k: Non-T2DM metabolic markers (3 files) ---
  list(panel = "k", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "glucose_ac_w0", var_w8 = "glucose_ac_w8",
       ylabel = "Fasting Glucose (mg/dL)", suffix = "glucose", ylim = c(60, 110)),
  list(panel = "k", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "insulin_w0", var_w8 = "insulin_w8",
       ylabel = "Fasting Insulin (μU/mL)", suffix = "insulin", ylim = c(0, 25)),
  list(panel = "k", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "homa_ir_w0", var_w8 = "homa_ir_w8",
       ylabel = "HOMA-IR", suffix = "homa_ir", ylim = c(0, 6)),
  # --- Fig2l: T2DM metabolic markers (3 files) ---
  list(panel = "l", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "glucose_ac_w0", var_w8 = "glucose_ac_w8",
       ylabel = "Fasting Glucose (mg/dL)", suffix = "glucose", ylim = c(55, 200)),
  list(panel = "l", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "insulin_w0", var_w8 = "insulin_w8",
       ylabel = "Fasting Insulin (μU/mL)", suffix = "insulin", ylim = c(0, 40)),
  list(panel = "l", dm_group = "T2DM", dm_tag = "T2D",
       var_w0 = "homa_ir_w0", var_w8 = "homa_ir_w8",
       ylabel = "HOMA-IR", suffix = "homa_ir", ylim = c(0, 10))
)

for (.s in .fig2_main_specs) {
  .p <- make_violin_single(
    df_complete_full,
    .s$var_w0, .s$var_w8,
    .s$ylabel,
    .s$dm_group,
    ylim = .s$ylim
  ) +
    labs(title = paste0("Fig2", .s$panel, ". ", .s$ylabel,
                        " — ", .s$dm_group, " (Baseline vs Week 8)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

  .fname <- paste0("Fig2", .s$panel, "_", .s$suffix, "_", .s$dm_tag, ".pdf")
  ggsave(.fname, plot = .p, width = 7, height = 5, bg = "white")
  cat("Saved:", .fname, "\n")
}
cat("Figures Fig2a–l saved.\n")


make_plot_dm_split <- function(df, col_w0, col_w8,
                               levels_order, fill_values,
                               legend_breaks, title, footnote) {
  # Use diabetes_label (derived from diabetes variable) for facet grouping
  df_plot <- df %>%
    filter(!is.na(diabetes_label), !is.na(.data[[col_w0]]) | !is.na(.data[[col_w8]])) %>%
    mutate(
      diabetes_label = factor(diabetes_label, levels = c("non-T2DM", "T2DM"))
    )

  bind_rows(
    df_plot %>% filter(!is.na(.data[[col_w0]])) %>%
      count(diabetes_label, group = .data[[col_w0]]) %>%
      mutate(timepoint = "Week 0"),
    df_plot %>% filter(!is.na(.data[[col_w8]])) %>%
      count(diabetes_label, group = .data[[col_w8]]) %>%
      mutate(timepoint = "Week 8")
  ) %>%
    mutate(
      timepoint = factor(timepoint, levels = c("Week 8", "Week 0")),
      group     = factor(group, levels = levels_order)
    ) %>%
    group_by(timepoint, diabetes_label) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup() %>%
    ggplot(aes(x = pct, y = timepoint, fill = group)) +
    geom_col(width = 0.55) +
    # Add % labels inside each bar segment (hide if segment too small)
    geom_text(
      aes(label = ifelse(pct >= 0.05,
                         scales::percent(pct, accuracy = 0.1), "")),
      position = position_stack(vjust = 0.5),
      color = "white", size = 3.5, fontface = "bold"
    ) +
    facet_wrap(~ diabetes_label, ncol = 2) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      breaks = seq(0, 1, by = 0.25),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_fill_manual(values = fill_values, breaks = legend_breaks) +
    labs(title = title, x = NULL, y = NULL,
         fill    = "Clinical criteria",
         caption = footnote) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold"),
      plot.caption       = element_text(hjust = 0, size = 8, color = "gray40"),
      legend.position    = "right",
      legend.title       = element_text(face = "bold", size = 11),
      strip.text         = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 11),
      panel.spacing.x    = unit(1, "cm")   # increase gap between non-T2DM / T2DM facets
    )
}

p_met <- make_plot_dm_split(
  df_met_complete, "met_gp", "met_gp_w8",
  levels_order  = c("Metabolic syndrome", "Normal"),
  fill_values   = c("Metabolic syndrome" = "#B22222", "Normal" = "#034078"),
  legend_breaks = c("Normal", "Metabolic syndrome"),
  title    = "Metabolic syndrome status",
  footnote = ""   # footnote_met
)

p_hlp <- make_plot_dm_split(
  df_hlp_complete, "hlp_gp", "hlp_gp_w8",
  levels_order  = c("Hyperlipidemia", "Normal"),
  fill_values   = c("Hyperlipidemia" = "#B22222", "Normal" = "#034078"),
  legend_breaks = c("Normal", "Hyperlipidemia"),
  title    = "Lipid status",
  footnote = ""   # footnote_hlp
)

figure_2m <- cowplot::plot_grid(p_met, p_hlp,
                                ncol = 1, align = "v", axis = "lr")

ggsave("Fig2m.pdf", plot = figure_2m,
       width = 12, height = 8, bg = "white")
cat("Saved: Fig2m.pdf\n")


# ==============================================================================
# Extended Data Figure 1 — Lipid Panel Changes (W0 vs W8)
# Metrics: Total Cholesterol (TC), HDL-C, LDL-C, Triglyceride
# Stratified by T2DM / non-T2DM
# Output: ExtDataFig1a–h_*.pdf (8 panel files)
# ==============================================================================

.extdata1_specs <- list(
  list(panel = "a", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "total_cholesterol_w0",              var_w8 = "total_cholesterol_w8",
       ylabel = "Total Cholesterol (mg/dL)", suffix = "tc",  ylim = c(100, 350)),
  list(panel = "b", dm_group = "T2DM",     dm_tag = "T2D",
       var_w0 = "total_cholesterol_w0",              var_w8 = "total_cholesterol_w8",
       ylabel = "Total Cholesterol (mg/dL)", suffix = "tc",  ylim = c(80, 330)),
  list(panel = "c", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "high_density_lipoprotein_w0",       var_w8 = "high_density_lipoprotein_w8",
       ylabel = "HDL-C (mg/dL)",             suffix = "hdl", ylim = c(25, 100)),
  list(panel = "d", dm_group = "T2DM",     dm_tag = "T2D",
       var_w0 = "high_density_lipoprotein_w0",       var_w8 = "high_density_lipoprotein_w8",
       ylabel = "HDL-C (mg/dL)",             suffix = "hdl", ylim = c(20, 95)),
  list(panel = "e", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "low_density_lipoprotein_cholesterol_w0",
       var_w8 = "low_density_lipoprotein_cholesterol_w8",
       ylabel = "LDL-C (mg/dL)",             suffix = "ldl", ylim = c(40, 270)),
  list(panel = "f", dm_group = "T2DM",     dm_tag = "T2D",
       var_w0 = "low_density_lipoprotein_cholesterol_w0",
       var_w8 = "low_density_lipoprotein_cholesterol_w8",
       ylabel = "LDL-C (mg/dL)",             suffix = "ldl", ylim = c(50, 270)),
  list(panel = "g", dm_group = "non-T2DM", dm_tag = "nonT2D",
       var_w0 = "triglyceride_w0",                   var_w8 = "triglyceride_w8",
       ylabel = "Triglyceride (mg/dL)",      suffix = "tg",  ylim = c(0, 350)),
  list(panel = "h", dm_group = "T2DM",     dm_tag = "T2D",
       var_w0 = "triglyceride_w0",                   var_w8 = "triglyceride_w8",
       ylabel = "Triglyceride (mg/dL)",      suffix = "tg",  ylim = c(0, 350))
)

for (.s in .extdata1_specs) {
  .p <- make_violin_single(
    df_complete_full,
    .s$var_w0, .s$var_w8,
    .s$ylabel,
    .s$dm_group,
    ylim = .s$ylim
  ) +
    labs(title = paste0("Extended Data Fig1", .s$panel, ". ", .s$ylabel,
                        " — ", .s$dm_group, " (Baseline vs Week 8)")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

  .fname <- paste0("ExtDataFig1", .s$panel, "_", .s$suffix, "_", .s$dm_tag, ".pdf")
  ggsave(.fname, plot = .p, width = 7, height = 5, bg = "white")
  cat("Saved:", .fname, "\n")
}
cat("Extended Data Figure 1 (lipid panel) saved.\n")


# ==============================================================================
# Figure 3 — Forest Plots: GLM Models (m_pct), Standardized β only
# Predictor of interest: % Change in Fasting Insulin (insulin_change_pct)
# Y variables and model formulas:
#   Fig3a (m_pct):
#     Y = weight_loss_pct ~ age + gender + diabetes + total_energy + insulin_change_pct
#   Fig3b (m_muscle_pct):
#     Y = muscle_change_pct ~ age + gender + diabetes + total_energy + insulin_change_pct
#   Fig3c (m_log_fat_pct):
#     Y = log(body_fat_mass_w8) ~ log(body_fat_mass_w0) + age + gender + diabetes
#         + total_energy + insulin_change_pct
# Note: Depends on m_pct, m_muscle_pct, m_log_fat_pct from Supplementary Table 4.
#       Raw β plots removed; standardized β only.
# Output: Fig3a_weight_loss.pdf, Fig3b_muscle.pdf, Fig3c_log_fat.pdf
# ==============================================================================

# ── Standardized β forest plots — % Δinsulin ─────────────────────────────────
p_m_pct_std         <- plot_forest_std(prepare_forest_std(m_pct,        labels_m_XY_pct),
                                        "Fig3a. Y = Weight Loss %  [Δ Insulin %]  — Standardized β")
p_m_muscle_pct_std  <- plot_forest_std(prepare_forest_std(m_muscle_pct,  labels_m_XY_pct),
                                        "Fig3b. Y = Muscle Change %  [Δ Insulin %]  — Standardized β")
p_m_log_fat_pct_std <- plot_forest_std(prepare_forest_std(m_log_fat_pct, labels_m_log_fat_pct),
                                        "Fig3c. Y = log(Fat Mass W8)  [Δ Insulin %]  — Standardized β")

ggsave("Fig3a_weight_loss.pdf", p_m_pct_std,        width=6, height=4,   bg="white")
ggsave("Fig3b_muscle.pdf",      p_m_muscle_pct_std,  width=6, height=4,   bg="white")
ggsave("Fig3c_log_fat.pdf",     p_m_log_fat_pct_std, width=6, height=4.5, bg="white")
cat("Figures Fig3a–c (standardized β, % Δinsulin) saved.\n")


# ==============================================================================
# Figure 4 + Extended Data Figure 2 — Exploratory Forest Plots:
#   Macronutrient Composition Effects
# Models: M1 (protein g/kg/day + total energy intake) and
#         M2 (carbohydrate % energy + fat % energy)
# Y variables and model formulas:
#   Fig4a/b (log insulin W8, M1/M2):
#     Y = log(insulin_w8) ~ age + gender + diabetes + bmi_w0 + log(insulin_w0)
#         + weight_loss_pct + [M1: protein_per_kg + total_energy |
#                              M2: carb_pct_energy + fat_pct_energy]
#   Fig4c/d (log HOMA-IR W8, M1/M2):
#     Y = log(homa_ir_w8) ~ age + gender + diabetes + bmi_w0 + log(homa_ir_w0)
#         + weight_loss_pct + [M1 / M2 dietary]
#   Fig4e/f (muscle_change_pct, M1/M2):
#     Y = muscle_change_pct ~ age + gender + diabetes + bmi_w0 + muscle_mass_w0
#         + weight_loss_pct + [M1 / M2 dietary]
#   Fig4g/h (log fat mass W8, M1/M2):
#     Y = log(body_fat_mass_w8) ~ age + gender + diabetes + log(body_fat_mass_w0)
#         + weight_loss_pct + [M1 / M2 dietary]
# Note: Depends on models fitted in Supplementary Table 5&6 section.
#       Weight loss models (W M1/M2) and abs-insulin panels removed.
# Output: Fig4a–h_*.pdf (8 files), ExtDataFig2a–d_*.pdf (4 files)
# ==============================================================================

# ── Shared labels and forest helper ──────────────────────────────────────────
lbl_demo   <- c("age"          = "Age",
                "gendermale"   = "Gender",
                "diabetesYes"  = "Diabetes")
lbl_wloss  <- c("weight_loss_pct" = "Weight Loss (%)")

# Dietary labels by model type
lbl_M1 <- c("protein_per_kg"  = "Protein (g/kg/day)",
             "total_energy"  = "Total energy intake (kcal/day)")
lbl_M2 <- c("carb_pct_energy" = "Carbohydrate (% energy)",
             "fat_pct_energy"  = "Fat (% energy)")

# Baseline covariates by Y
lbl_base_I       <- c("bmi_w0"               = "Baseline BMI",
                       "log(insulin_w0)"      = "log(Baseline Insulin)")
lbl_base_H       <- c("bmi_w0"               = "Baseline BMI",
                       "log(homa_ir_w0)"      = "log(Baseline HOMA-IR)")
lbl_base_Muorig  <- c("bmi_w0"               = "Baseline BMI",
                       "muscle_mass_w0"        = "Baseline Muscle (kg)")
lbl_base_Fa      <- c("log(body_fat_mass_w0)" = "log(Baseline Fat)")

# ── Helper: single forest panel ─────────────────────────────────
make_expl_forest <- function(model, var_labels, title) {
  plot_forest(prepare_forest(model, var_labels), title)
}

# ── Fig4a/b: Y = log(Insulin W8) — M1 and M2 ─────────────────────────────────
p_I_M1 <- make_expl_forest(mod_I_M1, c(lbl_demo, lbl_base_I, lbl_wloss, lbl_M1),
                            "Fig4a. Y = log(Insulin W8) | M1: Protein (g/kg) + Total energy intake")
p_I_M2 <- make_expl_forest(mod_I_M2, c(lbl_demo, lbl_base_I, lbl_wloss, lbl_M2),
                            "Fig4b. Y = log(Insulin W8) | M2: Carb% + Fat%")
ggsave("Fig4a_log_insulin_M1.pdf", p_I_M1, width=7, height=6, bg="white")
ggsave("Fig4b_log_insulin_M2.pdf", p_I_M2, width=7, height=6, bg="white")
cat("Fig4a/b log(Insulin W8) — M1/M2 saved.\n")

# ── Fig4c/d: Y = log(HOMA-IR W8) — M1 and M2 ────────────────────────────────
p_H_M1 <- make_expl_forest(mod_H_M1, c(lbl_demo, lbl_base_H, lbl_wloss, lbl_M1),
                            "Fig4c. Y = log(HOMA-IR W8) | M1: Protein (g/kg) + Total energy intake")
p_H_M2 <- make_expl_forest(mod_H_M2, c(lbl_demo, lbl_base_H, lbl_wloss, lbl_M2),
                            "Fig4d. Y = log(HOMA-IR W8) | M2: Carb% + Fat%")
ggsave("Fig4c_log_homa_ir_M1.pdf", p_H_M1, width=7, height=6, bg="white")
ggsave("Fig4d_log_homa_ir_M2.pdf", p_H_M2, width=7, height=6, bg="white")
cat("Fig4c/d log(HOMA-IR W8) — M1/M2 saved.\n")

# ── Fig4e/f: Y = log(Fat Mass W8) — M1 and M2 ────────────────────────────────
p_Fa_M1 <- make_expl_forest(mod_Fa_M1, c(lbl_demo, lbl_base_Fa, lbl_wloss, lbl_M1),
                            "Fig4e. Y = log(Fat W8) | M1: Protein (g/kg) + Total energy intake")
p_Fa_M2 <- make_expl_forest(mod_Fa_M2, c(lbl_demo, lbl_base_Fa, lbl_wloss, lbl_M2),
                            "Fig4f. Y = log(Fat W8) | M2: Carb% + Fat%")
ggsave("Fig4e_log_fat_M1.pdf", p_Fa_M1, width=7, height=6, bg="white")
ggsave("Fig4f_log_fat_M2.pdf", p_Fa_M2, width=7, height=6, bg="white")
cat("Fig4e/f log(Fat W8) — M1/M2 saved.\n")

# ── Fig4g/h: Y = Muscle Change % — M1 and M2 ─────────────────────────────────
p_Muorig_M1 <- make_expl_forest(mod_Mu_orig_M1, c(lbl_demo, lbl_base_Muorig, lbl_wloss, lbl_M1),
                                 "Fig4g. Y = Muscle Change % | M1: Protein (g/kg) + Total energy intake")
p_Muorig_M2 <- make_expl_forest(mod_Mu_orig_M2, c(lbl_demo, lbl_base_Muorig, lbl_wloss, lbl_M2),
                                 "Fig4h. Y = Muscle Change % | M2: Carb% + Fat%")
ggsave("Fig4g_muscle_pct_M1.pdf", p_Muorig_M1, width=7, height=6, bg="white")
ggsave("Fig4h_muscle_pct_M2.pdf", p_Muorig_M2, width=7, height=6, bg="white")
cat("Fig4g/h Muscle Change % — M1/M2 saved.\n")



# ── Insulin predictor label ───────────────────────────────────────────────────
lbl_ins_pct <- c("insulin_change_pct" = "Δ Fasting insulin (%)")

# ── Extended Data Fig2a/b: Y = Muscle Change % + Δinsulin (%) — M1 and M2 ────────
p_Muorig_InsPct_M1 <- make_expl_forest(
  mod_Mu_orig_InsPct_M1,
  c(lbl_demo, lbl_base_Muorig, lbl_wloss, lbl_ins_pct, lbl_M1),
  "ExtDataFig2a. Y = Muscle Change % + ΔInsulin (%) | M1: Protein (g/kg) + Total energy intake"
)
p_Muorig_InsPct_M2 <- make_expl_forest(
  mod_Mu_orig_InsPct_M2,
  c(lbl_demo, lbl_base_Muorig, lbl_wloss, lbl_ins_pct, lbl_M2),
  "ExtDataFig2b. Y = Muscle Change % + ΔInsulin (%) | M2: Carb% + Fat%"
)
ggsave("ExtDataFig2a_muscle_pct_InsPct_M1.pdf", p_Muorig_InsPct_M1, width=7, height=7, bg="white")
ggsave("ExtDataFig2b_muscle_pct_InsPct_M2.pdf", p_Muorig_InsPct_M2, width=7, height=7, bg="white")
cat("Extended Data Fig2a/b Muscle + ΔIns% — M1/M2 saved.\n")

# ── Extended Data Fig2c/d: Y = log(Fat W8) + Δinsulin (%) — M1 and M2 ──────────
p_Fa_InsPct_M1 <- make_expl_forest(
  mod_Fa_InsPct_M1,
  c(lbl_demo, lbl_base_Fa, lbl_wloss, lbl_ins_pct, lbl_M1),
  "ExtDataFig2c. Y = log(Fat W8) + ΔInsulin (%) | M1: Protein (g/kg) + Total energy intake"
)
p_Fa_InsPct_M2 <- make_expl_forest(
  mod_Fa_InsPct_M2,
  c(lbl_demo, lbl_base_Fa, lbl_wloss, lbl_ins_pct, lbl_M2),
  "ExtDataFig2d. Y = log(Fat W8) + ΔInsulin (%) | M2: Carb% + Fat%"
)
ggsave("ExtDataFig2c_log_fat_InsPct_M1.pdf", p_Fa_InsPct_M1, width=7, height=7, bg="white")
ggsave("ExtDataFig2d_log_fat_InsPct_M2.pdf", p_Fa_InsPct_M2, width=7, height=7, bg="white")
cat("Extended Data Fig2c/d log(Fat) + ΔIns% — M1/M2 saved.\n")


# ==============================================================================
# Figure 5 — 3D Body Shape Trajectory and Circumference Change
# Panels:
#   Fig5b: UMAP spotlight — Non-T2DM and T2DM highlighted (2 files)
#   Fig5c: Lollipop grid — median Δcircumference by DM × Sex (4 panels)
#   Fig5d: Lollipop — Spearman r, Δcircumference vs Δinsulin by DM group
#   Fig5e: Dumbbell plot — individual circumference W0→W8 vs Δinsulin (top 5 sites)
# Note: circ_base_vars, .circ_label_map, .anatomical_levels were already
#       defined in Section I and Supplementary Table 7&8 respectively;
#       the UMAP section below recomputes circ_base_vars from body_3d_wide
#       (= df_complete_full proxy) as a cross-check.
# Output: Fig5b_umap_spotlight_nont2d.pdf, Fig5b_umap_spotlight_t2d.pdf,
#         Fig5c_median_delta_circ_by_dm_sex.pdf, Fig5d_corr_circ_insulin.pdf,
#         Fig5e_dumbbell_circ_insulin.pdf
# ==============================================================================

# ── 3D circumference variables and UMAP data pipeline ────────────────────────
# df_3d, delta_cols, circ_base_vars — all pre-built in the Supplementary Table 9
# section above (requires date_w0, date_w8, insulin_change to be non-NA).
# circ_base_vars cross-check: body_3d_wide == df_complete_full (Section I, line 64).
cat("\n3D circumference variables identified:", length(circ_base_vars), "\n")
cat(paste(circ_base_vars, collapse = "\n"), "\n")

# ── UMAP embedding (used by Fig5b spotlight plots) ───────────────────────────
# W0 and W8 circumference profiles embedded JOINTLY in one UMAP space.
# Each client: grey dot (W0) → coloured dot/arrow (W8), by DM group.
# Colour scheme: non-T2DM = blue (#4E79A7), T2DM = red-orange (#E15759).
# Per-DM-group centroid arrows added (thick, same DM colour).
# NOTE on insulin_change: NEGATIVE values = improvement
#   (fasting insulin decreased → insulin_change = W8 − W0 < 0 = better)

# ── Exclude clients with any missing circumference value (W0 or W8) ──────────
# ~2% of clients have ≥1 NA circumference cell due to scanner detection failures.
# Exclusion is preferred over imputation for transparency; 2% loss is negligible.
.mat_w0_raw <- as.matrix(df_3d[, paste0(circ_base_vars, "_w0")])
.mat_w8_raw <- as.matrix(df_3d[, paste0(circ_base_vars, "_w8")])
.traj_ok    <- complete.cases(.mat_w0_raw) & complete.cases(.mat_w8_raw)

cat(sprintf("\nTrajectory UMAP — complete circumference data: %d / %d clients\n",
            sum(.traj_ok), nrow(df_3d)))
cat(sprintf("  Excluded due to ≥1 missing circ value: %d (%.1f%%)\n",
            sum(!.traj_ok), 100 * mean(!.traj_ok)))

df_traj   <- df_3d[.traj_ok, ]
.mat_w0ok <- .mat_w0_raw[.traj_ok, ]
.mat_w8ok <- .mat_w8_raw[.traj_ok, ]
n_traj    <- nrow(df_traj)

# ── Joint scaling: pool W0+W8 together, then scale ───────────────────────
# This preserves the magnitude relationship between timepoints in UMAP space.
.mat_pool <- rbind(.mat_w0ok, .mat_w8ok)
.pmeans   <- colMeans(.mat_pool, na.rm = TRUE)
.psds     <- apply(.mat_pool, 2, sd, na.rm = TRUE)
.psds[.psds == 0] <- 1                            # guard against zero-SD cols
.mat_sc   <- sweep(sweep(.mat_pool, 2, .pmeans, "-"), 2, .psds, "/")

# ── UMAP on the pooled (2n × p) matrix ───────────────────────────────────
if (n_traj < 3L)
  stop("Figure 5 UMAP requires ≥ 3 clients with complete circumference scans ",
       "(current n = ", n_traj, "). ",
       "This section is designed for the full dataset, not sample data.")

set.seed(42)
.emb <- umap(.mat_sc,
             n_neighbors = min(15L, nrow(.mat_sc) - 1L),
             min_dist    = 0.1,
             metric      = "euclidean",
             verbose     = FALSE)

# ── Build trajectory dataframe ────────────────────────────────────────────
df_traj_plot <- data.frame(
  client_id           = df_traj$client_id,
  diabetes_label      = factor(
    dplyr::recode(as.character(df_traj$diabetes_label),
                  "non-T2DM" = "Non-T2D", "T2DM" = "T2D"),
    levels = c("Non-T2D", "T2D")),
  Gender              = as.character(df_traj$Gender),
  fill_group          = df_traj$fill_group,
  insulin_change = df_traj$insulin_change,
  x0 = .emb[seq_len(n_traj),          1],   # W0 UMAP coordinates
  y0 = .emb[seq_len(n_traj),          2],
  x1 = .emb[seq_len(n_traj) + n_traj, 1],   # W8 UMAP coordinates
  y1 = .emb[seq_len(n_traj) + n_traj, 2],
  stringsAsFactors = FALSE
)

# ── Group-mean centroid arrows ────────────────────────────────────────────
# Overall centroid (thick black arrow, one per plot)
df_gm_overall <- data.frame(
  x0 = mean(df_traj_plot$x0), y0 = mean(df_traj_plot$y0),
  x1 = mean(df_traj_plot$x1), y1 = mean(df_traj_plot$y1)
)

# Per-DM-group centroids (used in faceted version)
df_gm_dm <- do.call(rbind, lapply(
  split(df_traj_plot, df_traj_plot$diabetes_label),
  function(d) data.frame(
    diabetes_label = unique(d$diabetes_label),
    x0 = mean(d$x0), y0 = mean(d$y0),
    x1 = mean(d$x1), y1 = mean(d$y1)
  )
))

# ── Mean Δ circumference per client (colour variable for arrows) ─────────
# mean_delta_circ = mean of all (W8−W0) circumference deltas per client
# Negative = body overall shrank (good), Positive = overall expanded
df_traj_plot$mean_delta_circ <- rowMeans(
  as.matrix(df_traj[, delta_cols]),
  na.rm = TRUE
)

# Symmetric colour scale limit (rounded up to nearest integer)
.circ_clim <- ceiling(max(abs(df_traj_plot$mean_delta_circ), na.rm = TRUE))

# Attach group-mean delta to df_gm_dm so centroid arrow can share gradient scale
df_gm_dm <- df_gm_dm %>%
  left_join(
    df_traj_plot %>%
      group_by(diabetes_label) %>%
      summarise(mean_delta_circ = mean(mean_delta_circ, na.rm = TRUE),
                .groups = "drop"),
    by = "diabetes_label"
  )

# Annotation stats
n_circ_dec   <- sum(df_traj_plot$mean_delta_circ < 0, na.rm = TRUE)
pct_circ_dec <- round(100 * n_circ_dec / n_traj)

# DM group colours used for dots and centroid arrows
.dm_dot_colors <- c("Non-T2D" = "#E69F00", "T2D" = "#CC79A7")

# ── KDE density grid for background contours ─────────────────────────────
.kde_grid <- MASS::kde2d(
  c(df_traj_plot$x0, df_traj_plot$x1),
  c(df_traj_plot$y0, df_traj_plot$y1),
  n = 80
)
df_kde    <- expand.grid(x = .kde_grid$x, y = .kde_grid$y)
df_kde$z  <- as.vector(.kde_grid$z)
n_kde_lev <- 6
.kde_brks <- quantile(df_kde$z[df_kde$z > 0],
                      probs = seq(0, 1, length.out = n_kde_lev + 1))


# ── Figure 5B: UMAP Spotlight plots ─────────────────────────────────────────
# Two plots: one T2D highlighted, one Non-T2D highlighted.
# Arrow colour = mean Δ circumference (highlighted group only).
# Background group rendered in grey at low alpha.

.make_umap_spotlight <- function(highlight_group, dot_color) {

  df_fg    <- df_traj_plot %>% filter(as.character(diabetes_label) == highlight_group)
  df_bg    <- df_traj_plot %>% filter(as.character(diabetes_label) != highlight_group)
  df_gm_fg <- df_gm_dm     %>% filter(as.character(diabetes_label) == highlight_group)

  n_fg      <- nrow(df_fg)
  n_dec_fg  <- sum(df_fg$mean_delta_circ < 0, na.rm = TRUE)
  pct_dec_fg <- round(100 * n_dec_fg / n_fg)

  ggplot() +
    # ① KDE contours (grey background density)
    geom_contour_filled(
      data        = df_kde,
      aes(x = x, y = y, z = z),
      breaks      = .kde_brks,
      alpha       = 0.40,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = colorRampPalette(c("grey97", "grey78"))(n_kde_lev),
      guide  = "none"
    ) +
    # ② Background group — grey arrows + dots (draw first, underneath)
    geom_segment(
      data      = df_bg,
      aes(x = x0, y = y0, xend = x1, yend = y1),
      arrow     = arrow(length = unit(0.10, "cm"), type = "closed"),
      colour    = "grey72", linewidth = 0.32, alpha = 0.38, lineend = "round"
    ) +
    geom_point(data = df_bg, aes(x = x0, y = y0),
               shape = 1,  size = 1.1, colour = "grey72", alpha = 0.38) +
    geom_point(data = df_bg, aes(x = x1, y = y1),
               shape = 16, size = 1.5, colour = "grey72", alpha = 0.38) +
    # ③ Highlighted group — arrows coloured by mean Δ circumference
    geom_segment(
      data      = df_fg,
      aes(x = x0, y = y0, xend = x1, yend = y1,
          colour = mean_delta_circ),
      arrow     = arrow(length = unit(0.13, "cm"), type = "closed"),
      linewidth = 0.45, alpha = 0.75, lineend = "round"
    ) +
    scale_colour_gradient2(
      name     = "Mean Δ Circ\n(cm)\nNegative = Decreased",
      low      = "#240046",
      mid      = "#e9ecef",
      high     = "#f35b04",
      midpoint = 0,
      limits   = c(-.circ_clim, .circ_clim),
      guide    = guide_colorbar(
        barheight      = unit(7, "cm"),
        barwidth       = unit(0.85, "cm"),
        title.position = "right",
        title.vjust    = 0.5,
        title.theme    = element_text(angle = 90, hjust = 0.5, size = 9)
      )
    ) +
    # ④ Highlighted group — dots (same gradient colour as the segments)
    geom_point(data = df_fg, aes(x = x0, y = y0, colour = mean_delta_circ),
               shape = 1,  size = 1.5, alpha = 0.65) +
    geom_point(data = df_fg, aes(x = x1, y = y1, colour = mean_delta_circ),
               shape = 16, size = 2.0, alpha = 0.90) +
    # ⑤ Group-mean centroid arrow (highlighted group) — same gradient as segments
    geom_segment(
      data      = df_gm_fg,
      aes(x = x0, y = y0, xend = x1, yend = y1, colour = mean_delta_circ),
      arrow     = arrow(length = unit(0.3, "cm"), type = "closed"),
      linewidth = 2.5
    ) +
    geom_point(data = df_gm_fg, aes(x = x0, y = y0, colour = mean_delta_circ),
               shape = 18, size = 6.5) +
    geom_point(data = df_gm_fg, aes(x = x1, y = y1, colour = mean_delta_circ),
               shape = 5,  size = 6.5, stroke = 1.8) +
    # ⑥ Annotation
    annotate(
      "label",
      x     = -Inf, y = -Inf,
      label = paste0(
        highlight_group, "  n = ", n_fg, "\n",
        n_dec_fg, " (", pct_dec_fg, "%) mean Δ circ < 0 cm"
      ),
      hjust = -0.05, vjust = -0.15,
      size  = 3.2, label.size = 0.35,
      fill  = "white", colour = "grey30"
    ) +
    labs(
      x        = "UMAP Dimension 1",
      y        = "UMAP Dimension 2",
      title    = paste0("Trajectory UMAP — ", highlight_group, " Highlighted"),
      subtitle = paste0(
        "8-Week Intervention  |  ◆ = group mean W0   ◇ = group mean W8\n"
        #"Grey = T2D group"
      )
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle    = element_text(hjust = 0.5, size = 9.5, colour = "grey35"),
      legend.position  = "right",
      legend.title     = element_text(size = 8.5),
      legend.text      = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey93")
    )
}

fig_umap_spot_nont2d <- .make_umap_spotlight("Non-T2D", .dm_dot_colors["Non-T2D"])
fig_umap_spot_t2d    <- .make_umap_spotlight("T2D",     .dm_dot_colors["T2D"])

# Save individually (Fig5b — two files, same panel)
ggsave("Fig5b_umap_spotlight_nont2d.pdf", fig_umap_spot_nont2d,
       width = 9.5, height = 8, bg = "white")
cat("Saved: Fig5b_umap_spotlight_nont2d.pdf\n")

ggsave("Fig5b_umap_spotlight_t2d.pdf", fig_umap_spot_t2d,
       width = 9.5, height = 8, bg = "white")
cat("Saved: Fig5b_umap_spotlight_t2d.pdf\n")


# ── Figure 5C: Lollipop Grid — Median Δcircumference by DM × Sex ─────────────
# 2-row (Male / Female) × 2-col (Non-T2D / T2D) grid

.delta_long6 <- df_3d %>%
  mutate(
    diabetes_label = factor(
      dplyr::recode(as.character(diabetes_label),
                    "non-T2DM" = "Non-T2D", "T2DM" = "T2D",
                    "non-T2D"  = "Non-T2D"),
      levels = c("Non-T2D", "T2D")
    ),
    Gender = factor(Gender, levels = c("Male", "Female"))
  ) %>%
  dplyr::select(diabetes_label, Gender, all_of(delta_cols)) %>%
  tidyr::pivot_longer(
    cols      = all_of(delta_cols),
    names_to  = "delta_var",
    values_to = "delta_cm"
  ) %>%
  mutate(
    raw_var = gsub("^delta_", "", delta_var),
    raw_var = gsub("_", ".", raw_var)
  ) %>%
  filter(!is.na(delta_cm), !is.na(diabetes_label), !is.na(Gender))

lollipop5c_data <- .delta_long6 %>%
  group_by(diabetes_label, Gender, raw_var) %>%
  summarise(
    median_delta = median(delta_cm, na.rm = TRUE),
    q25          = quantile(delta_cm, 0.25, na.rm = TRUE),
    q75          = quantile(delta_cm, 0.75, na.rm = TRUE),
    n            = sum(!is.na(delta_cm)),
    .groups      = "drop"
  ) %>%
  mutate(
    var_label = dplyr::recode(raw_var, !!!.circ_label_map, .default = raw_var),
    var_label = factor(var_label, levels = .anatomical_levels[
      .anatomical_levels %in% dplyr::recode(raw_var, !!!.circ_label_map, .default = raw_var)
    ]),
    # colour key: must match names in .gender_colors_by_dm ("non-T2DM Male" etc.)
    dm_orig    = dplyr::recode(as.character(diabetes_label),
                                "Non-T2D" = "non-T2DM", "T2D" = "T2DM"),
    fill_group = factor(paste(dm_orig, Gender),
                        levels = names(.gender_colors_by_dm))
  )

# Build column-strip labels with per-subgroup n
.strip_n6 <- lollipop5c_data %>%
  group_by(diabetes_label, Gender) %>%
  summarise(n_val = max(n), .groups = "drop")

.dm_strip_labels6 <- lollipop5c_data %>%
  distinct(diabetes_label) %>%
  rowwise() %>%
  mutate(
    n_male   = .strip_n6$n_val[.strip_n6$diabetes_label == diabetes_label &
                                  .strip_n6$Gender == "Male"],
    n_female = .strip_n6$n_val[.strip_n6$diabetes_label == diabetes_label &
                                  .strip_n6$Gender == "Female"],
    dm_strip = paste0(diabetes_label)
  ) %>%
  ungroup()

.dm_recode6 <- setNames(.dm_strip_labels6$dm_strip,
                         as.character(.dm_strip_labels6$diabetes_label))

lollipop5c_data <- lollipop5c_data %>%
  mutate(
    dm_strip = factor(
      dplyr::recode(as.character(diabetes_label), !!!.dm_recode6),
      levels = unname(.dm_recode6)
    )
  )

fig_lollipop5c <- ggplot(lollipop5c_data,
                        aes(x = median_delta, y = var_label, colour = fill_group)) +
  geom_vline(xintercept = 0, linetype = "solid",
             colour = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = q25, xmax = q75),
                 height = 0.35, linewidth = 0.7, alpha = 0.6) +
  geom_segment(aes(x = 0, xend = median_delta, yend = var_label),
               linewidth = 0.7) +
  geom_point(size = 3.5) +
  facet_grid(Gender ~ dm_strip) +
  scale_colour_manual(values = .gender_colors_by_dm, name = NULL) +
  scale_x_continuous(
    name   = "Δ Circumference (W8 − W0, cm)",
    expand = expansion(mult = c(0.20, 0.20))
  ) +
  labs(
    title    = "Median Change in Body Circumferences after 8-Week Intervention",
    subtitle = "Median (IQR)  |  Stratified by Diabetes Status × Sex",
    y        = NULL
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle      = element_text(hjust = 0.5, colour = "grey40", size = 10),
    axis.text.y        = element_text(size = 10, colour = "black"),
    axis.text.x        = element_text(size = 10),
    axis.ticks.y       = element_line(linewidth = 0.3),
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.4),
    strip.text.x       = element_text(face = "bold", size = 11),
    strip.text.y       = element_text(face = "bold", size = 11, angle = 0),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    legend.position    = "top",
    legend.text        = element_text(size = 11),
    panel.spacing.x    = unit(1.2, "lines"),
    panel.spacing.y    = unit(0.8, "lines")
  )

print(fig_lollipop5c)

.lollipop6_h <- max(10, length(levels(lollipop5c_data$var_label)) * 0.45 + 4)
ggsave("Fig5c_median_delta_circ_by_dm_sex.pdf", fig_lollipop5c,
       width = 16, height = .lollipop6_h, bg = "white")
cat("Saved: Fig5c_median_delta_circ_by_dm_sex.pdf\n")




# Clean labels and significance stars
lollipop1_data <- lollipop1_data %>%
  mutate(
    sig_label      = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    ),
    # Map raw variable name → clean display label via explicit lookup
    var_label      = dplyr::recode(variable, !!!.circ_label_map, .default = variable),
    # x position and hjust for significance stars (avoid overlapping with dot)
    sig_x          = r + ifelse(r >= 0,  0.04, -0.04),
    sig_hjust      = ifelse(r >= 0, 0, 1),
    diabetes_label = factor(
      dplyr::recode(as.character(diabetes_label),
                    "non-T2DM" = "Non-T2D", "T2DM" = "T2D"),
      levels = c("Non-T2D", "T2D")
    )
  )

# Apply anatomical ordering (only levels present in data are used)
.present   <- unique(as.character(lollipop1_data$var_label))
.in_order  <- .anatomical_levels[.anatomical_levels %in% .present]
.remainder <- setdiff(.present, .in_order)   # any unmapped vars go to bottom
.var_order <- c(.remainder, .in_order)

lollipop1_data <- lollipop1_data %>%
  mutate(var_label = factor(var_label, levels = .var_order))

# Panel-level sample sizes for strip annotation
.strip_n <- lollipop1_data %>%
  distinct(diabetes_label, n) %>%       # n is per variable; take one typical value
  group_by(diabetes_label) %>%
  summarise(n_min = min(n), n_max = max(n), .groups = "drop") %>%
  mutate(strip_lab = paste0(diabetes_label, "\n(n = ", n_min,
                            ifelse(n_min == n_max, "", paste0("–", n_max)), ")"))
.strip_labels <- setNames(.strip_n$strip_lab, .strip_n$diabetes_label)

lollipop1_data <- lollipop1_data %>%
  mutate(dm_strip = factor(
    dplyr::recode(as.character(diabetes_label), !!!.strip_labels),
    levels = .strip_labels
  ))

fig_lollipop1 <- ggplot(lollipop1_data,
                        aes(x = r, y = var_label, color = diabetes_label)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey60", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = r, yend = var_label),
               linewidth = 0.9) +
  geom_point(size = 3.5) +
  geom_text(aes(x = sig_x, label = sig_label, hjust = sig_hjust),
            size = 3.5, show.legend = FALSE) +
  scale_color_manual(
    name   = NULL,
    values = c("Non-T2D" = "#E69F00", "T2D" = "#CC79A7")
  ) +
  facet_wrap(~ dm_strip, ncol = 2, scales = "free_y") +
  labs(
    title    = "Spearman r: Δ Circumference vs Δ Fasting insulin",
    subtitle = "Δ = W8 − W0  |  * p<0.05  ** p<0.01  *** p<0.001 (BH-adjusted within group)",
    x        = "Spearman r",
    y        = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title         = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle      = element_text(hjust = 0.5, size = 9, color = "grey40"),
    legend.position    = "none",
    strip.text         = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

print(fig_lollipop1)

.lollipop1_h <- max(6, length(circ_base_vars) * 0.42 + 2)
ggsave("Fig5d_corr_circ_insulin.pdf", fig_lollipop1,
       width = 12, height = .lollipop1_h)
cat("Saved: Fig5d_corr_circ_insulin.pdf\n")

# ── Figure 5E: Dumbbell — circumference W0→W8 coloured by Δinsulin ───────────
# Requires .top_vars (from Fig5d / lollipop1_data above) and .circ_label_map,
# .anatomical_levels (defined in Supplementary Table 7&8 section).
# x = insulin_change (continuous, one point per client)
# Vertical segment: circumference W0 → W8 for that client
# Top 6 circumference variables by mean |r| across both DM groups
# Facet: var_label (row) × DM group (column)

# Select top 6 variables
.top_vars <- lollipop1_data %>%
  group_by(variable) %>%
  summarise(mean_abs_r = mean(abs(r)), .groups = "drop") %>%
  slice_max(mean_abs_r, n = 5, with_ties = FALSE) %>%
  pull(variable)

cat("\nTop 5 circumference variables by pooled |r|:\n")
cat(paste(.top_vars, collapse = "\n"), "\n")

# Build long dataset: one row per client × top variable, with W0, W8, Δ
.w0_cols <- paste0(.top_vars, "_w0")
.w8_cols <- paste0(.top_vars, "_w8")

df_dumbbell <- df_3d %>%
  dplyr::select(client_id, insulin_change, diabetes_label,
         all_of(.w0_cols), all_of(.w8_cols)) %>%
  # Pivot W0 and W8 columns together → long
  pivot_longer(
    cols      = c(all_of(.w0_cols), all_of(.w8_cols)),
    names_to  = "col_name",
    values_to = "circ_value"
  ) %>%
  mutate(
    timepoint = if_else(grepl("_w0$", col_name), "W0", "W8"),
    variable  = sub("_(w0|w8)$", "", col_name)
  ) %>%
  dplyr::select(-col_name) %>%
  pivot_wider(names_from = timepoint, values_from = circ_value) %>%
  filter(!is.na(W0), !is.na(W8)) %>%
  mutate(
    circ_change    = W8 - W0,
    direction      = factor(
      ifelse(circ_change <= 0, "Decreased ↓", "Increased ↑"),
      levels = c("Decreased ↓", "Increased ↑")
    ),
    # Map raw variable name → clean display label via explicit lookup
    var_label      = dplyr::recode(variable, !!!.circ_label_map, .default = variable),
    # Preserve anatomical order among the selected top vars
    var_label      = factor(var_label,
                            levels = .anatomical_levels[
                              .anatomical_levels %in%
                                dplyr::recode(.top_vars, !!!.circ_label_map,
                                              .default = .top_vars)]),
    diabetes_label = factor(
      dplyr::recode(as.character(diabetes_label),
                    "non-T2DM" = "Non-T2D", "T2DM" = "T2D",
                    "non-T2D"  = "Non-T2D"),
      levels = c("Non-T2D", "T2D"))
  )

cat("\ndf_dumbbell rows:", nrow(df_dumbbell),
    " (", n_distinct(df_dumbbell$client_id), "clients ×",
    n_distinct(df_dumbbell$variable), "vars)\n")

fig_lollipop3 <- ggplot(df_dumbbell, aes(x = insulin_change)) +
  # Vertical segment: W0 → W8 per client
  geom_segment(
    aes(xend = insulin_change, y = W0, yend = W8, color = direction),
    alpha = 0.50, linewidth = 0.45
  ) +
  # W0 endpoint (open circle, grey)
  geom_point(aes(y = W0),
             shape = 21, fill = "white", color = "grey55",
             size = 1.0, alpha = 0.65) +
  # W8 endpoint (filled, colored by direction)
  geom_point(aes(y = W8, color = direction),
             size = 1.3, alpha = 0.7) +
  # Loess trend: W0 baseline (dashed) and W8 post-intervention (solid)
  geom_smooth(aes(y = W0), method = "loess", formula = y ~ x,
              color = "#B09C85", linewidth = 1.0,
              se = TRUE, alpha = 0.15, linetype = "dashed") +
  geom_smooth(aes(y = W8), method = "loess", formula = y ~ x,
              color = "#D55E00", linewidth = 1.0,
              se = TRUE, alpha = 0.15, linetype = "solid") +
  scale_color_manual(
    name   = "Circumference Change",
    values = c("Decreased ↓" = "#240046", "Increased ↑" = "#f35b04")
  ) +
  facet_grid(var_label ~ diabetes_label, scales = "free_y", switch = "y") +
  labs(
    title    = "3D Circumference W0 → W8 by Δ Fasting insulin",
    subtitle = paste0(
      "Top 5 circumferences by |Spearman r|  |  ",
      "open = W0, filled = W8  |  loess: dashed = W0, solid = W8"
    ),
    x = "Δ Fasting insulin (μU/mL)",
    y = "Circumference (cm)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle     = element_text(hjust = 0.5, size = 8, color = "grey40"),
    strip.text.x      = element_text(face = "bold", size = 10),
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5),
    strip.placement   = "outside",
    legend.position   = "bottom",
    legend.text       = element_text(size = 9),
    panel.grid.minor  = element_blank(),
    panel.spacing.y   = unit(0.6, "lines")
  )

print(fig_lollipop3)

ggsave("Fig5e_dumbbell_circ_insulin.pdf", fig_lollipop3,
       width = 10, height = 14)
cat("Saved: Fig5e_dumbbell_circ_insulin.pdf\n")


# ==============================================================================
# Session Info — for reproducibility (required by Nature journals)
# ==============================================================================
# Prints and saves the full R session information:
#   R version, platform, locale, attached packages (with version numbers)
# Output: session_info.txt (saved to OUT_DIR)
# Include this file when submitting code for peer review.

cat("\n=== R Session Information ===\n")
print(sessionInfo())

writeLines(
  capture.output(sessionInfo()),
  con = file.path(OUT_DIR, "session_info.txt")
)
cat("session_info.txt saved.\n")
