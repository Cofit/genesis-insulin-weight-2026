# ==============================================================================
# generate_synthetic_data.R
#
# Generates sample_data.csv: a FULLY SYNTHETIC dataset whose only
# purpose is to let reviewers run analysis_code_submission.R end-to-end.
#
# - No real participant data were read, copied, perturbed or fitted.
# - Every distribution parameter below was set by hand from general clinical
#   plausibility (adults with overweight/obesity in an 8-week weight-management
#   programme); none was estimated from the study dataset.
# - Derived variables (BMI, HOMA-IR, deltas, % changes, energy, protein/kg,
#   clinical status groups) are computed by formula from generated values.
# - No associations between exposures (insulin change, protein intake) and
#   outcomes were built in. Results obtained from this file are meaningless
#   and must not be compared with the manuscript.
# - IDs (SYN_0001 ...) and dates (2030-2031) are artificial.
#
# Requires base R only. Run:  Rscript generate_synthetic_data.R
# ==============================================================================

set.seed(20261118)

N        <- 200
OUT_FILE <- "sample_data.csv"

# ── helpers ───────────────────────────────────────────────────────────────────
clip  <- function(x, lo, hi) pmin(pmax(x, lo), hi)
rnc   <- function(n, mean, sd, lo, hi) clip(rnorm(n, mean, sd), lo, hi)
rln   <- function(n, median, sdlog, lo, hi) clip(median * exp(rnorm(n, 0, sdlog)), lo, hi)
miss  <- function(x, keep) { x[!keep] <- NA; x }

# ── group structure (fixed counts so every T2DM x sex subgroup is populated) ─
n_dm   <- 70
dm     <- sample(c(rep(TRUE, n_dm), rep(FALSE, N - n_dm)))
male   <- logical(N)
male[dm]  <- sample(rep(c(TRUE, FALSE), c(round(n_dm * 0.3), n_dm - round(n_dm * 0.3))))
male[!dm] <- sample(rep(c(TRUE, FALSE), c(round((N - n_dm) * 0.25),
                                          (N - n_dm) - round((N - n_dm) * 0.25))))
m <- as.numeric(male)

# ── demographics & baseline body composition ─────────────────────────────────
age    <- ifelse(dm, rnc(N, 52, 10, 25, 75), rnc(N, 42, 11, 20, 70))
height <- ifelse(male, rnorm(N, 172, 6), rnorm(N, 159, 5.5))           # cm
bmi0   <- rnc(N, 30.5, 4.5, 24, 48) + ifelse(dm, 1, 0)
w0     <- round(bmi0 * (height / 100)^2, 1)
bmi0   <- w0 / (height / 100)^2

bfp0   <- clip(1.2 * bmi0 + 0.23 * age - 10.8 * m - 5.4 + rnorm(N, 0, 3), 15, 60)
fat0   <- round(w0 * bfp0 / 100, 1)
ffm0   <- w0 - fat0
mus0   <- round(ffm0 * rnorm(N, 0.55, 0.02), 1)
waist0 <- round(30 + 2.1 * bmi0 + 4 * m + rnorm(N, 0, 5), 1)

# ── baseline metabolic markers ────────────────────────────────────────────────
fpg0   <- round(ifelse(dm, rnc(N, 145, 30, 100, 300), rnc(N, 95, 9, 75, 125)))
hba0   <- round(ifelse(dm, clip(4.0 + 0.025 * fpg0 + rnorm(N, 0, 0.6), 6.0, 12),
                           clip(4.2 + 0.013 * fpg0 + rnorm(N, 0, 0.25), 4.6, 6.4)), 1)
g2h0   <- round(ifelse(dm, clip(fpg0 * rnorm(N, 1.6, 0.3), 140, 450),
                           clip(fpg0 * rnorm(N, 1.3, 0.25), 70, 199)))
ins0   <- round(rln(N, 12 * (bmi0 / 30)^1.5 * ifelse(dm, 1.2, 1), 0.45, 2, 80), 1)
ins1h0 <- round(ins0 * exp(rnorm(N, log(ifelse(dm, 4.5, 7)), 0.5)), 1)
ins2h0 <- round(ins0 * exp(rnorm(N, log(ifelse(dm, 5.0, 6)), 0.55)), 1)

tg0    <- round(rln(N, ifelse(dm, 150, 120), 0.45, 40, 800))
hdl0   <- round(clip(ifelse(male, rnorm(N, 43, 8), rnorm(N, 52, 10)) - 3 * dm, 25, 100))
ldl0   <- round(rnc(N, 125, 30, 40, 250))
tc0    <- round(ldl0 + hdl0 + tg0 / 5 + rnorm(N, 0, 5))
sbp0   <- round(rnc(N, 126 + 4 * dm, 14, 90, 190))
dbp0   <- round(clip(0.5 * sbp0 + rnorm(N, 17, 7), 55, 120))

# ── 8-week change (no exposure-outcome association is built in) ──────────────
ppk    <- rnc(N, 1.0, 0.25, 0.4, 1.8)                                   # protein g/kg
wlp    <- rnc(N, -7, 2.8, -16, 1)                                       # % weight change
w8     <- round(w0 * (1 + wlp / 100), 1)
dW     <- w8 - w0
dFat   <- dW * rnc(N, 0.75, 0.12, 0.4, 1.0)
fat8   <- round(fat0 + dFat, 1)
dFFM   <- dW - (fat8 - fat0)
mus8   <- round(mus0 + dFFM * rnorm(N, 0.55, 0.05) + rnorm(N, 0, 0.3), 1)
waist8 <- round(waist0 + 0.9 * dW + rnorm(N, 0, 1.5), 1)

fpg8   <- round(clip(fpg0 + ifelse(dm, rnorm(N, -18, 15), rnorm(N, -4, 7)), 70, 300))
hba8   <- round(clip(hba0 + ifelse(dm, rnorm(N, -0.7, 0.5), rnorm(N, -0.15, 0.15)), 4.5, 12), 1)
ins8   <- round(clip(ins0 * exp(rnorm(N, -0.2, 0.35)), 1.5, 80), 1)
g2h8   <- round(clip(g2h0 + rnorm(N, -15, 20), 60, 450))
tg8    <- round(clip(tg0 * exp(rnorm(N, -0.15, 0.25)), 35, 800))
hdl8   <- round(clip(hdl0 + rnorm(N, -1, 4), 25, 100))
ldl8   <- round(clip(ldl0 + rnorm(N, -8, 15), 35, 250))
tc8    <- round(ldl8 + hdl8 + tg8 / 5 + rnorm(N, 0, 5))
sbp8   <- round(clip(sbp0 + rnorm(N, -5, 9), 90, 190))
dbp8   <- round(clip(dbp0 + rnorm(N, -3, 6), 50, 120))

# ── diet during the programme ────────────────────────────────────────────────
kcal_t  <- ifelse(male, rnc(N, 1300, 170, 800, 1900), rnc(N, 1050, 150, 700, 1600))
prot_g  <- round(ppk * w0, 1)
carb_g  <- round(kcal_t * rnc(N, 18, 6, 5, 40) / 100 / 4, 1)
fat_g   <- round(pmax((kcal_t - 4 * prot_g - 4 * carb_g) / 9, 15), 1)
energy  <- round(4 * prot_g + 9 * fat_g + 4 * carb_g, 1)

# ── missingness by design (MCAR) ─────────────────────────────────────────────
bp_w0_ok   <- runif(N) > 0.15
bp_w8_ok   <- runif(N) > 0.15
ogtt_w0_ok <- runif(N) > 0.10
g2h_w8_ok  <- ogtt_w0_ok & runif(N) < 0.20
scan_w0_ok <- runif(N) < 0.90                 # 3D body scan done at week 0
scan_w8_ok <- scan_w0_ok & runif(N) < 0.70    # ... and at week 8

sbp0 <- miss(sbp0, bp_w0_ok); dbp0 <- miss(dbp0, bp_w0_ok)
sbp8 <- miss(sbp8, bp_w8_ok); dbp8 <- miss(dbp8, bp_w8_ok)
g2h0   <- miss(g2h0,   ogtt_w0_ok)
ins1h0 <- miss(ins1h0, ogtt_w0_ok)
ins2h0 <- miss(ins2h0, ogtt_w0_ok)
g2h8   <- miss(g2h8,   g2h_w8_ok)

# ── clinical status groups (definitions as in column_dictionary.pdf) ────────
dm_group <- function(fpg, a1c) ifelse(a1c >= 6.5 | fpg >= 126, "DM criteria",
                                ifelse(a1c >= 5.7 | fpg >= 100, "Pre-DM criteria",
                                       "Normal glycemia"))
met_group <- function(waist, tg, hdl, sbp, dbp, fpg) {
  k <- (waist >= ifelse(male, 90, 80)) + (tg >= 150) +
       (hdl < ifelse(male, 40, 50)) +
       (!is.na(sbp) & (sbp >= 130 | dbp >= 85)) + (fpg >= 100)
  ifelse(k >= 3, "Metabolic syndrome", "Normal")
}
hlp_group <- function(tc, ldl, tg, hdl) ifelse(tc >= 200 | ldl >= 130 | tg >= 150 |
                                               hdl < ifelse(male, 40, 50),
                                               "Hyperlipidemia", "Normal")

# ── 3D body scan anthropometrics (cm) ────────────────────────────────────────
r1 <- function(x) round(x, 1)
scan0 <- data.frame(
  Body.Height                = height + rnorm(N, 0, 0.8),
  Rise.Height                = 0.22 * height + rnorm(N, 0, 2),
  Scye.Depth.A               = 0.12 * height + rnorm(N, 0, 1.2),
  Halter.Length              = 0.47 * height + 0.3 * (bmi0 - 30) + rnorm(N, 0, 2.5),
  Vertical.Trunk             = 60 + 0.45 * height + 1.5 * bmi0 + rnorm(N, 0, 6),
  Chest.Circumference        = 50 + 1.6 * bmi0 + 6 * m + rnorm(N, 0, 4),
  Narrow.Waist.Circumference = waist0 - 12 + rnorm(N, 0, 3),
  Low.Hip.Circumference      = 45 + 2.0 * bmi0 - 4 * m + rnorm(N, 0, 4)
)
scan0$High.Hip.Circumference <- scan0$Low.Hip.Circumference - 3 + rnorm(N, 0, 3)
scan0$Belly.Circumference    <- waist0 - 3 + rnorm(N, 0, 4)
arm    <- 0.33 * height; inseam <- 0.44 * height
thigh  <- 20 + 1.3 * bmi0 - 3 * m;  knee  <- 22 + 0.55 * bmi0 + 2 * m
calf   <- 18 + 0.65 * bmi0 + 2 * m; ankle <- 14 + 0.25 * bmi0 + 1.5 * m
uarm   <- 8 + 0.8 * bmi0 + m;       wrist <- 10 + 0.15 * bmi0 + 1.5 * m
for (side in c("Left", "Right")) scan0[[paste0(side, ".Arm.Length")]]          <- arm    + rnorm(N, 0, 1.2)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Inseam.Length")]]       <- inseam + rnorm(N, 0, 1.5)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Thigh.Circumference")]] <- thigh  + rnorm(N, 0, 3)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Knee.Circumference")]]  <- knee  + rnorm(N, 0, 1.5)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Calf.Circumference")]]  <- calf  + rnorm(N, 0, 1.8)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Ankle.Circumference")]] <- ankle + rnorm(N, 0, 1.5)
scan0$Neck.Base.Circumference <- 22 + 0.55 * bmi0 + 5 * m + rnorm(N, 0, 1.8)
scan0$Collar.Circumference    <- scan0$Neck.Base.Circumference - 2 + rnorm(N, 0, 1)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Upper.Arm.Circumference")]] <- uarm  + rnorm(N, 0, 2)
for (side in c("Left", "Right")) scan0[[paste0(side, ".Wrist.Circumference")]]     <- wrist + rnorm(N, 0, 1)
scan0$Total.Rise.Length <- 40 + 1.6 * bmi0 + rnorm(N, 0, 4)
scan0$Pants.Rise.Length <- 0.6 * scan0$Total.Rise.Length + 12 + rnorm(N, 0, 3)

# Week-8 scan: circumferences shrink in proportion to weight change
# (k = sensitivity to % weight change); lengths only get measurement noise.
k_scale <- c(Chest = 0.5, Narrow.Waist = 0.6, Low.Hip = 0.45, High.Hip = 0.45,
             Belly = 0.6, Thigh = 0.4, Knee = 0.2, Calf = 0.3, Ankle = 0.15,
             Neck.Base = 0.2, Collar = 0.2, Upper.Arm = 0.4, Wrist = 0.1,
             Vertical.Trunk = 0.3, Total.Rise = 0.4)
scan8 <- scan0
for (v in names(scan0)) {
  key <- sub("^(Left|Right)\\.", "", sub("\\.(Circumference|Length)$", "", v))
  if (key %in% names(k_scale)) {
    scan8[[v]] <- scan0[[v]] * (1 + k_scale[[key]] * wlp / 100) + rnorm(N, 0, 0.8)
  } else {
    scan8[[v]] <- scan0[[v]] + rnorm(N, 0, 0.6)
  }
}
scan0[] <- lapply(scan0, function(x) r1(miss(x, scan_w0_ok)))
scan8[] <- lapply(scan8, function(x) r1(miss(x, scan_w8_ok)))
names(scan0) <- paste0(names(scan0), "_w0")
names(scan8) <- paste0(names(scan8), "_w8")

d0     <- as.Date("2030-01-01") + sample(0:729, N, replace = TRUE)
d8     <- d0 + round(rnc(N, 56, 3, 49, 63))
date0  <- ifelse(scan_w0_ok, format(d0), NA)
date8  <- ifelse(scan_w8_ok, format(d8), NA)

# ── assemble (column order matches the study dataset, minus birthday) ────────
bmi8 <- w8 / (height / 100)^2
homa0 <- fpg0 * ins0 / 405
homa8 <- fpg8 * ins8 / 405
r2 <- function(x) round(x, 2)
r3 <- function(x) round(x, 3)

df <- data.frame(
  client_id                              = sprintf("SYN_%04d", seq_len(N)),
  gender                                 = ifelse(male, "male", "female"),
  age                                    = r1(age),
  weight_w0                              = w0,
  bmi_w0                                 = r2(bmi0),
  body_fat_mass_w0                       = fat0,
  body_fat_mass_percentage_w0            = r2(fat0 / w0 * 100),
  weight_without_fat_w0                  = r1(ffm0),
  muscle_mass_w0                         = mus0,
  waist_circumference_w0                 = waist0,
  weight_w8                              = w8,
  bmi_w8                                 = r2(bmi8),
  body_fat_mass_w8                       = fat8,
  body_fat_mass_percentage_w8            = r2(fat8 / w8 * 100),
  muscle_mass_w8                         = mus8,
  waist_circumference_w8                 = waist8,
  glucose_ac_w0                          = fpg0,
  glucose_pc_2hr_w0                      = g2h0,
  insulin_w0                             = ins0,
  insulin_pc_1hr_w0                      = ins1h0,
  insulin_pc_2hr_w0                      = ins2h0,
  hba1c_w0                               = hba0,
  homa_ir_w0                             = r2(homa0),
  triglyceride_w0                        = tg0,
  total_cholesterol_w0                   = tc0,
  high_density_lipoprotein_w0            = hdl0,
  low_density_lipoprotein_cholesterol_w0 = ldl0,
  diastolic_blood_pressure_w0            = dbp0,
  systolic_blood_pressure_w0             = sbp0,
  glucose_ac_w8                          = fpg8,
  glucose_pc_2hr_w8                      = g2h8,
  insulin_w8                             = ins8,
  hba1c_w8                               = hba8,
  homa_ir_w8                             = r2(homa8),
  triglyceride_w8                        = tg8,
  total_cholesterol_w8                   = tc8,
  high_density_lipoprotein_w8            = hdl8,
  low_density_lipoprotein_cholesterol_w8 = ldl8,
  diastolic_blood_pressure_w8            = dbp8,
  systolic_blood_pressure_w8             = sbp8,
  diet_w8_calories                       = energy,
  diet_w8_protein                        = prot_g,
  diet_w8_fat                            = fat_g,
  diet_w8_carbohydrate                   = carb_g,
  weight_loss_pct                        = r3((w8 - w0) / w0 * 100),
  diabetes                               = ifelse(dm, "Yes", "No"),
  insulin_change                         = r1(ins8 - ins0),
  insulin_change_pct                     = r2((ins8 - ins0) / ins0 * 100),
  total_energy                           = energy,
  diabetes_label                         = ifelse(dm, "T2DM", "non-T2DM"),
  muscle_change_pct                      = r3((mus8 - mus0) / mus0 * 100),
  delta_muscle_mass                      = r1(mus8 - mus0),
  delta_weight                           = r1(w8 - w0),
  delta_bmi                              = r2(bmi8 - bmi0),
  delta_waist_circumference              = r1(waist8 - waist0),
  delta_body_fat_mass                    = r1(fat8 - fat0),
  delta_body_fat_mass_percentage         = r2(fat8 / w8 * 100 - fat0 / w0 * 100),
  delta_weight_without_fat               = r1((w8 - fat8) - ffm0),
  delta_hba1c                            = r1(hba8 - hba0),
  delta_glucose_ac                       = fpg8 - fpg0,
  delta_insulin                          = r1(ins8 - ins0),
  delta_homa_ir                          = r2(homa8 - homa0),
  delta_triglyceride                     = tg8 - tg0,
  delta_total_cholesterol                = tc8 - tc0,
  delta_high_density_lipoprotein         = hdl8 - hdl0,
  delta_low_density_lipoprotein_cholesterol = ldl8 - ldl0,
  delta_systolic_blood_pressure          = sbp8 - sbp0,
  delta_diastolic_blood_pressure         = dbp8 - dbp0,
  muscle_mass_percentage_w0              = r2(mus0 / w0 * 100),
  muscle_mass_percentage_w8              = r2(mus8 / w8 * 100),
  protein_g                              = prot_g,
  fat_g                                  = fat_g,
  carb_g                                 = carb_g,
  protein_per_kg                         = r3(prot_g / w0),
  date_w0                                = date0,
  check.names = FALSE, stringsAsFactors = FALSE
)
df <- cbind(df, scan0, date_w8 = date8, scan8)
df$dm_gp     <- dm_group(fpg0, hba0)
df$met_gp    <- met_group(waist0, tg0, hdl0, sbp0, dbp0, fpg0)
df$hlp_gp    <- hlp_group(tc0, ldl0, tg0, hdl0)
df$dm_gp_w8  <- dm_group(fpg8, hba8)
df$met_gp_w8 <- met_group(waist8, tg8, hdl8, sbp8, dbp8, fpg8)
df$hlp_gp_w8 <- hlp_group(tc8, ldl8, tg8, hdl8)

# Missing values are written as "NA" so read.csv() reads them as NA in
# character columns (e.g. date_w0 / date_w8) as well as numeric ones.
write.csv(df, OUT_FILE, row.names = FALSE, na = "NA")
cat("Wrote", OUT_FILE, ":", nrow(df), "rows x", ncol(df), "columns\n")
