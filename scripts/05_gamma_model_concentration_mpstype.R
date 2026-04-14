# ============================================================
# Title: Gamma GLM for MP Concentration by Polymer Type
# Author: Dr. Carlo Andrea Cossu
# Description:
# Fits a Gamma GLM using glmmTMB to assess how sampling site,
# sample type, animal species, and MP polymer type influence
# microplastic concentration at the polymer-type level.
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(glmmTMB)
library(broom)
library(dplyr)

# ---- 2. Load data ----
dat <- mps_by_type

# ---- 3. Data preprocessing ----
dat_filtered <- dat[
  sampling_site == "Timbavati reserve", sampling_site := "Greater Kruger National Park"
][
  sampling_site == "Kruger National Park", sampling_site := "Greater Kruger National Park"
][
  concentration_by_mps_type != 0
]

# ---- 4. Convert to factors ----
dat_filtered[, sample_type := factor(sample_type)]
dat_filtered[, sampling_site := factor(sampling_site)]
dat_filtered[, animal_species := factor(animal_species)]
dat_filtered[, mps_type := factor(mps_type)]

# ---- 5. Inspect levels ----
levels(dat_filtered$sample_type)
levels(dat_filtered$sampling_site)
levels(dat_filtered$animal_species)
levels(dat_filtered$mps_type)

# ---- 6. Set reference levels ----
dat_filtered[, sample_type := relevel(sample_type, ref = "EDTA blood")]
dat_filtered[, sampling_site := relevel(sampling_site, ref = "Karoo National Park")]
dat_filtered[, animal_species := relevel(animal_species, ref = "Eland")]
dat_filtered[, mps_type := relevel(mps_type, ref = "Polyester")]

# ---- 7. Exploratory analysis ----
hist(dat_filtered$concentration_by_mps_type,
     breaks = 30,
     main = "MP concentration by polymer type",
     xlab = "Concentration")

# ---- 8. Model fitting ----
model_gamma_mpstype <- glmmTMB(
  concentration_by_mps_type ~ sampling_site + sample_type + mps_type + animal_species,
  family = Gamma(link = "log"),
  data = dat_filtered
)

summary(model_gamma_mpstype)

# ---- 9. Extract fixed effects ----
exp(fixef(model_gamma_mpstype)$cond)

# ---- 10. Diagnostics ----
res <- residuals(model_gamma_mpstype)

# QQ plot
qqnorm(res)
qqline(res)

# Normality test (interpret cautiously for GLMs)
shapiro.test(res)

# Residuals vs fitted
plot(fitted(model_gamma_mpstype), res,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0)

par(mfrow = c(2,2))

# ---- 11. Build publication-ready table ----
model_table <- tidy(model_gamma_mpstype, conf.int = TRUE) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    exp_estimate = exp(estimate),
    exp_conf_low = exp(conf.low),
    exp_conf_high = exp(conf.high),

    `log estimate (CI)` = paste0(
      round(estimate, 2), " [",
      round(conf.low, 2), "–",
      round(conf.high, 2), "]"
    ),

    `exp estimate (CI)` = paste0(
      round(exp_estimate, 2), " [",
      round(exp_conf_low, 2), "–",
      round(exp_conf_high, 2), "]"
    ),

    p.value = signif(p.value, 2)
  ) %>%
  select(term, `log estimate (CI)`, `exp estimate (CI)`, p.value)

# ---- 12. Clean term labels ----
model_table <- model_table %>%
  mutate(term = recode(term,
    "(Intercept)" = "Intercept (EDTA blood + Karoo NP + Eland + Polyester)",

    # Sample type
    "sample_typeKidney" = "Kidney",

    # Sampling sites
    "sampling_siteAddo Elephant National Park" = "Addo Elephant NP",
    "sampling_siteGreater Kruger National Park" = "Greater Kruger NP",
    "sampling_siteMokala National Park" = "Mokala NP",

    # Animal species
    "animal_speciesAfrican buffalo" = "African buffalo",
    "animal_speciesAfrican elephant" = "African elephant",
    "animal_speciesHartebeest" = "Hartebeest",
    "animal_speciesHippo" = "Hippo",
    "animal_speciesSpotted hyaena" = "Spotted hyaena",
    "animal_speciesWarthog" = "Warthog",
    "animal_speciesWaterbuck" = "Waterbuck",

    # Polymer types
    "mps_typePolyamide" = "Polyamide",
    "mps_typePolypropylene" = "Polypropylene",
    "mps_typePolyethylene" = "Polyethylene",
    "mps_typePolyurethane" = "Polyurethane",
    "mps_typeStyrene-Isoprene" = "Styrene-Isoprene",
    "mps_typePolyethylene-co-Polypropylene" = "PE-PP copolymer",
    "mps_typePolyvinyl Chloride" = "Polyvinyl Chloride",
    "mps_typePolymethyl Methacrylate" = "Polymethyl Methacrylate",
    "mps_typePolystyrene" = "Polystyrene",
    .default = term
  ))

# ---- 13. Save outputs ----
write.csv(
  model_table,
  file = "outputs/tables/gamma_concentration_mpstype.csv",
  row.names = FALSE
)
