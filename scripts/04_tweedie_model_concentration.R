# ============================================================
# Title: Tweedie Model for MP Concentration (including zeros)
# Author: Dr. Carlo Andrea Cossu
# Description:
# Fits a Tweedie GLM using glmmTMB to model microplastic
# concentration including zero values. Produces diagnostics
# and publication-ready tables.
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(glmmTMB)
library(broom.mixed)
library(dplyr)

# ---- 2. Load data ----
dat <- unique(mps, by = c("sample_id"))

# ---- 3. Data preprocessing ----
dat_filtered <- dat[
  sampling_site == "Timbavati reserve", sampling_site := "Greater Kruger National Park"
][
  sampling_site == "Kruger National Park", sampling_site := "Greater Kruger National Park"
][
  sampling_site != "Mountain zebra national park"
][
  sample_type != "Serum"
][
  !animal_species %in% c("Plains zebra", "Gemsbok")
]

# ---- 4. Factor conversion ----
dat_filtered[, sample_type := factor(sample_type)]
dat_filtered[, sampling_site := factor(sampling_site)]
dat_filtered[, animal_species := factor(animal_species)]

# ---- 5. Set reference levels ----
dat_filtered[, sample_type := relevel(sample_type, ref = "EDTA blood")]
dat_filtered[, sampling_site := relevel(sampling_site, ref = "Karoo National Park")]
dat_filtered[, animal_species := relevel(animal_species, ref = "Eland")]

# ---- 6. Exploratory visualization ----
hist(dat_filtered$concentration,
     breaks = 30,
     main = "Distribution of MP concentration (including zeros)",
     xlab = "Concentration")

# ---- 7. Model fitting ----
model_tweedie <- glmmTMB(
  concentration ~ sample_type + sampling_site + animal_species,
  family = tweedie(link = "log"),
  data = dat_filtered
)

summary(model_tweedie)

# ---- 8. Model diagnostics ----
res <- residuals(model_tweedie)

# Normality (approximate, interpret cautiously for GLMs)
shapiro.test(res)

# QQ plot
qqnorm(res)
qqline(res)

# Residuals vs fitted
plot(fitted(model_tweedie), res,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0)

# ---- 9. Exponentiated coefficients ----
exp(fixef(model_tweedie)$cond)

# ---- 10. Create publication-ready table ----
model_table <- tidy(
  model_tweedie,
  effects = "fixed",
  component = "cond",
  conf.int = TRUE
) %>%
  filter(!is.na(estimate),
         !is.na(conf.low),
         !is.na(conf.high)) %>%
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

# ---- 11. Clean term labels ----
model_table <- model_table %>%
  mutate(term = recode(term,
    "(Intercept)" = "Intercept (EDTA blood + Karoo NP + Eland)",
    "sample_typeKidney" = "Kidney",
    "sampling_siteAddo Elephant National Park" = "Addo Elephant National Park",
    "sampling_siteGreater Kruger National Park" = "Greater Kruger National Park",
    "sampling_siteMokala National Park" = "Mokala National Park",
    "animal_speciesAfrican buffalo" = "African buffalo",
    "animal_speciesAfrican elephant" = "African elephant",
    "animal_speciesHartebeest" = "Hartebeest",
    "animal_speciesHippo" = "Hippo",
    "animal_speciesSpotted hyaena" = "Spotted hyaena",
    .default = term
  ))

# ---- 12. Save outputs ----
write.csv(
  model_table,
  file = "outputs/tables/tweedie_model_concentration.csv",
  row.names = FALSE
)

