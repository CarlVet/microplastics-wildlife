# ============================================================
# Title: Gamma GLM for MP Concentration (Sample-level, no zeros)
# Author: [Your Name]
# Description: Fits Gamma model using glmmTMB and produces
#              diagnostic plots and publication-ready tables
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(glmmTMB)
library(broom)
library(dplyr)

# ---- 2. Load data ----
dat <- unique(mps, by = c("sample_id"))

# ---- 3. Data preprocessing ----
dat_filtered <- dat[
  sampling_site == "Timbavati reserve", sampling_site := "Greater Kruger National Park"
][
  sampling_site == "Kruger National Park", sampling_site := "Greater Kruger National Park"
][
  concentration != 0
]

# Convert to factors
dat_filtered[, sample_type := factor(sample_type)]
dat_filtered[, sampling_site := factor(sampling_site)]
dat_filtered[, animal_species := factor(animal_species)]

# Set reference levels
dat_filtered[, sample_type := relevel(sample_type, ref = "EDTA blood")]
dat_filtered[, sampling_site := relevel(sampling_site, ref = "Karoo National Park")]
dat_filtered[, animal_species := relevel(animal_species, ref = "Eland")]

# ---- 4. Exploratory visualization ----
hist(dat_filtered$concentration, breaks = 30,
     main = "Distribution of MP concentration",
     xlab = "Concentration")

# ---- 5. Model fitting ----
model_gamma <- glmmTMB(
  concentration ~ sample_type + sampling_site + animal_species,
  family = Gamma(link = "log"),
  data = dat_filtered
)

summary(model_gamma)

# ---- 6. Model diagnostics ----
res <- residuals(model_gamma)

# QQ plot
qqnorm(res)
qqline(res)

# Normality test (optional)
shapiro.test(res)

# Residuals vs fitted
plot(fitted(model_gamma), res,
     xlab = "Fitted values",
     ylab = "Residuals")
abline(h = 0)

# ---- 7. Extract and format results ----
model_table <- tidy(model_gamma, conf.int = TRUE) %>%
  filter(!is.na(estimate)) %>%
  mutate(
    exp_estimate = exp(estimate),
    exp_conf_low = exp(conf.low),
    exp_conf_high = exp(conf.high),
    `log estimate (CI)` = paste0(round(estimate, 2), " [",
                                 round(conf.low, 2), "-",
                                 round(conf.high, 2), "]"),
    `exp estimate (CI)` = paste0(round(exp_estimate, 2), " [",
                                 round(exp_conf_low, 2), "-",
                                 round(exp_conf_high, 2), "]"),
    p.value = signif(p.value, 2)
  ) %>%
  select(term, `log estimate (CI)`, `exp estimate (CI)`, p.value)

# ---- 8. Clean term labels ----
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
    "animal_speciesSpotted hyaena" = "Spotted hyaena"
  ))

# ---- 9. Save outputs ----
write.csv(model_table, "outputs/gamma_model_concentration_table.csv", row.names = FALSE)

