# ============================================================
# Title: Gamma GLM for Microplastic Mean Feret Diameter
# Author: Dr. Carlo Andrea Cossu
# Description:
# Models variation in microplastic particle size (mean Feret
# diameter) using a Gamma GLM with log link.
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(glmmTMB)

# ---- 2. Load data ----
dat <- mps

# ---- 3. Data preprocessing ----
dat_filtered <- dat[
  sampling_site == "Timbavati reserve", sampling_site := "Greater Kruger National Park"
][
  sampling_site == "Kruger National Park", sampling_site := "Greater Kruger National Park"
][
  !is.na(m_feret)
]

# ---- 4. Quick structure check ----
xtabs(~ animal_species + sampling_site + mps_type, data = dat_filtered)

# ---- 5. Exploratory analysis ----
hist(dat_filtered$m_feret,
     breaks = 30,
     main = "Distribution of mean Feret diameter",
     xlab = "Feret diameter")

# ---- 6. Factor conversion ----
dat_filtered[, sample_type := factor(sample_type)]
dat_filtered[, sampling_site := factor(sampling_site)]
dat_filtered[, animal_species := factor(animal_species)]
dat_filtered[, mps_type := factor(mps_type)]

# ---- 7. Reference levels ----
dat_filtered[, sample_type := relevel(sample_type, ref = "EDTA blood")]
dat_filtered[, sampling_site := relevel(sampling_site, ref = "Greater Kruger National Park")]
dat_filtered[, animal_species := relevel(animal_species, ref = "Waterbuck")]
dat_filtered[, mps_type := relevel(mps_type, ref = "Styrene-Isoprene")]

# ---- 8. Model fitting ----
model_feret <- glmmTMB(
  m_feret ~ sample_type + animal_species + mps_type + sampling_site,
  family = Gamma(link = "log"),
  data = dat_filtered
)

summary(model_feret)

# ---- 9. Exponentiated coefficients ----
exp(fixef(model_feret)$cond)

# ---- 10. Diagnostics ----
res <- residuals(model_feret)

qqnorm(res)
qqline(res)

shapiro.test(res)
