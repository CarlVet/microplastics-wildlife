# ============================================================
# Title: Exploratory Data Analysis & Visualisation
# Author: Dr. Carlo Andrea Cossu
# Description:
# Generates descriptive statistics and visualisations of
# microplastic distribution across species, sampling sites,
# polymer types, and particle size. Includes barplots,
# boxplots, and interactive spatial mapping.
# ============================================================

# ---- 1. Load packages ----
library(data.table)
library(ggplot2)
library(ggtext)
library(dplyr)
library(leaflet)
library(sf)
library(htmltools)
library(glue)

# ============================================================
# SECTION A — MICROPLASTIC COUNTS BARPLOT
# ============================================================

# ---- Data aggregation ----
dat <- mps[, .(Count = .N),
           by = c("animal_species_label", "mps_type", "sample_type")]

dat <- dat[!mps_type == "-"]

# ---- Colour palette ----
mps_colours <- c(
  "Polyamide" = "#D7263D",
  "Polyester" = "pink",
  "Polyethylene" = "#0096FF",
  "Polyethylene-co-Polypropylene" = "#E2C044",
  "Polymethyl Methacrylate" = "#bf8bff",
  "Polypropylene" = "#4CBB17",
  "Polystyrene" = "blue",
  "Polyurethane" = "navy",
  "Polyvinyl Chloride" = "violet",
  "Styrene-Isoprene" = "black"
)

# ---- Barplot: MP counts ----
ggplot(dat,
       aes(x = animal_species_label,
           y = Count,
           fill = mps_type)) +

  geom_bar(stat = "identity") +

  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            size = 4,
            color = "white") +

  facet_grid(~sample_type, scales = "free_x") +

  scale_fill_manual(values = mps_colours) +

  theme_linedraw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title = element_blank()
  ) +

  ylab("MP count") +
  xlab("")

# ============================================================
# SECTION B — CONCENTRATION BY SAMPLING SITE, SAMPLE TYPE AND SPECIES (BOXPLOT)
# ============================================================

dat <- unique(mps, by = "sample_id")[animal_species == "African buffalo", animal_species := "Buffalo"
           ][animal_species == "African elephant", animal_species := "Elephant"
           ][animal_species == "Spotted hyaena", animal_species := "Hyaena"]

# ---- Species labels per site ----
species_labels <- dat[, .(
  species_list = paste(sort(unique(animal_species)), collapse = ", ")
), by = sampling_site]

dat <- merge(dat, species_labels, by = "sampling_site")

dat[, sampling_label := paste0(
  sampling_site,
  "<br><span style='font-size:8pt;'>(",
  species_list,
  ")</span>"
)]

# ---- Mean ± SD function ----
mean_sd_label <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  data.frame(y = m,
             label = paste0(round(m, 2), " (±", round(s, 2), ")"))
}

# ---- Boxplot ----
ggplot(dat,
       aes(x = sampling_label,
           y = concentration,
           colour = sample_type)) +

  geom_boxplot() +

  stat_summary(fun = mean,
               geom = "line",
               aes(group = 1),
               colour = "black",
               linewidth = 0.5) +

  stat_summary(fun = mean,
               geom = "point",
               colour = "black",
               size = 2) +

  stat_summary(fun.data = mean_sd_label,
               geom = "text",
               aes(label = after_stat(label)),
               vjust = -0.5,
               size = 3) +

  theme_linedraw() +
  theme(
    legend.position = "top",
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1)
  ) +

  ylab("Concentration (MPs/g)") +
  xlab("")

# ============================================================
# SECTION C — FERET DIAMETER BY ANIMAL SPECIES, SAMPLING SITE AND SAMPLE TYPE
# ============================================================

dat <- mps[mps_type != "-"]

dat <- dat[sampling_site == unique(dat$sampling_site)[5], sampling_site := "AENP"
][sampling_site == "Kruger National Park", sampling_site := "KNP"
][sampling_site == "Mokala National Park", sampling_site := "MokNP"
][sampling_site == "Karoo National Park", sampling_site := "KaNP"
][sampling_site == "Timbavati reserve", sampling_site := "TPNR"]

# ---- Site labels per species ----
site_labels <- dat[, .(
  site_list = paste(sort(unique(sampling_site)), collapse = ", ")
), by = animal_species]

dat <- merge(dat, site_labels, by = "animal_species")

dat[, species_label := paste0(
  animal_species,
  "<br><span style='font-size:8pt;'>(",
  site_list,
  ")</span>"
)]

# ---- Mean ± SD ----
mean_sd_label <- function(x) {
  m <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  data.frame(y = m,
             label = paste0(round(m, 2), " (±", round(s, 2), ")"))
}

# ---- Boxplot Feret ----
ggplot(dat,
       aes(x = species_label,
           y = m_feret,
           colour = sample_type)) +

  geom_boxplot() +

  stat_summary(fun = mean, geom = "line", aes(group = 1),
               colour = "black", linewidth = 0.5) +

  stat_summary(fun = mean, geom = "point",
               colour = "black", size = 2) +

  stat_summary(fun.data = mean_sd_label,
               geom = "text",
               aes(label = after_stat(label)),
               vjust = -0.5,
               size = 3) +

  scale_colour_manual(values = c(
    "EDTA blood" = "#F8766D",
    "Kidney" = "#00BA38"
  )) +

  theme_linedraw() +
  theme(
    legend.position = "top",
    axis.text.x = ggtext::element_markdown(angle = 90, vjust = 0.5, hjust = 1)
  ) +

  ylab("Mean Feret diameter") +
  xlab("")

# ============================================================
# SECTION D — INTERACTIVE MAP (LEAFLET)
# ============================================================

dat <- mps
dat <- dat[!mps_type == "-"]

# ---- Polymer factor order ----
dat$mps_type <- factor(dat$mps_type,
                       levels = names(mps_colours))

# ---- Mean concentration per park ----
park_concentration <- dat %>%
  group_by(sampling_site, sample_id) %>%
  summarise(
    mps = n(),
    weight = unique(weight),
    .groups = "drop"
  ) %>%
  mutate(conc = mps / weight) %>%
  group_by(sampling_site) %>%
  summarise(mean_concentration = mean(conc, na.rm = TRUE))

dat <- left_join(dat, park_concentration, by = "sampling_site")

# ---- Leaflet map ----
leaflet(options = leafletOptions(zoomControl = FALSE)) %>%
  addProviderTiles(providers$CartoDB.Voyager) %>%

  addPolygons(data = Karoo_boundaries, weight = 2, color = "grey") %>%
  addPolygons(data = Mokala_boundaries, weight = 2, color = "grey") %>%
  addPolygons(data = Addo_boundaries, weight = 2, color = "grey") %>%
  addPolygons(data = MZNP_boundaries, weight = 2, color = "grey") %>%
  addPolygons(data = Kruger_boundaries, weight = 2, color = "grey") %>%

  addCircleMarkers(
    data = dat,
    lng = ~longitude,
    lat = ~latitude,
    radius = 4,
    fillOpacity = 0.8,
    popup = ~paste0(
      "<b>Polymer:</b> ", mps_type, "<br>",
      "<b>Mean concentration:</b> ", round(mean_concentration, 2)
    )
  ) %>%

  addScaleBar(position = "bottomleft") %>%

  addLegend(
    position = "bottomright",
    colors = mps_colours,
    labels = names(mps_colours),
    title = "Polymer Type"
  )
