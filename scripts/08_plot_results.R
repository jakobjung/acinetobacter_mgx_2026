#!/usr/bin/env Rscript

# =============================================================================
# Plotting script for mSWEEP v3 results (species-filtered)
# Xiao et al. CRAB VAP dataset (PRJNA681291)
# =============================================================================

library(tidyverse)
library(patchwork)

# --- Load data ---------------------------------------------------------------

base_dir <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "..", "analysis")
# If not running in RStudio, set manually:
# base_dir <- "analysis"

abundances <- read_tsv(file.path(base_dir, "all_abundances_v3.tsv"),
                       col_names = c("SRR", "SC", "abundance"),
                       show_col_types = FALSE)

metadata <- read_tsv(file.path(base_dir, "sample_metadata.tsv"),
                     show_col_types = FALSE)

sample_info <- read_tsv(file.path(base_dir, "sample_info_v3.tsv"),
                        show_col_types = FALSE)

qc_meta <- read_tsv(file.path(base_dir, "xiao_qc_metadata.tsv"),
                    show_col_types = FALSE)

sc_label_map <- read_tsv(file.path(base_dir, "sc_labels.tsv"),
                         show_col_types = FALSE)

culture <- read_tsv(file.path(base_dir, "culture_st.tsv"),
                    show_col_types = FALSE)

# --- Merge metadata ----------------------------------------------------------

abun <- abundances %>%
  left_join(metadata, by = "SRR") %>%
  left_join(sample_info %>% select(SRR, aligned, total), by = "SRR")

# --- Define SC categories and colors ----------------------------------------

sc_labels <- setNames(sc_label_map$short_label, sc_label_map$SC)
major_scs <- sc_label_map$SC

# Color palette — focus on Xiao Oxford STs
sc_colors <- c(
  "Oxf ST208" = "#EE6677",
  "Oxf ST938" = "#228833",
  "Oxf ST369" = "#CCBB44",
  "Oxf ST195" = "#66CCEE",
  "Oxf ST1791" = "#AA3377",
  "Oxf ST547" = "#BBBBBB",
  "Oxf ST191" = "#EE8866",
  "Oxf ST2143" = "#44BB99",
  "Oxf ST540" = "#FFAABB",
  "Oxf ST136" = "#99DDFF",
  "Oxf ST381" = "#DDCC77",
  "Oxf ST436" = "#882255",
  "Oxf ST1333 (non-GC2)" = "#CC3311",
  "Oxf ST940" = "#4477AA",
  "Oxf ST940-2" = "#6699CC",
  "Pas ST151" = "#AAAAAA",
  "Pas ST2 (GC2)" = "#004488",
  "Oxf ST893" = "#997700",
  "Pas ST738" = "#559988",
  "Unknown" = "#CCCCCC",
  "Other (<1%)" = "#EEEEEE"
)

# --- Filter low-read samples -------------------------------------------------

MIN_READS <- 500  # minimum A. baumannii reads for reliable estimates

low_read_samples <- sample_info %>%
  filter(as.numeric(total) < MIN_READS) %>%
  pull(SRR)

cat(sprintf("Filtering %d samples with < %d A. baumannii reads\n",
            length(low_read_samples), MIN_READS))

abun_filtered <- abun %>%
  filter(!SRR %in% low_read_samples)

# --- Prepare abundance data for plotting -------------------------------------

abun_plot <- abun_filtered %>%
  mutate(SC_label = ifelse(SC %in% major_scs & abundance > 0.01,
                           sc_labels[SC],
                           "Other (<1%)")) %>%
  group_by(SRR, sample, group, SC_label) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  mutate(SC_label = factor(SC_label, levels = rev(names(sc_colors))))

# Order samples: CRAB-N first, then by number of SCs (most diverse last)
n_scs_per_sample <- abun_plot %>%
  filter(SC_label != "Other (<1%)" & abundance > 0.05) %>%
  group_by(sample, group) %>%
  summarise(n_scs = n(), .groups = "drop")

sample_order <- n_scs_per_sample %>%
  arrange(group, n_scs) %>%
  pull(sample)

# Add samples with 0 SCs (CRAB-N)
all_samples <- metadata$sample
missing <- setdiff(all_samples, sample_order)
sample_order <- c(missing, sample_order)

abun_plot <- abun_plot %>%
  mutate(sample = factor(sample, levels = sample_order))

# =============================================================================
# Plot 1: Stacked bar chart of SC abundances (v3 — species-filtered)
# =============================================================================
library(viridis)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "pink",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "khaki2",
   "#FB9A99", # lt pink
  "skyblue","gold1",
  "#CAB2D6", # lt purple

  "gray70","#FDBF6F", "steelblue")

# Calculate y-position for culture "X" marker (midpoint of the culture-reported SC bar)
# ggplot stacks bars in reverse factor level order, so sort accordingly
culture_markers <- abun_plot %>%
  filter(!is.na(sample)) %>%
  arrange(sample, desc(SC_label)) %>%
  group_by(sample, group) %>%
  mutate(
    ymax = cumsum(abundance),
    ymin = lag(ymax, default = 0),
    ymid = (ymin + ymax) / 2
  ) %>%
  ungroup() %>%
  inner_join(culture %>% select(sample, culture_SC), by = "sample") %>%
  mutate(culture_label = sc_labels[culture_SC]) %>%
  filter(SC_label == culture_label)

p1 <- ggplot(abun_plot %>% filter(!is.na(sample)),
             aes(x = sample, y = abundance, fill = SC_label)) +
  geom_bar(stat = "identity", width = 0.85) +
  geom_point(data = culture_markers, aes(x = sample, y = ymid),
             shape = 4, size = 3, stroke = 1.2, color = "white", inherit.aes = FALSE) +
  scale_fill_manual(values = c25, name = "Sequence cluster\n(Oxford ST)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "Relative A. baumannii SC abundance\n(mSWEEP, x = culture-reported ST)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(base_dir, "v3_plot1_sc_abundances.pdf"), p1, width = 14, height = 6)
ggsave(file.path(base_dir, "v3_plot1_sc_abundances.png"), p1, width = 14, height = 6, dpi = 300)

# =============================================================================
# Plot 2: A. baumannii read count per sample (Kraken2-classified)
# =============================================================================

read_data <- metadata %>%
  left_join(sample_info %>% select(SRR, aligned, total), by = "SRR") %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    total = as.numeric(total)
  )

p2 <- ggplot(read_data %>% filter(!is.na(sample)),
             aes(x = sample, y = total + 1)) +
  geom_bar(stat = "identity", width = 0.8, fill = "steelblue") +
  scale_y_log10(labels = scales::comma, expand = expansion(mult = c(0, 0.05))) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "A. baumannii reads\n(Kraken2-classified, log10)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(base_dir, "v3_plot2_abau_reads.pdf"), p2, width = 14, height = 4)
ggsave(file.path(base_dir, "v3_plot2_abau_reads.png"), p2, width = 14, height = 4, dpi = 300)

# =============================================================================
# Plot 3: Number of SCs per sample (co-infection complexity)
# =============================================================================

n_scs_plot <- abun_filtered %>%
  filter(abundance > 0.05) %>%
  group_by(SRR, sample, group) %>%
  summarise(n_scs = n(), .groups = "drop") %>%
  # add samples with 0 SCs
  right_join(metadata, by = c("SRR", "sample", "group")) %>%
  mutate(
    n_scs = replace_na(n_scs, 0),
    sample = factor(sample, levels = sample_order)
  )

# p3a: SC count for ALL samples (darkorange, for figure 1)
p3a <- ggplot(n_scs_plot %>% filter(!is.na(sample)),
              aes(x = sample, y = n_scs)) +
  geom_bar(stat = "identity", width = 0.8, fill = "steelblue") +
  scale_y_continuous(breaks = seq(0, 10, 1), expand = c(0, 0)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "Number of A. baumannii SCs\n(abundance > 5%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(base_dir, "v3_plot3a_sc_count_all.pdf"), p3a, width = 14, height = 4)
ggsave(file.path(base_dir, "v3_plot3a_sc_count_all.png"), p3a, width = 14, height = 4, dpi = 300)

# p3b: SC count for only samples with enough A. baumannii reads (steelblue, for figure 2)
n_scs_plot_filtered <- n_scs_plot %>%
  filter(!SRR %in% low_read_samples)

p3b <- ggplot(n_scs_plot_filtered %>% filter(!is.na(sample)),
              aes(x = sample, y = n_scs)) +
  geom_bar(stat = "identity", width = 0.8, fill = "steelblue") +
  scale_y_continuous(breaks = seq(0, 10, 1), expand = c(0, 0)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "Number of A. baumannii SCs\n(abundance > 5%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(base_dir, "v3_plot3b_sc_count_filtered.pdf"), p3b, width = 14, height = 4)
ggsave(file.path(base_dir, "v3_plot3b_sc_count_filtered.png"), p3b, width = 14, height = 4, dpi = 300)

# =============================================================================
# Plot 4: Heatmap of SC abundances
# =============================================================================

# Only include SCs that appear >1% in at least one sample
relevant_scs <- abun_filtered %>%
  filter(abundance > 0.01, SC %in% major_scs) %>%
  pull(SC) %>%
  unique()

heatmap_data <- abun_filtered %>%
  filter(SC %in% relevant_scs) %>%
  mutate(
    SC_label = sc_labels[SC],
    sample = factor(sample, levels = sample_order),
    abundance_pct = abundance * 100
  )

# Complete grid for all sample-SC combinations
all_combos <- expand.grid(
  sample = factor(sample_order, levels = sample_order),
  SC_label = unique(heatmap_data$SC_label)
)

heatmap_data <- heatmap_data %>%
  select(sample, SC_label, abundance_pct) %>%
  right_join(all_combos, by = c("sample", "SC_label")) %>%
  left_join(metadata %>% select(sample, group), by = "sample") %>%
  mutate(abundance_pct = replace_na(abundance_pct, 0))

# Order SCs by total abundance
sc_order <- heatmap_data %>%
  group_by(SC_label) %>%
  summarise(total = sum(abundance_pct), .groups = "drop") %>%
  arrange(total) %>%
  pull(SC_label)

heatmap_data$SC_label <- factor(heatmap_data$SC_label, levels = sc_order)

p4 <- ggplot(heatmap_data, aes(x = sample, y = SC_label, fill = abundance_pct)) +
  geom_tile(color = "black", linewidth = 0.01) +
  scale_fill_gradient(low = "white", high = "darkorange",
                      limits = c(0, 100), name = "Abundance (%)") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = "A. baumannii SC (Oxford ST)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 13)
  )

ggsave(file.path(base_dir, "heatmap.pdf"), p4, width = 14, height = 4)
ggsave(file.path(base_dir, "v3_plot4_heatmap.png"), p4, width = 14, height = 6, dpi = 300)

# =============================================================================
# Plot 5: Combined figure (plots 1 + 2 + 3 stacked)
# =============================================================================

p_combined <- p1 / p2 / p3 +
  plot_layout(heights = c(3, 1.5, 1.5)) +
  plot_annotation(
    title = "A. baumannii sub-lineage diversity in CRAB VAP tracheal aspirate metagenomics",
    subtitle = "Xiao et al. 2022 (PRJNA681291) — v3 species-filtered mSWEEP analysis",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )

ggsave(file.path(base_dir, "v3_plot5_combined.pdf"), p_combined, width = 14, height = 14)
ggsave(file.path(base_dir, "v3_plot5_combined.png"), p_combined, width = 14, height = 14, dpi = 300)

# =============================================================================
# Plot 6: Co-infection summary — culture ST vs mSWEEP dominant + minor SCs
# =============================================================================

coinfection_data <- abun_filtered %>%
  filter(abundance > 0.05, !is.na(sample)) %>%
  mutate(SC_label = sc_labels[SC]) %>%
  group_by(sample, group) %>%
  summarise(
    n_scs = n(),
    dominant = SC_label[which.max(abundance)],
    dominant_pct = max(abundance) * 100,
    all_scs = paste(sprintf("%s (%.0f%%)", SC_label, abundance * 100), collapse = " + "),
    .groups = "drop"
  ) %>%
  arrange(group, -n_scs)

# Print co-infection summary to console
cat("\n=== Co-infection summary (v3, SCs > 5%) ===\n\n")
cat("Samples with multiple SCs:\n")
coinfection_data %>%
  filter(n_scs > 1) %>%
  mutate(line = sprintf("  %s (%s): %s", sample, group, all_scs)) %>%
  pull(line) %>%
  cat(sep = "\n")

cat("\n\nSingle-strain samples:\n")
coinfection_data %>%
  filter(n_scs == 1) %>%
  mutate(line = sprintf("  %s (%s): %s", sample, group, all_scs)) %>%
  pull(line) %>%
  cat(sep = "\n")

cat("\n\nSamples with no A. baumannii:\n")
no_abau <- metadata %>%
  filter(!sample %in% coinfection_data$sample)
cat(paste("  ", no_abau$sample, "(", no_abau$group, ")"), sep = "\n")

# =============================================================================
# Plot 7: Kraken2 species composition
# =============================================================================

kraken_file <- file.path(base_dir, "kraken2_species.tsv")
if (file.exists(kraken_file)) {
  kraken <- read_tsv(kraken_file, show_col_types = FALSE) %>%
    left_join(metadata, by = "SRR")

  # Top 15 species across all samples, always include A. baumannii
  top_species <- kraken %>%
    group_by(species) %>%
    summarise(total_reads = sum(reads), .groups = "drop") %>%
    arrange(-total_reads) %>%
    head(15) %>%
    pull(species)

  # Ensure A. baumannii is always included
  if (!"Acinetobacter baumannii" %in% top_species) {
    top_species <- c("Acinetobacter baumannii", top_species[1:14])
  }


  kraken_colors_16 <- c("darkred",
    "darkorange",
    "steelblue",
    "darkgreen",
    "darkgrey",
    "purple",
    "brown",
    "cyan",
    "magenta",
    "goldenrod",
    "darkblue",
    "darkcyan",
    "darkmagenta",
    "darkgoldenrod",
    "darkslategray",
    "darkolivegreen",
    "lightgrey")  # for "Other"

  kraken_plot <- kraken %>%
    mutate(species_label = ifelse(species %in% top_species, species, "Other species (<0.5%)")) %>%
    group_by(SRR, sample, group, species_label) %>%
    summarise(pct = sum(pct), .groups = "drop") %>%
    # Add "Other species" row to fill to 100%
    group_by(SRR, sample, group) %>%
    mutate(shown_total = sum(pct)) %>%
    ungroup()

  # Add remainder as "Other species (<0.5%)" for each sample
  remainder <- kraken_plot %>%
    distinct(SRR, sample, group, shown_total) %>%
    mutate(
      species_label = "Other species (<0.5%)",
      pct = pmax(0, 100 - shown_total)
    ) %>%
    select(-shown_total)

  kraken_plot <- kraken_plot %>%
    select(-shown_total) %>%
    bind_rows(remainder) %>%
    group_by(SRR, sample, group, species_label) %>%
    summarise(pct = sum(pct), .groups = "drop") %>%
    mutate(sample = factor(sample, levels = sample_order))

  # Order species: Other first (bottom), Homo sapiens second, then all others alphabetically
  # but Klebsiella before Acinetobacter, and A. baumannii last (top)
  all_species_labels <- unique(kraken_plot$species_label)
  other_species <- sort(setdiff(all_species_labels,
                                c("Other species (<0.5%)", "Homo sapiens",
                                  "Acinetobacter baumannii")))
  # move Klebsiella species before Acinetobacter species
  kleb <- grep("Klebsiella", other_species, value = TRUE)
  acin <- grep("Acinetobacter", other_species, value = TRUE)
  rest <- setdiff(other_species, c(kleb, acin))
  other_species <- c(rest, kleb, acin)

  species_level_order <- c("Other species (<0.5%)", "Homo sapiens",
                           other_species, "Acinetobacter baumannii")
  kraken_plot <- kraken_plot %>%
    mutate(species_label = factor(species_label, levels = species_level_order))


  # Colors for species
  # Use distinct qualitative palette for other species, fixed colors for key ones
  qual_cols <- c("#228833", "#4477AA", "#AA3377", "#CCBB44", "#66CCEE",
                 "#EE6677", "#004488", "#997700", "#882255", "#559988",
                 "#DDCC77", "#CC6677", "#117733", "#332288")
  species_cols <- setNames(qual_cols[seq_along(other_species)], other_species)
  species_cols["Homo sapiens"] <- "lightgrey"
  species_cols["Other species (<0.5%)"] <- "darkgrey"
  species_cols["Acinetobacter baumannii"] <- "steelblue"
  if ("Klebsiella pneumoniae" %in% names(species_cols)) {
    species_cols["Klebsiella pneumoniae"] <- "darkred"
  }

  p7 <- ggplot(kraken_plot %>% filter(!is.na(sample)),
               aes(x = sample, y = pct, fill = species_label)) +
    geom_bar(stat = "identity", width = 0.85) +

    scale_fill_manual(values = species_cols, name = "Species (>0.5%)") +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    labs(
      x = NULL,
      y = "% of classified reads\n(Kraken2, species > 0.5%)"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.y = element_text(size = 12),
      legend.position = "right",
      legend.text = element_text(size = 9, face = "italic"),
      panel.grid.major.x = element_blank(),
      strip.text = element_text(face = "bold", size = 13)
    )

  ggsave(file.path(base_dir, "kraken_species.pdf"), p7, width = 16, height = 6)
  ggsave(file.path(base_dir, "v3_plot7_kraken_species.png"), p7, width = 16, height = 6, dpi = 300)
  cat("Kraken2 species plot saved.\n")
} else {
  cat("Kraken2 species file not found, skipping plot 7.\n")
}

# =============================================================================
# Combined Figure 1: Kraken species (A) + A. baumannii reads (B)
# =============================================================================

if (exists("p7")) {
  fig1 <- p7 / p2 / p3a +
    plot_layout(heights = c(3, 1.5, 1.5)) +
    plot_annotation(tag_levels = "A")

  ggsave(file.path(base_dir, "figure1_kraken_abau_reads.pdf"), fig1, width = 16, height = 14)
  ggsave(file.path(base_dir, "figure1_kraken_abau_reads.png"), fig1, width = 16, height = 14, dpi = 300)
  cat("Figure 1 saved.\n")
}

# =============================================================================
# Combined Figure 2: SC abundances (A) + SC count filtered (B)
# =============================================================================

fig2 <- p1 / p3b +
  plot_layout(heights = c(3, 1.5)) +
  plot_annotation(tag_levels = "A")

ggsave(file.path(base_dir, "figure2_sc_sc_count.pdf"), fig2, width = 16, height = 10)
ggsave(file.path(base_dir, "figure2_sc_sc_count.png"), fig2, width = 16, height = 10, dpi = 300)
cat("Figure 2 saved.\n")

cat("\nAll plots saved to", base_dir, "\n")
