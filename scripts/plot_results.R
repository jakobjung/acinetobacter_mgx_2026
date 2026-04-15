#!/usr/bin/env Rscript

# =============================================================================
# Plotting script for mSWEEP v2 results (Xiao et al. CRAB VAP dataset)
# =============================================================================

library(tidyverse)
library(patchwork)

# --- Load data ---------------------------------------------------------------

base_dir <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "..", "analysis")
# If not running in RStudio, set manually:
# base_dir <- "analysis"

abundances <- read_tsv(file.path(base_dir, "all_abundances.tsv"),
                       col_names = c("SRR", "SC", "abundance"),
                       show_col_types = FALSE)

metadata <- read_tsv(file.path(base_dir, "sample_metadata.tsv"),
                     show_col_types = FALSE)

sample_info <- read_tsv(file.path(base_dir, "sample_info.tsv"),
                        show_col_types = FALSE)

qc_meta <- read_tsv(file.path(base_dir, "xiao_qc_metadata.tsv"),
                    show_col_types = FALSE)

sc_label_map <- read_tsv(file.path(base_dir, "sc_labels.tsv"),
                         show_col_types = FALSE)

# --- Merge metadata ----------------------------------------------------------

abun <- abundances %>%
  left_join(metadata, by = "SRR") %>%
  left_join(sample_info %>% select(SRR, aligned, total), by = "SRR")

# --- Define SC categories and colors ----------------------------------------

# Build label lookup from sc_labels.tsv
sc_labels <- setNames(sc_label_map$short_label, sc_label_map$SC)
major_scs <- sc_label_map$SC

# Color palette (Paul Tol's qualitative + greys for minor SCs)
sc_colors <- c(
  "Oxf ST940" = "#4477AA",
  "Oxf ST940-2" = "#6699CC",
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
  "Pas ST151" = "#AAAAAA",
  "Pas ST2 (GC2)" = "#004488",
  "Oxf ST893" = "#997700",
  "Pas ST738" = "#559988",
  "Unknown" = "#CCCCCC",
  "Other (<1%)" = "#EEEEEE"
)

# --- Prepare abundance data for plotting -------------------------------------

abun_plot <- abun %>%
  mutate(SC_label = ifelse(SC %in% major_scs & abundance > 0.01,
                           sc_labels[SC],
                           "Other (<1%)")) %>%
  group_by(SRR, sample, group, SC_label) %>%
  summarise(abundance = sum(abundance), .groups = "drop") %>%
  mutate(SC_label = factor(SC_label, levels = rev(names(sc_colors))))

# Order samples by group then by ST940 abundance
sample_order <- abun_plot %>%
  filter(SC_label == "Oxf ST940") %>%
  arrange(group, -abundance) %>%
  pull(sample)

abun_plot <- abun_plot %>%
  mutate(sample = factor(sample, levels = sample_order))

# =============================================================================
# Plot 1: Stacked bar chart of SC abundances
# =============================================================================

p1 <- ggplot(abun_plot, aes(x = sample, y = abundance, fill = SC_label)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_fill_manual(values = sc_colors, name = "Sequence cluster\n(Oxford ST)") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "A. baumannii SC-level abundance per sample (mSWEEP v2)",
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "right",
    legend.text = element_text(size = 8),
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(base_dir, "plot1_sc_abundances.pdf"), p1, width = 14, height = 6)
ggsave(file.path(base_dir, "plot1_sc_abundances.png"), p1, width = 14, height = 6, dpi = 300)

# =============================================================================
# Plot 2: Alignment rate per sample (% reads mapping to A. baumannii)
# =============================================================================

align_data <- metadata %>%
  left_join(sample_info %>% select(SRR, aligned, total), by = "SRR") %>%
  mutate(
    pct_aligned = aligned / total * 100,
    sample = factor(sample, levels = sample_order)
  )

p2 <- ggplot(align_data, aes(x = sample, y = pct_aligned, fill = group)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("CRAB-I" = "#E63946", "CRAB-C" = "#457B9D", "CRAB-N" = "#A8DADC"),
                    name = "Group") +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "Percentage of reads pseudoaligned to A. baumannii reference panel",
    x = NULL,
    y = "% reads aligned"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(base_dir, "plot2_alignment_rates.pdf"), p2, width = 14, height = 4)
ggsave(file.path(base_dir, "plot2_alignment_rates.png"), p2, width = 14, height = 4, dpi = 300)

# =============================================================================
# Plot 3: Number of SCs per sample (co-infection complexity)
# =============================================================================

n_scs <- abun %>%
  filter(abundance > 0.05) %>%
  group_by(SRR, sample, group) %>%
  summarise(n_scs = n(), .groups = "drop") %>%
  mutate(sample = factor(sample, levels = sample_order))

p3 <- ggplot(n_scs, aes(x = sample, y = n_scs, fill = group)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("CRAB-I" = "#E63946", "CRAB-C" = "#457B9D", "CRAB-N" = "#A8DADC"),
                    name = "Group") +
  scale_y_continuous(breaks = seq(0, 10, 1), expand = c(0, 0)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "Number of A. baumannii SCs detected per sample (abundance > 5%)",
    x = NULL,
    y = "Number of SCs"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(base_dir, "plot3_sc_count.pdf"), p3, width = 14, height = 4)
ggsave(file.path(base_dir, "plot3_sc_count.png"), p3, width = 14, height = 4, dpi = 300)

# =============================================================================
# Plot 4: Host proportion (from paper QC data)
# =============================================================================

qc_plot <- qc_meta %>%
  mutate(
    sample = factor(sample, levels = sample_order),
    microbial_pct = (1 - host_proportion) * 100
  )

p4 <- ggplot(qc_plot, aes(x = sample, y = microbial_pct, fill = group)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = c("CRAB-I" = "#E63946", "CRAB-C" = "#457B9D", "CRAB-N" = "#A8DADC"),
                    name = "Group") +
  scale_y_continuous(expand = c(0, 0)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "Microbial read fraction after human read removal",
    x = NULL,
    y = "% non-human reads"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(base_dir, "plot4_host_proportion.pdf"), p4, width = 14, height = 4)
ggsave(file.path(base_dir, "plot4_host_proportion.png"), p4, width = 14, height = 4, dpi = 300)

# =============================================================================
# Plot 5: Heatmap of top SC abundances across samples
# =============================================================================

heatmap_data <- abun %>%
  filter(SC %in% major_scs, abundance > 0.001) %>%
  mutate(
    SC_label = sc_labels[SC],
    sample = factor(sample, levels = sample_order),
    abundance_pct = abundance * 100
  )

# add zero rows for missing SCs in each sample to ensure complete heatmap
all_combinations <- expand.grid(
    sample = levels(heatmap_data$sample),
    SC_label = factor(names(sc_colors), levels = names(sc_colors))
)

heatmap_data <- heatmap_data %>%
  select(sample, SC_label, abundance_pct) %>%
  right_join(all_combinations, by = c("sample", "SC_label")) %>%
  left_join(metadata %>% select(sample, group), by = "sample") %>%
  mutate(abundance_pct = replace_na(abundance_pct, 0))

# ok change order so that st940 is at the top, ST208 next, then the rest in descending order of overall abundance
sc_order <- heatmap_data %>%
  group_by(SC_label) %>%
  summarise(total_abundance = sum(abundance_pct), .groups = "drop") %>%
  arrange(total_abundance) %>%
  pull(SC_label)

heatmap_data$SC_label <- factor(heatmap_data$SC_label, levels = sc_order)

# also order samples as before
heatmap_data$sample <- factor(heatmap_data$sample, levels = sample_order)



library(viridis)
p5 <- ggplot(heatmap_data, aes(x = sample, y = SC_label, fill = abundance_pct)) +
  geom_tile(color = "black", linewidth = 0.01) +
  # scale_fill_viridis(option = "viridis", name = "Abundance (%)", na.value = "white") +
  scale_fill_gradient(low = "white", high = "darkorange", limits = c(0, 100)) +
  facet_grid(~ group, scales = "free_x", space = "free_x") +
  labs(
    title = "SC abundance heatmap across all samples",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )
# a

ggsave(file.path(base_dir, "plot5_heatmap.pdf"), p5, width = 14, height = 6)
ggsave(file.path(base_dir, "plot5_heatmap.png"), p5, width = 14, height = 6, dpi = 300)

# =============================================================================
# Plot 6: Combined figure (plots 1 + 2 + 3 stacked)
# =============================================================================

p_combined <- p1 / p2 / p3 +
  plot_layout(heights = c(3, 1.5, 1.5)) +
  plot_annotation(
    title = "A. baumannii sub-lineage diversity in CRAB VAP tracheal aspirate metagenomics",
    subtitle = "Xiao et al. 2022 (PRJNA681291) — mSWEEP v2 analysis",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "grey40")
    )
  )

ggsave(file.path(base_dir, "plot6_combined.pdf"), p_combined, width = 14, height = 14)
ggsave(file.path(base_dir, "plot6_combined.png"), p_combined, width = 14, height = 14, dpi = 300)

cat("All plots saved to", base_dir, "\n")
