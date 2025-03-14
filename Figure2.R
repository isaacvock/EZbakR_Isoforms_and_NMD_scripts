### PURPOSE OF THIS SCRIPT
## Reproduce Figure 2


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpointdensity)
library(ggh4x)
library(EZbakR)


##### User input #####

fastq2EZbakR_EZbdo <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/refseq/RefSeq_EZbakRFit_withgenewide.rds"
GRANDSLAM_EZbdo <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Comparisons/GRAND_SLAM/EZbakRFit_GRANDSLAM_11j.rds"

figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Paper/Figures/Figure_2/"


##### Functions used throughout #####

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


make_volcano_plot <- function(res, label_df = NULL,
                              nolabel = FALSE,
                              FDR = 0.05, L2FC_cut = 1,
                              ps = 0.5, lw = 0.5, ts = 2.5,
                              t_hj = .5, t_vj = -1.5,
                              lab_pshape = 21, lab_ps = 1, lab_pstroke = 0.5,
                              xlabel = "L2FC(kdeg)", ylabel = "-log10(padj)"){
  
  
  if(nolabel){
    
    res %>%
      mutate(
        sig = factor(case_when(
          padj < FDR & L2FC_kdeg < -abs(L2FC_cut) ~ "Stabilized",
          padj < FDR & L2FC_kdeg > abs(L2FC_cut) ~ "Destabilized",
          .default = "NS"
        ),
        levels = c("NS", "Stabilized", "Destabilized"))
      ) %>%
      arrange(sig) %>%
      ggplot(
        aes(x = L2FC_kdeg,
            y = -log10(padj),
            color = sig)
      ) + 
      geom_point(size = ps) + 
      theme_classic() + 
      xlab(xlabel) + 
      ylab(ylabel) +
      scale_color_manual(
        values = c("Stabilized" = "#56B4E9", 
                   "Destabilized" = "#e18228", 
                   "NS" = "grey")
      ) +
      geom_hline(
        yintercept = -log10(FDR),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      geom_vline(
        xintercept = -abs(L2FC_cut),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      geom_vline(
        xintercept = abs(L2FC_cut),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      theme(legend.position="none")
    
  }else{
    
    res %>%
      mutate(
        sig = factor(case_when(
          padj < FDR & L2FC_kdeg < -abs(L2FC_cut) ~ "Stabilized",
          padj < FDR & L2FC_kdeg > abs(L2FC_cut) ~ "Destabilized",
          .default = "NS"
        ),
        levels = c("NS", "Stabilized", "Destabilized"))
      ) %>%
      arrange(sig) %>%
      ggplot(
        aes(x = L2FC_kdeg,
            y = -log10(padj),
            color = sig)
      ) + 
      geom_point(size = ps) + 
      theme_classic() + 
      xlab(xlabel) + 
      ylab(ylabel) +
      scale_color_manual(
        values = c("Stabilized" = "#56B4E9", 
                   "Destabilized" = "#e18228", 
                   "NS" = "grey")
      ) +
      geom_hline(
        yintercept = -log10(FDR),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      geom_vline(
        xintercept = -abs(L2FC_cut),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      geom_vline(
        xintercept = abs(L2FC_cut),
        color = 'darkred',
        linetype = "dotted",
        linewidth = lw
      ) +
      theme(legend.position="none") +
      # Add label for "SRSF3"
      geom_text(data = label_df, aes(x = L2FC_kdeg, y = -log10(padj), label = gene_name), 
                color = "black", size = ts, hjust = t_hj, vjust = t_vj)  +  
      # Add black ring around the "SRSF3" point
      geom_point(data = label_df, aes(x = L2FC_kdeg, y = -log10(padj)), 
                 shape = lab_pshape, size = lab_ps, fill = NA, color = "black", 
                 stroke = lab_pstroke)  # Black ring around SRSF3
    
  }
  
  
}


# Method comparison ------------------------------------------------------------

#### Load data ####


gs_ez <- readRDS(GRANDSLAM_EZbdo)
ezbdo <- readRDS(fastq2EZbakR_EZbdo)


#### Get comparisons ####

gs_comp <- EZget(gs_ez, type = "comparisons")
ezgene_comp <- EZget(ezbdo, type = "comparisons",
                     features = "XF")
eztx_comp <- EZget(ezbdo, type = "comparisons",
                   features = c("transcript_id"),
                   exactMatch = FALSE)

genes_to_keep <- intersect(gs_comp$XF,
                           ezgene_comp$XF) %>%
  intersect(unique(eztx_comp$XF))

gs_comp <- gs_comp %>%
  filter(XF %in% genes_to_keep)
ezgene_comp <- ezgene_comp %>%
  filter(XF %in% genes_to_keep)
eztx_comp <- eztx_comp %>%
  filter(XF %in% genes_to_keep)


#### Calculate metrics ####

stabilized_genes <- unique(
  c(gs_comp$XF[gs_comp$difference < -log(2)],
    ezgene_comp$XF[ezgene_comp$difference < -log(2)],
    eztx_comp$XF[eztx_comp$difference < -log(2)])
)

stabilized_genes_padj <- unique(
  c(gs_comp$XF[gs_comp$difference < -log(2) & gs_comp$padj < 0.05],
    ezgene_comp$XF[ezgene_comp$difference < -log(2) & ezgene_comp$padj < 0.05],
    eztx_comp$XF[eztx_comp$difference < -log(2) & eztx_comp$padj < 0.05])
)


power_df <- tibble(
  Power = c(
    length(unique(gs_comp$XF[gs_comp$XF %in% stabilized_genes & gs_comp$difference < -log(2)]))/length(stabilized_genes),
    length(unique(ezgene_comp$XF[ezgene_comp$XF %in% stabilized_genes & ezgene_comp$difference < -log(2)]))/length(stabilized_genes),
    length(unique(gs_comp$XF[eztx_comp$XF %in% stabilized_genes & eztx_comp$difference < -log(2)]))/length(stabilized_genes)
  ),
  Power_padj = c(
    length(unique(gs_comp$XF[gs_comp$XF %in% stabilized_genes & gs_comp$difference < -log(2) & gs_comp$padj < 0.05]))/length(stabilized_genes_padj),
    length(unique(ezgene_comp$XF[ezgene_comp$XF %in% stabilized_genes & ezgene_comp$difference < -log(2) & ezgene_comp$padj < 0.05]))/length(stabilized_genes_padj),
    length(unique(gs_comp$XF[eztx_comp$XF %in% stabilized_genes & eztx_comp$difference < -log(2) & eztx_comp$padj < 0.05]))/length(stabilized_genes_padj)
  ),
  Method = factor(c(
    "GRAND-SLAM", "EZbakR (Gene)", "EZbakR (Transcript)"
  ),
  levels = c("EZbakR (Gene)", "GRAND-SLAM" , "EZbakR (Transcript)"))
)

power_df

gpow <- power_df %>%
  ggplot(aes(x = Method,
             y = Power_padj)) + 
  geom_bar(
    stat = "identity",
    fill = "darkgray",
    color = "black",
    linewidth = 0.1,
    width = 0.75
  ) +
  theme_classic() +
  xlab("Method") +
  ylab("Power") +
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10),
    axis.text.x = element_text(angle = 45, hjust=1))


setwd(figure_savedir)  
ggsave(
  filename = "Power_comparison.pdf",
  gpow,
  width = 2.75,
  height = 2.5
)

# Volcano plots ----------------------------------------------------------------

ezbdo <- readRDS(fastq2EZbakR_EZbdo)


##### Transcripts ######

res_df <- EZget(ezbdo,
                type = "comparisons",
                features = "transcript_id",
                exactMatch = FALSE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  )

label_data <- res_df %>%
  filter(XF == "SRSF3") %>%
  dplyr::mutate(
    gene_name = ifelse(difference < -1,
                       "SRSF3-PTC",
                       "SRSF3")
  )

volc_tx <- make_volcano_plot(res_df, label_data,
                             ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + 
  coord_cartesian(ylim = c(0, 400),
                  xlim = c(-4, 4))


##### Genes ######

g_res_df <- EZget(ezbdo,
                  type = "comparisons",
                  features = "XF",
                  exactMatch = TRUE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  ) %>%
  filter(
    XF %in% unique(res_df$XF)
  )

g_res_df %>%
  count(padj < 0.05 & L2FC_kdeg < -1)


g_res_df %>%
  count(padj < 0.05 & L2FC_kdeg > 1)


g_label_data <- g_res_df %>%
  filter(XF == "SRSF3") %>%
  dplyr::mutate(
    gene_name = "SRSF3"
  )

volc_gene <- make_volcano_plot(g_res_df, g_label_data,
                               ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + 
  coord_cartesian(ylim = c(0, 400),
                  xlim = c(-4, 4))

volc_gene_nolabel <- make_volcano_plot(g_res_df, nolabel = TRUE,
                                       ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + 
  coord_cartesian(ylim = c(0, 400),
                  xlim = c(-4, 4))


#### Save plots #####

setwd(figure_savedir)
ggsave(
  filename = "Gene_volcano_plot.pdf",
  plot = volc_gene,
  width = 2.5, 
  height = 2
)
ggsave(
  filename = "Isoform_volcano_plot.pdf",
  plot = volc_tx,
  width = 2.5, 
  height = 2
)


setwd(figure_savedir)
ggsave(
  filename = "Gene_volcano_plot.png",
  plot = volc_gene,
  width = 2.5, 
  height = 2
)
ggsave(
  filename = "Gene_volcano_plot_nolabel.png",
  plot = volc_gene_nolabel,
  width = 2.5, 
  height = 2
)
ggsave(
  filename = "Isoform_volcano_plot.png",
  plot = volc_tx,
  width = 2.5, 
  height = 2
)

# Gene vs. Transcript ----------------------------------------------------------


ezbdo <- readRDS(fastq2EZbakR_EZbdo)

res_df <- EZget(ezbdo,
                type = "comparisons",
                features = "transcript_id",
                exactMatch = FALSE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  )

g_res_df <- EZget(ezbdo,
                  type = "comparisons",
                  features = "XF",
                  exactMatch = TRUE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  )


g_res_df

res_df_XF <- res_df %>%
  dplyr::group_by(XF) %>%
  dplyr::summarise(
    max_tx_difference = min(L2FC_kdeg),
    max_tx_padj = unique(padj[L2FC_kdeg == min(L2FC_kdeg)])
  ) %>%
  dplyr::mutate(
    tx_sig = factor(case_when(
      max_tx_padj < 0.05 & max_tx_difference < -1 ~ "Stabilized",
      max_tx_padj < 0.05 & max_tx_difference > 1 ~ "Destabilized",
      .default = "NS"
    ),
    levels = c("NS", "Stabilized", "Destabilized"))
  )

res_df_XF %>%
  count(tx_sig)

compare_XF <- g_res_df %>%
  inner_join(
    res_df_XF,
    by = "XF"
  )

comp_label <- compare_XF %>%
  filter(XF == "SRSF3")

comp_plot <- compare_XF %>%
  arrange(tx_sig) %>%
  ggplot(
    aes(
      x = max_tx_difference,
      y = L2FC_kdeg
    )
  ) +
  geom_point(size = 0.2,
             aes(color = tx_sig)) + 
  scale_color_manual(
    values = c("Stabilized" = "#56B4E9", 
               "Destabilized" = "#e18228", 
               "NS" = "grey")
  ) +
  theme_classic() + 
  xlab("Transcript min L2FC(kdeg)") + 
  ylab("Gene L2FC(kdeg)") +
  geom_hline(
    yintercept = -1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_hline(
    yintercept = 1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = -1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = 1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  coord_cartesian(
    xlim = c(-4, 4),
    ylim = c(-4, 4)
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)
  ) +
  theme(legend.position="none") +
  geom_text(data = comp_label, aes(x = max_tx_difference, y = L2FC_kdeg, label = XF), 
            color = "black", size = 2.5, hjust = 0.5, vjust = -2.5)  +  
  geom_point(data = comp_label, aes(x = max_tx_difference, y = L2FC_kdeg), 
             shape = 21, size = 1, fill = NA, color = "black", 
             stroke = 0.5
  )

#### Save plots #####

setwd(figure_savedir)
ggsave(
  filename = "Comparison_plot.pdf",
  plot = comp_plot,
  width = 2.5, 
  height = 2
)
ggsave(
  filename = "Comparison_plot.png",
  plot = comp_plot,
  width = 2.5, 
  height = 2
)



# Gene vs. Transcript (GRAND-SLAM) ---------------------------------------------

gs_ez <- readRDS(GRANDSLAM_EZbdo)

ezbdo <- readRDS(fastq2EZbakR_EZbdo)


res_df <- EZget(ezbdo,
                type = "comparisons",
                features = "transcript_id",
                exactMatch = FALSE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  )

g_res_df <- EZget(gs_ez,
                  type = "comparisons",
                  features = "XF",
                  exactMatch = TRUE) %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  )


g_res_df

res_df_XF <- res_df %>%
  dplyr::group_by(XF) %>%
  dplyr::summarise(
    max_tx_difference = min(L2FC_kdeg),
    max_tx_padj = unique(padj[L2FC_kdeg == min(L2FC_kdeg)])
  ) %>%
  dplyr::mutate(
    tx_sig = factor(case_when(
      max_tx_padj < 0.05 & max_tx_difference < -1 ~ "Stabilized",
      max_tx_padj < 0.05 & max_tx_difference > 1 ~ "Destabilized",
      .default = "NS"
    ),
    levels = c("NS", "Stabilized", "Destabilized"))
  )

res_df_XF %>%
  count(tx_sig)

compare_XF <- g_res_df %>%
  inner_join(
    res_df_XF,
    by = "XF"
  )



comp_label <- compare_XF %>%
  filter(XF == "SRSF3")

comp_plot <- compare_XF %>%
  arrange(tx_sig) %>%
  ggplot(
    aes(
      x = max_tx_difference,
      y = L2FC_kdeg
    )
  ) +
  geom_point(size = 0.2,
             aes(color = tx_sig)) + 
  scale_color_manual(
    values = c("Stabilized" = "#56B4E9", 
               "Destabilized" = "#e18228", 
               "NS" = "grey")
  ) +
  theme_classic() + 
  xlab("Transcript min L2FC(kdeg)") + 
  ylab("Gene L2FC(kdeg)") +
  geom_hline(
    yintercept = -1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_hline(
    yintercept = 1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = -1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  geom_vline(
    xintercept = 1,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  coord_cartesian(
    xlim = c(-4, 4),
    ylim = c(-4, 4)
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)
  ) +
  theme(legend.position="none") +
  geom_text(data = comp_label, aes(x = max_tx_difference, y = L2FC_kdeg, label = XF), 
            color = "black", size = 2.5, hjust = 0.5, vjust = -2.5)  +  
  geom_point(data = comp_label, aes(x = max_tx_difference, y = L2FC_kdeg), 
             shape = 21, size = 1, fill = NA, color = "black", 
             stroke = 0.5
  )

#### Save plots #####

setwd(figure_savedir)
ggsave(
  filename = "Comparison_plot_GRANDSLAM.pdf",
  plot = comp_plot,
  width = 2.5, 
  height = 2
)
ggsave(
  filename = "Comparison_plot_GRANDSLAM.png",
  plot = comp_plot,
  width = 2.5, 
  height = 2
)




# SRSF3 kdeg bar plot ----------------------------------------------------------

ezbdo <- readRDS(fastq2EZbakR_EZbdo)

res_df <- EZget(ezbdo,
                type = "kinetics",
                features = "transcript_id",
                exactMatch = FALSE) %>%
  filter(
    sample %in% c("DMSO_8hr_1", "11j_8hr_1") &
      XF == "SRSF3"
  ) %>%
  mutate(
    name = case_when(
      transcript_id == "NR_036610.2" ~ "PTC",
      .default = "no PTC"
    ),
    level = "Isoform"
  ) %>%
  dplyr::select(
    sample, name, level, kdeg, log_kdeg, se_log_kdeg
  )

g_res_df <- EZget(ezbdo,
                  type = "kinetics",
                  features = "XF",
                  exactMatch = TRUE) %>%
  filter(
    sample %in% c("DMSO_8hr_1", "11j_8hr_1") &
      XF == "SRSF3"
  ) %>%
  mutate(
    name = "gene-avg",
    level = "Gene"
  ) %>%
  dplyr::select(
    sample, name, level, kdeg, log_kdeg, se_log_kdeg
  )

bar_df <- bind_rows(
  res_df,
  g_res_df
)


DMSO_bar <- bar_df %>%
  dplyr::mutate(
    name = factor(name,
                  levels = c("no PTC", "PTC", "gene-avg")),
    level = factor(level,
                   levels = c("Isoform", "Gene"))
  ) %>%
  dplyr::filter(sample == "DMSO_8hr_1") %>%
  mutate(
    thalf = log(2) / kdeg,
    thalf_min = log(2) / (kdeg - se_log_kdeg),
    thalf_max = log(2) / (kdeg + se_log_kdeg)
  ) %>%
  ggplot(aes(
    x = name,
    y = log(2) / kdeg
  )) + 
  geom_bar(
    stat = "identity",
    color = 'black',
    fill = 'darkgray',
    linewidth = 0.25
  ) + 
  geom_point(size = 0.2) +
  facet_grid(.~ level, scales = "free_x", space = "free_x", switch = "x") +
  geom_errorbar(
    aes(
      x = name,
      ymin = thalf_max,
      ymax = thalf_min
    ),
    width = 0.3,
    linewidth = 0.25
  ) +
  theme_classic() + 
  xlab("Feature") + 
  ylab("t1/2 (hrs)") +
  coord_cartesian(
    ylim = c(0, 4.5)
  ) +
  theme(
    axis.text=element_text(size=6), #change font size of axis text
    axis.title=element_text(size=8), #change font size of axis titles
    legend.text=element_text(size=6), #change font size of legend text
    legend.title=element_text(size=8),
    axis.title.x=element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


SMG1i_bar <- bar_df %>%
  dplyr::mutate(
    name = factor(name,
                  levels = c("no PTC", "PTC", "gene-avg")),
    level = factor(level,
                   levels = c("Isoform", "Gene"))
  ) %>%
  dplyr::filter(sample == "11j_8hr_1") %>%
  mutate(
    thalf = log(2) / kdeg,
    thalf_min = log(2) / (kdeg - se_log_kdeg),
    thalf_max = log(2) / (kdeg + se_log_kdeg)
  ) %>%
  ggplot(aes(
    x = name,
    y = log(2) / kdeg
  )) + 
  geom_bar(
    stat = "identity",
    color = 'black',
    fill = 'darkgray',
    linewidth = 0.25
  ) + 
  facet_grid(.~ level, scales = "free_x", space = "free_x", switch = "x") +
  geom_point(size = 0.2) +
  geom_errorbar(
    aes(
      x = name,
      ymin = thalf_max,
      ymax = thalf_min
    ),
    width = 0.3,
    linewidth = 0.25
  ) +
  theme_classic() + 
  xlab("Feature") + 
  ylab("t1/2 (hrs)") +
  coord_cartesian(
    ylim = c(0, 4.5)
  ) +
  theme(
    axis.text=element_text(size=6), #change font size of axis text
    axis.title=element_text(size=8), #change font size of axis titles
    legend.text=element_text(size=6), #change font size of legend text
    legend.title=element_text(size=8),
    axis.title.x=element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

DMSO_bar
SMG1i_bar


##### Save plots #####

setwd(figure_savedir)
ggsave(
  filename = "DMSO_half_lives.pdf",
  plot = DMSO_bar,
  width = 1.7, 
  height = 1
)
ggsave(
  filename = "SMG1i_half_lives.pdf",
  plot = SMG1i_bar,
  width = 1.7, 
  height = 1
)

