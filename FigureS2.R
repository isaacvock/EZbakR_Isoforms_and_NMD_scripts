### PURPOSE OF THIS SCRIPT
## Reproduce supplemental Figure 2


# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(MASS)
library(tidyr)
library(EZbakR)
library(DESeq2)



##### User input #####

EZbdo_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/refseq/RefSeq_EZbakRFit_withgenewide.rds"

figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Paper/Figures/Figure_2/"

##### FUNCTIONS USED THROUGHOUT #####

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


#' ps = point size for scatterplot
#' lw = linewidth for y = x line
compare_replicates <- function(df, x, y,
                               xlabel = gsub("_", " ", x),
                               ylabel = gsub("_", " ", y),
                               xscale = 1,
                               yscale = 1,
                               color_var = "density",
                               n = 100,
                               ps = 0.6,
                               lw = 0.8,
                               lw_lm = 0.5,
                               add_yeqx = TRUE,
                               add_regression = FALSE,
                               log_scale = FALSE){
  
  
  if(color_var == "density"){
    vars_to_keep <- c(x, y)
  }else{
    vars_to_keep <- c(x, y, color_var)
    df <- df %>%
      arrange(!!color_var)
  }
  
  
  if(!log_scale){
    
    df <- df %>%
      dplyr::select(!!vars_to_keep) %>%
      na.omit() %>%
      dplyr::mutate(
        density = get_density(
          x = !!dplyr::sym(x)*xscale,
          y = !!dplyr::sym(y)*yscale,
          n = n
        )
      )
    
    
    gcorr <- df  %>%
      ggplot(
        aes(x = !!dplyr::sym(x)*xscale,
            y = !!dplyr::sym(y)*yscale)
      )
    
  }else{
    
    df <- df %>%
      dplyr::select(!!vars_to_keep) %>%
      na.omit() %>%
      dplyr::mutate(
        density = get_density(
          x = log(!!dplyr::sym(x)*xscale),
          y = log(!!dplyr::sym(y)*yscale),
          n = n
        )
      )
    
    
    gcorr <- df  %>%
      ggplot(
        aes(x = log(!!dplyr::sym(x)*xscale),
            y = log(!!dplyr::sym(y)*yscale))
      )
    
  }
  
  
  gcorr <- gcorr + 
    geom_point(aes(color = !!dplyr::sym(color_var)),
               size = ps) + 
    theme_classic() +
    scale_color_viridis_c() +
    xlab(xlabel) + 
    ylab(ylabel)
  
  if(add_yeqx){
    
    gcorr <- gcorr +
      geom_abline(
        slope = 1,
        intercept = 0,
        color = 'darkred',
        linetype = 'dotted',
        linewidth = lw
      )
    
  }
  
  if(add_regression){
    
    gcorr <- gcorr +                           
      stat_smooth(method = "lm", 
                  formula = y ~ x, 
                  geom = "smooth",
                  linewidth = lw_lm) 
    
  }
  
  gcorr
  
  
}


# How much does kdeg explain TPM -----------------------------------------------

ez_refseq <- readRDS(EZbdo_path)

TPMs <- ez_refseq$readcounts$isoform_quant_rsem %>%
  dplyr::filter(grepl("DMSO", sample)) %>%
  dplyr::group_by(transcript_id) %>%
  summarise(
    avg_TPM = mean(TPM)
  )

avg_kdegs <- ez_refseq$averages$averages1

XF_to_TID <- avg_kdegs %>%
  dplyr::select(XF, transcript_id) %>%
  dplyr::distinct()


TPMs <- XF_to_TID %>%
  inner_join(
    TPMs,
    by = "transcript_id"
  )


TPMs <- TPMs %>%
  group_by(XF) %>%
  mutate(
    IsoPct = avg_TPM / sum(avg_TPM),
    niso = dplyr::n()
  )


combined_TPMs <- TPMs %>%
  inner_join(
    avg_kdegs,
    by = c("XF", "transcript_id")
  )



gTPM <- combined_TPMs %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    sd_treatmentDMSO_posterior < 0.15
  ) %>%
  dplyr::mutate(
    density = get_density(
      x = mean_treatmentDMSO,
      y = log(avg_TPM),
      n = 200
    )
  ) %>%
  ggplot(
    aes(
      x = mean_treatmentDMSO,
      y = log(avg_TPM),
      color = density
    )
  ) + 
  geom_point(
    size = 0.2
  ) + 
  scale_color_viridis_c() + 
  theme_classic() + 
  xlab("DMSO log(kdeg)") + 
  ylab("DMSO log(TPM)") +
  geom_smooth(method=lm,
              color = "darkred",
              linewidth = 0.3,
              formula = y ~ x)+ 
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")


gIsoPct <- combined_TPMs %>%
  dplyr::ungroup() %>%
  dplyr::filter(
    sd_treatmentDMSO_posterior < 0.15 &
      niso > 1
  ) %>%
  dplyr::mutate(
    density = get_density(
      x = mean_treatmentDMSO,
      y = EZbakR:::logit(IsoPct),
      n = 200
    )
  ) %>%
  ggplot(
    aes(
      x = mean_treatmentDMSO,
      y = EZbakR:::logit(IsoPct),
    )
  ) + 
  geom_point(aes(color = density),
             size = 0.2) + 
  scale_color_viridis_c() + 
  theme_classic() + 
  xlab("DMSO log(kdeg)") + 
  ylab("DMSO logit(IsoPct)") +
  geom_smooth(method=lm,
              color = "darkred",
              linewidth = 0.3,
              formula = y ~ x)+ 
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



### calculate R^2's
TPMfit <- lm(log(avg_TPM) ~ mean_treatmentDMSO, 
             data = combined_TPMs %>%
               dplyr::filter(
                 sd_treatmentDMSO_posterior < 0.15
               )) %>%
  summary()

IsoPctfit <- lm(EZbakR:::logit(IsoPct) ~ mean_treatmentDMSO, 
                data = combined_TPMs %>%
                  dplyr::filter(
                    sd_treatmentDMSO_posterior < 0.15 &
                      niso > 1
                  )) %>%
  summary()



Rdf <- tibble(
  R2 = c(TPMfit$r.squared,
         IsoPctfit$r.squared),
  Variable = factor(
    c("TPM", "IsoPct"),
    levels = c("TPM", "IsoPct")
  )
)


gRs <- Rdf %>%
  ggplot(
    aes(x = Variable,
        y = R2)
  ) +
  geom_bar(
    stat = "identity",
    color = "black",
    fill = "darkgray",
    linewidth = 0.3,
    width = 0.5
  ) + 
  theme_classic() + 
  xlab("Quantity") + 
  ylab("Fraction explaind by kdeg") + 
  coord_cartesian(
    ylim = c(0, 1)
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none")



setwd(figure_savedir)
ggsave(
  filename = "TPM_vs_kdeg.pdf",
  plot = gTPM,
  width = 2.5,
  height = 2.1
)
ggsave(
  filename = "IsoPct_vs_kdeg.pdf",
  plot = gIsoPct,
  width = 2.5,
  height = 2.1
)
ggsave(
  filename = "Percent_kdeg_explains.pdf",
  plot = gRs,
  width = 2.1,
  height = 2.1
)




# DE vs. DK --------------------------------------------------------------------

##### Load data #####

ez_refseq <- readRDS(EZbdo_path)

DK <- EZget(ez_refseq,
            type = "comparisons",
            features = "transcript_id",
            exactMatch = FALSE)

##### Run DESeq2 #####

isocnt <- ez_refseq$readcounts$isoform_quant_rsem

samps <- unique(isocnt$sample)
replace_list <- vector(mode = "list",
                       length = length(samps))
for(i in seq_along(replace_list)){
  
  replace_list[[i]] <- list(0)
  
}

names(replace_list) <- samps

isocnt_mat <- isocnt %>%
  dplyr::select(sample, transcript_id, expected_count) %>%
  pivot_wider(names_from = "sample",
              values_from = "expected_count") %>%
  replace_na(replace = replace_list)

isocnt_matrix <- isocnt_mat %>%
  dplyr::select(-transcript_id) %>%
  as.matrix() %>%
  round()
rownames(isocnt_matrix) <- isocnt_mat$transcript_id


metadf_deseq <- data.frame(
  condition = factor(c("DMSO", "DMSO", "DMSO", "11j",
                       "11j", "11j", "DMSO", "11j"))
)


dds <- DESeqDataSetFromMatrix(
  countData = isocnt_matrix,
  colData = metadf_deseq,
  design = ~ condition
)

dds <- DESeq(dds)
resultsNames(dds)
res <- lfcShrink(dds, coef="condition_DMSO_vs_11j")
res_uns <- results(dds, name = "condition_DMSO_vs_11j")

res_df <- as_tibble(res)
res_df$transcript_id <- rownames(res)

res_uns_df <- as_tibble(res_uns)
res_uns_df$transcript_id <- rownames(res_uns)


##### Compare DEseq2 to EZbakR #####

compare_df <- DK %>%
  inner_join(res_df %>%
               dplyr::rename(DE_padj = padj),
             by = "transcript_id")


gcomp <- compare_df %>%
  mutate(
    density = get_density(
      x = difference * log2(exp(1)),
      y = -log2FoldChange,
      n = 200
    )
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log2FoldChange,
        color = density)
  ) + 
  geom_point(size = 0.1) + 
  theme_classic() +
  scale_color_viridis_c() + 
  xlab("L2FC(kdeg)") + 
  ylab("L2FC(RNA)") +
  geom_abline(
    slope = -1,
    intercept = 0,
    color = 'darkred',
    linetype = "dotted",
    linewidth = 0.3
  ) +
  theme(
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none") +
  coord_cartesian(
    xlim = c(-5, 5),
    ylim = c(-5, 5)
  )




compare_df %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  ) %>%
  mutate(
    conclusion = factor(case_when(
      log2FoldChange < -1 & difference < -log(2) ~ "SandU",
      log2FoldChange > -1 & difference < -log(2) ~ "SandNU",
      log2FoldChange < -1 & difference > -1 ~ "ITE",
      .default = "NSandNU",
    ),
    levels = c("SandU", "SandNU", "ITE", "NSandNU")),
    mechanism = case_when(
      log2FoldChange < 0 ~ (L2FC_kdeg + 0.5)*(-log2FoldChange),
      abs(log2FoldChange) < 0.5 & abs(L2FC_kdeg) < 0.5 ~ 0,
      .default = (L2FC_kdeg - 0.5)*(-log2FoldChange)
    )
  ) %>%
  mutate(
    mechanism = ifelse(
      mechanism < -5, -max(mechanism), mechanism
    )
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log2FoldChange,
        color = mechanism)
  ) + 
  geom_point() + 
  theme_classic() +
  #scale_color_viridis_c() + 
  scale_color_gradient2(mid = "gray90") +
  # scale_color_gradient(low = "darkblue",
  #                      high = "darkred") + 
  # scale_color_manual(
  #   values = c("darkcyan", "darkblue", "darkred", "darkgray")
  # ) +
  xlab("L2FC(kdeg)") + 
  ylab("L2FC(RNA)")



compare_df %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  ) %>%
  mutate(
    conclusion = factor(case_when(
      log2FoldChange < -1 & L2FC_kdeg < -1 ~ "SandU",
      log2FoldChange > -0.25 & L2FC_kdeg < -1 ~ "SandNU",
      log2FoldChange < -1 & L2FC_kdeg > -0.25 ~ "ITE",
      .default = "NSandNU",
    ),
    levels = c("SandU", "SandNU", "ITE", "NSandNU")),
    mechanism = case_when(
      log2FoldChange < 0 ~ (L2FC_kdeg + 0.5)*(-log2FoldChange),
      abs(log2FoldChange) < 0.5 & abs(L2FC_kdeg) < 0.5 ~ 0,
      .default = (L2FC_kdeg - 0.5)*(-log2FoldChange)
    )
  ) %>%
  dplyr::count(conclusion)


gconclusion <- compare_df %>%
  mutate(
    L2FC_kdeg = difference * log2(exp(1))
  ) %>%
  mutate(
    conclusion = factor(case_when(
      log2FoldChange < -1 & L2FC_kdeg < -1 ~ "SandU",
      log2FoldChange > -0.25 & L2FC_kdeg < -1 ~ "SandNU",
      log2FoldChange < -1 & L2FC_kdeg > -0.25 ~ "ITE",
      .default = "NSandNU",
    ),
    levels = c("SandU", "SandNU", "ITE", "NSandNU")),
    mechanism = case_when(
      log2FoldChange < 0 ~ (L2FC_kdeg + 0.5)*(-log2FoldChange),
      abs(log2FoldChange) < 0.5 & abs(L2FC_kdeg) < 0.5 ~ 0,
      .default = (L2FC_kdeg - 0.5)*(-log2FoldChange)
    )
  ) %>%
  mutate(
    mechanism = ifelse(
      mechanism < -5, -max(mechanism), mechanism
    )
  ) %>%
  ggplot(
    aes(x = difference * log2(exp(1)),
        y = -log2FoldChange,
        color = conclusion)
  ) + 
  geom_point(size = 0.2) + 
  theme_classic() +
  scale_color_manual(
    values = c("darkcyan", "darkblue", "darkred", "gray90")
  ) +
  xlab("L2FC(kdeg)") + 
  ylab("L2FC(RNA)") +
  theme(
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) + #change font size of legend title
  theme(legend.position = "none") +
  coord_cartesian(
    xlim = c(-5, 5),
    ylim = c(-5, 5)
  )



setwd(figure_savedir)
ggsave(
  filename = "L2FCkdeg_vs_L2FCRNA.pdf",
  plot = gcomp,
  width = 3.25,
  height = 2.75
)
ggsave(
  filename = "L2FCkdeg_vs_L2FCRNA_conclusions.pdf",
  plot = gconclusion,
  width = 3.25,
  height = 2.75
)
