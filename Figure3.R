### PURPOSE OF THIS SCRIPT
## Reproduce Figure 3


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(tidyr)
library(MASS)
library(readr)
library(data.table)
library(rtracklayer)
library(devtools)
library(glmnet)
library(forcats)
library(EZbakR)


##### User input #####

# JCC score files
jcc_SRonly_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Annotation_gamut/JCCs/JCC_SRonly_trimmed_JCC.rds"
jcc_Mix_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Annotation_gamut/JCCs/JCC_mix_trimmed_JCC.rds"
jcc_Ensembl_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Annotation_gamut/JCCs/ensembl_fixed_genescores.csv"
jcc_RefSeq_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Annotation_gamut/JCCs/refseq_genescores.csv"


# Annotations
  # "Trimmed" because that means AnnotationCleaner output with unsupported flag
Mix_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/Annotations/mix_trimmed.gtf"
SRonly_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/Annotations/SRonly_trimmed.gtf"
RefSeq_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/Annotations/refseq_trimmed.gtf"
Ensembl_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/Annotations/ensembl_corrected_trimmed.gtf"


# RSEM quantifications
rsem_Mix <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/RSEM/mix_trimmed/11j_8hr_1.isoforms.results"
rsem_SRonly <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/RSEM/SRonly_then_trimmed/11j_8hr_1.isoforms.results"
rsem_RefSeq <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/RSEM/RefSeq/11j_8hr_1.isoforms.results"
rsem_Ensembl <-  "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/RSEM/Ensembl_fixed/11j_8hr_1.isoforms.results"


# EZbakRFits
ezbdo_Mix <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/mix_trimmed/Mix_trimmed_EZbakRFit.rds"
ezbdo_SRonly <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/SRonly_trimmed/SRonly_trimmed_EZbakRFit.rds"
ezbdo_RefSeq <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/refseq/RefSeq_EZbakRFit.rds"
ezbdo_Ensembl <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses/refseq/RefSeq_EZbakRFit.rds"


# factR2 transcripts info
factR2_Mix <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/Annotations/factR2/mix_trimmed/mix_trimmed_factR2_transcript.tsv"
factR2_SRonly <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/Annotations/factR2/SRonly_trimmed/SRonly_trimmed_factR2_transcript.tsv"
factR2_RefSeq <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/Annotations/factR2/refseq/hg38_refseq_factR2_transcript.tsv"
factR2_Ensembl <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/Annotations/factR2/ensembl/Hs_ensembl_113 _factR2_transcript.tsv"


# Isoform features for LASSO regression
feature_table_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/ML_features/NMD_tables/Filtered_kdeg_feature_table.csv"



# Path to function_library.R
function_library_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Scripts/ML/function_library.R"

# Directory to save figures in
figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Paper/Figures/Figure_3/"

##### Load functions used throughout #####

source(function_library_path)


# JCC CDF ----------------------------------------------------------------------


jcc_sr_T <- readRDS(jcc_SRonly_path)
jcc_mix_T <- readRDS(jcc_Mix_path)
jcc_ensembl <- read_csv(jcc_Ensembl_path)
jcc_refseq <- read_csv(jcc_RefSeq_path)

jcc_comp <- bind_rows(
  list(
    jcc_sr_T$geneScores %>%
      dplyr::select(jccscore) %>%
      na.omit() %>%
      dplyr::filter(jccscore > 0) %>%
      mutate(assembly = "SRonly_T"),
    jcc_mix_T$geneScores %>%
      dplyr::select(jccscore) %>%
      na.omit() %>%
      dplyr::filter(jccscore > 0) %>%
      mutate(assembly = "mix_T"),
    jcc_ensembl %>%
      dplyr::select(jccscore) %>%
      na.omit() %>%
      dplyr::filter(jccscore > 0) %>%
      mutate(assembly = "Ensembl"),
    jcc_refseq %>%
      dplyr::select(jccscore) %>%
      na.omit() %>%
      dplyr::filter(jccscore > 0) %>%
      mutate(assembly = "RefSeq")
  )
)


jcc_comp %>%
  dplyr::group_by(assembly) %>%
  summarise(fraction_bad = sum(jccscore > 0.4) / n())


jcc_cdf <- jcc_comp %>%
  mutate(
    assembly = factor(assembly,
                      levels = c("Ensembl", "RefSeq", "SRonly_T", "mix_T"))
  ) %>%
  filter(jccscore < 1) %>%
  ggplot(aes(x = jccscore, color = assembly)) +
  stat_ecdf(geom = "line") + 
  theme_classic() + 
  # scale_color_viridis_d() +
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  xlab("JCC") +
  ylab("eCDF") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10))


jcc_cdf_zoom <- jcc_cdf + 
  coord_cartesian(
    xlim = c(0, 0.5),
    ylim = c(0, 1)
  )

jcc_cdf_zoom



setwd(figure_savedir)
ggsave(
  filename = "JCC_CDF.pdf",
  plot = jcc_cdf_zoom,
  width = 2.5, 
  height = 2
)



# Problematic barplot ----------------------------------------------------------

# Function
assess_problematic <- function(gtf, rsem,
                               assembly){
  
  expressed <- rsem[TPM > 1] %>%
    dplyr::select(transcript_id, gene_id)
  
  result <- mcols(gtf) %>%
    as_tibble() %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::inner_join(expressed,
                      by = c("transcript_id", "gene_id")) %>%
    dplyr::summarise(
      fraction_problematic = sum(problematic == "TRUE") / n(),
      isoform_count = n()
    )
  
  return(result %>%
           mutate(assembly = assembly))
}


# Need to filter out isoforms with low TPM (< 1)
mix <- rtracklayer::import(Mix_path)
SRonly <- rtracklayer::import(SRonly_path)
refseq <- rtracklayer::import(RefSeq_path)
ensembl <- rtracklayer::import(Ensembl_path)

rsem_mix_T <- fread(rsem_Mix)
rsem_SRonly_T <- fread(rsem_SRonly)
rsem_refseq <- fread(rsem_RefSeq)
rsem_ensembl <-  fread(rsem_Ensembl)


compare_df <- bind_rows(
  list(
    assess_problematic(mix,
                       rsem_mix_T,
                       "Mix"),
    assess_problematic(SRonly,
                       rsem_SRonly_T,
                       "SRonly"),
    assess_problematic(refseq,
                       rsem_refseq,
                       "RefSeq"),
    assess_problematic(ensembl,
                       rsem_ensembl,
                       "Ensembl")
  )
)


prob_bar <- compare_df %>%
  mutate(
    assembly = factor(assembly,
                      levels = c("Ensembl", "RefSeq", "SRonly", "Mix"))
  ) %>%
  ggplot(
    aes(x = assembly,
        y = fraction_problematic,
        fill = assembly)
  ) +
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.6,
    linewidth = 0.3
  ) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  coord_cartesian(
    ylim = c(0, 0.3)
  ) +
  scale_fill_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  xlab("Assembly") + 
  ylab("Fraction problematic") +
  theme(
    legend.position = "none"
  )

prob_bar


setwd(figure_savedir)
ggsave(
  filename = "Problematic_barplot.pdf",
  plot = prob_bar,
  width = 2.25, 
  height = 2
)


# Number of high confidence NMD isoforms ---------------------------------------


### Annotations with problematic flags for filtering isoforms
mix_T <- rtracklayer::import(Mix_path)
SRonly_T <- rtracklayer::import(SRonly_path)
refseq <- rtracklayer::import(RefSeq_path)
ensembl <- rtracklayer::import(Ensembl_path)

### JCC files for filtering out problematic genes
jcc_sr_T <- readRDS(jcc_SRonly_path)$geneScores
jcc_mix_T <- readRDS(jcc_Mix_path)$geneScores
jcc_ensembl <- read_csv(jcc_Ensembl_path)
jcc_refseq <- read_csv(jcc_RefSeq_path)

### EZbakR conclusions for each
ez_mix_T <- readRDS(ezbdo_Mix) 
ez_mix_T <- ez_mix_T %>% EZget(type = "comparisons", features = "transcript_id",exactMatch = FALSE)
ez_sr_T <- readRDS(ezbdo_SRonly)
ez_sr_T <- ez_sr_T %>% EZget(type = "comparisons", features = "transcript_id",exactMatch = FALSE)
ez_refseq <- readRDS(ezbdo_RefSeq) %>%
  EZget(type = "comparisons", features = "transcript_id",exactMatch = FALSE)
ez_ensembl <- readRDS(ezbdo_Ensembl)


### Figure out what genes have decent data
good_genes_sr_T <- jcc_sr_T %>%
  filter(jccscore > 0 & jccscore < 0.4) %>%
  dplyr::select(gene) %>%
  unlist() %>% unname()

good_genes_mix_T <- jcc_mix_T %>%
  filter(jccscore > 0 & jccscore < 0.4) %>%
  dplyr::select(gene) %>%
  unlist() %>% unname()

good_genes_refseq <- jcc_refseq %>%
  filter(jccscore > 0 & jccscore < 0.4) %>%
  dplyr::select(gene) %>%
  unlist() %>% unname()

good_genes_ensembl <- jcc_ensembl %>%
  filter(jccscore > 0 & jccscore < 0.4) %>%
  dplyr::select(gene) %>%
  unlist() %>% unname()


### Figure out which isoforms are legit

good_isoforms_sr_t <- SRonly_T %>%
  as_tibble() %>%
  dplyr::filter(type == "transcript" & problematic != "TRUE") %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

good_isoforms_mix_t <- mix_T %>%
  as_tibble() %>%
  dplyr::filter(type == "transcript" & problematic != "TRUE") %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

good_isoforms_refseq <- refseq %>%
  as_tibble() %>%
  dplyr::filter(type == "transcript" & problematic != "TRUE") %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()


good_isoforms_ensembl <- ensembl %>%
  as_tibble() %>%
  dplyr::filter(type == "transcript" & problematic != "TRUE") %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

### Figure out which transcripts are high confidence NMD targets in each
NMD_targets_sr_T <- ez_sr_T %>%
  filter(XF %in% good_genes_sr_T &
           difference < -log(2) &
           padj < 0.01 &
           transcript_id %in% good_isoforms_sr_t) %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

NMD_targets_mix_T <- ez_mix_T %>%
  filter(XF %in% good_genes_mix_T &
           difference < -log(2) &
           padj < 0.01 &
           transcript_id %in% good_isoforms_mix_t) %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

NMD_targets_refseq <- ez_refseq %>%
  filter(XF %in% good_genes_refseq &
           difference < -log(2) &
           padj < 0.01 &
           transcript_id %in% good_isoforms_refseq) %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()

NMD_targets_ensembl <- ez_ensembl %>%
  EZget(type = "comparisons") %>%
  filter(XF %in% good_genes_ensembl &
           difference < -log(2) &
           padj < 0.01 &
           transcript_id %in% good_isoforms_ensembl) %>%
  dplyr::select(transcript_id) %>%
  unlist() %>% unname()



### Get counts

NMDcnt_df <- tibble(
  num_NMD = c(
    length(NMD_targets_ensembl),
    length(NMD_targets_refseq),
    length(NMD_targets_mix_T),
    length(NMD_targets_sr_T)
  ),
  assembly = factor(
    c("Ensembl", "RefSeq", "Mix", "SRonly"),
    levels = c("Ensembl", "RefSeq", "SRonly", "Mix")
  )
)


NMD_bar <- NMDcnt_df %>%
  ggplot(
    aes(x = assembly,
        y = num_NMD,
        fill = assembly)
  ) +
  geom_bar(
    stat = "identity",
    color = "black",
    width = 0.6,
    linewidth = 0.3
  ) + 
  scale_fill_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5)) +
  theme(
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12)) +
  xlab("Assembly") + 
  ylab("# high conf. stabilized") +
  theme(legend.position = "none")

NMD_bar


setwd(figure_savedir)
ggsave(
  filename = "NMD_barplot.pdf",
  plot = NMD_bar,
  width = 3.25,
  height = 2.75
)

NMDcnt_df



# Kinetic parameter plot like that made by Justin and Bobby --------------------


##### Load EZbakR analyses #####

avg_mix <- readRDS(ezbdo_Mix) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_sr <- readRDS(ezbdo_SRonly) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_refseq <- readRDS(ezbdo_RefSeq) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_ensembl <- readRDS(ezbdo_Ensembl) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)


comp_mix <- readRDS(ezbdo_Mix) %>%
  EZget(type = "comparisons",
        features = "transcript_id",
        exactMatch = FALSE)

comp_sr <- readRDS(ezbdo_SRonly) %>%
  EZget(type = "comparisons",
        features = "transcript_id",
        exactMatch = FALSE)

comp_refseq <- readRDS(ezbdo_RefSeq) %>%
  EZget(type = "comparisons",
        features = "transcript_id",
        exactMatch = FALSE)

comp_ensembl <- readRDS(ezbdo_Ensembl) %>%
  EZget(type = "comparisons",
        features = "transcript_id",
        exactMatch = FALSE)



##### Load factR2 annotations #####

factr2_mix <- fread(factR2_Mix)
factr2_sronly <- fread(factR2_SRonly)
factr2_refseq <- fread(factR2_RefSeq)
factr2_ensembl <- fread(factR2_Ensembl)

### Process

# Average log(kdeg)'s
mix_PTC <- factr2_mix %>%
  filter(cds == "yes") %>%
  dplyr::select(transcript_id, is_NMD)

sronly_PTC <- factr2_sronly %>%
  filter(cds == "yes") %>%
  dplyr::select(transcript_id, is_NMD)


refseq_PTC <- factr2_refseq %>%
  filter(cds == "yes") %>%
  dplyr::select(transcript_id, is_NMD)

ensembl_PTC <- factr2_ensembl %>%
  filter(cds == "yes") %>%
  dplyr::select(transcript_id, is_NMD)



### Combine PTC with kinetic parameter estimates

combined_mix <- avg_mix %>%
  dplyr::inner_join(
    mix_PTC,
    by = "transcript_id"
  )


combined_sr <- avg_sr %>%
  dplyr::inner_join(
    sronly_PTC,
    by = "transcript_id"
  )


combined_refseq <- avg_refseq %>%
  dplyr::inner_join(
    refseq_PTC,
    by = "transcript_id"
  )

combined_ensembl <- avg_ensembl %>%
  dplyr::inner_join(
    ensembl_PTC,
    by = "transcript_id"
  )


##### Make density plots


all_analyses <- bind_rows(list(
  combined_mix %>%
    mutate(annotation = "Mix"),
  combined_sr %>%
    mutate(annotation = "SRonly"),
  combined_refseq %>%
    mutate(annotation = "RefSeq"),
  combined_ensembl %>%
    mutate(annotation = "Ensembl")
)
) %>%
  mutate(
    annotation = factor(annotation,
                        levels = c("Ensembl", "RefSeq", "SRonly", "Mix")),
    LFC_kdeg = mean_treatment11j - mean_treatmentDMSO
  )


all_analyses %>%
  ggplot(
    aes(x = mean_treatmentDMSO,
        color = annotation)
  ) +
  geom_density() + 
  theme_classic() + 
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +  xlab("DMSO log(kdeg)") + 
  ylab("Density") +
  coord_cartesian(
    xlim = c(-5, 1)
  )


all_analyses %>%
  ggplot(
    aes(x = mean_treatment11j,
        color = annotation)
  ) +
  geom_density() + 
  theme_classic() + 
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +  xlab("SGM1i log(kdeg)") + 
  ylab("Density") +
  coord_cartesian(
    xlim = c(-5, 1)
  )


all_analyses %>%
  ggplot(
    aes(x = log(2)/exp(mean_treatment11j),
        color = annotation,
        fill = annotation)
  ) +
  geom_histogram(
    alpha = 0.1,
    position = "identity",
    binwidth = 0.2,
    aes(y = after_stat(count / sum(count)))
  ) + 
  theme_classic() + 
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +  scale_fill_viridis_d() + 
  xlab("SGM1i log(kdeg)") + 
  ylab("Density") +
  coord_cartesian(xlim = c(0, 30))



all_analyses %>%
  filter(
    annotation == "Ensembl"
  ) %>%
  ggplot(
    aes(x = LFC_kdeg*log2(exp(1)),
        color = is_NMD,
        fill = is_NMD)
  ) +
  geom_density(
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c("darkgray",
               "darkred")
  ) +
  scale_fill_manual(
    values = c("darkgray",
               "darkred")
  ) +
  theme_classic()  +
  xlab("L2FC(kdeg)") +
  ylab("Density") +
  coord_cartesian(
    xlim = c(-4, 3)
  )



all_analyses %>%
  filter(
    annotation == "RefSeq"
  ) %>%
  ggplot(
    aes(x = LFC_kdeg*log2(exp(1)),
        color = is_NMD,
        fill = is_NMD)
  ) +
  geom_density(
    alpha = 0.2
  ) +
  scale_color_manual(
    values = c("darkgray",
               "darkred")
  ) +
  scale_fill_manual(
    values = c("darkgray",
               "darkred")
  ) +
  theme_classic() +
  xlab("L2FC(kdeg)") +
  ylab("Density")



L2FC_violins <- all_analyses %>%
  ggplot(
    aes(x = LFC_kdeg*log2(exp(1)),
        y = annotation,
        color = annotation,
        fill = annotation)
  ) +
  geom_violin(
    color = 'black',
    linewidth = 0.1
  ) +
  geom_boxplot(
    color = "black",
    fill = "gray90",
    outlier.fill = NA,
    outlier.color = NA,
    notch = TRUE,
    width = 0.25,
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  scale_fill_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  facet_grid(
    rows = vars(is_NMD)
  ) +
  theme_classic() +
  xlab("L2FC(kdeg)") +
  ylab("Annotation") +
  theme(
    legend.position = "none",
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10)
  ) +
  coord_cartesian(
    xlim = c(-3, 1.25)
  )



DMSO_violins <- all_analyses %>%
  ggplot(
    aes(x = mean_treatmentDMSO,
        y = annotation,
        color = annotation,
        fill = annotation)
  ) +
  geom_violin(
    color = 'black',
    linewidth = 0.1
  ) +
  geom_boxplot(
    color = "black",
    fill = "gray90",
    outlier.fill = NA,
    outlier.color = NA,
    notch = TRUE,
    width = 0.25,
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  scale_fill_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  facet_grid(
    rows = vars(is_NMD)
  ) +
  theme_classic() +
  xlab("DMSO log(kdeg)") +
  ylab("Annotation") +
  theme(
    legend.position = "none",
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10)
  ) +
  coord_cartesian(
    xlim = c(-5, 1)
  )



SMG1i_violins <- all_analyses %>%
  ggplot(
    aes(x = mean_treatment11j,
        y = annotation,
        color = annotation,
        fill = annotation)
  ) +
  geom_violin(
    color = 'black',
    linewidth = 0.1
  ) +
  geom_boxplot(
    color = "black",
    fill = "gray90",
    outlier.fill = NA,
    outlier.color = NA,
    notch = TRUE,
    width = 0.25,
    linewidth = 0.25
  ) +
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  scale_fill_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  facet_grid(
    rows = vars(is_NMD)
  ) +
  theme_classic() +
  xlab("SGM1i log(kdeg)") +
  ylab("Annotation") +
  theme(
    legend.position = "none",
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10)
  ) +
  coord_cartesian(
    xlim = c(-5, 1)
  )



setwd(figure_savedir)
ggsave(
  filename = "L2FC_violins.pdf",
  plot = L2FC_violins,
  width = 3,
  height = 3
)
ggsave(
  filename = "DMSO_violins.pdf",
  plot = DMSO_violins,
  width = 3,
  height = 3
)
ggsave(
  filename = "SMG1i_violins.pdf",
  plot = SMG1i_violins,
  width = 3,
  height = 3
)



# Isoform range bar plot -------------------------------------------------------
### Actually decided to plot log(kdeg) instead

assess_range <- function(avgs){
  
  ranges <- avgs %>%
    dplyr::group_by(XF) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::summarise(
      DMSO_range = max(mean_treatmentDMSO) - min(mean_treatmentDMSO),
      DMSO_range_se = sqrt(sum(sd_treatmentDMSO_posterior^2)),
      SMG1i_range = max(mean_treatment11j) - min(mean_treatment11j),
      SMG1i_range_se = sqrt(sum(sd_treatment11j_posterior^2)),
      DMSO_avg_log_readcnt = mean(coverage_treatmentDMSO),
      SGM1i_avg_log_readcnt = mean(coverage_treatment11j)
    ) %>%
    dplyr::mutate(
      range_difference = SMG1i_range - DMSO_range,
      range_difference_se = sqrt(DMSO_range_se^2 + SMG1i_range_se^2)
    )  %>%
    dplyr::mutate(
      DMSO_range_pval = 2*pnorm(-abs(DMSO_range/DMSO_range_se)),
      SMG1i_range_pval = 2*pnorm(-abs(SMG1i_range/SMG1i_range_se)),
      range_difference_pval = 2*pnorm(-abs(range_difference/range_difference_se))
    ) %>%
    dplyr::mutate(
      DMSO_range_padj = p.adjust(DMSO_range_pval, method = "BH"),
      SMG1i_range_padj = p.adjust(SMG1i_range_pval, method = "BH"),
      range_difference_padj = p.adjust(range_difference_pval, method = "BH")
    )
  
  
  return(ranges)
  
  
}


setwd("C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/Annotation_gamut_analyses")
avg_mix <- readRDS(ezbdo_Mix) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_sr <- readRDS(ezbdo_SRonly) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_refseq <- readRDS(ezbdo_RefSeq) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)

avg_ensembl <- readRDS(ezbdo_Ensembl) %>%
  EZget(type = "averages",
        features = "transcript_id",
        exactMatch = FALSE)


avg_ensembl %>%
  filter(XF == "ENSG00000007392") %>%
  dplyr::mutate(
    half_life = log(2) / exp(mean_treatment11j)
  ) %>%
  dplyr::select(
    transcript_id, half_life
  )

avg_sr %>%
  filter(XF == "MSTRG.162") %>%
  dplyr::mutate(
    half_life = log(2) / exp(mean_treatment11j)
  ) %>%
  dplyr::select(
    transcript_id, half_life
  )

avg_refseq %>%
  filter(XF == "LUC7L") %>%
  dplyr::mutate(
    half_life = log(2) / exp(mean_treatment11j)
  ) %>%
  dplyr::select(
    transcript_id, half_life
  )

avg_mix %>%
  filter(XF == "MSTRG.4369") %>%
  dplyr::mutate(
    half_life = log(2) / exp(mean_treatment11j)
  ) %>%
  dplyr::select(
    transcript_id, half_life
  )



range_mix <- assess_range(avg_mix)
range_sr <- assess_range(avg_sr)
range_refseq <- assess_range(avg_refseq)
range_ensembl <- assess_range(avg_ensembl)


range_df <- tibble(
  SMG1i_ranges = c(
    avg_mix$mean_treatment11j,
    avg_sr$mean_treatment11j,
    avg_refseq$mean_treatment11j,
    avg_ensembl$mean_treatment11j
  ),
  assembly = factor(rep(
    c("Mix", "SRonly", "RefSeq", "Ensembl"),
    times = c(
      nrow(avg_mix),
      nrow(avg_sr),
      nrow(avg_refseq),
      nrow(avg_ensembl)
    )
  ),
  levels = c("Ensembl",
             "RefSeq",
             "SRonly",
             "Mix"))
)


IDR_dens <- range_df %>%
  ggplot(
    aes(x = SMG1i_ranges,
        color = assembly)
  ) + 
  geom_density(linewidth = 0.35) + 
  theme_classic() +
  scale_color_manual(
    values = c("#004488", "#6699CC", "#44AA99", "#DDCC77")
  ) +
  xlab("log(kdeg) after SMG1i") +
  ylab("Density") +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) +
  coord_cartesian(
    xlim = c(-5, 1)
  )

IDR_dens


setwd(figure_savedir)
ggsave(
  filename = "lkdeg_density.pdf",
  plot = IDR_dens,
  width = 2.25, 
  height = 2
)



# Stability modeling -----------------------------------------------------------

##### FUNCTIONS #####

fit_lasso <- function(X, y, ID = 1){
  
  
  cv_RNAdeg <- cv.glmnet(X, 
                         y, alpha = 1)
  
  best_lambda  <- cv_RNAdeg$lambda.min
  best_lambda 
  # 0.00287
  
  lasso_RNAdeg  <- glmnet(X,
                          y,
                          alpha = 1,
                          lambda = best_lambda)
  
  
  
  coefs_lasso_RNAdeg <- coef(lasso_RNAdeg)
  
  vars <- coefs_lasso_RNAdeg@Dimnames[[1]][coefs_lasso_RNAdeg@i + 1]
  vals <- coefs_lasso_RNAdeg@x
  
  return(
    tibble(value = vals,
           variable = factor(vars, levels = vars[order(-abs(vals))]),
           ID = ID)
  )
}


##### LOAD DATA #####


kdeg_features <- read_csv(feature_table_path) %>%
  mutate(
    DRACH_count = CDS_DRACH_count + threeputr_DRACH_count + fiveputr_DRACH_count,
    miRNAseed_count = fiveputr_miRNAseed_count + threeputr_miRNAseed_count + CDS_miRNAseed_count
  ) %>%
  dplyr::select(-CDS_DRACH_count, -threeputr_DRACH_count, -fiveputr_DRACH_count,
                -CDS_miRNAseed_count, -threeputr_miRNAseed_count, -fiveputr_miRNAseed_count)

ftokeep <- colnames(kdeg_features)
ftokeep <- ftokeep[!(ftokeep %in% c("mean_treatmentDMSO", "mean_treatment11j", "transcript_id",
                                    "minus1_AA", "minus2_AA"))]

##### Engineer features ######

kdeg_features <- engineer_features(kdeg_features ,
                                   features = ftokeep,
                                   outcome = "mean_treatmentDMSO",
                                   additional_info_to_keep = c("transcript_id"))



##### LASSO #####

X_RNAdeg <- kdeg_features %>%
  dplyr::select(!!ftokeep) %>%
  data.matrix()

y_RNAdeg <- kdeg_features %>%
  dplyr::select(mean_treatmentDMSO) %>%
  data.matrix()


lasso_results <- vector(
  mode = "list",
  length = 100
)

for(i in 1:100){
  
  rows <- sample(1:nrow(X_RNAdeg), size = nrow(X_RNAdeg), replace = TRUE)
  
  X_sub <- X_RNAdeg[rows,]
  y_sub <- y_RNAdeg[rows,]
  
  lasso_results[[i]] <- fit_lasso(X_sub, y_sub, ID = i)
  
}

bind_rows(lasso_results)

LASSO_coefplot <- bind_rows(lasso_results) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    variable = fct_recode(variable,
                          "# of exons" = "total_exons",
                          "CDS length" = "cds_length",
                          "# of downstream EJs" = "num_of_downEJs" ,
                          "# of DRACH motifs" = "DRACH_count",
                          "# of m6A sites" = "m6A_cnts_CLIPseq",
                          "5' UTR length" = "utr_5_length",
                          "3' UTR length" = "threepUTR_length",
                          "# of exons" = "total_exons",
                          "Codon optimality" = "log_CAI",
                          "miRNA seed count" = "miRNAseed_count"
    )
  ) %>%
  group_by(variable) %>%
  dplyr::summarise(
    value_avg = mean(value),
    value_sd = sd(value)
  ) %>%
  ggplot(aes(x = variable, y = value_avg)) + 
  geom_point(size = 0.6) +
  geom_errorbar(
    aes(ymin = value_avg - value_sd*3,
        ymax = value_avg + value_sd*3),
    linewidth = 0.2,
    width = 0.8
  ) +
  theme_classic() + 
  xlab("Feature") + 
  ylab("LASSO coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(
    yintercept = 0,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + 
  coord_cartesian(
    ylim = c(-0.65, 0.65)
  )



LASSO_barplot <- bind_rows(lasso_results) %>%
  filter(variable != "(Intercept)") %>%
  mutate(
    variable = fct_recode(variable,
                          "# of exons" = "total_exons",
                          "CDS length" = "cds_length",
                          "# of downstream EJs" = "num_of_downEJs" ,
                          "# of DRACH motifs" = "DRACH_count",
                          "# of m6A sites" = "m6A_cnts_CLIPseq",
                          "5' UTR length" = "utr_5_length",
                          "3' UTR length" = "threepUTR_length",
                          "# of exons" = "total_exons",
                          "Codon optimality" = "log_CAI",
                          "miRNA seed count" = "miRNAseed_count"
    )
  ) %>%
  ggplot(aes(x = variable, y = value)) + 
  geom_point(size = 0.6) +
  geom_boxplot(outlier.colour = NA,
               outlier.fill = NA) +
  theme_classic() + 
  xlab("Feature") + 
  ylab("LASSO coefficient") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_hline(
    yintercept = 0,
    color = 'darkred',
    linewidth = 0.5,
    linetype = "dotted"
  ) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + 
  coord_cartesian(
    ylim = c(-0.65, 0.65)
  )

LASSO_barplot


setwd(figure_savedir)
ggsave("LASSO_coefficients_mix_trimmed.pdf",
       plot = LASSO_coefplot,
       width = 3.5,
       height = 3)



