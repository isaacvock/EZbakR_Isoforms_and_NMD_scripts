### PURPOSE OF THIS SCRIPT
## Reproduce Supplemental Figure 1


# Load dependencies ------------------------------------------------------------

library(dplyr)
library(data.table)
library(ggplot2)
library(MASS)
library(tidyr)
library(EZbakR)
library(readr)



##### User input #####

# PE data
NRsim_EZbdo_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/EZbakRFits/NRsim_PE_validation_nopremRNA_refseq/NRsim_refseq_nocB.rds"
ground_truth_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Simulations/PE_validation_nopreRNA_refseq/generate_transcript_kinetics/kinetics.csv"

# Different pnew data
NRsim_pnewrange_EZbdo_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/EZbakRFits/NRsim_pnew_range/NRsim_pnew_range_nocB.rds"
ground_truth_pnewrange_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Simulations/PEval_pnew_range/generate_transcript_kinetics/kinetics.csv"

# SE data
NRsim_SE_EZbdo_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/EZbakRFits/NRsim_SE_validation/NRsim_SE_validation_nocB.rds"
ground_truth_SE_path <- "C:/Users/isaac/Yale University/Simon Lab – RNA - Documents/IWV/Hogg_lab/Simulations/SE_validation_refseq/generate_transcript_kinetics/kinetics.csv"

# Directory to save final figures to
figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Paper/Figures/Figure_1/"

##### Functions used throughout #####

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


assess_accuracy <- function(ests,
                            sample_name,
                            min_reads = 100,
                            density_res = 100,
                            lw = 0.8,
                            ps = 0.5,
                            xvar = "fn",
                            yvar = "fraction_highTC",
                            lt = "dotted",
                            color_var = "density"){
  
  ests <- ests %>%
    dplyr::filter(sample == sample_name &
                    n > min_reads) %>%
    ungroup()
  
  if(color_var == "density" & !("density" %in% colnames(ests))){
    
    ests <- ests %>%
      dplyr::mutate(
        density = get_density(
          x = !!dplyr::sym(xvar),
          y = !!dplyr::sym(yvar),
          n = density_res
        )
      )
    
  }
  
  
  gacc <- ests %>%
    ggplot(aes(x = !!dplyr::sym(xvar),
               y = !!dplyr::sym(yvar),
               color = !!dplyr::sym(color_var))) + 
    geom_point(size = ps) +
    theme_classic() + 
    scale_color_viridis_c() + 
    geom_abline(slope = 1,
                intercept = 0,
                color = 'darkred',
                linewidth = lw,
                linetype = "dotted")
  
  return(gacc)
  
}




# PE read length ---------------------------------------------------------------

ezbdo <- readRDS(NRsim_EZbdo_path)
gt <- read_csv(ground_truth_path)

fn_ests <- ezbdo$fractions$isoforms %>%
  inner_join(gt, by = "transcript_id")

# 25 million reads, PE-150
g150 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                          filter(se_logit_fraction_highTC < 0.5), 
                        sample_name = "sim_5",
                        lw = 0.3,
                        ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")

# PE-100
g100 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                          filter(se_logit_fraction_highTC < 0.5), 
                        sample_name = "sim_2",
                        lw = 0.3,
                        ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


# PE-50
g50 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                         filter(se_logit_fraction_highTC < 0.5), 
                       sample_name = "sim_7",
                       lw = 0.3,
                       ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


### Save figures
setwd(figure_savedir)
ggsave(filename = "NRsim_validation_PE150.pdf",
       plot = g150,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_PE100.pdf",
       plot = g100,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_PE50.pdf",
       plot = g50,
       width = 2,
       height = 1.67)



# Pnew -------------------------------------------------------------------------


ezbdo <- readRDS(NRsim_pnewrange_EZbdo_path)

gt <- read_csv(ground_truth_pnewrange_path)

fn_ests <- ezbdo$fractions$isoforms %>%
  inner_join(gt, by = "transcript_id")

# 10 million reads, PE-100
g10 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                         filter(se_logit_fraction_highTC < 0.5), 
                       sample_name = "sim_1",
                       lw = 0.3,
                       ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


g05 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                         filter(se_logit_fraction_highTC < 0.5), 
                       sample_name = "sim_2",
                       lw = 0.3,
                       ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


g01 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                         filter(se_logit_fraction_highTC < 0.5), 
                       sample_name = "sim_4",
                       lw = 0.3,
                       ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


### Save figures
setwd(figure_savedir)
ggsave(filename = "NRsim_validation_pnew10perc.pdf",
       plot = g10,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_pnew5perc.pdf",
       plot = g05,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_pnew1perc.pdf",
       plot = g01,
       width = 2,
       height = 1.67)




# SE read length ---------------------------------------------------------------



ezbdo <- readRDS(NRsim_SE_EZbdo_path)

gt <- read_csv(ground_truth_SE_path)

fn_ests <- ezbdo$fractions$isoforms %>%
  inner_join(gt, by = "transcript_id")

# SE-50
gSE50 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                           filter(se_logit_fraction_highTC < 0.5), 
                         sample_name = "sim_5",
                         lw = 0.3,
                         ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")

# SE-100
gSE100 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                            filter(se_logit_fraction_highTC < 0.5), 
                          sample_name = "sim_7",
                          lw = 0.3,
                          ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


# SE-150
gSE150 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                            filter(se_logit_fraction_highTC < 0.5), 
                          sample_name = "sim_8",
                          lw = 0.3,
                          ps = 0.2) +
  theme(
    axis.text=element_text(size=8), #change font size of axis text
    axis.title=element_text(size=10), #change font size of axis titles
    legend.text=element_text(size=8), #change font size of legend text
    legend.title=element_text(size=10)) + #change font size of legend title
  theme(legend.position = "none") +
  xlab("True fn") + 
  ylab("Est. fn")


### Save figures
setwd(figure_savedir)
ggsave(filename = "NRsim_validation_SE50.pdf",
       plot = gSE50,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_SE100.pdf",
       plot = gSE100,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_SE150.pdf",
       plot = gSE150,
       width = 2,
       height = 1.67)
