### PURPOSE OF THIS SCRIPT
## Reproduce figure 1B

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


# Figure 1B --------------------------------------------------------------------

ezbdo <- readRDS(NRsim_EZbdo_path)
gt <- read_csv(ground_truth_path)

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


g50 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
                         filter(se_logit_fraction_highTC < 0.5), 
                       sample_name = "sim_3",
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


g100 <- assess_accuracy(fn_ests %>% filter(!grepl("_", location)) %>%
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
ggsave(filename = "NRsim_validation_10mil_reads.pdf",
       plot = g10,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_50mil_reads.pdf",
       plot = g50,
       width = 2,
       height = 1.67)
ggsave(filename = "NRsim_validation_100mil_reads.pdf",
       plot = g100,
       width = 2,
       height = 1.67)
