### PURPOSE OF THIS SCRIPT
## Reproduce Figure S4


# Load dependencies ------------------------------------------------------------

library(data.table)
library(dplyr)
library(ggplot2)
library(readr)
library(devtools)
library(MASS)
library(tidyr)
devtools::load_all("C:/Users/isaac/Documents/Simon_Lab/EZbakR/")
library(rtracklayer)
library(GenomicFeatures)
library(tximport)
library(stringr)
library(glmnet)
library(forcats)


##### User input #####


# kdeg estimates
EZbakR_estimates <- "C:/Users/isaac/Box/TimeLapse/Annotation_gamut/DataTables/RNAdeg_data_model_features.csv"
ActD_estimates <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/ML_features/ActD_kdeg_ests.csv"


# Isoform features for LASSO regression
feature_table_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Data/ML_features/NMD_tables/Filtered_kdeg_feature_table.csv"

# Path to function_library.R
function_library_path <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Kinetics/Scripts/ML/function_library.R"

# Directory to save figures in
figure_savedir <- "C:/Users/isaac/Documents/Simon_Lab/Isoform_Paper/Figures/Figure_3/"


##### FUNCTIONS USED THROUGHOUT #####

source(function_library_path)

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


# Transcription inhibition comparison ------------------------------------------

##### Compare with EZbakR estimates #####

ez_ests <- fread(EZbakR_estimates)
kd_ests <- fread(ActD_estimates)


compare_ests <- ez_ests %>%
  as_tibble() %>%
  dplyr::inner_join(kd_ests,
                    by = "transcript_id") %>%
  dplyr::mutate(log_threepUTR_length = log(threepUTR_length + 1),
                avg_log_DMSO_TPM = log(avg_TPM_DMSO),
                log_fivepUTR_length = log(str_length(fiveputr_seq) + 1),
                log_CDS_length = log(str_length(CDS_seq) + 1))

gtxi_vs_EZ <- compare_replicates(compare_ests,
                                 "log_kdeg_DMSO",
                                 "log_kdeg",
                                 add_regression = TRUE,
                                 add_yeqx = TRUE,
                                 lw = 0.5) +
  theme(legend.position="none") +
  theme(
    axis.text=element_text(size=10), #change font size of axis text
    axis.title=element_text(size=12), #change font size of axis titles
    legend.text=element_text(size=10), #change font size of legend text
    legend.title=element_text(size=12))

setwd(figure_savedir)
ggsave(
  filename = "ActD_vs_EZbakR_kdeg.pdf",
  plot = gtxi_vs_EZ,
  width = 3.25,
  height = 2.75
)


# Transcription inhibition kdeg feature modeling -------------------------------

kd_ests <- fread(ActD_estimates)


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

kdeg_features <- kdeg_features %>%
  inner_join(kd_ests %>%
               dplyr::rename(log_kdeg_ActD = log_kdeg) %>%
               dplyr::select(transcript_id, log_kdeg_ActD),
             by = "transcript_id")

kdeg_features <- engineer_features(kdeg_features ,
                                   features = ftokeep,
                                   outcome = "log_kdeg_ActD",
                                   additional_info_to_keep = c("transcript_id"))



##### LASSO #####

X_RNAdeg <- kdeg_features %>%
  dplyr::select(!!ftokeep) %>%
  data.matrix()

y_RNAdeg <- kdeg_features %>%
  dplyr::select(log_kdeg_ActD) %>%
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

glasso_actd <- bind_rows(lasso_results)  %>%
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


setwd(figure_savedir)
ggsave("LASSO_coefficients_ActD.pdf",
       plot = glasso_actd,
       width = 3.5,
       height = 3)


