#### Organized Olink HT QC Pipeline #####
## Packages Used
library(tidyverse)
library(arrow)
library(tools)


Olink_Intensity_Norm <- function(data_olink) {
  # calculate ExtNPX based on counts in the Olink parquet file
  data_olink_cal <-
    data_olink %>%
    dplyr::group_by(SampleID, PlateID, Block) %>%
    mutate(LogProtExp = ExtNPX + log2(1e5)) %>% 
    ungroup()
  
  
  return(data_olink_cal) # return the output
}

pre_comb <- function(data_olink){
  
  ht_nc_vals <- data_olink %>% 
    filter(SampleType == "NEGATIVE_CONTROL") %>% 
    group_by(Assay, OlinkID) %>% 
    summarise(median_nc = median(LogProtExp),
              iqr_nc = quantile(LogProtExp, 0.75), .groups = 'drop')
  
  # Calculating Plate Control coefficient of variance
  ht_pc_vals <- data_olink %>% 
    filter(SampleType == "PLATE_CONTROL") %>% 
    group_by(Assay, OlinkID) %>% 
    summarise(pc_cv = 100*sd(LogProtExp)/mean(LogProtExp), .groups = 'drop') %>% 
    mutate(high_var_assay = case_when(pc_cv > 20 ~ "High Variance",
                                      T ~ "Pass")) %>% 
    select(-pc_cv)
  
  # This is the "sample level" qc, calculates ith sample in jth assay that needs to be replaced with zero or LLOQ
  # also labels those values in a new column - sample_level_qc
  ht_scaled_npx_sample <- data_olink %>% 
    filter(SampleType == "SAMPLE") %>% 
    left_join(., ht_nc_vals, by = c("Assay", "OlinkID")) %>% 
    left_join(., ht_pc_vals, by = c("Assay", "OlinkID")) %>% 
    mutate(sample_level_qc = case_when(LogProtExp < median_nc ~ "Below LLOD",
                                       LogProtExp < iqr_nc ~ "Below LLOQ",
                                       T ~ "Pass")) %>% 
    mutate(LogProtExp = case_when(LogProtExp < median_nc ~ 0,
                                  LogProtExp < iqr_nc ~ iqr_nc,
                                  T ~ LogProtExp)) 
  

  return(ht_scaled_npx_sample)
  
  
}

post_comb <- function(data_olink){
  # specifying number of samples present in the total combined dataset
   n_samples <- data_olink %>% 
    group_by(SampleID) %>%
    summarise() %>% 
    nrow()
  
  
  # Assay level QC - if 50% of samples are below LLOQ, labeled as semi-continuous
  # if 75% of samples are below LLOD, labeled as categorical
  # Test to adjust how I calculate categorical, semi-continuous, or continuous
  ht_scaled_npx_assay <- data_olink %>% 
    group_by(Assay, OlinkID, sample_level_qc) %>%
    summarise(percentage = 100*n()/n_samples, .groups = 'drop') %>%
    filter(sample_level_qc %in% c("Below LLOD", "Below LLOQ")) %>% 
    pivot_wider(., names_from = sample_level_qc, values_from = percentage) %>% 
    mutate(`Below LLOD` = case_when(is.na(`Below LLOD`) == T ~ 0,
                                    T ~ `Below LLOD`),
           `Below LLOQ` = case_when(is.na(`Below LLOQ`) == T ~ 0,
                                    T ~ `Below LLOQ`)) %>% 
    mutate(`Below LLOQ` = `Below LLOQ` + `Below LLOD`) %>% 
    mutate(assay_level_qc = case_when(`Below LLOD` > 75 ~ "Categorical",
                                      `Below LLOQ` > 50 ~ "Semi-Continuous",
                                      T ~ "Continuous")) %>% 
    select(-c(`Below LLOD`, `Below LLOQ`))
  
  # For samples that didn't include the assay_level_qc call, we can determine that there are 0 samples in that assay
  # that were below the LLOD, and that assay is good quality and should be called continuous
  # This data only contains the patient samples themselves
  ht_scaled_npx_qc <- data_olink %>% 
    left_join(., ht_scaled_npx_assay, by = c("OlinkID", "Assay")) %>% 
    mutate(assay_level_qc = case_when(is.na(assay_level_qc) == T ~ "Continuous",
                                      T ~ assay_level_qc))  %>% 
    select(-c(median_nc, iqr_nc))
}


### HERE IS WHERE THE ISSUE IS ####
# Input single parquet file or list of parquet files that need to be qc'd
olink_qc_run <- function(file){ 
  
  if(length(file) > 1){
    message("Multiple Olink HT runs detected")
    colname_list <- lapply(file, function(y) arrow::open_dataset(y)$schema$names)
    
    
    bool <- length(unique(lapply(colname_list, function(x) sort(toupper(x))))) == 1
    
    if(bool == FALSE) warning("File columns do not match. Check parquet files for differences")
    
    single_dat <- lapply(file, function(y) arrow::read_parquet(y))
    
    filter_dat <- lapply(single_dat, function(y) y %>% filter(is.na(ExtNPX) == FALSE))
    
    scaled_dat <- lapply(filter_dat, function(x) Olink_Intensity_Norm(data_olink = x))

    # Calculating LLOD and LLOQ based on negative control values
    
    pre_comb_dat <- lapply(scaled_dat, function(x) pre_comb(data_olink = x))
    
    lapply(pre_comb_dat, function(x) arrow::write_parquet(x, paste0("Run_Specific_QC", xname, ".parquet")))
    
    return(pre_comb_dat)
    
  } else if (length(file) == 1){
    message("1 Olink HT run detected")
    dat <- read_parquet(file)
    
    filter_dat <- dat %>% filter(is.na(ExtNPX) == FALSE)
    
    scaled_dat <- Olink_Intensity_Norm(data_olink = filter_dat)
    
    pre_comb_dat <- pre_comb(data_olink = scaled_dat)
    
    return(pre_comb_dat)
    
      
  } else (stop('length of file is less than 1'))
  }
  


olink_qc_project <- function(file){
  
  comb_dat <- bind_rows(pre_comb_dat)
  
  post_comb_dat <- post_comb(comb_dat)
  
  return(post_comb_dat)
}





dat <- olink_qc_run(file = c("\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_1-2_NPX_2024-02-07.parquet", 
                         "\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_3-4_NPX_2024-02-07.parquet", 
                         "\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_5-6_NPX_2024-02-12.parquet"))








