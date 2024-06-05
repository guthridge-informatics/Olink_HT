### Level 1 and Level 2 QC for Olink HT Explore Data ###
require(dplyr)
require(arrow)
require(tidyr)
require(stringr)
require(umap)
require(ggplot2)
require(gridExtra)
require(dbscan)
require(ggrepel)
require(readxl)


##### Subfunctions ###################################################################################
## CHANGE MADE - REMOVED 8:12 DESIGNATION
# Misc function to read in multiple files into a list and name them based on the file name without file extension
Multifile_read <- function(directory, file_extension) {
  files_run <- list.files(directory) %>% # list all files within the given directory
    str_extract(string = ., 
                pattern = paste0(".*", file_extension)) %>%
    na.omit %>% # only extract files with the provided file extension
    as.character() 
     # omit empty rows
  if (file_extension == "parquet") {
    data <- lapply(X = files_run, FUN = read_parquet)
  } else if (file_extension == "xlsx"){
    data <-
      lapply(
        X = files_run,
        FUN = function(x)
          read_excel(path = x, sheet = 1) %>% select(c(4, 6))
      )
  }
  names(data) = str_extract_all(string = files_run, pattern = paste0("[^\\.", file_extension, "]+"), simplify = TRUE)
  return(data)
}

# Misc function to write multiple files into a list and name them based on the file name without file extension
Multifile_write <- function(data, file_extension) {
  names = names(data) # get the names of the list of the data
  if (file_extension == "parquet") {
    lapply(
      X = names,
      FUN = function(x)
        arrow::write_parquet(data[[x]], sink = paste0("SDP/Level_1/", x, ".parquet"))
    )
  } else if (file_extension == "csv"){
    lapply(
      X = names,
      FUN = function(x)
        write.csv(data[[x]], file = paste0("SDP/Level_1/", x, ".csv"))
    )
  }
}

####################################################################################################################

## Pre Level 1
# This reader function scan the directory for parquet files and read them into a list format
Olink_Reader <- function(directory){
  # Defining the directory all Olink data is stored in
  setwd(directory) # Need to set the working directory to such or else subsequent code will not work
  
  # make sub-directory for SDP hiearchical file structure
  dirs = c(
    root = paste0(directory,"/SDP"), # root directory for SDP
    lvl1 = paste0(directory, "/SDP/Level_1"), # level 1 directory for SDP
    lvl2 = paste0(directory, "/SDP/Level_2"), # level 2 directory for SDP
    code = paste0(directory, "/SDP/code")
  )
  lapply(X = dirs, function(i) dir.create(path = i))
  
  # Extract all .parquet files within the directory and read into a list of data frames
  data <- Multifile_read(directory, "parquet")
  # Extract all .csv files within the directory and read into a list of data frames that corresponding to the sample manifest of each OLINK Explore HT run
  manifest <- Multifile_read(directory, "xlsx")
  
  # Return parquet files
  return(list(data = data, manifest = manifest))
}

## Level 1 QC
Olink_lvl1 <- function(olink_files, proj_names){
  names = names(olink_files$data) # get the names of all the files aka list names
  # store in the standard control names in the SampleID column
  ctrls = c("PC1", "PC2", "PC3", "PC4", "PC5", "NC1", "NC2", "SC1", "SC2", "SC3")
  manifest_filtered <-
    lapply(
      olink_files$manifest,
      FUN = function(x)
        dplyr::filter(x, Project %in% proj_names)
    ) # filter the manifest file by the project names first
  Multifile_write(data = manifest_filtered, file_extension = "csv") # write the SMDs
  
  # filter the current raw data by only selecting the sample only pertaining to the project names
  data_filtered <-
    lapply(
      names,
      FUN = function(x)
        dplyr::filter(
          olink_files$data[[x]],
          SampleID %in% c(ctrls, manifest_filtered[[x]]$sample_id)
        ) %>%
        left_join(olink_files$manifest[[x]], by = dplyr::join_by("SampleID" == "sample_id"))
    ) # filter only the sample IDs belongs to the projects
  names(data_filtered) = names
  Multifile_write(data = data_filtered, file_extension = "parquet") # write the Level 1 parquet
  return(list(data = data_filtered, manifest = manifest_filtered))
}

## Batch Correction Process (For Multi-Plate Studies only)
Olink_QC <- function(data) {
  table_assay <- 
    data %>%
    dplyr::filter(AssayType == "assay") %>% # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, OlinkID, AssayQC) %>% # group by PlateID and OlinkID to get the split among the Assay QC
    tally %>% # tabulate the counts
    ungroup %>% # this is to allow the ungrouping of the tibble and return to original tibble without grouping
    dplyr::group_by(PlateID, AssayQC) %>% # grouping to allow tabulation of assays
    tally
  
  table_sample <- 
    data %>%
    dplyr::filter(SampleType == "SAMPLE") %>% # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, SampleID, SampleQC) %>% # group by PlateID and OlinkID to get the split among the Assay QC
    tally %>% # tabulate the counts
    mutate(Frequency = prop.table(n)) # calculate the frequencies within each sample
  
  return(list(AssayOlinkQC = table_assay, SampleOlinkQC = table_sample))
}

# Calculate the median and variance of each plate controls on a 384-well plate for each run
ctrl_ref <- function(data) {
  ref <-
    data %>%
    dplyr::filter(SampleType == "PLATE_CONTROL" & AssayType == "assay") %>% # filter down to just the plate controls and the assay
    dplyr::group_by(OlinkID) %>% # grouped it by just the OlinkID
    summarise(Median = median(na.omit(ExtNPX)), # calculate the median
              Variance = var(na.omit(ExtNPX))) # calculate variance
  return(ref)
}

# This function batch correct the 
median_correction <- function(data, Meds){
  # calculate the medians for each 384-well plate 
  ref_med <- 
    Reduce(function(x, y) left_join(x, y, by = "OlinkID"), Meds) %>% # left join all the run median and variance per each run
    dplyr::select(contains("Median")) %>%  # filter only the median column names
    apply(1, mean) %>% # calculate the mean of all medians to scale it to it
    data.frame(OlinkID = Meds[[1]]$OlinkID, ReferenceMedian = .) # make a resulting data frame that contains the OlinkID and the reference median
  
  Meds_correction <- 
    lapply(Meds, FUN = function(x) left_join(x, ref_med)) %>%
    lapply(FUN = function(x) mutate(x, Correction = Median - ReferenceMedian))
  
  data_correction <-
    mapply(function(x, y) {left_join(x, y, by = "OlinkID")}, 
           x = data, 
           y = Meds_correction, 
           SIMPLIFY = FALSE) %>%
    lapply(FUN = function(x) mutate(x, ExtNPX_Corrected = ExtNPX - Correction))
  return(data_correction)  
}

# This is the main routing function for the batch correction and normalization
batch_correction <- function(data, method = "median"){
  if (method == "median"){
    Meds <- lapply(data, FUN = ctrl_ref)
    data_corrected <- median_correction(data, Meds)
  }
  return(data_corrected)
}

## Level 2 QC Part 1
Olink_lvl2_prep <- function(data){
  data <- 
    data %>%
    mutate(LogProtExp = ExtNPX_Corrected + log2(1e5)) %>% # transform into LogProExp
    mutate(LogProtExp_Raw = ExtNPX_Corrected + log2(1e5)) # transform into LogProExp
  
  ht_nc_vals <- data %>% 
    filter(SampleType == "NEGATIVE_CONTROL") %>% 
    group_by(Assay, OlinkID) %>% 
    summarise(median_nc = median(na.omit(LogProtExp)),
              iqr_nc = as.numeric(quantile(na.omit(LogProtExp), 0.75)), .groups = 'drop') 
  
  # Calculating Plate Control coefficient of variance
  ht_pc_vals <- data %>% 
    filter(SampleType == "PLATE_CONTROL") %>% 
    group_by(Assay, OlinkID) %>% 
    summarise(pc_cv = 100*sd(LogProtExp)/mean(LogProtExp), .groups = 'drop') %>% 
    mutate(high_var_assay = case_when(pc_cv > 20 ~ "High Variance",
                                      T ~ "Pass")) %>% 
    select(-pc_cv)
  
  # This is the "sample level" qc, calculates ith sample in jth assay that needs to be replaced with zero or LLOQ
  # also labels those values in a new column - sample_level_qc
  ht_scaled_npx_sample <- data %>% 
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

## Level 2 QC Part 2
Olink_lvl2 <- function(data){
  # concatenate all the data together
  
  data <- Reduce(function(x, y) rbind(x, y), data)
  
  # specifying number of samples present in the total combined dataset
  n_samples <- data %>% 
    group_by(SampleID) %>%
    summarise() %>% 
    nrow()
  
  
  # Assay level QC - if 50% of samples are below LLOQ, labeled as semi-continuous
  # if 75% of samples are below LLOD, labeled as categorical
  # Test to adjust how I calculate categorical, semi-continuous, or continuous
  ht_scaled_npx_assay <- data %>% 
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
  
  data <- left_join(data, dplyr::select(ht_scaled_npx_assay, OlinkID, assay_level_qc), by = "OlinkID")
  
  

  arrow::write_parquet(data, sink = paste0("SDP/Level_2/Level 2 SDP.parquet")) # write the parquet files of the level 2
  return(data_lvl2 = data)
}

##### BATCH CORRECTION AND GRAPHING FUNCTIONS 

# This function tallies the number of failed and passed assay and samples per each 96-well plate
Olink_QC <- function(data) {
  table_assay <- 
    data %>%
    dplyr::filter(AssayType == "assay") %>% # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, OlinkID, AssayQC) %>% # group by PlateID and OlinkID to get the split among the Assay QC
    tally %>% # tabulate the counts
    ungroup %>% # this is to allow the ungrouping of the tibble and return to orignal tibble without grouping
    dplyr::group_by(PlateID, AssayQC) %>% # grouping to allow tabulation of assays
    tally
  
  table_sample <- 
    data %>%
    dplyr::filter(SampleType == "SAMPLE") %>% # filter for only assays instead of all extension, plate controls
    dplyr::group_by(PlateID, SampleID, SampleQC) %>% # group by PlateID and OlinkID to get the split among the Assay QC
    tally %>% # tabulate the counts
    mutate(Frequency = prop.table(n)) # calculate the frequencies within each sample
  
  return(list(AssayOlinkQC = table_assay, SampleOlinkQC = table_sample))
}


ctrl_ref <- function(data) {
  ref <-
    data %>%
    dplyr::filter(SampleType == "PLATE_CONTROL" & AssayType == "assay") %>% # filter down to just the plate controls and the assay
    dplyr::group_by(OlinkID) %>% # grouped it by just the OlinkID
    summarise(Median = median(na.omit(ExtNPX)), # calculate the median
              Variance = var(na.omit(ExtNPX))) # calculate variance
  return(ref)
}

# This function batch correct the 
median_correction <- function(data, Meds){
  # calculate the medians for each 384-well plate 
  ref_med <- 
    Reduce(function(x, y) left_join(x, y, by = "OlinkID"), Meds) %>% # left join all the run median and variance per each run
    dplyr::select(contains("Median")) %>%  # filter only the median column names
    apply(1, mean) %>% # calculate the mean of all medians to scale it to it
    data.frame(OlinkID = Meds[[1]]$OlinkID, ReferenceMedian = .) # make a resulting data frame that contains the OlinkID and the reference median
  
  Meds_correction <- 
    lapply(Meds, FUN = function(x) left_join(x, ref_med)) %>%
    lapply(FUN = function(x) mutate(x, Correction = Median - ReferenceMedian))
  
  data_correction <-
    mapply(function(x, y) {left_join(x, y, by = "OlinkID")}, 
           x = data, 
           y = Meds_correction, 
           SIMPLIFY = FALSE) %>%
    lapply(FUN = function(x) mutate(x, ExtNPX_Corrected = ExtNPX - Correction))
  return(data_correction)  
}

# This is the main routing function for the batch correction and normalization
batch_correction <- function(data, method = "median"){
  if (method == "median"){
    Meds <- lapply(data, FUN = ctrl_ref)
    data_corrected <- median_correction(data, Meds)
  }
  return(data_corrected)
}

# this function automatically graphcs the before and post normalization and check for umap of each run

normalization_check <- function(data_corrected, pt.size = 0.5){
  data_corrected_combined <- 
    data_corrected %>% 
    dplyr::filter(AssayType == "assay" & SampleType == "SAMPLE") %>%
    mutate(ProteinID = paste0(Assay, "_", OlinkID)) %>%
    dplyr::select(SampleID, PlateID, ProteinID, ExtNPX, ExtNPX_Corrected, LogProtExp_Raw) %>%
    pivot_wider(names_from = ProteinID, values_from = c(ExtNPX, ExtNPX_Corrected, LogProtExp_Raw))
  
  plot1 <-
    data_corrected_combined %>%
    na.omit %>%
    dplyr::select(contains("ExtNPX")) %>%
    UMAP_groups(groups = na.omit(data_corrected_combined)$PlateID, n_neighbors = 30, pt.size = pt.size) + ggtitle("Raw ExtNXP")
  
  plot2 <-
    data_corrected_combined %>%
    na.omit %>%
    dplyr::select(contains("ExtNPX_Corrected")) %>%
    UMAP_groups(groups = na.omit(data_corrected_combined)$PlateID, n_neighbors = 30, pt.size = pt.size) + ggtitle("Batch-corrected ExtNXP")
  
  plot3 <-
    data_corrected_combined %>%
    na.omit %>%
    dplyr::select(contains("LogProtExp_Raw")) %>%
    UMAP_groups(groups = na.omit(data_corrected_combined)$PlateID, pt.size = pt.size, n_neighbors = 30, min_dist = 1) + ggtitle("Batch-corrected LogProExp")
  
  return(list(plot1 = plot1, plot2 = plot2, plot3 = plot3))
}




###### UMAP OPTIMAZATION #################
## This module is made to better optimize the UMAP visualization for any omics
require(gridExtra) # needed for the grid plot arrangement
require(pathviewr) # needed for find_curve_elbow for PCA post hoc analysis

# This function uses variance stabilization tranformation to detect variables/analytes with above expected variances
vst_to_pca <- function(data, exclude = ""){
  # data excludes the categorical columns that will not be used for the vst process
  data <- dplyr::select(data, -which(names(data) %in% exclude))
  # calculate the vst matrix of mean and variance
  if (class(data)[1] == "tbl_df") {
    mean <- apply(X = data, MARGIN = 2, FUN = function(x) mean(unlist(x)))
    variance <- apply(X = data, MARGIN = 2, FUN = function(x) var(unlist(x)))
  } else {
    mean <- apply(X = data, MARGIN = 2, mean)
    variance <- apply(X = data, MARGIN = 2, var)
  }
  # create the data_vst data frame that contains mean and variances
  data_vst <- data.frame(Assay = colnames(data), 
                         Mean = mean,
                         Variance = variance)
  
  # plot the scatter plots of the mean and variance of each protein
  fit_lowess <- loess(Variance ~ Mean, data_vst, span = 0.2)
  
  # graph the loess function in the plot of mean vs variance
  plot <-
    ggplot(data = data_vst, aes(x = Mean, y = Variance)) +
    geom_point() +
    geom_smooth(method = loess, formula = y ~ x, lty = 2, method.args = list(span = 0.2)) +
    theme_classic()
  
  # add predicted variance 
  data_vst$PredictedVar <- predict(fit_lowess, data_vst$Mean)
  
  # add difference score above the predicted means
  data_vst$Diff <- data_vst$Variance - data_vst$PredictedVar
  return(list(data_vst = data_vst, plot = plot))
}

# This function automatically find the best PCA dimensional cut-off for the subsequent UMAP
pca_to_umap <- function(data){
  # calculate the principal components
  pca <- prcomp(na.omit(data), scale. = TRUE, center = TRUE)
  # get PCA importance from the principal component analysis
  res_pca <- summary(pca)$importance %>% t() %>% data.frame %>% mutate(PCs = c(1:nrow(.)))
  n_pca <- find_curve_elbow(data_frame = res_pca[, c("PCs", "Proportion.of.Variance")], plot_curve = TRUE)
  pca <- 
    pca$x[, c(1:n_pca)] %>% 
    data.frame()
  return(pca)
}

# Optimization function for n_neighbors by setting the spread = 10, min_dist = 0.1 as default by can be changed as needed
optimize_n_neighbor <- function(data, groups, spread = 10, min_dist = 0.1, min = 5, max = 16, step = 1) {
  plots <- list()
  counter <- 1
  if (is.na(groups)) {
    groups = "black"
  }
  for (n_neighbors in seq(min, max, step)) {
    umap_cyto <- umap(data, spread = spread, min_dist = min_dist, n_neighbors = n_neighbors, random_state = 123)
    dat_umap <- umap_cyto$layout %>% data.frame()
    colnames(dat_umap) <- c("UMAP1", "UMAP2")
    # graph the umap plot
    plot_0 <- ggplot(data = dat_umap, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(aes(color = groups), size = 4) +
      theme_classic() + 
      ggtitle(paste0("n_neighbor = ", n_neighbors)) +
      theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12))
    plots[[counter]] <- plot_0
    counter = counter + 1
  }
  grid.arrange(grobs = plots, ncol = 4)
}

# Optimization function for spread after n_neighbor parameter has been fixed, by default this function uses min_dist of 0.1
optimize_spread <- function(data, groups, n_neighbor, min_dist = 0.1, min = 1, max=15, step = 1){
  # spread
  plots <- list()
  counter <- 1
  if (is.na(groups)) {
    groups = "black"
  }
  for (spread in seq(min, max, step)) {
    umap_cyto <- umap(data, spread = spread, min_dist = min_dist, n_neighbors = n_neighbor, random_state = 123)
    dat_umap <- umap_cyto$layout %>% data.frame()
    colnames(dat_umap) <- c("UMAP1", "UMAP2")
    # graph the umap plot
    plot_0 <- ggplot(data = dat_umap, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(aes(color = groups), size = 4) +
      theme_classic() + 
      ggtitle(paste0("spread = ", spread)) +
      theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12))
    plots[[counter]] <- plot_0
    counter = counter + 1
  }
  grid.arrange(grobs = plots, ncol = 4)
}

# The final optimization function that is to fix the min_dist after n_neighbor and spread parameters have been fixed
optimize_min_dist <- function(data, groups, spread, n_neighbor, min = 0.01, max = 0.5, step = 0.03){
  plots <- list()
  counter <- 1
  if (is.na(groups)) {
    groups = "black"
  }
  for (min_dist in seq(min, max, step)) {
    umap_cyto <- umap(data, spread = spread, min_dist = min_dist, n_neighbors = n_neighbors, random_state = 123)
    dat_umap <- umap_cyto$layout %>% data.frame()
    colnames(dat_umap) <- c("UMAP1", "UMAP2")
    # graph the umap plot
    plot_0 <- ggplot(data = dat_umap, aes(x = UMAP1, y = UMAP2)) + 
      geom_point(aes(color = groups), size = 4) +
      theme_classic() + 
      ggtitle(paste0("min_dst = ", min_dist)) +
      theme(legend.text = element_text(size = 12), axis.text = element_text(size = 12))
    plots[[counter]] <- plot_0
    counter = counter + 1
  }
  grid.arrange(grobs = plots, ncol = 4)
}

# graph the clean reworked umap
UMAP_snn <- function(data, k = NA, eps = 7, minPts = 10, arrow_size = 0.1, pt.size = 4, arrowtip_size = 2, cols = NULL, label = FALSE, label.size = 15, spread = 5, min_dist = 0.25, n_neighbors = NA) {
  if (is.na(n_neighbors)){
    n_neighbors = nrow(data)/10
  }
  umap_cyto <- umap(data, spread = spread, min_dist = min_dist, n_neighbors = n_neighbors, random_state = 123)
  dat_umap <- umap_cyto$layout %>% data.frame()
  colnames(dat_umap) <- c("UMAP1", "UMAP2")
  
  # sNN clustering the UMAP coordinates
  if (is.na(k)) {
    k = nrow(dat_umap)/10
  }
  groups <- as.factor(sNNclust(dat_umap, k = k, eps = eps, minPts = minPts)$cluster)
  
  # graph the umap plot
  p <- ggplot(data = dat_umap, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(color = groups), size = 4) +
    theme_void() 
  # y-range
  yrange = layer_scales(p)$y$range$range
  # x-range
  xrange = layer_scales(p)$x$range$range
  p <- p + 
    guides(fill = guide_legend(title = paste0("k = ", k))) +
    theme(legend.title = element_blank()) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1], yend = yrange[1] + (yrange[2]-yrange[1])*arrow_size, size = 0.8, arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] + 0.2, y = yrange[1], label = "UMAP1"), hjust = 0, vjust = 1, size = 4) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1] + (xrange[2]-xrange[1])*arrow_size, yend = yrange[1], size = 0.8,  arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] - 0.6, y = yrange[1], label = "UMAP2"), angle = 90, hjust = 0, vjust = 1, size = 4) +
    guides(colour = guide_legend(ncol = 1, label.theme = element_text(face = "bold", size = label.size),override.aes = list(size = 6))) +
    ggtitle(label = "")
  return(list(plot = p, data_clust = cbind(dat_umap, Cluster = groups)))
}

# graph the clean reworked umap
UMAP_groups <- function(data, groups, eps = 7, minPts = 10, arrow_size = 0.1, pt.size = 0.5, arrowtip_size = 2, cols = NULL, label = FALSE, label.size = 15, spread = 5, min_dist = 0.25, n_neighbors = NA) {
  if (is.na(n_neighbors)){
    n_neighbors = nrow(data)/10
  }
  umap_cyto <- umap(data, spread = spread, min_dist = min_dist, n_neighbors = n_neighbors, random_state = 123)
  dat_umap <- umap_cyto$layout %>% data.frame()
  colnames(dat_umap) <- c("UMAP1", "UMAP2")
  
  # graph the umap plot
  p <- ggplot(data = dat_umap, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(color = groups), size = 4) +
    theme_void() 
  # y-range
  yrange = layer_scales(p)$y$range$range
  # x-range
  xrange = layer_scales(p)$x$range$range
  p <- p + 
    theme(legend.title = element_blank()) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1], yend = yrange[1] + (yrange[2]-yrange[1])*arrow_size, size = 0.8, arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] + (xrange[2] - xrange[1])*0.01, y = yrange[1], label = "UMAP1"), hjust = 0, vjust = 1, size = 4) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1] + (xrange[2]-xrange[1])*arrow_size, yend = yrange[1], size = 0.8,  arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] - (xrange[2] - xrange[1])*0.02, y = yrange[1] + (xrange[2] - xrange[1])*0.01, label = "UMAP2"), angle = 90, hjust = 0, vjust = 1, size = 4) +
    guides(colour = guide_legend(ncol = 1, label.theme = element_text(face = "bold", size = label.size),override.aes = list(size = 6))) +
    ggtitle(label = "")
  return(p)
}

# graphing pretty umap graph
graph_UMAP <- function(data_umap, groups, arrow_size = 0.1, pt.size = 0.5, arrowtip_size = 2, cols = NULL, label = FALSE, label.size = 15) {
  # graph the umap plot
  p <- ggplot(data = data_umap, aes(x = UMAP1, y = UMAP2)) + 
    geom_point(aes(color = as.factor(groups)), size = 0.5) +
    theme_void() 
  # y-range
  yrange = layer_scales(p)$y$range$range
  # x-range
  xrange = layer_scales(p)$x$range$range
  p <- p + 
    theme(legend.title = element_blank()) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1], yend = yrange[1] + (yrange[2]-yrange[1])*arrow_size, size = 0.8, arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] + (xrange[2] - xrange[1])*0.01, y = yrange[1], label = "UMAP1"), hjust = 0, vjust = 1, size = 4) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1] + (xrange[2]-xrange[1])*arrow_size, yend = yrange[1], size = 0.8,  arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] - (xrange[2] - xrange[1])*0.02, y = yrange[1] + (xrange[2] - xrange[1])*0.01, label = "UMAP2"), angle = 90, hjust = 0, vjust = 1, size = 4) +
    guides(colour = guide_legend(ncol = 1, label.theme = element_text(face = "bold", size = label.size),override.aes = list(size = 6))) +
    ggtitle(label = "")
  return(p)
}








