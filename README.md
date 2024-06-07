# Olink_HT Usage Instructions
Instructions for Olink HT Level 2 QC using functions from the organized_olink_qc_functions.R script


Set the working directory to one containing both raw parquet and internal manifest files  
1. ``` setwd(desired_directory) ```

Initialize the Standard Data Package folder (specified by the Data Standards Document) and check that manifests and parquet files contain the same names  
2. ```  file_list <-  Olink_Reader(".")   ```

Split out the raw parquet files by desired "PROJECT" present in the internal manifest files and write split out files to Level_1 subfolder within SDP  
3. ```  files_lvl1 <-  Olink_lvl1(olink_files = file_list, proj_names = "PROJECT")   ```

If more than 1 file is present, perform median batch correction between runs 
4.  ```  files_batch_corrected <-  batch_correction(data = files_lvl1$data, method = "median")   ```

Prepare the data for Level 2 (specified by Data Standards Document) QC by evaluating the lower limits of detection and quantification per run
5.  ``` files_lvl2_prep <- lapply(files_batch_corrected, Olink_lvl2_prep)  ```

Perform Level 2 QC 
6.  ``` files_lvl2 <- Olink_lvl2(data = files_lvl2_prep)  ```

Evaluate the data for batch effects or other technical issues using UMAPS of the raw Extension NPX, batch corrected Extension NPX, and batch corrected LogProtExp (specified by Data Standards Document)
7.  ``` plots <- normalization_check(files_lvl2)  ```

In our experience, we have found that addition of a step removing technical high variance analytes (~600) minimizes the presence of a false "batch effect" seen in UMAPs of multi-file projects  
8.
 ```
   high_var_analytes <- files_lvl2 %>%
                              filter(high_var_assay == "High Variance") %>%
                              group_by(OlinkID) %>% 
                              summarise(., .groups = 'drop')
        low_var_lvl2 <- files_lvl2 %>%
                           filter(!OlinkID %in% high_var_analytes$OlinkID)
        plots_low_var <- normalization_check(low_var_lvl2)
 ```
