#### Organized Olink HT QC Pipeline #####

#' Olink_Intensity_Norm
#'
#' @param data_olink
#'
#' @importFrom dplyr filter group_by mutate left_join join_by select group_split bind_rows
#' @importFrom arrows write_parquet
#' @importFrom lubridate as_date
#' @importFrom stringr str_extract
#' @return
#' @export
#'
#' @examples

project_split <- function(data_olink, manifest, level){
  
  sample_key <- manifest |>
    dplyr::filter(project != "Bridge") |> 
    dplyr::filter(sample_type == "SAMPLE") |> 
    dplyr::select(sample_id, project) 
  
  bridge_key <- manifest |> 
    dplyr::filter(project == "Bridge") |> 
    dplyr::select(sample_id, project)
  
  control_dat <- data_olink |> 
    dplyr::filter(SampleType != "SAMPLE") |> 
    dplyr::mutate(project = "CONTROL")
  
  sample_df_list <- data_olink |> 
    dplyr::filter(SampleID %in% sample_key$sample_id) |> 
    dplyr::left_join( y = sample_key,
                      by = dplyr::join_by("SampleID" == "sample_id")) |> 
    dplyr::group_by(project) |> 
    dplyr::group_split()
  
  
  df_list <- lapply(sample_df_list, function(x) dplyr::bind_rows(x, control_dat))
  
  
  
  
  if(nrow(bridge_key) >= 1){
    message("Writing level 1 parquet...")
    
    bridge_dat <- data_olink |> 
      dplyr::filter(SampleID %in% bridge_key$sample_id) |> 
      dplyr::mutate(project = "Bridge")
    
    bridge_controls <- dplyr::bind_rows(bridge_dat, control_dat)
    
    arrow::write_parquet(x = bridge_controls , 
                         sink = paste0(unique(bridge_controls[bridge_controls$project != "CONTROL",]$project),
                                       "_",
                                       stringr::str_extract(unique(bridge_controls$file_name),
                                                            "\\d+(-\\d+)_NPX.*$"
                                       )
                         )
    )
    
    lapply(df_list, 
           function(x) {
             arrow::write_parquet(x = x,
                                  sink = 
                                    paste0(level, "_",
                                           unique(x[x$project != "CONTROL",]$project),
                                           "_",
                                           stringr::str_extract(
                                             unique(x$file_name),
                                             "\\d+(-\\d+)_NPX.*$"
                                           )
                                    )
             )
           }
    )
    
    
  } else{
    message("Writing Level 2 parquet...")
    lapply(df_list, 
           function(x) 
             arrow::write_parquet(x = x,
                                  sink = 
                                    paste0(level, "_",
                                           unique(x[x$project != "CONTROL",]$project),
                                           "_",
                                           stringr::str_extract(
                                             unique(x$file_name),
                                             "\\d+(-\\d+)_NPX.*$"
                                           )
                                    )
             )
    )
    
    
  }
}






Olink_Intensity_Norm <- function(data_olink) {
  # calculate ExtNPX based on counts in the Olink parquet file
  data_olink |>
    dplyr::group_by(SampleID, PlateID, Block) |>
    dplyr::mutate(LogProtExp = ExtNPX + log2(1e5)) #|>
    dplyr::ungroup()
}

#' pre_comb
#'
#' @param data_olink
#'
#' @importFrom dplyr filter group_by summarise join_by mutate case_when
#' @importFrom stats median quantile
#'
#' @return
#' @export
#'
#' @examples
pre_comb <- function(data_olink) {
  ht_nc_vals <- data_olink |>
    dplyr::filter(SampleType == "NEGATIVE_CONTROL") |>
    dplyr::group_by(Assay, OlinkID) |>
    dplyr::summarise(
      median_nc = median(LogProtExp),
      iqr_nc = quantile(LogProtExp, 0.75),
      .groups = 'drop'
    ) |> 
    select(-Assay)
  
  # Calculating Plate Control coefficient of variance
  ht_pc_vals <- data_olink |>
    dplyr::filter(SampleType == "PLATE_CONTROL") |>
    dplyr::group_by(Assay, OlinkID) |>
    dplyr::summarise(pc_cv = 100*sd(LogProtExp)/mean(LogProtExp), .groups = 'drop') |>
    dplyr::mutate(
      high_var_assay =
        dplyr::case_when(
          pc_cv > 20 ~ "High Variance",
          T ~ "Pass"
        )
    ) |>
    dplyr::select(-c(pc_cv, Assay))
  
  # This is the "sample level" qc, calculates ith sample in jth assay that needs to be replaced with zero or LLOQ
  # also labels those values in a new column - sample_level_qc
  ht_scaled_npx_sample <- data_olink |>
    dplyr::filter(SampleType == "SAMPLE") |>
    dplyr::left_join(y = ht_nc_vals, by = dplyr::join_by("OlinkID" == "OlinkID")) |>
    dplyr::left_join(y = ht_pc_vals, by = dplyr::join_by("OlinkID" == "OlinkID")) |>
    dplyr::mutate(
      sample_level_qc =
        dplyr::case_when(
          LogProtExp < median_nc ~ "Below LLOD",
          LogProtExp < iqr_nc ~ "Below LLOQ",
          T ~ "Pass"
        )
    ) |>
    dplyr::mutate(
      LogProtExp =
        dplyr::case_when(
          LogProtExp < median_nc ~ 0,
          LogProtExp < iqr_nc ~ iqr_nc,
          T ~ LogProtExp
        )
    )
  
  ht_scaled_npx_sample
}

#' post_comb
#'
#' @param data_olink
#'
#' @importFrom dplyr group_by summarise filter mutate case_when select left_join join_by pull
#' @importFrom tidyr pivot_wider
#' @return
#' @export
#'
#' @examples
post_comb <- function(data_olink) {
  # specifying number of samples present in the total combined dataset
  n_samples <-
    data_olink |>
    dplyr::pull(SampleID) |>
    unique() |>
    length()
  
  
  # Assay level QC - if 50% of samples are below LLOQ, labeled as semi-continuous
  # if 75% of samples are below LLOD, labeled as categorical
  # Test to adjust how I calculate categorical, semi-continuous, or continuous
  ht_scaled_npx_assay <-
    data_olink |>
    dplyr::group_by(Assay, OlinkID, sample_level_qc) |>
    dplyr::summarise(percentage = 100*dplyr::n()/n_samples, .groups = 'drop') |>
    dplyr::filter(sample_level_qc %in% c("Below LLOD", "Below LLOQ")) |>
    tidyr::pivot_wider(names_from = sample_level_qc, values_from = percentage) |>
    dplyr::mutate(
      `Below LLOD` =
        dplyr::case_when(
          is.na(`Below LLOD`) == T ~ 0,
          T ~ `Below LLOD`
        ),
      `Below LLOQ` =
        dplyr::case_when(
          is.na(`Below LLOQ`) == T ~ 0,
          T ~ `Below LLOQ`
        ),
      `Below LLOQ` = `Below LLOQ` + `Below LLOD`,
      assay_level_qc =
        dplyr::case_when(
          `Below LLOD` > 75 ~ "Categorical",
          `Below LLOQ` > 50 ~ "Semi-Continuous",
          T ~ "Continuous"
        )
    ) |>
    dplyr::select(-c(`Below LLOD`, `Below LLOQ`, Assay))
  
  # For samples that didn't include the assay_level_qc call,
  # we can determine that there are 0 samples in that assay
  # that were below the LLOD, and that assay is good quality and should be
  # called continuous. This data only contains the patient samples themselves
  ht_scaled_npx_qc <-
    data_olink |>
    dplyr::left_join(
      y = ht_scaled_npx_assay,
      by = dplyr::join_by("OlinkID" == "OlinkID")
    ) |>
    dplyr::mutate(
      assay_level_qc =
        dplyr::case_when(
          is.na(assay_level_qc) == T ~ "Continuous",
          T ~ assay_level_qc
        )
    ) |>
    dplyr::select(-c(median_nc, iqr_nc))
  
  ht_scaled_npx_qc
}

#'
#'
#' @param npx_file
#'
#' @importFrom arrow open_dataset read_parquet write_parquet
#' 
#' @importFrom purrr imap
#' @importFrom readxl read_excel cell_cols
#' @importFrom janitor clean_names
#' 

olink_qc_level1 <- function(npx_file, manifest_file){
  
  message(stringr::str_glue("Found {length(npx_file)} parquet files..."))
  if (length(npx_file) > 1) {
    message("Multiple Olink HT runs detected")
    colname_list <- lapply(
      X = npx_file,
      FUN = \(y) arrow::open_dataset(y)$schema$names
    )
    colname_list_len <- length(unique(lapply(colname_list, \(x) sort(toupper(x)))))
    if (colname_list_len != 1) {
      warning("File columns do not match. Check parquet files for differences")
    }
    
    message("Reading files...")
    single_dat <-
      purrr::imap(
        .x = npx_file,
        .f = \(x, i) {
          arrow::read_parquet(x) |> dplyr::mutate(file_name = basename(npx_file[[i]]))
        }
      )
    
    
    message("Separating projects and writing raw files to disk...")
    manifest_read <- lapply(X = manifest_file,
                            FUN = function(x) readxl::read_excel(x, range = readxl::cell_cols("A:F")) |> 
                              janitor::clean_names())
    manifest_ref <- dplyr::bind_rows(manifest_read)
    lapply(single_dat, function(x) project_split(x, manifest = manifest_ref, level = "Level 1"))
    
  } else if (length(npx_file) == 1) {
    message("Single Olink HT run detected")
    single_dat <- arrow::read_parquet(npx_file) |> dplyr::mutate(file_name = basename(npx_file[[i]]))
    
    
    message("Separating projects and writing raw files to disk...")
    manifest_read <- readxl::read_excel(manifest_file, 
                                        range = readxl::cell_cols("A:F")) |> 
      janitor::clean_names()
    
    project_split(single_dat, manifest = manifest_read, level = "Level_1")
  }
  
  
  
}




# Input single parquet file or list of parquet npx_files that need to be qc'd
#' olink_qc_run
#'
#' @param npx_file
#'
#' @importFrom arrow open_dataset read_parquet write_parquet
#' @importFrom dplyr filter `%>%`
#' @importFrom purrr imap
#' @importFrom readxl read_excel
#' @importFrom janitor clean_names
#' @importFrom dplyr bind_rows
#'
#' @return
#' @export
#'
#' @examples
olink_qc_level2 <- function(npx_file){
  
  message(stringr::str_glue("Found {length(npx_file)} parquet files..."))
  if (length(npx_file) > 1) {
    message("Multiple Olink HT runs detected")
    colname_list <- lapply(
      X = npx_file,
      FUN = \(y) arrow::open_dataset(y)$schema$names
    )
    colname_list_len <- length(unique(lapply(colname_list, \(x) sort(toupper(x)))))
    if (colname_list_len != 1) {
      warning("File columns do not match. Check parquet files for differences")
    }
    
    # Reading in output from the previous function
    
    message("Reading in files...")
    single_dat <-
      lapply(
        X = npx_file, 
        FUN = arrow::read_parquet
      ) 
    
    message("Filtering...")
    sample_fail <- unlist(lapply(X = single_dat, 
                                 FUN = function(y) y |> 
                                   dplyr::filter(
                                     SampleQC == "FAIL") |> 
                                   dplyr::group_by(SampleID) |> 
                                   dplyr::summarise(.data = _, .groups = 'drop') |>
                                   dplyr::pull(SampleID)))
    
    
    filter_dat <-
      lapply(
        X = single_dat,
        FUN = \(y) dplyr::filter(y, !SampleID %in% sample_fail)
      )
    
    message(stringr::str_glue("Samples ", paste0(sample_fail, collapse = ', '), " were removed"))
    message("Scaling...")
    scaled_dat <-
      lapply(
        X = filter_dat,
        FUN = Olink_Intensity_Norm
      )
    
    # Calculating LLOD and LLOQ based on negative control values
    message("Calculating LLOD and LLOQ values...")
    pre_comb_dat <-
      lapply(
        X = scaled_dat,
        FUN = pre_comb
      )
    
    message("Combining multiple runs...")
    comb_dat <- dplyr::bind_rows(pre_comb_dat)
    
    message("Calculating assay level qc...")
    post_comb_dat <- post_comb(comb_dat)
    
    
    
    
    file_entry <- stringr::str_extract(unique(post_comb_dat$file_name),
                                       "\\d+[-]\\d+[-]\\d+")
    
    file_date <- lubridate::as_date(file_entry)
    
    
    # controls may be present as well
    if(length(unique(comb_dat$project)) == 1){
      
      message("Writing QC'ed files to disk...")
      arrow::write_parquet(x = post_comb_dat |> dplyr::select(-c(file_name, project)), 
                           sink = paste0("Level_2_",
                                         unique(post_comb_dat$project),
                                         "_Runs_",
                                         as.character(min(file_date)),
                                         "_",
                                         as.character(max(file_date)),
                                         ".parquet"))
    } else{
      
      warning(stringr::str_glue("{length(unique(comb_dat$project))} unique projects have been combined"))
      
      message("Writing QC'ed files to disk...")
      project_name <- paste(unique(post_comb_dat$project), collapse = "_")
      
      arrow::write_parquet(x = post_comb_dat |> dplyr::select(-c(file_name, project)), 
                           sink = paste0("Level_2_",
                                         project_name,
                                         "_Runs_",
                                         as.character(min(file_date)),
                                         "_",
                                         as.character(max(file_date)),
                                         ".parquet"))
      
    }
    
    
  } else if (length(npx_file) == 1) {
    message("Single Olink HT run detected")
    dat <- arrow::read_parquet(npx_file)
    message("Filtering...")
    filter_dat <-
      dplyr::filter(
        dat,
        is.na(ExtNPX) == FALSE
      )
    message("Scaling...")
    scaled_dat <-
      Olink_Intensity_Norm(data_olink = filter_dat)
    
    message("Calculating LLOD and LLOQ...")
    pre_comb_dat <-
      pre_comb(data_olink = scaled_dat)
    
    message("Calculating assay_level_qc")
    post_comb_dat <- post_comb(pre_comb_dat)
    
    
    message("Writing files to disk...")
    arrow::write_parquet(x = post_comb_dat |> dplyr::select(-c(file_name, project)),
                         sink = 
                           paste0("Level 2", "_",
                                  unique(post_comb_dat[post_comb_dat$project != "CONTROL",]$project),
                                  "_",
                                  stringr::str_extract(
                                    unique(post_comb_dat$file_name),
                                    "\\d+(-\\d+)_NPX.*$"
                                  )
                           )
    )
    
  } else (
    stop('length of file is less than 1')
  )
  
  post_comb_dat
}



# base_dir <- c(wherever all the files are)
# file_list <- paste0(base_dir, list.files(base_dir, recursive = TRUE))

# dat <- olink_qc_level1(npx_file = parquet_file_list, manifest_file = file_list)
# 
# 
# 
# dat2 <- olink_qc_level2(c("Level 1_AMP SLE_1-2_NPX_2024-02-07.parquet",
#                            "Level 1_AMP SLE_3-4_NPX_2024-02-07.parquet",
#                            "Level 1_AMP SLE_5-6_NPX_2024-02-12.parquet",
#                            "Bridge_1-2_NPX_2024-02-07.parquet"))




