#### Organized Olink HT QC Pipeline #####

#' Olink_Intensity_Norm
#'
#' @param data_olink
#'
#' @importFrom dplyr group_by mutate ungroup
#'
#' @return
#' @export
#'
#' @examples
Olink_Intensity_Norm <- function(data_olink) {
  # calculate ExtNPX based on counts in the Olink parquet file
  data_olink |>
    dplyr::group_by(SampleID, PlateID, Block) |>
    dplyr::mutate(LogProtExp = ExtNPX + log2(1e5)) |>
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
      .groups = "drop"
    )

  # Calculating Plate Control coefficient of variance
  ht_pc_vals <- data_olink |>
    dplyr::filter(SampleType == "PLATE_CONTROL") |>
    dplyr::group_by(Assay, OlinkID) |>
    dplyr::summarise(pc_cv = 100 * sd(LogProtExp) / mean(LogProtExp), .groups = "drop") |>
    dplyr::mutate(
      high_var_assay =
        dplyr::case_when(
          pc_cv > 20 ~ "High Variance",
          T ~ "Pass"
        )
    ) |>
    dplyr::select(-pc_cv)

  # This is the "sample level" qc, calculates ith sample in jth assay that needs to be replaced with zero or LLOQ
  # also labels those values in a new column - sample_level_qc
  ht_scaled_npx_sample <- data_olink |>
    dplyr::filter(SampleType == "SAMPLE") |>
    dplyr::left_join(y = ht_nc_vals, by = dplyr::join_by("Assay" == "OlinkID")) |>
    dplyr::left_join(y = ht_pc_vals, by = dplyr::join_by("Assay" == "OlinkID")) |>
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
    dplyr::summarise(percentage = 100 * n() / n_samples, .groups = "drop") |>
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
    dplyr::select(-c(`Below LLOD`, `Below LLOQ`))

  # For samples that didn't include the assay_level_qc call,
  # we can determine that there are 0 samples in that assay
  # that were below the LLOD, and that assay is good quality and should be
  # called continuous. This data only contains the patient samples themselves
  ht_scaled_npx_qc <-
    data_olink |>
    dplyr::left_join(
      y = ht_scaled_npx_assay,
      by = dplyr::join_by("OlinkID" == "Assay")
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


### HERE IS WHERE THE ISSUE IS ####
# Input single parquet file or list of parquet npx_files that need to be qc'd
#' olink_qc_run
#'
#' @param npx_file
#'
#' @importFrom arrow open_dataset read_parquet write_parquet
#' @importFrom dplyr filter
#' @importFrom purrr imap
#'
#' @return
#' @export
#'
#' @examples
olink_qc_run <- function(npx_file) {
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
      lapply(
        X = npx_file,
        FUN = arrow::read_parquet
      )
    message("Filtering...")
    filter_dat <-
      lapply(
        X = single_dat,
        FUN = \(y) dplyr::filter(y, !is.na(ExtNPX))
      )
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
    message("Writing QC'ed files to disk...")
    purrr::imap(
      .x = pre_comb_dat,
      .f = \(x, i) {
        arrow::write_parquet(
          x = x,
          sink = paste0("Run_Specific_QC_", stringr::str_split_i(string = npx_file[[i]], pattern = "/", i = -1) |> stringr::str_remove(pattern = ".parquet"), ".parquet")
        )
      }
    )
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
    message("")
    pre_comb_dat <-
      pre_comb(data_olink = scaled_dat)
  } else {
    (
      stop("length of file is less than 1")
    )
  }

  pre_comb_dat
}



#' olink_qc_project
#'
#' @param npx_file
#'
#' @importFrom dplyr bind_rows
#'
#' @return
#' @export
#'
#' @examples
olink_qc_project <- function(npx_file) {
  comb_dat <- dplyr::bind_rows(pre_comb_dat)
  post_comb(comb_dat)
}

# dat <- olink_qc_run(npx_file = c("\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_1-2_NPX_2024-02-07.parquet",
#                          "\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_3-4_NPX_2024-02-07.parquet",
#                          "\\\\data/ADI/Phenotyping_Core/01_Instrumentation/Olink_HT/01_Data/2024-02-07_AMP_SLE/AMP_SLE_5-6_NPX_2024-02-12.parquet"))
