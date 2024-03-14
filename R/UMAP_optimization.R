## This module is made to better optimize the UMAP visualization for any omics

#' optimized_pca
#'
#' @description
#' Determine the optimal number of principal components for a given set of data
#' and return the principal component values for just those PCs
#'
#'
#' @param data A wide
#'
#' @importFrom irlba prcomp_irlba
#' @importFrom dplyr mutate row_number select
#' @importFrom tibble as_tibble column_to_rownames
#' @importFrom pathviewr find_curve_elbow
#' @importFrom tidyselect num_range
#' @importFrom rlang `%||%`
#'
#' @return
#' @export
#'
#' @examples
optimized_pca <- function(data, num_pcs=10) {

  num_pcs <- num_pcs %||% ncol(data)

  pca_res <- irlba::prcomp_irlba(data, n = num_pcs, scale. = TRUE, center = TRUE)

  # get PCA importance from the principal component analysis
  res_pca <- summary(pca_res)[["importance"]] |>
    t() |>
    tibble::as_tibble() |>
    dplyr::mutate(PCs = dplyr::row_number())
  n_pcs <-
    pathviewr::find_curve_elbow(
      data_frame =
        dplyr::select(
          .data = res_pca,
          PCs,
          `Proportion of Variance`
          ),
      plot_curve = TRUE
      )
  pca_res |>
    purrr::pluck("x") |>
    tibble::as_tibble() |>
    dplyr::select(
      tidyselect::num_range(
        prefix = "PC",
        range = seq(n_pcs)
        )
    ) |>
    dplyr::mutate(sample = rownames(data)) |>
    tibble::column_to_rownames(sample)
}

# TODO: combine the three optimize functions and just switch the title/optimizing
# ranges since the code is nearly identical

#' optimize_n_neighbor
#'
#' @description
#' Optimization function for n_neighbors by setting the spread = 10,
#' min_dist = 0.1 as default by can be changed as needed
#'
#'
#' @param data
#' @param groups
#' @param spread
#' @param min_dist
#' @param min
#' @param max
#' @param step
#'
#' @importFrom rlang enquo
#' @importFrom dplyr select rename
#' @importFrom purrr map pluck
#' @importFrom tibble as_tibble
#' @importFrom umap umap
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom cowplot plot_grid
#'
#' @return
#' @export
#'
#' @examples
optimize_n_neighbor <- function(
    data,
    groups,
    spread = 10,
    min_dist = 0.1,
    min = 5,
    max = 16,
    step = 1
    ) {

  # we can handle both quoted and unquoted variables for "groups"
  groups <- rlang::enquo(groups)

  npx_df <-
    dplyr::select(
      .data = data,
      1,
      where(is.numeric)
    ) |>
    tibble::column_to_rownames(var = colnames(data)[[1]])
  md <- dplyr::select(.data = data, 1, {{groups}})

  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(x){

    umap_cyto <-
      umap::umap(
        d = npx_df,
        spread = spread,
        min_dist = min_dist,
        n_neighbors = x,
        seed = 123
      ) |>
      purrr::pluck("layout") |>
      tibble::as_tibble(rownames = colnames(data)[[1]]) |>
      dplyr::rename(UMAP1 = "V1", UMAP2 = "V2")

    umap_cyto |>
      dplyr::left_join(md) |>
      ggplot2::ggplot(
        mapping =
          ggplot2::aes(
            x = UMAP1,
            y = UMAP2,
            color = {{groups}}
          )
      ) +
      ggplot2::geom_point(size = 4) +
      ggplot2::theme_classic() +
      ggplot2::ggtitle(paste0("n_neighbors = ", x)) +
      ggplot2::theme(
        legend.text = ggplot2::element_text(size = 12),
        axis.text = ggplot2::element_text(size = 12)
      )
  })
  cowplot::plot_grid(plotlist = plots, ncol = 4)
}

#
#' optimize_spread
#'
#' @description
#' Find the optimal value to pass to umap for the spread parameter.  Run after
#' optimizating the value for n_neighbor.
#'
#'
#' @param data
#' @param groups
#' @param n_neighbor
#' @param min_dist
#' @param min
#' @param max
#' @param step
#'
#' @importFrom rlang enquo
#' @importFrom dplyr select rename
#' @importFrom purrr map pluck
#' @importFrom tibble as_tibble
#' @importFrom umap umap
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom cowplot plot_grid
#'
#' @return
#' @export
#'
#' @examples
optimize_spread <- function(
    data,
    groups,
    n_neighbor,
    min_dist = 0.1,
    min = 1,
    max = 15,
    step = 1
    ) {
  # spread

  groups <- rlang::enquo(groups)

  npx_df <-
    dplyr::select(
      .data = data,
      1,
      where(is.numeric)
    ) |>
    tibble::column_to_rownames(var = colnames(data)[[1]])
  md <- dplyr::select(.data = data, 1, {{groups}})

  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(spread){
      umap_cyto <-
        umap::umap(
          d = npx_df,
          spread = spread,
          min_dist = min_dist,
          n_neighbors = n_neighbor,
          seed = 123
        ) |>
        purrr::pluck("layout") |>
        tibble::as_tibble(rownames = colnames(data)[[1]]) |>
        dplyr::rename(UMAP1 = "V1", UMAP2 = "V2")

      umap_cyto |>
        dplyr::left_join(md) |>
        ggplot2::ggplot(
          mapping =
            ggplot2::aes(
              x = UMAP1,
              y = UMAP2,
              color = {{groups}}
            )
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(paste0("spread = ", spread)) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12)
        )
      }
    )
cowplot::plot_grid(plotlist = plots, ncol = 4)
}


# The final optimization function that is to fix the min_dist after n_neighbor and spread parameters have been fixed
#' Title
#'
#' @param data
#' @param groups
#' @param spread
#' @param n_neighbor
#' @param min
#' @param max
#' @param step
#'
#' @importFrom rlang enquo
#' @importFrom dplyr select rename
#' @importFrom purrr map pluck
#' @importFrom tibble as_tibble
#' @importFrom umap umap
#' @importFrom ggplot2 ggplot aes geom_point theme_classic ggtitle theme element_text
#' @importFrom cowplot plot_grid
#'
#' @return
#' @export
#'
#' @examples
optimize_min_dist <- function(
    data,
    groups,
    spread,
    n_neighbor,
    min = 0.01,
    max = 0.5,
    step = 0.03
    ) {

  groups <- rlang::enquo(groups)

  npx_df <-
    dplyr::select(
      .data = data,
      1,
      where(is.numeric)
    ) |>
    tibble::column_to_rownames(var = colnames(data)[[1]])
  md <- dplyr::select(.data = data, 1, {{groups}})

  plots <- purrr::map(
    .x = seq(min, max, step),
    .f = \(min_dist){
      umap_cyto <-
        umap::umap(
          d = npx_df,
          spread = spread,
          min_dist = min_dist,
          n_neighbors = n_neighbor,
          seed = 123
        ) |>
        purrr::pluck("layout") |>
        tibble::as_tibble(rownames = colnames(data)[[1]]) |>
        dplyr::rename(UMAP1 = "V1", UMAP2 = "V2")

      umap_cyto |>
        dplyr::left_join(md) |>
        ggplot2::ggplot(
          mapping =
            ggplot2::aes(
              x = UMAP1,
              y = UMAP2,
              color = {{groups}}
            )
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::theme_classic() +
        ggplot2::ggtitle(paste0("min_dist = ", min_dist)) +
        ggplot2::theme(
          legend.text = ggplot2::element_text(size = 12),
          axis.text = ggplot2::element_text(size = 12)
        )
    }
  )
  cowplot::plot_grid(plotlist = plots, ncol = 4)
}


# graph the clean reworked umap
#' Title
#'
#' @param data
#' @param groups
#' @param arrow_size
#' @param pt.size
#' @param arrowtip_size
#' @param cols
#' @param label
#' @param label.size
#'
#' @return
#' @export
#'
#' @examples
pub_UMAP <- function(
    data,
    groups,
    spread = 5,
    min_dist = 0.25,
    n_neighbors = 9,
    random_state = 123,
    arrow_size = 0.1,
    pt.size = 3,
    arrowtip_size = 2,
    cols = NULL,
    label = FALSE,
    label.size = 15
    ) {

  groups <- rlang::enquo(groups)

  npx_df <-
    dplyr::select(
      .data = data,
      1,
      where(is.numeric)
      ) |>
    tibble::column_to_rownames(var = colnames(data)[[1]])

  md <- dplyr::select(.data = data, 1, {{groups}})

  umap_cyto <- umap::umap(
    npx_df,
    spread = spread,
    min_dist = min_dist,
    n_neighbors = n_neighbors,
    random_state = random_state
    ) |>
  purrr::pluck("layout") |>
  tibble::as_tibble(rownames = colnames(data)[[1]]) |>
  dplyr::rename(UMAP1 = "V1", UMAP2 = "V2")

  message("Stuff!")
  p <-
    umap_cyto |>
    dplyr::left_join(md) |>
    ggplot2::ggplot(
      mapping =
        ggplot2::aes(
          x = UMAP1,
          y = UMAP2,
          color = {{groups}}
        )
    ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0("min_dist = ", min_dist)) +
    ggplot2::theme(
      legend.text = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 12)
    ) +
    ggplot2::theme_void()

  yrange <- ggplot2::layer_scales(p)[["y"]][["range"]][["range"]]
  xrange <- ggplot2::layer_scales(p)[["x"]][["range"]][["range"]]

  p + ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::geom_segment(
      x = xrange[1],
      y = yrange[1],
      xend = xrange[1],
      yend = yrange[1] + (yrange[2] - yrange[1]) * arrow_size,
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(
          arrowtip_size,
          "mm"
          ),
        type = "closed"
        )
      ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = xrange[1] + 0.2,
        y = yrange[1],
        label = "UMAP1"
        ),
      hjust = 0,
      vjust = 1,
      size = 4
      ) +
    ggplot2::geom_segment(
      x = xrange[1],
      y = yrange[1],
      xend = xrange[1] + (xrange[2] - xrange[1]) * arrow_size,
      yend = yrange[1],
      size = 0.8,
      arrow = grid::arrow(
        length = grid::unit(
          arrowtip_size,
          "mm"
          ),
        type = "closed"
        )
      ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = xrange[1] - 0.6,
        y = yrange[1],
        label = "UMAP2"
        ),
      angle = 90,
      hjust = 0,
      vjust = 1,
      size = 4
      ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = 1,
        label.theme = ggplot2::element_text(
          face = "bold",
          size = label.size
          ),
        override.aes = list(size = 6)
        )
      ) +
    ggplot2::ggtitle(label = "")
}
