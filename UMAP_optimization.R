## This module is made to better optimize the UMAP visualization for any omics
require(gridExtra) # needed for the grid plot arrangement
require(pathviewr) # needed for find_curve_elbow for PCA post hoc analysis

# This function automatically find the best PCA dimensional cut-off for the subsequent UMAP
pca_to_umap <- function(data){
  # calculate the principal components
  pca <- prcomp(data, scale. = TRUE, center = TRUE)
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
pub_UMAP <- function(data, groups, arrow_size = 0.1, pt.size = 3, arrowtip_size = 2, cols = NULL, label = FALSE, label.size = 15) {
  umap_cyto <- umap(data, spread = 5, min_dist = 0.25, n_neighbors = 9, random_state = 123)
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
    geom_text(aes(x = xrange[1] + 0.2, y = yrange[1], label = "UMAP1"), hjust = 0, vjust = 1, size = 4) +
    geom_segment(x = xrange[1], y = yrange[1], xend = xrange[1] + (xrange[2]-xrange[1])*arrow_size, yend = yrange[1], size = 0.8,  arrow = arrow(length = unit(arrowtip_size,"mm"), type = "closed")) +
    geom_text(aes(x = xrange[1] - 0.6, y = yrange[1], label = "UMAP2"), angle = 90, hjust = 0, vjust = 1, size = 4) +
    guides(colour = guide_legend(ncol = 1, label.theme = element_text(face = "bold", size = label.size),override.aes = list(size = 6))) +
    ggtitle(label = "")
  return(p)
}