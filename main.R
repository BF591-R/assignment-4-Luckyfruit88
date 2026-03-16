library('tidyverse')
library('RColorBrewer')

#' Read the expression data "csv" file as a dataframe, not tibble
#'
#' @param filename (str): the path of the file to read
#' @param delimiter (str): generalize the function so it can read in data with
#'   your choice of delimiter
#'
#' @return A dataframe containing the example intensity data with rows as probes
#'   and columns as samples
#' @export
#'
#' @examples
read_data <- function(intensity_data, delimiter) {
    return(NULL)
}

###answer#####

read_data <- function(filename, delimiter) {
  df <- read.table(
    file = filename,
    sep = delimiter,
    header = TRUE,
    row.names = 1,
    check.names = FALSE
  )
  return(df)
}


#' Define a function to calculate the proportion of variance explained by each PC
#'
#' @param pca_results (obj): the results returned by `prcomp()`
#'
#' @return A vector containing the values of the variance explained by each PC
#' @export
#'
#' @examples
calculate_variance_explained <- function(pca_results) {
    return(NULL)
}

####answer#####
calculate_variance_explained <- function(pca_results){
  pc_variance <- pca_results$sdev^2
  
  prop_variance_explained <- pc_variance / sum(pc_variance)
  
  return(prop_variance_explained)
}


#' Define a function that takes in the variance values and the PCA results to
#' make a tibble with PC names, variance explained by each PC, and the
#' cumulative sum of variance explained. These columns should be named 
#' "principal_components", "variance_explained", and "cumulative", respectively.
#' 
#'
#'
#' @param pca_ve (vector): the vector generated in the previous function with
#'   the variance explained values
#' @param pca_results (object): the results returned by `prcomp()`
#' @return A tibble that contains the names of the PCs, the individual variance
#'   explained, and the cumulative variance explained with names described above
#' @export
#' @examples 
make_variance_tibble <- function(pca_ve, pca_results) {
    return(NULL)
}

#####answer###

make_variance_tibble <- function(pca_ve, pca_results){
  pc_names <- colnames (pca_results$rotation)
  cumulative_variance <- cumsum (pca_ve)
  variance_tibble <- tibble::tibble(
    principal_components = pc_names,
    variance_explained = pca_ve,
    cumulative = cumulative_variance
  )
  return(variance_tibble)
}


#' Define a function to create a biplot of PC1 vs. PC2 labeled by
#' SixSubTypesClassification
#'
#' @param metadata (str): The path to the proj_metadata.csv file
#' @param pca_results (obj): The results returned by `prcomp()`
#'
#' @return A ggplot consisting of a scatter plot of PC1 vs PC2 labeled by
#'   SixSubTypesClassification found in the metadata
#' @export
#'
#' @examples
make_biplot <- function(metadata, pca_results) {
    return(NULL)
}

####answer#######

make_biplot <- function(metadata, pca_results) {
  meta_df <- readr::read_csv(metadata, show_col_types = FALSE)
  
  pca_df <- as.data.frame(pca_results$x)
  
  pca_df$geo_accession <- rownames(pca_df)
  plot_data <- dplyr::inner_join(pca_df, meta_df, by = "geo_accession")
  
  biplot <- ggplot2::ggplot(
    plot_data, 
    ggplot2::aes(x = PC1, y = PC2, color = SixSubtypesClassification)
  ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA Biplot of PC1 vs PC2",
      x = "PC1",
      y = "PC2",
      color = "Subtype Classification"
    )
  
  return(biplot)
}


#' Define a function to return a list of probeids filtered by signifiance
#'
#' @param diff_exp_tibble (tibble): A tibble containing the differential expression results
#' @param fdr_threshold (float): an appropriate FDR threshold, we will use a
#'   value of .01. This is the column "padj" in the tibble.
#'
#' @return A list with the names of the probeids passing the fdr_threshold
#' @export
#'
#' @examples
list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
    return(NULL)
}

####answer#####

list_significant_probes <- function(diff_exp_tibble, fdr_threshold) {
  sig_probes <- diff_exp_tibble %>%
    dplyr::filter(!is.na(padj))%>%
    dplyr::filter(padj < fdr_threshold)%>%
    dplyr::pull(probeid)
  
  return(sig_probes)
}

#' Define a function that uses the list of significant probeids to return a
#' matrix with the intensity values for only those probeids.
#' @param intensity (dataframe): The dataframe of intensity data generated in
#'   part 1
#' @param sig_ids_list (list/vector): The list of differentially expressed
#'   probes generated in part 6
#'
#' @return A `matrix()` of the probe intensities for probes in the list of
#'   significant probes by FDR determined in the previous function.
#'
#' @export
#'
#' @examples
return_de_intensity <- function(intensity, sig_ids_list) {
    return(NULL)
}

#######answer####

return_de_intensity <- function(intensity, sig_ids_list) {
  sig_ids <- unlist (sig_ids_list)
  
  subset_intensity <- intensity [rownames(intensity) %in% sig_ids, ]
  
  intensity_matrix <- as.matrix(subset_intensity)
  return(intensity_matrix)
}

#' Define a function that takes the intensity values for significant probes and
#' creates a color-blind friendly heatmap
#'
#' @param de_intensity (matrix): The matrix of intensity values for significant
#'   differentially expressed probes returned in part 7
#' @param num_colors (int): The number of colors in a specificed RColorBrewer
#'   palette
#' @param palette (str): The name of the chosen RColorBrewer palette
#'
#' @return A heatmap displaying the intensity values for the differentially
#'   expressed probes
#' @export
#'
#' @examples
plot_heatmap <- function(de_intensity, num_colors, palette) {
    return(NULL)
}

#######answer######

plot_heatmap <- function(de_intensity, num_colors, palette) {
  heatmap_colors <- RColorBrewer::brewer.pal(n = num_colors, name = palette)
  
  color_func <- colorRampPalette(heatmap_colors)(256)
  
  hm <- heatmap(
    x = de_intensity,
    col = color_func,
    scale = "row",
    main = "Differentially Expressed Probes Heatmap",
    xlab = "Samples",
    ylab = "Significant Probes",
    margins = c(10, 10)  
  )
  return(hm)
}                                    

