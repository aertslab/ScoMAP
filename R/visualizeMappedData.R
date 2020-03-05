#' VM_loom
#'
#' Generates a minimal loom file to visualize data on SCope (http://scope.aertslab.org).
#' @param VM Virtual map data frame, after running mapCells()
#' @param omics_mat Matrix containing features as rows (e.g. genes, regions) and real cells as columns.
#' @param loom_path Path where the loom file will be generated
#' @param genome Genome version to which data belongs
#'
#' @examples
#' 
#' VM_loom(VM, dgem, loom_path='VM.loom', genome='dm6')
#' 
#' @export

VM_loom <- function(VM, omics_mat ,loom_path, genome){
  # Check
  if (!'RM_assignment' %in% colnames(VM)){
    stop('Please, run mapCells() first.')
  }
  if(sum(VM$RM_assignment %in% colnames(omics_mat)) != nrow(VM)){
    stop('Are the cells asisgned to this virtual map on the omics feature matrix?')
  }
  # Get profiles from the real cells mapped
  omics_mat <- omics_mat[,VM$RM_assignment]
  colnames(omics_mat) <- rownames(VM)
  # Get virtual map coordinates
  default.coordinates <- VM[, c('x', 'y')]
  default.coordinates.name <- 'Virtual Map'
  # Get metadata
  metadata <- VM[,c('cluster_annot', 'is.landmark', 'RM_cell_type')]
  metadata[is.na(metadata[,'is.landmark']),'is.landmark'] <- ''
  colnames(metadata) <- c('Spatial_domain', 'Landmark', 'Cell_type')
  # Build minimal loom
  build_loom(
    file.name=loom_path,
    dgem= omics_mat,
    title="Virtual map",
    genome=genome, 
    default.embedding=default.coordinates,
    default.embedding.name=default.coordinates.name 
  )
  # Add metadata
  loom <- open_loom(loom_path)
  Spatial_domain <- as.vector(unlist(metadata[,'Spatial_domain']))
  names(Spatial_domain) <- rownames(metadata)
  add_col_attr(loom=loom, key = 'Spatial_domain', value=Spatial_domain, as.annotation = T)
  Landmark <- as.vector(unlist(metadata[,'Landmark']))
  names(Landmark) <- rownames(metadata)
  add_col_attr(loom=loom, key = 'Landmark', value=Landmark, as.annotation = T)
  Cell_type <- as.vector(unlist(metadata[,'Cell_type']))
  names(Cell_type) <- rownames(Cell_type)
  add_col_attr(loom=loom, key = 'Cell_type', value=Cell_type, as.annotation = T)
  rm(loom)
}

#' PlotVMFeatures
#'
#' Plot omics data into the virtual template.
#' @param VM Virtual map data frame, after running mapCells()
#' @param omics_mat Matrix containing features as rows (e.g. genes, regions) and real cells as columns.
#' @param features Vector with the names of features to be plotted on the virtual template.  
#' Features will be assigned to the corresponding RGB channel depending on the order. 
#' @param thr Threshold (between 0 and 1) to color cells with low feature value in grey. Default: 0.
#'
#' @examples
#' 
#' PlotVMFeatures(VM, dgem, features=c('ato', 'hth', 'Optix'))
#' 
#' @export

PlotVMFeatures <- function(VM, omics_mat, features, thr=0){
  if (!'RM_assignment' %in% colnames(VM)){
    stop('Please, run mapCells() first.')
  }
  if(sum(VM$RM_assignment %in% colnames(omics_mat)) != nrow(VM)){
    stop('Are the cells asisgned to this virtual map on the omics feature matrix?')
  }
  if(sum(features %in% rownames(omics_mat)) < length(features)){
    stop('Features are not in the omics matrix.')
  }
  
  if (length(features) == 1){
    red <- features[1]
    geneMat <- omics_mat[red,,drop=FALSE]
    modCols <- list(red=red)
  } else if (length(features) == 2) {
    red <- features[1]
    green <- features[2]
    geneMat <- omics_mat[c(red, green),,drop=FALSE]
    modCols <- list(red=red, green=green)
  } else if (length(features) == 3) {
    red <- features[1]
    green <- features[2]
    blue <- features[3]
    geneMat <- omics_mat[c(red, green, blue),,drop=FALSE]
    modCols <- list(red=red, green=green, blue=blue)
  } else {
    stop('A minimum of 1 and maximum of 3 the features can be provided.')
  }
  coordinates <- VM[,c('x', 'y')]
  geneMat <- t(apply(geneMat, 1, function(x) (x-min(x))/(max(x)-min(x))))
  offColor <- "#c0c0c030" # Transparent light gray
  modCols <- lapply(modCols, function(x) sapply(x, function(gene) rownames(geneMat)[gene]))
  cellColChan <- sapply(modCols, function(modsCol) apply(as.matrix(geneMat[names(modsCol),]), 1, mean))
  if (length(features) == 1){
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], 0, 0, alpha=.8))
  } else if (length(features) == 2) {
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], 0, alpha=.8))
  } else if (length(features) == 3) {
    cellCol <- apply(cellColChan, 1, function(x) rgb(x["red"], x["green"], x["blue"], alpha=.8))
  }
  cellCol <- cellCol[VM$RM_assignment]
  names(cellCol) <- rownames(VM)
  cellCol[as.vector(which(rowSums(cellColChan) < thr))] <- offColor
  plot(coordinates, col=cellCol[order(cellCol[VM$RM_assignment])], pch=16, axes=FALSE, frame.plot=FALSE, ann=FALSE, cex=0.5)
  legend("topright", legend =  features, fill=c('red', 'green', 'blue'), cex=0.5)
}
