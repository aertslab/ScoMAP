#' readTemplate
#'
#' Read virtual template (in jpeg) into R. The spatially located regions must be coloured differently,
#' and no borders should be used.
#' @param pathToImage Path to the jpeg or png image. The number of virtual cells will be determined by the 
#' size of the picture.
#' @param k Number of spatial domains in the image
#' @param bgThr Threshold (between 0 and 1) above which all the color channels (RGB) have to be above for the
#' pixel to be considered background. Default:0.93
#' @param neighbours Number of neighbours to be used to refine color cluster assignment. It must be the same length 
#' as the number of refinement rounds. Default: c(20, 5).
#' @param refineRound Number of cluster refinement rounds. Default:2
#' @param plot Boolean indicating whether to return the computed template colored by spatial sections. 
#' Default: TRUE.
#'
#' @return A virtual map data.frame containing the virtual cell coordinates, assigned spatial cluster (original and refined),
#' and cluster color
#'
#' @examples
#' virtualTemplate <- readTemplate(pathToJpeg, k=8, bgThr=0.93, neighbours=c(20,5), refineRound=2)
#'
#' @import jpeg
#' @import png
#' @import data.table
#' 
#' @export

readTemplate<- function(pathToImage,
                        k,
                        bgThr=0.93,
                        neighbours=c(20,5),
                        refineRound=2,
                        plot=TRUE){
  # Check
  if(length(neighbours) < refineRound){
    stop('The number of neighbours has not been indicated for all refinement rounds')
  }
  
  # Read image into R
  if (strsplit(pathToImage, split = "\\.")[[1]][2] == 'jpeg' | strsplit(pathToImage, split = "\\.")[[1]][2] == 'jpg'){
    VM <- readJPEG(pathToImage)
  } else if (strsplit(pathToImage, split = "\\.")[[1]][2] == 'png'){
    VM <- readPNG(pathToImage)
  } else {
    stop('Only jpeg or png formats are accepted')
  }
  
  # Obtain the dimension
  VMDm <- dim(VM)
  # Assign RGB channels to data frame
  VMRGB <- data.frame(
    x = rep(1:VMDm[2], each = VMDm[1]),
    y = rep(VMDm[1]:1, VMDm[2]),
    R = as.vector(VM[,,1]),
    G = as.vector(VM[,,2]),
    B = as.vector(VM[,,3]),
    color = rgb(VM[,,1], VM[,,2], VM[,,3], maxColorValue=1),
    row.names=paste0(rep(1:VMDm[2], each = VMDm[1]), '_', rep(VMDm[1]:1, VMDm[2])),
                     stringsAsFactors=FALSE
  )
  # Remove white (and close to white), pixels
  if (sum(VMRGB[,'R'] > bgThr & VMRGB[,'G'] > bgThr & VMRGB[,'B'] > bgThr) > 0){
    VMRGB <- VMRGB[-which(VMRGB[,'R'] > bgThr & VMRGB[,'G'] > bgThr & VMRGB[,'B'] > bgThr),]
  }
  
  # Color clustering
  col_d <- dist(VMRGB[, c("R", "G", "B")], method = "euclidean") # distance matrix
  suppressWarnings(fit <- hclust(col_d, method="ward.D"))
  clusters <- cutree(fit, k=k)
  VMRGB$cluster <- clusters
  centroids <- sapply(1:k, function(i) rgb(t(colMeans(VMRGB[which(VMRGB$cluster == i), c('R', 'G', 'B')])), maxColorValue=1))
  names(centroids) <- 1:k
  VMRGB$cluster_color <- centroids[clusters]
  
  # Refine clusters
  pos_d <- as.matrix(dist(VMRGB[, c('x', 'y')], method = "euclidean")) 
  cluster_refined <- sapply(colnames(pos_d), function (virtual_cell)names(which.max(table(VMRGB$cluster[order(pos_d[, virtual_cell])][1:neighbours[1]]))))
  VMRGB$cluster_refined <- cluster_refined
  VMRGB$cluster_refined_color <- centroids[VMRGB$cluster_refined]
  
  # Further refinement
  if (refineRound > 1){
    for (i in 2:refineRound){
      VMRGB$cluster_refined <- sapply(colnames(pos_d), function (virtual_cell)names(which.max(table(VMRGB$cluster_refined[order(pos_d[, virtual_cell])][1:neighbours[i]]))))
      VMRGB$cluster_refined_color <- centroids[VMRGB$cluster_refined]
    }
  }
  
  if (plot==TRUE){
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(VMRGB[,c('x', 'y')], col=VMRGB$cluster_refined_color, pch=16, axes=FALSE, xlab='', ylab='')
    legend("topright", inset=c(-0.2,0), legend =  names(centroids), fill=centroids, cex=0.5)
  }
  
  return(as.data.frame(VMRGB))
}

#' readTemplateBySection
#'
#' Read separate sections from virtual template (in jpeg) into R. One image has to be provided per spatial cluster.
#' 
#' @param pathToJpegFolder Path to the folder containin jpeg images. The number of virtual cells will be determined by the 
#' size of the picture.
#' @param plot Boolean indicating whether to return the computed template colored by spatial sections. 
#' Default: TRUE.
#'
#' @return A virtual map data.frame containing the virtual cell coordinates, assigned spatial cluster (based on the file name),
#' and cluster color
#'
#' @examples
#' virtualTemplate <- readTemplateBySection(pathToJpeg)
#'
#' @import jpeg
#' @import data.table
#' @export

readTemplateBySection <- function(pathToJpegFolder,
                        plot=TRUE){
  sectionPath <- paste0(pathToJpegFolder, list.files(pathToJpegFolder))
  sectionName <- gsub('.jpg', '', list.files(pathToJpegFolder))
  sectionColor <- .distinctColorPalette(length(sectionName))
  VM <- lapply(1:length(sectionPath), function (i) .readTemplateSection(sectionPath[i], i, sectionColor[i], sectionName[i]))
  VM <- rbindlist(VM)
  VM <- VM[!duplicated(paste0(as.vector(unlist(VM[,1])), '_', as.vector(unlist(VM[,2])))),]
  rownames(VM) <- paste0(as.vector(unlist(VM[,1])), '_', as.vector(unlist(VM[,2])))
  
  if (plot==TRUE){
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(VM[,c('x', 'y')], col=VM$cluster_refined_color, pch=16, axes=FALSE, xlab='', ylab='')
    legend("topright", inset=c(-0.2,0), legend =  paste0(sectionName, '(', 1:length(sectionPath), ')'), fill=sectionColor, cex=0.5)
  }
  return(as.data.frame(VM))
}

# Helper function
.readTemplateSection <- function(pathToJpeg,
                                 cluster,
                                 cluster_color, 
                                 cluster_annot){
  VM <- readJPEG(pathToJpeg)
  # Obtain the dimension
  VMDm <- dim(VM)
  # Assign RGB channels to data frame
  VMRGB <- data.frame(
    x = rep(1:VMDm[2], each = VMDm[1]),
    y = rep(VMDm[1]:1, VMDm[2]),
    R = as.vector(VM[,,1]),
    G = as.vector(VM[,,2]),
    B = as.vector(VM[,,3]),
    color = rgb(VM[,,1], VM[,,2], VM[,,3], maxColorValue=1),
    row.names=paste0(rep(1:VMDm[2], each = VMDm[1]), '_', rep(VMDm[1]:1, VMDm[2]))
  )
  kClusters <- 3
  kMeans <- kmeans(VMRGB[, c("R", "G", "B")], centers = kClusters)
  clusters <- kMeans$cluster
  centroids <- sapply(1:3, function(i) rgb(t(colMeans(VMRGB[which(clusters == i), c('R', 'G', 'B')])), maxColorValue=1))
  bg <- names(which(clusters == which(centroids == max(centroids))))
  VMRGB <- VMRGB[-which(rownames(VMRGB) %in% bg),]
  
  VMRGB$cluster <- rep(cluster, nrow(VMRGB)) 
  VMRGB$cluster_color <- rep(cluster_color, nrow(VMRGB)) 
  VMRGB$cluster_refined <- rep(cluster, nrow(VMRGB)) 
  VMRGB$cluster_refined_color <- rep(cluster_color, nrow(VMRGB))
  VMRGB$cluster_annot <- rep(cluster_annot, nrow(VMRGB))
  rownames(VMRGB) <- paste0(as.vector(unlist(VMRGB[,1])), '_', as.vector(unlist(VMRGB[,2])))
  
  return(as.data.frame(VMRGB))
}

# Helper function
.distinctColorPalette <-function(k) {
  set.seed(123)
  if(packageVersion("scales") >= '1.1.0'){
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=85)(2e3))))
  } else {
    ColorSpace <- t(unique(col2rgb(scales::hue_pal(l=60:100)(2e3))))
  }
  km <- kmeans(ColorSpace, k, iter.max=20)
  colors <- rgb(round(km$centers), maxColorValue=255)
  return(colors)
}

  
#' addClusterAnnotation
#'
#' Add cluster annotation to a VM data frame
#' @param VM Virtual map data frame
#' @param clusterAnnotation Vector with cluster names as values and cluster numbers as name
#'
#' @return A virtual map data.frame containing a cluster_annot column
#'
#' @examples
#' VM <- addClusterAnnotation(VM, clusterAnnotation)
#'
#' @export
addClusterAnnotation <- function(VM,
                           clusterAnnotation){
  VM$cluster_annot <- clusterAnnotation[VM$cluster_refined]
  return(as.data.frame(VM))
}

#' intercalateCells
#'
#' Intercalates two or more cell types within the same spatial cluster
#' @param VM Virtual map data frame
#' @param clusterAnnotation Vector with cluster names as values and cluster numbers as name
#'
#' @return A virtual map data.frame containing a cluster_annot column
#'
#' @examples
#' VM <- addClusterAnnotation(VM, clusterAnnotation)
#'
#' @export
intercalateCells <- function(VM,
                             targetAnnot,
                             subclusters){
  if (is.null(VM$cluster_annot)){
    stop('Please, provide cluster annotation with addClusterAnnotation().')
  }
  subclusters <- rep(subclusters, sum(VM$cluster_annot == targetAnnot))[1:sum(VM$cluster_annot == targetAnnot)]
  VM$cluster_annot[which(VM$cluster_annot == targetAnnot)] <- subclusters
  return(as.data.frame(VM))
}

intercalateCells <- function(VM,
                             targetAnnot,
                             subclusters){
  if (is.null(VM$cluster_annot)){
    stop('Please, provide cluster annotation with addClusterAnnotation().')
  }
  subclusters <- rep(subclusters, sum(VM$cluster_annot == targetAnnot))[1:sum(VM$cluster_annot == targetAnnot)][sample(1:sum(VM$cluster_annot == targetAnnot))]
  VM$cluster_annot[which(VM$cluster_annot == targetAnnot)] <- subclusters
  return(as.data.frame(VM))
}

#' plotAnnotatedVM
#'
#' Plots VM based on cluster annotation
#' @param VM Virtual map data frame
#' @param colVars Vector with selected colors and cluster annotation as names. Default=NULL.
#' @param inset Indicates the positon out of the plot the legend will be. Default: c(-0.2,0)
#' @param show.landmark Show landmark cells in red. Default: FALSE.
#' @examples
#' VM <- plotAnnotatedVM(VM, clusterAnnotation)
#'
#' @export
plotAnnotatedVM <- function(VM, colVar=NULL, inset=c(-0.2,0), show.landmark=FALSE){
  if (is.null(colVar)){
    sectionName <- unique(VM$cluster_annot)
    sectionColor <- .distinctColorPalette(length(sectionName))
    names(sectionColor) <- sectionName
    color <- sectionColor[VM$cluster_annot]
    if (show.landmark == TRUE){
      color[!is.na(VM$is.landmark)] <- 'red'
    }
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(VM[,c('x', 'y')], col=color, pch=16, axes=FALSE, xlab='', ylab='')
    legend("topright", inset=c(-0.2,0), legend =  sectionName, fill=sectionColor, cex=0.5)
  } else {
    sectionColor <- colVar
    color <- sectionColor[VM$cluster_annot]
    if (show.landmark == TRUE){
      color[!is.na(VM$is.landmark)] <- 'red'
    }
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(VM[,c('x', 'y')], col=color, pch=16, axes=FALSE, xlab='', ylab='')
    legend("topright", inset=inset, legend =  names(sectionColor)[which(names(sectionColor) %in% unique(VM$cluster_annot))], fill=sectionColor, cex=0.5)
  }
}

#' addExternalCluster
#'
#' Plots VM based on cluster annotation
#' @param VM Virtual map data frame
#' @param x X position to draw external cluster
#' @param y Y position to draw external cluster
#' @param number_cell Minimal number of cells that the external cluster should have. Default:40
#' @param color Color given to the new external cluster
#' @param cluster_annot Annotation of the new external cluster
#' 
#' @return A virtual map data.frame containing additional rows with the new cluster information
#'
#' @examples
#' VM <- addExternalCluster(VM, x=x, y=y, number_cells=40, color='black', 'Cluster_X')
#'
#' @export

addExternalCluster <- function(VM,
                                x,
                                y,
                                number_cells=40,
                                color,
                                cluster_annot){
  ExternalVM <- cbind(x,y)
  while(nrow(ExternalVM) < number_cells){
    additional_cells <- ExternalVM[1,]
    for (row in 1:nrow(ExternalVM)){
      x <- ExternalVM[row, 'x']
      y <- ExternalVM[row, 'y']
      additional_cells <- rbind(additional_cells,
                                cbind(x, y+1),
                                cbind(x, y-1),
                                cbind(x+1, y),
                                cbind(x-1, y),
                                cbind(x+1, y+1),
                                cbind(x+1, y-1),
                                cbind(x-1, y+1),
                                cbind(x-1, y-1))
    }
    ExternalVM <- rbind(ExternalVM, additional_cells)
    ExternalVM <- ExternalVM[!duplicated(paste0(ExternalVM[,1], '_', ExternalVM[,2])),]
  }
  
  ExternalVM <- as.data.frame(ExternalVM)
  ExternalVM$R <- rep(col2rgb(color)[1], nrow(ExternalVM))
  ExternalVM$G <- rep(col2rgb(color)[2], nrow(ExternalVM))
  ExternalVM$B <- rep(col2rgb(color)[3], nrow(ExternalVM))
  ExternalVM$color <- rep(rgb(col2rgb(color)[1], col2rgb(color)[2], col2rgb(color)[3], maxColorValue=255), nrow(ExternalVM))
  ExternalVM$cluster <- rep(length(unique(VM$cluster))+1, nrow(ExternalVM))
  ExternalVM$cluster_refined <- rep(length(unique(VM$cluster))+1, nrow(ExternalVM))
  ExternalVM$cluster_color <- rep(color, nrow(ExternalVM))
  ExternalVM$cluster_refined_color <- rep(color, nrow(ExternalVM))
  ExternalVM$cluster_annot <- rep(cluster_annot, nrow(ExternalVM))
  VM <- rbind(VM, ExternalVM)
  rownames(VM) <- paste0(as.vector(unlist(VM[,1])), '_', as.vector(unlist(VM[,2])))
  return(as.data.frame(VM))
}


#' selectLandmark
#'
#' Select landmark on the virtual map
#' @param VM Virtual map data frame
#' @param reference_group Group in which landmark will be located. The reference group must be included in cluster_annot
#' @param type Whether the landmark will be the central position on the cluster ('centroid'), the
#' left-most cells ('vertical_line') or the upper-most cells ('horizontal-line')
#' @param landmark_name Name assigned to the landmark
#' 
#' @return A virtual map data.frame containing an additional column called is.landmark with the 
#' corresponding values
#'
#' @examples
#' VM <- selectLandmark(VM, reference_group='Antenna_A3_Arista', type='centroid', landmark_name='Antenna')
#' VM <- selectLandmark(VM, reference_group='MF_Morphogenetic_Furrow', type='vertical_line', landmark_name='Eye')
#'
#' @export
selectLandmark <- function(VM, 
                           reference_group,
                           type,
                           landmark_name){
  if (!reference_group %in% VM$cluster_annot){
    stop('The reference_group must be a spatial domain specified in cluster_annot') 
  }
  
  ref_group_coor <- VM[which(VM$cluster_annot == reference_group),]
  if (type=='centroid'){
    landmark <- paste0(as.integer(mean(as.vector(unlist(ref_group_coor[,1])))), '_', as.integer(mean(as.vector(unlist(ref_group_coor[,2])))))
  } 
  else if (type=='vertical_line'){
    y_coord <- as.vector(unlist(unique(ref_group_coor[,2])))
    x_coord <- sapply(1:length(y_coord), function(i) min(as.vector(unlist(ref_group_coor[which(ref_group_coor[,'y'] == y_coord[i]), 'x']))))
    landmark <- paste0(x_coord, '_', y_coord)
  } else if (type=='horizontal_line'){
    x_coord <- as.vector(unlist(unique(ref_group_coor[,1])))
    y_coord <- sapply(1:length(x_coord), function(i) max(as.vector(unlist(ref_group_coor[which(ref_group_coor[,'x'] == x_coord[i]), 'y']))))
    landmark <- paste0(x_coord, '_', y_coord)
  }
  
  if(!'is.landmark' %in% colnames(VM)){
    VM$is.landmark <- rep(NA, nrow(VM))
  }
  
  if (landmark %in% rownames(VM)){
    VM$is.landmark[which(rownames(VM) %in% landmark)] <- landmark_name
  }
  else {
    ref_group_coor$is.landmark <- rep(NA, nrow(ref_group_coor))
    VM <- rbind(VM, ref_group_coor[1,,drop=FALSE])
    rownames(VM)[nrow(VM)] <- landmark
    VM[landmark, 'x'] <- strsplit(landmark, '_')[[1]][1]
    VM[landmark, 'y'] <- strsplit(landmark, '_')[[1]][2]
    VM$is.landmark[which(rownames(VM) %in% landmark)] <- landmark_name
  }
  
  return(as.data.frame(VM))
}
