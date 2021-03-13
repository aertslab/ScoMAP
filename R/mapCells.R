#' mapCells
#'
#' Maps cells on a real cells into a virtual map.
#' @param VM Virtual map data frame
#' @param RM Real map data frame. It must contain the following columns:
#' - 'Cell_type': Annotated cell type
#' - 'Landmark': Reference landmark to use. It must be named as the landmarks in the virtual map, or 
#' 'None' is the cells will be randomly mapped on the matching spatial cluster.
#' - 'Pseudotime' : Pseudotime value for the cell. Landmark cell type must be included in the trajectory.
#' - 'PosLandmark' : Relative position to the landmark. If it is set to'Before', it means that the landmark is later in pseudotime,
#' if it is set to 'After' it means that the landmark is located earlier in pseudotime.
#'-  'Spatial_cluster': Spatial cluster to which it will be mapped. These clusters must aggree with VM$cluster_annot.
#' @param target_clusters Spatial clusters to map. If NULL, all clusters will be mapped. Default: NULL.
#' @param nr_bin Number of bins in which the pseudotime and spatial map (as distance from landmark) ordering are divided for mapping.
#' Default: 10. 
#' @param trajectory_list A list with dataframes with Pseudotime and Poslandmark columns inferred from spatial data
#' @param seed Seed (as assignment is random). Default: 123.
#'
#' @return A virtual map data.frame containing the assigned real cell and its information.
#'
#' @examples
#' VM_RNA <- mapCells(VM_RNA, RM_RNA, nr_bin=10)
#'
#' @export


mapCells <- function(VM,
                     RM,
                     target_clusters=NULL,
                     nr_bin=10,
                     seed=123,
                     trajectory_list=list()){
  if(is.null(target_clusters)){
    target_clusters <- unique(RM$Spatial_cluster)[!is.na(unique(RM$Spatial_cluster))]
  }
  for (target_cluster in target_clusters){
    VM <- .mapCells_int(VM, RM, target_cluster, nr_bin, seed, trajectory_list)
  }
  return(VM)
}

# Helper function
.mapCells_int <- function(VM,
         RM,
         target_cluster,
         nr_bin=10,
         seed=123,
         trajectory_list=list()){
  # Check
  if(!target_cluster %in% RM$Spatial_cluster){
    stop('The specified cluster is not defined in the real map (`RM$Spatial_cluster`).')
  }
  
  if(!target_cluster %in% VM$cluster_annot){
    stop('The specified cluster is not defined in the virtual map (`VM$cluster_annot`).')
  }
  
  # Set seed (random sampling)
  set.seed(seed)
  # Prepare data for target type
  target_RM <- RM[which(RM$Spatial_cluster == target_cluster),]
  target_VM <- VM[which(VM$cluster_annot==target_cluster),]
  # Take landmark
  target_landmark_RM <- unique(target_RM$Landmark)
  
  if (target_landmark_RM == 'None'){
    print(paste0('The target cluster ', target_cluster, ' will be ramdomly mapped.'))
    if (nrow(target_RM) > nrow(target_VM)){
      cell_assigment <- rownames(target_RM)[sample(nrow(target_VM))]
    } else {
      if (nrow(target_VM)-nrow(target_RM) > nrow(target_RM)){
        cell_assigment <- sample(c(sample(rownames(target_RM), nrow(target_RM)), sample(rownames(target_RM), nrow(target_VM)-nrow(target_RM), replace=TRUE)))
      } else {
        cell_assigment <- sample(c(sample(rownames(target_RM), nrow(target_RM)), sample(rownames(target_RM), nrow(target_VM)-nrow(target_RM))))
      }
    }
    names(cell_assigment) <- rownames(target_VM)
    target_RM$PseudotimeRank <- rep('0', nrow(target_RM))
    distance2landmark <- rep('0', nrow(target_VM))
    names(distance2landmark) <- rownames(target_VM)
  } else {
    if (target_landmark_RM == 'Trajectory'){
      print(paste0('The target cluster ', target_cluster, ' will be mapped to a trajectory.'))
      distance2landmark <- as.vector(unlist(trajectory_list[[target_cluster]][rownames(target_VM),'Pseudotime']))
    } else {
      print(paste0('The target cluster ', target_cluster, ' will be mapped to a spatial axis.'))
      # Change names now if not mapping to trajectory
      rownames(target_VM) <- paste0(as.vector(unlist(target_VM[,'x'])), '_', as.vector(unlist(target_VM[,'y'])))
      # Check
      if(!target_landmark_RM %in% VM$is.landmark){
        stop('The landmark specified in the real map (`RM$Landmark`) does not exist in the virtual map (`VM$is.landmark`)')
      }
      # Calculate distance of virtual cell to landmark
      landmark <- VM[which(VM$is.landmark == target_landmark_RM), c('x', 'y')]
      if (nrow(landmark) == 1){
        distance2landmark <- sapply(1:nrow(target_VM), function (i) sqrt(sum((as.numeric(target_VM[i,c('x','y')])-as.numeric(landmark))^2)))
      } else if (nrow(landmark) > 1){
        distance2landmark <- as.vector(unlist(lapply(1:nrow(target_VM), function (i) min(sqrt(rowSums((rbind(target_VM[i,c('x','y')][rep(1, nrow(landmark)), ])-landmark)^2))))))
      }
    }
    names(distance2landmark) <- rownames(target_VM)
    distance2landmark <- distance2landmark[order(distance2landmark)]
    VM_bin <- split(distance2landmark, factor(sort(rank(distance2landmark)%%nr_bin)))
    # Make bins on pseudotime
    pseudotime_order <- as.numeric(target_RM$Pseudotime)
    names(pseudotime_order) <- rownames(target_RM)
    if (target_landmark_RM == 'Trajectory'){
      if (unique(as.vector(unlist(target_RM$PosLandmark))) != unique(as.vector(unlist(trajectory_list[[target_cluster]]$PosLandmark)))){
        print(paste('Direction spatial trajectory:', unique(as.vector(unlist(target_RM$PosLandmark)))))
        print(paste('Direction omics trajectory:', unique(as.vector(unlist(trajectory_list[[target_cluster]]$PosLandmark)))))
        pseudotime_order <- -pseudotime_order
      }
    } else {
      if (unique(as.vector(unlist(target_RM$PosLandmark))) == 'Before'){
        pseudotime_order <- -pseudotime_order
      }
    }
    pseudotime_order <- rank(pseudotime_order)
    pseudotime_order <- pseudotime_order[order(pseudotime_order)]
    target_RM[names(pseudotime_order), 'PseudotimeRank'] <- pseudotime_order
    RM_bin <- split(pseudotime_order, factor(sort(rank(pseudotime_order)%%nr_bin)))
    denom_VM <- max(lengths(VM_bin))
    denom_RM <- max(lengths(RM_bin))
    # Map cells
    if (denom_VM < denom_RM){
      if (sum(lengths(RM_bin) > denom_VM) >= nr_bin){
        cell_assigment <- as.vector(unlist(lapply(1:length(VM_bin), function(i) sample(names(RM_bin[[i]]), length(VM_bin[[i]])))))
      } else{
        print(length(RM_bin))
        print(length(VM_bin))
        cell_assigment <- as.vector(unlist(lapply(1:length(VM_bin), function(i) sample(names(RM_bin[[i]]), length(VM_bin[[i]]), replace=TRUE))))
      }
    } else {
      cell_assigment <- as.vector(unlist(sapply(1:length(VM_bin), function(i)  sample(c(names(RM_bin[[i]]), sample(names(RM_bin[[i]]), length(VM_bin[[i]])-length(RM_bin[[i]]), replace=TRUE))))))
    } 
    names(cell_assigment) <- as.vector(unlist(lapply(VM_bin, names)))
  }
  
  if(!'RM_assignment' %in% colnames(VM)){
    VM$RM_assignment <- rep(NA, nrow(VM))
    VM$RM_cell_type <- rep(NA, nrow(VM))
    VM$RM_landmark <- rep(NA, nrow(VM))
    VM$RM_pseudotime <- rep(0, nrow(VM))
    VM$RM_posLandmark <- rep(NA, nrow(VM))
    VM$RM_spatial_cluster <- rep(NA, nrow(VM))
    VM$DistanceFromLandmark- rep(0, nrow(VM))
    VM$PseudotimeRank- rep(0, nrow(VM))
  }
  VM[names(cell_assigment), 'RM_assignment'] <- cell_assigment
  VM[names(cell_assigment), 'RM_cell_type'] <- as.vector(unlist(target_RM[cell_assigment,'Cell_type']))
  VM[names(cell_assigment), 'RM_landmark'] <- as.vector(unlist(target_RM[cell_assigment,'Landmark']))
  VM[names(cell_assigment), 'RM_pseudotime'] <- as.vector(unlist(target_RM[cell_assigment,'Pseudotime']))
  VM[names(cell_assigment), 'RM_posLandmark'] <- as.vector(unlist(target_RM[cell_assigment,'PosLandmark']))
  VM[names(cell_assigment), 'RM_spatial_cluster'] <- as.vector(unlist(target_RM[cell_assigment,'Spatial_cluster']))
  VM[names(distance2landmark), 'DistanceFromLandmark'] <- distance2landmark
  VM[names(cell_assigment), 'PseudotimeRank'] <- as.vector(unlist(target_RM[cell_assigment,'PseudotimeRank']))
  
  return(as.data.frame(VM))
}

#' getVirtualFeatureMatrix
#'
#' Get the feature values on the virtual cells
#' @param VM Virtual map data frame, after running mapCells()
#' @param omics_mat Matrix containing features as rows (e.g. genes, regions) and real cells as columns.
#'
#' @return Matrix containing virtual cells as columns and features as rows
#'
#' @examples
#' 
#' VM_DGEM <- getVirtualFeatureMatrix(VM_RNA, DGEM)
#' 
#' @export
getVirtualFeatureMatrix <- function(VM,
                                    omics_mat){
  if (!'RM_assignment' %in% colnames(VM)){
    stop('Please, run mapCells() first.')
  }
  if(sum(VM$RM_assignment %in% colnames(omics_mat)) != nrow(VM)){
    stop('Are the cells asisgned to this virtual map on the omics feature matrix?')
  }
  # Get profiles from the real cells mapped
  omics_mat <- omics_mat[,VM$RM_assignment]
  colnames(omics_mat) <- rownames(VM)
  return(omics_mat)
}
