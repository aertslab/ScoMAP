#' getSearchSpace
#'
#' Get regions per gene in which to look for potential enhancers to be linked
#' @param txdb Txdb object matching with the genome assembly used for the analysis
#' @param org.db Org.Db objet for the corresponding species
#' @param genes Genes for which enhancer-to-gene links want to be inferred
#' @param extend Space around the TSS that must be considered (first value, upstream; second, downstream TSS).
#' In addition, intronic regions in the gene will be considered. Default=c(50000, 50000)
#'
#' @return Genomic ranges object containing the regions considered for each gene (as metadata)
#'
#' @examples
#' searchSpace <- getSearchSpace(txdb, org.db, rownames(DGEM), extend=c(50000, 50000))
#'
#' @import AnnotationDbi
#' @import GenomicRanges
#'
#' @export

getSearchSpace <- function(txdb, org.db, genes, extend=c(50000, 50000)) {
  # Check up
  if(!is(txdb,'TxDb')){
    stop('txdb has to be an object of class TxDb')
  }
  if(!is(org.db,'OrgDb')){
    stop('org.db has to be an object of class OrgDb')
  }
  # Get search space around TSS
  # Genes to ensemble dict
  ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = genes, columns="ENSEMBL", keytype="SYMBOL")
  if (sum(is.na(ENS2SYMBOL[,2])) > 0){ENS2SYMBOL <- ENS2SYMBOL[-which(is.na(ENS2SYMBOL[,2])),]}
  # Select genes in the list
  filter_list <- list()
  filter_list[['GENEID']] <- ENS2SYMBOL[,2]
  # Take promoter coordinates for the specific genes
  TSS <- promoters(txdb, upstream =extend[1], downstream=extend[2], filter=filter_list, columns=c("GENEID"))
  # Annotate to symbol
  ENS2SYMBOL_VECTOR <- as.vector(ENS2SYMBOL[,1])
  names(ENS2SYMBOL_VECTOR) <- ENS2SYMBOL[,2]
  elementMetadata(TSS)$SYMBOL <- ENS2SYMBOL_VECTOR[unlist(as.vector(elementMetadata(TSS)$GENEID))]
  elementMetadata(TSS) <- elementMetadata(TSS)[ , -which(colnames(elementMetadata(TSS)) %in% c('GENEID', 'width'))]
  colnames(elementMetadata(TSS)) <- 'SYMBOL'
  elementMetadata(TSS)$RegionType <- rep('Extended Promoter', length(unlist(as.vector(elementMetadata(TSS)$SYMBOL))))
  
  
  # Get introns
  introns <- intronicParts(txdb, linked.to.single.gene.only=TRUE)
  # Get Flybase to symbol name
  FB2SYMBOL <- select(org.Dm.eg.db, keys = unlist(as.vector(elementMetadata(introns)$gene_id)), columns="SYMBOL", keytype="FLYBASE")
  if (sum(is.na(FB2SYMBOL[,1])) > 0){FB2SYMBOL <- FB2SYMBOL[-which(is.na(FB2SYMBOL[,1])),]}
  FB2SYMBOL_VECTOR <- as.vector(FB2SYMBOL[,2])
  names(FB2SYMBOL_VECTOR) <- FB2SYMBOL[,1]
  elementMetadata(introns)$SYMBOL <- FB2SYMBOL_VECTOR[unlist(as.vector(elementMetadata(introns)$gene_id))]
  elementMetadata(introns) <- elementMetadata(introns)[ , -which(colnames(elementMetadata(introns)) %in% c('gene_id', 'tx_name', 'tx_id'))]
  colnames(elementMetadata(introns)) <- 'SYMBOL'
  # Subset for selected genes
  introns <- introns[which(elementMetadata(introns)$SYMBOL %in% genes),]
  elementMetadata(introns)$RegionType <- rep('Intron', length(unlist(as.vector(elementMetadata(introns)$SYMBOL))))
  
  # Merge regions
  searchSpace <- c(TSS, introns)
  return(searchSpace)
}

#' enhancerToGene
#'
#' Link enhancers to target genes
#' @param VM_RNA_mat Matrix containing genes as rows, virtual cells as columns and gene expression as values (see getVirtualFeatureMatrix())
#' @param VM_RNA_mat Matrix containing regions as rows, virtual cells as columns and cisTopic's region accessibility probabilities as values (recommended, see getVirtualFeatureMatrix())
#' @param searchSpace Search space GenomicRanges object as obtained from getSearchSpace()
#' @param method Whether to use pearson correlation between region accessibility and gene expression ('Correlation')
#' or a random forest ('RF') model to infer enhancer-to-gene links.
#' @param minoverlap Minimum overlap between the candidate regulatory regions and the search space. Default: 0.4.
#' @param nCores How many cores to use if method='RF'. Default: 1
#' @param nTrees How many trees to use per RF model if method='RF'. Default: 1000
#'
#' @return A list containing a data frame for each gene, in which the RF importance (if method='RF') or the
#' correlation values (if method='Correlation') are given.
#'
#' @examples
#' RF_links <- enhancerToGene(VM_DGEM, VM_PRTM, searchSpace, method='RF') 
#' Cor_links <- enhancerToGene(VM_DGEM, VM_PRTM, searchSpace, method='Correlation') 
#'
#' @import parallel
#' @import doSNOW
#' @import Matrix
#' @importFrom plyr llply
#'
#' @export

enhancerToGene <- function(VM_RNA_mat,
                           VM_ATAC_mat,
                           searchSpace,
                           method,
                           minOverlap=0.4,
                           nCores = 1,
                           nTrees=1000
){
  #Get region Granges
  regions <- rownames(VM_ATAC_mat)
  GR_regions <- .regionVectorToGenomicRanges(regions)
  region2gene <- .getOverlapRegionsFromGR_regions(GR_regions, searchSpace)
  region2gene_split <- split(region2gene, region2gene$gene)
  gene_names <- names(region2gene_split) 
  region2gene_list <-  llply(1:length(region2gene_split), function(i) as.vector(unlist(region2gene_split[[i]][,1])))
  names(region2gene_list) <- gene_names 
  
  region2gene_list <- region2gene_list[names(region2gene_list) %in% rownames(VM_RNA_mat)]
  gene_names<- names(region2gene_list)
  region2gene_list <-  llply(1:length(region2gene_list), function(i) region2gene_list[[i]][region2gene_list[[i]] %in% rownames(VM_ATAC_mat)])
  names(region2gene_list) <- gene_names 
  region2gene_list <- llply(1:length(region2gene_list), function (i) rbind(VM_RNA_mat[(names(region2gene_list)[i]),,drop=FALSE], VM_ATAC_mat[region2gene_list[[i]],]))
  names(region2gene_list) <- gene_names 
  
  if(method == 'RF'){
    if (nCores > 1){
      cl <- makeCluster(nCores, type = "SOCK")
      registerDoSNOW(cl)
      clusterEvalQ(cl, library(GENIE3))
      clusterExport(cl, c('region2gene_list'), envir=environment())
      opts <- list(preschedule=TRUE)
      clusterSetRNGStream(cl, 123)
      output <- llply(1:length(region2gene_list), function(i) as.data.frame(GENIE3(as.matrix(region2gene_list[[i]]), treeMethod = "RF", K = "sqrt", nTrees = nTrees,
                                                                regulators = rownames(region2gene_list[[i]])[-1], targets = rownames(region2gene_list[[i]])[1],
                                                                nCores = nCores, verbose = FALSE) , .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
    } else {
      output <- llply(1:length(region2gene_list), function(i) as.data.frame(GENIE3(as.matrix(region2gene_list[[i]]), treeMethod = "RF", K = "sqrt", nTrees = nTrees,
                                                                     regulators = rownames(region2gene_list[[i]])[-1], targets = rownames(region2gene_list[[i]])[1],
                                                                     nCores = nCores, verbose = FALSE) , .parallel = FALSE, .paropts = list(.options.snow=opts), .inform=FALSE))
    }
    output <- lapply(1:length(output), function (i) {colnames(output[[i]]) <- 'RF_importance'; output[[i]]})
    names(output) <- gene_names
  } else if (method == 'Correlation'){
    output <- lapply(1:length(region2gene_list), function(i) as.data.frame(apply(region2gene_list[[i]][2:nrow(region2gene_list[[i]]),] , 1 , cor , y = region2gene_list[[i]][1,])))
    output <- lapply(1:length(output), function (i) {colnames(output[[i]]) <- 'Correlation'; output[[i]]})
    names(output) <- names(region2gene_list)
  }
  return(output)
}
  
.regionVectorToGenomicRanges <- function(regionVector){
  chr <-  sapply(strsplit(regionVector, split = ":"), "[", 1)
  coord <-  sapply(strsplit(regionVector, split = ":"), "[", 2)
  start <- as.numeric(sapply(strsplit(coord, split = "-"), "[", 1))
  end <- as.numeric(sapply(strsplit(coord, split = "-"), "[", 2))
  dataframe <- as.data.frame(cbind(chr, start, end))
  colnames(dataframe) <- c('seqnames', 'start', 'end')
  rownames(dataframe) <- regionVector
  Gr <- makeGRangesFromDataFrame(dataframe, keep.extra.columns = TRUE) 
  return(Gr)
}

.getOverlapRegionsFromGR_regions <- function(
  GR_regions,
  searchSpace,
  minOverlap=0.4,
  overlapping=TRUE,
  ...)
{
  dbRegionsOverlap <- findOverlaps(searchSpace, GR_regions, type='any', select="all", ignore.strand=TRUE)
  
  if(minOverlap>0){
    overlaps <- pintersect(searchSpace[queryHits(dbRegionsOverlap)], GR_regions[subjectHits(dbRegionsOverlap)])
    percentOverlapGR_regions <- width(overlaps) / width(searchSpace[queryHits(dbRegionsOverlap)])
    percentOverlapRegions <- width(overlaps) / width(GR_regions[subjectHits(dbRegionsOverlap)])
    maxOverlap <- apply(cbind(percentOverlapGR_regions, percentOverlapRegions), 1, max)
    dbRegionsOverlap <- dbRegionsOverlap[which(maxOverlap > minOverlap)]
    maxOverlap <- maxOverlap[which(maxOverlap > minOverlap)]
  }
  
  selectedRegions <- searchSpace[queryHits(dbRegionsOverlap)]
  symbol <- as.data.frame(selectedRegions)$SYMBOL
  selectedRegions <- paste(as.vector(seqnames(selectedRegions)), ':', as.vector(start(selectedRegions)), '-', as.vector(end(selectedRegions)), sep='')
  selectedGR_regions <- names(GR_regions[subjectHits(dbRegionsOverlap)])
  
  selectedMapped <- data.frame(selectedGR_regions, selectedRegions, symbol, maxOverlap, row.names=NULL)
  colnames(selectedMapped) <- c('dataRegion', 'searchSpace', 'gene', 'maxOverlap')
  
  return(selectedMapped)
}

#' plotLinks
#'
#' Plot enhancer-to-gene links (To be enhanced)
#' @param links Data frames list containing correlation or RF scores for each region in each gene as returned by enhancerToGene().
#' @param annot Annotation data frame, as required by cicero (Pliner et al., 2019)
#' @param txdb Txdb object matching with the genome assembly used for the analysis
#' @param org.db Org.Db objet for the corresponding species
#' @param gene Gene for which enhancer-to-gene links wants to be plotted
#' @param chr Chromosome name of the genomic window that wants to be plotted
#' @param start Start position of the genomic window that wants to be plotted
#' @param end End position of the genomic window that wants to be plotted
#' @param cutoff Value below which links will not be shown. Default: -1.
#'
#' @examples
#' plotLinks(RF_links, dm6_annot, TxDb.Dmelanogaster.UCSC.dm6.ensGene, org.Dm.eg.db, gene='dac', 'chr2L', 16470000, 16490000)
#'
#' @import AnnotationDbi
#' @import GenomicRanges
#'
#' @export


plotLinks <- function(links,
                      annot,
                      txdb,
                      org.db,
                      gene, 
                      chr,
                      start,
                      end,
                      cutoff=-1){
  # Check up
  if(! "cicero" %in% installed.packages()){
    stop('Please, install cicero: \n BiocManager::install("cicero")')
  } else {
    require(cicero)
  }
  
  # Check up
  if(!is(txdb,'TxDb')){
    stop('txdb has to be an object of class TxDb')
  }
  if(!is(org.db,'OrgDb')){
    stop('org.db has to be an object of class OrgDb')
  }
  # Get search space around TSS
  # Genes to ensemble dict
  ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = gene, columns="ENSEMBL", keytype="SYMBOL")
  if (sum(is.na(ENS2SYMBOL[,2])) > 0){ENS2SYMBOL <- ENS2SYMBOL[-which(is.na(ENS2SYMBOL[,2])),]}
  # Select genes in the list
  filter_list <- list()
  filter_list[['GENEID']] <- ENS2SYMBOL[,2]
  # Take promoter coordinates for the specific genes
  TSS <- promoters(txdb, upstream = 0, downstream= 0, filter=filter_list, columns=c("GENEID"))
  TSS <- paste(seqnames(TSS), start(TSS), end(TSS), sep='_')
  # Form conns data frame
  linksgene <- links[[gene]]
  regions <- gsub(':', '_', rownames(linksgene))
  regions <- gsub('-', '_', rownames(linksgene))
  conns <- as.data.frame(cbind(regions, rep(TSS, length(regions)), linksgene[,1], rep('black', length(regions))))
  names(conns) <-  c('Peak1', 'Peak2', 'coaccess', 'color')
  conns[,3] <- as.numeric(as.vector(unlist(conns[,3])))
  cicero::plot_connections(conns, chr, start, end,
                   alpha_by_coaccess = TRUE,
                   gene_model = annot, 
                   connection_color='color',
                   connection_color_legend=F,
                   gene_model_color='blue',
                   coaccess_cutoff = cutoff, 
                   connection_width = .5, 
                   collapseTranscripts = "longest",
                   peak_color = "black")
}