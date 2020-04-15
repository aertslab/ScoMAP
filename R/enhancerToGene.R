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
  if (taxonomyId(org.db) == 9606){
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = genes, columns="ENTREZID", keytype="SYMBOL")
  } else {
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = genes, columns="ENSEMBL", keytype="SYMBOL")
  }
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
  introns <- intronicParts(txdb, linked.to.single.gene.only=FALSE)
  # Get Ens to symbol name
  if (taxonomyId(org.db) == 9606){
    ENS2SYMBOLL <- AnnotationDbi::select(org.db, keys = unlist(as.vector(elementMetadata(introns)$gene_id)), columns="SYMBOL", keytype="ENTREZID")
  } else {
    ENS2SYMBOLL <- AnnotationDbi::select(org.db, keys = unlist(as.vector(elementMetadata(introns)$gene_id)), columns="SYMBOL", keytype="ENSEMBL")
  }
  if (sum(is.na(ENS2SYMBOLL[,1])) > 0){ENS2SYMBOLL <- ENS2SYMBOLL[-which(is.na(ENS2SYMBOLL[,1])),]}
  ENS2SYMBOLL_VECTOR <- as.vector(ENS2SYMBOLL[,2])
  names(ENS2SYMBOLL_VECTOR) <- ENS2SYMBOLL[,1]
  introns <- S4Vectors::expand(introns, "gene_id")
  elementMetadata(introns)$SYMBOL <- ENS2SYMBOLL_VECTOR[unlist(as.vector(elementMetadata(introns)$gene_id))]
  elementMetadata(introns) <- elementMetadata(introns)[ , -which(colnames(elementMetadata(introns)) %in% c('gene_id', 'tx_name', 'tx_id'))]
  colnames(elementMetadata(introns)) <- 'SYMBOL'
  # Subset for selected genes
  introns <- introns[which(elementMetadata(introns)$SYMBOL %in% genes),]
  elementMetadata(introns)$RegionType <- rep('Intron', length(unlist(as.vector(elementMetadata(introns)$SYMBOL))))
  
  # Merge regions
  searchSpace <- c(TSS, introns)
  names(searchSpace) <- NULL
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
  region2gene_list <-  llply(region2gene_split, function(x) unique(as.vector(unlist(x[,1]))))
  
  region2gene_list <- region2gene_list[names(region2gene_list) %in% rownames(VM_RNA_mat)]
  region2gene_list <-  llply(region2gene_list, function(x) x[x %in% rownames(VM_ATAC_mat)])
  gene_names <- names(region2gene_list) 
  region2gene_list <- llply(1:length(region2gene_list), function (i) rbind(VM_RNA_mat[(names(region2gene_list)[i]),,drop=FALSE], VM_ATAC_mat[region2gene_list[[i]],]))
  names(region2gene_list) <- gene_names
  
  if(method == 'RF'){
    if(! "GENIE3" %in% installed.packages()){
      stop('Please, install GENIE3: \n BiocManager::install("GENIE3")')
    } else {
      require(GENIE3)
    }
    if (nCores > 1){
      cl <- makeCluster(nCores, type = "SOCK")
      registerDoSNOW(cl)
      clusterEvalQ(cl, library(GENIE3))
      clusterExport(cl, c('region2gene_list'), envir=environment())
      opts <- list(preschedule=TRUE)
      clusterSetRNGStream(cl, 123)
      output <- llply(region2gene_list, function(x) as.data.frame(tryCatch(GENIE3(as.matrix(x), treeMethod = "RF", K = "sqrt", nTrees = nTrees,
                                                                regulators = rownames(x)[-1], targets = rownames(x)[1],
          nCores = 1, verbose = FALSE), error=function(e) NULL), .parallel = TRUE, .paropts = list(.options.snow=opts), .inform=FALSE))
    } else {
        output <- llply(region2gene_list, function(x) as.data.frame(tryCatch(GENIE3(as.matrix(x), treeMethod = "RF", K = "sqrt", nTrees = nTrees,
                                                                     regulators = rownames(x)[-1], targets = rownames(x)[1],
            nCores = 1, verbose = FALSE), error=function(e) NULL), .parallel = FALSE, .inform=FALSE, .progress = "text"))
    }
    if (sum(sapply(output, is.null)) > 0){
        output <- output[-which(sapply(output, is.null))]
    }
    if (0 %in% sapply(output, nrow)){
      output <- output[-which(sapply(output, nrow) == 0)]
    }
    output <- lapply(output, function (x) {colnames(x) <- 'RF_importance'; x})
  } else if (method == 'Correlation'){
    output <- lapply(region2gene_list, function (x) t(as.data.frame(t(apply(x[2:nrow(x),,drop=FALSE], 1 , cor , y = x[1,])))))
    output <- lapply(output, function (x) {colnames(x) <- 'Correlation'; x})
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
#' @param RF_links Data frames list containing RF scores for each region in each gene as returned by enhancerToGene().
#' @param Cor_links Data frames list containing correlation scores for each region in each gene as returned by enhancerToGene(). If both RF and correlation links are provided, the height of the links will represent the RF importance and the color whether the correlation is positive or negative. If only RF is provided, links will be colored black; and if only correlation links are provided the height of the lnk will indicate the absolute correlation value and the color whether it is positive or negative.
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


plotLinks <- function(RF_links=NULL,
                      Cor_links=NULL,
                      annot,
                      txdb,
                      org.db,
                      gene, 
                      chr,
                      start,
                      end,
                      cutoff=0){
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
  if (taxonomyId(org.db) == 9606){
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = gene, columns="ENTREZID", keytype="SYMBOL")
  } else {
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = gene, columns="ENSEMBL", keytype="SYMBOL")
  }
  if (sum(is.na(ENS2SYMBOL[,2])) > 0){ENS2SYMBOL <- ENS2SYMBOL[-which(is.na(ENS2SYMBOL[,2])),]}
  # Select genes in the list
  filter_list <- list()
  filter_list[['GENEID']] <- ENS2SYMBOL[,2]
  # Take promoter coordinates for the specific genes
  TSS <- promoters(txdb, upstream = 0, downstream= 0, filter=filter_list, columns=c("GENEID"))
  TSS <- paste(seqnames(TSS), start(TSS), end(TSS), sep='_')
  
  if (!is.null(RF_links)){
    # Form conns data frame
    linksgene <- RF_links[[gene]]
    regions <- gsub(':', '_', rownames(linksgene))
    regions <- gsub('-', '_', rownames(linksgene))
    conns <- as.data.frame(cbind(regions, rep(TSS, length(regions)), linksgene[,1], rep('black', length(regions))))
    if(!is.null(Cor_links)){
      corlinksgene <- as.vector(unlist(Cor_links[[gene]]))
      mypal <- colorRampPalette(c('#FF0000', '#228B22'))(10)
      color <- .map2color(corlinksgene,mypal)
      conns <- as.data.frame(cbind(regions, rep(TSS, length(regions)), linksgene[,1], color))
    }
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
  } else if (!is.null(Cor_links)) {
    # Form conns data frame
    linksgene <- Cor_links[[gene]]
    regions <- gsub(':', '_', rownames(linksgene))
    regions <- gsub('-', '_', rownames(linksgene))
    mypal <- colorRampPalette(c('#FF0000', '#228B22'))(10)
    color <- .map2color(linksgene,mypal)
    linksgene <- abs(linksgene)
    conns <- as.data.frame(cbind(regions, rep(TSS, length(regions)), linksgene[,1], color))
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

}

#' pruneLinks
#'
#' Prune enhancer-to-gene links 
#' @param RF_links Data frames list containing RF scores for each region in each gene as returned by enhancerToGene().
#' @param Cor_links Data frames list containing correlation scores for each region in each gene as returned by enhancerToGene(). If both RF and correlation links are provided, the height of the links will represent the RF importance and the color whether the correlation is positive or negative. If only RF is provided, links will be colored black; and if only correlation links are provided the height of the lnk will indicate the absolute correlation value and the color whether it is positive or negative.
#' @param cor_prob Probability threshold on the fitted distribution above which positive and negative correlation will be taken.
#' @param cor_thr A vector containing the lower and upper thresholds for the correlation links. 
#'
#' @return A list containing slots with the pruned RF links ('RF_links') and/or correlation links ('Cor_links')
#' @examples
#' pruneLinks(RF_links, Cor_links)
#'
#' @import plyr
#'
#' @export


pruneLinks <- function(RF_links=NULL,
                      Cor_links=NULL,
                      cor_prob=0.01,
                      cor_thr=NULL){
  
  if(is.null(RF_links) & is.null(Cor_links)){
    stop('Please, provide at least either RF or correlation links.')
  }
  
  if (!is.null(RF_links)){
    if(! "Binarize" %in% installed.packages()){
      stop('Please, install Binarize: \n install.packages("Binarize")')
    } else {
      require(Binarize)
    }
    RF_links_tmp <- llply(RF_links, as.matrix)
    RF_links_tmp <- llply(RF_links_tmp, function (x)  x[complete.cases(x),])
    RF_links_2_enh <- RF_links_tmp[lengths(RF_links_tmp) <= 2]
    RF_links_tmp <- RF_links_tmp[lengths(RF_links_tmp) > 2]
    RF_links_tmp <- llply(RF_links_tmp, function (x) x[which(x >= binarize.BASC(x)@threshold)])
    if(length(RF_links_2_enh) > 0){
      RF_links_tmp <- c(RF_links_tmp, RF_links_2_enh)
    }
    RF_links_tmp <- llply(RF_links_tmp, as.data.frame)
    RF_links_tmp <- llply(RF_links_tmp, function (x) {colnames(x) <- 'RF_importance'; x}) 
  } 
  
  if(!is.null(Cor_links)){
    if(! "fitdistrplus" %in% installed.packages()){
      stop('Please, install fitdistrplus: \n install.packages("fitdistrplus")')
    } else {
      require(fitdistrplus)
    }
    Cor_links_tmp <- llply(Cor_links, as.matrix)
    fit <- suppressWarnings(fitdistrplus::fitdist(unlist(Cor_links_tmp)[!is.na(unlist(Cor_links_tmp))], "norm", method='mme')) 
    cutofflow <- as.numeric(unlist(quantile(fit, probs = cor_prob))[1])
    cutoffup <- as.numeric(unlist(quantile(fit, probs = 1-cor_prob))[1])
    if (!is.null(cor_thr)){
      cutofflow <- cor_thr[1] 
      cutoffup <- cor_thr[2] 
    }
    print(paste0('Low cutoff: ', cutofflow, '; Upper cutoff:', cutoffup))
    Cor_links_tmp <- llply(Cor_links_tmp, function (x) x[which(x > cutoffup | x < cutofflow),])
    Cor_links_tmp <- llply(Cor_links_tmp, as.data.frame)
    Cor_links_tmp <- llply(Cor_links_tmp, function (x) {colnames(x) <- 'Correlation'; x})
  }
  
  if (!is.null(RF_links) & !is.null(Cor_links)){
    # Check-up
    RF_links_tmp <- RF_links_tmp[names(RF_links_tmp) %in% names(Cor_links_tmp)]
    Cor_links_tmp <- Cor_links_tmp[names(RF_links_tmp)]
    RF_thr_enh <- llply(RF_links_tmp, rownames)
    Cor_links_enh <- llply(Cor_links_tmp, rownames)
    thr_enh <- llply(1:length(RF_thr_enh), function(i) c(RF_thr_enh[[i]], Cor_links_enh[[i]]))
    RF_links <- RF_links[names(RF_links_tmp)]
    RF_links_tmp <- llply(1:length(thr_enh), function(i) RF_links[[i]][rownames(RF_links[[i]]) %in% thr_enh[[i]],,drop=FALSE])
    names(RF_links_tmp) <- names(RF_links)
    Cor_links <- Cor_links[names(Cor_links_tmp)]
    Cor_links_tmp <- llply(1:length(thr_enh), function(i) Cor_links[[i]][rownames(Cor_links[[i]]) %in% thr_enh[[i]],,drop=FALSE])
    names(Cor_links_tmp) <- names(Cor_links)
  }
  
  prunedLinks <- list()
  if (!is.null(RF_links_tmp)){
    prunedLinks[['RF_links']] <- RF_links_tmp
  }
  if (!is.null(Cor_links_tmp)){
    prunedLinks[['Cor_links']] <- Cor_links_tmp
  }
  return(prunedLinks)
}


#' exportBB
#'
#' Export links to bigInteract format for visualization in UCSC.
#' @param RF_links Data frames list containing RF scores for each region in each gene as returned by enhancerToGene().
#' @param Cor_links Data frames list containing correlation scores for each region in each gene as returned by enhancerToGene(). If both RF and correlation links are provided, the height of the links will represent the RF importance and the color whether the correlation is positive or negative. If only RF is provided, links will be colored black; and if only correlation links are provided the height of the lnk will indicate the absolute correlation value and the color whether it is positive or negative.
#' @param annot Annotation data frame, as required by cicero (Pliner et al., 2019)
#' @param txdb Txdb object matching with the genome assembly used for the analysis
#' @param org.db Org.Db objet for the corresponding species
#' @param standardized Whether link scores (RF or correlation based) should be standardized per gene.
#' This is recommended for visualization in UCSC, as scores must be between 0-1000. Default=TRUE.
#' @param save_path Path to save bb file.
#' 
#' @return A dataframe containing links in bigInteract format. For more information, please visit
#' https://genome.ucsc.edu/goldenPath/help/interact.html. Scores values will
#' depend on whether RF links and or correlation links are provided, if both are provided
#' RF importances will be used for the score and correlations will be used for the color.
#'
#' @examples
#' exportBB(RF_links, Cor_links)
#'
#' @import plyr
#' @import AnnotationDbi
#' @import data.table
#'
#' @export


exportBB <- function(RF_links=NULL,
                      Cor_links=NULL,
                      txdb,
                      org.db,
                      standardized=TRUE,
                      save_path=NULL){
  
  # Check up
  if(!is(txdb,'TxDb')){
    stop('txdb has to be an object of class TxDb')
  }
  if(!is(org.db,'OrgDb')){
    stop('org.db has to be an object of class OrgDb')
  }
  genes <- names(RF_links)
  # Get search space around TSS
  # Genes to ensemble dict
  if (taxonomyId(org.db) == 9606){
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = genes, columns="ENTREZID", keytype="SYMBOL")
  } else {
    ENS2SYMBOL <- AnnotationDbi::select(org.db, keys = genes, columns="ENSEMBL", keytype="SYMBOL")
  }
  if (sum(is.na(ENS2SYMBOL[,2])) > 0){ENS2SYMBOL <- ENS2SYMBOL[-which(is.na(ENS2SYMBOL[,2])),]}
  # Select genes in the list
  filter_list <- list()
  filter_list[['GENEID']] <- ENS2SYMBOL[,2]
  # Take promoter coordinates for the specific genes
  TSS <- promoters(txdb, upstream = 0, downstream= 0, filter=filter_list, columns=c("GENEID"))
  ENS2SYMBOL_VECTOR <- as.vector(ENS2SYMBOL[,1])
  names(ENS2SYMBOL_VECTOR) <- ENS2SYMBOL[,2]
  elementMetadata(TSS)$GENEID <- ENS2SYMBOL_VECTOR[unlist(as.vector(elementMetadata(TSS)$GENEID))]
  TSS <- as.data.frame(TSS)
  TSS <- split(TSS, TSS$GENEID)
  TSS <- llply(TSS, function (x) x[1,])
  TSS <- llply(TSS, function (x) paste0(x[,1], ':', x[,2], '-', x[,3]))
  
  if (!is.null(RF_links) & !is.null(Cor_links)){
    RF_links <- RF_links[names(RF_links) %in% names(Cor_links)]
    if (standardized == TRUE){
      RF_links <- lapply(RF_links, function(x) {x[,1] <- (as.numeric(x[,1])-min(as.numeric(x[,1])))/(max(as.numeric(x[,1]))-min(as.numeric(x[,1])));x})
      RF_links <- lapply(RF_links, function(x) {x[x[,1] == 'NaN',1] <- 1;x})
    }
    Cor_links <- Cor_links[names(RF_links)]
    TSS <- TSS[names(RF_links)]
    BB <- as.data.frame(data.table::rbindlist(llply(1:length(RF_links), function(i) cbind(rownames(RF_links[[i]]), RF_links[[i]], Cor_links[[i]], rep(names(RF_links)[i], nrow(RF_links[[i]])), rep(TSS[[i]], nrow(RF_links[[i]]))))))
    BB[,3] <- as.numeric(BB[,3])
    color <- BB[!is.na(BB[,3]),3]
    mypal <- colorRampPalette(c('#FF0000', '#228B22'))(10)
    BB[!is.na(BB[,3]),3] <- .map2color(color,mypal)
    if(sum(is.na(BB[,3])) > 0){
      BB[which(is.na(BB[,3]))] <- 'black'
    } 
  }
  else if (!is.null(RF_link) & is.null(Cor_links)){
    if (standardized == TRUE){
      RF_links <- lapply(RF_links, function(x) {x[,1] <- (as.numeric(x[,1])-min(as.numeric(x[,1])))/(max(as.numeric(x[,1]))-min(as.numeric(x[,1])));x})
      RF_links <- lapply(RF_links, function(x) {x[x[,1] == 'NaN',1] <- 1;x})
    }
    TSS <- TSS[names(RF_links)]
    BB <- as.data.frame(data.table::rbindlist(llply(1:length(RF_links), function(i) cbind(rownames(RF_links[[i]]), RF_links[[i]], rep('black', nrow(RF_links[[i]])), rep(names(RF_links)[i], nrow(RF_links[[i]])), rep(TSS[[i]], nrow(RF_links[[i]]))))))
  }
  else if (is.null(RF_link) & !is.null(Cor_links)){
    if (standardized == TRUE){
      Cor_links <- lapply(Cor_links, function(x) {x[,1] <- (abs(as.numeric(x[,1]))-abs(min(as.numeric(x[,1]))))/(max(abs(as.numeric(x[,1])))-min(abs(as.numeric(x[,1]))));x})
      Cor_links <- lapply(Cor_links, function(x) {x[x[,1] == 'NaN',1] <- 1;x})
    }
    TSS <- TSS[names(Cor_links)]
    BB <- as.data.frame(data.table::rbindlist(llply(1:length(Cor_links), function(i) as.data.frame(cbind(rownames(Cor_links[[i]]), Cor_links[[i]], Cor_links[[i]], rep(names(Cor_links)[i], nrow(Cor_links[[i]])), rep(TSS[[i]], nrow(Cor_links[[i]])))))))
    BB[,3] <- as.numeric(BB[,3])
    color <- BB[!is.na(BB[,3]),3]
    mypal <- colorRampPalette(c('#FF0000', '#228B22'))(10)
    BB[!is.na(BB[,3]),3] <- .map2color(color,mypal)
    if(sum(is.na(BB[,3])) > 0){
      BB[which(is.na(BB[,3]))] <- 'black'
    } 
  } else {
    stop('Please, provide at least either correlation or RF links.')
  }
    BB[,2] <- as.numeric(BB[,2])
    colnames(BB) <- c('Enhancer', 'Score', 'Color', 'Gene', 'TSS')
    BB <- as.matrix(BB)
    Enhancer_seqnames <- sapply(strsplit(BB[,1], split = ":"), "[", 1)
    Enhancer_coord <- sapply(strsplit(BB[,1], split = ":"), "[", 2)
    Enhancer_start <- round((as.numeric(sapply(strsplit(Enhancer_coord, split = "-"), "[", 1))+as.numeric(sapply(strsplit(Enhancer_coord, split = "-"), "[", 2)))/2)
    Enhancer_end <- Enhancer_start+1
    
    TSS_seqnames <- sapply(strsplit(BB[,5], split = ":"), "[", 1)
    TSS_coord <- sapply(strsplit(BB[,5], split = ":"), "[", 2)
    TSS_start <- round((as.numeric(sapply(strsplit(TSS_coord, split = "-"), "[", 1))+as.numeric(sapply(strsplit(TSS_coord, split = "-"), "[", 2)))/2)
    TSS_end <- TSS_start+1
    
    BB <- cbind(Enhancer_seqnames, Enhancer_start, TSS_end, BB[,'Gene'], round(as.numeric(BB[,'Score'])*1000),
                round(as.numeric(BB[,'Score'])*10), BB[,'Gene'], BB[, 'Color'], Enhancer_seqnames, Enhancer_start,
                Enhancer_end, BB[,'Enhancer'], rep('.', nrow(BB)), TSS_seqnames, TSS_start, 
                TSS_end, BB[, 'TSS'], rep('.', nrow(BB)))
    
    BB[which(as.numeric(BB[,3]) < as.numeric(BB[,2])), c(2,3)] <- BB[which(as.numeric(BB[,3]) < as.numeric(BB[,2])), c(3,2)]
    BB <- as.data.frame(BB)
    colnames(BB) <- c('chrom', 'chromStart', 'chromEnd', 'name', 'score', 'value', 'exp', 'color',
                      'sourceChrom', 'sourceStart', 'sourceEnd', 'sourceName', 'sourceStrand',
                      'targetChrom', 'targetStart', 'targetEnd', 'targetName', 'targetStrand')
    if (!is.null(file)){
      write.table(BB, file=save_path, row.names=FALSE, col.names = FALSE, quote=FALSE,  sep = "\t", eol = "\n")
    }
    return(BB)
}

# Helper fuction
.map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}
