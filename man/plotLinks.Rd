% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enhancerToGene.R
\name{plotLinks}
\alias{plotLinks}
\title{plotLinks}
\usage{
plotLinks(
  RF_links = NULL,
  Cor_links = NULL,
  annot,
  txdb,
  org.db,
  gene,
  chr,
  start,
  end,
  cutoff = 0
)
}
\arguments{
\item{RF_links}{Data frames list containing RF scores for each region in each gene as returned by enhancerToGene().}

\item{Cor_links}{Data frames list containing correlation scores for each region in each gene as returned by enhancerToGene(). If both RF and correlation links are provided, the height of the links will represent the RF importance and the color whether the correlation is positive or negative. If only RF is provided, links will be colored black; and if only correlation links are provided the height of the lnk will indicate the absolute correlation value and the color whether it is positive or negative.}

\item{annot}{Annotation data frame, as required by cicero (Pliner et al., 2019)}

\item{txdb}{Txdb object matching with the genome assembly used for the analysis}

\item{org.db}{Org.Db objet for the corresponding species}

\item{gene}{Gene for which enhancer-to-gene links wants to be plotted}

\item{chr}{Chromosome name of the genomic window that wants to be plotted}

\item{start}{Start position of the genomic window that wants to be plotted}

\item{end}{End position of the genomic window that wants to be plotted}

\item{cutoff}{Value below which links will not be shown. Default: -1.}
}
\description{
Plot enhancer-to-gene links (To be enhanced)
}
\examples{
plotLinks(RF_links, dm6_annot, TxDb.Dmelanogaster.UCSC.dm6.ensGene, org.Dm.eg.db, gene='dac', 'chr2L', 16470000, 16490000)

}
