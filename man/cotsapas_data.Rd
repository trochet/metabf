\name{cotsapas}
\alias{cotsapas}
\docType{data}
\title{Summary Statistics from \dQuote{Pervasive Sharing of Genetic Effects in Autoimmune Disease}}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{cotsapas}
\format{
  A data frame with the summary statistics for seven autoimmune diseases across 107 single nucleotide polymorphisms. The columns are described below.
  \describe{
    \item{\code{SNP}}{The rsid of the single nucleotide polymorphism (SNP).}
    \item{\code{CHR}}{The chromosome where the SNP is located.}
    \item{\code{POS}}{The basepair position of the SNP on the chromosome.}
    \item{\code{major_al}}{The major allele of the SNP.}
    \item{\code{minor_al}}{The minor allele of the SNP.}
    \item{\code{freq}}{The frequency of the minor allele.}
    \item{\code{RA.}}{Values with this prefix pertain to data from a genome-wide associiation study of rheumatoid arthritis.}
    \item{\code{PS.}}{Values with this prefix pertain to data from a genome-wide associiation study of psoriasis.}
    \item{\code{MS.}}{Values with this prefix pertain to data from a genome-wide associiation study of multiple sclerosis.}
    \item{\code{SLE.}}{Values with this prefix pertain to data from a genome-wide associiation study of systemic lupus erythematosus.}
    \item{\code{CD.}}{Values with this prefix pertain to data from a genome-wide associiation study of Crohn's disease.}
    \item{\code{CeD.}}{Values with this prefix pertain to data from a genome-wide associiation study of coeliac disease.}
    \item{\code{T1D.}}{Values with this prefix pertain to data from a genome-wide associiation study of type 1 diabetes.}
    \item{\code{.beta}}{Values with this suffix show the effect size estimate (or log of the odds ratio) of the SNP on the disease.}
    \item{\code{.SE}}{Values with this suffix show the standard error of the effect size estimate of the SNP on the disease.}
    \item{\code{.P}}{Values with this suffix show the \emph{p}-value of the effect of SNP on the disease.}
  }
}
\details{These data come from studies collated by Cotsapas \emph{et al.} (see References below). The original paper listed Z-scores and \emph{p}-values. We have used 1000 Genomes and the individual study sizes to estimate the standard errors and from that, the effect sizes. The listed SNP positions have also been update to match GRCh37.}

\source{Dataset S1 from \dQuote{Pervasive Sharing of Genetic Effects in Autoimmune Disease}. (See reference below)}

\references{Cotsapas C, Voight BF, Rossin E, Lage K, Neale BM, Wallace C, \emph{et al.} (2011) Pervasive Sharing of Genetic Effects in Autoimmune Disease. \emph{PLoS Genetics} 7(8): e1002254. doi:10.1371/journal.pgen.1002254

The 1000 Genomes Projet Consortium. (2015) A global reference for human genetic variation. \emph{Nature} \strong{526}, 68--74. doi:10.1038/nature15393}
\examples{
cotsapas[1:10,]
names(cotsapas)

#To obtain all the effect size estimates
betas <- cotsapas[,grep("\\\\.beta", names(cotsapas))]}
\keyword{datasets}
