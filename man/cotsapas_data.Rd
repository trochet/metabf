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
    \item{\code{RA.}}{values with this prefix pertain to data from a genome-wide associiation study of rheumatoid arthritis.}
    \item{\code{PS.}}{values with this prefix pertain to data from a genome-wide associiation study of psoriasis.}
    \item{\code{MS.}}{values with this prefix pertain to data from a genome-wide associiation study of multiple sclerosis.}
    \item{\code{SLE.}}{values with this prefix pertain to data from a genome-wide associiation study of systemic lupus erythematosus.}
    \item{\code{CD.}}{values with this prefix pertain to data from a genome-wide associiation study of Crohn's disease.}
    \item{\code{CeD.}}{values with this prefix pertain to data from a genome-wide associiation study of coeliac disease.}
    \item{\code{T1D.}}{values with this prefix pertain to data from a genome-wide associiation study of type 1 diabetes.}
    \item{\code{.beta}}{values with this suffix show the effect size estimate (or log of the odds ratio) of the SNP on the disease.}
    \item{\code{.SE}}{values with this suffix show the standard error of the effect size estimate of the SNP on the disease.}
    \item{\code{.P}}{values with this suffix show the \emph{p}-value of the effect of SNP on the disease.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{\dQuote{Pervasive Sharing of Genetic Effects in Autoimmune Disease}. Cotsapas \emph{et. al.} }
\examples{
cotsapas[1:10,]}
\keyword{datasets}
