\name{meta.abf}
\alias{meta.abf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Meta-Analysis with an Approximate Bayes Factor}

\description{This function performs a meta-analysis on summary statistics from multiple genome-wide association studies. It takes as its inputs the effect size estimates (\code{betas}) and associated standard errors (\code{ses}), and for some prior on true effect size and on the relationship between the effect size across studies, returns a Bayes factor that shows the evidence in favour of the alternative model against the null model of no true effect.}
\usage{
meta.abf(betas, ses, prior.sigma, prior.cor = "indep",
prior.rho = NA, cryptic.cor = NA, log = FALSE,
log10 = FALSE, tolerance = 1e-1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{This is a vector, matrix, or data frame of observed effect sizes of a single SNP in a set of studies.  If \code{betas} is a vector, then this assumes a single SNP and each element of the vector corresponds to a study. If \code{betas} is a matrix, then the rows correspond to SNPs and the columns to studies.}
  \item{ses}{This is a vector, matrix, or data frame of standard errors corresponding to those in \code{betas}. Its class and dimensions should be identical to those of \code{betas}.}
  \item{prior.sigma}{This is the prior on true effect sizes for each SNP in each study.  It can be a single value that applies to all studies, set individually for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis) or set for each study and SNP (i.e. a matrix of same dimension as \code{betas} and \code{ses}).}
  \item{prior.cor}{This is the prior on the relationship between true effect sizes in each study.  It is a square matrix whose row and column numbers are the same as the number of studies.  Its elements are the pairwise correlation coefficients between true effect sizes of the studies.

  It can also take values \dQuote{indep} (independent effects, where the values are 0 uniformly), \dQuote{fixed} (fixed effects, where the values are 1 uniformly), \dQuote{correlated} (uniform correlated effects, which requires the prior.rho parameter to be set).  You may also simply input the desired matrix directly.

  If \code{betas} and \code{ses} are matrices, the same prior.cor will be applied to every row (representing every SNP).}
  \item{prior.rho}{If \code{prior.cor} is set to \dQuote{correlated}, then this is the value of the off-diagonal elements when the pair-wise correlation coefficients are all assumed to be the same value (for instance, all 0.5 or all 1).  If \code{prior.cor} takes any other value, \code{prior.rho} is ignored.}
  \item{cryptic.cor}{This is a square matrix whose row and coumn numbers are the same as the number of studies.  The off-diagonal elements are the pairwise correlation coefficients between the studies.  If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.}
  \item{log}{When \code{TRUE}, the calculation will return the log of the Bayes factors. Cannot be \code{TRUE} if \code{log10 = TRUE}.}
  \item{log10}{When \code{TRUE}, the calculation will return the log10 of the Bayes factors. Cannot be \code{TRUE} if \code{log = TRUE}.}
  \item{tolerance}{This is the relative tolerance to detect zero singular values.  Used in an internal function during the calculation of the approximate Bayes factor.  This should probably never be altered, but is there as a last resort if the calculation is not returning the values that are expected.}
}

\value{Returns the approximate Bayes factor for the data under the given priors against the null. For a single SNP, this is a single value. For a matrix with rows of SNPs, this will return a numeric vector.}

\note{Missing data should be marked \code{NA}. Where these are included, the calculation will be performed with the subset of the data for which there is information.}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{Trochet, Holly and Pirinen, Matti}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{exh.abf}}, \code{\link{wakefield.abf}}}
\examples{
#Obtain the matrices of betas and standard errors
betas <- cotsapas[,grep("\\\\.beta", names(cotsapas))]
ses <- cotsapas[,grep("\\\\.SE", names(cotsapas))]

betas[1:10,]
ses[1:10,]

#Perform an independent effects meta-analysis
meta.abf(betas, ses, prior.sigma = 0.2, prior.cor = "indep", log10 = TRUE)

#Single SNP analysis
i <- 1
meta.abf(betas[i,], ses[i,], prior.sigma = 0.2, prior.cor = "correlated", prior.rho = 0.5, log10 = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ htest }
