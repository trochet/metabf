\name{exh.abf}
\alias{exh.abf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Subset-Exhaustive Approximate Bayes Factors}
\description{This function takes the summary data at a single SNP across a set of genome-wide association studies and user-set priors to calculate approximate Bayes factors for every possible subset of the data.}
\usage{
exh.abf(betas, ses, prior.sigma, prior.cor = "indep",
prior.rho = NA, cryptic.cor = NA, log = FALSE,
log10 = FALSE, study.names = NULL, na.rm = FALSE,
tolerance=1e-1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{betas}{This is a numeric vector of observed effect sizes of a single SNP in a set of studies.  Each element of the vector is assumed to correspond to a study.}
  \item{ses}{This is a numeric vector of standard errors corresponding to those in betas.  It should have the same length as betas.}
  \item{prior.sigma}{This is the prior on the true effect sizes for the SNP in each study.  It can be a flat value (i.e. the same for all studies) or set for each study (i.e. a vector whose length is equal to the number of studies in the meta-analysis).}
  \item{prior.cor}{This is the prior on the relationship between true effect sizes in each study.  It is a square matrix whose row and column numbers are the same as the number of studies.  Its elements are the pairwise correlation coefficients between true effect sizes of the studies.

  It can take values \dQuote{indep} (independent effects, where the values are 0 uniformly), \dQuote{fixed} (fixed effects, where the values are 1 uniformly), \dQuote{correlated} (correlated effects, which requires the prior.rho parameter to be set).  You may also simply input the desired matrix directly.}
  \item{prior.rho}{If \code{prior.cor} is set to \dQuote{correlated}, then this is the value of the off-diagonal elements.  This can be either a single value for all pair-wise correlation coefficients or it can be the upper triangle of the desired matrix.  If \code{prior.cor} takes any other value, \code{prior.rho} is ignored.}
  \item{cryptic.cor}{This is a square matrix whose row and coumn numbers are the same as the number of studies.  The off-diagonal elements are the pairwise correlation coefficients between the studies.  If the studies in the meta-analysis are not independent of each other, it may be necessary to set this parameter so that the covariance in null effects is accounted for.}
  \item{log}{When \code{TRUE}, the calculation will return the log of the Bayes factors. Cannot be \code{TRUE} if \code{log10 = TRUE}.}
  \item{log10}{When \code{TRUE}, the calculation will return the log10 of the Bayes factors. Cannot be \code{TRUE} if \code{log = TRUE}.}
  \item{study.names}{Optional parameter that, if set, labels the output the columns with the study names.}
  \item{na.rm}{By default, the output will include a column of NAs for the study with the missing data. Changing this parameter to TRUE removes this column.}
  \item{tolerance}{This is the relative tolerance to detect zero singular values.  Used in an internal function during the calculation of the approximate Bayes factor.  This should probably never be altered, but is there as a last resort if the calculation is not returning the values that are expected.}
}
\details{To work properly, the values in \code{betas} must reflect the effect of the same allele in the same SNP across all studies.}
\value{Returns a data frame where the first \emph{n} columns correspond to the studies in \code{betas} and \code{ses}. If \code{study.names = TRUE}, then these will be labelled with the user-provided study names.

Each line is a binary vector where \code{1} corresponds to the studies with true effects in the model and \code{0} corresponds to the studies without. For example, the string that is uniformly \code{0} is the null model, and the string that is uniformly \code{1} is the model where every study shows a true effect.

The final column, labelled \emph{ABF} indicates the approximate Bayes factor calculated for each of these models against the null.}
\references{
%% ~put references to the literature/web site here ~
}
\author{Trochet, Holly and Pirinen, Matti}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{meta.abf}},\code{\link{wakefield.abf}}}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{ htest }
