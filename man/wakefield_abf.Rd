\name{wakefield.abf}
\alias{wakefield.abf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Wakefield's Approximate Bayes Factor}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
wakefield.abf(beta,se,prior.sigma,log=FALSE,log10=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{beta}{A single value or numeric vector. It should represent the observed effect size of a SNP on a trait from a genome-wide association study.}
  \item{se}{A single value or numeric vector that corresponds to the value(s) in \code{beta}. It should represent the standard error of the effect size of a SNP on a trait from a genome-wide association study.}
  \item{prior.sigma}{A single value or numeric vector.  This is the prior on true effect sizes.  If \code{beta} and \code{se} are vectors, this \code{prior.sigma} is permitted to be a single value, in which case, it will be applied over all values of \code{beta} and \code{se}.}
  \item{log}{When \code{TRUE}, the calculation will return the log of the Bayes factors. Cannot be \code{TRUE} if \code{log10 = TRUE}.}
  \item{log10}{When \code{TRUE}, the calculation will return the log10 of the Bayes factors. Cannot be \code{TRUE} if \code{log = TRUE}.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Wakefield, Jon. Bayes Factors for Genome-Wide Association Studies: Comparison with P-values. \emph{Genetic Epidemiology.} 33:79--76. 2009.}
\author{Trochet, Holly}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{meta.abf}},\code{\link{exh.abf}}}
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
\keyword{ univariate }
\keyword{ htest }
