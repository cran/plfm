\name{print.plfm}
\alias{print.plfm}
\title{Printing plfm objects}
\description{Printing method for probabilistic latent feature analysis objects.}
\usage{
\method{print}{plfm}(x,\dots)
}
\arguments{
\item{x}{Probabilistic latent feature analysis object returned by \code{\link{plfm}}.}
\item{\dots}{Further arguments are ignored.}
}
\details{The printing method for probabilistic latent feature analysis objects displays 
	(1) the parameters used to call the \code{\link{plfm}} function, 
	(2) information on the descriptive fit of the model (i.e. correlation between observed and expected frequencies, 
	and proportion of the variance in the observed frequencies accounted for by the model), 
	and (3) the estimated object- and attribute parameters.
} 
%%\seealso{\code{\link{print.summary.plfm}}}
\examples{
## example print.plfm(plfm(...))
}

