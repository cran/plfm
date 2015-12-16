\name{print.LCplfm}
\alias{print.LCplfm}
\title{Printing LCplfm objects}
\description{Printing method for latent class probabilistic feature analysis objects.}
\usage{
\method{print}{LCplfm}(x,\dots)
}
\arguments{
\item{x}{Latent class probabilistic  feature analysis object returned by \code{\link{LCplfm}}.}
\item{\dots}{Further arguments are ignored.}
}
\details{The printing method for latent class probabilistic  feature analysis objects displays 
	(1) the parameters used to call the \code{\link{LCplfm}} function, 
        (2) the estimated object-, attribute- and class size parameters,
        (3) the estimated standard errors of object-, attribute- and class size parameters,
        (4) fit measures  
} 
%%\seealso{\code{\link{print.summary.LCplfm}}}
\examples{
\dontrun{
# example 
print(LCplfm(data=anger$data,F=2, T=2, M=1))
}
}

