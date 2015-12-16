\name{print.stepLCplfm}
\alias{print.stepLCplfm}
\title{Printing \code{\link{stepLCplfm}} objects}
\description{Printing method for a series of latent class probabilistic latent feature analysis objects.}
\usage{
\method{print}{stepLCplfm}(x,\dots)
}
\arguments{
\item{x}{Latent class probabilistic latent feature analysis object returned by \code{\link{stepLCplfm}}.}
\item{\dots}{Further arguments are ignored.}
}
\details{The printing method for \code{\link{stepLCplfm}} displays summary tables about the fit of models with different numbers of features and different numbers of latent classes. 
         Two tables are printed which summarize the fit of models with different numbers of features/classes in terms of (1) information criteria (AIC, BIC,...), 
          and (2) the descriptive fit of the model to the \emph{JXK} 
         frequency table (correlation between observed and expected frequencies, and variance accounted for by the model).
	
} 
\examples{
# example print.stepLCplfm(stepLCplfm(...))
}

