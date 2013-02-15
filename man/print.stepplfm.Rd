\name{print.stepplfm}
\alias{print.stepplfm}
\title{Printing \code{\link{stepplfm}} objects}
\description{Printing method for a series of probabilistic latent feature analysis objects.}
\usage{
\method{print}{stepplfm}(x,\dots)
}
\arguments{
\item{x}{probabilistic latent feature analysis object returned by \code{\link{stepplfm}}.}
\item{\dots}{Further arguments are ignored.}
}
\details{The printing method for \code{\link{stepplfm}} displays summary tables about the fit of models with different numbers of features. 
         For each mapping rule, three tables are printed which summarize the fit of models with different numbers of features in terms of (1) information criteria (AIC, BIC,...), 
         (2) the statistical fit of the model to the \emph{JXK} frequency table (Chi-square value, df and corresponding p-value), and (3) the descriptive fit of the model to the \emph{JXK} 
         frequency table (correlation between observed and expected frequencies, and variance accounted for by the model).
	
} 
\examples{
## example print.stepplfm(stepplfm(...))
}

