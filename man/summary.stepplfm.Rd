\name{summary.stepplfm}
\alias{summary.stepplfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for \code{\link{stepplfm}} objects}
\description{
  The function \code{summary.stepplfm} summarizes the fit measures of a series of \code{\link{plfm}} objects. 
}
\usage{\method{summary}{stepplfm}(object, \dots)}
\arguments{
  \item{object}{Probabilistic latent feature analysis object returned by \code{\link{stepplfm}}}
  \item{\dots}{Further arguments are ignored}
}
\details{
  The summary of the \code{\link{stepplfm}} generates a table of fit measures (information criteria, Chi-square fit on the \emph{J X K} table, and descriptive fit measures) 
  for models with different numbers of features and/or mapping rules that are involved in the analysis.
}
\author{Michel Meulders}
\seealso{
\code{\link{stepplfm}}, \code{\link{print.stepplfm}}, \code{\link{print.summary.stepplfm}}
}


