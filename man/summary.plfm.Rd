\name{summary.plfm}
\alias{summary.plfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing probabilistic latent feature analysis}
\description{
  The function \code{summary.plfm} summarizes the main output of \code{\link{plfm}} including estimates and standard errors for 
  object- and attribute parameters, model selection criteria, and goodness-of-fit measures.
}
\usage{\method{summary}{plfm}(object, \dots)}
\arguments{
  \item{object}{Probabilistic latent feature analysis object returned by \code{\link{plfm}}}
  \item{\dots}{Further arguments are ignored}
}
\details{
The summary of probabilistic latent feature analysis objects displays: 
\enumerate{
\item The parameters used to call the \code{\link{plfm}} function.
\item The value of the loglikelihood, the deviance, the logarithm of the posterior density, the information criteria AIC and BIC.
\item The result of a Pearson chi-square goodness-of-fit test on the \emph{J X K} table.
\item Information on the descriptive fit of the model (i.e. correlation between observed and expected frequencies.
      and proportion of the variance in the observed frequencies accounted for by the model).
\item The estimated object- and attribute parameters.
\item Asymptotic standard errors of the object- and attribute parameters.
}
}
\value{
  \item{call}{Parameters used to call the function.}
  \item{informationcriteria}{List of information criteria that can be used for model selection.}
  \item{chisquaretest}{Pearson Chi-square test to evaluate the statistical goodness-of-fit of the model on the \emph{J X K} object by attribute table of association frequencies.}
  \item{descriptivefit}{A list of measures to evaluate the descriptive goodness-of-fit of the model on the \emph{J X K} object by attribute table of association frequencies.}
  \item{objpar}{A \emph{J X F} matrix of estimated object parameters.}
  \item{SE.objpar}{A \emph{J X F} matrix of estimated standard errors of object parameters.}
  \item{attpar}{A \emph{K X F} matrix of estimated attribute parameters.}
  \item{SE.attpar}{A \emph{K X F} matrix of estimated standard errors of attribute parameters}
}
\author{Michel Meulders}
\seealso{
\code{\link{plfm}}, \code{\link{print.plfm}}, \code{\link{print.summary.plfm}}
}
\examples{
## Perceptual analysis of associations between car models and car attributes

##load car data
data(car)

##compute the disjunctive model with 4 features
carf4<-plfm(maprule="disj",freq1=car$freq1,freqtot=car$freqtot,F=4,M=1)

## display a summary of the results
summary(carf4)
}
