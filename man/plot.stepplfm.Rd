\name{plot.stepplfm}
\alias{plot.stepplfm}
\title{Plot fit of \code{\link{stepplfm}} objects}
\description{Plot method to visualize the fit of probabilistic latent feature analysis objects with different numbers of features.}
\usage{
\method{plot}{stepplfm}(x,which="BIC",\dots)
}
\arguments{
\item{x}{List of probabilistic latent feature analysis objects returned by \code{\link{stepplfm}}.}
\item{which}{Fit criterion for which models with different numbers of features are compared. 
             The argument \code{which} can take the following values: 
             "AIC", "BIC", "Deviance", "Chisquare", "Correlation", "VAF"}
\item{\dots}{Further arguments are ignored.}
}
\examples{

\dontrun{
## example 1:Perceptual analysis of associations between car models and car attributes

##load car data
data(car)


## compute 5 runs of disjunctive and conjunctive models with 1 up to 4 features
car.lst<-stepplfm(minF=1,maxF=4,maprule="disj/conj",freq1=car$freq1,
                  freqtot=car$freqtot,M=5)

## visualize the fit of models with different mapping rules 
## and a different number of features

par(pty="s")
par(mfrow=c(2,2))
plot(car.lst,which="BIC")
plot(car.lst,which="AIC")
plot(car.lst,which="VAF")

}

## example 2: analysis on determinants of anger-related behavior

## load anger data
data(anger)

## compute 1 run of disjunctive and conjunctive models with 1 up to 3 features
anger.lst<-stepplfm(minF=1,maxF=3,maprule="disj/conj",freq1=anger$freq1,
                    freqtot=anger$freqtot,M=1)

## visualize the fit of models with different mapping rules 
## and a different number of features

par(pty="s")
par(mfrow=c(2,2))
plot(anger.lst,which="BIC")
plot(anger.lst,which="AIC")
plot(anger.lst,which="VAF")


}

