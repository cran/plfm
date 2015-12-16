\name{plot.stepLCplfm}
\alias{plot.stepLCplfm}
\title{Plot fit of \code{\link{stepLCplfm}} objects}
\description{Plot method to visualize the fit of latent class probabilistic feature models with different numbers of features/classes.}
\usage{
\method{plot}{stepLCplfm}(x,which="BIC",\dots)
}
\arguments{
\item{x}{List of latent class probabilistic latent feature analysis objects returned by \code{\link{stepLCplfm}}.}
\item{which}{Fit criterion for which models with different numbers of features are compared. 
             The argument \code{which} can take the following values: 
             "AIC", "BIC", "Deviance", "Correlation", "VAF"}
\item{\dots}{Further arguments are ignored.}
}
\examples{

\dontrun{
# example 1: analysis on determinants of anger-related behavior

# load anger data
data(anger)

# compute 5 runs of disjunctive latent class probabilistic feature models
# with 1 up to 3 features and with 1 up to 2 latent classes
# assume constant situation classification per person 
# and class-specific situation parameters (i.e. model=1) 

anger.lst<-stepLCplfm(minF=1,maxF=3,minT=1,maxT=2,data=anger$data,
                      maprule="disj",M=5,emcrit1=1e-3,emcrit2=1e-8,model=1)


# visualize BIC of fitted models 

par(pty="s")
plot(anger.lst)

# print overview fit measures for all estimated models

anger.lst

# print model with 3 features and 1 latent class

anger.lst[[3,1]]
}

\dontrun{
# example 2:Perceptual analysis of associations between car models and car attributes

# load car data
data(car)


# compute 5 runs of disjunctive models with 4 features and 1 up to 3 latent classes
# assume constant attribute classification per respondent 
# and class-specific car parameters (i.e. model 4)

car.lst<-stepLCplfm(minF=4,maxF=4,minT=1,maxT=3,data=car$data3w,
                      maprule="disj",M=5,emcrit1=1e-3,emcrit2=1e-8,model=4,printrun=TRUE)


# visualize BIC of fitted models
plot(car.lst)

# print overview of fitmeasures for all fitted models
car.lst
}
}