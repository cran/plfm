\name{stepLCplfm}
\alias{stepLCplfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Latent class probabilistic latent feature analysis of three-way three-mode binary data}
\description{The function \code{\link{stepLCplfm}} subsequently applies the \code{\link{LCplfm}} function to fit disjunctive and/or conjunctive models with \emph{minF} up to \emph{maxF} 
latent features and \emph{minT} to \emph{maxT} latent classes. The results of the estimated models are stored in a list with \emph{F X T} components.}  
\usage{
	stepLCplfm(minF=1,maxF=3,minT=1,maxT=3,
                   data,maprule="disj",M=5,emcrit1=1e-3,emcrit2=1e-8,
                   model=1,delta=0.0001,printrun=FALSE)
}
\arguments{
  \item{minF}{Minimum number of latent features included in the model.}
  \item{maxF}{Maximum number of latent features included in the model.}
  \item{minT}{Minimum number of latent classes included in the model.}
  \item{maxT}{Maximum number of latent classes included in the model.}
  \item{data}{A \emph{I X J X K} data array of binary observations. Observation \emph{(i,j,k)} \emph{(i=1,..,I; j=1,..,J; k=1,..,K)}  indicates 
               whether object \emph{j} is associated to attribute \emph{k} according to rater \emph{i}.} 
  \item{maprule}{Fit disjunctive models (\code{maprule}="disj") or conjunctive models (\code{maprule}="conj").}
  \item{M}{The number of exploratory runs of the EM algorithm using random starting points for each model.}
  \item{emcrit1}{Convergence criterion to be used for the estimation of candidate models in the exploration step.}
  \item{emcrit2}{Convergence criterion to be used for the estimation of the best model in the final analysis.}
  \item{model}{The type of dependency and heterogeneity assumption included in the model. \code{model}=1, \code{model}=2, \code{model}=3 represent models with a constant 
               object-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters. \code{model}=4, \code{model}=5, \code{model}=6 represent models with a constant 
               attribute-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters.}
  \item{delta}{The precision used to compute standard errors of the model parameters with the method of finite differences.}
  \item{printrun}{\code{printrun}=TRUE prints the analysis type (disjunctive or conjunctive), the number of features (\emph{F}), the number of latent classes (\emph{T}) 
                  and the number of the run to the output screen, whereas \code{printrun}=FALSE suppresses the printing.}
}

\details{The results of subsequent \code{\link{LCplfm}} analyses are stored in a matrix of lists with \code{(maxF-minF+1,maxT-minT+1)} components.} 
         

\author{
Michel Meulders
}
%%\note{
%%  ~~further notes~~
%%}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
