\name{plot.LCplfm}
\alias{plot.LCplfm}
\title{plot parameters in \code{\link{LCplfm}} object}
\description{Plot method to visualize the parameters  of latent class probabilistic feature models with different numbers of features/classes.}
\usage{
\method{plot}{LCplfm}(x, feature=1, class=0, element="object", cexsymb=1, cexlabel=1,
                      minpositionlabel = -1, positionlabel = -0.8, xlegend = "topright", 
                      ylegend=NULL, x.intersplegend=1, y.intersplegend=1, \dots)
}
\arguments{
\item{x}{Latent class probabilistic feature model object returned by \code{\link{LCplfm}}.}
\item{feature}{Latent feature for which parameters are visualized.}
\item{class}{Latent class for which parameters are visualized. When the model contains class-specific object- or attribute parameters, \code{class=0} means that parameters of 
             all classes are included in the plot.}
\item{element}{Object parameters are plotted if \code{element}="object" and attribute parameters are plotted if \code{element}="attribute".}
\item{cexsymb}{Size of symbol used for plotting points.}
\item{cexlabel}{Size of object- or attribute labels in plot.}
\item{minpositionlabel}{Value smaller than 0 that defines space for plotting object- or attribute labels.}
\item{positionlabel}{Value between \code{minpositionlabel} and 0 to align object- or attribute labels.}
\item{xlegend, ylegend}{The x and y co-ordinates to be used to position the legend. They can be specified by keyword or in 
      any way which is accepted by xy.coords: See "Details" of legend.}
\item{x.intersplegend}{Character interspacing factor for horizontal (x) spacing in legend.}
\item{y.intersplegend}{Character interspacing factor for vertical (y) line distances in legend.}
\item{\dots}{Further arguments are ignored.}
}
\examples{

\dontrun{
# example 1: analysis on determinants of anger-related behavior

# load anger data
data(anger)

# compute 5 runs of disjunctive latent class probabilistic feature model
# with 4 features and 2 latent classes
# assume constant situation classification per person 
# and class-specific situation parameters (i.e. model=1) 

anger.m1<-LCplfm(data=anger$data,F=4,T=2,maprule="disj",
                  M=5,emcrit1=1e-3,emcrit2=1e-8,model=1)

# visualize object and attribute parameters 
# of both classes per feature in one figure

par(mfrow=c(2,2),pty="s")
plot(anger.m1,element="attribute",feature=1, main="Feature 1", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m1,element="attribute",feature=2, main="Feature 2", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m1,element="attribute",feature=3, main="Feature 3", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m1,element="attribute",feature=4, main="Feature 4", 
     minpositionlabel=-2, positionlabel=-1)

par(mfrow=c(2,2),pty="s")
plot(anger.m1,element="object",feature=1,main="Feature 1", 
     minpositionlabel=-1.5, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m1,element="object",feature=2,main="Feature 2", 
     minpositionlabel=-1.5, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m1,element="object",feature=3,main="Feature 3", 
     minpositionlabel=-1.5, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m1,element="object",feature=4,main="Feature 4", 
     minpositionlabel=-1.5, positionlabel=-1, y.intersplegend=0.7)

# compute 5 runs of disjunctive latent class probabilistic feature model
# with 2 features and 2 latent classes
# assume constant situation classification per person 
# and class-specific situation and behavior parameters (i.e. model=3) 

anger.m3<-LCplfm(data=anger$data,F=2,T=2,maprule="disj",
                  M=5,emcrit1=1e-3,emcrit2=1e-8,model=3)

# visualize object and attribute parameters of feature 1,2 
# for class 1
par(mfrow=c(2,2)) 
plot(anger.m3,element="attribute",feature=1, class=1,main="Feature 1, class 1", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m3,element="attribute",feature=2, class=1,main="Feature 2, class 1", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m3,element="object",feature=1, class=1,main="Feature 1, class 1", 
     minpositionlabel=-2, positionlabel=-1)
plot(anger.m3,element="object",feature=2, class=1,main="Feature 2, class 1", 
     minpositionlabel=-2, positionlabel=-1)


# visualize object and attribute parameters of feature 1,2 
# for class 2
par(mfrow=c(2,2))
plot(anger.m3,element="attribute",feature=1, class=2,main="Feature 1, class 2", 
     minpositionlabel=-1.7, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m3,element="attribute",feature=2, class=2,main="Feature 2, class 2", 
     minpositionlabel=-1.7, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m3,element="object",feature=1, class=2,main="Feature 1, class 2", 
     minpositionlabel=-1.7, positionlabel=-1, y.intersplegend=0.7)
plot(anger.m3,element="object",feature=2, class=2,main="Feature 2, class 2", 
     minpositionlabel=-1.7, positionlabel=-1, y.intersplegend=0.7)

}
}