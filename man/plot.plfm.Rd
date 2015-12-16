\name{plot.plfm}
\alias{plot.plfm}
\title{plot parameters in \code{\link{plfm}} object}
\description{Plot method to visualize the parameters  of probabilistic feature models.}
\usage{
\method{plot}{plfm}(x,feature=1,element="object",cexsymb=1,cexlabel=1,\dots)
}
\arguments{
\item{x}{Probabilistic feature model object returned by \code{\link{plfm}}.}
\item{feature}{Latent feature for which parameters are visualized.}
\item{element}{Object parameters are plotted if \code{element}="object" and attribute parameters are plotted if \code{element}="attribute".}
\item{cexsymb}{Size of symbol used for plotting points.}
\item{cexlabel}{Size of object- or attribute labels in plot.}
\item{\dots}{Further arguments are ignored.}
}
\examples{

# examples


\dontrun{
# example 1:Perceptual analysis of associations between car models and car attributes

# load car data
data(car)

# compute 1 run of a disjunctive model with 4 features
# use components of a data frame as input
cardisj4<-plfm(datatype="dataframe",data=car$datalongformat,object=objectlabel,
               attribute=attributelabel,rating=rating,maprule="disj",F=4,M=1)

# plot car and attribute parameters per feature
par(mfrow=c(1,2))
plot(cardisj4,feature=1,element="object",main="Car parameters Feature 1")
plot(cardisj4,feature=1,element="attribute",main="Attribute parameters Feature 1")

par(mfrow=c(1,2))
plot(cardisj4,feature=2,element="object",main="Car parameters Feature 2")
plot(cardisj4,feature=2,element="attribute",main="Attribute parameters Feature 2")

par(mfrow=c(1,2))
plot(cardisj4,feature=3,element="object",main="Car parameters Feature 3")
plot(cardisj4,feature=3,element="attribute",main="Attribute parameters Feature 3")

par(mfrow=c(1,2))
plot(cardisj4,feature=4,element="object",main="Car parameters Feature 4")
plot(cardisj4,feature=4,element="attribute",main="Attribute parameters Feature 4")
}



par(mfrow=c(1,2))

# example 2: analysis on determinants of anger-related behavior

# load anger data
data(anger)

# compute 1 run of a disjunctive model with 4 features
# use frequency data as input
angerdisj2<-plfm(maprule="disj",freq1=anger$freq1,freqtot=anger$freqtot,F=2,M=1)

# plot situation and behavior parameters
par(mfrow=c(2,2))
for (f in 1:2){
plot(angerdisj2,feature=f,element="object",main=paste("Situation parameters Feature",f,sep=" "))}
for (f in 1:2){
plot(angerdisj2,feature=f,element="attribute",main=paste("Behavior parameters Feature",f,sep=" "))}

}