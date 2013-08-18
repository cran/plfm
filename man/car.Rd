\name{car}
\docType{data}
\alias{car}
\title{Ratings of associations between car models and car attributes}
\description{ The data describe the ratings of 78 respondents about the association between each of 14 car models and each of 27 car attributes.}
\usage{data(car)}
\format{The data consist of a list of 3 objects:
\enumerate{
\item data: A data frame that consists of 6 components:  Each row of the data frame  describes the outcome of a binary rater judgement
            about the association between a certain car and a certain attribute. The components \emph{IDobject} and \emph{objectlabel} contain an ID and 
             label for the car models, the components \emph{IDattribute} and \emph{attributelabel} contain an ID and  label for the attributes, the component \emph{IDrater} 
            contains a rater ID, and the component \emph{rating} contains the binary  ratings (1 if the car model has the attribute according to the rater, and 0 otherwise).  
\item freq1: A 14 X 27 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates how many of 78 respondents 
 indicate an association between car model \emph{j} and attribute \emph{k}.
\item freqtot: A 14 X 27 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the total number of respondents 
 who judged the car-attribute pair \emph{(j,k)}.
}
}
\source{Van Gysel, E. (2011). \emph{Perceptuele analyse van automodellen met probabilistische feature modellen.} [translation from Dutch: Perceptual analysis of car models with probabilistic feature models] 
        Master thesis. Hogeschool-Universiteit Brussel.}
\references{Meulders, M. (2013). An R Package for Probabilistic Latent Feature Analysis of Two-Way Two-Mode Frequencies. \emph{Journal of Statistical Software, 54(14)}, 1-29. 
   URL http://www.jstatsoft.org/v54/i14/.
}
\keyword{datasets}