\name{car2}
\docType{data}
\alias{car2}
\title{Judgements on associations between car models and car attributes}
\description{ The data consist of the binary judgements of 147 respondents about the association between each of 12 car models and each of 23 car attributes.}
\usage{data(car2)}
\format{The data consist of a list of 5 objects:
\enumerate{
\item data3w: A 147 X 12 X 23 array of binary judgements. The observation in cell \emph{(i,j,k)} equals 1 if respondent \emph{i} indicates that car \emph{j} has 
         attribute \emph{k}, and 0 otherwise.   
\item freq1: A 12 X 23 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates how many of 147 respondents 
 indicate an association between car model \emph{j} and attribute \emph{k}.
\item freqtot: A 12 X 23 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the total number of respondents 
 who judged the car-attribute pair \emph{(j,k)}.
}
}
\source{M\"ahler, R. (2014). \emph{Analyse van perceptie en preferentie van middelgrote wagens.} [translation from Dutch: Analysis of perception and preference for midsize cars.] 
        Master thesis. KU Leuven.} 
\keyword{datasets}