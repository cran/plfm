\name{hostility}
\docType{data}
\alias{hostility}
\title{self-reported hostile behavior in frustrating situations}
\description{
 	The data consist of the judgments of 316 first-year psychology students who indicated on a three point scale the extent to which they would display 
	each of 4 hostile behaviors in each of 14 frustrating situations (0= you do not display this response in this situation, 
        1= you display this response to a limited extent in this situation, 2= you display this response to a strong extent in this situation). 
}
\usage{data(hostility)}
\format{The data consist of a list of 6 objects:
\enumerate{
\item data: A 316 X 14 X 4 array of dichotomized judgements (0 versus 1 or 2). The observation in cell \emph{(i,j,k)} equals 1 if person \emph{i} 
 would display behavior \emph{k} in situation \emph{j} to a limited or strong extent and 0 if person \emph{i} 
 would not display behavior \emph{k} in situation \emph{j}.
\item freq1: A 14 X 4 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the number of respondents 
 who indicate that they would display behavior \emph{k} in situation \emph{j}.
\item freqtot: A 14 X 4 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the total number of respondents 
 who judged the situation-response pair \emph{(j,k)}.
\item situation: A vector with descriptions of the situations.
\item rowlabels: A vector of labels for the situations.
\item columnlabels: A vector of labels for the anger-related behaviors.
}
}
\source{Vansteelandt, K. (1999). A formal model for the competency-demand hypothesis. \emph{European Journal of Personality, 13}, 429-442.
}
\references{Vansteelandt, K. and Van Mechelen, I. (1998). Individual differences in situation-behavior profiles: A triple typology model. 
\emph{Journal of Personality and Social Psychology, 75}, 751-765. 
}
\keyword{datasets}