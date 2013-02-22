\name{anger}
\docType{data}
\alias{anger}
\title{Situational determinants of anger-related behavior}
\description{
 	The raw data consist of the binary judgments of 101 first-year psychology students who indicated whether or not they would display 
	each of 8 anger-related behaviors when being angry at someone in each of 6 situations. 
	The 8 behaviors consist of 4 pairs of reactions that reflect a particular strategy to deal with situations in which 
	one is angry at someone, namely, (1) fighting (fly off the handle, quarrel), 
	(2) fleeing (leave, avoid), (3) emotional sharing (pour out one's heart, tell one's story), 
	and (4) making up (make up, clear up the matter). 
	The six situations are constructed from two factors  with three levels: 
	(1) the extent to which one likes the  instigator of anger (like, dislike, unfamiliar), 
	and (2) the status of the instigator of anger (higher, lower, equal). 
	Each situation is presented as one level of a factor, 
	without specifying a level for the other factor.
}
\usage{data(anger)}
\format{The data consist of a list of 5 objects:
\enumerate{
\item freq1: A 6 X 8 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates how many of 101 respondents 
 would display reaction  \emph{k} in situation \emph{j}.
\item freqtot: A 6 X 8 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the total number of respondents 
 who judged the situation-response pair \emph{(j,k)}.
\item rowlabels: A vector of labels for the situations.
\item columnlabels: A vector of labels for the anger-related reactions.
\item data: A 101 X 6 X 8 array of binary (0/1) values. the value in cell \emph{(i,j,k)} equals 1 if
 person \emph{i} would display behavior \emph{k} in situation \emph{j}, and 0 otherwise.
}
}
\source{Meulders, M., De Boeck, P., Kuppens, P., and Van Mechelen, I. (2002). Constrained latent
class analysis of three-way three-mode data. \emph{Journal of Classification, 19}, 277-302.
}
\references{
Kuppens, P., Van Mechelen, I., and Meulders, M. (2004). Every cloud has a silver lining: 
   Interpersonal and individual differences determinants of anger-related behaviors. \emph{Personality and Social Psychology Bulletin, 30}, 1550-1564. 


Meulders, M. (2011). \emph{An R package for probabilistic Latent feature analysis of two-way two-mode frequencies.} Manuscript submitted for publication.

Vermunt, J. K. (2007). A hierarchical mixture model for clustering three-way data
sets. \emph{Computational Statistics and Data Analysis, 51}, 5368-5376.

}
\keyword{datasets}