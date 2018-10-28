\name{anger2}
\docType{data}
\alias{anger2}
\title{Situational determinants of anger-related behavior}
\description{
 	The raw data consist of the binary judgments of 115 first-year psychology students who indicated whether or not they would display 
	each of 14 anger-related behaviors when being angry at someone in each of 9 situations. 
	The 14 behaviors consist of 7 pairs of reactions that reflect a particular strategy to deal with situations in which 
	one is angry at someone:
\enumerate{
\item Anger-out: (a) You flew off the handle, (b) You started a fight
\item Avoidance: (a) You avoided a confrontation, (b) You went out of the other's way
\item Social sharing (a) You unburdened your heart to others, (b) You told others what had happened 
\item Assertive behavior: (a) You said what was bothering you in a direct and sober way, (b) You calmly explained what was bothering you 
\item Indirect behavior (a) You showed something was bothering you without saying anything, (b) You started to sulk 
\item Anger-in: (a) You suppressed your anger, (b) You bottled up your anger 
\item Reconciliation (a) You reconciled, (b) You talked things out}

	The six situations are constructed by crossing the levels of two factors  with three levels: 
	(1) the extent to which one likes the  instigator of anger (like, unfamilar, dislike), 
	and (2) the status of the instigator of anger (lower status, equal status, higher status)
}
\usage{data(anger2)}
\format{The data consist of a list of 5 objects:
\enumerate{
\item data: A 115 X 9 X 14 matrix of binary observations (0/1). The observation in cell \emph{(i,j,k)} equals 1 if person \emph{i} would display behavior
\emph{k} in situation \emph{j} and 0 otherwise.
\item freq1: A 9 X 14 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the number of respondents 
 who indicate that they would display behavior \emph{k} in situation \emph{j}.
\item freqtot: A 9 X 14 matrix of frequencies. The frequency in cell \emph{(j,k)} indicates the total number of respondents 
 who judged the situation-response pair \emph{(j,k)}.
\item rowlabels: A vector of labels for the situations.
\item columnlabels: A vector of labels for the anger-related behaviors.
}
}
\source{Kuppens, P., Van Mechelen, I., and Meulders, M. (2004). Every cloud has a silver lining: 
   Interpersonal and individual differences determinants of anger-related behaviors. \emph{Personality and Social Psychology Bulletin, 30}, 1550-1564.}
\references{Meulders, M. and De Bruecker, P. (2018). Latent class probabilistic latent feature analysis of three-way three-mode binary data. 
  \emph{Journal of Statistical Software, 87(1), 1-45}, 1-29.
}
\keyword{datasets}