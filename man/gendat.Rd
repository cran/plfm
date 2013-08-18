\name{gendat}
\alias{gendat}
\title{Data generation}
\description{Computation of association probabilities and data generation for disjunctive or conjunctive probabilistic latent feature models.}
\usage{gendat(maprule="disj", N, objpar, attpar)}
\arguments{
  \item{maprule}{Disjunctive (\code{maprule="disj"}) or conjunctive (\code{maprule="conj"}) 
	mapping rule of the probabilistic latent feature model.}
  \item{N}{Number of replications for which binary associations are generated.}
  \item{objpar}{True objectparameters. As object parameters are probabilities they should be between 0 and 1.}
  \item{attpar}{True attributeparameters. As attribute parameters are probabilities they should be between 0 and 1.}
}
\details{
 The function \code{gendat} computes for all pairs of \emph{J} objects and \emph{K} attributes association probabilities and it generates 
 association frequencies (i.e. the number of replications \code{N} for which an object is associated to an attribute), 
 according to a disjunctive or a conjunctive probabilistic latent feature model. In addition, the function computes a matrix with in each cell 
 the total number of replications \code{N}. 
 If the requested number of replications \code{N} equals 0, 
 the function only computes association probabilities and does not generate new data. 
 
 To compute association probabilities the function \emph{gendat} uses a \emph{J X F} matrix of object parameters and a \emph{K X F} matrix 
 of attribute parameters as input. The \emph{F} object parameters of object \emph{j} represent, for each of \emph{F} features, 
 the probability that object \emph{j} has feature \emph{f}. 
 Similarly, the \emph{F} attribute parameters of attribute \emph{k} reflect, for each of \emph{F} features, 
 the probability that attribute \emph{k} is linked to feature \emph{f}. 

 According to the \emph{disjunctive} probabilistic latent feature model, object \emph{j} is associated
 to attribute \emph{k} if the object and the attribute have at least one feature in common. 
 More specifically, the association probability in cell \emph{(j,k)} for the disjunctive model 
 can be computed as: 

\deqn{p(j,k)=1-\prod_f(1-objpar[j,f]*attpar[k,f]).}
 
 According to the \emph{conjunctive} probabilistic latent feature model, object \emph{j} and attribute \emph{k} 
 are associated if object \emph{j} has all the features that are linked to attribute \emph{k}. 
 For the conjunctive model the association probability in cell \emph{(j,k)} is computed as:

\deqn{p(j,k)=\prod_f(1-(1-objpar[j,f])*attpar[k,f]).}
}
\value{
\item{call}{Parameters used to call the function.}
\item{prob1}{\emph{J X K} matrix of association probabilities.}
\item{freq1}{\emph{J X K} matrix of association frequencies.}
\item{freqtot}{\emph{J X K} matrix with number of replications.}
}
\references{
Maris, E., De Boeck, P., and Van Mechelen, I. (1996). Probability matrix decomposition models. \emph{Psychometrika, 61}, 7-29.

Meulders, M., De Boeck, P., Van Mechelen, I., Gelman, A., and Maris, E. (2001). Bayesian inference with probability matrix decomposition models. 
\emph{Journal of Educational and Behavioral Statistics, 26}, 153-179.

Meulders, M., De Boeck, P., Van Mechelen, I., & Gelman, A. (2005). Probabilistic feature analysis of facial perception of emotions. 
\emph{Applied Statistics, 54}, 781-793.
}
\author{Michel Meulders}
\seealso{\code{\link{plfm}}}
\examples{
## define constants
J<-20
K<-15
F<-2

## generate true parameters
set.seed(43565)
objectparameters<-matrix(runif(J*F),nrow=J)
attributeparameters<-matrix(runif(K*F),nrow=K)

## compute association probabilities for a conjunctive model
probconj<-gendat(maprule="conj",N=0,
             objpar=objectparameters,attpar=attributeparameters)

## generate data for a disjunctive model using N=200 replications
gdat<-gendat(maprule="disj",N=200,
             objpar=objectparameters,attpar=attributeparameters)
}
