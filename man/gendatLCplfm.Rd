\name{gendatLCplfm}
\alias{gendatLCplfm}
\title{Data generation}
\description{Data generation for disjunctive or conjunctive latent class probabilistic latent feature models.}
\usage{gendatLCplfm(N,objpar,attpar,sizepar,maprule="disj",model=1)}
\arguments{
  \item{N}{Number of replications (e.g. persons) for which binary object-attribute associations are generated.}
  \item{objpar}{True objectparameters. If \code{model}=1, \code{model}=3, \code{model}=4, or \code{model}=6 \code{objpar} is a \emph{J X F X T} array, 
                if  \code{model}=2 or \code{model}=5 \code{objpar} is a \emph{J X F} matrix. As object parameters are probabilities they should be between 0 and 1.}
  \item{attpar}{True attributeparameters. If \code{model}=2, \code{model}=3, \code{model}=5, or \code{model}=6 \code{attpar} is a \emph{K X F X T} array, 
                if  \code{model}=1 or \code{model}=4 \code{attpar} is a \emph{K X F} matrix. As attribute parameters are probabilities they should be between 0 and 1.}
  \item{sizepar}{A \emph{T}-vector of true class size parameters.}
  \item{maprule}{Disjunctive (\code{maprule="disj"}) or conjunctive (\code{maprule="conj"}) 
	        mapping rule of the latent class probabilistic latent feature model.}
  \item{model}{The type of dependency and heterogeneity assumption included in the model. \code{model}=1, \code{model}=2, \code{model}=3 represent models with a constant 
               object-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters. \code{model}=4, \code{model}=5, \code{model}=6 represent models with a constant 
               attribute-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters.}
}
\details{
 The function \code{gendatLCplfm} generates binary object-attribute associations for \emph{N} replications 
  according to a disjunctive or a conjunctive latent class probabilistic latent feature model of a specific model type. 
  In addition, the function computes the \emph{J X K} matrix of marginal object-attribute association probabilities and
  a  \emph{J X K X T} array of class-specific object-attribute association probabilities.                      
  To compute association probabilities the function \emph{gendatLCplfm} uses a vector of class size parameters (\code{sizepar}) a matrix or array of 
  object parameters (\code{objpar}) and a matrix or array of true  attribute parameters (\code{attpar}) as input. 
 
 According to the \emph{disjunctive} probabilistic latent feature model, object \emph{j} is associated
 to attribute \emph{k} if the object and the attribute have at least one feature in common. 
 More specifically, for \code{model=1} the class-specific object-attribute association probability in cell \emph{(j,k)} for the disjunctive model 
 can be computed as: 

\deqn{p(j,k|t)=1-\prod_f(1-objpar[j,f,t]*attpar[k,f]).}
 
 According to the \emph{conjunctive} probabilistic latent feature model, object \emph{j} and attribute \emph{k} 
 are associated if object \emph{j} has all the features that are linked to attribute \emph{k}. 

 In particular, for \code{model=1},the class-specific object-attribute association probability in cell \emph{(j,k)} for the conjunctive model can be computed as:

\deqn{p(j,k|t)=\prod_f(1-(1-objpar[j,f,t])*attpar[k,f]).}

The marginal object-attribute association probability can be computed as follows:

\deqn{p(j,k)=\sum_t sizepar[t]*p(j,k|t).}
}
\value{
\item{call}{Parameters used to call the function.}
\item{data}{\emph{I X J X K} matrix of association probabilities.}
\item{class}{\emph{I}-vector that contains latent class membership of each replication.}
\item{condprob.JKT}{\emph{J X K X T} array of class-specific conditional object-attribute association probabilities.}
\item{margprob.JK}{\emph{J X K} matrix  of marginal object-attribute association probabilities.}
}
\references{
Meulders, M., Tuerlinckx, F., and Vanpaemel, W. (2013). Constrained multilevel latent class models for the analysis of three-way three-mode binary data. 
\emph{Journal of Classification, 30 (3)}, 306-337.
}
\author{Michel Meulders}
\seealso{\code{\link{LCplfm}}}
\examples{
\dontrun{
# define constants
I<-500
J<-10
K<-8
F<-2
T<-2

# model 1

# generate true parameters
objpar<-array(runif(J*F*T),c(J,F,T))
attpar<-matrix(runif(K*F),c(K,F))
sizepar<-rep(1/T,T)
# generate data
d<-gendatLCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=1)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=1,maprule="conj")


# model 2

# generate true parameters
objpar<-matrix(runif(J*F),nrow=J)
attpar<-array(runif(K*F*T),c(K,F,T))
sizepar<-rep(1/T,T)
# generate data
d<-gendatLCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=2)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=2,maprule="conj")

# model 3

# generate true parameters
objpar<-array(runif(J*F*T),c(J,F,T))
attpar<-array(runif(K*F*T),c(K,F,T))
sizepar<-rep(1/T,T)
# generate data
d<-gendat.LCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=3)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=3,maprule="conj")

# model 4

# generate true parameters
objpar<-array(runif(J*F*T),c(J,F,T))
attpar<-matrix(runif(K*F),c(K,F))
sizepar<-rep(1/T,T)
# generate data
d<-gendat.LCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=4)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=4,maprule="conj")

# model 5

# generate true parameters
objpar<-matrix(runif(J*F),nrow=J)
attpar<-array(runif(K*F*T),c(K,F,T))
sizepar<-rep(1/T,T)
# generate data
d<-gendat.LCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=5)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=5,maprule="conj")


# model 6
# generate true parameters
objpar<-array(runif(J*F*T),c(J,F,T))
attpar<-array(runif(K*F*T),c(K,F,T))
sizepar<-rep(1/T,T)
# generate data
d<-gendatLCplfm(N=I,objpar=objpar,attpar=attpar,sizepar=sizepar,maprule="conj",model=6)
# estimate parameters of true model
res<-LCplfm(data=d$data,F=2,T=2,model=6,maprule="conj")
}
}
