\name{bayesplfm}
\alias{bayesplfm}
\title{Bayesian analysis of probabilistic latent feature models for two-way two-mode frequency data}
\description{Computation of a sample of the posterior distribution for disjunctive or conjunctive probabilistic latent feature models with \emph{F} features.}
\usage{
bayesplfm(data,object,attribute,rating,freq1,freqtot,F,
          Nchains=2,Nburnin=0,maxNiter=4000,
          Nstep=1000,Rhatcrit=1.2,maprule="disj",datatype="freq",
          start.bayes="best",fitted.plfm=NULL)

}
\arguments{
  \item{data}{ A data frame that consists of three components: the variables 
               \code{object}, \code{attribute} and \code{rating}. Each row of the data frame describes the outcome of a binary rater judgement
               about the association between a certain object and a certain attribute.} 
  \item{object}{The name of the \code{object} component in the data frame \code{data}. The values of the vector \code{data$object} should be (non-missing) numeric or character values.}
  \item{attribute}{The name of the \code{attribute} component in the data frame \code{data}. The values of the vector \code{data$attribute} should be (non-missing) numeric or character values.}
  \item{rating}{The name of the \code{rating} component in the data frame \code{data}. The elements of the vector \code{data$rating} should be the numeric values 0 (no association) or 1 (association), 
                or should be specified as missing (NA).}
  \item{freq1}{A \emph{J X K} matrix of observed association frequencies.}
  \item{freqtot}{A \emph{J X K} matrix with the total number of binary ratings in each cell \emph{(j,k)}. If the total number of ratings is the same for all cells of the matrix 
                 it is sufficient to enter a single numeric value rather than a matrix. For instance, if \emph{N} raters have judged \emph{J X K} associations, one may specify \code{freqtot}\emph{=N}}
  \item{F}{The number of latent features included in the model.}  
  \item{Nchains}{The number of Markov-chains that are simulated using a data-augmented Gibbs sampling algorithm.}
  \item{Nburnin}{The number of burn-in iterations.}
  \item{maxNiter}{The maximum number of iterations that will be computed for each chain.}
  \item{Nstep}{The convergence of the chains to the true posterior will be checked for each parameter after c*\code{Nstep} iterations with c=1,2,...
               The convergence will only be checked when \code{Nchains}>1.}
  \item{Rhatcrit}{The estimation procedure will be stopped if the Rhat convergence diagnostic is smaller than \code{Rhatcrit} 
                  for each object- and attribute parameter. By default \code{Rhatcrit}=1.2.}  
  \item{maprule}{Disjunctive (\emph{maprule="disj"}) or conjunctive (\emph{maprule="conj"}) mapping rule of the probabilistic latent feature model.}
  \item{datatype}{The type of data used as input. When \code{datatype}="freq" one should specify frequency data \code{freq1} and \code{freqtot}, and when \code{datatype}="dataframe" one should 
                  specify the name of the data frame \code{data}, and its components, \code{object}, \code{attribute} and \code{rating}.}
  \item{start.bayes}{This argument can be used to define the type of starting point for the Bayesian analysis. If \code{start.bayes}="best" the best solution of a \code{\link{plfm}} analysis 
                     is used as the starting point for the Bayesian analysis, and if \code{start.bayes} = "fitted.plfm", the
                     starting point is read from the (\code{\link{plfm}}) object assigned to the argument 
                     \code{fitted.plfm}. If \code{start.bayes}="random", a random starting point is used for the Bayesian analysis.
                     }
  \item{fitted.plfm}{The name of the \code{\link{plfm}} object that contains posterior mode estimates for the specified model.} 


}
\details{
  The function \code{\link{bayesplfm}} can be used to compute a sample of the posterior 
  distribution of disjunctive or conjunctive probabilistic latent feature models with a particular number of features 
  using a data-augmented Gibbs sampling algorithm 
  (Meulders, De Boeck, Van Mechelen, Gelman, and Maris, 2001; Meulders, De Boeck, Van Mechelen, and Gelman, 2005; Meulders, 2013).

  By specifying the parameter \code{Nchains} the function can be used to compute one single chain, or multiple chains. 
  When only one chain is computed, no convergence measure is reported. When more than one chain is computed, for each parameter, 
  convergence to the true posterior distribution is assessed using the Rhat convergence diagnostic proposed by Gelman and Rubin (1992).
  
  When using \code{\link{bayesplfm}} for Bayesian analysis the same starting point will be used for each simulated chain. The reason for using the same
  starting point for each of the chains is that the posterior distribution of probabilistic feature models with \emph{F>2} is always multimodal 
  (local maxima may exist, and one may switch feature labels), so that the aim of the Bayesian analysis is to compute a sample in the neigbourhood 
  of one specific posterior mode.  It is recommended to use the best posterior mode obtained
  with the \code{\link{plfm}} function as a starting point for the Bayesian analysis (use \code{start.bayes}="best", or specify \code{start.bayes}="fitted.plfm" and 
  \code{fitted.plfm}=object) with "object" being a \code{\link{plfm}} object that contains posterior mode estimates for the specified model. As an alternative to using the \code{plfm()}, 
   function one may use random  starting points for the Bayesian analysis (\code{start.bayes}="random") to explore the posterior distribution.

  The function \code{bayesplfm()} will converge well if the  distinct posterior modes are well-separated and if the different chains only visit the same mode during the estimation process. 
  However, if the posterior distribution is multimodal, it may fail to converge  if the Gibbs sampler starts visiting different posterior modes within
  one chain, or if different chains sample from distinct posterior modes. 
 
}
\value{
\item{call}{Parameters used to call the function.}
\item{sample.objpar}{A \emph{J X F X Niter X Nchains} array with parameter values for the object parameters. 
	The matrix \code{sample.objpar[,,i,c]} contains the draw of the object parameters 
	in iteration \emph{i} of chain \emph{c}. Note: when \code{Nchains}=1 the chain length \emph{Niter} equals \code{maxNiter}, 
        and when \code{Nchains}>1 the chain length \emph{Niter} equals the number of iterations required to obtain convergence.}
\item{sample.attpar}{A \emph{K X F X Niter X Nchains} array with parameter values for the attribute parameters. 
	              The matrix \code{sample.attpar[,,i,c]} contains the draw of the attribute parameters 
	              in iteration \emph{i} of chain \emph{c}. Note: when \code{Nchains}=1 the chain length \emph{Niter} equals \code{maxNiter}, 
                      and when \code{Nchains}>1 the chain length \emph{Niter} equals the number of iterations required to obtain convergence.}
\item{pmean.objpar}{A \emph{J X F} matrix with the posterior mean of the object parameters computed on all iterations and chains in the sample.}
\item{pmean.attpar}{A \emph{K X F} matrix with the posterior mean of the attribute parameters computed on all iterations and chains in the sample.}
\item{p95.objpar}{A \emph{3 X J X F} array which contains for each object parameter the percentiles 2.5, 50 and 97.5.}
\item{p95.attpar}{A \emph{3 X K X F} array which contains for each attribute parameter the percentiles 2.5, 50 and 97.5.}
\item{Rhat.objpar}{A \emph{J X F} matrix with Rhat convergence values for the object parameters.}
\item{Rhat.attpar}{A \emph{K X F} matrix with Rhat convergence values for the attribute parameters.}
\item{fitmeasures}{A list with two measures of descriptive fit on the \emph{J X K} table: (1) the correlation between observed and expected frequencies, 
                   and (2) the proportion of the variance in the observed frequencies accounted for by the model. 
                   The association probabilities and corresponding expected frequencies are computed using the posterior mean of the parameters.}
\item{convstat}{The number of object-and attribute parameters that do not meet the convergence criterion.}
}
\references{

Gelman, A., and Rubin, D. B. (1992). Inference from iterative simulation using multiple
sequences. \emph{Statistical Science, 7} , 457-472.

Meulders, M., De Boeck, P., Van Mechelen, I., Gelman, A., and Maris, E. (2001). Bayesian inference with probability matrix decomposition models. 
\emph{Journal of Educational and Behavioral Statistics, 26}, 153-179.


Meulders, M., De Boeck, P., Van Mechelen, I., and Gelman, A. (2005). Probabilistic feature analysis of facial perception of emotions. 
\emph{Applied Statistics, 54}, 781-793.

Meulders, M. and De Bruecker, P. (2018). Latent class probabilistic latent feature analysis of three-way three-mode binary data. 
\emph{Journal of Statistical Software, 87(1)}, 1-45.

Meulders, M. (2013). An R Package for Probabilistic Latent Feature Analysis of Two-Way Two-Mode Frequencies. \emph{Journal of Statistical Software, 54(14)}, 1-29. 
   URL http://www.jstatsoft.org/v54/i14/.

}
\author{
Michel Meulders
}
\seealso{
\code{\link{plfm}}, \code{\link{summary.bayesplfm}},\code{\link{print.summary.bayesplfm}}
}
\examples{

\dontrun{
## example 1: Bayesian analysis using data generated under the model

## define number of objects
J<-10
## define number of attributes
K<-10
## define number of features
F<-2

## generate true parameters
set.seed(43565)
objectparameters<-matrix(runif(J*F),nrow=J)
attributeparameters<-matrix(runif(K*F),nrow=K)

## generate data for conjunctive model using N=100 replications
gdat<-gendat(maprule="conj",N=100,
              objpar=objectparameters,attpar=attributeparameters)

## Use stepplfm to compute posterior mode(s) for 1 up to 3 features 

conj.lst<-stepplfm(minF=1,maxF=3,maprule="conj",freq1=gdat$freq1,freqtot=100,M=5)


## Compute a sample of the posterior distribution 
## for the conjunctive model with two features
## use the posterior mode obtained with stepplfm as starting point
conjbayes2<-bayesplfm(maprule="conj",freq1=gdat$freq1,freqtot=100,F=2,
                      maxNiter=3000,Nburnin=0,Nstep=1000,Nchains=2,
                      start.bayes="fitted.plfm",fitted.plfm=conj.lst[[2]])

}

\dontrun{
## example 2: Bayesian analysis of situational determinants of anger-related behavior

## load data
data(anger)

## Compute one chain of 500 iterations (including 250 burn-in iterations) 
## for the disjunctive model with two features
## use a random starting point

bayesangerdisj2a<-bayesplfm(maprule="disj",freq1=anger$freq1,freqtot=anger$freqtot,F=2,
                      maxNiter=500,Nstep=500,Nburnin=250,Nchains=1,start.bayes="random")

##print a summary of the output 
summary(bayesangerdisj2a)


## Compute a sample of the posterior distribution 
## for the disjunctive model with two features
## compute starting points with plfm
## run 2 chains with a maximum length of 10000 iterations
## compute convergence after each 1000 iterations

bayesangerdisj2b<-bayesplfm(maprule="disj",freq1=anger$freq1,freqtot=anger$freqtot,F=2,
                      maxNiter=10000,Nburnin=0,Nstep=1000,Nchains=2,start.bayes="best")


## print the output of the disjunctive 2-feature model for the anger data
print(bayesangerdisj2b)


## print a summary of the output of the disjunctive 2-feature model 
##for the anger data
summary(bayesangerdisj2b)
}


}
