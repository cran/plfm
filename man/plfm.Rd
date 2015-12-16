\name{plfm}
\alias{plfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Probabilistic latent feature analysis of two-way two-mode frequency data}
\description{Computation of parameter estimates, standard errors, criteria for model selection, and goodness-of-fit criteria for 
 disjunctive or conjunctive probabilistic latent feature models with \emph{F} features.}
\usage{
plfm(data,object,attribute,rating,freq1,freqtot,F,
     datatype="freq",maprule="disj",M=5,emcrit1=1e-2,
     emcrit2=1e-10,printrun=TRUE)
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
  \item{datatype}{The type of data used as input. When \code{datatype}="freq" one should specify frequency data \code{freq1} and \code{freqtot}, and when \code{datatype}="dataframe" one should 
                  specify the name of the data frame \code{data}, and its components, \code{object}, \code{attribute} and \code{rating}.}
  \item{maprule}{Disjunctive (\code{maprule}="disj") or conjunctive (\code{maprule}="conj") mapping rule of the probabilistic latent feature model.}
  \item{M}{The number of times a particular model is estimated using random starting points.}
  \item{emcrit1}{Convergence criterion which indicates when the estimation algorithm should switch from Expectation-Maximization (EM) steps to EM+Newton-Rhapson steps.}
  \item{emcrit2}{Convergence criterion which indicates final convergence to a local maximum.}
  \item{printrun}{\code{printrun}=TRUE prints the analysis type (disjunctive or conjunctive), the number of features (\emph{F}) and the number of the run to the output screen, whereas 
                 \code{printrun}=FALSE suppresses the printing.}
}
\details{
\emph{Estimation}

	The function \code{\link{plfm}} uses an accelerated EM-algorithm to locate the posterior mode(s) of the probabilistic latent feature model. 
	The algorithm starts with a series of Expectation-Maximization (EM) steps until the difference between subsequent values of the 
	logarithm of the posterior density becomes smaller than the convergence criterion \code{emcrit1}, and then switches to an accelerated algorithm which 
	consists of EM + Newton-Rhapson steps. The accelerated algorithm stops when  the difference between subsequent values of the 
	logarithm of the posterior density becomes smaller than the convergence criterion \code{emcrit2}.  
	Computational details about the implementation of the EM-steps for PLFMs are described in Maris, De Boeck, and Van Mechelen (1996). 
	The general scheme of the accelerated algorithm is described in Louis (1982) and Tanner (1996). Computational details about implementing the accelerated algorithm 
	for PLFMs are described in Meulders (2013).

	When using the function \code{\link{plfm}} to estimate a particular PLFM (i.e. with a certain number of latent features and a disjunctive or conjunctive mapping rule), 
	one may locate the distinct posterior mode(s) by running the algorithm \code{M} times using random starting points. The estimated object-and attribute parameters of 
	each run are stored in the objects \code{objpar.runs} and \code{attpar.runs} of the output list. Next, a number of additional statistics 
	(estimated object- and attribute parameters, asymptotic standard errors of object- and attribute parameters, model selecion criteria and goodness-of-fit measures) 
	are computed for the best model (i.e. the model among \code{M} runs with the highest posterior density).

\emph{Model selection criteria and goodness-of-fit measures}

        To choose among models with different numbers of features, or with different mapping rules, 
	one may use information criteria such as the Akaike Information Criterion (AIC, Akaike, 1973, 1974),
	or the Schwarz Bayesian Information Criterion (BIC, Schwarz, 1978). AIC and BIC are computed as \emph{-2*loglikelihood+k*Npar}.
	For AIC \emph{k} equals 2 and for BIC \emph{k} equals \emph{log(N)}, with \emph{N} the observed number of replications 
	for which object-attribute associations are collected. \emph{Npar} represents the number of model parameters; 
	for probabilistic latent feature models this equals \emph{(J+K)F}. Models with the lowest value for AIC or BIC should be selected.
	
	To assess the statistical fit of the probabilistic feature model one may use a Pearson chi-square measure 
	on the \emph{J X K} frequency table to evaluate whether predicted 
	frequencies deviate significantly from observed frequencies (see Meulders et al., 2001). In addition, one may assess the descriptive fit of the model 
	using the correlation between observed and expected frequencies in the \emph{J X K} table, 
	and the proportion of the variance in the observed frequencies accounted for by the model (VAF)
	(i.e. the squared correlation between observed and expected frequencies).

	The model selection criteria AIC and BIC, the results of the Pearson goodness-of fit test, 
	and the descriptive fit measures (correlation observed and expected frequencies, and VAF) are stored in the object \code{fitmeasures} of the output list

}
\value{
\item{call}{Parameters used to call the function.}
\item{objpar}{A \emph{J X F} matrix of object parameters.}
\item{attpar}{A \emph{K X F} matrix of attribute parameters.}
\item{fitmeasures}{A list of model selection criteria and goodness-of-fit criteria for the model with the highest posterior density.}
\item{logpost.runs}{A list with the logarithm of the posterior density for each of the \emph{M} computed models.}
\item{objpar.runs}{A \emph{M X J X F} array which contains the object parameters for each of the \emph{M} computed models.}
\item{attpar.runs}{A \emph{M X K X F} array which contains the attribute parameters for each of the \emph{M} computed models.}
\item{bestsolution}{An index which indicates the model with the highest posterior density among each of the \emph{M} computed models.}
\item{gradient.objpar}{A \emph{J X F} gradient matrix for the object parameters in the best solution.}
\item{gradient.attpar}{A \emph{K X F} gradient matrix for the attribute parameters in the best solution.}
\item{SE.objpar}{A \emph{J X F} matrix of asymptotic standard errors for the object parameters in the best solution.}
\item{SE.attpar}{A \emph{K X F} matrix of asymptotic standard errors for the attribute parameters in the best solution.}
\item{prob1}{A \emph{J X K} matrix of expected association probabilities for the best solution.}
}
\references{

Akaike, H. (1973). Information theory and an extension of the maximum likelihood
principle. In B. N. Petrov and F. Csaki (Eds.), \emph{Second international symposium on
information theory} (p. 271-283). Budapest: Academiai Kiado.

Akaike, H. (1974). A new look at the statistical model identification. \emph{IEEE Transactions
on Automatic Control, 19}, 716-723.

Candel, M. J. J. M., and Maris, E. (1997). Perceptual analysis of two-way two-mode
frequency data: probability matrix decomposition and two alternatives.
\emph{International Journal of Research in Marketing, 14}, 321-339.

Louis, T. A. (1982). Finding observed information using the em algorithm. \emph{Journal of the
Royal Statistical Society, Series B, 44}, 98-130.


Maris, E., De Boeck, P., and Van Mechelen, I. (1996). Probability matrix decomposition models. \emph{Psychometrika, 61}, 7-29.

Meulders, M., De Boeck, P., and Van Mechelen, I. (2001). Probability matrix decomposition
models and main-effects generalized linear models for the analysis of replicated
binary associations. \emph{Computational Statistics and Data Analysis, 38}, 217-233.

Meulders, M., De Boeck, P., Van Mechelen, I., Gelman, A., and Maris, E. (2001). Bayesian inference with probability matrix decomposition models. 
\emph{Journal of Educational and Behavioral Statistics, 26}, 153-179.

Meulders, M., De Boeck, P., Van Mechelen, I., and Gelman, A. (2005). Probabilistic feature analysis of facial perception of emotions. 
\emph{Applied Statistics, 54}, 781-793.

Meulders, M. (2013). An R Package for Probabilistic Latent Feature Analysis of Two-Way Two-Mode Frequencies. \emph{Journal of Statistical Software, 54(14)}, 1-29. 
   URL http://www.jstatsoft.org/v54/i14/.

Tanner, M. A. (1996). \emph{Tools for statistical inference: Methods for the exploration of
posterior distributions and likelihood functions} (Third ed.). New York:
Springer-Verlag.

Schwarz, G. (1978). Estimating the dimensions of a model. \emph{Annals of Statistics, 6}, 461-464.


}
\author{
Michel Meulders
}
%%\note{
%%  ~~further notes~~
%%}
\seealso{\code{\link{gendat}}, \code{\link{print.plfm}}, 
          \code{\link{summary.plfm}}, \code{\link{print.summary.plfm}}}
\examples{

\dontrun{
## example 1: Analysis of data generated under the model

# define constants
J<-20
K<-15
F<-2

# generate true parameters
set.seed(43565)
objectparameters<-matrix(runif(J*F),nrow=J)
attributeparameters<-matrix(runif(K*F),nrow=K)

# generate data for disjunctive model using N=200 replications
gdat<-gendat(maprule="disj",N=200,
              objpar=objectparameters,attpar=attributeparameters)

# Estimate a disjunctive probabilistic latent feature model with 2 features
# suppress printing the type of analysis to the output screen
disj2<-plfm(maprule="disj",freq1=gdat$freq1,freqtot=200,F=2,M=1,printrun=FALSE)
}


\dontrun{
# example 2:Perceptual analysis of associations between car models and car attributes

# load car data
data(car)

# compute 1 run of a disjunctive model with 4 features
# use components of a data frame as input
cardisj4<-plfm(datatype="dataframe",data=car$datalongformat,object=objectlabel,
               attribute=attributelabel,rating=rating,maprule="disj",F=4,M=1)

# print the output of a disjunctive 4-feature model  
# for data on the perception of car models
print (cardisj4)


# print a summary of the output of a disjunctive 4-feature model  
# for data on the perception of car models
sumcardisj4<-summary(cardisj4)
sumcardisj4
}

\dontrun{
# example 3: analysis on determinants of anger-related behavior

# load anger data
data(anger)

# compute 1 run of a disjunctive model with 2 features
# use frequency data as input
angerdisj2<-plfm(maprule="disj",freq1=anger$freq1,freqtot=anger$freqtot,F=2,M=1)


# print the output of a disjunctive 2-feature model 
# for data on the situational determinants of anger-related behaviors
print (angerdisj2)


# print a summary of the output of a disjunctive 2-feature model 
# for data on the situational determinants of anger-related behaviors
sumangerdisj2<-summary(angerdisj2)
sumangerdisj2
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
