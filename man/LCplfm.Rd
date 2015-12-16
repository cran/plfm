\name{LCplfm}
\alias{LCplfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Latent class probabilistic feature analysis of three-way three-mode binary data}
\description{Computation of parameter estimates, standard errors, criteria for model selection, and measures of descriptive fit for 
 disjunctive or conjunctive latent class probabilistic  feature models.}
\usage{  
LCplfm(data,F=2,T=2,M=5,maprule="disj",emcrit1=1e-3,emcrit2=1e-8,
       model=1,start.objectparameters=NULL,start.attributeparameters=NULL,
       start.sizeparameters=NULL,delta=0.0001,printrun=FALSE,
       update.objectparameters=NULL,update.attributeparameters=NULL)
}
\arguments{
  \item{data}{A \emph{I X J X K} data array of binary observations. Observation \emph{(i,j,k)} \emph{(i=1,..,I; j=1,..,J; k=1,..,K)}  indicates 
               whether object \emph{j} is associated to attribute \emph{k} according to rater \emph{i}.} 
  \item{F}{The number of latent features included in the model.}
  \item{T}{The number of latent classes included in the model.}
  \item{M}{The number of times a particular model is estimated using random starting points.}
  \item{maprule}{Disjunctive (\code{maprule}="disj") or conjunctive (\code{maprule}="conj") mapping rule of the probabilistic latent feature model.}
  \item{emcrit1}{Convergence criterion to be used for the estimation of candidate models.}
  \item{emcrit2}{Convergence criterion to be used for the estimation of the best model.}
  \item{model}{The type of dependency and heterogeneity assumption included in the model. \code{model}=1, \code{model}=2, \code{model}=3 represent models with a constant 
               object-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters. \code{model}=4, \code{model}=5, \code{model}=6 represent models with a constant 
               attribute-feature classification per person and with, respectively, class-specific object parameters, class-specific attribute parameters, 
               and class-specific object- and attribute parameters.}
  \item{start.objectparameters}{An array of object parameters to be used as starting value for each run.
                                The size of the array equals \emph{J x F x T x M}  when \code{model = 1,4,3,6} and 
                                   \emph{J x F x M} when \code{model = 2,5}.  
                                 If \code{start.objectparameters=NULL} randomly generated object parameters 
                                 are used as starting values.}
  \item{start.attributeparameters}{An  array of attribute parameters to be used as starting value for each run. 
                                   The size of the array equals \emph{K x F x T x M}  when \code{model = 2,5,3,6} 
                                   and \emph{K x F x M} when \code{model = 1,3}. 
                                   If \code{start.attributeparameters=NULL} randomly generated attribute parameters 
                                   are used as starting values.}
  \item{start.sizeparameters}{A \emph{T x M} matrix of latent class size parameters to be used as starting value for each run. 
                              If \code{start.sizeparameters=NULL} randomly  generated class size parameters are used 
                              as starting values.}
  \item{delta}{The precision used to compute standard errors of the parameters with the method of finite differences.}
  \item{printrun}{\code{printrun}=TRUE prints the analysis type (disjunctive or conjunctive), the number of features (\emph{F}), the number of latent classes (\emph{T}) 
                  and the number of the run to the output screen, whereas \code{printrun}=FALSE suppresses the printing.}
  \item{update.objectparameters}{A binary valued array that indicates for each object parameter whether it has to be 	estimated 
                                 from the data or constrained to the starting value. A value of 1 means that the corresponding object parameter is estimated and 
                                 a value of 0 means that the corresponding object parameter is constrained to the starting value provided by the user.
                                 The size of the array equals \emph{J x F x T}  when \code{model = 1,4,3,6} and \emph{J x F} when \code{model = 2,5}.  
                                 If \code{update.objectparameters} \code{= NULL} all object parameters are estimated from the data.}
  \item{update.attributeparameters}{A binary valued array that indicates for each attribute parameter whether it has to be estimated 
                                 from the data or constrained to the starting value. A value of 1 means that the corresponding attribute parameter is estimated and 
                                 a value of 0 means that the corresponding attribute parameter is constrained to the starting value provided by the user.
                                 The size of the array equals \emph{K x F x T}  when \code{model = 2,5,3,6} and \emph{K x F} when \code{model = 1,3}.  
                                 If \code{update.attributeparameters} \code{= NULL} all attribute parameters are estimated from the data.}

}
\details{
\emph{Estimation}
      The estimation algorithm includes two steps. In a first exploratory step an EM algorithm is used to conduct \code{M} runs using random starting points. Each exploratory run is 
      terminated if the convergence criterium (i.e., the sum of absolute differences between parameter values in subsequent iterations) is smaller than \code{emcrit1}. In a second step, 
      the best solution  among the \code{M} runs (i.e., with the highest posterior density) is used as the starting point of the EM algorithm for conducting a final analysis. The final 
      analysis is terminated if the convergence criterion \code{emcrit2} is smaller than the convergence criterion. 

\emph{Model selection criteria and goodness-of-fit measures}

        To choose among models with different numbers of features, or with different mapping rules, 
	one may use information criteria such as the Akaike Information Criterion (AIC, Akaike, 1973, 1974),
	or the Schwarz Bayesian Information Criterion (BIC, Schwarz, 1978). AIC and BIC are computed as \emph{-2*loglikelihood+k*Npar}.
	For AIC \emph{k} equals 2 and for BIC \emph{k} equals \emph{log(N)}, with \emph{N} the observed number of replications (\emph{I}) 
	for which object-attribute associations are collected. \emph{Npar} represents the number of model parameters.
	Models with the lowest value for AIC or BIC should be selected.
	
	The descriptive fit of the model is assessed with the correlation between observed and expected frequencies in the \emph{J X K} table, 
	and the proportion of the variance in the observed frequencies accounted for by the model (VAF)
	(i.e. the squared correlation between observed and expected frequencies).

	The model selection criteria AIC and BIC, and the descriptive fit measures (correlation observed and expected frequencies, and VAF) 
        are stored in the object \code{fitmeasures} of the output list

}
\value{
\item{call}{Parameters used to call the function.}
\item{logpost.runs}{A list with the logarithm of the posterior density for each of the \emph{M} computed models.}
\item{best}{An index which indicates the model with the highest posterior density among each of the \emph{M} computed models.}
\item{objpar}{Estimated object parameters for the best model.}
\item{attpar}{Estimated attribute parameters for the best model.}
\item{sizepar}{A vector of \code{T} class size parameters for the best model.}
\item{SE.objpar}{Estimated standard errors for the object parameters of the best model.}
\item{SE.attpar}{Estimated standard errors for the attribute parameters of the best model.}
\item{SE.sizepar}{Estimated standard errors for the class size parameters  of the best model.}
\item{gradient.objpar}{Gradient of the object parameters for the best model.}
\item{gradient.attpar}{Gradient of the attribute parameters for the best model.}
\item{gradient.sizepar}{Gradient of the class size parameters for the best model.}
\item{fitmeasures}{A list of model selection criteria and goodness-of-fit criteria for the model with the highest posterior density.}
\item{postprob}{A \emph{I X T} matrix of posterior probabilities for the best model.}
\item{margprob.JK}{A \emph{J X K} matrix of marginal object-attribute association probabilities.}
\item{condprob.JKT}{A \emph{J X K X T} array of conditional object-attribute association probabilities (i.e., probability of object-attribute association given latent class membership).}
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

Maris, E. (1999). Estimating multiple classification latent class models. \emph{Psychometrika, 64}, 187-212.

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

Meulders, M., Tuerlinckx, F., and Vanpaemel, W. (2013). Constrained multilevel latent class models for the analysis of three-way three-mode binary data. 
\emph{Journal of Classification, 30 (3)}, 306-337.

Tanner, M. A. (1996). \emph{Tools for statistical inference: Methods for the exploration of
posterior distributions and likelihood functions} (Third ed.). New York:
Springer-Verlag.

Tatsuoka, K. (1984). \emph{Analysis of errors in fraction addition and subtraction problems}. Final Report for NIE-G-81-0002, University of Illinois, Urbana-Champaign. 

Schwarz, G. (1978). Estimating the dimensions of a model. \emph{Annals of Statistics, 6}, 461-464.


}
\author{
Michel Meulders and Philippe De Bruecker
}
%%\note{
%%  ~~further notes~~
%%}
\seealso{\code{\link{print.LCplfm}}}
\examples{

\dontrun{

# example 1: analysis on determinants of anger-related behavior

# load anger data
data(anger)

# estimate a disjunctive LCplfm model with F=2 and T=2 
# assume constant situation-feature classification
# and class-specific situation parameters (i.e. model 1)
# use 10 exploratory runs with random starting points 
anger.LCplfm<-LCplfm(data=anger$data,F=2, T=2, M=10)

# print the output of the model 
print (anger.LCplfm)


# estimate a disjunctive LCplfm model with F=4 and T=2
# assume constant situation-feature classifications
# and class-specific situation parameters (i.e. model 1)
# use 20 exploratory runs with random starting points (M=20)
# constrain parameters of subsequent behavior pairs to "load"
# on only one feature

# specify which attribute parameters have to be estimated from the data
update.attribute<-matrix(rep(0,8*4),ncol=4)
update.attribute[1:2,1]<-c(1,1)
update.attribute[3:4,2]<-c(1,1)
update.attribute[5:6,3]<-c(1,1)
update.attribute[7:8,4]<-c(1,1)

# specify starting values for attribute parameters in each of M=20 runs
# for parameters with update.attribute==0 starting values are constrained to 1e-6
# for parameters with update.attribute==1 starting values are sampled from a unif(0,1)
start.attribute<-array(runif(8*4*20),c(8,4,20))
start.attribute[update.attribute\%o\%rep(1,20)==0]<-1e-6 

# estimate the constrained model
anger.LCplfm.constr<-LCplfm(data=anger$data,F=4, T=2, M=20, 
                     update.attributeparameters=update.attribute,
                     start.attributeparameters=start.attribute)

# estimate a disjunctive LCplfm model with F=4 and T=2
# assume constant situation-feature classifications
# class-specific situation and bahavior parameters (i.e. model 3)
# use 20 exploratory runs with random starting points (M=20)
# constrain parameters of subsequent behavior pairs to "load"
# on only one feature

# specify which attribute parameters have to be estimated from the data
 
update.attribute<-matrix(rep(0,8*4),ncol=4)
update.attribute[1:2,1]<-c(1,1)
update.attribute[3:4,2]<-c(1,1)
update.attribute[5:6,3]<-c(1,1)
update.attribute[7:8,4]<-c(1,1)
update.attribute<-update.attribute\%o\%rep(1,2)

# specify starting values for attribute parameters in each of M=20 runs
# for parameters with update.attribute==0 starting values are constrained to 1e-6
# for parameters with update.attribute==1 starting values are sampled from a unif(0,1)
start.attribute<-array(runif(8*4*2*20),c(8,4,2,20))
start.attribute[update.attribute\%o\%rep(1,20)==0]<-1e-6 

# estimate the constrained model
anger.LCplfm.m3.constr<-LCplfm(data=anger$data,F=4, T=2, M=20, model=3, 
                     update.attributeparameters=update.attribute,
                     start.attributeparameters=start.attribute)

}

\dontrun{
# example 2: analysis of car perception data

# load car data
data(car)

# estimate a disjunctive LCplfm with F=3 and T=2
# assume constant attribute-feature classification
# and class-specific car parameters (i.e. model 4)
# use 10 exploratory runs with random starting points 
car.LCplfm<-LCplfm(data=car$data3w,F=3, T=2, M=10,model=4)

# print the output of the model 
print(car.LCplfm)
}

\dontrun{

# example 3: estimation of multiple classification latent class 
# model (Maris, 1999) for cognitive diagnosis


# load subtraction data
library(CDM)
data(fraction.subtraction.data)
data(fraction.subtraction.qmatrix)


# create three-way data as input for LCplfm
I<-536
J<-1
K<-20
data3w<-array(c(as.matrix(fraction.subtraction.data)),c(I,J,K))

# add item labels

itemlabel<-c("5/3 - 3/4", 
"3/4 - 3/8", 
"5/6 - 1/9",
"3 1/2 - 2 3/2", 
"4 3/5 - 3 4/10", 
"6/7 - 4/7", 
"3 - 2 1/5", 
"2/3 - 2/3", 
"3 7/8 - 2", 
"4 4/12 - 2 7/12", 
"4 1/3 - 2 4/3", 
"1 1/8 - 1/8", 
"3 3/8 - 2 5/6", 
"3 4/5 - 3 2/5", 
"2 - 1/3", 
"4 5/7 - 1 4/7", 
"7 3/5 - 4/5", 
"4 1/10 - 2 8/10", 
"4 - 1 4/3", 
"4 1/3 - 1 5/3") 

dimnames(data3w)[[3]]<-itemlabel

# estimate multiple classification latent class model (Maris, 1999)

set.seed(537982)
subtract.m1.lst<-stepLCplfm(data3w,minF=3,maxF=5,minT=1,maxT=3,model=1,M=20,maprule="conj")


# print BIC values
sumar<-summary(subtract.m1.lst)
as.matrix(sort(sumar[,5]))

# print output best model
subtract.m1.lst[[5,2]]

# correlation between extracted skills and qmatrix
round(cor(fraction.subtraction.qmatrix,subtract.m1.lst[[5,2]]$attpar),2)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
