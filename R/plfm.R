##################################################################
### function gendat
##################################################################

gendat<-function(maprule="disj",N,objpar,attpar){
## function probability disjunctive model

pdisj<-function(F,objpar,attpar){
prob<-matrix(rep(1,J*K),nrow=J)
for (f in 1:F) {prob<-prob*(1-objpar[,f]%o%attpar[,f])}
pdisj<-1-prob
}

## function probability conjunctive model

pconj<-function(F,objpar,attpar){
prob<-matrix(rep(1,J*K),nrow=J)
for (f in 1:F) {prob<-prob*(1-(1-objpar[,f])%o%attpar[,f])}
pconj<-prob
}

J<-dim(objpar)[1]
F<-dim(objpar)[2]
K<-dim(attpar)[1]


if ((dim(objpar)[2]!=dim(attpar)[2])) print("Number of columns in matrix of object parameters and matrix of attribute parameters should be the same")
else if (sum((objpar<0)|(objpar>1)|is.nan(objpar))>0) print("Object parameters should all be non-missing and between 0 and 1")
else if (sum((attpar<0)|(attpar>1)|is.nan(attpar))>0) print("Attribute parameters should all be between 0 and 1")
else if ((maprule!="disj")&(maprule!="conj")) print("Incorrect specification of mapping rule")
else if (N<0) print("The number of replications should be positive (N>0)")
else if (N==0)
{
if (maprule=="disj") prob1<-pdisj(F,objpar,attpar)
else if (maprule=="conj") prob1<-pconj(F,objpar,attpar)
gendat<-list(call=match.call(),prob1=prob1)
}
else if (N>0){
if (maprule=="disj") prob1<-pdisj(F,objpar,attpar)
else if (maprule=="conj") prob1<-pconj(F,objpar,attpar)
freqtot<-matrix(rep(N,J*K),nrow=J)
freq1<-matrix(rbinom(J*K,N,prob1),nrow=J)
gendat<-list(call=match.call(),prob1=prob1,freq1=freq1,freqtot=freqtot)
}
}


##################################################################
### function PLFM
##################################################################


plfm<-function(data,object,attribute,rating,freq1,freqtot,F,datatype="freq",maprule="disj",M=5,emcrit1=1e-2,emcrit2=1e-10,printrun=TRUE){



#### functions


pi<-function(J,K,F,sigma,rho){
prob<-matrix(rep(1,J*K),nrow=J)
for (f in 1:F) {prob<-prob*(1-sigma[,f]%o%rho[,f])}
pi<-1-prob
}

logpost<-function(p1,sigma,rho,freq1,freq0){
logpost<-sum(log(sigma)+log(1-sigma))+sum(log(rho)+log(1-rho))+sum(freq1*log(p1)+freq0*log(1-p1))
}

loglik<-function(p1,freq1,freq0){
loglik<-sum(freq1*log(p1)+freq0*log(1-p1))
}

computeHessian<-function(sigma.o,rho.o,psi,nsigma.o,nrho.o,nfreq1,nfreq0,np1.o,freq1mat,pi1mat,J,K,F,NPsigma,NPtot){
  ## Hsigma
  
  Hsigma<-1/sigma.o**2+1/((1-sigma.o)**2)+apply((nrho.o/(1-psi))**2*(nfreq1*((1-np1.o)/np1.o)**2+nfreq0),c(1,3),sum)  

  ## Hrho
  
  Hrho<-1/rho.o**2+1/((1-rho.o)**2)+apply((nsigma.o/(1-psi))**2*(nfreq1*((1-np1.o)/np1.o)**2+nfreq0),c(2,3),sum)

  ## Hsigma.common
  
  temp1<-c(t(rho.o))%o%rep(1,F)
  temp2<-matrix(c(aperm((rep(1,F)%o%rho.o),c(3,1,2))),nrow=F*K,byrow=TRUE)
  temp3<-aperm(array(c(t(temp1*temp2)),c(F,F,K)),c(3,1,2))
  temp4<-rep(1,J)%o%temp3
  temp5<-c(aperm(psi,c(3,1,2)))%o%rep(1,F)
  temp6<-matrix(c(aperm(psi%o%rep(1,F),c(3,4,1,2))),nrow=J*K*F,byrow=TRUE)
  temp7<-(1-temp5)*(1-temp6)
  temp8<-aperm(array(t(temp7),c(F,F,J,K)),c(3,4,1,2))
  Hsigma.common<-apply(freq1mat*temp4*(1-pi1mat)/((pi1mat**2)*temp8),c(1,3,4),sum)


  ## Hrho.common
  
  temp1<-c(t(sigma.o))%o%rep(1,F)
  temp2<-matrix(c(aperm((rep(1,F)%o%sigma.o),c(3,1,2))),nrow=F*J,byrow=TRUE)
  temp3<-aperm(array(c(t(temp1*temp2)),c(F,F,J)),c(3,1,2))
  temp4<-rep(1,K)%o%temp3
  temp4<-aperm(temp4,c(2,1,3,4))
  Hrho.common<-apply(freq1mat*temp4*(1-pi1mat)/((pi1mat**2)*temp8),c(2,3,4),sum)


  ## Hsigmarho.common
  
  Hsigmarho.common<-nfreq1*((1-np1.o)/(1-psi))*(-(1/np1.o)+(psi*(1-np1.o)/((1-psi)*np1.o**2)))+nfreq0/((1-psi)**2)

  ## Hsigmarho.distinct

  temp9<-(rep(1,J)%o%rho.o)%o%rep(1,F)
  temp10<-rep(1,K)%o%(aperm(array(t(c(t(sigma.o))%o%rep(1,F)),c(F,F,J)),c(3,1,2)))
  temp11<-aperm(temp10,c(2,1,3,4))
  temp12<-aperm(temp9*temp11,c(1,2,4,3))
  temp13<-aperm(sigma.o%o%rep(1,K),c(1,3,2))
  temp14<-rep(1,J)%o%rho.o
  temp15<-(temp13*temp14)%o%(rep(1,F))
  temp16<-aperm(temp15,c(1,2,4,3))
  temp17<-(1-temp15)*(1-temp16)
  Hsigmarho.distinct<-aperm((freq1mat*temp12*(1-pi1mat))/((pi1mat**2)*temp17),c(1,2,4,3))

  #########################
  ##complete Hessian matrix
  #########################

  ## deriv(sigma,sigma)

  Hes.sigma<-array(rep(0,J*J*F*F),c(J,F,J,F))

  for (j1 in 1:J){
    for (f1 in 1:F){
      Hes.sigma[j1,f1,j1,f1]<-Hsigma[j1,f1]}}

  for (j in 1:J){
    for (f1 in 1:F){
      for (f2 in 1:F){
    
      if (f1!=f2) Hes.sigma[j,f1,j,f2]<-Hsigma.common[j,f1,f2]}}}

  Hesmat.sigma<-matrix(aperm(Hes.sigma,c(2,1,4,3)),nrow=J*F)

  ## deriv(rho,rho)

  Hes.rho<-array(rep(0,K*K*F*F),c(K,F,K,F))

  for (k in 1:K){
    for (f1 in 1:F){
     Hes.rho[k,f1,k,f1]<-Hrho[k,f1]}}

  for (k in 1:K){
    for (f1 in 1:F){
      for (f2 in 1:F){
        if (f1!=f2) Hes.rho[k,f1,k,f2]<-Hrho.common[k,f1,f2]}}}

  Hesmat.rho<-matrix(aperm(Hes.rho,c(2,1,4,3)),nrow=K*F)

  ##deriv(sigma,rho)

  Hes.sigmarho<-array(rep(0,J*F*K*F),c(J,F,K,F))

  for (j in 1:J){
   for (k in 1:K){
    for (f1 in 1:F){
     Hes.sigmarho[j,f1,k,f1]<-Hsigmarho.common[j,k,f1]}}}

  for (j in 1:J){
    for (k in 1:K){
      for (f1 in 1:F){
        for (f2 in 1:F){
         if (f1!=f2) Hes.sigmarho[j,f1,k,f2]<-Hsigmarho.distinct[j,k,f1,f2]}}}}

   Hesmat.sigmarho<-matrix(c(aperm(Hes.sigmarho,c(2,1,4,3))),nrow=J*F)

   ##deriv(rho,sigma) 

   Hessian<-matrix(rep(0,(NPtot)**2),nrow=NPtot)
   Hessian[1:NPsigma,1:NPsigma]<-Hesmat.sigma
   Hessian[1:NPsigma,(NPsigma+1):NPtot]<-Hesmat.sigmarho  
   Hessian[((NPsigma+1):NPtot),1:NPsigma]<-t(Hesmat.sigmarho)
   Hessian[(NPsigma+1):NPtot,(NPsigma+1):NPtot]<-Hesmat.rho
   computeHessian<-Hessian

}



## parameter estimation

if (datatype=="freq") {
 if (length(freqtot)==1) freqtot<-matrix(rep(freqtot,((dim(freq1)[1])*(dim(freq1)[2]))),nrow=dim(freq1)[1])
}

if (datatype=="dataframe"){
agg<-table(data$object,data$attribute,data$rating)
freq1<-agg[,,"1"]
freqtot<-agg[,,"0"]+agg[,,"1"]}


if (sum(is.nan(freq1)|is.na(freq1)|freq1<0)){print("The elements of the matrix freq1 should be non-missing and larger than or equal to 0")}
else if (sum(is.nan(freqtot)|is.na(freqtot)|freqtot<1)){print("The elements of the matrix freqtot should be non-missing and larger than 0")}
else if (sum(abs(dim(freqtot)-dim(freq1)))>0){print("The matrices freq1 and freqtot should have the same dimensionality")}
else if (sum(freq1>freqtot)>0) {print("The corresponding values of the matrix freqtot should be larger than or equal to the elements in the matrix freq1")}
else if ((maprule!="disj") & (maprule!="conj")) {print("The mapping rule is not correctly specified")}
else if (M<1) {print("The number of requested runs should be at least one (M>=1)")}
else if (((dim(freq1)[1])*(dim(freq1)[2]))<F*sum(dim(freq1))) {print("The number of model parameters (J+K)*F should be smaller than the number of observations J*K used to fit the model")}
else if (emcrit1<0 | emcrit2<0) {print("emcrit1 and emcrit2 should be small positive numbers")}


else{

## initialisation

N<-max(freqtot)
J<-dim(freq1)[1]
K<-dim(freq1)[2]
freq0<-freqtot-freq1

if (is.null(rownames(freq1))=="FALSE") rowlabels<-rownames(freq1)
else rowlabels<-paste("R_", seq(1,J),sep="")
if (is.null(colnames(freq1))=="FALSE") columnlabels<-colnames(freq1)
else columnlabels<-paste("C_", seq(1,K),sep="")

if (maprule=="conj"){
freq0<-freq1
freq1<-freqtot-freq0
}


NPsigma<-J*F
NPrho<-K*F
NPtot<-J*F+K*F

featurelabels<-paste("F",seq(1,F),sep="")
runlabels<-paste("RUN",seq(1,M),sep="")

nfreq1<-freq1%o%rep(1,F)
nfreq0<-freq0%o%rep(1,F)
nfreqtot<-freqtot%o%rep(1,F)
freq1mat<-(freq1%o%rep(1,F))%o%rep(1,F)

margobj<-apply(freqtot,1,sum)%o%rep(1,F)
margatt<-apply(freqtot,2,sum)%o%rep(1,F)

psi<-array(rep(0,J*K*F),c(J,K,F))
ksi<-array(rep(0,J*K*F),c(J,K,F))
eps<-array(rep(0,J*K*F),c(J,K,F))
gamma<-array(rep(0,J*K*F),c(J,K,F))
alfa<-array(rep(0,J*K*F),c(J,K,F))

sigma.runs<-array(rep(0,J*F*M),c(J,F,M))
rho.runs<-array(rep(0,K*F*M),c(K,F,M))
logpost.runs<-rep(0,M)

dimnames(sigma.runs)[[1]]<-rowlabels
dimnames(sigma.runs)[[2]]<-featurelabels
dimnames(sigma.runs)[[3]]<-runlabels
dimnames(rho.runs)[[1]]<-columnlabels
dimnames(rho.runs)[[2]]<-featurelabels
dimnames(rho.runs)[[3]]<-runlabels
names(logpost.runs)<-runlabels

for (run in 1:M){

## set convergence criterion to switch from EM to EM+NR
emcritX<-emcrit1

## initialisation starting values

sigma.o<-matrix(runif(J*F),nrow=J,byrow=TRUE)
sigma.n<-matrix(runif(J*F),nrow=J,byrow=TRUE)
rho.o<-matrix(runif(K*F),nrow=K,byrow=TRUE)
rho.n<-matrix(runif(K*F),nrow=K,byrow=TRUE)
theta.o<-c(t(sigma.o),t(rho.o))
theta.n<-c(t(sigma.n),t(rho.n))
p1.o<-pi(J,K,F,sigma.o,rho.o)
p1.n<-pi(J,K,F,sigma.n,rho.n)
logpost.o<-logpost(p1.o,sigma.o,rho.o,freq1,freq0)
logpost.n<-logpost(p1.n,sigma.n,rho.n,freq1,freq0)


while ((abs(logpost.o-logpost.n))>emcrit2){
  
  flagNR<-0
  sigma.o<-sigma.n
  rho.o<-rho.n
  nrho.o<-rep(1,J)%o%rho.o
  nsigma.o<-aperm(sigma.o%o%rep(1,K),c(1,3,2))
  theta.o<-c(t(sigma.o),t(rho.o))

  p1.o<-pi(J,K,F,sigma.o,rho.o)
  np1.o<-p1.o%o%rep(1,F)
  pi1mat<-np1.o%o%rep(1,F)

  logpost.o<-logpost(p1.o,sigma.o,rho.o,freq1,freq0)

  for (f in 1:F){psi[,,f]<-(sigma.o[,f]%o%rho.o[,f])}
  for (f in 1:F){ksi[,,f]<-(sigma.o[,f]%o%(1-rho.o[,f]))}
  for (f in 1:F){eps[,,f]<-((1-sigma.o[,f])%o%rho.o[,f])}

  gamma<-ksi/(1-psi)
  alfa<-eps/(1-psi)
  sigma.n<-(1+apply(nfreq1*(1-gamma)*psi*(1/np1.o)+nfreqtot*gamma,c(1,3),sum))/(2+margobj)
  rho.n<-(1+apply(nfreq1*(1-alfa)*psi*(1/np1.o)+nfreqtot*alfa,c(2,3),sum))/(2+margatt)
  p1.n<-pi(J,K,F,sigma.n,rho.n)
  logpost.n<-logpost(p1.n,sigma.n,rho.n,freq1,freq0)

 while ((flagNR==0) & (abs(logpost.o-logpost.n)<emcritX) & (abs(logpost.o-logpost.n)>emcrit2)){

  
  ##compute EM step
  


  sigma.o<-sigma.n
  rho.o<-rho.n
  nrho.o<-rep(1,J)%o%rho.o
  nsigma.o<-aperm(sigma.o%o%rep(1,K),c(1,3,2))
  theta.o<-c(t(sigma.o),t(rho.o))

  p1.o<-pi(J,K,F,sigma.o,rho.o)
  np1.o<-p1.o%o%rep(1,F)
  pi1mat<-np1.o%o%rep(1,F)

  logpost.o<-logpost(p1.o,sigma.o,rho.o,freq1,freq0)

  for (f in 1:F){psi[,,f]<-(sigma.o[,f]%o%rho.o[,f])}
  for (f in 1:F){ksi[,,f]<-(sigma.o[,f]%o%(1-rho.o[,f]))}
  for (f in 1:F){eps[,,f]<-((1-sigma.o[,f])%o%rho.o[,f])}

  gamma<-ksi/(1-psi)
  alfa<-eps/(1-psi)

  sigma.EM<-(1+apply(nfreq1*(1-gamma)*psi*(1/np1.o)+nfreqtot*gamma,c(1,3),sum))/(2+margobj)
  rho.EM<-(1+apply(nfreq1*(1-alfa)*psi*(1/np1.o)+nfreqtot*alfa,c(2,3),sum))/(2+margatt)
  theta.EM<-c(t(sigma.EM),t(rho.EM))
  
  
  
  
  ##compute minus second derivative augmented posterior
  

  ES1<-apply(nfreq1*(1/np1.o)*(psi+ksi*(np1.o-psi)/(1-psi))+nfreq0*ksi/(1-psi),c(1,3),sum)
  ET1<-apply(nfreq1*(1/np1.o)*(psi+eps*(np1.o-psi)/(1-psi))+nfreq0*eps/(1-psi),c(2,3),sum)
  secderivsigma<-(1+ES1)/(sigma.o**2)+(1+margobj-ES1)/((1-sigma.o)**2) 
  secderivrho<-(1+ET1)/(rho.o**2)+(1+margatt-ET1)/((1-rho.o)**2)
  secderiv<-diag(J*F+K*F)*c(t(secderivsigma),t(secderivrho))

  

  ##compute Hessian


  completeHessian<-computeHessian(sigma.o,rho.o,psi,nsigma.o,nrho.o,nfreq1,nfreq0,np1.o,freq1mat,pi1mat,J,K,F,NPsigma,NPtot)
  

  ##compute acceleration step

  
  theta.n<-theta.o+solve(completeHessian)%*%(secderiv%*%(theta.EM-theta.o))
  sigma.n<-matrix(theta.n[1:NPsigma],nrow=J,byrow=TRUE)
  rho.n<-matrix(theta.n[(NPsigma+1):NPtot],nrow=K,byrow=TRUE)
  
  ## check convergence Newton-Rhapson step
  if (sum(ifelse(((sigma.n>0) & (sigma.n<1)),0,1))+sum(ifelse(((rho.n>0) & (rho.n<1)),0,1))>0) 
  {
   sigma.n<-sigma.EM
   rho.n<-rho.EM
   
   
   emcritX<-emcritX/10
   flagNR<-1
   ##print(c("change EMcrit",emcritX,emcrit1))
  }
  
  
  p1.n<-pi(J,K,F,sigma.n,rho.n)
  logpost.n<-logpost(p1.n,sigma.n,rho.n,freq1,freq0)
  
  
}## end loop EM+NR

}##end loop EM



 sigma.runs[,,run]<-sigma.n
 rho.runs[,,run]<-rho.n
 logpost.runs[run]<-logpost.n

 if (printrun=="TRUE"){
 	if (maprule=="disj") {print(paste("DISJUNCTIVE ANALYSIS  ","F=",F,"  RUN=",run,sep=""),quote="FALSE")} 
 	else if (maprule=="conj") {print(paste("CONJUNCTIVE ANALYSIS  ","F=",F,"  RUN=",run,sep=""),quote="FALSE")} 
 }

}## end M



### select best solution


## select best solution and compute fit measures
best<-order(logpost.runs,decreasing=TRUE)[1]
sigma.best<-as.matrix(sigma.runs[,,best])
rownames(sigma.best)<-rowlabels
colnames(sigma.best)<-featurelabels
rho.best<-as.matrix(rho.runs[,,best])
rownames(rho.best)<-columnlabels
colnames(rho.best)<-featurelabels
p1.best<-pi(J,K,F,sigma.best,rho.best)
rownames(p1.best)<-rowlabels
colnames(p1.best)<-columnlabels
np1.best<-p1.best%o%rep(1,F)
pi1mat.best<-np1.best%o%rep(1,F)
for (f in 1:F){psi[,,f]<-(sigma.best[,f]%o%rho.best[,f])}
nrho.best<-rep(1,J)%o%rho.best
nsigma.best<-aperm(sigma.best%o%rep(1,K),c(1,3,2))



## computation gradient final solution


## gradient observed posterior
gradsigma<-1/sigma.best-1/(1-sigma.best)+apply(((nrho.best)/(1-psi))*(nfreq1*(1-np1.best)/np1.best-nfreq0),c(1,3),sum)
gradrho<-1/rho.best-1/(1-rho.best)+apply(((nsigma.best)/(1-psi))*(nfreq1*(1-np1.best)/np1.best-nfreq0),c(2,3),sum)


##compute Hessian best solution


completeHessian<-computeHessian(sigma.best,rho.best,psi,nsigma.best,nrho.best,nfreq1,nfreq0,np1.best,freq1mat,pi1mat.best,J,K,F,NPsigma,NPtot)
SEsigma<-1/sqrt(matrix(diag(completeHessian)[1:NPsigma],nrow=J,byrow=TRUE))
rownames(SEsigma)<-rowlabels
colnames(SEsigma)<-featurelabels

SErho<-1/sqrt(matrix(diag(completeHessian)[(NPsigma+1):NPtot],nrow=K,byrow=TRUE))
rownames(SErho)<-columnlabels
colnames(SErho)<-featurelabels



## compute fit measures



loglik.best<-loglik(p1.best,freq1,freq0)
logpost.best<-logpost(p1.best,sigma.best,rho.best,freq1,freq0)
correlation<-cor(c(p1.best*freqtot),c(freq1))
VAF<-correlation**2
deviance<--2*loglik.best
AIC<-deviance+2*NPtot
BIC<-deviance+NPtot*log(N)
chisquare.fit<-sum(((freq1-freqtot*p1.best)**2/(freqtot*p1.best*(1-p1.best))))
df<-(J*K)-(J+K)*F
pchi<-1-pchisq(chisquare.fit,df)
fitmeasures<-matrix(c(loglik.best,logpost.best,deviance,AIC,BIC,chisquare.fit,df,pchi,correlation,VAF),ncol=1)
rownames(fitmeasures)<-c("Log Likelihood","Log Posterior","Deviance","AIC","BIC","Pearson Chi-square","df","p-value","Correlation observed and expected frequencies","VAF observed frequencies")

if (maprule=="disj")
plfm<-list(call=match.call(),objpar=sigma.best,attpar=rho.best,fitmeasures=fitmeasures,logpost.runs=logpost.runs,objpar.runs=sigma.runs,attpar.runs=rho.runs,bestsolution=best,gradient.objpar=gradsigma,gradient.attpar=gradrho,SE.objpar=SEsigma,SE.attpar=SErho,prob1=p1.best)

else if (maprule=="conj")

plfm<-list(call=match.call(),objpar=1-sigma.best,attpar=rho.best,fitmeasures=fitmeasures,logpost.runs=logpost.runs,objpar.runs=1-sigma.runs,attpar.runs=rho.runs,bestsolution=best,gradient.objpar=gradsigma,gradient.attpar=gradrho,SE.objpar=SEsigma,SE.attpar=SErho,prob1=1-p1.best)


class(plfm)<-"plfm"
plfm
}
}



############################
## function print.plfm
############################

print.plfm<-function(x, ...)
{
cat("Call:\n")
print(x$call)

cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
descript<-as.matrix(x$fitmeasures[9:10,])
colnames(descript)<-""
print(descript,digits=3)

cat("\nESTIMATE OBJECT PARAMETERS:\n")
strout<-as.matrix(capture.output(round(x$objpar,2)))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\nESTIMATE ATTRIBUTE PARAMETERS:\n")
strout<-as.matrix(capture.output(round(x$attpar,2)))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")
}


###############################
## function summary.plfm
###############################

summary.plfm <- function(object, ...)
{

crit<-as.matrix(object$fitmeasures[1:5,])
colnames(crit)<-""
crit<-round(crit)

chi<-as.matrix(object$fitmeasures[6:8,])
colnames(chi)<-""
chi<-round(chi,3)

descript<-as.matrix(object$fitmeasures[9:10,])
colnames(descript)<-""
descript<-round(descript,3)


featurelabels<-colnames(object$objpar)
F<-length(featurelabels)
estobjectpar<-round(object$objpar,2)
colnames(estobjectpar)<-paste("Est(",featurelabels,")",sep="")
seobjectpar<-round(object$SE.objpar,3)
colnames(seobjectpar)<-paste("SE(",featurelabels,")",sep="")

estattpar<-round(object$attpar,2)
colnames(estattpar)<-paste("Est(",featurelabels,")",sep="")
seattpar<-round(object$SE.attpar,3)
colnames(seattpar)<-paste("SE(",featurelabels,")",sep="")

result<-list(call=object$call,
             informationcriteria=crit,
             chisquaretest=chi,
             descriptivefit=descript,
             objpar=estobjectpar,
             SE.objpar=seobjectpar,
             attpar=estattpar,
             SE.attpar=seattpar)
class(result)<-"summary.plfm"
result
}

##############################################
## function print.summary.plfm
##############################################

print.summary.plfm<-function(x, ...)
{
cat("Call:\n")
print(x$call)

cat("\nINFORMATION CRITERIA:\n")
print(x$informationcriteria)

cat("\nPEARSON CHI SQUARE TEST OBJECT X ATTRIBUTE TABLE:\n")
print(x$chisquaretest)

cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
print(x$descriptivefit)

cat("\nESTIMATE OBJECT PARAMETERS:\n")
strout<-as.matrix(capture.output(x$objpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\nSTANDARD ERROR OBJECT PARAMETERS:\n")
strout<-as.matrix(capture.output(x$SE.objpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$SE.objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\nESTIMATE ATTRIBUTE PARAMETERS:\n")
strout<-as.matrix(capture.output(x$attpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\nSTANDARD ERROR ATTRIBUTE PARAMETERS:\n")
strout<-as.matrix(capture.output(x$SE.attpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$SE.attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")
}


##############################
## function stepplfm
##############################

stepplfm<-function(minF,maxF,data,object,attribute,rating,freq1,freqtot,datatype="freq",maprule="disj",M=5,emcrit1=1e-2,emcrit2=1e-10,printrun=TRUE){

tempcall<-match.call()
tempcall<-tempcall[c(-2,-3)]


if (maprule=="disj"){
      stepplfm<-vector(mode="list",maxF-minF+1)        
	for (f in (1:(maxF-minF+1))){
       stepplfm[[f]]<-plfm(maprule="disj",data=data,object=object,attribute=attribute,rating=rating,freq1=freq1,freqtot=freqtot,datatype=datatype,F=f,M=M,emcrit1=emcrit1,emcrit2=emcrit2,printrun=printrun)
       tempcall$F<-as.numeric(f)
       stepplfm[[f]]$call<-tempcall
      }

 }     

else if (maprule=="conj"){
      stepplfm<-vector(mode="list",maxF-minF+1)          
	for (f in (1:(maxF-minF+1))){
       stepplfm[[f]]<-plfm(maprule="conj",data=data,object=object,attribute=attribute,rating=rating,freq1=freq1,freqtot=freqtot,datatype=datatype,F=f,M=M,emcrit1=emcrit1,emcrit2=emcrit2,printrun=printrun)
       tempcall$F<-as.numeric(f)
       stepplfm[[f]]$call<-tempcall

      }

 }  
else if (maprule=="disj/conj"){
	
      stepdisj<-vector(mode="list",maxF-minF+1)
      stepconj<-vector(mode="list",maxF-minF+1)          
	for (f in (1:(maxF-minF+1))){
       stepdisj[[f]]<-plfm(maprule="disj",data=data,object=object,attribute=attribute,rating=rating,freq1=freq1,freqtot=freqtot,datatype=datatype,F=f,M=M,emcrit1=emcrit1,emcrit2=emcrit2,printrun=printrun)
       tempcall$F<-as.numeric(f)
      stepdisj[[f]]$call<-tempcall
      }
      for (f in (1:(maxF-minF+1))){
       stepconj[[f]]<-plfm(maprule="conj",data=data,object=object,attribute=attribute,rating=rating,freq1=freq1,freqtot=freqtot,datatype=datatype,F=f,M=M,emcrit1=emcrit1,emcrit2=emcrit2,printrun=printrun)
       tempcall$F<-as.numeric(f)
      stepconj[[f]]$call<-tempcall

      }

	stepplfm<-list(disj=stepdisj,conj=stepconj)
}
class(stepplfm)<-"stepplfm"
stepplfm
}



############################
## function print.stepplfm
############################

print.stepplfm<-function(x, ...)
{
if ((length(x$disj)==0)&(length(x$conj)==0)){
	nfeat<-length(x)
	minfeat<-dim(x[[1]]$objpar)[2]
	maxfeat<-dim(x[[nfeat]]$objpar)[2]
	fitm<-t(sapply(x,function(obj) obj$fitmeasures))
	colnames(fitm)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	rownames(fitm)<-paste("F=",seq(minfeat:maxfeat),sep="")

      if (x[[1]]$call$maprule=="disj"){
       	cat("\n*****************************") 
      	cat("\nDISJUNCTIVE MODEL\n")
      	cat("******************************\n") 
    		}
     else if (x[[1]]$call$maprule=="conj"){
		cat("\n*********************************") 
     		cat("\nCONJUNCTIVE MODEL\n")
     		cat("*********************************\n") 
		}
	cat("\nINFORMATION CRITERIA:\n")
      cat("\n")
	print(round(fitm[,c(1:5)],0))


	cat("\nPEARSON CHI SQUARE TEST OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitm[,c(6:8)],3))

	cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitm[,c(9,10)],3))
	}
else if ((length(x$disj)>0)&(length(x$conj)>0)){
	nfeat<-length(x$disj)
      minfeat<-dim(x$disj[[1]]$objpar)[2]
      maxfeat<-dim(x$disj[[nfeat]]$objpar)[2]
      fitdisj<-t(sapply(x$disj,function(obj) obj$fitmeasures))
      fitconj<-t(sapply(x$conj,function(obj) obj$fitmeasures))
      colnames(fitdisj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	rownames(fitdisj)<-paste("F=",seq(minfeat:maxfeat),sep="")
      colnames(fitconj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	rownames(fitconj)<-paste("F=",seq(minfeat:maxfeat),sep="")
      
      cat("\n*****************************") 
      cat("\nDISJUNCTIVE MODEL\n")
      cat("******************************\n") 
      cat("\nINFORMATION CRITERIA:\n")
      cat("\n")
	print(round(fitdisj[,c(1:5)],0))


	cat("\nPEARSON CHI SQUARE TEST OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitdisj[,c(6:8)],3))

	cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitdisj[,c(9,10)],3))

      cat("\n*********************************") 
      cat("\nCONJUNCTIVE MODEL\n")
      cat("*********************************\n") 
      cat("\nINFORMATION CRITERIA:\n")
      cat("\n")
	print(round(fitconj[,c(1:5)],0))


	cat("\nPEARSON CHI SQUARE TEST OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitconj[,c(6:8)],3))

	cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
      cat("\n")
	print(round(fitconj[,c(9,10)],3))


}
}




################################
## function plot.stepplfm
################################

plot.stepplfm<-function(x,which="BIC",...){

if ((length(x$disj)==0)&(length(x$conj)==0)){
	nfeat<-length(x)
      if (nfeat==1){print("The fit of the model is only displayed if F>1",quote="FALSE")}
      else{ 
		minfeat<-dim(x[[1]]$objpar)[2]
		maxfeat<-dim(x[[nfeat]]$objpar)[2]
		fitm<-t(sapply(x,function(obj) obj$fitmeasures))
		colnames(fitm)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
		rownames(fitm)<-paste("F=",seq(minfeat:maxfeat),sep="")
            plot(seq(minfeat:maxfeat),fitm[,which],xlab="Number of features",ylab=which,type="b",xaxp=c(minfeat,maxfeat,maxfeat-minfeat),...)


      	
	}
}
else if ((length(x$disj)>0)&(length(x$conj)>0)){
	nfeat<-length(x$disj)
      if (nfeat==1){print("The fit of the model is only displayed if F>1",quote="FALSE")}
      else{ 
      	minfeat<-dim(x$disj[[1]]$objpar)[2]
      	maxfeat<-dim(x$disj[[nfeat]]$objpar)[2]
      	fitdisj<-t(sapply(x$disj,function(obj) obj$fitmeasures))
      	fitconj<-t(sapply(x$conj,function(obj) obj$fitmeasures))
		colnames(fitdisj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
		rownames(fitdisj)<-paste("F=",seq(minfeat:maxfeat),sep="")
		colnames(fitconj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
		rownames(fitconj)<-paste("F=",seq(minfeat:maxfeat),sep="")
            comb<-rbind(fitdisj,fitconj)
            minstat<-apply(comb,2,min)
            maxstat<-apply(comb,2,max)

		plot(seq(minfeat:maxfeat),fitdisj[,which],xlab="Number of features",ylab=which,type="b",xaxp=c(minfeat,maxfeat,maxfeat-minfeat),ylim=c(minstat[which],maxstat[which]),...)
		lines(seq(minfeat:maxfeat),fitconj[,which],lty=2,type="b")
		if ((which=="BIC")|(which=="AIC")|(which=="Deviance")|(which=="Chisquare")){
            legend("topright",c("Disjunctive","Conjunctive"),lty=c(1,2),border=NULL,bty="n")}
            else if ((which=="Correlation")|(which=="VAF")){
     	      legend("bottomright",c("Disjunctive","Conjunctive"),lty=c(1,2),border=NULL,bty="n")}
 

      }

}
}



############################
## function summary.stepplfm
############################

summary.stepplfm<-function(object, ...)
{
 if ((length(object$disj)==0)&(length(object$conj)==0)){
	nfeat<-length(object)
	minfeat<-dim(object[[1]]$objpar)[2]
	maxfeat<-dim(object[[nfeat]]$objpar)[2]
	fitm<-t(sapply(object,function(x) x$fitmeasures))
	colnames(fitm)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	

      if (object[[1]]$call$maprule=="disj"){
       rownames(fitm)<-paste("DISJUNCTIVE F=",seq(minfeat:maxfeat),sep="")
    		}
     else if (object[[1]]$call$maprule=="conj"){
       rownames(fitm)<-paste("CONJUNCTIVE F=",seq(minfeat:maxfeat),sep="")

		}

      result<-fitm
      class(result)<-"summary.stepplfm"
      result 
	}

else if ((length(object$disj)>0)|(length(object$conj)>0)){
      nfeat<-length(object$disj)
      minfeat<-dim(object$disj[[1]]$objpar)[2]
      maxfeat<-dim(object$disj[[nfeat]]$objpar)[2]
      fitdisj<-t(sapply(object$disj,function(x) x$fitmeasures))
      fitconj<-t(sapply(object$conj,function(x) x$fitmeasures))
      colnames(fitdisj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	rownames(fitdisj)<-paste("DISJUNCTIVE F=",seq(minfeat:maxfeat),sep="")
      colnames(fitconj)<-c("LogLik","LogPost","Deviance","AIC","BIC","Chisquare","df","p-value","Correlation","VAF")
	rownames(fitconj)<-paste("CONJUNCTIVE F=",seq(minfeat:maxfeat),sep="")
      result<-rbind(fitdisj,fitconj)      
      class(result)<-"summary.stepplfm"
      result 

}
}




#####################################
## function print.summary.stepplfm
#####################################

print.summary.stepplfm<-function(x,digits=2,...)
{
print(x[],digits=digits,...)
}



##################################################################
### function bayesplfm
##################################################################

bayesplfm<-function(data,object,attribute,rating,freq1,freqtot,F,Nchains=2,Nburnin=0,maxNiter=4000,Nstep=1000,Rhatcrit=1.2,maprule="disj",datatype="freq",start.bayes="best",fitted.plfm=NULL){



## functions


pi<-function(J,K,F,sigma,rho){
prob<-matrix(rep(1,J*K),nrow=J)
for (f in 1:F) {prob<-prob*(1-sigma[,f]%o%rho[,f])}
pi<-1-prob
}

post95<-function(x){post95<-quantile(x,c(0.025,0.5,0.975))}

logit<-function(x){logit<-log(x/(1-x))}

Rhat<-function(x){
nrow<-dim(x)[1]
ncol<-dim(x)[2]
gsample<-apply(x,2,mean)
varsample<-apply(x,2,var)
between<-nrow*var(gsample)
within<-mean(varsample)
varW<-var(varsample)/ncol
varB<-2*(between^2)/(ncol-1)
kwadvec<-gsample^2
cov1<-cov(varsample,kwadvec)
cov2<-cov(varsample,gsample)
muhat<-mean(gsample)
covWB<-(cov1-2*muhat*cov2)*(nrow/ncol)
postvar<-((nrow-1)/nrow)*within+(1/nrow)*(1+(1/ncol))*between
varpostvar<-((nrow-1)^2*varW +(1+(1/ncol))^2*varB +2*(nrow-1)*(1+(1/ncol))*covWB)/(nrow^2)
df<-2*(postvar^2)/varpostvar
Rhat<-sqrt((postvar/within)*((df+3)/(df+1)))
}

computedraw<-function(J,K,F,C,binmat,sigma,rho){


## compute P(X_jki=x_jki|d_ijk,theta)


probd.x0<-matrix(rep(1,C*K),nrow=C)
for (f in 1:F){probd.x0<-probd.x0*(1-binmat[,f]%o%rho[,f])}
nprobd.x0<-rep(1,J)%o%t(probd.x0)
nprobd.x1<-1-nprobd.x0

margprob.sigma<-exp(binmat%*%t(log(sigma))+(1-binmat)%*%log(1-t(sigma)))
nmargprob.sigma<-aperm(margprob.sigma%o%rep(1,K),c(2,3,1))

probx.d1<-nprobd.x1*nmargprob.sigma
probx.d0<-nprobd.x0*nmargprob.sigma

## normalize
norm1<-(apply(probx.d1,c(1,2),sum))%o%rep(1,C)
norm0<-(apply(probx.d0,c(1,2),sum))%o%rep(1,C)
probx.d1<-probx.d1/norm1
probx.d0<-probx.d0/norm0

## draw latent data objects
latx1<-array(rep(0,J*K*C),c(J,K,C))
latx0<-array(rep(0,J*K*C),c(J,K,C))
for (j in 1:J){
 for (k in 1:K){
  if (freq1[j,k]>0) latx1[j,k,]<-rmultinom(1,freq1[j,k],probx.d1[j,k,]) 
  if (freq0[j,k]>0) latx0[j,k,]<-rmultinom(1,freq0[j,k],probx.d0[j,k,]) 
}}

statobj1<-apply(latx1+latx0,c(1,3),sum)%*%binmat
statobj0<-apply(latx1+latx0,c(1,3),sum)%*%(1-binmat)

## draw objectparameters
sigma.n<-matrix(rbeta(J*F,1+statobj1,1+statobj0),nrow=J)


## compute P(Y_kji=y_kji|X_jki=x_jki,d_ijk,theta)


probd.xy0<-matrix(rep(1,C*C),nrow=C)
for (f in 1:F){probd.xy0<-probd.xy0*(1-binmat[,f]%o%binmat[,f])}
nprobd.xy0<-rep(1,K)%o%probd.xy0 
nprobd.xy1<-1-nprobd.xy0

margprob.rho<-t(exp(binmat%*%t(log(rho))+(1-binmat)%*%log(1-t(rho))))
nmargprob.rho<-aperm(margprob.rho%o%rep(1,C),c(1,3,2))

proby.xd1<-nprobd.xy1*nmargprob.rho
proby.xd0<-nprobd.xy0*nmargprob.rho

## normalize

norm1<-(apply(proby.xd1,c(1,2),sum))%o%rep(1,C)
norm0<-(apply(proby.xd0,c(1,2),sum))%o%rep(1,C)

proby.xd1<-proby.xd1/norm1
proby.xd0<-proby.xd0/norm0

proby.xd1[is.na(proby.xd1)]<-0
proby.xd0[is.na(proby.xd0)]<-0

Slatx1<-apply(latx1,c(2,3),sum)
Slatx0<-apply(latx0,c(2,3),sum)

## draw latent data attributes

laty1<-array(rep(0,K*C*C),c(K,C,C))
laty0<-array(rep(0,K*C*C),c(K,C,C))
for (k in 1:K){
 for (c in 1:C){
  if (Slatx1[k,c]>0) laty1[k,c,]<-rmultinom(1,Slatx1[k,c],proby.xd1[k,c,])
  if (Slatx0[k,c]>0) laty0[k,c,]<-rmultinom(1,Slatx0[k,c],proby.xd0[k,c,])
}}

statatt1<-(apply(laty1+laty0,c(1,3),sum))%*%binmat
statatt0<-(apply(laty1+laty0,c(1,3),sum))%*%(1-binmat)

## draw attribute parameters
rho.n<-matrix(rbeta(K*F,1+statatt1,1+statatt0),nrow=K)

computedraw<-list(sigma.n=sigma.n,rho.n=rho.n)
}

## derive freq1 and freqtot  data from input

if (datatype=="freq") {
 if (length(freqtot)==1) freqtot<-matrix(rep(freqtot,((dim(freq1)[1])*(dim(freq1)[2]))),nrow=dim(freq1)[1])
}

if (datatype=="dataframe"){
agg<-table(data$object,data$attribute,data$rating)
freq1<-agg[,,"1"]
freqtot<-agg[,,"0"]+agg[,,"1"]}



## check input parameters, define starting values

if (sum(is.nan(freq1)|is.na(freq1)|freq1<0)){print("The elements of the matrix freq1 should be non-missing and larger than or equal to 0")}
else if (sum(is.nan(freqtot)|is.na(freqtot)|freqtot<1)){print("The elements of the matrix freqtot should be non-missing and larger than 0")}
else if (sum(abs(dim(freqtot)-dim(freq1)))>0){print("The matrices freq1 and freqtot should have the same dimensionality")}
else if (sum(freq1>freqtot)>0) {print("The corresponding values of the matrix freqtot should be larger or equal to the elements in the matrix freq1")}
else if ((maprule!="disj") & (maprule!="conj")) {print("The mapping rule is not correctly specified")}
else if (((dim(freq1)[1])*(dim(freq1)[2]))<F*sum(dim(freq1))) {print("The number of model parameters (J+K)*F should be smaller than the number of observations J*K used to fit the model")}
else if (F<1){print("The requested number of features should be larger than zero")}
else if (Nburnin<0){print("The requested number of burn-in iterations should be larger than or equal to zero")}
else if (Nchains<1){print("The requested number of chains should be larger than or equal to 1")}
else if (Nstep>maxNiter) {print("maxNiter should be larger than Nstep")}
else if (Rhatcrit<1.05){print("The convergence criterion Rhatcrit should not be smaller than 1.05")}


else {


J<-dim(freq1)[1]
K<-dim(freq1)[2]
C<-2**F
featurelabels<-paste("F",seq(1,F),sep="")
if (is.null(rownames(freq1))=="FALSE") rowlabels<-rownames(freq1)
else rowlabels<-paste("R_", seq(1,J),sep="")
if (is.null(colnames(freq1))=="FALSE") columnlabels<-colnames(freq1)
else columnlabels<-paste("C_", seq(1,K),sep="")

## define start values

if (!(is.null(fitted.plfm)) & (start.bayes=="fitted.plfm")){ 
     		if  ( (!(is.list(fitted.plfm)))) {print("fitted.plfm should be a list with components objpar and attpar")}
	else if  ( (is.list(fitted.plfm)) & (length(fitted.plfm$objpar)==0)) {print("The list fitted.plfm should contain the component objpar")}
	else if  ( (is.list(fitted.plfm)) & (length(fitted.plfm$attpar)==0)) {print("The list fitted.plfm should contain the component attpar")}
	else if  ( (is.list(fitted.plfm)) & (sum(is.nan(fitted.plfm$objpar)|is.na(fitted.plfm$objpar)|fitted.plfm$objpar<0|fitted.plfm$objpar>1)>0)) {print("The values in the matrix fitted.plfm$objpar should be non-missing and between 0 and 1")}
	else if  ( (is.list(fitted.plfm)) & (sum(is.nan(fitted.plfm$attpar)|is.na(fitted.plfm$attpar)|fitted.plfm$attpar<0|fitted.plfm$attpar>1)>0)) {print("The values in the matrix fitted.plfm$attpar should be non-missing and between 0 and 1")}
	else if  ( (is.list(fitted.plfm)) & (dim(fitted.plfm$objpar)[1]!=dim(freq1)[1])){print("The matrix fitted.plfm$objpar should have the same number of rows as the matrix freq1")}
	else if  ( (is.list(fitted.plfm)) & (dim(fitted.plfm$attpar)[1]!=dim(freq1)[2])){print("The matrix number of rows of the matrix fitted.plfm$attpar should equal the number of columns of the matrix freq1")}
	else if  ( (is.list(fitted.plfm)) & (dim(fitted.plfm$objpar)[2]!=F)){print(paste("The matrix fitted.plfm$objpar should have ", F," columns",sep=""))}
	else if ( (is.list(fitted.plfm)) & (dim(fitted.plfm$attpar)[2]!=F)){print(paste("The matrix fitted.plfm$attpar should have ", F," columns",sep=""))}
	else {
  		start.objpar<-fitted.plfm$objpar
  		start.attpar<-fitted.plfm$attpar}
}

else if (start.bayes=="best") 
 {
      postmode<-plfm(maprule=maprule,object=object,attribute=attribute,rating=rating,freq1=freq1,freqtot=freqtot,datatype=datatype,F=F,M=20,emcrit1=1e-2,emcrit2=1e-10,printrun="FALSE")
      start.objpar<-postmode$objpar
      start.attpar<-postmode$attpar
 }

else if (start.bayes=="random")
 {
  start.objpar<-matrix(runif(J*F),nrow=J)
  start.attpar<-matrix(runif(K*F),nrow=K)
 }

initial.objpar<-array(rep(0,J,F,Nchains),c(J,F,Nchains))
for (chain in 1:Nchains){initial.objpar[,,chain]<-start.objpar}
initial.attpar<-array(rep(0,K,F,Nchains),c(K,F,Nchains))
for (chain in 1:Nchains){initial.attpar[,,chain]<-start.attpar}


freq0<-freqtot-freq1
if (maprule=="conj"){
freq0<-freq1
freq1<-freqtot-freq0
initial.objpar<-1-initial.objpar
}



## generate binary matrices
binmat<-matrix(rep(0,C*F),nrow=C)
for (i in 1:(C-1)){binmat[i+1,]<-c(digitsBase(i,base=2,F))}


## define output objects

sample.sigma<-array(rep(0,J*F*maxNiter*Nchains),c(J,F,maxNiter,Nchains))
dimnames(sample.sigma)[[1]]<-rowlabels
dimnames(sample.sigma)[[2]]<-featurelabels

sample.rho<-array(rep(0,K*F*maxNiter*Nchains),c(K,F,maxNiter,Nchains))
dimnames(sample.rho)[[1]]<-columnlabels
dimnames(sample.rho)[[2]]<-featurelabels


for (chain in 1:Nchains){


## define initial values

sigma.o<-as.matrix(initial.objpar[,,chain])
rho.o<-as.matrix(initial.attpar[,,chain])


for (iter in 1:Nburnin){

##if (iter%%100==0) print(c("iteration ",iter),justify="left",quote="FALSE")
newdraw<-computedraw(J,K,F,C,binmat,sigma.o,rho.o)
sigma.o<-newdraw$sigma.n
rho.o<-newdraw$rho.n
}

for (iter in 1:Nstep){
##if (iter%%1000==0) print(c("chain ",chain,"iteration ",iter),justify="left",quote="FALSE")
newdraw<-computedraw(J,K,F,C,binmat,sigma.o,rho.o)
sigma.o<-newdraw$sigma.n
rho.o<-newdraw$rho.n
sample.sigma[,,iter,chain]<-newdraw$sigma.n 
sample.rho[,,iter,chain]<-newdraw$rho.n
}
}

## permute features to avoid label switching

if ((Nchains>1) & (F>1)){
temp<-abind(sample.sigma[,,(1:Nstep),],sample.rho[,,(1:Nstep),],along=1)
gtemp<-apply(temp,c(1,2,4),mean)
new<-temp
for (chain in 2:Nchains){
new[,,,chain]<-temp[,c(apply(cor(gtemp[,,1],gtemp[,,chain]),1,order)[F,]),,chain]}
sample.sigma[,,(1:Nstep),]<-new[(1:J),,,]
sample.rho[,,(1:Nstep),]<-new[((J+1):(J+K)),,,]
}

## compute convergence


if (Nchains==1)
{convstat<-(J+K)*F
print(c("chainlength",Nstep),justify="left",quote="FALSE")}

if (Nchains>1){
Rhat.sigma<-matrix(rep(0,J*F),nrow=J)
Rhat.rho<-matrix(rep(0,K*F),nrow=K)
rownames(Rhat.sigma)<-rowlabels
colnames(Rhat.sigma)<-featurelabels
rownames(Rhat.rho)<-columnlabels
colnames(Rhat.rho)<-featurelabels


for (j in 1:J){
for (f in 1:F){
mat<-logit(sample.sigma[j,f,(1:Nstep),])
Rhat.sigma[j,f]<-Rhat(mat)
}}

for (k in 1:K){
for (f in 1:F){
mat<-logit(sample.rho[k,f,(1:Nstep),])
Rhat.rho[k,f]<-Rhat(mat)    
}}

convstat<-sum(Rhat.sigma>Rhatcrit)+sum(Rhat.rho>Rhatcrit)
print(c("chainlength",Nstep,"# parameters to converge",convstat),justify="left",quote="FALSE")
}




stepnumber<-1

while ((convstat>0) & (((stepnumber+1)*Nstep)<=maxNiter)){



## define output step

for (chain in 1:Nchains){

## define initial values

sigma.o<-sample.sigma[,,(stepnumber*Nstep),chain]
rho.o<-sample.rho[,,(stepnumber*Nstep),chain]


for (iter in 1:Nstep){
newdraw<-computedraw(J,K,F,C,binmat,sigma.o,rho.o)
sigma.o<-newdraw$sigma.n
rho.o<-newdraw$rho.n
sample.sigma[,,(stepnumber*Nstep+iter),chain]<-newdraw$sigma.n 
sample.rho[,,(stepnumber*Nstep+iter),chain]<-newdraw$rho.n
}
}

stepnumber<-stepnumber+1

## permute features to avoid label switching

if ((Nchains>1) & (F>1)){
temp<-abind(sample.sigma[,,(1:(stepnumber*Nstep)),],sample.rho[,,(1:(stepnumber*Nstep)),],along=1)
gtemp<-apply(temp,c(1,2,4),mean)
new<-temp
for (chain in 2:Nchains){
new[,,,chain]<-temp[,c(apply(cor(gtemp[,,1],gtemp[,,chain]),1,order)[F,]),,chain]}
sample.sigma[,,(1:(stepnumber*Nstep)),]<-new[(1:J),,,]
sample.rho[,,(1:(stepnumber*Nstep)),]<-new[((J+1):(J+K)),,,]
}



## compute convergence


if (Nchains==1)
{print(c("chainlength",(stepnumber*Nstep)),justify="left",quote="FALSE")}

if (Nchains>1){
Rhat.sigma<-matrix(rep(0,J*F),nrow=J)
Rhat.rho<-matrix(rep(0,K*F),nrow=K)
rownames(Rhat.sigma)<-rowlabels
colnames(Rhat.sigma)<-featurelabels
rownames(Rhat.rho)<-columnlabels
colnames(Rhat.rho)<-featurelabels


for (j in 1:J){
for (f in 1:F){
mat<-logit(sample.sigma[j,f,(1:(stepnumber*Nstep)),])
Rhat.sigma[j,f]<-Rhat(mat)
}}

for (k in 1:K){
for (f in 1:F){
mat<-logit(sample.rho[k,f,(1:(stepnumber*Nstep)),])
Rhat.rho[k,f]<-Rhat(mat)    
}}

convstat<-sum(Rhat.sigma>Rhatcrit)+sum(Rhat.rho>Rhatcrit)
print(c("chainlength",(stepnumber*Nstep),"# parameters to converge",convstat),justify="left",quote="FALSE")

}
}




##compute posterior mean and posterior intervals
if (F==1)
{
pmean.sigma<-as.matrix(apply(sample.sigma[,,(1:(stepnumber*Nstep)),],1,mean))
pmean.rho<-as.matrix(apply(sample.rho[,,(1:(stepnumber*Nstep)),],1,mean))
p95.sigma<-as.matrix(apply(sample.sigma[,,(1:(stepnumber*Nstep)),],1,post95))
p95.rho<-as.matrix(apply(sample.rho[,,(1:(stepnumber*Nstep)),],1,post95))
}

if (F>1)
{
pmean.sigma<-apply(sample.sigma[,,(1:(stepnumber*Nstep)),],c(1,2),mean)
pmean.rho<-apply(sample.rho[,,(1:(stepnumber*Nstep)),],c(1,2),mean)
p95.sigma<-apply(sample.sigma[,,(1:(stepnumber*Nstep)),],c(1,2),post95)
p95.rho<-apply(sample.rho[,,(1:(stepnumber*Nstep)),],c(1,2),post95)
}


## compute fitmeasures


p1.mean<-pi(J,K,F,pmean.sigma,pmean.rho)
correlation<-cor(c(p1.mean*freqtot),c(freq1))
VAF<-correlation**2
fitmeasures<-as.matrix(c(correlation,VAF))
rownames(fitmeasures)<-c("Correlation observed and expected frequencies","VAF observed frequencies")
colnames(fitmeasures)<-""

if ((maprule=="disj") & (Nchains>1))
bayesplfm<-list(call=match.call(),sample.objpar=sample.sigma[,,(1:(stepnumber*Nstep)),],sample.attpar=sample.rho[,,(1:(stepnumber*Nstep)),],pmean.objpar=pmean.sigma,pmean.attpar=pmean.rho,p95.objpar=p95.sigma,p95.attpar=p95.rho,Rhat.objpar=Rhat.sigma,Rhat.attpar=Rhat.rho,fitmeasures=fitmeasures,convstat=convstat)
else if ((maprule=="disj") & (Nchains==1))
bayesplfm<-list(call=match.call(),sample.objpar=sample.sigma[,,(1:(stepnumber*Nstep)),],sample.attpar=sample.rho[,,(1:(stepnumber*Nstep)),],pmean.objpar=pmean.sigma,pmean.attpar=pmean.rho,p95.objpar=p95.sigma,p95.attpar=p95.rho,fitmeasures=fitmeasures)
else if ((maprule=="conj") & (Nchains>1))
bayesplfm<-list(call=match.call(),sample.objpar=1-sample.sigma[,,(1:(stepnumber*Nstep)),],sample.attpar=sample.rho[,,(1:(stepnumber*Nstep)),],pmean.objpar=1-pmean.sigma,pmean.attpar=pmean.rho,p95.objpar=1-p95.sigma,p95.attpar=p95.rho,Rhat.objpar=Rhat.sigma,Rhat.attpar=Rhat.rho,fitmeasures=fitmeasures,convstat=convstat)
else if ((maprule=="conj") & (Nchains==1))
bayesplfm<-list(call=match.call(),sample.objpar=1-sample.sigma[,,(1:(stepnumber*Nstep)),],sample.attpar=sample.rho[,,(1:(stepnumber*Nstep)),],pmean.objpar=1-pmean.sigma,pmean.attpar=pmean.rho,p95.objpar=1-p95.sigma,p95.attpar=p95.rho,fitmeasures=fitmeasures)

class(bayesplfm)<-"bayesplfm"
bayesplfm
}
}





############################
## function print.bayesplfm
############################


print.bayesplfm <- function(x, ...)
{
J<-dim(x$pmean.objpar)[1]
K<-dim(x$pmean.attpar)[1]
F<-dim(x$pmean.objpar)[2]
Nchains<-dim(x$sample.objpar)[4]
Npar<-(J+K)*F

cat("CALL:\n")
print(x$call)
if (is.na(Nchains)=="FALSE"){
cat("\nNUMBER OF PARAMETERS THAT DO NOT MEET CONVERGENCE CRITERION:\n")
convstat<-matrix(c(Npar,x$convstat),ncol=1)
colnames(convstat)<-""
rownames(convstat)<-c("total number of parameters","number of parameters without convergence")
print(convstat)
}

cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
print(round(x$fitmeasures,3))

cat("\nPOSTERIOR MEAN OBJECTPARAMETERS:\n")
strout<-as.matrix(capture.output(round(x$pmean.objpar,2)))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$pmean.objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\nPOSTERIOR MEAN ATTRIBUTEPARAMETERS:\n")
strout<-as.matrix(capture.output(round(x$pmean.attpar,2)))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$pmean.attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")
}


###############################
## function summary.bayesplfm
###############################

summary.bayesplfm <- function(object, ...)
{
J<-dim(object$pmean.objpar)[1]
K<-dim(object$pmean.attpar)[1]
F<-dim(object$pmean.objpar)[2]
Nchains<-dim(object$sample.objpar)[4]

descript<-round(object$fitmeasures,3)
estobjpar<-round(object$pmean.objpar,2)
estattpar<-round(object$pmean.attpar,2)
p95objpar<-matrix(rep("",J*F),nrow=J)
rownames(p95objpar)<-rownames(object$pmean.objpar)
colnames(p95objpar)<-colnames(object$pmean.objpar)
for (f in 1:F) p95objpar[,f]<-paste("[",round(object$p95.objpar[1,,f],3),";",round(object$p95.objpar[3,,f],3),"]",sep="")
p95attpar<-matrix(rep("",K*F),nrow=K)
rownames(p95attpar)<-rownames(object$pmean.attpar)
colnames(p95attpar)<-colnames(object$pmean.attpar)
for (f in 1:F) p95attpar[,f]<-paste("[",round(object$p95.attpar[1,,f],3),";",round(object$p95.attpar[3,,f],3),"]",sep="")
if (is.na(Nchains)){
result<-list(call=object$call,
             descriptivefit=descript,
             objpar=estobjpar,
             attpar=estattpar,
             p95objpar=p95objpar,
             p95attpar=p95attpar)
}
else if (Nchains>1){
Rhatobjpar<-round(object$Rhat.objpar,3)
Rhatattpar<-round(object$Rhat.attpar,3)
result<-list(call=object$call,
             descriptivefit=descript,
             objpar=estobjpar,
             attpar=estattpar,
             p95objpar=p95objpar,
             p95attpar=p95attpar,
             Rhatobjpar=Rhatobjpar,
             Rhatattpar=Rhatattpar)
}
class(result)<-"summary.bayesplfm"
result
}

####################################
## function print.summary.bayesplfm
####################################



print.summary.bayesplfm<-function(x, ...)
{
cat("Call:\n")
print(x$call)

cat("\nDESCRIPTIVE FIT OBJECT X ATTRIBUTE TABLE:\n")
print(x$descriptivefit)

cat("\nPOSTERIOR MEAN OBJECTPARAMETERS:\n")
strout<-as.matrix(capture.output(x$objpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\n95% POSTERIOR INTERVAL OBJECTPARAMETERS:\n")
strout<-as.matrix(capture.output(x$p95objpar))
strout<-gsub("\"", "", strout,fixed="TRUE")
strout<-gsub("[0.", "[.", strout,fixed="TRUE")
strout<-gsub(";0.", ";.", strout,fixed="TRUE")
rownames(strout)<-rep("",(1+dim(x$p95objpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")


if (is.null(x$Rhatobjpar)=="FALSE")
{cat("\n RHAT CONVERGENCE OBJECTPARAMETERS:\n")
print(x$Rhatobjpar)}

cat("\nPOSTERIOR MEAN ATTRIBUTEPARAMETERS:\n")
strout<-as.matrix(capture.output(x$attpar))
strout<- gsub(" 0", "  ", strout)
rownames(strout)<-rep("",(1+dim(x$attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

cat("\n95% POSTERIOR INTERVAL ATTRIBUTEPARAMETERS:\n")
strout<-as.matrix(capture.output(x$p95attpar))
strout<-gsub("\"", "", strout,fixed="TRUE")
strout<-gsub("[0.", "[.", strout,fixed="TRUE")
strout<-gsub(";0.", ";.", strout,fixed="TRUE")
rownames(strout)<-rep("",(1+dim(x$p95attpar)[1]))
colnames(strout)<-rep("",1)
print(strout,quote="FALSE")

if (is.null(x$Rhatattpar)=="FALSE")
{cat("\n RHAT CONVERGENCE ATTRIBUTEPARAMETERS:\n")
print(x$Rhatattpar)}
}



