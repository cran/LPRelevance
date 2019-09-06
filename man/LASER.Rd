\name{LASER}
\alias{LASER}
\alias{g2l.sampler}
\title{
Generates Artificial RELevance Samples.
}
\description{
This function generates the artificial relevance samples (LASER) by nonparametrically modeling the conditional density \eqn{f(z|X=x)}.
}
\usage{
LASER(nsample = length(z), X, z, X.target, m = c(6, 8), centering = 'LP',
	coef.smooth='BIC', parallel = FALSE)
}
\arguments{
  \item{nsample}{Number of relevance samples to generate for each case.}
  \item{X}{A \eqn{n}-by-\eqn{d} matrix of covariate values}
  \item{z}{A length \eqn{n} vector containing observations of z values.}
  \item{X.target}{A \eqn{k}-by-\eqn{d} matrix providing k sets of target points for which the LASERs are required.}
  \item{m}{An ordered pair. First number indicates how many LP-nonparametric basis to construct for each \eqn{X}, second number indicates how many to construct for \eqn{z}. Default: \code{m=c(6,8)}}
  \item{centering}{Centering method for z with respect to covariates, default is \code{LP}, which uses LP regression; \code{lm} uses simple linear regression, and \code{spline} uses spline smoothing. Set to \code{NULL} to disable. Note: 
\code{spline} option only works for univariate \code{X}.}
  \item{coef.smooth}{Specifies the method to use for LP coefficient smoothing (AIC or BIC). Uses BIC by default.}
  \item{parallel}{Use parallel computing for obtaining the relevance samples, mainly used for very huge \code{nsample}, default if FALSE.}
}

\value{
A list containing the following items:
   \item{data}{The relevance sample points generated for \code{X.target}.}
   \item{LPcoef}{The LP coefficient values for \eqn{z} given \eqn{\boldsymbol{x}}.}
}
\references{
Mukhopadhyay, S., and Wang, K (2019) "On The Problem of Relevance in Statistical Inference". Technical Report. 
}
\author{
Mukhopadhyay, S. and Wang, K.

Maintainer: Kaijun Wang <kaijunwang.19@gmail.com>
}
\examples{

##Toy data:
set.seed(20)
x<-rep(c(1:10),rep(50,10))
sig<-sapply(x,FUN=function(x){min(x,5)})
z<-rnorm(500,0,sig)
out<-rnorm(10,10,.1)
x<-c(x,rep(1,10))
z<-c(z,out);z<-z/sd(z)

##LASER samples at x=1
sample.x<-LASER(X=x,z=z,X.target=1,m=c(4,8),centering=NULL,parallel=FALSE)$data
hist(sample.x,50)

\donttest{
data(RAZOR)
X<-RAZOR$x
z<-RAZOR$stat
sample.x30<-LASER(X,z,X.target=30,m=c(4,8))$data
hist(sample.x30,50)
}
}

\keyword{ Main Functions }
