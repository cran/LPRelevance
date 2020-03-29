\name{LASER}
\alias{LASER}
\alias{g2l.sampler}
\alias{z.lp.center}
\title{
Generates Artificial RELevance Samples.
}
\description{
This function generates the artificial relevance samples (LASER). These are "sharpened" z-samples manufactured by the relevance-function \eqn{d_{Y|X=x}(u)}.
}
\usage{
LASER( X,z, X.target, m=c(4,6), nsample=length(z), lp.reg.method='lm',
       coef.smooth='BIC', centering=TRUE,parallel=FALSE,...)
}
\arguments{
  \item{X}{A \eqn{n}-by-\eqn{d} matrix of covariate values}
  \item{z}{A length \eqn{n} vector containing observations of z values.}
  \item{X.target}{A \eqn{k}-by-\eqn{d} matrix providing k sets of target points for which the LASERs are required.}
  \item{m}{An ordered pair. First number indicates how many LP-nonparametric basis to construct for each \eqn{X}, second number indicates how many to construct for \eqn{z}. Default: \code{m=c(4,6)}}
  \item{nsample}{Number of relevance samples to generate for each case.}
  \item{lp.reg.method}{Method for estimating the relevance function and its conditional LP-Fourier coefficients. We currently support thee options: lm (inbuilt with subset selection), glmnet, and knn.}
  \item{centering}{Whether to perform regression-adjustment to center the data, default is TRUE.}
  \item{coef.smooth}{Specifies the method to use for LP coefficient smoothing (AIC or BIC). Uses BIC by default.}
  \item{parallel}{Use parallel computing for obtaining the relevance samples, mainly used for very huge \code{nsample}, default if FALSE.}
  \item{...}{Extra parameters to pass to other functions. Currently only supports the arguments for \code{knn()}.}
}

\value{
A list containing the following items:
   \item{data}{The relevance sample points generated for \code{X.target}.}
   \item{LPcoef}{The LP coefficient values for \eqn{z} given \eqn{\boldsymbol{x}}.}
}
\references{
Mukhopadhyay, S., and Wang, K (2020) "On The Problem of Relevance in Statistical Inference". Technical Report. 
}
\author{
Mukhopadhyay, S. and Wang, K.

Maintainer: Kaijun Wang <kwang2@fredhutch.org>
}
\examples{
\donttest{
data(funnel)
X<-funnel$x
z<-funnel$stat
z.laser.x30<-LASER(X,z,X.target=30,m=c(4,8))$data
hist(z.laser.x30,50)
}
}

\keyword{ Main Functions }
