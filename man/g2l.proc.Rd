\name{g2l.proc}
\alias{g2l.proc}
\alias{g2l.infer}
\alias{g2l.infer.boot}
\alias{fdr.thresh}
\alias{get_bh_threshold}
\alias{getNullProb}
\alias{Predict.LP.poly}

\title{
Procedures for global and local inference.
}
\description{
This function performs customized fdr analyses tailored to each individual cases.
}
\usage{
g2l.proc(X, z, X.target = NULL, z.target = NULL, m = c(6, 8), alpha = 0.05, 
	niter = NULL, nsample = length(z), approx.method = "direct",
	ngrid = 2000, centering = 'LP', coef.smooth='BIC', 
	fdr.method = "locfdr", plot = TRUE, rel.null = "custom", 
	locfdr.df = 10, fdr.th.fixed = NULL, parallel = TRUE)
}
\arguments{
  \item{X}{A \eqn{n}-by-\eqn{d} matrix of covariate values}
  \item{z}{A length \eqn{n} vector containing observations of z values.}
  \item{X.target}{A \eqn{k}-by-\eqn{d} matrix providing \eqn{k} sets of covariates for target cases to investigate. Set to NULL to investigate all cases and provide global inference results.}
  \item{z.target}{A vector of length \eqn{k}, providing the target \eqn{z} values to investigate}
  \item{m}{An ordered pair. First number indicates how many LP-nonparametric basis to construct for each \eqn{X}, second number indicates how many to construct for \eqn{z}. Default: \code{m=c(6,8)}.}
  \item{alpha}{Confidence level for determining signals.}
  \item{niter}{Number of iterations to use for each target case, each time a new set of relevance samples will be generated for analysis, and the resulting fdr curves are aggregated together by taking the mean values. Set to \code{NULL} to disable.}
  \item{nsample}{Number of relevance samples generated for each case.}
  \item{approx.method}{Method used to approximate customized fdr curve, default is "direct".When set to "indirect", the customized fdr is computed by modifying pooled fdr using comparison density.}
  \item{ngrid}{Number of gridpoints to use for computing customized fdr curve.}
  \item{centering}{Centering method for z with respect to covariates, default is \code{LP}, which uses LP regression; \code{lm} uses simple linear regression, and \code{spline} uses spline smoothing. Set to \code{NULL} to disable}
  \item{coef.smooth}{Specifies the method to use for LP coefficient smoothing (AIC or BIC). Uses BIC by default.}
  \item{fdr.method}{Method for controlling false discoveries (either "locfdr" or "BH"), default choice is "locfdr".}
  \item{plot}{Whether to include plots in the results, default is \code{TRUE}.}
  \item{rel.null}{Indicates how null hypothesis behaves with respect to \eqn{X}: case-relevant (\code{custom}) or fixed standard normal (\code{th}). By default the hypothesis changes with different cases of \eqn{X}.}
  \item{locfdr.df}{Degrees of freedom to use for \code{locfdr()}}
  \item{fdr.th.fixed}{Use fixed fdr threshold for finding signals. Default set to \code{NULL}, which finds different thresholds for different cases. }
  \item{parallel}{Use parallel computing for obtaining the relevance samples, mainly used for very huge \code{nsample}, default if FALSE.}
}

\value{
A list containing the following items:
	\item{macro}{Available when \code{X.target} set to \code{NULL}, contains the following items:}
		\item{result}{A list of global inference results:}
		\item{X}{Matrix of covariates, same as input \code{X}.}
		\item{z}{Vector of observations, same as input \code{z}.}
		\item{probnull}{A vector of length \eqn{n}, indicating how likely the observed z belongs to local null.}
		\item{signal}{A binary vector of length \eqn{n}, discoveries are indicated by \eqn{1}.}
		\item{}{}	
	\item{plots}{A list of plots for global inference:}
		\item{signal_x}{A plot of signals discovered, marked in red}
		\item{dps_xz}{A scatterplot of z on x, colored based on the discovery propensity scores, only available when \code{fdr.method = "locfdr"}.}
		\item{dps_x}{A scatterplot of discovery propensity scores on x, only available when \code{fdr.method = "locfdr"}.}
	\item{}{}
	\item{micro}{Available when \code{X.target} are provided with values, contains the following items:}
		\item{result}{Customized estimates for null probabilities for target \eqn{X} and \eqn{z}}
		\item{global}{Pooled global estimates for null probabilities for target \eqn{X} and \eqn{z}}
		\item{plots}{Customized fdr plots for the target cases.}
		\item{}{}
	\item{\code{m.lp}}{Same as input \code{m}}
}

\references{
Mukhopadhyay, S., and Wang, K (2019) "On The Problem of Relevance in Statistical Inference". Technical Report. 
}
\author{
Mukhopadhyay, S. and Wang, K.

Maintainer: Kaijun Wang <kaijunwang.19@gmail.com>
}

\examples{
##Toy data
set.seed(20)
x<-rep(c(1:10),rep(50,10))
sig<-sapply(x,FUN=function(x){min(x,5)})
z<-rnorm(500,0,sig)
out<-rnorm(10,10,.1)
x<-c(x,rep(1,10))
z<-c(z,out);z<-z/sd(z)

##macro inference
g2l_macro<-g2l.proc(x,z,niter=NULL,alpha=.05, nsample=500, centering=NULL,
	fdr.method = 'locfdr',parallel=FALSE)
g2l_macro$macro$plots

#micro-inference at x=1,z=2.3:
x.target=1
z.target=2.3
g2l_micro<-g2l.proc(x,z,x.target,z.target,niter = 10,m=c(4,8),alpha=.05,
	centering=NULL,parallel=FALSE)
g2l_micro$micro$result
g2l_micro$micro$global


\donttest{
data(RAZOR)
X<-RAZOR$x
z<-RAZOR$z
##macro-inference using locfdr and LASER:
g2l_macro<-g2l.proc(X,z,m=c(4,8),niter=NULL,alpha=.05,
fdr.method = 'locfdr',parallel=FALSE)
g2l_macro$macro$plots

##micro-inference on point (30,4.09), using 10 iterations:
X.target=30
z.target=4.09
g2l_micro<-g2l.proc(X,z,X.target,z.target,niter = 10,m=c(4,8),alpha=.05,parallel=FALSE)
g2l_micro$micro$result
g2l_micro$micro$global
g2l_micro$micro$plots
}
}

\keyword{ Main Functions }