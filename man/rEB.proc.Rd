\name{rEB.proc}
\alias{rEB.proc}
\alias{Finite.rEB}
\alias{LP.post.conv}

\title{
Relevance-Integrated Empirical Bayes Inference
}
\description{
Performs custom-tailored empirical Bayes inference with relevant samples.
}
\usage{
rEB.proc(X, z, X.target, z.target, m = c(6, 8), niter = NULL, centering = 'LP', 
	coef.smooth='BIC',nsample = length(z), theta.set.prior = NULL, 
	theta.set.post = NULL, LP.type = "L2", g.method = "DL", 
	sd0 = NULL, m.EB = 8, parallel = FALSE, avg.method = "mean", 
	post.curve = "HPD", post.alpha = 0.8, color = "red")
}
\arguments{
  \item{X}{A \eqn{n}-by-\eqn{d} matrix of covariate values}
  \item{z}{A length \eqn{n} vector containing observations of target random variable.}
  \item{X.target}{A length \eqn{d} vector providing the set of covariates for the target case. }
  \item{z.target}{the target \eqn{z} to investigate}
  \item{m}{An ordered pair. First number indicates how many LP-nonparametric basis to construct for each \eqn{X}, second number indicates how many to construct for \eqn{z}.}
  \item{niter}{Number of iterations to use for Finite Bayes, set to \code{NULL} to disable.}
    \item{centering}{Centering method for z with respect to covariates, default is \code{LP}, which uses LP regression; \code{lm} uses simple linear regression, and \code{spline} uses spline smoothing. Set to \code{NULL} to disable}
  \item{coef.smooth}{Specifies the method to use for LP coefficient smoothing (AIC or BIC). Uses BIC by default.}
  \item{nsample}{Number of relevance samples generated for the target case.}
  \item{theta.set.prior}{This indicates the set of grid points to compute prior density.}
  \item{theta.set.post}{This indicates the set of grid points to compute posterior density.}
  \item{LP.type}{User selects either "L2" for LP-orthogonal series representation of comparison density \eqn{d} or "MaxEnt" for the maximum entropy representation. Default is L2.}
  \item{g.method}{Determines the method to find \eqn{\tau^2}: "DL" uses Dersimonian and Lard technique,"SJ" uses Sidik-Jonkman}
  \item{sd0}{Fixed standard error for \eqn{z|\theta}. Default is NULL, the standard error will be calculated from data.}
  \item{m.EB}{The truncation point reflecting the concentration of true nonparametric prior density \eqn{\pi} around known prior distribution \eqn{g}}
  \item{parallel}{Use parallel computing for obtaining the relevance samples, mainly used for very huge \code{nsample}, default if FALSE.}
  \item{avg.method}{For Finite Bayes, this specifies how the results from different iterations are aggregated. ("\code{mean}" or "\code{median}".)}
  \item{post.curve}{For plotting, this specifies what to show on posterior curve. "\code{HPD}" provides HPD interval, "\code{band}" gives confidence band.}
  \item{post.alpha}{Confidence level to use when plotting posterior confidence band, or the alpha level for HPD interval.}
  \item{color}{The color of the plots.}
}
\value{
A list containing the following items:
	\item{result}{contains the results for prior and posterior density:}
		\item{prior}{Prior results:}
		\item{g.par}{Parameters for \eqn{g}.} 
		\item{LP.coef}{reports the LP coefficient values for \eqn{z} given \eqn{\boldsymbol{x}}.}
		\item{}{}
		\item{posterior}{Posterior results: }
		\item{post.mean}{Posterior mean for \eqn{\pi(\theta|\boldsymbol{x})}.}
		\item{post.mean.sd}{Standard error for the posterior mean. Only available for Finite Bayes.}
		\item{HPD.interval}{The HPD interval for posterior \eqn{\pi(\theta|\boldsymbol{x})}.}
		\item{post.alpha}{same as input \code{post.alpha}.}
		\item{}{}
	\item{plots}{The plots for prior and posterior density.}
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
x<-rep(c(1:5),rep(50,5))
sig<-sapply(x,FUN=function(x){min(x,5)})
z<-rnorm(250,0,sig)
out<-rnorm(5,5,.1)
x<-c(x,rep(1,5))
z<-c(z,out);z<-z/sd(z)

x.target=1
z.target=2.3
rEB.out<-rEB.proc(x,z,x.target,z.target,m=c(4,8),nsample=50,
	coef.smooth='AIC', centering=NULL,m.EB=4,parallel=FALSE)
rEB.out$plots$rEB.prior
rEB.out$plots$rEB.post

\donttest{
data(RAZOR)
X<-RAZOR$x
z<-RAZOR$stat
X.target=30
z.target=4.09
rEB.out<-rEB.proc(X,z,X.target,z.target,m=c(4,8),
	theta.set.prior=seq(-2,2,length.out=200),
	theta.set.post=seq(-2,5,length.out=200),
	centering=TRUE,m.EB=6,parallel=FALSE)
rEB.out$plots$rEB.post
rEB.out$plots$rEB.prior
}
}
\keyword{ Main Functions }
