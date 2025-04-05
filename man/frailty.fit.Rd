\name{frailty.fit}
\alias{frailty.fit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fitted different shared frailty models
}
\description{
frailty.fit computes the maximum likelihood estimates based on the EM algorithm for the
shared gamma, inverse gaussian, weighted Lindley, Birnbaum-Saunders, truncated normal,
mixture of inverse gaussian and mixture of Birbaum-Saunders frailty models. 
}
\usage{
frailty.fit(formula, data, dist.frail="gamma", dist = "np", prec = 1e-04, 
        max.iter = 1000, part=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
A formula that contains on the left hand side an object of the type Surv and on the right hand side a 
+cluster(id) statement, possibly with the covariates definition.
}
  \item{data}{
A data.frame in which the formula argument can be evaluated
}
  \item{dist.frail}{
the distribution assumed for the frailty. Supported values: gamma (GA also is valid), IG (inverse gaussian), WL
(weighted Lindley), BS (Birnbaum-Saunders), TN (truncated normal), MIG (mixture of IG), MBS (mixture of BS) and
GE (generalized exponential).
}
  \item{dist}{
the distribution assumed for the basal model. Supported values: weibull, pe (piecewise exponential),
exponential and np (non-parametric).
}
  \item{prec}{
The convergence tolerance for parameters.  
}
  \item{max.iter}{
The maximum number of iterations.
}
  \item{part}{
partition time (only for piecewise exponential distribution).
}
}
\details{
For the weibull, exponential and piecewise exponential distributions as the basal model, the M1-step is performed 
using the optim function. For the non-parametric case, the M1-step is based on the coxph function 
from the survival package.
}
\value{
an object of class "extrafrail" is returned. The object returned for this functions is a list containing the following
components:
\item{coefficients}{A named vector of coefficients}
\item{se}{A named vector of the standard errors for the estimated coefficients.}
\item{t}{The vector of times.}
\item{delta}{The failure indicators.}
\item{id}{A variable indicating the cluster which belongs each observation.}
\item{x}{The regressor matrix based on cov.formula (without intercept term).}
\item{dist}{The distribution assumed for the basal model.}
\item{dist.frail}{The distribution assumed for the frailty variable.}
\item{tau}{The Kendall's tau coefficient.}
\item{logLik}{The log-likelihood function (only when the Weibull model is specified for the basal distribution).}
\item{Lambda0}{The observed times and the associated cumulative hazard function (only when the non-parametric
option is specified for the basal distribution)}
\item{part}{the partition time (only for piecewise exponential model).}
}
\references{
Gallardo, D.I., Bourguignon, M. (2022) The shared weighted Lindley frailty model for
cluster failure time data. Biometrical journal, 67, e70044.

Gallardo, D.I., Bourguignon, M., Romeo, J. (2024) Birnbaum-Saunders frailty regression 
models for clustered survival data. Statistics and Computing, 34, 141.
}
\author{
Diego Gallardo, Marcelo Bourguignon and John Santibanez.
}
\examples{
\donttest{
require(survival)
#require(frailtyHL)
data(rats, package="frailtyHL")
#Fit for WL frailty model
fit.WL <- frailty.fit(survival::Surv(time, status)~ rx+ survival::cluster(litter), 
dist.frail="WL", data = rats)
summary(fit.WL)
#Fit for gamma frailty model
fit.GA <- frailty.fit(survival::Surv(time, status) ~ rx + survival::cluster(litter), 
dist.frail="gamma", data = rats)
summary(fit.GA)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
