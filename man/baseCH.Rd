\name{baseCH}
\alias{baseCH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Computes the baseline cumulative hazard function. 
}
\description{
Provides the baseline cumulative hazard function (\eqn{\Lambda_0}) for an object with extrafrail class.
}
\usage{
baseCH(t, fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{
the vector of times for which the baseline cumulative hazard function should be computed.
}
  \item{fit}{
an object with extrafrail class. 
}
}
\details{
Provides the baseline cumulative hazard function. When the baseline distribution is assumed as
the Weibull model, this function is \eqn{\Lambda_0(t)=\lambda t^{\rho}}. For the piecewise
exponential model, this function is \eqn{\Lambda_0(t)=\sum_{j=1}^L \lambda_j \nabla_j(t)},
where \eqn{\nabla_j(t)=0, } if \eqn{t<a_{j-1}}, \eqn{\nabla_j(t)=t-a_{j-1}, } if \eqn{a_{j-1}\leq t < a_{j}}
and \eqn{\nabla_j(t)=a_j-a_{j-1}, } if \eqn{t\geq a_{j}}, with \eqn{a=(a_0=0, a_1, \ldots, a_{j-1}),} the 
corresponding partition time.
}
\value{
a vector with the same length that t, including the baseline cumulative hazard function related to t.
}
\references{
Gallardo, D.I., Bourguignon, M. (2022) The multivariate weighted Lindley frailty model for
cluster failure time data. Submitted.
}
\author{
Diego Gallardo and Marcelo Bourguignon.
}
\examples{
\donttest{
#require(frailtypack)
require(survival)
data(rats, package="frailtyHL")
#Example for WL frailty model
fit.WL <- frailty.fit(survival::Surv(time, status) ~ rx + survival::cluster(litter), 
dist.frail="WL", data = rats)
baseCH(c(80,90,100),fit.WL)
}
}
