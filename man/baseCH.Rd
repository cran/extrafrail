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
the Weibull model, this function is \eqn{\Lambda_0(t)=\lambda t^{\rho}}. 
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
data(readmission, package="frailtypack")
#Example for WL frailty model
fit.WL <- frailty.fit(survival::Surv(time, event) ~ dukes + survival::cluster(id), 
dist.frail="WL", data = readmission)
baseCH(c(500,1000,1500),fit.WL)
}
}
