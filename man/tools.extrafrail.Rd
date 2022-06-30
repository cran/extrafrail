\name{tools.extrafrail}
\alias{summary.extrafrail}
\alias{print.extrafrail}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print a summary for a object of the "extrafrail" class.
}
\description{
Summarizes the results for a object of the "extrafrail" class.
}
\usage{
\method{summary}{extrafrail}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of the "extrafrail" class.
}
  \item{\dots}{
for extra arguments.}
}
\details{
Supported frailty models are:
- weighted frailty model 
}
\value{
A complete summary for the coefficients extracted from a "extrafrail" object.  
}
\references{
Gallardo and Bourguignon (2022).
}
\author{
Diego Gallardo and Marcelo Bourguignon.}
\examples{
\donttest{
#require(frailtypack)
require(survival)
data(readmission, package="frailtypack")
fit <- frailtyWL(survival::Surv(time, event) ~ dukes + survival::cluster(id), 
data = readmission)
summary(fit)
}
}
