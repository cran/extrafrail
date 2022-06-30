\name{rWL}
\alias{rWL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generated random variables from the weighted Lindley distribution. 
}
\description{
Generated random variables from the weighted Lindley distribution with mean 1. 
}
\usage{
rWL(n, theta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of observations. If length(n) > 1, the length is taken to be the number required.
}
  \item{theta}{
variance of the variable.
}
}
\details{
The weighted Lindley distribution has probability density function 
\deqn{
f(z;\theta)=\frac{\theta}{2\Gamma(\theta)}a_{\theta}^{-b_{\theta}-1}z^{b_{\theta}-1}(1+z)\exp\left(-\frac{z}{a_{\theta}}\right), \quad z, \theta>0,
}
where \eqn{a_{\theta}=\frac{\theta(\theta+4)}{2(\theta+2)}} and \eqn{b_{\theta}=\frac{4}{\theta(\theta+4)}}. Under
this parametrization, E(Z)=1 and Var(Z)=\eqn{\theta}.
}
\value{
a vector of length n with the generated values.
}
\references{
Gallardo, D.I., Bourguignon, M. (2022) The multivariate weighted Lindley frailty model for
cluster failure time data. Submitted.
}
\author{
Diego Gallardo and Marcelo Bourguignon.
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
rWL(10, theta=0.5)
}

