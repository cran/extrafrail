rWL <-
function(n, theta=1)
{
if (is.null(n)) stop("sample size must be specified")
if (is.null(theta)) stop("theta must be specified")
if (round(n) != n | n <= 0) stop("sample size must be a positive integer")
if (theta<=0) stop("theta must be positive")
phi=4/(theta*(theta+4))
alpha=sqrt(phi*(phi+1))
y=rbinom(n, size=1, prob=alpha/(alpha+phi))
x1=rgamma(n, shape=phi, rate=alpha)
x2=rgamma(n, shape=phi+1, rate=alpha)
y*x1+(1-y)*x2
}
