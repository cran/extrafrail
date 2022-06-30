baseCH <-
function(t, fit)
{
if(!inherits(fit,"extrafrail")) stop("fit must be an object of extrafrail class")
if(!is.numeric(t)) stop("t must be a numeric vector")
t<-as.vector(t)
t[which(t<0)]=0
if(fit$dist=="np")
{
Lambda0<-fit$Lambda0
ind.min<-function(t0, time) min(which(time>=t0))
tt=Lambda0[,1]
t[which(t>max(tt))]=max(tt)
ss=Lambda0[,2][unlist(lapply(t, ind.min, time=tt))]
}
if(fit$dist=="weibull")
{
lambda<-fit$coefficients[length(fit$coefficients)-2:1]
ss=lambda[1]*t^lambda[2]
}
ss
}
