frailty.fit <-
function(formula, data, dist.frail="gamma", dist="np", prec=1e-4, max.iter=1000)
{
if (!inherits(formula, "formula")) {
        if (inherits(formula, "data.frame")) 
            warning("You gave a data.frame instead of a formula.")
        stop("formula is not an object of type formula")
    }
    if (!inherits(data, "data.frame")) {
        if (inherits(data, "formula")) 
            warning("You gave a formula instead of a data.frame.")
        stop("data is not an object of type data.frame")
    }
if (missing(formula)) 
        stop("Missing formula")
if (missing(data)) 
        stop("Missing data")
cluster <- function(x) x
mf <- model.frame(formula, data)
    pos_cluster <- grep("cluster", names(mf))
    if (length(pos_cluster) != 1) 
        stop("misspecified or non-specified cluster")
cluster <- mf[[pos_cluster]]
 Y <- mf[[1]]
    if (!inherits(Y, "Surv")) 
        stop("left hand side not a survival object")
    X1 <- model.matrix(formula, data)
    pos_cluster_X1 <- grep("cluster", colnames(X1))
    x <- X1[, -c(1, pos_cluster_X1), drop = FALSE]
    t<-Y[, 1]; delta<-Y[, 2]
if (is.null(dist)) stop("dist must be specified")
if (is.null(prec)) stop("prec must be specified")
if (is.null(max.iter)) stop("max.iter must be specified")
t<-c(t); delta<-c(delta); cluster<-c(cluster); ind=cluster
max.iter<-round(max.iter)
if (!any(dist.frail == c("gamma","IG","WL","BS"))) stop("frailty distribution is not recognized")
if (!any(dist == c("np","weibull"))) stop("distribution for time-to-event is not recognized")
if(length(t)!=length(delta)) stop("t and delta don't have the same length")
if(length(t)!=length(cluster)) stop("t and cluster don't have the same length")
if(prec>1) stop("prec is too high")
if(max.iter<=0) stop("max.iter at least 1")
H.base<-function(t, lambda, rho, dist){
if(dist=="weibull") Lambda0<-lambda*t^rho
Lambda0}
h.base<-function(t, lambda, rho, dist, log=TRUE){
if(dist=="weibull") log.lambda0<-log(lambda)+log(rho)+(rho-1)*log(t)
log.lambda0}
frailtyGA <- function(formula, data, dist="np", prec=1e-4, max.iter=1000)
{
prof.GA<-function(theta, z, logz)
{
ll=-(1/theta)*log(theta)-lgamma(1/theta)+(1/theta-1)*logz-z/theta
-sum(ll)
}
fail.cluster<-function(delta, indice){
sum.fail<-function(ind, delta){sum(delta[which(indice==ind)])}
unlist(lapply(1:max(indice), sum.fail, delta=delta))}
tau.GA<-function(theta){theta/(theta+2)}
if(dist=="weibull")
{
observed.llike.0.ga.dist<-function(eta, t, delta, ind, dist)
{
rho=eta[1]
lambda=eta[2]
theta=eta[3]
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
a=r+1/theta
b=fail.cluster(Lambda0, ind)+1/theta
P1=-lgamma(1/theta)-(1/theta)*log(theta)+lgamma(a)-a*log(b)
P2=delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.ga.dist<-function(eta, t, delta, x, ind,dist)
{
beta=eta[1:ncol(x)]
rho=eta[ncol(x)+1]
lambda=eta[ncol(x)+2]
theta=eta[ncol(x)+3]
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
a=r+1/theta
b=fail.cluster(Lambda0*exp(x%*%beta), ind)+1/theta
P1=-lgamma(1/theta)-(1/theta)*log(theta)+lgamma(a)-a*log(b)
P2=delta*x%*%beta+delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.ga.dist.Q1.0<-function(eta, theta, z, t, delta, x, ind, dist){
if(dist=="weibull") {rho=exp(eta[1]);lambda=exp(eta[2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*log.lambda0-z[ind]*Lambda0
-sum(P1)}
observed.llike.ga.dist.Q1<-function(eta, theta, z, t, delta, x, ind, dist){
beta=eta[1:ncol(x)]
if(dist=="weibull") {rho=exp(eta[ncol(x)+1]);lambda=exp(eta[ncol(x)+2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*x%*%beta+delta*log.lambda0-z[ind]*Lambda0*exp(x%*%beta)
-sum(P1)}
m=length(t); theta.last=0.5
z.last=rep(1, m); logz.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1; rho.last=1; lambda.last=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
a.aux=r+1/theta.last
b.aux=fail.cluster(Lambda0.last, cluster)+1/theta.last
z.new=a.aux/b.aux
logz.new=digamma(a.aux)-log(b.aux)
aux.1=optim(c(log(rho.last), log(lambda.last)), observed.llike.ga.dist.Q1.0,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
if(dist=="weibull") {rho.new=exp(aux.1$par[1]);lambda.new=exp(aux.1$par[2])}
theta.new=optim(theta.last, prof.GA, z=z.new, logz=logz.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(rho.last,lambda.last,theta.last)- c(rho.new,lambda.new,theta.new)))
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.ga.dist, x0=c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.0.ga.dist(c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
para=c(rho.new,lambda.new,theta.new)
names(para)=names(se)=c("rho","lambda","theta")
}
if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=survreg(Surv(t, delta)~x, dist = "weibull")
beta.last=-coef(cox.aux)[-1]/cox.aux$scale
lambda.last=exp(-coef(cox.aux)[1]/cox.aux$scale)
rho.last=1/cox.aux$scale
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
a.aux=r+1/theta.last
b.aux=fail.cluster(Lambda0.last*exp(x%*%beta.last), cluster)+1/theta.last
z.new=a.aux/b.aux
logz.new=digamma(a.aux)-log(b.aux)
aux.1=optim(c(beta.last, log(rho.last), log(lambda.last)), observed.llike.ga.dist.Q1,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
beta.new=aux.1$par[1:ncol(x)]
if(dist=="weibull") {rho.new=exp(aux.1$par[ncol(x)+1]);lambda.new=exp(aux.1$par[ncol(x)+2])}
theta.new=optim(theta.last, prof.GA, z=z.new, logz=logz.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(beta.last,rho.last,lambda.last,theta.last)- c(beta.new,rho.new,lambda.new,theta.new)))
beta.last=beta.new
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.ga.dist, x0=c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.ga.dist(c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
para=c(beta.last,rho.last,lambda.last,theta.last)
names(para)=names(se)=c(colnames(x),"rho","lambda","theta")
}

object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"GA"
object.out$tau<-tau.GA(theta.last)
object.out$logLik<-llike.obs
}

if(dist=="np")
{
observed.llike.0.ga<-function(eta, t, delta, ind, cox.aux)
{
theta=eta
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
a=r+1/theta
b=fail.cluster(Lambda0, ind)+1/theta
P1=-lgamma(1/theta)-(1/theta)*log(theta)+lgamma(a)-a*log(b)
-sum(P1)
}
observed.llike.ga<-function(eta, t, delta, x, ind, cox.aux)
{
theta=eta[length(eta)]
beta=eta[-length(eta)]
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
a=r+1/theta
b=fail.cluster(Lambda0*exp(x%*%beta), ind)+1/theta
P1=-lgamma(1/theta)-(1/theta)*log(theta)+lgamma(a)-a*log(b)
P2=delta*x%*%beta
-(sum(P1)+sum(P2))
}
cumhazard.basal<-function(t, coxph.object){
ind.min<-function(t0, time)
{min(which(time>=t0))}
bb=basehaz(coxph.object)
tt=bb$time
bb$hazard[unlist(lapply(t, ind.min, time=tt))]}
m=length(t); theta.last=0.5
z.last=rep(1, m); logz.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~offset(log(z.last[cluster])))
Lambda0.new=cumhazard.basal(t, cox.aux)
a.aux=r+1/theta.last
b.aux=fail.cluster(Lambda0.new, cluster)+1/theta.last
z.new=a.aux/b.aux
logz.new=digamma(a.aux)-log(b.aux)
theta.new=optim(theta.last, prof.GA, z=z.new, logz=logz.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(theta.last- theta.new))
theta.last=theta.new
z.last=z.new
logz.last=logz.new
i=i+1
}
aux.se=hessian(observed.llike.0.ga, x0=c(theta.last), t=t, delta=delta, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(theta.new)
names(para)=names(se)=c("theta")
}

if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=coxph(Surv(t, delta)~x)
beta.last=coef(cox.aux)
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~x+offset(log(z.last[cluster])))
beta.new=coef(cox.aux)
Lambda0.new=cumhazard.basal(t, cox.aux)
a.aux=r+1/theta.last
b.aux=fail.cluster(Lambda0.new*exp(x%*%beta.new), cluster)+1/theta.last
z.new=a.aux/b.aux
logz.new=digamma(a.aux)-log(b.aux)
theta.new=optim(theta.last, prof.GA, z=z.new, logz=logz.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(beta.last, theta.last)-c(beta.new, theta.new)))
beta.last=beta.new
theta.last=theta.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.ga, x0=c(beta.last, theta.last), t=t, delta=delta, x=x, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(beta.new, theta.new)
names(para)=names(se)=c(colnames(x),"theta")
}

bb=basehaz(cox.aux)
Lambda0=cbind(bb$time, bb$hazard)
colnames(Lambda0)=c("time","hazard")
object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$Lambda0<-Lambda0
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"GA"
object.out$tau<-tau.GA(theta.last)
}
object.out
}
frailtyWL <- function(formula, data, dist="np", prec=1e-4, max.iter=1000)
{
prof.WL<-function(theta, z, kappa)
{
phi=4/(theta*(theta+4)); alpha=sqrt(phi*(phi+1))
ll=(phi+1)*log(alpha)-log(alpha+phi)-lgamma(phi)+(phi-1)*kappa-alpha*z
-sum(ll)
}
fail.cluster<-function(delta, indice){
sum.fail<-function(ind, delta){sum(delta[which(indice==ind)])}
unlist(lapply(1:max(indice), sum.fail, delta=delta))}
tau.WL<-function(theta){
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
aux.int<-function(x,theta){
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
x*(1+a*x)^(-2*b-4)*(1+theta*x/2)*(1+theta*(x+1)/(theta+2))}
4*a*(1+b)*integrate(aux.int, lower=0, upper=Inf, theta=theta)$value-1}
if(dist=="weibull")
{
observed.llike.0.dist<-function(eta, t, delta, ind,dist)
{
rho=eta[1]
lambda=eta[2]
theta=eta[3]
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0, ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
P2=delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.dist<-function(eta, t, delta, x, ind,dist)
{
beta=eta[1:ncol(x)]
rho=eta[ncol(x)+1]
lambda=eta[ncol(x)+2]
theta=eta[ncol(x)+3]
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0*exp(x%*%beta), ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
P2=delta*x%*%beta+delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.wl.dist.Q1.0<-function(eta, theta, z, t, delta, x, ind,dist){
if(dist=="weibull"){rho=exp(eta[1]);lambda=exp(eta[2])}
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*log.lambda0-z[ind]*Lambda0
-sum(P1)}
observed.llike.wl.dist.Q1<-function(eta, theta, z, t, delta, x, ind,dist){
beta=eta[1:ncol(x)]
if(dist=="weibull") {rho=exp(eta[ncol(x)+1]); lambda=exp(eta[ncol(x)+2])}
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=lambda*t^rho
log.lambda0=log(lambda)+log(rho)+(rho-1)*log(t)
P1=delta*x%*%beta+delta*log.lambda0-z[ind]*Lambda0*exp(x%*%beta)
-sum(P1)}
m=length(t); theta.last=0.5
z.last=rep(1, m); kappa.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1; rho.last=1; lambda.last=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.last, cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
aux.1=optim(c(log(rho.last), log(lambda.last)), observed.llike.wl.dist.Q1.0,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
if(dist=="weibull") {rho.new=exp(aux.1$par[1]); lambda.new=exp(aux.1$par[2])}
theta.new=optim(theta.last, prof.WL, z=z.new, kappa=kappa.new, method="Brent", lower=0.0000001, upper=20)$par
dif=max(abs(c(rho.last,lambda.last,theta.last)- c(rho.new,lambda.new,theta.new)))
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
kappa.last=kappa.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.dist, x0=c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.0.dist(c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
para=c(rho.new,lambda.new,theta.new)
names(para)=names(se)=c("rho","lambda","theta")
}
if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=survreg(Surv(t, delta)~x, dist = "weibull")
beta.last=-coef(cox.aux)[-1]/cox.aux$scale
lambda.last=exp(-coef(cox.aux)[1]/cox.aux$scale)
rho.last=1/cox.aux$scale
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.last*exp(x%*%beta.last), cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
aux.1=optim(c(beta.last, log(rho.last), log(lambda.last)), observed.llike.wl.dist.Q1,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
beta.new=aux.1$par[1:ncol(x)]
if(dist=="weibull") {rho.new=exp(aux.1$par[ncol(x)+1]);lambda.new=exp(aux.1$par[ncol(x)+2])}
theta.new=optim(theta.last, prof.WL, z=z.new, kappa=kappa.new, method="Brent", lower=0.0000001, upper=20)$par
dif=max(abs(c(beta.last,rho.last,lambda.last,theta.last)- c(beta.new,rho.new,lambda.new,theta.new)))
beta.last=beta.new
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
kappa.last=kappa.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.dist, x0=c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.dist(c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
para=c(beta.last,rho.last,lambda.last,theta.last)
names(para)=names(se)=c(colnames(x),"rho","lambda","theta")
}

object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"WL"
object.out$tau<-tau.WL(theta.last)
object.out$logLik<-llike.obs
}

if(dist=="np")
{
observed.llike.0<-function(eta, t, delta, ind, cox.aux)
{
theta=eta
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0, ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
-sum(P1)
}
observed.llike<-function(eta, t, delta, x, ind, cox.aux)
{
theta=eta[length(eta)]
beta=eta[-length(eta)]
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0*exp(x%*%beta), ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
P2=delta*x%*%beta
-(sum(P1)+sum(P2))
}
cumhazard.basal<-function(t, coxph.object){
ind.min<-function(t0, time)
{min(which(time>=t0))}
bb=basehaz(coxph.object)
tt=bb$time
bb$hazard[unlist(lapply(t, ind.min, time=tt))]}
m=length(t); theta.last=0.5
z.last=rep(1, m); kappa.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~offset(log(z.last[cluster])))
Lambda0.new=cumhazard.basal(t, cox.aux)
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.new, cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
theta.new=optim(theta.last, prof.WL, z=z.new, kappa=kappa.new, method="Brent", lower=0.0000001, upper=20)$par
dif=max(abs(theta.last- theta.new))
theta.last=theta.new
kappa.last=kappa.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0, x0=c(theta.last), t=t, delta=delta, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(theta.new)
names(para)=names(se)=c("theta")
}

if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=coxph(Surv(t, delta)~x)
beta.last=coef(cox.aux)
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~x+offset(log(z.last[cluster])))

beta.new=coef(cox.aux)
Lambda0.new=cumhazard.basal(t, cox.aux)
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.new*exp(x%*%beta.new), cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
theta.new=optim(theta.last, prof.WL, z=z.new, kappa=kappa.new, method="Brent", lower=0.0000001, upper=20)$par
dif=max(abs(c(beta.last, theta.last)-c(beta.new, theta.new)))
beta.last=beta.new
theta.last=theta.new
kappa.last=kappa.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike, x0=c(beta.last, theta.last), t=t, delta=delta, x=x, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(beta.new, theta.new)
names(para)=names(se)=c(colnames(x),"theta")
}

bb=basehaz(cox.aux)
Lambda0=cbind(bb$time, bb$hazard)
colnames(Lambda0)=c("time","hazard")
object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$Lambda0<-Lambda0
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"WL"
object.out$tau<-tau.WL(theta.last)
}
object.out
}





summary.extrafrail<- function(object, ...)
{
asterisk<-function (x) 
{
    if (x > 0.1) {
        ast = " "}
    else {
        if (x > 0.05) {
            ast = "."}
        else {
            if (x > 0.01) {
                ast = "*"}
            else {
                if (x > 0.001) {
                  ast = "**"}
                else {
                  {
                    ast = "***"}
                }
            }
        }
    }
    return(ast)
}
 tt<-cbind(object$coefficients,object$se,exp(object$coefficients),object$coefficients/object$se,
pnorm(abs(object$coefficients/object$se),lower.tail=FALSE))
    ast = sapply(tt[,5],FUN=asterisk)
    tt = data.frame(round(tt, 5), ast)
 colnames(tt)<-c("coef","s.e.","exp(coef)","z value","Pr(>|z|)","")
uni<-all(as.numeric(names(table(table(object$id))))==1)
bi<-all(as.numeric(names(table(table(object$id))))==2)
if(object$dist.frail=="WL")
{
if(uni) cat("Univariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(bi) cat("Bivariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(!uni & !bi) cat("Multivariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
}
if(object$dist.frail=="BS")
{
if(uni) cat("Univariate frailty Birnbaum-Saunders model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(bi) cat("Bivariate frailty Birnbaum-Saunders model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(!uni & !bi) cat("Multivariate Birnbaum-Saunders model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
}
if(length(object$coefficients)>1)
      {
if(object$dist=="np")
{
        cat("-------------------------------------------------------------------------\n")
 cat("Regression Coefficients\n")
       print(tt[-nrow(tt),,drop=FALSE])
       cat("---\n")
       cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
       cat("-------------------------------------------------------------------------\n")
 cat("Frailty variance\n")
       print(tt[nrow(tt),1:2,drop=FALSE])
cat("-------------------------------------------------------------------------\n")
}
if(object$dist=="weibull")
{
if(nrow(tt)>3)
{
        cat("-------------------------------------------------------------------------\n")
 cat("Regression Coefficients\n")
       print(tt[-c(nrow(tt)-0:2),,drop=FALSE])
       cat("---\n")
       cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
}
       cat("-------------------------------------------------------------------------\n")
 cat("Parameters of baseline distribution\n")
       print(tt[nrow(tt)-2:1,1:2,drop=FALSE])
       cat("-------------------------------------------------------------------------\n")
 cat("Frailty variance\n")
       print(tt[nrow(tt),1:2,drop=FALSE])
cat("-------------------------------------------------------------------------\n")
}
}
if(length(object$coefficients)==1)
      {cat("-------------------------------------------------------------------------\n")
 cat("Frailty variance\n")
       print(tt[nrow(tt),1:2,drop=FALSE])
cat("-------------------------------------------------------------------------\n")
}
cat(paste("Kendall's tau:",round(object$tau,4)))
cat("\n---\n")

}
frailtyIG <- function(formula, data, dist="np", prec=1e-4, max.iter=1000)
{
prof.IG<-function(theta, z, z1)
{
ll=-(1/2)*log(theta)+1/theta-(1/2)*(z/theta+z1/theta)
-sum(ll)
}
fail.cluster<-function(delta, indice){
sum.fail<-function(ind, delta){sum(delta[which(indice==ind)])}
unlist(lapply(1:max(indice), sum.fail, delta=delta))}
tau.IG<-function(theta){0.5-1/theta+(2/theta^2)*exp(2/theta)*expint_E1(2/theta)}
if(dist=="weibull")
{
observed.llike.0.ig.dist<-function(eta, t, delta, ind, dist)
{
rho=eta[1]
lambda=eta[2]
theta=eta[3]
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0, ind)+1/theta
b=1/theta
P1=0.5*log(2)+1/theta-0.5*log(pi*theta)+log(besselK(sqrt(a*b),nu=p))-(p/2)*(log(a)-log(b))
P2=delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.ig.dist<-function(eta, t, delta, x, ind, dist)
{
beta=eta[1:ncol(x)]
rho=eta[ncol(x)+1]
lambda=eta[ncol(x)+2]
theta=eta[ncol(x)+3]
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0*exp(x%*%beta), ind)+1/theta
b=1/theta
P1=0.5*log(2)+1/theta-0.5*log(pi*theta)+log(besselK(sqrt(a*b),nu=p))-(p/2)*(log(a)-log(b))
P2=delta*x%*%beta+delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.ig.dist.Q1.0<-function(eta, theta, z, t, delta, x, ind, dist){
if(dist=="weibull") {rho=exp(eta[1]);lambda=exp(eta[2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*log.lambda0-z[ind]*Lambda0
-sum(P1)}
observed.llike.ig.dist.Q1<-function(eta, theta, z, t, delta, x, ind, dist){
beta=eta[1:ncol(x)]
if(dist=="weibull") {rho=exp(eta[ncol(x)+1]);lambda=exp(eta[ncol(x)+2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*x%*%beta+delta*log.lambda0-z[ind]*Lambda0*exp(x%*%beta)
-sum(P1)}
m=length(t); theta.last=0.5
z.last=rep(1, m); z1.last=rep(1, m)
r=fail.cluster(delta, cluster)
dif=10; i=1; rho.last=1; lambda.last=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
p.aux=r-1/2
a.aux=2*fail.cluster(Lambda0.last, ind)+1/theta.last
b.aux=1/theta.last
z.new=sqrt(b.aux/a.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)
z1.new=sqrt(a.aux/b.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)-2*p.aux/b.aux
#cox.aux=survreg(Surv(t, delta)~offset(log(z.last[cluster])), dist = "weibull")
#lambda.new=exp(-coef(cox.aux)/cox.aux$scale)
#rho.new=1/cox.aux$scale
aux.1=optim(c(log(rho.last), log(lambda.last)), observed.llike.ig.dist.Q1.0,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
if(dist=="weibull") {rho.new=exp(aux.1$par[1]);lambda.new=exp(aux.1$par[2])}
theta.new=optim(theta.last, prof.IG, z=z.new, z1=z1.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(rho.last,lambda.last,theta.last)- c(rho.new,lambda.new,theta.new)))
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.ig.dist, x0=c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.0.ig.dist(c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster, dist=dist)
para=c(rho.new,lambda.new,theta.new)
names(para)=names(se)=c("rho","lambda","theta")
}
if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=survreg(Surv(t, delta)~x, dist = "weibull")
beta.last=-coef(cox.aux)[-1]/cox.aux$scale
lambda.last=exp(-coef(cox.aux)[1]/cox.aux$scale)
rho.last=1/cox.aux$scale
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
p.aux=r-1/2
a.aux=2*fail.cluster(Lambda0.last*exp(x%*%beta.last), ind)+1/theta.last
b.aux=1/theta.last
z.new=sqrt(b.aux/a.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)
z1.new=sqrt(a.aux/b.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)-2*p.aux/b.aux
#cox.aux=survreg(Surv(t, delta)~x+offset(log(z.last[cluster])), dist = "weibull")
#beta.new=-coef(cox.aux)[-1]/cox.aux$scale
#lambda.new=exp(-coef(cox.aux)[1]/cox.aux$scale)
#rho.new=1/cox.aux$scale
aux.1=optim(c(beta.last, log(rho.last), log(lambda.last)), observed.llike.ig.dist.Q1,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
beta.new=aux.1$par[1:ncol(x)]
if(dist=="weibull") {rho.new=exp(aux.1$par[ncol(x)+1]);lambda.new=exp(aux.1$par[ncol(x)+2])}
theta.new=optim(theta.last, prof.IG, z=z.new, z1=z1.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(beta.last,rho.last,lambda.last,theta.last)- c(beta.new,rho.new,lambda.new,theta.new)))
beta.last=beta.new
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.ig.dist, x0=c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.ig.dist(c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster, dist=dist)
para=c(beta.last,rho.last,lambda.last,theta.last)
names(para)=names(se)=c(colnames(x),"rho","lambda","theta")
}

object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"IG"
object.out$tau<-tau.IG(theta.last)
object.out$logLik<-llike.obs
}

if(dist=="np")
{
observed.llike.0.ig<-function(eta, t, delta, ind, cox.aux)
{
theta=eta
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0, ind)+1/theta
b=1/theta
P1=0.5*log(2)+1/theta-0.5*log(pi*theta)+log(besselK(sqrt(a*b),nu=p))-(p/2)*(log(a)-log(b))
-sum(P1)
}
observed.llike.ig<-function(eta, t, delta, x, ind, cox.aux)
{
theta=eta[length(eta)]
beta=eta[-length(eta)]
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0*exp(x%*%beta), ind)+1/theta
b=1/theta
P1=0.5*log(2)+1/theta-0.5*log(pi*theta)+log(besselK(sqrt(a*b),nu=p))-(p/2)*(log(a)-log(b))
P2=delta*x%*%beta
-(sum(P1)+sum(P2))
}
cumhazard.basal<-function(t, coxph.object){
ind.min<-function(t0, time)
{min(which(time>=t0))}
bb=basehaz(coxph.object)
tt=bb$time
bb$hazard[unlist(lapply(t, ind.min, time=tt))]}
m=length(t); theta.last=0.5
z.last=rep(1, m); logz.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~offset(log(z.last[cluster])))
Lambda0.new=cumhazard.basal(t, cox.aux)
p.aux=r-1/2
a.aux=2*fail.cluster(Lambda0.new, ind)+1/theta.last
b.aux=1/theta.last
z.new=sqrt(b.aux/a.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)
z1.new=sqrt(a.aux/b.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)-2*p.aux/b.aux
theta.new=optim(theta.last, prof.IG, z=z.new, z1=z1.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(theta.last- theta.new))
theta.last=theta.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.ig, x0=c(theta.last), t=t, delta=delta, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(theta.new)
names(para)=names(se)=c("theta")
}

if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=coxph(Surv(t, delta)~x)
beta.last=coef(cox.aux)
while(i<=max.iter & dif>prec)
{
cox.aux=coxph(Surv(t, delta)~x+offset(log(z.last[cluster])))
beta.new=coef(cox.aux)
Lambda0.new=cumhazard.basal(t, cox.aux)
p.aux=r-1/2
a.aux=2*fail.cluster(Lambda0.new*exp(x%*%beta.new), ind)+1/theta.last
b.aux=1/theta.last
z.new=sqrt(b.aux/a.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)
z1.new=sqrt(a.aux/b.aux)*besselK(sqrt(a.aux*b.aux),nu=p.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p.aux)-2*p.aux/b.aux
theta.new=optim(theta.last, prof.IG, z=z.new, z1=z1.new, method="Brent", lower=0.0001, upper=20)$par
dif=max(abs(c(beta.last, theta.last)-c(beta.new, theta.new)))
beta.last=beta.new
theta.last=theta.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.ig, x0=c(beta.last, theta.last), t=t, delta=delta, x=x, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(beta.new, theta.new)
names(para)=names(se)=c(colnames(x),"theta")
}

bb=basehaz(cox.aux)
Lambda0=cbind(bb$time, bb$hazard)
colnames(Lambda0)=c("time","hazard")
object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$Lambda0<-Lambda0
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"IG"
object.out$tau<-tau.IG(theta.last)
}
object.out
}

frailtyBS <- function(formula, data, dist="np", prec=1e-4, max.iter=1000)
{
if (!inherits(formula, "formula")) {
        if (inherits(formula, "data.frame")) 
            warning("You gave a data.frame instead of a formula.")
        stop("formula is not an object of type formula")
    }
    if (!inherits(data, "data.frame")) {
        if (inherits(data, "formula")) 
            warning("You gave a formula instead of a data.frame.")
        stop("data is not an object of type data.frame")
    }
if (missing(formula)) 
        stop("Missing formula")
if (missing(data)) 
        stop("Missing data")
cluster <- function(x) x
mf <- model.frame(formula, data)
    pos_cluster <- grep("cluster", names(mf))
    if (length(pos_cluster) != 1) 
        stop("misspecified or non-specified cluster")
cluster <- mf[[pos_cluster]]
 Y <- mf[[1]]
    if (!inherits(Y, "Surv")) 
        stop("left hand side not a survival object")
    X1 <- model.matrix(formula, data)
    pos_cluster_X1 <- grep("cluster", colnames(X1))
    x <- X1[, -c(1, pos_cluster_X1), drop = FALSE]
    t<-Y[, 1]; delta<-Y[, 2]
if (is.null(dist)) stop("dist must be specified")
if (is.null(prec)) stop("prec must be specified")
if (is.null(max.iter)) stop("max.iter must be specified")
t<-c(t); delta<-c(delta); cluster<-c(cluster); ind=cluster
max.iter<-round(max.iter)
if (!any(dist == c("np","weibull"))) stop("distribution is not recognized")
if(length(t)!=length(delta)) stop("t and delta don't have the same length")
if(length(t)!=length(cluster)) stop("t and cluster don't have the same length")
if(prec>1) stop("prec is too high")
if(max.iter<=0) stop("max.iter at least 1")
prof.BS<-function(theta, z, z1, y)
{
phi=(1+sqrt(3*theta+1)-theta)/theta
ll=(1/2)*(phi+log1p(phi))+(1-y)*(log(phi)-log1p(phi))-(1/4)*(z*(1+phi)+z1*phi^2/(1+phi))
-sum(ll)
}
fail.cluster<-function(delta, indice){
sum.fail<-function(ind, delta){sum(delta[which(indice==ind)])}
unlist(lapply(1:max(indice), sum.fail, delta=delta))}
tau.BS<-function(theta){
aux.int<-function(x, theta){
phi<-(1+sqrt(3*theta+1)-theta)/theta
ll<--log(2)+log1p((1+4*x/(phi+1))^(-1/2))+(phi/2)*(1-(1+4*x/(phi+1))^(1/2))
L.z<-exp(ll)
p1<-exp((phi/2)*(1-(1+4*x/(phi+1))^(1/2)))
p2<-(phi*(phi^2+phi*(4*x+3)+2)*sqrt(1+4*x/(phi+1))+phi^3+phi^2*(4*x-5)-2*phi*(12*x+5)-4)*sqrt(1+4*x/(phi+1))+4*(2*(phi^2+phi*(4*x+3)+2)*sqrt(1+4*x/(phi+1))+phi*(phi+4*x+1))
L2.z<-p1*p2/(2*(phi+4*x+1)^3)
x*L.z*L2.z}
4*integrate(aux.int, lower=0, upper=Inf, theta=theta)$value-1}
if(dist=="weibull")
{
observed.llike.0.bs.dist<-function(eta, t, delta, ind)
{
rho=eta[1]
lambda=eta[2]
theta=eta[3]
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0, ind)+(phi+1)/2
b=phi^2/(2*(phi+1))
P1=phi/2+(1/2)*log1p(phi)-log(2)-0.5*log(pi)-(p/2)*(log(a)-log(b))+log(sqrt(b/a)*besselK(sqrt(a*b),nu=p+1)+(phi/(1+phi))*besselK(sqrt(a*b),nu=p))
P2=delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.bs.dist<-function(eta, t, delta, x, ind)
{
beta=eta[1:ncol(x)]
rho=eta[ncol(x)+1]
lambda=eta[ncol(x)+2]
theta=eta[ncol(x)+3]
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0*exp(x%*%beta), ind)+(phi+1)/2
b=phi^2/(2*(phi+1))
P1=phi/2+(1/2)*log1p(phi)-log(2)-0.5*log(pi)-(p/2)*(log(a)-log(b))+log(sqrt(b/a)*besselK(sqrt(a*b),nu=p+1)+(phi/(1+phi))*besselK(sqrt(a*b),nu=p))
P2=delta*x%*%beta+delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.bs.dist.Q1.0<-function(eta, theta, z, t, delta, x, ind, dist){
if(dist=="weibull") {rho=exp(eta[1]);lambda=exp(eta[2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*log.lambda0-z[ind]*Lambda0
-sum(P1)}
observed.llike.bs.dist.Q1<-function(eta, theta, z, t, delta, x, ind, dist){
beta=eta[1:ncol(x)]
if(dist=="weibull") {rho=exp(eta[ncol(x)+1]);lambda=exp(eta[ncol(x)+2])}
Lambda0=H.base(t, lambda, rho, dist)
log.lambda0=h.base(t, lambda, rho, dist)
P1=delta*x%*%beta+delta*log.lambda0-z[ind]*Lambda0*exp(x%*%beta)
-sum(P1)}
m=length(t); theta.last=0.5
z.last=rep(1, m); z1.last=rep(1, m); y.last=rep(0.5, m)
r=fail.cluster(delta, cluster)
dif=10; i=1; rho.last=1; lambda.last=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
phi.last=(1+sqrt(3*theta.last+1)-theta.last)/theta.last
p0.aux=r-1/2
p1.aux=r+1-1/2
a.aux=2*fail.cluster(Lambda0.last, cluster)+(phi.last+1)/2
b.aux=phi.last^2/(2*(phi.last+1))
xi0.aux=(phi.last/(1+phi.last))*besselK(sqrt(a.aux*b.aux), nu=p0.aux)/(a.aux/b.aux)^(p0.aux/2)
xi1.aux=besselK(sqrt(a.aux*b.aux), nu=p1.aux)/(a.aux/b.aux)^(p1.aux/2)
nu.aux=xi1.aux/(xi0.aux+xi1.aux)
y.new=nu.aux
z.new=sqrt(b.aux/a.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))
z1.new=sqrt(a.aux/b.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))-(2/b.aux)*((1-nu.aux)*p0.aux+nu.aux*p1.aux)
aux.1=optim(c(log(rho.last), log(lambda.last)), observed.llike.bs.dist.Q1.0,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
if(dist=="weibull") {rho.new=exp(aux.1$par[1]);lambda.new=exp(aux.1$par[2])}
theta.new=optim(theta.last, prof.BS, z=z.new, z1=z1.new, y=y.new, method="Brent", lower=0.0001, upper=4.9999)$par
dif=max(abs(c(rho.last,lambda.last,theta.last)- c(rho.new,lambda.new,theta.new)))
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.bs.dist, x0=c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.0.bs.dist(c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster)
para=c(rho.new,lambda.new,theta.new)
names(para)=names(se)=c("rho","lambda","theta")
}
if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=survreg(Surv(t, delta)~x, dist = "weibull")
beta.last=-coef(cox.aux)[-1]/cox.aux$scale
lambda.last=exp(-coef(cox.aux)[1]/cox.aux$scale)
rho.last=1/cox.aux$scale
while(i<=max.iter & dif>prec)
{
phi.last=(1+sqrt(3*theta.last+1)-theta.last)/theta.last
Lambda0.last=H.base(t, lambda.last, rho.last, dist)
p0.aux=r-1/2
p1.aux=r+1-1/2
a.aux=2*fail.cluster(Lambda0.last*exp(x%*%beta.last), cluster)+(phi.last+1)/2
b.aux=phi.last^2/(2*(phi.last+1))
xi0.aux=(phi.last/(1+phi.last))*besselK(sqrt(a.aux*b.aux), nu=p0.aux)/(a.aux/b.aux)^(p0.aux/2)
xi1.aux=besselK(sqrt(a.aux*b.aux), nu=p1.aux)/(a.aux/b.aux)^(p1.aux/2)
nu.aux=xi1.aux/(xi0.aux+xi1.aux)
y.new=nu.aux
z.new=sqrt(b.aux/a.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))
z1.new=sqrt(a.aux/b.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))-(2/b.aux)*((1-nu.aux)*p0.aux+nu.aux*p1.aux)
aux.1=optim(c(beta.last, log(rho.last), log(lambda.last)), observed.llike.bs.dist.Q1,
method="BFGS", theta=theta.last, z=z.new, t=t, delta=delta, x=x, ind=ind, dist=dist)
beta.new=aux.1$par[1:ncol(x)]
if(dist=="weibull") {rho.new=exp(aux.1$par[ncol(x)+1]);lambda.new=exp(aux.1$par[ncol(x)+2])}
theta.new=optim(theta.last, prof.BS, z=z.new, z1=z1.new, y=y.new, method="Brent", lower=0.0001, upper=4.9999)$par
dif=max(abs(c(beta.last,rho.last,lambda.last,theta.last)- c(beta.new,rho.new,lambda.new,theta.new)))
beta.last=beta.new
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.bs.dist, x0=c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.bs.dist(c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster)
para=c(beta.last,rho.last,lambda.last,theta.last)
names(para)=names(se)=c(colnames(x),"rho","lambda","theta")
}

object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"BS"
object.out$tau<-tau.BS(theta.last)
object.out$logLik<-llike.obs
}

if(dist=="np")
{
observed.llike.0.bs<-function(eta, t, delta, ind, cox.aux)
{
theta=eta
phi=(1+sqrt(3*theta+1)-theta)/theta
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0, ind)+(phi+1)/2
b=phi^2/(2*(phi+1))
P1=phi/2+(1/2)*log1p(phi)-log(2)-0.5*log(pi)-(p/2)*(log(a)-log(b))+log(sqrt(b/a)*besselK(sqrt(a*b),nu=p+1)+(phi/(1+phi))*besselK(sqrt(a*b),nu=p))
-sum(P1)
}
observed.llike.bs<-function(eta, t, delta, x, ind, cox.aux)
{
theta=eta[length(eta)]
phi=(1+sqrt(3*theta+1)-theta)/theta
beta=eta[-length(eta)]
Lambda0=cumhazard.basal(t, cox.aux)
r=fail.cluster(delta, ind)
p=r-1/2
a=2*fail.cluster(Lambda0*exp(x%*%beta), ind)+(phi+1)/2
b=phi^2/(2*(phi+1))
P1=phi/2+(1/2)*log1p(phi)-log(2)-0.5*log(pi)-(p/2)*(log(a)-log(b))+log(sqrt(b/a)*besselK(sqrt(a*b),nu=p+1)+(phi/(1+phi))*besselK(sqrt(a*b),nu=p))
P2=delta*x%*%beta
-(sum(P1)+sum(P2))
}
cumhazard.basal<-function(t, coxph.object){
ind.min<-function(t0, time)
{min(which(time>=t0))}
bb=basehaz(coxph.object)
tt=bb$time
bb$hazard[unlist(lapply(t, ind.min, time=tt))]}
m=length(t); theta.last=0.5
z.last=rep(1, m)
r=fail.cluster(delta, cluster)
dif=10; i=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
phi.last=(1+sqrt(3*theta.last+1)-theta.last)/theta.last
cox.aux=coxph(Surv(t, delta)~offset(log(z.last[cluster])))
Lambda0.new=cumhazard.basal(t, cox.aux)
p0.aux=r-1/2
p1.aux=r+1-1/2
a.aux=2*fail.cluster(Lambda0.new, cluster)+(phi.last+1)/2
b.aux=phi.last^2/(2*(phi.last+1))
xi0.aux=(phi.last/(1+phi.last))*besselK(sqrt(a.aux*b.aux), nu=p0.aux)/(a.aux/b.aux)^(p0.aux/2)
xi1.aux=besselK(sqrt(a.aux*b.aux), nu=p1.aux)/(a.aux/b.aux)^(p1.aux/2)
nu.aux=xi1.aux/(xi0.aux+xi1.aux)
y.new=nu.aux
z.new=sqrt(b.aux/a.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))
z1.new=sqrt(a.aux/b.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))-(2/b.aux)*((1-nu.aux)*p0.aux+nu.aux*p1.aux)
theta.new=optim(theta.last, prof.BS, z=z.new, z1=z1.new, y=y.new, method="Brent", lower=0.0001, upper=4.9999)$par
dif=max(abs(theta.last- theta.new))
theta.last=theta.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.bs, x0=c(theta.last), t=t, delta=delta, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(theta.new)
names(para)=names(se)=c("theta")
}

if(ncol(x)>0){
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
cox.aux=coxph(Surv(t, delta)~x)
beta.last=coef(cox.aux)
while(i<=max.iter & dif>prec)
{
phi.last=(1+sqrt(3*theta.last+1)-theta.last)/theta.last
cox.aux=coxph(Surv(t, delta)~x+offset(log(z.last[cluster])))
beta.new=coef(cox.aux)
Lambda0.new=cumhazard.basal(t, cox.aux)
p0.aux=r-1/2
p1.aux=r+1-1/2
a.aux=2*fail.cluster(Lambda0.new*exp(x%*%beta.new), cluster)+(phi.last+1)/2
b.aux=phi.last^2/(2*(phi.last+1))
xi0.aux=(phi.last/(1+phi.last))*besselK(sqrt(a.aux*b.aux), nu=p0.aux)/(a.aux/b.aux)^(p0.aux/2)
xi1.aux=besselK(sqrt(a.aux*b.aux), nu=p1.aux)/(a.aux/b.aux)^(p1.aux/2)
nu.aux=xi1.aux/(xi0.aux+xi1.aux)
y.new=nu.aux
z.new=sqrt(b.aux/a.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))
z1.new=sqrt(a.aux/b.aux)*((1-nu.aux)*besselK(sqrt(a.aux*b.aux),nu=p0.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p0.aux)+nu.aux*besselK(sqrt(a.aux*b.aux),nu=p1.aux+1)/besselK(sqrt(a.aux*b.aux),nu=p1.aux))-(2/b.aux)*((1-nu.aux)*p0.aux+nu.aux*p1.aux)
theta.new=optim(theta.last, prof.BS, z=z.new, z1=z1.new, y=y.new, method="Brent", lower=0.0001, upper=4.9999)$par
dif=max(abs(c(beta.last, theta.last)-c(beta.new, theta.new)))
beta.last=beta.new
theta.last=theta.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.bs, x0=c(beta.last, theta.last), t=t, delta=delta, x=x, ind=cluster, cox.aux=cox.aux)
se=sqrt(diag(solve(aux.se)))
para=c(beta.new, theta.new)
names(para)=names(se)=c(colnames(x),"theta")
}

bb=basehaz(cox.aux)
Lambda0=cbind(bb$time, bb$hazard)
colnames(Lambda0)=c("time","hazard")
object.out<-list(coefficients=para, se=se, z=z.new)
class(object.out)<-"extrafrail"
object.out$t<-t
object.out$delta<-delta
object.out$id<-cluster
object.out$Lambda0<-Lambda0
object.out$x<-x
object.out$dist<-dist
object.out$dist.frail<-"BS"
object.out$tau<-tau.BS(theta.last)
}
object.out
}
if(dist.frail=="gamma") val<-frailtyGA(formula, data, dist=dist, prec=prec, max.iter=max.iter)
if(dist.frail=="IG")    val<-frailtyIG(formula, data, dist=dist, prec=prec, max.iter=max.iter)
if(dist.frail=="WL")    val<-frailtyWL(formula, data, dist=dist, prec=prec, max.iter=max.iter)
if(dist.frail=="BS")    val<-frailtyBS(formula, data, dist=dist, prec=prec, max.iter=max.iter)
val
}
