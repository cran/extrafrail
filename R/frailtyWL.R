frailtyWL <-
function(formula, data, dist="np", prec=1e-4, max.iter=1000)
{
if (!inherits(formula, "formula")) {
        if (inherits(formula, "data.frame")) 
            warning("You gave a data.frame instead of a formula.\n                                                Argument order has changed; now it's emfrail(formula, data, etc..).")
        stop("formula is not an object of type formula")
    }
    if (!inherits(data, "data.frame")) {
        if (inherits(data, "formula")) 
            warning("You gave a formula instead of a data.frame.\n                                            Argument order has changed; now it's emfrail(formula, data, etc..).")
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
#if(is.null(t)) stop("t must have at least one element")
#if(is.null(delta)) stop("delta must have at least one element")
#if(is.null(cluster)) stop("cluster must have at least one element")
if (is.null(dist)) stop("dist must be specified")
if (is.null(prec)) stop("prec must be specified")
if (is.null(max.iter)) stop("max.iter must be specified")
t<-c(t); delta<-c(delta); cluster<-c(cluster)
max.iter<-round(max.iter)
if (!any(dist == c("np","weibull"))) stop("distribution is not recognized")
if(length(t)!=length(delta)) stop("t and delta don't have the same length")
if(length(t)!=length(cluster)) stop("t and cluster don't have the same length")
if(prec>1) stop("prec is too high")
if(max.iter<=0) stop("max.iter at least 1")
#x=NULL
#if(cov.formula!=~1)
#{
#x=model.matrix(cov.formula)[,-1,drop=FALSE]
#}
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
observed.llike.0.we<-function(eta, t, delta, ind)
{
rho=eta[1]
lambda=eta[2]
theta=eta[3]
Lambda0=lambda*t^rho
log.lambda0=log(lambda)+log(rho)+(rho-1)*log(t)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0, ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
P2=delta*log.lambda0
-(sum(P1)+sum(P2))
}
observed.llike.we<-function(eta, t, delta, x, ind)
{
beta=eta[1:ncol(x)]
rho=eta[ncol(x)+1]
lambda=eta[ncol(x)+2]
theta=eta[ncol(x)+3]
Lambda0=lambda*t^rho
log.lambda0=log(lambda)+log(rho)+(rho-1)*log(t)
r=fail.cluster(delta, ind)
a=theta*(theta+4)/(2*(theta+2))
b=4/(theta*(theta+4))
a.i=1/(fail.cluster(Lambda0*exp(x%*%beta), ind)+1/a)
b.i=r+b
P1=log(theta)-log(2)-(b+1)*log(a)-lgamma(b)+lgamma(b.i)+b.i*log(a.i)+log1p(a.i*b.i)
P2=delta*x%*%beta+delta*log.lambda0
-(sum(P1)+sum(P2))
}
m=length(t); theta.last=0.5
z.last=rep(1, m); kappa.last=rep(0, m)
r=fail.cluster(delta, cluster)
dif=10; i=1; rho.last=1; lambda.last=1
if(ncol(x)==0){
while(i<=max.iter & dif>prec)
{
cox.aux=survreg(Surv(t, delta)~offset(log(z.last[cluster])), dist = "weibull")
lambda.new=exp(-coef(cox.aux)/cox.aux$scale)
rho.new=1/cox.aux$scale
Lambda0.new=lambda.new*t^rho.new
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.new, cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
theta.new=optim(theta.last, prof.WL, z=z.new, kappa=kappa.new, method="Brent", lower=0.0000001, upper=20)$par
dif=max(abs(c(rho.last,lambda.last,theta.last)- c(rho.new,lambda.new,theta.new)))
theta.last=theta.new
rho.last=rho.new
lambda.last=lambda.new
kappa.last=kappa.new
z.last=z.new
i=i+1
}
aux.se=hessian(observed.llike.0.we, x0=c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.0.we(c(rho.last,lambda.last,theta.last), t=t, delta=delta, ind=cluster)
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
cox.aux=survreg(Surv(t, delta)~x+offset(log(z.last[cluster])), dist = "weibull")
beta.new=-coef(cox.aux)[-1]/cox.aux$scale
lambda.new=exp(-coef(cox.aux)[1]/cox.aux$scale)
rho.new=1/cox.aux$scale
Lambda0.new=lambda.new*t^rho.new
phi.aux=4/(theta.last*(theta.last+4))+r
alpha.aux=2*(theta.last+2)/(theta.last*(theta.last+4))+fail.cluster(Lambda0.new, cluster)
z.new=phi.aux*(alpha.aux+phi.aux+1)/(alpha.aux*(alpha.aux+phi.aux))
kappa.new=-alpha.aux/(phi.aux*(alpha.aux+phi.aux))+digamma(phi.aux+1)-log(alpha.aux)
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
aux.se=hessian(observed.llike.we, x0=c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster)
se=sqrt(diag(solve(aux.se)))
llike.obs=-observed.llike.we(c(beta.last,rho.last,lambda.last,theta.last), t=t, x=x, delta=delta, ind=cluster)
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
object.out$tau<-tau.WL(theta.last)
}
object.out
}
