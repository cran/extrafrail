summary.extrafrail <-
function(object, ...)
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
    ast = sapply(tt[,4],FUN=asterisk)
    tt = data.frame(round(tt, 5), ast)
 colnames(tt)<-c("coef","s.e.","exp(coef)","z value","Pr(>|z|)","")
uni<-all(as.numeric(names(table(table(object$id))))==1)
bi<-all(as.numeric(names(table(table(object$id))))==2)
if(uni) cat("Univariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(bi) cat("Bivariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
if(!uni & !bi) cat("Multivariate frailty weighted Lindley model with\n",ifelse(object$dist=="np",
"non-parametric",object$dist)," survival function\n",sep="")
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
