print.extrafrail <-
function(x, digits = max(3L, getOption("digits") - 3L), ...) 
{
uni<-all(as.numeric(names(table(table(x$id))))==1)
bi<-all(as.numeric(names(table(table(x$id))))==2)
if(uni) cat("Univariate frailty weighted Lindley model with\n",ifelse(x$dist=="np",
"non-parametric",x$dist)," survival function\n",sep="")
if(bi) cat("Bivariate frailty weighted Lindley model with\n",ifelse(x$dist=="np",
"non-parametric",x$dist)," survival function\n",sep="")
if(!uni & !bi) cat("Multivariate frailty weighted Lindley model with\n",ifelse(x$dist=="np",
"non-parametric",x$dist)," survival function\n",sep="")
    cat("\n")
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
            quote = FALSE)
    invisible(x)
}
