
print.extrafrail<- function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    uni <- all(as.numeric(names(table(table(x$id)))) == 1)
    bi <- all(as.numeric(names(table(table(x$id)))) == 2)
    if (x$dist.frail == "WL") {
        if (uni) 
            cat("Univariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
    }
    if (x$dist.frail == "GA") {
        if (uni) 
            cat("Univariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "IG") {
        if (uni) 
            cat("Univariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
    }
    if (x$dist.frail == "BS") {
        if (uni) 
            cat("Univariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist=="pe", "piecewise exponential", x$dist)), 
                " survival function\n", sep = "")
    }
    cat("\n")
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
        quote = FALSE)
    invisible(x)
}

