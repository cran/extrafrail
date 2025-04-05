print.extrafrail<-function (x, digits = max(3L, getOption("digits") - 3L), ...) 
{
    uni <- all(as.numeric(names(table(table(x$id)))) == 1)
    bi <- all(as.numeric(names(table(table(x$id)))) == 2)
    if (x$dist.frail == "WL") {
        if (uni) 
            cat("Univariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate weighted Lindley frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "GA") {
        if (uni) 
            cat("Univariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist == "pe", 
                "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist == "pe", 
                "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate gamma frailty model with\n", ifelse(x$dist == 
                "np", "non-parametric", ifelse(x$dist == "pe", 
                "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "IG") {
        if (uni) 
            cat("Univariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate inverse gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "BS") {
        if (uni) 
            cat("Univariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "TN") {
        if (uni) 
            cat("Univariate truncated normal frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate truncated normal frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate truncated normal frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "MIG") {
        if (uni) 
            cat("Univariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "MBS") {
        if (uni) 
            cat("Univariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    if (x$dist.frail == "GE") {
        if (uni) 
            cat("Univariate generalized exponential frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (bi) 
            cat("Bivariate generalized exponential frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
        if (!uni & !bi) 
            cat("Multivariate generalized exponential frailty model with\n", 
                ifelse(x$dist == "np", "non-parametric", ifelse(x$dist == 
                  "pe", "piecewise exponential", x$dist)), " survival function\n", 
                sep = "")
    }
    cat("\n")
    cat("Coefficients:\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2L, 
        quote = FALSE)
    invisible(x)
}