

summary.extrafrail <- function (object, ...) 
{
    asterisk <- function(x) {
        if (x > 0.1) {
            ast = " "
        }
        else {
            if (x > 0.05) {
                ast = "."
            }
            else {
                if (x > 0.01) {
                  ast = "*"
                }
                else {
                  if (x > 0.001) {
                    ast = "**"
                  }
                  else {
                    {
                      ast = "***"
                    }
                  }
                }
            }
        }
        return(ast)
    }
    tt <- cbind(object$coefficients, object$se, exp(object$coefficients), 
        object$coefficients/object$se, pnorm(abs(object$coefficients/object$se), 
            lower.tail = FALSE))
    ast = sapply(tt[, 5], FUN = asterisk)
    tt = data.frame(round(tt, 5), ast)
    colnames(tt) <- c("coef", "s.e.", "exp(coef)", "z value", 
        "Pr(>|z|)", "")
    uni <- all(as.numeric(names(table(table(object$id)))) == 
        1)
    bi <- all(as.numeric(names(table(table(object$id)))) == 2)
    if (object$dist.frail == "WL") {
        if (uni) 
            cat("Univariate weighted Lindley frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate weighted Lindley frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate weighted Lindley frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "TN") {
        if (uni) 
            cat("Univariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "GA") {
        if (uni) 
            cat("Univariate gamma frailty model with\n", ifelse(object$dist == 
                "np", "non-parametric", ifelse(object$dist == 
                "pe", "piecewise exponential", object$dist)), 
                " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate gamma frailty model with\n", ifelse(object$dist == 
                "np", "non-parametric", ifelse(object$dist == 
                "pe", "piecewise exponential", object$dist)), 
                " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate gamma frailty model with\n", ifelse(object$dist == 
                "np", "non-parametric", ifelse(object$dist == 
                "pe", "piecewise exponential", object$dist)), 
                " survival function\n", sep = "")
    }
    if (object$dist.frail == "IG") {
        if (uni) 
            cat("Univariate inverse gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate inverse gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate inverse gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "BS") {
        if (uni) 
            cat("Univariate Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "TN") {
        if (uni) 
            cat("Univariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate truncated normal frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "MIG") {
        if (uni) 
            cat("Univariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate mixture of inverse Gaussian frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "MBS") {
        if (uni) 
            cat("Univariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate mixture of Birnbaum-Saunders frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (object$dist.frail == "GE") {
        if (uni) 
            cat("Univariate generalized exponential frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (bi) 
            cat("Bivariate generalized exponential frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
        if (!uni & !bi) 
            cat("Multivariate generalized exponential frailty model with\n", 
                ifelse(object$dist == "np", "non-parametric", 
                  ifelse(object$dist == "pe", "piecewise exponential", 
                    object$dist)), " survival function\n", sep = "")
    }
    if (length(object$coefficients) > 1) {
        if (object$dist == "np") {
            cat("-------------------------------------------------------------------------\n")
            cat("Regression Coefficients\n")
            print(tt[-nrow(tt), , drop = FALSE])
            cat("---\n")
            cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
            cat("-------------------------------------------------------------------------\n")
            cat("Frailty variance\n")
            print(tt[nrow(tt), 1:2, drop = FALSE])
            cat("-------------------------------------------------------------------------\n")
        }
        if (object$dist != "np") {
            if (object$dist == "weibull") {
                if (nrow(tt) > 3) {
                  cat("-------------------------------------------------------------------------\n")
                  cat("Regression Coefficients\n")
                  print(tt[-c(nrow(tt) - 0:2), , drop = FALSE])
                  cat("---\n")
                  cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
                }
            }
            if (object$dist == "exponential") {
                if (nrow(tt) > 2) {
                  cat("-------------------------------------------------------------------------\n")
                  cat("Regression Coefficients\n")
                  print(tt[-c(nrow(tt) - 0:1), , drop = FALSE])
                  cat("---\n")
                  cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
                }
            }
            if (object$dist == "pe") {
                if (nrow(tt) > length(object$part) + 1) {
                  cat("-------------------------------------------------------------------------\n")
                  cat("Regression Coefficients\n")
                  print(tt[-c(nrow(tt) - 0:(length(object$part))), 
                    , drop = FALSE])
                  cat("---\n")
                  cat("Signif. codes:  0 \"***\" 0.001 \"**\" 0.01 \"*\" 0.05 \".\" 0.1 \" \" 1\n")
                }
            }
            cat("-------------------------------------------------------------------------\n")
            cat("Parameters of baseline distribution\n")
            ss = ifelse(object$dist == "pe", length(object$par), 
                ifelse(object$dist == "exponential", 1, 2))
            print(tt[nrow(tt) - ss:1, 1:2, drop = FALSE])
            cat("-------------------------------------------------------------------------\n")
            cat("Frailty variance\n")
            print(tt[nrow(tt), 1:2, drop = FALSE])
            cat("-------------------------------------------------------------------------\n")
        }
    }
    if (length(object$coefficients) == 1) {
        cat("-------------------------------------------------------------------------\n")
        cat("Frailty variance\n")
        print(tt[nrow(tt), 1:2, drop = FALSE])
        cat("-------------------------------------------------------------------------\n")
    }
    cat(paste("Kendall's tau:", round(object$tau, 4)))
    cat("\n---\n")
}
    