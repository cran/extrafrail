frailty.fit<- function (formula, data, dist.frail = "gamma", dist = "np", prec = 1e-04, 
    max.iter = 1000, part = NULL) 
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
    t <- Y[, 1]
    delta <- Y[, 2]
    if (is.null(dist)) 
        stop("dist must be specified")
    if (is.null(prec)) 
        stop("prec must be specified")
    if (is.null(max.iter)) 
        stop("max.iter must be specified")
    t <- c(t)
    delta <- c(delta)
    cluster <- c(cluster)
    ind = cluster
    max.iter <- round(max.iter)
    if (!any(dist.frail == c("gamma", "GA", "IG", "WL", "BS", 
        "TN", "MIG","MBS"))) 
        stop("frailty distribution is not recognized")
    if (dist.frail == "gamma") 
        dist.frail = "GA"
    if (!any(dist == c("np", "weibull", "pe", "exponential"))) 
        stop("distribution for time-to-event is not recognized")
    if (length(t) != length(delta)) 
        stop("t and delta don't have the same length")
    if (length(t) != length(cluster)) 
        stop("t and cluster don't have the same length")
    if (prec > 1) 
        stop("prec is too high")
    if (max.iter <= 0) 
        stop("max.iter at least 1")
    H.base <- function(t, lambda, rho = 1, dist, part = NULL) {
        if (dist == "weibull") {
            Lambda0 <- lambda * t^rho
        }
        if (dist == "pe") {
            Lambda0 <- -ppexp(t, rate = lambda, t = part, lower.tail = FALSE, 
                log.p = TRUE)
        }
        Lambda0
    }
    h.base <- function(t, lambda, rho = 1, dist, log = TRUE, 
        part = NULL) {
        if (dist == "weibull") {
            log.lambda0 <- log(lambda) + log(rho) + (rho - 1) * 
                log(t)
        }
        if (dist == "pe") {
            log.lambda0 <- dpexp(t, rate = lambda, t = part, 
                log = TRUE) - ppexp(t, rate = lambda, t = part, 
                lower.tail = FALSE, log.p = TRUE)
        }
        log.lambda0
    }
    frailtyGA <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
        prof.GA <- function(theta, z, logz) {
            ll = -(1/theta) * log(theta) - lgamma(1/theta) + 
                (1/theta - 1) * logz - z/theta
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        tau.GA <- function(theta) {
            theta/(theta + 2)
        }
        if (dist == "weibull") {
            observed.llike.0.ga.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                theta = eta[3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0, ind) + 1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ga.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                theta = eta[ncol(x) + 3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind) + 
                  1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ga.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.ga.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            logz.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.last, cluster) + 
                    1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.ga.dist.Q1.0, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(rho.last, lambda.last, theta.last) - 
                    c(rho.new, lambda.new, theta.new)))
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ga.dist, x0 = c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.ga.dist(c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                para = c(rho.new, lambda.new, theta.new)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) + 1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.ga.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    theta.last) - c(beta.new, rho.new, lambda.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ga.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.ga.dist(c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                para = c(beta.last, rho.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "GA"
            object.out$tau <- tau.GA(theta.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.ga.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                lambda = eta[1:length(part)]
                theta = eta[length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0, ind) + 1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ga.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                theta = eta[ncol(x) + length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind) + 
                  1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ga.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.ga.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            logz.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
	    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.last, cluster) + 
                    1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  aux.1 = optim(log(lambda.last), observed.llike.ga.dist.Q1.0, 
                    method = "BFGS", theta = theta.last, z = z.new, 
                    t = t, delta = delta, x = x, ind = ind, dist = dist, 
                    part = part)
                  lambda.new = exp(aux.1$par)
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(lambda.last, theta.last) - 
                    c(lambda.new, theta.new)))
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ga.dist, x0 = c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.ga.dist(c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                para = c(lambda.new, theta.new)
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) + 1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.ga.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part = part)
                  beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, lambda.last, theta.last) - 
                    c(beta.new, lambda.new, theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ga.dist, x0 = c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.ga.dist(c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                para = c(beta.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "GA"
            object.out$tau <- tau.GA(theta.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0.ga <- function(eta, t, delta, ind, 
                cox.aux) {
                theta = eta
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0, ind) + 1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                -sum(P1)
            }
            observed.llike.ga <- function(eta, t, delta, x, ind, 
                cox.aux) {
                theta = eta[length(eta)]
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                a = r + 1/theta
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind) + 
                  1/theta
                P1 = -lgamma(1/theta) - (1/theta) * log(theta) + 
                  lgamma(a) - a * log(b)
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            logz.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(z.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.new, cluster) + 
                    1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(theta.last - theta.new))
                  theta.last = theta.new
                  z.last = z.new
                  logz.last = logz.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ga, x0 = c(theta.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(theta.new)
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(z.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  a.aux = r + 1/theta.last
                  b.aux = fail.cluster(Lambda0.new * exp(x %*% 
                    beta.new), cluster) + 1/theta.last
                  z.new = a.aux/b.aux
                  logz.new = digamma(a.aux) - log(b.aux)
                  theta.new = optim(theta.last, prof.GA, z = z.new, 
                    logz = logz.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, theta.last) - c(beta.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ga, x0 = c(beta.last, 
                  theta.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(beta.new, theta.new)
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "GA"
            object.out$tau <- tau.GA(theta.last)
        }
        object.out
    }
    frailtyWL <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
        prof.WL <- function(theta, z, kappa) {
            phi = 4/(theta * (theta + 4))
            alpha = sqrt(phi * (phi + 1))
            ll = (phi + 1) * log(alpha) - log(alpha + phi) - 
                lgamma(phi) + (phi - 1) * kappa - alpha * z
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        tau.WL <- function(theta) {
            a = theta * (theta + 4)/(2 * (theta + 2))
            b = 4/(theta * (theta + 4))
            aux.int <- function(x, theta) {
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                x * (1 + a * x)^(-2 * b - 4) * (1 + theta * x/2) * 
                  (1 + theta * (x + 1)/(theta + 2))
            }
            4 * a * (1 + b) * integrate(aux.int, lower = 0, upper = Inf, 
                theta = theta)$value - 1
        }
        if (dist == "weibull") {
            observed.llike.0.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                theta = eta[3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0, ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.dist <- function(eta, t, delta, x, 
                ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                theta = eta[ncol(x) + 3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.wl.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.wl.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = lambda * t^rho
                log.lambda0 = log(lambda) + log(rho) + (rho - 
                  1) * log(t)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            kappa.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.last, 
                    cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.wl.dist.Q1.0, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(c(rho.last, lambda.last, theta.last) - 
                    c(rho.new, lambda.new, theta.new)))
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.dist, x0 = c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.dist(c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                para = c(rho.new, lambda.new, theta.new)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.last * 
                    exp(x %*% beta.last), cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.wl.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    theta.last) - c(beta.new, rho.new, lambda.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.dist(c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                para = c(beta.last, rho.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "WL"
            object.out$tau <- tau.WL(theta.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                lambda = eta[1:length(part)]
                theta = eta[length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0, ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.dist <- function(eta, t, delta, x, 
                ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                theta = eta[ncol(x) + length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.wl.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.wl.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            kappa.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
	    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.last, 
                    cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  aux.1 = optim(log(lambda.last), observed.llike.wl.dist.Q1.0, 
                    method = "BFGS", theta = theta.last, z = z.new, 
                    t = t, delta = delta, x = x, ind = ind, dist = dist, 
                    part = part)
                  lambda.new = exp(aux.1$par[1:length(part)])
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(c(lambda.last, theta.last) - 
                    c(lambda.new, theta.new)))
                  theta.last = theta.new
                  lambda.last = lambda.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.dist, x0 = c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.dist(c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                para = c(lambda.new, theta.new)
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(exp(-coef(cox.aux)[1]/cox.aux$scale), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.last * 
                    exp(x %*% beta.last), cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.wl.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part = part)
                  beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, lambda.last, theta.last) - 
                    c(beta.new, lambda.new, theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  lambda.last = lambda.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.dist, x0 = c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.dist(c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                para = c(beta.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "WL"
            object.out$tau <- tau.WL(theta.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0 <- function(eta, t, delta, ind, 
                cox.aux) {
                theta = eta
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0, ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                -sum(P1)
            }
            observed.llike <- function(eta, t, delta, x, ind, 
                cox.aux) {
                theta = eta[length(eta)]
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                a = theta * (theta + 4)/(2 * (theta + 2))
                b = 4/(theta * (theta + 4))
                a.i = 1/(fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/a)
                b.i = r + b
                P1 = log(theta) - log(2) - (b + 1) * log(a) - 
                  lgamma(b) + lgamma(b.i) + b.i * log(a.i) + 
                  log1p(a.i * b.i)
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            kappa.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(z.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.new, 
                    cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(theta.last - theta.new))
                  theta.last = theta.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0, x0 = c(theta.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(theta.new)
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(z.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  phi.aux = 4/(theta.last * (theta.last + 4)) + 
                    r
                  alpha.aux = 2 * (theta.last + 2)/(theta.last * 
                    (theta.last + 4)) + fail.cluster(Lambda0.new * 
                    exp(x %*% beta.new), cluster)
                  z.new = phi.aux * (alpha.aux + phi.aux + 1)/(alpha.aux * 
                    (alpha.aux + phi.aux))
                  kappa.new = -alpha.aux/(phi.aux * (alpha.aux + 
                    phi.aux)) + digamma(phi.aux + 1) - log(alpha.aux)
                  theta.new = optim(theta.last, prof.WL, z = z.new, 
                    kappa = kappa.new, method = "Brent", lower = 1e-07, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, theta.last) - c(beta.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  kappa.last = kappa.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike, x0 = c(beta.last, 
                  theta.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(beta.new, theta.new)
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "WL"
            object.out$tau <- tau.WL(theta.last)
        }
        object.out
    }
    frailtyIG <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
        prof.IG <- function(theta, z, z1) {
            ll = -(1/2) * log(theta) + 1/theta - (1/2) * (z/theta + 
                z1/theta)
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        tau.IG <- function(theta) {
            0.5 - 1/theta + (2/theta^2) * exp(2/theta) * expint_E1(2/theta)
        }
        if (dist == "weibull") {
            observed.llike.0.ig.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                theta = eta[3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ig.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                theta = eta[ncol(x) + 3]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ig.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.ig.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            z1.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last, ind) + 
                    1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.ig.dist.Q1.0, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(rho.last, lambda.last, theta.last) - 
                    c(rho.new, lambda.new, theta.new)))
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ig.dist, x0 = c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.ig.dist(c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                para = c(rho.new, lambda.new, theta.new)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), ind) + 1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.ig.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    theta.last) - c(beta.new, rho.new, lambda.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ig.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.ig.dist(c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster, dist = dist)
                para = c(beta.last, rho.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "IG"
            object.out$tau <- tau.IG(theta.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.ig.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                lambda = eta[1:length(part)]
                theta = eta[length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ig.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                theta = eta[ncol(x) + length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.ig.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.ig.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            z1.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
	    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last, ind) + 
                    1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  aux.1 = optim(log(lambda.last), observed.llike.ig.dist.Q1.0, 
                    method = "BFGS", theta = theta.last, z = z.new, 
                    t = t, delta = delta, x = x, ind = ind, dist = dist, 
                    part = part)
                  lambda.new = exp(aux.1$par[1:length(part)])
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(lambda.last, theta.last) - 
                    c(lambda.new, theta.new)))
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ig.dist, x0 = c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.ig.dist(c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                para = c(lambda.new, theta.new)
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), ind) + 1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.ig.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part = part)
                  beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, lambda.last, theta.last) - 
                    c(beta.new, lambda.new, theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ig.dist, x0 = c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.ig.dist(c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                para = c(beta.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "IG"
            object.out$tau <- tau.IG(theta.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0.ig <- function(eta, t, delta, ind, 
                cox.aux) {
                theta = eta
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                -sum(P1)
            }
            observed.llike.ig <- function(eta, t, delta, x, ind, 
                cox.aux) {
                theta = eta[length(eta)]
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + 1/theta
                b = 1/theta
                P1 = 0.5 * log(2) + 1/theta - 0.5 * log(pi * 
                  theta) + log(besselK(sqrt(a * b), nu = p)) - 
                  (p/2) * (log(a) - log(b))
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            logz.last = rep(0, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(z.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.new, ind) + 
                    1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(theta.last - theta.new))
                  theta.last = theta.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.ig, x0 = c(theta.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(theta.new)
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(z.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  p.aux = r - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.new * exp(x %*% 
                    beta.new), ind) + 1/theta.last
                  b.aux = 1/theta.last
                  z.new = sqrt(b.aux/a.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux)
                  z1.new = sqrt(a.aux/b.aux) * besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux + 1)/besselK(sqrt(a.aux * 
                    b.aux), nu = p.aux) - 2 * p.aux/b.aux
                  theta.new = optim(theta.last, prof.IG, z = z.new, 
                    z1 = z1.new, method = "Brent", lower = 1e-04, 
                    upper = 20)$par
                  dif = max(abs(c(beta.last, theta.last) - c(beta.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.ig, x0 = c(beta.last, 
                  theta.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(beta.new, theta.new)
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "IG"
            object.out$tau <- tau.IG(theta.last)
        }
        object.out
    }
    frailtyBS <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
        prof.BS <- function(theta, z, z1, y) {
            phi = (1 + sqrt(3 * theta + 1) - theta)/theta
            ll = (1/2) * (phi + log1p(phi)) + (1 - y) * (log(phi) - 
                log1p(phi)) - (1/4) * (z * (1 + phi) + z1 * phi^2/(1 + 
                phi))
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        tau.BS <- function(theta) {
            aux.int <- function(x, theta) {
                phi <- (1 + sqrt(3 * theta + 1) - theta)/theta
                ll <- -log(2) + log1p((1 + 4 * x/(phi + 1))^(-1/2)) + 
                  (phi/2) * (1 - (1 + 4 * x/(phi + 1))^(1/2))
                L.z <- exp(ll)
                p1 <- exp((phi/2) * (1 - (1 + 4 * x/(phi + 1))^(1/2)))
                p2 <- (phi * (phi^2 + phi * (4 * x + 3) + 2) * 
                  sqrt(1 + 4 * x/(phi + 1)) + phi^3 + phi^2 * 
                  (4 * x - 5) - 2 * phi * (12 * x + 5) - 4) * 
                  sqrt(1 + 4 * x/(phi + 1)) + 4 * (2 * (phi^2 + 
                  phi * (4 * x + 3) + 2) * sqrt(1 + 4 * x/(phi + 
                  1)) + phi * (phi + 4 * x + 1))
                L2.z <- p1 * p2/(2 * (phi + 4 * x + 1)^3)
                x * L.z * L2.z
            }
            4 * integrate(aux.int, lower = 0, upper = Inf, theta = theta)$value - 
                1
        }
        if (dist == "weibull") {
            observed.llike.0.bs.dist <- function(eta, t, delta, 
                ind, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                theta = eta[3]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.bs.dist <- function(eta, t, delta, 
                x, ind, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                theta = eta[ncol(x) + 3]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.bs.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.bs.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            z1.last = rep(1, m)
            y.last = rep(0.5, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last, cluster) + 
                    (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.bs.dist.Q1.0, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(c(rho.last, lambda.last, theta.last) - 
                    c(rho.new, lambda.new, theta.new)))
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.bs.dist, x0 = c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.bs.dist(c(rho.last, 
                  lambda.last, theta.last), t = t, delta = delta, 
                  ind = cluster)
                para = c(rho.new, lambda.new, theta.new)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) + (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.bs.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    theta.last) - c(beta.new, rho.new, lambda.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.bs.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.bs.dist(c(beta.last, 
                  rho.last, lambda.last, theta.last), t = t, 
                  x = x, delta = delta, ind = cluster)
                para = c(beta.last, rho.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "BS"
            object.out$tau <- tau.BS(theta.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.bs.dist <- function(eta, t, delta, 
                ind, part = NULL) {
                lambda = eta[1:length(part)]
                theta = eta[length(part) + 1]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.bs.dist <- function(eta, t, delta, 
                x, ind, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                theta = eta[ncol(x) + length(part) + 1]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.bs.dist.Q1.0 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.bs.dist.Q1 <- function(eta, theta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            z1.last = rep(1, m)
            y.last = rep(0.5, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
	    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last, cluster) + 
                    (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  aux.1 = optim(log(lambda.last), observed.llike.bs.dist.Q1.0, 
                    method = "BFGS", theta = theta.last, z = z.new, 
                    t = t, delta = delta, x = x, ind = ind, dist = dist, 
                    part = part)
                  lambda.new = exp(aux.1$par[1:length(part)])
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(c(lambda.last, theta.last) - 
                    c(lambda.new, theta.new)))
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.bs.dist, x0 = c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.bs.dist(c(lambda.last, 
                  theta.last), t = t, delta = delta, ind = cluster, 
                  part = part)
                para = c(lambda.new, theta.new)
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) + (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.bs.dist.Q1, method = "BFGS", 
                    theta = theta.last, z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part = part)
                  beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(c(beta.last, lambda.last, theta.last) - 
                    c(beta.new, lambda.new, theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  lambda.last = lambda.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.bs.dist, x0 = c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.bs.dist(c(beta.last, 
                  lambda.last, theta.last), t = t, x = x, delta = delta, 
                  ind = cluster, part = part)
                para = c(beta.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "BS"
            object.out$tau <- tau.BS(theta.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0.bs <- function(eta, t, delta, ind, 
                cox.aux) {
                theta = eta
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0, ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                -sum(P1)
            }
            observed.llike.bs <- function(eta, t, delta, x, ind, 
                cox.aux) {
                theta = eta[length(eta)]
                phi = (1 + sqrt(3 * theta + 1) - theta)/theta
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                p = r - 1/2
                a = 2 * fail.cluster(Lambda0 * exp(x %*% beta), 
                  ind) + (phi + 1)/2
                b = phi^2/(2 * (phi + 1))
                P1 = phi/2 + (1/2) * log1p(phi) - log(2) - 0.5 * 
                  log(pi) - (p/2) * (log(a) - log(b)) + log(sqrt(b/a) * 
                  besselK(sqrt(a * b), nu = p + 1) + (phi/(1 + 
                  phi)) * besselK(sqrt(a * b), nu = p))
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            m = length(t)
            theta.last = 0.5
            z.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(z.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.new, cluster) + 
                    (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(theta.last - theta.new))
                  theta.last = theta.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.bs, x0 = c(theta.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(theta.new)
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  phi.last = (1 + sqrt(3 * theta.last + 1) - 
                    theta.last)/theta.last
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(z.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  p0.aux = r - 1/2
                  p1.aux = r + 1 - 1/2
                  a.aux = 2 * fail.cluster(Lambda0.new * exp(x %*% 
                    beta.new), cluster) + (phi.last + 1)/2
                  b.aux = phi.last^2/(2 * (phi.last + 1))
                  xi0.aux = (phi.last/(1 + phi.last)) * besselK(sqrt(a.aux * 
                    b.aux), nu = p0.aux)/(a.aux/b.aux)^(p0.aux/2)
                  xi1.aux = besselK(sqrt(a.aux * b.aux), nu = p1.aux)/(a.aux/b.aux)^(p1.aux/2)
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  z.new = sqrt(b.aux/a.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux))
                  z1.new = sqrt(a.aux/b.aux) * ((1 - nu.aux) * 
                    besselK(sqrt(a.aux * b.aux), nu = p0.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p0.aux) + 
                    nu.aux * besselK(sqrt(a.aux * b.aux), nu = p1.aux + 
                      1)/besselK(sqrt(a.aux * b.aux), nu = p1.aux)) - 
                    (2/b.aux) * ((1 - nu.aux) * p0.aux + nu.aux * 
                      p1.aux)
                  theta.new = optim(theta.last, prof.BS, z = z.new, 
                    z1 = z1.new, y = y.new, method = "Brent", 
                    lower = 1e-04, upper = 4.9999)$par
                  dif = max(abs(c(beta.last, theta.last) - c(beta.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  z.last = z.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.bs, x0 = c(beta.last, 
                  theta.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                para = c(beta.new, theta.new)
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = z.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "BS"
            object.out$tau <- tau.BS(theta.last)
        }
        object.out
    }
    frailtyTN <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
        ratio <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, 
            log.p = TRUE))
        prof.TN <- function(nu, z, z2) {
            ll = log1p(nu + ratio(nu) - 1) - 0.5 * (z2 * (nu + 
                ratio(nu))^2 - 2 * z * (nu + ratio(nu)) * nu + 
                nu^2) - pnorm(nu, log.p = TRUE)
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        vari.tn <- function(nu) {
            (nu + ratio(nu))^(-2) - ratio(nu) * (nu + ratio(nu))^(-1)
        }
        C.k.aux <- function(x, k, ri, bi, nu) {
            exp((ri + k) * log(x) - x * bi - x^2 * 0.5 * (nu + 
                ratio(nu))^2)
        }
        C.k <- function(k, ri, bi, nu) {
            integrate(C.k.aux, lower = 0, upper = 10, k = k, 
                ri = ri, bi = bi, nu = nu, abs.tol = 1e-15)$value
        }
        dg <- function(x) {
            -2 * (x + ratio(x))^(-3) * (1 - ratio(x) * (x + ratio(x))) + 
                ratio(x) + ratio(x) * (x + ratio(x))^(-2) * (1 - 
                ratio(x) * (x + ratio(x)))
        }
        tau.TN.aux <- function(x, theta) {
            eq.vari.tn = function(nu, theta) {
                vari.tn(nu) - theta
            }
            nu = uniroot(eq.vari.tn, lower = -20, upper = 20, 
                theta = theta, extendInt = "yes")$root
            b = nu + exp(dnorm(nu, log = TRUE) - pnorm(nu, log.p = TRUE))
            L = exp(pnorm(nu - x/b, log.p = TRUE) - pnorm(nu, 
                log.p = TRUE) + (x/b) * (-nu + x/(2 * b)))
            L2 = L/b^2 * ((nu - x/b) * (exp(dnorm(nu - x/b, log = TRUE) - 
                pnorm(nu - x/b, log.p = TRUE)) + nu - x/b) + 
                1)
            pmax(0, x * L^2/b^2 * ((nu - x/b) * (exp(dnorm(nu - 
                x/b, log = TRUE) - pnorm(nu - x/b, log.p = TRUE)) + 
                nu - x/b) + 1))
        }
        L.deriv <- function(s, nu, n) {
            gamma = nu + ratio(nu)
            alpha = nu - s/gamma
            L0 = exp(pnorm(alpha, log.p = TRUE) - pnorm(nu, log.p = TRUE) + 
                (s/gamma) * (0.5 * s/gamma - nu))
            L1 = -(L0/gamma) * (alpha + ratio(alpha))
            L2 = (L0/gamma^2) * (alpha * (ratio(alpha) + alpha) + 
                1)
            if (n == 0) 
                L = L0
            if (n == 1) 
                L = L1
            if (n == 2) 
                L = L2
            if (n > 2) {
                L = c(L1, L2)
                for (i in 3:n) {
                  L[i] = (i - 1) * L[i - 2]/gamma^2 - alpha * 
                    L[i - 1]/gamma
                }
                L = L[n]
            }
            L
        }
        tau.TN = function(theta) 4 * integrate(tau.TN.aux, lower = 0, 
            upper = 1000, theta = theta, stop.on.error = FALSE)$value - 
            1
        if (dist == "weibull") {
            observed.llike.0.tn.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                nu = eta[3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0, ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.tn.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                nu = eta[ncol(x) + 3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.tn.dist.Q1.0 <- function(eta, z, t, 
                delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.tn.dist.Q1 <- function(eta, z, t, 
                delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            nu.last = 2
            theta.last = vari.tn(nu.last)
            z.last = rep(1, m)
            z2.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  b.aux = fail.cluster(Lambda0.last, cluster) - 
                    nu.last * (nu.last + ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.tn.dist.Q1.0, method = "BFGS", 
                    z = z.new, t = t, delta = delta, x = x, ind = ind, 
                    dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = max(-40, 
                      nu.last - 3), upper = min(40, nu.last + 
                      3))$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(c(rho.last, lambda.last, theta.last) - 
                    c(rho.new, lambda.new, theta.new)))
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.tn.dist, x0 = c(rho.last, 
                  lambda.last, nu.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                llike.obs = -observed.llike.0.tn.dist(c(rho.last, 
                  lambda.last, nu.last), t = t, delta = delta, 
                  ind = cluster, dist = dist)
                para = c(rho.new, lambda.new, theta.new)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
                  b.aux = fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) - nu.last * (nu.last + 
                    ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.tn.dist.Q1, method = "BFGS", 
                    z = z.new, t = t, delta = delta, x = x, ind = ind, 
                    dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = -40, 
                    upper = 40)$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    theta.last) - c(beta.new, rho.new, lambda.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.tn.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, nu.last), t = t, x = x, 
                  delta = delta, ind = cluster, dist = dist)
                se = sqrt(diag(solve(aux.se)))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                llike.obs = -observed.llike.tn.dist(c(beta.last, 
                  rho.last, lambda.last, nu.last), t = t, x = x, 
                  delta = delta, ind = cluster, dist = dist)
                para = c(beta.new, rho.new, lambda.new, theta.new)
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new, z2 = z2.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "TN"
            object.out$tau <- tau.TN(theta.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.tn.dist <- function(eta, t, delta, 
                ind, dist, part = NULL) {
                lambda = eta[1:length(part)]
                nu = eta[length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0, ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.tn.dist <- function(eta, t, delta, 
                x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                nu = eta[ncol(x) + length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.tn.dist.Q1.0 <- function(eta, z, t, 
                delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.tn.dist.Q1 <- function(eta, z, t, 
                delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            m = length(t)
            nu.last = 2
            theta.last = vari.tn(nu.last)
            z.last = rep(1, m)
            z2.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
	    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  b.aux = fail.cluster(Lambda0.last, cluster) - 
                    nu.last * (nu.last + ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  aux.1 = optim(log(lambda.last), observed.llike.tn.dist.Q1.0, 
                    method = "BFGS", z = z.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part = part)
                  lambda.new = exp(aux.1$par)
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = -40, 
                    upper = 40)$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(c(lambda.last, theta.last) - 
                    c(lambda.new, theta.new)))
                  theta.last = theta.new
                  lambda.last = lambda.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.tn.dist, x0 = c(lambda.last, 
                  nu.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.tn.dist(c(lambda.last, 
                  nu.last), t = t, delta = delta, ind = cluster, 
                  dist = dist, part = part)
                se = abs(c(lambda.last, theta.last))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                llike.obs = lambda.last
                para = c(lambda.new, theta.new)
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
                  b.aux = fail.cluster(Lambda0.last * exp(x %*% 
                    beta.last), cluster) - nu.last * (nu.last + 
                    ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.tn.dist.Q1, method = "BFGS", 
                    z = z.new, t = t, delta = delta, x = x, ind = ind, 
                    dist = dist, part = part)
                  beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = -40, 
                    upper = 40)$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(c(beta.last, lambda.last, theta.last) - 
                    c(beta.new, lambda.new, theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  lambda.last = lambda.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.tn.dist, x0 = c(beta.last, 
                  lambda.last, nu.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                se = sqrt(diag(solve(aux.se)))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                llike.obs = -observed.llike.tn.dist(c(beta.last, 
                  lambda.last, nu.last), t = t, x = x, delta = delta, 
                  ind = cluster, dist = dist, part = part)
                para = c(beta.last, lambda.last, theta.last)
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = z.new, z2 = z2.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "TN"
            object.out$tau <- tau.TN(theta.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0.tn <- function(eta, t, delta, ind, 
                cox.aux) {
                nu = eta
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0, ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                -sum(P1)
            }
            observed.llike.tn <- function(eta, t, delta, x, ind, 
                cox.aux) {
                nu = eta[length(eta)]
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
                r = fail.cluster(delta, ind)
                b = fail.cluster(Lambda0 * exp(x %*% beta), ind)
                P1 = log((-1)^r * mapply(L.deriv, s = b, nu = nu, 
                  n = r))
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            m = length(t)
            nu.last = 2
            theta.last = vari.tn(nu.last)
            z.last = rep(1, m)
            z2.last = rep(1, m)
            r = fail.cluster(delta, cluster)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(z.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  b.aux = fail.cluster(Lambda0.new, cluster) - 
                    nu.last * (nu.last + ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = -40, 
                    upper = 40)$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(theta.last - theta.new))
                  theta.last = theta.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.tn, x0 = c(nu.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                para = c(theta.new)
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(z.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
                  b.aux = fail.cluster(Lambda0.new * exp(x %*% 
                    beta.last), cluster) - nu.last * (nu.last + 
                    ratio(nu.last))
                  C0 = mapply(C.k, k = 0, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C1 = mapply(C.k, k = 1, ri = r, bi = b.aux, 
                    nu = nu.last)
                  C2 = mapply(C.k, k = 2, ri = r, bi = b.aux, 
                    nu = nu.last)
                  z.new = C1/C0
                  z2.new = C2/C0
                  nu.new = optim(nu.last, prof.TN, z = z.new, 
                    z2 = z2.new, method = "Brent", lower = -40, 
                    upper = 40)$par
                  theta.new = vari.tn(nu.new)
                  dif = max(abs(c(beta.last, theta.last) - c(beta.new, 
                    theta.new)))
                  beta.last = beta.new
                  theta.last = theta.new
                  nu.last = nu.new
                  z.last = z.new
                  z2.last = z2.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.tn, x0 = c(beta.last, 
                  nu.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
                se[length(se)] = se[length(se)] * abs(dg(nu.last))
                para = c(beta.new, theta.new)
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = z.new, z2 = z2.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "TN"
            object.out$tau <- tau.TN(theta.last)
        }
        object.out
    }
frailtyMIG <- function(formula, data, dist = "np", prec = 1e-04, 
        max.iter = 1000, part = NULL) {
	Var_U <- function(m) {m^3*(1-m)+(1-m^2)^3/m^2}
        prof.MIG <- function(m, u1, y, a.0, a.1, b, p, nu.aux) {            
		kappa = (1-nu.aux)*m^(-2)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*(m+1)^(-2)*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
		kappa1 = (m+1)^(-1)*nu.aux+m^(-1)*(1-nu.aux)
            ll = log(m)-0.5*log1p(-m)-0.5*(m^2/(1-m))*(u1+kappa-2*kappa1)+
			y*log1p(-m)+(1-y)*log(m)
            -sum(ll)
        }
        fail.cluster <- function(delta, indice) {
            sum.fail <- function(ind, delta) {
                sum(delta[which(indice == ind)])
            }
            unlist(lapply(1:max(indice), sum.fail, delta = delta))
        }
        tau.MIG <- function(m) {
            aux.int <- function(x, m) {
		a1<-(m/(1-m))*(1-sqrt(1+2*(1-m)*x))
		a2<-(m^2/(1-m^2))*(1-sqrt(1+2*(1+m)^2*(1-m)*x/m^2))
                L.z <- m*exp(a1)+(1-m)*exp(a2)
			p1<-m^2*(exp(a1)*(1-m)*(1+2*(1-m)*x)^(-3/2)+(1+2*(1-m)*x)^(-1)*m*exp(a1))
			p2<-(1-m^2)*(exp(a2)*(1+m)^2*(1-m)*m^(-2)*(1+2*(1+m)^2*(1-m)*x/m^2)^(-3/2)+(1+2*(1+m)^2*(1-m)*x/m^2)^(-1)*exp(a2)*(1+m))
                L2.z <- p1+p2
                x * L.z * L2.z
            }
            4 * integrate(aux.int, lower = 0, upper = Inf, m=m)$value - 1
        }
        if (dist == "weibull") {
            observed.llike.0.mig.dist <- function(eta, t, delta, 
                ind, part = NULL) {
                rho = eta[1]
                lambda = eta[2]
                m = eta[3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0.last, cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.mig.dist <- function(eta, t, delta, 
                x, ind, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = eta[ncol(x) + 1]
                lambda = eta[ncol(x) + 2]
                m = eta[ncol(x) + 3]
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0.last*exp(x %*% beta), cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.mig.dist.Q1.0 <- function(eta, 
                z, t, delta, x, ind, dist, part = NULL) {
                rho = exp(eta[1])
                lambda = exp(eta[2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.mig.dist.Q1 <- function(eta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                rho = exp(eta[ncol(x) + 1])
                lambda = exp(eta[ncol(x) + 2])
                Lambda0 = H.base(t, lambda, rho, dist)
                log.lambda0 = h.base(t, lambda, rho, dist)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            mm = length(t)
            m.last = 0.75
            u.last = rep(1, mm)
            u1.last = rep(1, mm)
            y.last = rep(0.5, mm)
            dif = 10
            i = 1
            rho.last = 1
            lambda.last = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.last, cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                    observed.llike.mig.dist.Q1.0, method = "BFGS", 
                    z = u.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  rho.new = exp(aux.1$par[1])
                  lambda.new = exp(aux.1$par[2])
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(c(rho.last, lambda.last, m.last) - 
                    c(rho.new, lambda.new, m.new)))
                  m.last = m.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.mig.dist, x0 = c(rho.last, 
                  lambda.last, m.last), t = t, delta = delta, 
                  ind = cluster)
                se = sqrt(diag(solve(aux.se)))
                llike.obs = -observed.llike.0.mig.dist(c(rho.last, 
                  lambda.last, m.last), t = t, delta = delta, 
                  ind = cluster)
                para = c(rho.new, lambda.new, Var_U(m.new))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                names(para) = names(se) = c("rho", "lambda", 
                  "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull", data=data)
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
                rho.last = 1/cox.aux$scale
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda.last, rho.last, 
                    dist)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.last*exp(x %*% beta.last), cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                    observed.llike.mig.dist.Q1, method = "BFGS",
                    z = u.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist)
                  beta.new = aux.1$par[1:ncol(x)]
                  rho.new = exp(aux.1$par[ncol(x) + 1])
                  lambda.new = exp(aux.1$par[ncol(x) + 2])
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(c(beta.last, rho.last, lambda.last, 
                    m.last) - c(beta.new, rho.new, lambda.new, 
                    m.new)))
                  beta.last = beta.new
                  m.last = m.new
                  rho.last = rho.new
                  lambda.last = lambda.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.mig.dist, x0 = c(beta.last, 
                  rho.last, lambda.last, m.last), t = t, 
                  x = x, delta = delta, ind = cluster)
                se = sqrt(diag(solve(aux.se)))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                llike.obs = -observed.llike.mig.dist(c(beta.last, 
                  rho.last, lambda.last, m.last), t = t, 
                  x = x, delta = delta, ind = cluster)
                para = c(beta.new, rho.new, lambda.new, Var_U(m.new))
                names(para) = names(se) = c(colnames(x), "rho", 
                  "lambda", "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = u.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "MIG"
            object.out$tau <- tau.MIG(m.last)
            object.out$logLik <- llike.obs
        }
        if (dist == "pe" | dist == "exponential") {
            dist.aux = "pe"
            if (dist == "exponential") {
                dist = "pe"
                part = 0
                dist.aux = "exponential"
            }
            observed.llike.0.mig.dist <- function(eta, t, delta, 
                ind, part = NULL) {
                lambda = eta[1:length(part)]
                m = eta[length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0.last, cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                P2 = delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.mig.dist <- function(eta, t, delta, 
                x, ind, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = eta[ncol(x) + 1:length(part)]
                m = eta[ncol(x) + length(part) + 1]
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0.last*exp(x %*% beta), cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                P2 = delta * x %*% beta + delta * log.lambda0
                -(sum(P1) + sum(P2))
            }
            observed.llike.mig.dist.Q1.0 <- function(eta, 
                z, t, delta, x, ind, dist, part = NULL) {
                lambda = exp(eta[1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * log.lambda0 - z[ind] * Lambda0
                -sum(P1)
            }
            observed.llike.mig.dist.Q1 <- function(eta, 
                z, t, delta, x, ind, dist, part = NULL) {
                beta = eta[1:ncol(x)]
                lambda = exp(eta[ncol(x) + 1:length(part)])
                Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                  part = part)
                log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                  part = part)
                P1 = delta * x %*% beta + delta * log.lambda0 - 
                  z[ind] * Lambda0 * exp(x %*% beta)
                -sum(P1)
            }
            mm = length(t)
            m.last = 0.75
            u.last = rep(1, mm)
            u1.last = rep(1, mm)
            y.last = rep(0.5, mm)
            dif = 10
            i = 1
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.last, cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  aux.1 = optim(c(log(lambda.last)), 
                    observed.llike.mig.dist.Q1.0, method = "BFGS", 
                    z = u.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part=part)
			lambda.new = exp(aux.1$par)
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(c(lambda.last, m.last) - 
                    c(lambda.new, m.new)))
                  m.last = m.new
                  lambda.last = lambda.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.mig.dist, x0 = c(lambda.last, 
                  m.last), t = t, delta = delta, ind = cluster, 
                  part = part)
                se = sqrt(diag(solve(aux.se)))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                 llike.obs = -observed.llike.0.mig.dist(c(lambda.last, 
                  m.last), t = t, delta = delta, ind = cluster, 
                  part = part)
                para = c(lambda.new, Var_U(m.new))
                names(para) = names(se) = c(paste("lambda", 1:length(part), 
                  sep = ""), "theta")
               if (dist.aux == "exponential") 
                  names(para) = names(se) = c("lambda", "theta")
            }
            if (ncol(x) > 0) {
                cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
                beta.last = -coef(cox.aux)[-1]/cox.aux$scale
                lambda.last = rep(1/mean(t[which(delta==1)]), 
                  length(part))
                while (i <= max.iter & dif > prec) {
                  Lambda0.last = H.base(t, lambda = lambda.last, 
                    dist = dist, part = part)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.last*exp(x %*% beta.last), cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  aux.1 = optim(c(beta.last, log(lambda.last)), 
                    observed.llike.mig.dist.Q1, method = "BFGS",
                    z = u.new, t = t, delta = delta, 
                    x = x, ind = ind, dist = dist, part=part)
			beta.new = aux.1$par[1:ncol(x)]
                  lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(c(beta.last, lambda.last, 
                    m.last) - c(beta.new, lambda.new, 
                    m.new)))
                  beta.last = beta.new
                  m.last = m.new
                  lambda.last = lambda.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.mig.dist, x0 = c(beta.last, 
                  lambda.last, m.last), t = t, x = x, delta = delta, 
                  ind = cluster, part = part)
                se = sqrt(diag(solve(aux.se)))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                llike.obs = -observed.llike.mig.dist(c(beta.last, 
                  lambda.last, m.last), t = t, x = x, delta = delta, 
                  ind = cluster, part = part)
                para = c(beta.last, lambda.last, Var_U(m.last))
                names(para) = names(se) = c(colnames(x), paste("lambda", 
                  1:length(part), sep = ""), "theta")
                if (dist.aux == "exponential") 
                  names(para) = names(se) = c(colnames(x), "lambda", 
                    "theta")
            }
            object.out <- list(coefficients = para, se = se, 
                z = u.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$x <- x
            object.out$dist <- dist.aux
            object.out$dist.frail <- "MIG"
            object.out$tau <- tau.MIG(m.last)
            object.out$logLik <- llike.obs
            if (dist.aux == "pe") 
                object.out$part <- part
        }
        if (dist == "np") {
            observed.llike.0.mig <- function(eta, t, delta, ind, 
                cox.aux) {
                m = eta
                Lambda0 = cumhazard.basal(t, cox.aux)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0, cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                -sum(P1)
            }
            observed.llike.mig <- function(eta, t, delta, x, ind, 
                cox.aux) {
                m = eta[length(eta)]
                beta = eta[-length(eta)]
                Lambda0 = cumhazard.basal(t, cox.aux)
		    r = fail.cluster(delta, ind)
		    A = fail.cluster(Lambda0*exp(x %*% beta), cluster)
		    a.0=2*A+1/(1-m)
		    a.1=2*A+m^2/((m+1)^2*(1-m))
		    b=m^2/(1-m)
		    p=r-1/2
		    aux1=m*exp(m/(1-m))*2*besselK(sqrt(a.0*b),nu=p)*(a.0/b)^(-p/2)
		    aux2=(1-m)*exp(m^2/(1-m^2))*2*besselK(sqrt(a.1*b),nu=p)*(a.1/b)^(-p/2)
		    P1 = log(m)-0.5*(log(2*pi)+log1p(-m))+log(aux1+aux2)
                P2 = delta * x %*% beta
                -(sum(P1) + sum(P2))
            }
            cumhazard.basal <- function(t, coxph.object) {
                ind.min <- function(t0, time) {
                  min(which(time >= t0))
                }
                bb = basehaz(coxph.object)
                tt = bb$time
                bb$hazard[unlist(lapply(t, ind.min, time = tt))]
            }
            mm = length(t)
            m.last = 0.75
            u.last = rep(1, mm)
            u1.last = rep(1, mm)
            y.last = rep(0.5, mm)
            dif = 10
            i = 1
            if (ncol(x) == 0) {
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ offset(log(u.last[cluster])))
                  Lambda0.new = cumhazard.basal(t, cox.aux)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.new, cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(m.last - m.new))
                  m.last = m.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.0.mig, x0 = c(m.last), 
                  t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                para = c(Var_U(m.new))
                names(para) = names(se) = c("theta")
            }
            if (ncol(x) > 0) {
                cox.aux = coxph(Surv(t, delta) ~ x)
                beta.last = coef(cox.aux)
                while (i <= max.iter & dif > prec) {
                  cox.aux = coxph(Surv(t, delta) ~ x + offset(log(u.last[cluster])))
                  beta.new = coef(cox.aux)
                  Lambda0.new = cumhazard.basal(t, cox.aux)
			r = fail.cluster(delta, ind)
			A = fail.cluster(Lambda0.new*exp(x %*% beta.last), cluster)
			a.0=2*A+1/(1-m.last)
			a.1=2*A+m.last^2/((m.last+1)^2*(1-m.last))
			b=m.last^2/(1-m.last)
			p=r-1/2
                  xi0.aux = m.last*besselK(sqrt(a.0*b), nu=p)*a.0^(-p/2)*exp(m.last/(1-m.last))
                  xi1.aux = (1-m.last)*besselK(sqrt(a.1*b), nu=p)*a.1^(-p/2)*exp(m.last^2/((m.last+1)*(1-m.last)))
                  nu.aux = xi1.aux/(xi0.aux + xi1.aux)
                  y.new = nu.aux
                  u.new = (1-nu.aux)*sqrt(b)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(a.0)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(b)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(a.1)*besselK(sqrt(a.1*b), nu=p))
                  u1.new = (1-nu.aux)*sqrt(a.0)*besselK(sqrt(a.0*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.0*b), nu=p))+nu.aux*sqrt(a.1)*besselK(sqrt(a.1*b), nu=p+1)/(sqrt(b)*besselK(sqrt(a.1*b), nu=p))-2*p/b
                  m.new = optimize(prof.MIG, c(0.0001,0.9999), 
			u1 = u1.new, y = y.new, a.0=a.0, a.1=a.1, b=b, p=p, nu.aux=nu.aux)$minimum
                  dif = max(abs(c(beta.last, m.last) - c(beta.new, m.new)))
                  beta.last = beta.new
                  m.last = m.new
                  y.last = y.new
			u.last = u.new
			u1.last = u1.new
                  i = i + 1
                }
                aux.se = hessian(observed.llike.mig, x0 = c(beta.last, 
                  m.last), t = t, delta = delta, x = x, ind = cluster, 
                  cox.aux = cox.aux)
                se = sqrt(diag(solve(aux.se)))
			se[length(se)] = se[length(se)]*abs((8*m.new^6-3*m.new^5-6*m.new^4+2)/m.new^3)
                para = c(beta.new, Var_U(m.new))
                names(para) = names(se) = c(colnames(x), "theta")
            }
            bb = basehaz(cox.aux)
            Lambda0 = cbind(bb$time, bb$hazard)
            colnames(Lambda0) = c("time", "hazard")
            object.out <- list(coefficients = para, se = se, 
                z = u.new)
            class(object.out) <- "extrafrail"
            object.out$t <- t
            object.out$delta <- delta
            object.out$id <- cluster
            object.out$Lambda0 <- Lambda0
            object.out$x <- x
            object.out$dist <- dist
            object.out$dist.frail <- "MIG"
            object.out$tau <- tau.MIG(m.last)
        }
        object.out
    }

frailtyMBS <- function(formula, data, dist = "np", prec = 1e-04, 
                       max.iter = 1000, part = NULL) {
  Var_U <- function(m) {(2*m^2/(1-m)+5)*(m^2/(1-m)+1)^(-2)*(m^4+(1+m)^2*(1-m)^2)}
  d.Var_U<-function(m){(-12*m^8+48*m^7-100*m^6+132*m^5-90*m^4-28*m^3+84*m^2-36*m)/(m^2-m+1)^3}
  prof.MBS <- function(m, u, y, v, a.0, a.1, nu.0, nu.1, b.0, b.1, p.0, p.1) {            
    gamma=m^2/(2*m^2+1-m)
    ##kappa es E((V_i-1/2)*log(m+Y_i))
    kappa=E.g(m, type=4, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
    ##kappa1 es E(Ui/(m+Yi))
    kappa1=E.g(m, type=2, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
    ##kappa2 es E((m+Yi)/Ui)
    kappa2=E.g(m, type=3, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
    ll = y*log1p(-m)+(1-y)*log(m)+kappa+(1/2-v)*log1p((1-m)/m^2)-log(besselK(0.5*m^2/(1-m), nu=1/2))-(m^2/(1-m)+1)*kappa1/4-m^4*(1-m)^(-1)*(m^2+1-m)^(-1)*kappa2/4
    -sum(ll)
  }
  a<-function(y, a.0, a.1){ifelse(y==0, a.0, a.1)}
  b<-function(y, b.0, b.1){ifelse(y==0, b.0, b.1)}
  p<-function(y, p.0, p.1){ifelse(y==0, p.0, p.1)}
  nu<-function(y, nu.0, nu.1){ifelse(y==0, nu.0, nu.1)}
  g<-function(y, v, m, a.0, a.1, b.0, b.1, p.0, p.1, type=1)
  {
    ##type=1 for Esp(U|...), type=2 for E(U/(m+Y)|...), type=3 for E((m+Y)/U|...) and type=4 for E((V-1/2)*log(m+Y)|...)
    if(type==1){gg=sqrt(b(y,b.0, b.1)/a(y,a.0, a.1))*besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1)+1)/besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1))}
    if(type==2){gg=(m+y)^(-1)*sqrt(b(y,b.0, b.1)/a(y,a.0, a.1))*besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1)+1)/besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1))}
    if(type==3){gg=(m+y)*sqrt(a(y,a.0, a.1)/b(y,b.0, b.1))*besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1)+1)/besselK(sqrt(a(y,a.0, a.1)*b(y,b.0, b.1)),nu=p(v,p.0, p.1))-2*(m+y)*p(v,p.0, p.1)/b(y,b.0, b.1)}
    if(type==4){gg=(v-1/2)*log(m+y)}
    gg
  }
  E.g<-function(m, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
  {
    zeros=rep(0, length(a.0))
    ones=rep(1, length(a.0))
    g(zeros,zeros,m,a.0, a.1, b.0, b.1, p.0, p.1,type)*(1-gamma)+g(zeros,ones,m,a.0, a.1, b.0, b.1, p.0, p.1,type)*gamma+(g(ones,zeros,m,a.0, a.1, b.0, b.1, p.0, p.1,type)-g(zeros,zeros,m,a.0, a.1, b.0, b.1, p.0, p.1,type))*nu(zeros,nu.0,nu.1)*(1-gamma)+(g(ones,ones,m,a.0, a.1, b.0, b.1, p.0, p.1,type)-g(zeros,ones,m,a.0, a.1, b.0, b.1, p.0, p.1,type))*nu(ones,nu.0,nu.1)*gamma
  }
  fail.cluster <- function(delta, indice) {
    sum.fail <- function(ind, delta) {
      sum(delta[which(indice == ind)])
    }
    unlist(lapply(1:max(indice), sum.fail, delta = delta))
  }
  ###
  ##Hay que modificar este tau
  ###
  tau.MBS <- function(m) {
    
    
    L_MBS_deriv <- function(t, m, n) {
      
      a <- (m^2 / (1 - m) + 1) / (2 * m)
      b <- (m * (m^2 / (1 - m))^2) / (2 * (m^2 / (1 - m) + 1))
      c <- (m^2 / (1 - m) + 1) / (2 * (m + 1))
      d <- ((m + 1) * (m^2 / (1 - m))^2) / (2 * (m^2 / (1 - m) + 1))
      
      # Define the modified Bessel functions
      K <- function(nu, x) {
        return(besselK(x, nu))
      }
      sqrt_term_m <- sqrt(m^2 / (1 - m) + 1)
      sqrt_term_m1 <- sqrt((m^2 / (1 - m)) + 1 + 4 * m * t)
      sqrt_term_m1_m1 <- sqrt((m^2 / (1 - m)) + 1 + 4 * (m + 1) * t)
      L_m <- exp((m^2 / (1 - m) * (sqrt_term_m - sqrt_term_m1)) / (2 * sqrt_term_m))
      L_m1 <- exp((m^2 / (1 - m) * (sqrt_term_m - sqrt_term_m1_m1)) / (2 * sqrt_term_m))
      term1 <- (1 / 2) * m * L_m * (1 + sqrt_term_m / sqrt_term_m1)
      term2 <- (1 / 2) * (1 - m) * L_m1 * (1 + sqrt_term_m / sqrt_term_m1_m1)
      L0 <- term1 + term2
      
      term1_deriv <- -(m / 2) *(( (a / b)^(1/4) / (sqrt((a + 2 * t) / b))^(3/2) * K(3/2, sqrt(b * (a + 2 * t))) / K(1/2, sqrt(a * b))) +
                                  ((a / b)^(-1/4) / (sqrt((a + 2 * t) / b))^(1/2) * K(1/2, sqrt(b * (a + 2 * t))) / K(-1/2, sqrt(a * b))))
      term2_deriv <- -(1 - m) / 2 * (( (c / d)^(1/4) / (sqrt((c + 2 * t) / d))^(3/2) * K(3/2, sqrt(d * (c + 2 * t))) / K(1/2, sqrt(c * d))) +
                                       ((c / d)^(-1/4) / (sqrt((c + 2 * t) / d))^(1/2) * K(1/2, sqrt(d * (c + 2 * t))) / K(-1/2, sqrt(c * d))))
      L1 <- term1_deriv + term2_deriv 
      
      
      term1_deriv2 <-(m / 2)*(( (a / b)^(1/4) / (sqrt((a + 2 * t) / b))^(5/2) * K(5/2, sqrt(b * (a + 2 * t))) / K(1/2, sqrt(a * b))) +
                                ((a / b)^(-1/4) / (sqrt((a + 2 * t) / b))^(3/2) * K(3/2, sqrt(b * (a + 2 * t))) / K(-1/2, sqrt(a * b))))
      term2_deriv2 <- (1 - m) / 2 * (( (c / d)^(1/4) / (sqrt((c + 2 * t) / d))^(5/2) * K(5/2, sqrt(d * (c + 2 * t))) / K(1/2, sqrt(c * d))) +
                                       ((c / d)^(-1/4) / (sqrt((c + 2 * t) / d))^(3/2) * K(3/2, sqrt(d * (c + 2 * t))) / K(-1/2, sqrt(c * d))))
      L2 <- term1_deriv2 + term2_deriv2
      if (n == 0) L = L0
      if (n == 1) L = L1
      if (n == 2) L = L2
      L
    }
    aux.int <- function(x, m) {
      x*L_MBS_deriv(x, m, 0)*L_MBS_deriv(x, m, 2)
    }
    4 * integrate(aux.int, lower = 0, upper = Inf, m=m)$value - 1
  }
  
  if (dist == "weibull") {
    observed.llike.0.mbs.dist <- function(eta, t, delta, 
                                          ind, part = NULL) {
      rho = eta[1]
      lambda = eta[2]
      m = eta[3]
      Lambda0 = H.base(t, lambda, rho, dist)
      log.lambda0 = h.base(t, lambda, rho, dist)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0.last, cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      P2 = delta * log.lambda0
      -(sum(P1) + sum(P2))
    }
    observed.llike.mbs.dist <- function(eta, t, delta, 
                                        x, ind, part = NULL) {
      beta = eta[1:ncol(x)]
      rho = eta[ncol(x) + 1]
      lambda = eta[ncol(x) + 2]
      m = eta[ncol(x) + 3]
      Lambda0 = H.base(t, lambda, rho, dist)
      log.lambda0 = h.base(t, lambda, rho, dist)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0.last*exp(x %*% beta), cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      P2 = delta * x %*% beta + delta * log.lambda0
      -(sum(P1) + sum(P2))
    }
    observed.llike.mbs.dist.Q1.0 <- function(eta, 
                                             z, t, delta, x, ind, dist, part = NULL) {
      rho = exp(eta[1])
      lambda = exp(eta[2])
      Lambda0 = H.base(t, lambda, rho, dist)
      log.lambda0 = h.base(t, lambda, rho, dist)
      P1 = delta * log.lambda0 - z[ind] * Lambda0
      -sum(P1)
    }
    observed.llike.mbs.dist.Q1 <- function(eta, 
                                           z, t, delta, x, ind, dist, part = NULL) {
      beta = eta[1:ncol(x)]
      rho = exp(eta[ncol(x) + 1])
      lambda = exp(eta[ncol(x) + 2])
      Lambda0 = H.base(t, lambda, rho, dist)
      log.lambda0 = h.base(t, lambda, rho, dist)
      P1 = delta * x %*% beta + delta * log.lambda0 - 
        z[ind] * Lambda0 * exp(x %*% beta)
      -sum(P1)
    }
    mm = length(t)
    m.last = 0.5
    u.last = rep(1, mm)
    v.last = rep(0.5, mm)
    y.last = rep(0.5, mm)
    dif = 10
    i = 1
    rho.last = 1
    lambda.last = 1
    if (ncol(x) == 0) {
      while (i <= max.iter & dif > prec) {
        Lambda0.last = H.base(t, lambda.last, rho.last, 
                              dist)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.last, cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(0-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(0+1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(0-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        aux.1 = optim(c(log(rho.last), log(lambda.last)), 
                      observed.llike.mbs.dist.Q1.0, method = "BFGS", 
                      z = u.new, t = t, delta = delta, 
                      x = x, ind = ind, dist = dist)
        rho.new = exp(aux.1$par[1])
        lambda.new = exp(aux.1$par[2])
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(c(rho.last, lambda.last, m.last) - 
                        c(rho.new, lambda.new, m.new)))
        m.last = m.new
        rho.last = rho.new
        lambda.last = lambda.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      m.last=m.new=max(m.new, 0.001)
      aux.se = hessian(observed.llike.0.mbs.dist, x0 = c(rho.last, 
                                                         lambda.last, m.last), t = t, delta = delta, 
                       ind = cluster)
      se = sqrt(diag(solve(aux.se)))
      llike.obs = -observed.llike.0.mbs.dist(c(rho.last, 
                                               lambda.last, m.last), t = t, delta = delta, 
                                             ind = cluster)
      para = c(rho.new, lambda.new, Var_U(m.new))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      names(para) = names(se) = c("rho", "lambda", 
                                  "theta")
    }
    if (ncol(x) > 0) {
      cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull", data=data)
      beta.last = -coef(cox.aux)[-1]/cox.aux$scale
      lambda.last = exp(-coef(cox.aux)[1]/cox.aux$scale)
      rho.last = 1/cox.aux$scale
      while (i <= max.iter & dif > prec) {
        Lambda0.last = H.base(t, lambda.last, rho.last, 
                              dist)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.last*exp(x %*% beta.last), cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        aux.1 = optim(c(beta.last, log(rho.last), log(lambda.last)), 
                      observed.llike.mbs.dist.Q1, method = "BFGS",
                      z = u.new, t = t, delta = delta, 
                      x = x, ind = ind, dist = dist)
        beta.new = aux.1$par[1:ncol(x)]
        rho.new = exp(aux.1$par[ncol(x) + 1])
        lambda.new = exp(aux.1$par[ncol(x) + 2])
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(c(beta.last, rho.last, lambda.last, 
                        m.last) - c(beta.new, rho.new, lambda.new, 
                                    m.new)))
        beta.last = beta.new
        m.last = m.new
        rho.last = rho.new
        lambda.last = lambda.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      aux.se = hessian(observed.llike.mbs.dist, x0 = c(beta.last, 
                                                       rho.last, lambda.last, m.last), t = t, 
                       x = x, delta = delta, ind = cluster)
      se = sqrt(diag(solve(aux.se)))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      llike.obs = -observed.llike.mbs.dist(c(beta.last, 
                                             rho.last, lambda.last, m.last), t = t, 
                                           x = x, delta = delta, ind = cluster)
      para = c(beta.new, rho.new, lambda.new, Var_U(m.new))
      names(para) = names(se) = c(colnames(x), "rho", 
                                  "lambda", "theta")
    }
    m.last=m.new=max(m.new, 0.001)
    object.out <- list(coefficients = para, se = se, 
                       z = u.new)
    class(object.out) <- "extrafrail"
    object.out$t <- t
    object.out$delta <- delta
    object.out$id <- cluster
    object.out$x <- x
    object.out$dist <- dist
    object.out$dist.frail <- "MBS"
    object.out$tau <- tau.MBS(m.last)
    object.out$logLik <- llike.obs
  }
  if (dist == "pe" | dist == "exponential") {
    dist.aux = "pe"
    if (dist == "exponential") {
      dist = "pe"
      part = 0
      dist.aux = "exponential"
    }
    observed.llike.0.mbs.dist <- function(eta, t, delta, 
                                          ind, part = NULL) {
      lambda = eta[1:length(part)]
      m = eta[length(part) + 1]
      Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                       part = part)
      log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                           part = part)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0.last, cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      P2 = delta * log.lambda0
      -(sum(P1) + sum(P2))
    }
    observed.llike.mbs.dist <- function(eta, t, delta, 
                                        x, ind, part = NULL) {
      beta = eta[1:ncol(x)]
      lambda = eta[ncol(x) + 1:length(part)]
      m = eta[ncol(x) + length(part) + 1]
     Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                       part = part)
      log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                           part = part)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0.last*exp(x %*% beta), cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      P2 = delta * x %*% beta + delta * log.lambda0
      -(sum(P1) + sum(P2))
    }
    observed.llike.mbs.dist.Q1.0 <- function(eta, 
                                             z, t, delta, x, ind, dist, part = NULL) {
      lambda = exp(eta[1:length(part)])
      Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                       part = part)
      log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                           part = part)
      P1 = delta * log.lambda0 - z[ind] * Lambda0
      -sum(P1)
    }
    observed.llike.mbs.dist.Q1 <- function(eta, 
                                           z, t, delta, x, ind, dist, part = NULL) {
      beta = eta[1:ncol(x)]
      lambda = exp(eta[ncol(x) + 1:length(part)])
      Lambda0 = H.base(t, lambda = lambda, dist = dist, 
                       part = part)
      log.lambda0 = h.base(t, lambda = lambda, dist = dist, 
                           part = part)
      P1 = delta * x %*% beta + delta * log.lambda0 - 
        z[ind] * Lambda0 * exp(x %*% beta)
      -sum(P1)
    }
    mm = length(t)
    m.last = 0.5
    u.last = rep(1, mm)
    v.last = rep(0.5, mm)
    y.last = rep(0.5, mm)
    dif = 10
    i = 1
    lambda.last = rep(1/mean(t[which(delta==1)]), 
                      length(part))
    if (ncol(x) == 0) {
      while (i <= max.iter & dif > prec) {
        Lambda0.last = H.base(t, lambda = lambda.last, 
                              dist = dist, part = part)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.last, cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        aux.1 = optim(c(log(lambda.last)), 
                      observed.llike.mbs.dist.Q1.0, method = "BFGS", 
                      z = u.new, t = t, delta = delta, 
                      x = x, ind = ind, dist = dist, part=part)
        lambda.new = exp(aux.1$par)
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(c(lambda.last, m.last) - 
                        c(lambda.new, m.new)))
        m.last = m.new
        lambda.last = lambda.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      m.last=m.new=max(m.new, 0.001)
      aux.se = hessian(observed.llike.0.mbs.dist, x0 = c(lambda.last, 
                                                         m.last), t = t, delta = delta, ind = cluster, 
                       part = part)
      se = sqrt(diag(solve(aux.se)))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      llike.obs = -observed.llike.0.mbs.dist(c(lambda.last, 
                                               m.last), t = t, delta = delta, ind = cluster, 
                                             part = part)
      para = c(lambda.new, Var_U(m.new))
      names(para) = names(se) = c(paste("lambda", 1:length(part), 
                                        sep = ""), "theta")
      if (dist.aux == "exponential") 
        names(para) = names(se) = c("lambda", "theta")
    }
    if (ncol(x) > 0) {
      cox.aux = survreg(Surv(t, delta) ~ x, dist = "weibull")
      beta.last = -coef(cox.aux)[-1]/cox.aux$scale
      lambda.last = rep(sum(delta)/sum(t), 
                        length(part))
	
      while (i <= max.iter & dif > prec) {
        Lambda0.last = H.base(t, lambda = lambda.last, 
                              dist = dist, part = part)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.last*exp(x %*% beta.last), cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        aux.1 = optim(c(beta.last, log(lambda.last)), 
                      observed.llike.mbs.dist.Q1, method = "BFGS",
                      z = u.new, t = t, delta = delta, 
                      x = x, ind = ind, dist = dist, part=part)
        beta.new = aux.1$par[1:ncol(x)]
        lambda.new = exp(aux.1$par[ncol(x) + 1:length(part)])
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(c(beta.last, lambda.last, 
                        m.last) - c(beta.new, lambda.new, 
                                    m.new)))
        beta.last = beta.new
        m.last = m.new
        lambda.last = lambda.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      m.last=m.new=max(m.new, 0.001)
      aux.se = hessian(observed.llike.mbs.dist, x0 = c(beta.last, 
                                                       lambda.last, m.last), t = t, x = x, delta = delta, 
                       ind = cluster, part = part)
      se = sqrt(diag(solve(aux.se)))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      llike.obs = -observed.llike.mbs.dist(c(beta.last, 
                                             lambda.last, m.last), t = t, x = x, delta = delta, 
                                           ind = cluster, part = part)
      para = c(beta.last, lambda.last, Var_U(m.last))
      names(para) = names(se) = c(colnames(x), paste("lambda", 
                                                     1:length(part), sep = ""), "theta")
      if (dist.aux == "exponential") 
        names(para) = names(se) = c(colnames(x), "lambda", 
                                    "theta")
    }
    object.out <- list(coefficients = para, se = se, 
                       z = u.new)
    class(object.out) <- "extrafrail"
    object.out$t <- t
    object.out$delta <- delta
    object.out$id <- cluster
    object.out$x <- x
    object.out$dist <- dist.aux
    object.out$dist.frail <- "MBS"
    object.out$tau <- tau.MBS(m.last)
    object.out$logLik <- llike.obs
    if (dist.aux == "pe") 
      object.out$part <- part
  }
  if (dist == "np") {
    observed.llike.0.mbs <- function(eta, t, delta, ind, 
                                     cox.aux) {
      m = eta
      Lambda0 = cumhazard.basal(t, cox.aux)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0, cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      -sum(P1)
    }
    observed.llike.mbs <- function(eta, t, delta, x, ind, 
                                   cox.aux) {
      m = eta[length(eta)]
      beta = eta[-length(eta)]
      Lambda0 = cumhazard.basal(t, cox.aux)
      r = fail.cluster(delta, ind)
      A = fail.cluster(Lambda0*exp(x %*% beta), cluster)
      theta=m^2/(1-m)
      p=r-1/2
      a.m=2*A+(theta+1)/(2*m)
      b.m=m*theta^2/(2*(theta+1))
      a.m1=2*A+(theta+1)/(2*(m+1))
      b.m1=(m+1)*theta^2/(2*(theta+1))
      aux1=m*(sqrt(pi*m))^(-1)*(b.m/a.m)^(p/2)*(m*theta*(theta+1)^(-1)*besselK(sqrt(a.m*b.m),nu=p)+sqrt(b.m/a.m)*besselK(sqrt(a.m*b.m),nu=p+1))
      aux2=(1-m)*(sqrt(pi*(m+1)))^(-1)*(b.m1/a.m1)^(p/2)*((m+1)*theta*(theta+1)^(-1)*besselK(sqrt(a.m1*b.m1),nu=p)+sqrt(b.m1/a.m1)*besselK(sqrt(a.m1*b.m1),nu=p+1))
      P1 = (theta/2)-log(2)+(1/2)*log1p(theta)+log(aux1+aux2)
      P2 = delta * x %*% beta
      -(sum(P1) + sum(P2))
    }
    cumhazard.basal <- function(t, coxph.object) {
      ind.min <- function(t0, time) {
        min(which(time >= t0))
      }
      bb = basehaz(coxph.object)
      tt = bb$time
      bb$hazard[unlist(lapply(t, ind.min, time = tt))]
    }
    mm = length(t)
    m.last = 0.5
    u.last = rep(1, mm)
    v.last = rep(0.5, mm)
    y.last = rep(0.5, mm)
    dif = 10
    i = 1
    if (ncol(x) == 0) {
      while (i <= max.iter & dif > prec) {
        cox.aux = coxph(Surv(t, delta) ~ offset(log(u.last[cluster])))
        Lambda0.new = cumhazard.basal(t, cox.aux)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.new, cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(m.last - m.new))
        m.last = m.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      m.last=m.new=max(m.new, 0.001)
      aux.se = hessian(observed.llike.0.mbs, x0 = c(m.last), 
                       t = t, delta = delta, ind = cluster, cox.aux = cox.aux)
      se = sqrt(diag(solve(aux.se)))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      para = c(Var_U(m.new))
      names(para) = names(se) = c("theta")
    }
    if (ncol(x) > 0) {
      cox.aux = coxph(Surv(t, delta) ~ x)
      beta.last = coef(cox.aux)
      while (i <= max.iter & dif > prec) {
        cox.aux = coxph(Surv(t, delta) ~ x + offset(log(u.last[cluster])))
        beta.new = coef(cox.aux)
        Lambda0.new = cumhazard.basal(t, cox.aux)
        r = fail.cluster(delta, ind)
        A = fail.cluster(Lambda0.new*exp(x %*% beta.last), cluster)
        gamma=m.last^2/(2*m.last^2+1-m.last)
        a.0=2*A+(m.last^2/(1-m.last)+1)/(2*m.last)
        a.1=2*A+(m.last^2/(1-m.last)+1)/(2*(m.last+1))
        b.0=(m.last)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        b.1=(m.last+1)*(m.last^2/(1-m.last))^2/(2*(m.last^2/(1-m.last)+1))
        p.0=r+1/2
        p.1=r+1/2-1
        nu.0=besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.0)*(a.0/b.0)^(-p.0/2)*m.last^(1/2)+besselK(sqrt(a.1*b.1), nu=p.0)*(a.1/b.1)^(-p.0/2)*(m.last+1)^(-1/2)*(1-m.last))
        nu.1=besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last)/(besselK(sqrt(a.0*b.0), nu=p.1)*(a.0/b.0)^(-p.1/2)*m.last^(1+1/2)+besselK(sqrt(a.1*b.1), nu=p.1)*(a.1/b.1)^(-p.1/2)*(m.last+1)^(1-1/2)*(1-m.last))
        u.new=E.g(m.last, type=1, a.0, a.1, b.0, b.1, p.0, p.1, nu.0, nu.1, gamma)
        y.new=nu.0*(1-gamma)+nu.1*gamma
        v.new=gamma
        m.new = optimize(prof.MBS, c(0.0001,0.9999), 
                         u = u.new, y = y.new, v = v.new, a.0=a.0, a.1=a.1, b.0=b.0, b.1=b.1, p.0=p.0, p.1=p.1, nu.0=nu.0, nu.1=nu.1)$minimum
        dif = max(abs(c(beta.last, m.last) - c(beta.new, m.new)))
        beta.last = beta.new
        m.last = m.new
        y.last = y.new
        u.last = u.new
        v.last = v.new
        i = i + 1
      }
      m.last=m.new=max(m.new, 0.001)
      aux.se = hessian(observed.llike.mbs, x0 = c(beta.last, 
		m.last), t = t, delta = delta, x = x, ind = cluster,
		cox.aux = cox.aux)
      se = sqrt(diag(solve(aux.se)))
      se[length(se)] = se[length(se)]*abs(d.Var_U(m.new))
      para = c(beta.new, Var_U(m.new))
      names(para) = names(se) = c(colnames(x), "theta")
    }
    bb = basehaz(cox.aux)
    Lambda0 = cbind(bb$time, bb$hazard)
    colnames(Lambda0) = c("time", "hazard")
    object.out <- list(coefficients = para, se = se, 
                       z = u.new)
    class(object.out) <- "extrafrail"
    object.out$t <- t
    object.out$delta <- delta
    object.out$id <- cluster
    object.out$Lambda0 <- Lambda0
    object.out$x <- x
    object.out$dist <- dist
    object.out$dist.frail <- "MBS"
    object.out$tau <- tau.MBS(m.last)
  }
  object.out
}
    if (dist.frail == "GA") 
        val <- frailtyGA(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "IG") 
        val <- frailtyIG(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "WL") 
        val <- frailtyWL(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "BS") 
        val <- frailtyBS(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "TN") 
        val <- frailtyTN(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "MIG") 
        val <- frailtyMIG(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    if (dist.frail == "MBS") 
        val <- frailtyMBS(formula, data, dist = dist, prec = prec, 
            max.iter = max.iter, part = part)
    val
}


