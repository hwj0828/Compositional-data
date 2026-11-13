rm(list = ls())
library(abess)
library(MASS)
library(SparseDOSSA2)

gene.data <- function(n, para){
  p <- para$dim
  nfact <- para$nfact
  beta <- para$beta
  delta <- para$delta
  alpha <- para$alpha
  # u <- matrix(rnorm(n*ncol(alpha)), n, ncol(alpha))
  # h <- u%*%t(alpha) + mvrnorm(n, mu=rep(0,p), Sigma=diag(1 + (1:p)/5))
  # G <- diag(p) - matrix(1, p, p)/p
  # z <- h%*%G
  # y <- rbinom(n, 1, plogis(h%*%beta + u%*%delta))
  
  
  template <- stool.template(p)
  features <- names(template$EM_fit$fit$pi0)
  
  u <- matrix(rnorm(n*nfact), nrow=n, ncol=nfact)
  colnames(u) <- if(nfact == 1L) "u" else paste0("u", seq_len(nfact))
  spike_df <- data.frame(metadata_datum=rep(seq_len(nfact), each=p),
                         feature_spiked=rep(features, times=nfact), associated_property="abundance",
                         effect_size=as.numeric(alpha), stringsAsFactors=FALSE)
  sim <- SparseDOSSA2(template=template, n_sample=n, new_features=FALSE,
                      metadata_matrix=u, spike_metadata=spike_df, verbose=FALSE)
  # Each absolute abundance is a Poisson intensity on the same fixed count scale.
  lambda <- t(sim$simulated_matrices$a_spiked[features, , drop=FALSE])
  counts <- matrix(rpois(length(lambda), lambda), nrow=n, ncol=p, dimnames=dimnames(lambda))
  # Use one pseudocount for the observed log counts and their centered log ratios.
  h <- log1p(counts)
  z <- sweep(h, 1L, rowMeans(h), "-")
  y <- rbinom(n, 1, plogis(-10.3 + h%*%beta + u%*%delta))
  cbind(y,z,h,u)
}

# Use the same p most prevalent Stool taxa in every replicate.
stool.template <- function(p){
  stool <- utils::getFromNamespace("Stool", "SparseDOSSA2")
  pars <- stool$EM_fit$fit
  if(p > length(pars$pi0)) stop("p exceeds the number of taxa in the Stool template.")
  keep <- order(pars$pi0)[seq_len(p)]
  Sigma <- pars$Sigma[keep, keep, drop=FALSE]
  fixed <- list(pi0=pars$pi0[keep], mu=pars$mu[keep], sigma=pars$sigma[keep],
                Sigma=Sigma, Omega=solve(Sigma))
  list(EM_fit=list(fit=fixed), depth_fit=stool$depth_fit)
}


# x.rotate is x.res %*% ptilda, with dimx-1 columns.
# One call to stats::factanal estimates loadings (nfact by dimx) and uniqueness.
# Identification uses the row-deletion condition on [1, alpha].
# This covariance estimator and the finite-sample projected MLE can differ.
factor_analysis <- function(x.rotate, ptilda, nfact) {
  x.rotate <- as.matrix(x.rotate)
  ptilda <- as.matrix(ptilda)
  dimx <- ncol(x.rotate) + 1L
  reference <- dimx
  min_uniqueness <- 1e-4
  
  
  A <- cov(x.rotate) * (1 - 1 / nrow(x.rotate))
  if (nrow(x.rotate) < dimx || !isTRUE(rcond(A) > .Machine$double.eps)) {
    stop("The single fit requires a full-rank projected covariance A (n >= dimx).")
  }
  
  # x.rotate %*% t(contrast) gives x.res[, j] - x.res[, reference].
  contrast <- sweep(ptilda, 2L, ptilda[reference, ], "-")
  A_diff <- contrast %*% A %*% t(contrast)
  
  # Add the covariance of one shared independent auxiliary variable.
  tau2 <- mean(diag(A_diff)[-reference])
  A_aug <- (A_diff + t(A_diff)) / 2 + matrix(tau2, dimx, dimx)
  lower <- min(1e-6, min_uniqueness / max(diag(A_aug), tau2 + min_uniqueness))
  result <- factanal(
    covmat = A_aug, factors = nfact + 1L, rotation = "none",
    control = list(lower = lower, opt = list(maxit = 500L, factr = 1e5))
  )
  
  # Undo standardization and recover the original reference error variance.
  loadings <- t(sweep(unclass(result$loadings), 1L, sqrt(diag(A_aug)), "*"))
  uniqueness <- diag(A_aug) * unname(result$uniquenesses)
  ref_fraction <- unname(result$uniquenesses[reference])
  if (1 - ref_fraction <= sqrt(.Machine$double.eps) ||
      sum(loadings[, reference]^2) <= .Machine$double.eps * tau2 ||
      qr(t(loadings))$rank < nfact + 1L) {
    stop("The fitted loadings cannot separate the common reference factor.")
  }
  uniqueness[reference] <- tau2 * ref_fraction / (1 - ref_fraction)
  uniqueness <- pmax(uniqueness, min_uniqueness)
  
  # Remove the auxiliary factor to recover alpha_j - alpha_reference.
  rotation <- qr.Q(qr(matrix(loadings[, reference], ncol = 1L)), complete = TRUE)
  loadings <- crossprod(rotation[, -1L, drop = FALSE], loadings)
  loadings[, reference] <- 0
  
  # Enforce the two original equality constraints.
  loadings <- sweep(loadings, 1L, rowMeans(loadings), "-")
  Sigma_inv <- diag(1 / uniqueness, nrow = dimx)
  M <- loadings %*% Sigma_inv %*% t(loadings)
  loadings <- crossprod(eigen(M, symmetric = TRUE)$vectors, loadings)
  
  list(loadings = loadings, uniqueness = uniqueness)
}

# Gauss-Hermite nodes and weights for a standard normal variable
normal.quadrature <- function(k){
  J <- matrix(0, k, k)
  j <- seq_len(k-1L)
  J[cbind(j, j+1L)] <- sqrt(j)
  ev <- eigen(J + t(J), symmetric=TRUE)
  list(nodes=ev$values, weights=ev$vectors[1, ]^2)
}

# x has only the selected columns. S is Var((U,Q) | Z,v), or Var(U | H,v)
# for the comparison using observed H. With S omitted this fits the plug-in logit.
logistic_reg <- function(y,x,hatU,total,v=NULL,initial_para,dimx,S=NULL){
  x <- as.matrix(x)
  hatU <- as.matrix(hatU)
  nfact <- ncol(hatU)
  if(!is.null(total)) x <- x + as.numeric(total)/dimx
  regressor <- cbind(1, hatU, v, x)
  active <- ncol(regressor) - ncol(x) + seq_len(ncol(x))
  response.sign <- 2*as.numeric(y)-1
  
  # Map the regression coefficients to (delta, sum(beta)/dimx).
  nlatent <- nfact + as.integer(!is.null(total))
  A <- matrix(0, nlatent, ncol(regressor))
  A[cbind(seq_len(nfact), 1L+seq_len(nfact))] <- 1
  if(!is.null(total)) A[nlatent, active] <- 1/dimx
  if(is.null(S)) S <- matrix(0, nlatent, nlatent)
  K <- crossprod(A, S%*%A)
  
  loss <- function(para, rule){
    mu <- response.sign * drop(regressor%*%para)
    Kpara <- drop(K%*%para)
    sd <- sqrt(max(0, sum(para*Kpara)))
    eta <- outer(mu, sd*rule$nodes, "+")
    logp <- plogis(eta, log.p=TRUE)
    shift <- logp[, which.max(rule$nodes)]
    relative <- exp(logp - shift)
    prob <- drop(relative%*%rule$weights)
    value <- -mean(shift + log(prob))
    
    # Differentiate the same quadrature rule, including its parameter-dependent sd.
    deriv <- relative * (-expm1(logp))
    score.mu <- drop(deriv%*%rule$weights)/prob
    gradient <- -drop(crossprod(regressor, response.sign*score.mu))/length(y)
    if(sd > 0){
      score.sd <- drop(deriv%*%(rule$weights*rule$nodes))/prob
      gradient <- gradient - mean(score.sd)*Kpara/sd
    }
    attr(value, "gradient") <- gradient
    value
  }
  
  initial_para <- as.numeric(initial_para)
  if(length(initial_para) != ncol(regressor)) stop("Wrong length for initial_para.")
  # Check the likelihood and score again with twice as many nodes.
  tryCatch({
    for(k in c(24L,48L,96L)){
      fit <- nlm(loss, initial_para, rule=normal.quadrature(k), gradtol=1e-8,
                 iterlim=150, check.analyticals=FALSE)
      check <- loss(fit$estimate, normal.quadrature(2L*k))
      if(is.finite(check) && abs(as.numeric(check)-fit$minimum) < 1e-8 &&
         max(abs(attr(check,"gradient"))) < 2e-6) return(fit$estimate)
      initial_para <- fit$estimate
    }
    stop("The integrated logistic fit did not pass its likelihood and score checks.")
  }, error = function(e) {
    stop(errorCondition(conditionMessage(e), class = "integrated_fit_error",
                        call = conditionCall(e)))
  })
}

esti.compositional <- function(dta=NULL, subset=NULL, y,x,v=NULL, nfact){
  if(is.null(dta)){
    y <- as.matrix(y)
    x <- as.matrix(x)
    if(!is.null(v)) v <- as.matrix(v)
  }else{
    y <- as.matrix(dta[, y, drop=FALSE])
    x <- as.matrix(dta[, x, drop=FALSE])
    if(!is.null(v)) v <- as.matrix(dta[, v, drop=FALSE])
  }
  if(!is.null(subset)){
    y <- y[subset, , drop=FALSE]
    x <- x[subset, , drop=FALSE]
    if(!is.null(v)) v <- v[subset, , drop=FALSE]
  }
  n <- nrow(x)
  dimx <- ncol(x)
  y <- as.numeric(y)
  x.res <- if(is.null(v)) scale(x, scale=FALSE) else lm(x~v)$residuals
  
  # Factor fit using contrasts of the compositional coordinates.
  G <- diag(dimx) - matrix(1, dimx, dimx)/dimx
  Ptilda <- eigen(G, symmetric=TRUE)$vectors[, -dimx, drop=FALSE]
  x.rotate <- x.res%*%Ptilda
  result.fact <- factor_analysis(x.rotate, Ptilda, nfact)
  halpha <- result.fact$loadings
  uniqueness <- result.fact$uniqueness
  
  # Conditional means and covariance of (U,Q) given the observed contrasts.
  Hvar <- crossprod(halpha) + diag(uniqueness, nrow=dimx)
  sigma_x <- crossprod(Ptilda, Hvar%*%Ptilda)
  uv_z <- halpha%*%Ptilda
  S <- diag(nfact)
  uv_z <- rbind(uv_z, crossprod(rowSums(Hvar), Ptilda))
  cov.uq <- rowSums(halpha)
  S <- rbind(cbind(S, cov.uq), c(cov.uq, sum(Hvar)))
  slopes <- solve(sigma_x, t(uv_z))
  cond.mean <- x.rotate%*%slopes
  S <- S - uv_z%*%slopes
  S <- (S+t(S))/2
  hatU <- cond.mean[, seq_len(nfact), drop=FALSE]
  total <- cond.mean[, nfact+1L]
  x.obs <- x.rotate%*%t(Ptilda)
  hatH <- x.obs + total/dimx
  
  # Keep the binary y and all adjustment variables. Residualize only hatH.
  nuisance <- cbind(1, hatU, v)
  lm.H <- lm.fit(nuisance, hatH)
  if(lm.H$rank < ncol(nuisance)) stop("The adjustment variables are collinear.")
  k <- ncol(nuisance)-1L
  max.size <- min(floor((dimx - nfact - 1L)/2),
                  n - lm.H$rank - 1L)
  if(max.size < 0L) stop("Insufficient residual degrees of freedom.")
  regressor <- cbind(nuisance[, -1, drop=FALSE], lm.H$residuals)
  colnames(regressor) <- NULL
  abess_fit <- abess(regressor, y, family="binomial", tune.type="bic",
                     always.include=seq_len(k), support.size=k+(0:max.size))
  b <- as.numeric(coef(abess_fit, support.size=abess_fit$best.size))
  hbeta <- tail(b, dimx)
  x_index <- which(hbeta != 0)
  if(!length(x_index)) return(numeric(dimx))
  
  # Convert nuisance slopes back to the original hatH coordinates before refitting.
  initial <- head(b, k+1L) - drop(lm.H$coefficients[, x_index, drop=FALSE]%*%hbeta[x_index])
  estimate <- logistic_reg(y, x.obs[, x_index, drop=FALSE], hatU, total, v,
                           c(initial, hbeta[x_index]), dimx, S)
  hbeta[x_index] <- tail(estimate, length(x_index))
  return(as.numeric(hbeta))
}

esti.confounder <- function(dta=NULL, subset=NULL, y,x,v=NULL, nfact){
  if(is.null(dta)){
    y <- as.matrix(y)
    x <- as.matrix(x)
    if(!is.null(v)) v <- as.matrix(v)
  }else{
    y <- as.matrix(dta[, y, drop=FALSE])
    x <- as.matrix(dta[, x, drop=FALSE])
    if(!is.null(v)) v <- as.matrix(dta[, v, drop=FALSE])
  }
  if(!is.null(subset)){
    y <- y[subset, , drop=FALSE]
    x <- x[subset, , drop=FALSE]
    if(!is.null(v)) v <- v[subset, , drop=FALSE]
  }
  n <- nrow(x)
  dimx <- ncol(x)
  y <- as.numeric(y)
  x.res <- if(is.null(v)) scale(x, scale=FALSE) else lm(x~v)$residuals
  
  # Factor fit using all coordinates of the observed H.
  Ptilda <- diag(dimx)
  x.rotate <- x.res
  result.fact <- tryCatch(factanal(x.res, factors=nfact, rotation="none"),
                          error=function(e){warning(conditionMessage(e)); NULL})
  if(is.null(result.fact)) return(rep(NA_real_, dimx))
  sdx <- sqrt(diag(var(x.res)))
  halpha <- t(sweep(unclass(result.fact$loadings), 1, sdx, "*"))
  uniqueness <- sdx^2 * as.numeric(result.fact$uniquenesses)
  
  # Conditional means and covariance of U given the observed H.
  Hvar <- crossprod(halpha) + diag(uniqueness, nrow=dimx)
  sigma_x <- crossprod(Ptilda, Hvar%*%Ptilda)
  uv_z <- halpha%*%Ptilda
  S <- diag(nfact)
  slopes <- solve(sigma_x, t(uv_z))
  cond.mean <- x.rotate%*%slopes
  S <- S - uv_z%*%slopes
  S <- (S+t(S))/2
  hatU <- cond.mean[, seq_len(nfact), drop=FALSE]
  total <- NULL
  x.obs <- x.rotate%*%t(Ptilda)
  hatH <- x.obs
  
  # Keep the binary y and all adjustment variables. Residualize only hatH.
  nuisance <- cbind(1, hatU, v)
  lm.H <- lm.fit(nuisance, hatH)
  if(lm.H$rank < ncol(nuisance)) stop("The adjustment variables are collinear.")
  k <- ncol(nuisance)-1L
  max.size <- min(floor((dimx - nfact)/2), n - lm.H$rank - 1L)
  if(max.size < 0L) stop("Insufficient residual degrees of freedom.")
  regressor <- cbind(nuisance[, -1, drop=FALSE], lm.H$residuals)
  colnames(regressor) <- NULL
  abess_fit <- abess(regressor, y, family="binomial", tune.type="bic",
                     always.include=seq_len(k), support.size=k+(0:max.size))
  b <- as.numeric(coef(abess_fit, support.size=abess_fit$best.size))
  hbeta <- tail(b, dimx)
  x_index <- which(hbeta != 0)
  if(!length(x_index)) return(numeric(dimx))
  
  # Convert nuisance slopes back to the original hatH coordinates before refitting.
  initial <- head(b, k+1L) - drop(lm.H$coefficients[, x_index, drop=FALSE]%*%hbeta[x_index])
  estimate <- logistic_reg(y, x.obs[, x_index, drop=FALSE], hatU, total, v,
                           c(initial, hbeta[x_index]), dimx, S)
  hbeta[x_index] <- tail(estimate, length(x_index))
  return(as.numeric(hbeta))
}


# During simulation, a failed integrated refit marks only this method's row.
fit_or_na <- function(expr, method, i, seed, p) {
  tryCatch(expr, integrated_fit_error = function(e) {
    message(sprintf("Replicate %d (seed %d), %s: %s", i, seed, method, conditionMessage(e)))
    rep(NA_real_, p)
  })
}

n <- 500
p <- 10
nfact <- 1

set.seed(10)
alpha <- matrix(rnorm(p*nfact), nrow=p, ncol=nfact)
alpha <- sweep(alpha, 2, colMeans(alpha), "-")
para <- list(dim = p, nfact = nfact, beta = c(c(1, 2), rep(0, p-2)),
             delta = rep(1, nfact), alpha = alpha)

nsims <- 100
result1 <- matrix(NA, nrow = nsims, ncol = p)
result2 <- matrix(NA, nrow = nsims, ncol = p)
result3 <- matrix(NA, nrow = nsims, ncol = p)
for (i in 1:nsims) {
  # abess resets the RNG, so fix each simulated dataset separately.
  seed <- 20261100L + i
  set.seed(seed)
  data <- gene.data(n, para)
  result1[i,] <- fit_or_na(
    esti.compositional(dta=data, y=1, x=(2:(p+1)), nfact = nfact),
    "esti.compositional", i, seed, p)
  colnames(data)[2:(p+1)] <- paste0("V", seq_len(p))
  abess_fit <- abess(data[, 2:(p+1)], data[,1], family="binomial", tune.type = "bic")
  result2[i,]<- coef(abess_fit, support.size =
                       abess_fit$best.size)[-1]
  result3[i,] <- fit_or_na(
    esti.confounder(dta=data, y=1, x=((2+p):(1+2*p)), nfact = nfact),
    "esti.confounder", i, seed, p)
}

fit_ok <- lapply(list(compositional=result1, abess=result2, confounder=result3),
                 function(x) rowSums(!is.finite(x)) == 0L)
failed_reps <- lapply(fit_ok, function(ok) which(!ok))
print(data.frame(successful = nsims - lengths(failed_reps), failed = lengths(failed_reps)))

bias1 <- sweep(result1, 2, para$beta, "-")
n_ok <- sum(fit_ok$compositional)
mc.summary <- data.frame(truth=para$beta, mean=NA_real_, bias=NA_real_, mcse=NA_real_)
if (n_ok > 0L) {
  estimates <- result1[fit_ok$compositional, , drop=FALSE]
  mc.summary$mean <- colMeans(estimates)
  mc.summary$bias <- colMeans(bias1[fit_ok$compositional, , drop=FALSE])
  mc.summary$mcse <- apply(estimates, 2L, sd) / sqrt(n_ok)
}
print(round(mc.summary, 4))
group_names <- expression(beta[1], beta[2], beta[3], beta[4], beta[5],
                          beta[6], beta[7], beta[8], beta[9], beta[10])
plot_fits <- function(x, ok, label) {
  # heading <- sprintf("%s (%d/%d fits)", label, sum(ok), length(ok))
  heading <- sprintf("%s", label)
  if (!any(ok)) {
    plot.new()
    title(main = heading)
    text(0.5, 0.5, "No successful fits")
    return(invisible(NULL))
  }
  boxplot(x[ok, , drop=FALSE], names=group_names, outline=FALSE, main=heading)
}
plot_fits(bias1,   fit_ok$compositional, "Compositional estimation error")
plot_fits(result1, fit_ok$compositional, "Decipher estimates")
plot_fits(result2, fit_ok$abess,         "LASSO estimates")
plot_fits(result3, fit_ok$confounder,    "Confounder estimates")

# data <- gene.data(n, para)
# result1 <- esti.compositional(dta=data, y=1, x=(2:(p+1)), nfact = nfact)
# result2 <- esti.confounder(dta=data, y=1, x=((2+p):(1+2*p)), nfact = nfact)
# 
# result1
# result2

