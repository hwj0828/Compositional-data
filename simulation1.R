rm(list = ls())
library(abess)
library(MASS)
library(nloptr)
library(SparseDOSSA2)

gene.data <- function(n, para){
  p <- para$dim
  nfact <- para$nfact
  beta <- para$beta
  delta <- para$delta
  alpha <- para$alpha
  
  # u <- rnorm(n)
  # h <- u%*%t(alpha) + mvrnorm(n, mu=rep(0,p), Sigma = diag(1 + (1:p)/5))
  # G <- diag(p) - matrix(1, nrow = p, ncol = p)/p
  # z <- h%*%G
  # y <- h%*%beta + u%*%delta + rnorm(n)
  
  template <- stool.template(p)
  features <- names(template$EM_fit$fit$pi0)
  
  u <- matrix(rnorm(n*nfact), nrow=n, ncol=nfact)
  colnames(u) <- if(nfact == 1L) "u" else paste0("u", seq_len(nfact))
  spike_df <- data.frame(metadata_datum=rep(seq_len(nfact), each=p),
                         feature_spiked=rep(features, times=nfact), associated_property="abundance",
                         effect_size=as.numeric(alpha), stringsAsFactors=FALSE)
  sim <- SparseDOSSA2::SparseDOSSA2(template=template, n_sample=n, new_features=FALSE,
                                    metadata_matrix=u, spike_metadata=spike_df, verbose=FALSE)
  # Each absolute abundance is a Poisson intensity on the same fixed count scale.
  lambda <- t(sim$simulated_matrices$a_spiked[features, , drop=FALSE])
  counts <- matrix(rpois(length(lambda), lambda), nrow=n, ncol=p, dimnames=dimnames(lambda))
  # Use one pseudocount for the observed log counts and their centered log ratios.
  h <- log1p(counts)
  z <- sweep(h, 1L, rowMeans(h), "-")
  y <- h%*%beta + u%*%delta + rnorm(n)
  
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

esti.compositional <- function(dta=NULL, subset=NULL, y,x,v=NULL, nfact){
  if(is.null(dta)){
    if(is.null(subset)) subset <- seq_len(nrow(as.matrix(x)))
    y <- as.matrix(y)[subset, , drop=FALSE]
    x <- as.matrix(x)[subset, , drop=FALSE]
    if(!is.null(v)) v <- as.matrix(v)[subset, , drop=FALSE]
  }else{
    if(is.null(subset)) subset <- seq_len(nrow(dta))
    y <- as.matrix(dta[subset, y, drop=FALSE])
    x <- as.matrix(dta[subset, x, drop=FALSE])
    if(!is.null(v)) v <- as.matrix(dta[subset, v, drop=FALSE])
  }
  
  n <- nrow(x)
  dimx <- ncol(x)
  if(ncol(y) != 1 || nrow(y) != n || !is.numeric(y) || any(!is.finite(y))){
    stop("y must be a finite numeric response with one value per observation.")
  }
  
  # factor analysis
  if(is.null(v)){
    x.res <- x
  }else{
    x.res <- lm(x~v)$residuals
  }
  
  G <- diag(dimx) - matrix(1, nrow=dimx, ncol=dimx)/dimx
  P <- eigen(G)$vectors
  Ptilda <- P[, -dimx, drop=FALSE]
  x.rotate <- x.res%*%Ptilda
  
  result.fact <- factor_analysis(x.rotate, Ptilda, nfact)
  halpha <- result.fact$loadings
  uniqueness <- result.fact$uniqueness
  
  # estimate U, total and H
  sigma_h <- t(halpha)%*%halpha + diag(uniqueness, nrow=dimx)
  sigma_x <- t(Ptilda)%*%(sigma_h%*%Ptilda)
  hatU <- x.rotate%*%solve(sigma_x, t(Ptilda)%*%t(halpha))
  total <- x.rotate%*%solve(sigma_x, t(Ptilda)%*%(sigma_h%*%rep(1, dimx)))
  # the coefficient on total in hatH %*% beta is sum(beta)/dimx
  hatH <- x.rotate%*%t(Ptilda) + matrix(total/dimx, n, dimx)
  
  # residualize y and hatH on hatU and any observed v
  regressor <- cbind(hatU, v)
  lm.res <- lm(cbind(y, hatH) ~ regressor)
  y.res <- lm.res$residuals[, 1]
  hatH.res <- lm.res$residuals[, -1, drop=FALSE]
  colnames(hatH.res) <- NULL
  
  # select the nonzero coefficients of beta
  max.size <- min(floor((dimx - nfact - 1)/2), n - lm.res$rank - 1)
  if(max.size < 0){
    stop("Insufficient residual degrees of freedom for subset selection.")
  }
  abess_fit <- abess(hatH.res, y.res, family="gaussian", tune.type="bic",
                            support.size=0:max.size, fit.intercept=FALSE)
  hbeta <- coef(abess_fit, support.size=abess_fit$best.size)[-1]
  return(as.numeric(hbeta))
}



esti.confounder <- function(dta=NULL, subset=NULL, y,x,v=NULL, nfact){
  if(is.null(dta)){
    if(is.null(subset)) subset <- 1:nrow(as.matrix(x))
    y <- as.matrix(y)[subset,]
    x <- as.matrix(x)[subset,]
  }else{
    if(is.null(subset)) subset <- 1:nrow(as.matrix(dta))
    y <- as.matrix(dta[subset,y])
    x <- as.matrix(dta[subset,x])
  }
  
  n <- nrow(as.matrix(x))
  dimx <- ncol(as.matrix(x))
  nfact <- nfact
  
  # factor analysis
  if(is.null(v)){
    x.res <- x
    # ols Y~X
    # lm.y <-lm(y~x) 
  }else{
    if(is.null(dta)){
      v <- v[subset,]
    }else{
      v <- dta[subset, v]
    }
    x.res <- lm(x~v)$res
    # ols Y~X+V
    # lm.y <-lm(y~x+v) 
  }
  res.factor <- try(factanal(x.res, factors=nfact))
  # Check if there was an error
  if (inherits(res.factor, "try-error")) {
    return(rep(NA, dimx))
  }
  efa.fit <- as.matrix(res.factor$loadings)
  halpha <- efa.fit*sqrt(diag(var(x.res)))
  hgamma <- solve(var(x.res))%*%halpha
  hatU <- x.res%*%hgamma
  
  if(is.null(v)){
    regressor <- cbind(hatU,x)
    colnames(regressor) <- NULL
    abess_fit <- abess(regressor, y, family="gaussian", tune.type = "bic",always.include=c(1:(nfact)))
    hbeta = coef(abess_fit, support.size =
                   abess_fit$support.size[which.min(abess_fit$tune.value[which(abess_fit$support.size <= 
                                                                                 (dimx - nfact)/2 + nfact+1)])])[-c(1:(nfact+1))]
  } else{
    abess_fit <- abess(cbind(hatU,v,x), y, family="gaussian", tune.type = "bic",always.include=c(1:(nfact+ncol(v))))
    hbeta = coef(abess_fit, support.size =
                   abess_fit$support.size[which.min(abess_fit$tune.value[which(abess_fit$support.size <= 
                                                                                 (dimx - nfact)/2 + nfact+ncol(v)+1)])])[-c(1:(nfact+ncol(v)+1))]
  }
  
  return(c(hbeta))
}


n <- 200
p <- 10
nfact <- 1

set.seed(8)
alpha <- matrix(rnorm(p*nfact), nrow=p, ncol=nfact)
alpha <- sweep(alpha, 2, colMeans(alpha), "-")
para <- list(dim = p, nfact = nfact, beta = c(c(1, 2), rep(0, p-2)),
             delta = rep(2, nfact), alpha = alpha)


nsims <- 100
result1 <- matrix(NA, nrow = nsims, ncol = p)
result2 <- matrix(NA, nrow = nsims, ncol = p)
result3 <- matrix(NA, nrow = nsims, ncol = p)
for (i in 1:nsims) {
  set.seed(20261100 + i)
  data <- gene.data(n, para)
  result1[i,] <- esti.compositional(dta=data, y=1, x=(2:(p+1)), nfact = nfact)
  colnames(data)[2:(p+1)] <- paste0("V", seq_len(p))
  abess_fit <- abess(data[, 2:(p+1)], data[,1], family="gaussian", tune.type = "bic")
  result2[i,]<- coef(abess_fit, support.size =
                    abess_fit$best.size)[-1]
  # result2[i,] <- esti.confounder(dta=data, y=1, x=(2:(p+1)), nfact = nfact)
  result3[i,] <- esti.confounder(dta=data, y=1, x=((2+p):(1+2*p)), nfact = nfact)
}

# bias <- result1 - rep(1,nsims)%*%t(para$beta)
# bias2 <- result2 - rep(1,nsims)%*%t(para$beta)
# bias3 <- result3 - rep(1,nsims)%*%t(para$beta)

group_names <- expression(beta[1], beta[2], beta[3], beta[4], beta[5]
                          ,beta[6], beta[7], beta[8], beta[9], beta[10]
                          )
boxplot(result1, names = group_names, outline = FALSE)
boxplot(result2, names = group_names, outline = FALSE)
boxplot(result3, names = group_names, outline = FALSE)
# boxplot(result1, names = group_names, outline = FALSE)

# data <- gene.data(n, para)
# result1 <- esti.compositional(dta=data, y=1, x=(2:(p+1)), nfact = nfact)
# result2 <- esti.confounder(dta=data, y=1, x=((2+p):(1+2*p)), nfact = nfact)
# 
# result1
# result2
