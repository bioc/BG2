binomial_PQL <- function(Y,X_sig1=NULL, Beta, Z, Alpha){
  if(is.null(X_sig1)){
    exp_value <- exp(Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))
  }else{
    exp_value <- exp(X_sig1%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))
  }
  mu <- exp_value/(1+exp_value)
  V <- drop(exp_value/(1+exp_value)^2)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(exp_value)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

poisson_PQL <- function(Y,X_sig1=NULL, Beta, Z, Alpha, replicates){
  if(is.null(replicates)){
    replicates = 1
  }

  if(is.null(X_sig1)){
    mu <- exp(Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))*replicates
  }else{
    mu <- exp(X_sig1%*%Beta+ as.matrix(rowSums(mapply(function(x,y){x%*%y}, Z,Alpha)), ncol=1))*replicates
  }

  V <- drop(mu)
  inv_V_vector <- 1/V
  y_star <- inv_V_vector*(Y-mu)+log(mu)-log(replicates)
  return(list(mu=mu, V=V, inv_V_vector=inv_V_vector, y_star=y_star))
}

PQL <- function(Y, kinship, Z, X=NULL, Xs=NULL, indices_X=NULL, family, replicates=NULL, postprob=NULL){

  # kinship (covariance matrix for random effects): list
  # Z (design matrix for random effects): list
  # X: SNPs from screening
  # Xs: other covariates

  # family & default link function
  # binomial	(link = "logit")
  # poisson	(link = "log")

  # replicates: num of measurements for each individual (for poisson distribution)

  if(!is.null(X) | !is.null(Xs)) {
    if(!is.null(X)){
      X_sig <- X
      indices <- indices_X

      if(ncol(X_sig)>1){
        if(Matrix::rankMatrix(X_sig)[1] < ncol(X_sig)){
          dropped_cols <- caret::findLinearCombos(X_sig)$remove
          X_sig <- X_sig[,-dropped_cols]
          indices <- indices[-dropped_cols]
        }

      }

      if(ncol(X_sig)>1){
        total.p <- ncol(X_sig)
        redundancy.snp.indices <- c()
        for(snpi in 1:(total.p-1)){
          for(snpj in (snpi+1):total.p){
            if(!(snpi %in% redundancy.snp.indices)){
              if(stats::cor.test(X_sig[,snpi],X_sig[,snpj],method = "spearman",exact=FALSE,alternative = "greater")$p.value < 0.05){
                if(postprob[indices[snpi]] < postprob[indices[snpj]]){
                  redundancy.snp.indices <- c(redundancy.snp.indices,snpi)
                }else{
                  redundancy.snp.indices <- c(redundancy.snp.indices,snpj)
                }
              }
            }
          }
        }
        X_sig <- X_sig[,!((1:total.p) %in% redundancy.snp.indices)]
        indices <- indices[!((1:total.p) %in% redundancy.snp.indices)]

      }

      if(is.null(Xs)){
        X_sig1 <- cbind(1,X_sig)
      }else{
        X_sig1 <- cbind(1,Xs,X_sig)
      }
    }else{
      X_sig1 <- cbind(1,Xs)
    }

    glmfit <- stats::glm(Y~X_sig1[,-1], family = family)
    Beta <- glmfit$coefficients
    n <- length(Y)
    n_rf <- length(kinship)
    Alpha <- list()
    Kappa <- list()
    Kappa_temp <- list()
    for(i_rf in 1:n_rf){
      Alpha[[i_rf]] <- matrix(rep(0, nrow(kinship[[i_rf]])), ncol = 1)
      Kappa[[i_rf]] <- 0
      Kappa_temp[[i_rf]] <- 0
    }
    Beta_temp <- matrix(0, nrow = ncol(X_sig1), ncol = 1)


    while(sum(mapply(function(x,y){abs(x-y)>0.0001},Kappa, Kappa_temp))>=1 | max(abs(Beta-Beta_temp))>0.0001){

      Kappa_temp <- Kappa
      Beta_temp <- Beta


      if(family=="binomial"){
        par_est <- binomial_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha)
      }
      if(family=="poisson"){
        par_est <- poisson_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha, replicates=replicates)
      }

      mu <- par_est$mu
      V <- par_est$V
      inv_V_vector <- par_est$inv_V_vector
      y_star <- par_est$y_star


      logl <- function(x){
        # log posterior of kappa
        H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, x, kinship, Z)), ncol=n, nrow=n)
        diag(H) <- diag(H)+inv_V_vector
        eign <- eigen(H, symmetric = TRUE)
        P <- eign$vectors
        D <- eign$values
        D_inv <- 1/D
        D_inv[!is.finite(D_inv)] <- 0
        X.significant_1.tilde <- t(P)%*%X_sig1  # t(P) %*% (1 X)

        re.tilde <- t(P)%*%y_star-X.significant_1.tilde%*%Beta

        l <- (+0.5*sum(log(D_inv))   #log(det(H))
              -0.5*determinant(t(X.significant_1.tilde)%*%(D_inv*X.significant_1.tilde), logarithm=TRUE)$modulus[1] #t(X) %*% inv_H %*% X = t(X.tilde) %*% D_inv %*% X.tilde
              -0.5*sum(re.tilde^2*D_inv)) # (y-Xbeta)invH(y-Xbeta)
        return(l)
      }

      ts <- stats::optim(par=rep(0.1,n_rf), fn=logl, lower = rep(0,n_rf), upper = rep(50,n_rf), method = "L-BFGS-B", control = list(fnscale=-1))$par
      Kappa <- as.list(ts)

      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, kinship, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      D <- eign$values
      P <- eign$vectors
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0

      X.significant_1.tilde <- t(P)%*%X_sig1
      y.tilde <- t(P)%*%y_star
      #update beta and alpha
      xdx <- t(X.significant_1.tilde)%*%(D_inv*X.significant_1.tilde)
      eign <- eigen(xdx, symmetric = TRUE)
      D_xdx <- eign$values
      P_xdx <- eign$vectors
      D_inv_xdx <- 1/D_xdx
      D_inv_xdx[!is.finite(D_inv_xdx)] <- 0
      inv_xdx <- P_xdx%*%(D_inv_xdx*t(P_xdx))
      Beta <- inv_xdx%*%t(X.significant_1.tilde)%*%matrix(D_inv*y.tilde,ncol = 1)
      re.tilde <- y.tilde-X.significant_1.tilde%*%Beta   #t(P)%*%(y_star-XBeta)

      Alpha <- mapply(function(x,y,z){x*y%*%t(z)%*%P%*%matrix(D_inv*re.tilde,ncol=1)}, Kappa, kinship, Z, SIMPLIFY = FALSE)
      #Kappa*kinship%*%t(Z)%*%inv_H%*%(y_star-XBeta)

    }


    if(family=="binomial"){
      par_est <- binomial_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha)
    }
    if(family=="poisson"){
      par_est <- poisson_PQL(Y=Y,X_sig1=X_sig1, Beta=Beta, Z=Z, Alpha=Alpha, replicates=replicates)
    }

    mu <- par_est$mu
    V <- par_est$V
    inv_V_vector <- par_est$inv_V_vector
    y_star <- par_est$y_star

    return(list(y_star=y_star, kappa=Kappa, H=H, P=P, D_inv=D_inv, beta=Beta, inv_v=inv_V_vector, X_sig1=X_sig1))

  }else{

    glmfit <- stats::glm(Y~1, family = family)
    Beta <- glmfit$coefficients
    n <- length(Y)
    n_rf <- length(kinship)
    Alpha <- list()
    Kappa <- list()
    Kappa_temp <- list()
    for(i_rf in 1:n_rf){
      Alpha[[i_rf]] <- matrix(rep(0, nrow(kinship[[i_rf]])), ncol = 1)
      Kappa[[i_rf]] <- 0
      Kappa_temp[[i_rf]] <- 0
    }
    Beta_temp <- 0

    while(sum(mapply(function(x,y){abs(x-y)>0.0001},Kappa, Kappa_temp))>=1 | abs(Beta-Beta_temp)>0.0001){
      Kappa_temp <- Kappa
      Beta_temp <- Beta

      if(family=="binomial"){
        par_est <- binomial_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha)
      }
      if(family=="poisson"){
        par_est <- poisson_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha, replicates=replicates)
      }

      mu <- par_est$mu
      V <- par_est$V
      inv_V_vector <- par_est$inv_V_vector
      y_star <- par_est$y_star

      logl <- function(x){
        # log posterior of kappa
        H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, x, kinship, Z)), ncol=n, nrow=n)
        diag(H) <- diag(H)+inv_V_vector
        eign <- eigen(H, symmetric = TRUE)
        P <- eign$vectors
        D <- eign$values
        D_inv <- 1/D
        D_inv[!is.finite(D_inv)] <- 0
        one.tilde <- colSums(P)    # t(P) %*% 1
        #inv_H <- P%*%(D_inv*t(P))
        re.tilde <- t(P)%*%(y_star-Beta)

        l <- (+0.5*sum(log(D_inv)) #log(det(H))
              -0.5*log(sum(one.tilde^2*D_inv)) # t(1) %*% inv H %*% 1 = t(1) %*% P %*% inv D %*% t(P) %*% 1
              -0.5*sum(re.tilde^2*D_inv))  #t(re) %*% inv H %*% re = t(re.tilde) %*% inv D %*% tilde
        return(l)
      }

      ts <- stats::optim(par=rep(0.1,n_rf), fn=logl, lower = rep(0,n_rf), upper = rep(30,n_rf), method = "L-BFGS-B", control = list(fnscale=-1))$par
      Kappa <- as.list(ts)

      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, kinship, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      D <- eign$values
      P <- eign$vectors
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0

      X.tilde <- colSums(P)
      y.tilde <- t(P)%*%y_star
      #update beta and alpha
      xdx <- sum(X.tilde^2*D_inv)
      inv_xdx <- 1/xdx
      Beta <- inv_xdx*sum(X.tilde*D_inv*y.tilde)

      re.tilde <- y.tilde-X.tilde*Beta   #t(P)%*%(y_star-XBeta)
      Alpha <- mapply(function(x,y,z){x*y%*%t(z)%*%P%*%matrix(D_inv*re.tilde,ncol=1)}, Kappa, kinship, Z, SIMPLIFY = FALSE)
      #Kappa*kinship%*%t(Z)%*%inv_H%*%(y_star-XBeta)
    }

    if(family=="binomial"){
      par_est <- binomial_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha)
    }
    if(family=="poisson"){
      par_est <- poisson_PQL(Y=Y,X_sig1=NULL, Beta=Beta, Z=Z, Alpha=Alpha, replicates=replicates)
    }

    mu <- par_est$mu
    V <- par_est$V
    inv_V_vector <- par_est$inv_V_vector
    y_star <- par_est$y_star

    return(list(y_star=y_star, kappa=Kappa, H=H, P=P, D_inv=D_inv, beta=Beta, inv_v=inv_V_vector))
  }



}

P3D <- function(Y, SNPs, Xs=NULL, kinship, Z, family, replicates = NULL){
  n = length(Y)
  PQL_est <- PQL(Y=Y, Xs=Xs, kinship=kinship, Z=Z, family=family, replicates=replicates, postprob=NULL)

  y_star <- PQL_est$y_star
  Kappa <- PQL_est$kappa
  Beta <- PQL_est$beta
  inv_V_vector <- PQL_est$inv_v
  H <- PQL_est$H
  P <- PQL_est$P
  D_inv <- PQL_est$D_inv

  if(is.null(Xs)){

    A <- H-sum(P%*%(D_inv*t(P)))
    eign <- eigen(A, symmetric = TRUE)
    P <- eign$vectors
    y.tilde = t(P) %*% (y_star-Beta)
    D <- eign$values
    D_inv <- D^(-1)
    #D_inv[!is.finite(D_inv)] <- 0
    D_inv[n] <- 0
  }else{
    X_1 <-  PQL_est$X_sig1
    A <- H-X_1%*%solve(t(X_1)%*%P%*%(D_inv*t(P))%*%X_1)%*%t(X_1)
    eign <- eigen(A, symmetric = TRUE)
    P <- eign$vectors
    y.tilde = t(P) %*% (y_star-X_1%*%Beta)
    D <- eign$values
    D_inv <- D^(-1)
    D_inv[(n-ncol(X_1)+1):n] <- 0
  }

  X.tilde = t(P) %*% SNPs
  # Estimate beta with beta.hat and compute var(beta.hat)
  xj.t.xj = apply(X.tilde*(D_inv*X.tilde),2,sum)
  xj.t.y = t(X.tilde) %*% (D_inv*y.tilde)

  beta.hat = xj.t.y / xj.t.xj
  var.beta.hat = 1 / xj.t.xj
  t.statistic <- beta.hat / sqrt(var.beta.hat)
  pvalues <- 2*stats::pnorm(abs(t.statistic), mean=0, sd=1, lower.tail=FALSE)
  return_dat <- cbind(beta.hat,var.beta.hat,pvalues)
  return_dat <- as.data.frame(return_dat)
  colnames(return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values")
  return(return_dat)
}

log_marginal_likelihood_null <- function(y.tilde,D_inv){
  # This function computes the log marginal likelihood for the case
  # when there is no regressor in the model. y.tilde ~ N(0,D)
  n <- length(y.tilde)
  return(-0.5*n*log(2*pi)+0.5*sum(log(D_inv))-0.5*sum(D_inv*y.tilde^2))
}

log.marginal.likelihood.type1 <- function(k, x.tilde_m, y.tilde, D_inv, ydinvy, dinvy, C.k_m, tau, nboot)
{

  n <- length(y.tilde)     # sample size
  enum <- nboot
  #k <- ncol(x.tilde_m)
  if(is.matrix(C.k_m))
  {
    inv_C.k = solve(C.k_m)
    beta.k.hat = inv_C.k %*% t(x.tilde_m) %*% dinvy
    R.k <- ydinvy - t(dinvy)%*%x.tilde_m%*%beta.k.hat
    E2 <- MASS::mvrnorm(enum, mu = beta.k.hat, Sigma = inv_C.k)
    E1 <- t((t(E2)-drop(beta.k.hat))*sqrt(1+tau*n))

  } else {
    inv_C.k = 1/C.k_m
    beta.k.hat = inv_C.k %*% t(x.tilde_m) %*% dinvy
    R.k <- ydinvy - t(dinvy)%*%x.tilde_m%*%beta.k.hat
    E1 <- matrix(stats::rnorm(enum*length(beta.k.hat), mean = 0, sd = sqrt((1+tau*n)*inv_C.k)), ncol = length(beta.k.hat))
    E2 <- t((t(E1)+drop(beta.k.hat))/sqrt(1+tau*n))

  }

  E2 <- mean(apply(E2^2, 1, prod))
  E1 <- mean(apply(E1^2, 1, prod))

  return( -0.5*n*log(2*pi) - 0.5*k*log(1+tau*n) + 0.5*sum(log(D_inv))
          -0.5*R.k +log(E2) - log(E1))
}

BG2_terminal <- function(Y, Xs=NULL, SNPs,Z,kinship, P3D_return_dat, Tau = 1 ,family, FDR.threshold = 0.95,replicates=NULL,
                         maxiterations = 4000, runs_til_stop = 500,nboot = 2000){
  n= length(Y)
  #estimate tau and pi0
  beta.hat <- P3D_return_dat$Beta_Hat
  var.beta.hat <- P3D_return_dat$Var_Beta_Hat

  if(Tau == 1){
    #uniform dist for tau
    log.marginal.likelihood_tau = function(param)
    {
      pi0 = exp(param[1]) / (1+exp(param[1]))
      tau = exp(param[2])
      return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                       (1-pi0)*(2*pi*var.beta.hat)^(-0.5) * (n*tau+1)^(-1.5) *
                       (1+n*tau*beta.hat^2/((n*tau+1)*var.beta.hat)) *
                       exp(-beta.hat^2 / (2*var.beta.hat*(n*tau+1)))
      )))

    }

    param.start <- c(2,-2)
    result <- stats::optim(param.start, fn=log.marginal.likelihood_tau, lower = c(1,-20), method = "L-BFGS-B", hessian=TRUE, control = list(fnscale=-1))

    pi0.hat = exp(result$par[1]) / (1+exp(result$par[1]))
    tau.hat = exp(result$par[2])
  }else if(Tau == 2){
    #inverse gamma prior for tau
    log.marginal.likelihood_tau = function(param)
    {
      pi0 = exp(param[1]) / (1+exp(param[1]))
      tau = exp(param[2])
      return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                       (1-pi0)*(2*pi*var.beta.hat)^(-0.5) * (n*tau+1)^(-1.5) *
                       (1+n*tau*beta.hat^2/((n*tau+1)*var.beta.hat)) *
                       exp(-beta.hat^2 / (2*var.beta.hat*(n*tau+1)))
      )) +stats::dgamma(x=1/tau, shape = 0.55/0.022+1, scale = 1/0.55,log = TRUE))

    }

    param.start <- c(2,-2)
    result <- stats::optim(param.start, fn=log.marginal.likelihood_tau, lower = c(1,-20), method = "L-BFGS-B", hessian=TRUE, control = list(fnscale=-1))

    pi0.hat = exp(result$par[1]) / (1+exp(result$par[1]))
    tau.hat = exp(result$par[2])
  }else{
    log.marginal.likelihood_tau = function(param)
    {
      pi0 = exp(param[1]) / (1+exp(param[1]))
      return(sum(log(pi0*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
                       (1-pi0)*(2*pi*var.beta.hat)^(-0.5) * (n*Tau+1)^(-1.5) *
                       (1+n*Tau*beta.hat^2/((n*Tau+1)*var.beta.hat)) *
                       exp(-beta.hat^2 / (2*var.beta.hat*(n*Tau+1)))
      )))

    }

    param.start <- c(2)
    result <- stats::optim(param.start, fn=log.marginal.likelihood_tau, method = "L-BFGS-B", hessian=TRUE, control = list(fnscale=-1))

    pi0.hat = exp(result$par[1]) / (1+exp(result$par[1]))
    tau.hat = Tau
  }

  if( (pi0.hat > 0.6) & (tau.hat < 0.5) & (tau.hat > 10^(-10))){

    # Compute posterior probability of beta_j different than 0:
    numerator <- (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
      (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
      exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))
    denominator <- pi0.hat*stats::dnorm(beta.hat,mean= 0,sd=sqrt(var.beta.hat)) +
      (1-pi0.hat)*(2*pi*var.beta.hat)^(-0.5) * (n*tau.hat+1)^(-1.5) *
      (1+n*tau.hat*beta.hat^2/((n*tau.hat+1)*var.beta.hat)) *
      exp(-beta.hat^2 / (2*var.beta.hat*(n*tau.hat+1)))

    postprob <- numerator / denominator

    order.postprob <- order(postprob, decreasing=TRUE)
    postprob.ordered <- postprob[order.postprob]

    FDR.Bayes <- cumsum(postprob.ordered) / 1:ncol(SNPs)
    if(sum(FDR.Bayes > FDR.threshold) == 0){
      P3D_return_dat <- cbind(P3D_return_dat,postprob,FALSE)
      P3D_return_dat <- as.data.frame(P3D_return_dat)
      colnames(P3D_return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values","PostProb","Significant")
    }else{
      P3D_return_dat <- cbind(P3D_return_dat,postprob,postprob >= postprob.ordered[max(which(FDR.Bayes > FDR.threshold))])
      P3D_return_dat <- as.data.frame(P3D_return_dat)
      colnames(P3D_return_dat) <- c("Beta_Hat","Var_Beta_Hat","P_Values","PostProb","Significant")
    }

    if(sum(P3D_return_dat$Significant) > 0){
      indices_X <- which(P3D_return_dat$Significant)
      X <- SNPs[, indices_X,drop = FALSE]


      pi0 = pi0.hat
      tau = tau.hat

      PQL_est <- PQL(Y, kinship, Z = Z, X=X, Xs=Xs, indices_X=indices_X, family, replicates=replicates, postprob = postprob)

      y_star <- PQL_est$y_star
      Kappa <- PQL_est$kappa
      Beta <- PQL_est$beta
      inv_V_vector <- PQL_est$inv_v

      H <- matrix(rowSums(mapply(function(x,y,z){x*z%*%y%*%t(z)}, Kappa, kinship, Z)), ncol=n, nrow=n)
      diag(H) <- diag(H)+inv_V_vector
      eign <- eigen(H, symmetric = TRUE)
      P <- eign$vectors
      if(is.null(Xs)){
        y.tilde = t(P) %*% (y_star-Beta[1])

      }else{
        y.tilde = t(P) %*% (y_star-cbind(1,Xs)%*%Beta[1:(ncol(Xs)+1)])
      }
      X.tilde = t(P) %*% X
      D <- eign$values

      #y.tilda~N(X.tilde*Beta,D)  # full model, no intercept in X.tilde
      # elements with all k regressors
      D_inv <- 1/D
      D_inv[!is.finite(D_inv)] <- 0
      dinvy <- D_inv*y.tilde
      ydinvy <- sum(D_inv*y.tilde^2)
      xdinvx <- t(X.tilde)%*%(D_inv*X.tilde)
      C.k <- xdinvx*(1+1/tau/n)

      total.p <- ncol(X.tilde)

      if(total.p < 16){
        # Do full model search
        total.models <- 2^total.p
        log.unnormalized.posterior.probability <- rep(NA, total.models)
        log.unnormalized.posterior.probability[1] <- total.p * log(pi0) + log_marginal_likelihood_null(y.tilde,D_inv)
        dat <- rep(list(0:1), total.p)
        dat <- as.matrix(expand.grid(dat))
        for (i in 1:(total.models-1)){
          model <- unname(which(dat[i + 1,] == 1))
          k <- length(model)

          Xsub <- X.tilde[,model,drop = FALSE]
          if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
            dropped_cols <- caret::findLinearCombos(Xsub)$remove
            model <- model[-dropped_cols]
          }

          x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))
          C.k_m <- C.k[model, model]
          log.unnormalized.posterior.probability[i+1] <-
            k*log(1-pi0) + (total.p-k)*log(pi0) +
            log.marginal.likelihood.type1(k, x.tilde_m, y.tilde, D_inv, ydinvy, dinvy, C.k_m, tau, nboot)

        }
        log.unnormalized.posterior.probability <- log.unnormalized.posterior.probability - max(log.unnormalized.posterior.probability)
        unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
        posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)

      }else {
        # Do model search with genetic algorithm
        fitness_ftn <- function(string){
          if(sum(string) == 0){
            return(total.p * log(pi0) + log_marginal_likelihood_null(y.tilde,D_inv))
          }else{
            model <- which(string==1)
            k <- length(model)

            Xsub <- X.tilde[,model,drop = FALSE]
            if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
              dropped_cols <- caret::findLinearCombos(Xsub)$remove
              model <- model[-dropped_cols]
            }

            x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))
            C.k_m <- C.k[model, model]

            return(k*log(1-pi0) + (total.p-k)*log(pi0) +
                     log.marginal.likelihood.type1(k, x.tilde_m, y.tilde, D_inv, ydinvy, dinvy, C.k_m, tau, nboot)
            )
          }
        }

        if(total.p > 99){
          suggestedsol <- diag(total.p)
          tmp_log.unnormalized.posterior.probability <- vector()
          for(i in 1:total.p){
            model <- which(suggestedsol[i,]==1)
            k <- length(model)

            Xsub <- X.tilde[,model,drop = FALSE]
            if(Matrix::rankMatrix(Xsub)[1] < ncol(Xsub)){
              dropped_cols <- caret::findLinearCombos(Xsub)$remove
              model <- model[-dropped_cols]
            }

            x.tilde_m <- matrix(X.tilde[,model], ncol = length(model))
            C.k_m <- C.k[model, model]
            tmp_log.unnormalized.posterior.probability[i] <- (k*log(1-pi0) + (total.p-k)*log(pi0) +
                                                                log.marginal.likelihood.type1(k, x.tilde_m, y.tilde, D_inv, ydinvy, dinvy, C.k_m, tau, nboot))

          }
          suggestedsol <- rbind(0,suggestedsol[order(tmp_log.unnormalized.posterior.probability,decreasing = TRUE)[1:99],])
        }else{
          suggestedsol <- rbind(0,diag(total.p))
        }

        # maxiterations = 4000
        # runs_til_stop = 1000

        fitness_ftn <- memoise::memoise(fitness_ftn)
        ans <- GA::ga("binary", fitness = fitness_ftn, nBits = total.p,maxiter = maxiterations,popSize = 100,
                      elitism = min(c(10,2^total.p)),run = runs_til_stop,suggestions = suggestedsol,monitor = FALSE)
        memoise::forget(fitness_ftn)
        dat <- ans@population
        dupes <- duplicated(dat)
        dat <- dat[!dupes,]
        ans@fitness <- ans@fitness[!dupes]
        log.unnormalized.posterior.probability <- ans@fitness - max(ans@fitness)
        unnormalized.posterior.probability <- exp(log.unnormalized.posterior.probability)
        posterior.probability <- unnormalized.posterior.probability/sum(unnormalized.posterior.probability)
      }

      inclusion_prb <- unname((t(dat)%*%posterior.probability)/sum(posterior.probability))

      model <- dat[which.max(posterior.probability),]
      model_dat <- cbind(indices_X,model,inclusion_prb)
      model_dat <- as.data.frame(model_dat)
      colnames(model_dat) <- c("SNPs","BestModel","Inclusion_Prob")
      return(list(prescreen = P3D_return_dat,postprob=postprob,modelselection = model_dat,pi_0_hat = pi0.hat, tau_hat = tau.hat))

    }else{
      return(list(prescreen = P3D_return_dat,postprob=postprob,modelselection = "No significant in prescreen1",pi_0_hat = pi0.hat, tau_hat = tau.hat))
    }
  }else{
    return(list(prescreen = P3D_return_dat,modelselection = "No significant in prescreen1",pi_0_hat = pi0.hat, tau_hat = tau.hat))
  }

}

#' Performs BG2 as described in the manuscript, Xu, Williams, and Ferreira BG2: Bayesian variable selection in generalized linear mixed models with non-local priors for non-Gaussian GWAS data, Bioinformatics, Submitted.
#'
#'
#'
#' @param Y The observed phenotypes, count or binary.
#' @param SNPs The SNP matrix, where each column represents a single SNP encoded as the numeric coding 0, 1, 2. This is entered as a matrix object.
#' @param FDR_Nominal The nominal false discovery rate for which SNPs are selected from in the screening step.
#' @param Fixed A matrix of fixed covariates to control for. Do not include the intercept. The value is defaulted at NULL implying no fixed covariates.
#' @param family Specify if the response is count ("poisson") or binary ("bernoulli").
#' @param Covariance A list of covariance matrices that are the covariance matrices of the random effects. This matches the list of design matrices in Z.
#' @param Z A list of matrices specifying the design matrix of each random effect of interest.
#' @param replicates If family = "poisson", the replicates of each ecotype, can be a vector or a number if the number of replicates is the same for each ecotype. If family = "binomial", replicates = NULL.
#' @param Tau Specifying either a fixed value for the dispersion parameter of the nonlocal prior (0.022 and 0.348 are used in the paper). Or specify a prior for tau, either uniform or Inverse Gamma centered at ...
#' @param maxiterations The maximum iterations the genetic algorithm in the model selection step iterates for.
#' Defaulted at 400 which is the value used in the BG2 paper simulation studies.
#' @param runs_til_stop The number of iterations at the same best model before the genetic algorithm in the model selection step converges.
#' Defaulted at 40 which is the value used in the BG2 paper simulation studies.
#' @return The column indices of SNPs that were in the best model identified by BG2.
#' @examples
#' library(BG2)
#' data("Y_poisson");data("SNPs");data("kinship")
#' n <- length(Y_poisson)
#' covariance <- list()
#' covariance[[1]] <- kinship
#' covariance[[2]] <- diag(1, nrow = n, ncol = n)
#'
#' set.seed(1330)
#' output_poisson <- BG2(Y=Y_poisson, SNPs=SNPs, Fixed = NULL,
#'                    Covariance=covariance, Z=NULL, family="poisson",
#'                    replicates=4, Tau="uniform",FDR_Nominal = 0.05,
#'                    maxiterations = 4000, runs_til_stop = 400)
#'
#' @export
BG2 <- function(Y, SNPs, FDR_Nominal = 0.05, Fixed=NULL, family=c("poisson","bernoulli"),Covariance, Z=NULL,replicates=NULL, Tau="uniform", maxiterations = 4000, runs_til_stop = 400){
  #Y: observations

  #SNPs: all SNP variables

  #Fixed: fixed covariates

  #kinship: a list of covariance matrices of random effects

  #Z: a list of design matrices of random effects

  #family: "bernoulli" or "poisson"

  #replicates: if family = "poisson", the relipcates of each ecotype, can be a vector or a number.
  #            if family = "bernoulli", replicates = NULL

  #Tau: scale parameter of non-local prior
  # Tau = “uniform”: use uniform prior for tau
  # Tau = “IG”: use inverse gamma prior for tau
  # Tau = 0.022: fixed tau at 0.022
  # Tau = 0.348: fixed tau at 0.348

  family <- match.arg(family)

  if(sum(family %in% c("poisson","bernoulli")) == 0){
    stop("family must be either poisson or bernoulli")
  }
  if(FDR_Nominal > 1 | FDR_Nominal < 0){
    stop("FDR_Nominal has to be between 0 and 1")
  }
  if(!is.numeric(Y)){
    stop("Y has to be numeric")
  }
  if(!is.matrix(SNPs)){
    stop("SNPs has to be a matrix object")
  }
  if(!is.numeric(SNPs)){
    stop("SNPs has to contain numeric values")
  }
  if(maxiterations-floor(maxiterations)!=0){
    stop("maxiterations has to be a integer")
  }
  if(runs_til_stop-floor(runs_til_stop)!=0){
    stop("runs_til_stop has to be a integer")
  }
  if(maxiterations < runs_til_stop){
    stop("maxiterations has to be larger than runs_til_stop")
  }

  if(Tau == "uniform"){
    Tau <- 1
  }
  if(Tau == "IG"){
    Tau <- 2
  }
  if(family == "bernoulli"){
    family <- "binomial"
  }
  kinship <- Covariance
  SNPs <- scale(SNPs)
  Xs <- Fixed
  FDR.threshold <- 1 - FDR_Nominal
  nboot <- 1000

  n <- length(Y)
  if(is.null(Z)){
    n_rf <- length(kinship)
    Z <- list()
    for(rf in 1:n_rf){
      Z[[rf]] <- diag(1, ncol = n, nrow = n)
    }
  }


  P3D_return_dat <- P3D(Y=Y, SNPs=SNPs, Xs=Xs, kinship=kinship, Z=Z, family=family, replicates = replicates)

  tmp <- BG2_terminal(Y,Xs=Xs, SNPs = SNPs,Z = Z,kinship = kinship, P3D_return_dat= P3D_return_dat, Tau = Tau ,family = family,replicates = replicates,
                      FDR.threshold = FDR.threshold, maxiterations = maxiterations, runs_til_stop = runs_til_stop, nboot = nboot)


  if(!is.character(tmp$modelselection)){
    indices <- tmp$modelselection$SNPs[tmp$modelselection$BestModel == 1]
    n_indices <- length(indices)
    for(i in 1:n_indices){
      screen_indices <- tmp$modelselection$SNPs
      for(j in 1:length(screen_indices)){
        if(identical(all.equal(stats::cor(SNPs[,indices[i]], SNPs[,screen_indices[j]]),1),TRUE)){
          indices <- c(indices,screen_indices[j])
        }
      }
    }
    return(unique(indices)[order(unique(indices))])

  }else{
    return("No significant SNP")
  }

}

#' A. Thaliana Kinship matrix
#'
#' This is a kinship matrix from the TAIR9 genotype information for 328 A. Thaliana Ecotypes from the paper
#' Components of Root Architecture Remodeling in Response to Salt Stress (Julkowska et al. Genetic Components of Root Architecture Remodeling in Response to Salt Stress, The Plant Cell, Volume 29, Issue 12, December 2017, Pages 3198–3213). The kinship matrix was computed using all SNPs with minor allele frequency
#' greater than 0.01.
#'
#' @format ## `kinship`
#' A matrix with 328 rows and 328 columns corresponding to the 328 ecotypes.
"kinship"

#' A. Thaliana Genotype matrix
#'
#' This is a matrix with 328 observations and 9,000 SNPs. Each row is contains 9,000 SNPs from a single A. Thaliana ecotype in the paper
#' Components of Root Architecture Remodeling in Response to Salt Stress (Julkowska et al. Genetic Components of Root Architecture Remodeling in Response to Salt Stress, The Plant Cell, Volume 29, Issue 12, December 2017, Pages 3198–3213).
#'
#' @format ## `SNPs`
#' A matrix with 328 observations and 9,000 SNPs.
"SNPs"

#' A. Thaliana Simulated Phenotype matrix
#'
#' This is a phenotype matrix simulated from 9,000 SNPs. SNPs at positions 450, 1350, 2250, 3150, 4050,
#' 4950, 5850, 6750, 7650, and 8550 have nonzero coefficients. Further, the data was simulated under a poisson mixed effects model
#' with both a kinship random effect and an overdispersion random effect. The data was simulated using the kinship random effect provided
#' in the package.
#'
#' @format ## `Y_poisson`
#' A vector with 328 observations corresponding to the 328 ecotypes.
"Y_poisson"

#' A. Thaliana Simulated Phenotype matrix
#'
#' This is a phenotype matrix simulated from 9,000 SNPs. SNPs at positions 450, 1350, 2250, 3150, 4050,
#' 4950, 5850, 6750, 7650, and 8550 have nonzero coefficients. Further, the data was simulated under a binary mixed effects model
#' with only a kinship random effect. The data was simulated using the kinship random effect provided in the package.
#'
#' @format ## `Y_binary`
#' A vector with 328 observations corresponding to the 328 ecotypes.
"Y_binary"
