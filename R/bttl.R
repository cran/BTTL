
bttl <- function(y, M, c_s){

  bt_simple <- function(A, y){
    loglik = 0
    # assume y is a single list
    n <- nrow(y)   # number of rankers for group k
    alpha <- A
    for(i in 1:n){
      ranking = y[i, ]
      lambda = exp(alpha)
      lsum = sum(exp(alpha[ranking]))
      # plackettluce log likelihood, where objects are evaluated in order of rank
      for(j in 1:length(ranking)){
        loglik = loglik + log(lsum) - alpha[as.numeric(ranking[j])]
        lsum = lsum - lambda[as.numeric(ranking[j])] # remove already ranked item from denominator in lik
      }
    }
    loglik = loglik + 10*abs(sum(alpha)) # sum to zero restriction
    return(loglik)
  }

  loglik_bt <- function(alpha, y){
    loglik = 0
    n <- nrow(y)
    for(i in 1:n){
      ranking = y[i, ]
      lambda = exp(alpha)
      lsum = sum(exp(alpha[ranking]))
      # plackettluce log likelihood, where objects are evaluated in order of rank
      for(j in 1:length(ranking)){
        loglik = loglik + log(lsum) - alpha[as.numeric(ranking[j])]
        lsum = lsum - lambda[as.numeric(ranking[j])] # remove already ranked item from denominator in lik
      }
    }
    return(loglik)
  }

  oracle_bt <- function(y, M, informative){

    oracle_bt_step_1 <- function(U, y, informative){
      u <- U
      loglik = 0

      if(length(informative) > 0){
        y_informative <- vector("list", length(informative))
        for(s in 1:length(informative)){
          y_informative[[s]] <- y[[informative[s]]]
        }
        y_informative <- do.call(rbind, y_informative)
        y_informative <- rbind(y[[1]], y_informative)
      } else{
        y_informative <- y[[1]]
      }
      n <- nrow(y_informative)
      for(i in 1:n){
        ranking = y_informative[i, ]
        lambda = exp(u)
        lsum = sum(lambda)
        # plackettluce log likelihood, where objects are evaluated in order of rank
        for(j in 1:length(ranking)){
          loglik = loglik + log(lsum) - u[as.numeric(ranking[j])]
          lsum = lsum - lambda[as.numeric(ranking[j])] # remove already ranked item from denominator in lik
        }
      }
      # normalize
      loglik = loglik + 10*abs(sum(u)) # sum to zero restriction
      return(loglik)
    }

    oracle_bt_step_2 <- function(delta, y, u_est){
      loglik = 0
      p <- length(y)
      # columns 1 is target data
      n_0 = nrow(y[[1]])   # number of rankers for group k
      for(i in 1:n_0){
        ranking = y[[1]][i, ]
        lambda = exp(u_est + delta)
        lsum = sum(exp((u_est + delta)[ranking]))
        # plackettluce log likelihood, where objects are evaluated in order of rank
        for(j in 1:length(ranking)){
          loglik = loglik + log(lsum) - (u_est + delta)[as.numeric(ranking[j])]
          lsum = lsum - lambda[as.numeric(ranking[j])] # remove already ranked item from denominator in lik
        }
      }
      # penalize delta
      loglik = loglik + 50*sum((delta)^2) # l2 penalty

      # normalize
      loglik = loglik + 10*abs(sum(u_est + delta)) # sum to zero restriction
      #exp(alpha_est + delta) <- exp(alpha_est + delta)/sum(exp(alpha_est + delta))
      return(loglik)
    }

    U <- matrix(rep(0, M), ncol = 1)
    result <- optim(par = U, oracle_bt_step_1, method = "BFGS", y = y, informative = informative)
    u_est <- result$par

    delta <- matrix(rep(0, M), ncol = 1)
    result <- optim(par = delta, oracle_bt_step_2, method = "BFGS", u_est = u_est, y = y)
    delta_est <- result$par

    alpha_est <- u_est+delta_est
    return(alpha_est)
  }

  # first element of y is the primary variable
  # sample rows for splitting of primary data
  y_0_split_idx <- vector("list", 3)
  y_0_split_idx[[1]] <- sample(1:nrow(y[[1]]), ceiling(nrow(y[[1]])/3), replace = FALSE)
  y_0_split_idx[[2]] <- sample(c(1:nrow(y[[1]]))[-y_0_split_idx[[1]]], floor(nrow(y[[1]])/3), replace = FALSE)
  y_0_split_idx[[3]] <- c(1:nrow(y[[1]]))[-c(y_0_split_idx[[1]], y_0_split_idx[[2]])]

  # split primary data
  y_0_split <- vector("list", 3)
  for(i in 1:3){
    y_0_split[[i]] <- y[[1]][y_0_split_idx[[i]],]
  }
  # empty list of secondary attribute worth estimates
  alpha_s_list <- vector("list", 3)
  loglik_s_list <- vector("list", 3)
  for(q in 1:3){
    alpha_s_list[[q]] <- vector("list", length(y)-1)
    loglik_s_list[[q]] <- vector("list", length(y)-1)
  }

  # empty list of secondary attribute worth estimates
  alpha_0_list <- vector("list", 3)
  loglik_0 <- numeric(3)

  # compute estimates
  for(q in 1:3){
    A <- matrix(rep(0, M), ncol = 1)
    result <- optim(par = A, bt_simple, method = "BFGS", y = y[[1]][-y_0_split_idx[[q]],])
    alpha_0_list[[q]] <- result$par

    for(s in 2:length(y)){
      A <- matrix(rep(0, M), ncol = 1)
      result <- optim(par = A, bt_simple, method = "BFGS", y = rbind(y[[1]][-y_0_split_idx[[q]],], y[[s]]))
      alpha_s_list[[q]][[s-1]] <- result$par
      loglik_s_list[[q]][[s-1]] <- loglik_bt(alpha_s_list[[q]][[s-1]], y[[1]][y_0_split_idx[[q]],])
    }
    loglik_0[q] <- loglik_bt(alpha_0_list[[q]], y[[1]][y_0_split_idx[[q]],])
  }

  loglik_s <- vector("list", length(y)-1)
  for(q in 1:3){
    for(s in 1:(length(y)-1)){
      loglik_s[[s]][q] <- loglik_s_list[[q]][[s]]
    }
  }
  loglik_0_mu <- mean(loglik_0)

  loglik_s_mu <- numeric(3)
  for(s in 1:(length(y)-1)){
    loglik_s_mu[s] <- mean(loglik_s[[s]])
  }

  sigma <- numeric(3)
  for(q in 1:3){
    sigma[q] <- (loglik_0[q] - loglik_0_mu)^2
  }
  sigma <- sqrt(mean(sigma))

  S_hat <- numeric(length(y)-1)
  for(s in 1:(length(y)-1)){
    S_hat[s] <- abs(loglik_s_mu[s] - loglik_0_mu) <= c_s*(max(sigma, floor(nrow(y[[1]])/3)*0.01))
  }
  S_hat <- c(2:length(y))[S_hat == 1]

  alpha_est <- oracle_bt(y, M, informative = S_hat)

  return(list(alpha_est = alpha_est, S_hat = S_hat))
}

