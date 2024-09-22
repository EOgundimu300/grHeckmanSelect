

################################################################################
# Function for group variable selection in FIML estimator of Heckman model.  
# Group Lasso, Adaptive LASSO, SCAD and MCP are implemented. 
# The code can also be used for individual penalization where each variable
# is taken to belong to its own group in a highly computationally efficient way.
###################################################################################

#' Group variable selection in nonignorable missing data using LASSO, SCAD and MCP.
#'
#' @param W A matrix of covariates in selection equation (intercept is not included).
#' @param X A matrix of covariates in outcome equation (intercept is not included).
#' @param s Binary outcome for the selection equation (0/1 and false/true allowed).
#' @param y Continuous outcome for the outcome equation (0 for NA).
#' @param group_sel A vector describing the grouping of the coefficients in the selection equation. 
#' It is best if group is a vector of consecutive integers. If there are coefficients to be included 
#' in the model without being penalized, assign them to group 0.
#' @param group_out A vector describing the grouping of the coefficients in the outcome equation. 
#' It is best if group is a vector of consecutive integers. If there are coefficients to be included 
#' in the model without being penalized, assign them to group 0.
#' @param nlambda The number of lambda values. Default is 100.
#' @param lambda A user supplied sequence of lambda values.
#' Typically, this is left unspecified, and the function automatically 
#' computes a grid of lambda values.
#' @param penalty The penalty to be applied to the model, one of grLasso, grSCAD, or grMCP.
#' @param lambda.min The smallest value for lambda, as a fraction of lambda.max. Default is .001.
#' @param log.lambda When TRUE compute the grid values of lambda on log scale (default) or linear scale otherwise.

#' @param eps Convergence threshhold. Default is 1e-4
#' @param max_iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param gamma Tuning parameter of the group MCP/SCAD penalty. Default is 3 for MCP and 4 for SCAD
#' @param group_multiplier A vector of values representing multiplicative factors by which each group's penalty is to be multiplied.
#' The default is the square root of group size.
#' @param init_strat Use MLE to initialize the group descent algorithm otherwise use zeros.
#' @param penalty.factors Allows for the use of weighted versions of grLasso, grSCAD, or grMCP.


#' @return class grHeckSelect containing optimal penalized coefficients, lambda values,
#' bic values etc.
#' @export
#'
#' @examples
#' #grHeckSelect(W=W, X=X, s=s, y=y, group_sel, group_out, penalty="grLasso")

#'

grHeckSelect <- function(W, X, s, y, group_sel, group_out, penalty="grLasso",
                        nlambda=100,lambda,lambda.min=0.001,log.lambda=TRUE,
                        eps=1e-4,max_iter=10000,
                        gamma=ifelse(penalty=="grSCAD",4,3), group_multiplier,
                        init_strat = "MLE",penalty.factors=NULL){
  require(sampleSelection)
  group <- c(0,group_sel,0,group_out,0,0)
  
  if(is.null(colnames(W))) colnames(W) <- paste("W",1:ncol(W),sep='')
  if(is.null(colnames(X))) colnames(X) <- paste("X",1:ncol(X),sep='')
  joint_dat <- data.frame(cbind(W,X))
  s_form <- as.formula(paste("s ~",paste(colnames(W),collapse=" + ")))
  y_form <- as.formula(paste("y ~",paste(colnames(X),collapse=" + ")))
  nongroup_fit <- selection(s_form, y_form, data = joint_dat)
  est_params <- nongroup_fit$estimate
  pS <- ncol(W) + 1
  pO <- ncol(X) + 1
  est_params[(pS+1):(pS+pO)] <- est_params[(pS+1):(pS+pO)]/est_params[pS+pO+1]
  est_params[pS+pO+1] = -log(est_params[pS+pO+1])
  est_params[pS+pO+2] = atanh(est_params[pS+pO+2])
  hess <- -1*hesslik(est_params,
                     wS=cbind(1,W[s==0,]),
                     wO=cbind(1,W[s==1,]),
                     xO=cbind(1,X[s==1,]),
                     yO=y[s==1])
  hesseig <- eigen(hess)
  pseudoX <- hesseig$vectors %*% diag((hesseig$values) ** (1/2)) %*% t(hesseig$vectors)
  pseudoY <- pseudoX %*% est_params
  
  if(penalty=="grLasso"){
    penalty_fun <- lasso_pen
  } else if(penalty=="grMCP"){
    penalty_fun <- mcp_pen
  } else if(penalty=="grSCAD"){
    penalty_fun <- scad_pen
  } else{
    print("Unrecognized penalty - using grLasso")
    penalty_fun <- lasso_pen
  }
  
  if(is.null(penalty.factors)){
    penalty.factors <- rep(1,length(est_params))
  }
  
  scale <- apply(pseudoX, 2, function(x){sqrt(sum(x**2)/nrow(pseudoX))})
  XX <- pseudoX %*% diag(1/scale)
  
  g_order <- order(group)
  g_order_inv <- match(1:length(group),g_order)
  
  XX_ord <- XX[,g_order]
  
  g <- group[g_order]
  
  XX_orth <- orthogonalize(XX_ord,g)
  K <- as.integer(table(g))
  K1 <- cumsum(K)
  K0 = K1[1]
  
  if(missing(group_multiplier)) group_mult = sqrt(K)[-1]
  
  if(missing(lambda)){
    fity <- glm(pseudoY~XX_orth[, g==0]-1 , family="gaussian")
    r <- fity$residuals
    zmax <- maxgrad(XX_orth, r, K1, group_mult)/nrow(XX_orth)
    
    lambda.max = zmax
    
    if (log.lambda) {
      if (lambda.min==0) {
        lambda <- c(exp(seq(log(lambda.max), log(0.001*lambda.max), length=nlambda-1)), 0)
      } else {
        lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
      }
    } else {
      if (lambda.min==0) {
        lambda <- c(seq(lambda.max, 0.001*lambda.max, length = nlambda-1), 0)
      } else {
        lambda <- seq(lambda.max, lambda.min*lambda.max, length = nlambda)
      }
    }
  }
  
  # if(init_strat=="MLE"){
  #   init_b <- est_params
  #   init_b <- init_b*scale
  #   init_b <- init_b[g_order]
  #   init_b <- orthogonalize_init(init_b, XX_orth, g, intercept=FALSE)
  # } else{
  #   init_b <- rep(0,ncol(XX_orth))
  # }
  
  fity2 <- glm(pseudoY~XX_orth-1, family="gaussian")
  init_b <- fity2$coefficients
  
  
  fit_out <- my_gdfit(XX_orth, pseudoY, penalty_fun, K1, K0, lambda, gamma, eps,
                      max_iter, group_mult, init_b, penalty.factors)
  iterations <- fit_out$iterations
  beta <- fit_out$beta
  beta <- unorthogonalize(beta, XX_orth, g, intercept=FALSE)
  beta <- beta[g_order_inv,]
  beta <- apply(beta,2,function(x){return(x/scale)})
  
  loss <- apply(beta,2,loglik_new,
                wS=cbind(1,W[s==0,]),
                wO=cbind(1,W[s==1,]),
                xO=cbind(1,X[s==1,]),
                yO=y[s==1])
  df <- apply(beta,2,function(x){return(sum(x!=0))}) - 4
  bic <- df*log(nrow(X)) - 2*loss
  
  rownames(beta) <- rownames(nongroup_fit$hess)
  
  beta[pS+pO+2,] = tanh(beta[pS+pO+2,])
  beta[pS+pO+1,] = exp(-beta[pS+pO+1,])
  beta[(pS+1):(pS+pO),] = beta[(pS+1):(pS+pO),]*beta[pS+pO+1,]
  
  best_idx <- which(bic == min(bic))
  opt_params <- beta[,best_idx]
  opt_lambda <- lambda[best_idx]
  opt_bic <- min(bic)
  
  out <- list(params=beta,
              group=group,
              lambda=lambda,
              df=df,
              loss=loss,
              bic=bic,
              penalty=penalty,
              n=nrow(X),
              iter=iterations,
              opt_params=opt_params,
              opt_lambda=opt_lambda,
              opt_bic=opt_bic)
  
              class(out) <- "grHeckSelect"

  return(out)
}


#grlasso <- grHeckSelect(W=W, X=X, s=s, y=y, group_sel, group_out, penalty="grLasso")



#########################################################################
# Coefficients from object of class grHeckSelect
#########################################################################

#' Title Coefficients from objects returned by grHeckSelect function
#'
#' @param object a fitted object of class inheriting from "grHeckSelect"
#' @param ...
#'
#' @return vector of coefficients for the selection and outcome models
#' @export
#'
#' @examples
#' #coef(grHeckSelect(W=W, X=X, s=s, y=y, group_sel, group_out, penalty="grLasso"))



coef.grHeckSelect<- function(object,...){
 x <- object$opt_params
return(x)
}




my_gdfit <- function(X, y, penalty_fun, K1, K0, lambda, gamma, eps, 
                     max_iter, group_mult, init_b, penalty.factors){
  
  n <- nrow(X)
  sdy <- sqrt(sum(y**2)/n)
  tol <- eps*sdy
  output_b <- matrix(data=NA,ncol=length(lambda),nrow=nrow(X))
  iter_list <- rep(0,length(lambda))
  new_init_b <- init_b
  for(l in 1:length(lambda)){
    
    b <- new_init_b
    new_b <- b
    r <- y - X %*% b
    
    lam <- lambda[l]
    converged <- FALSE
    total_iter <- 0
    while(total_iter < max_iter){
      total_iter = total_iter + 1
      maxchange = 0
      
      # This commented bit is the co-ordinate descent for unpenalized parameters that performs very poorly for some reason
      # for(j in 1:K0){
      #   shift = 1/n * sum(X[,j] * r)
      #   new_b[j] = shift + b[j]
      #   r = r - shift * X[,j]
      #   if(abs(shift) > maxchange){
      #     maxchange = abs(shift)
      #   }
      # }
      
      for(g in 1:(length(K1) - 1)){
        j_cols <- (K1[g]+1):K1[g+1]
        lam_temp <- lam * group_mult[g]
        lam_temp = lam_temp * penalty.factors[g]
        list_temp <- my_gd(b, X, r, j_cols, n, penalty_fun, lam_temp, gamma)
        new_b[j_cols] = list_temp$b[j_cols]
        r = list_temp$r
        if(list_temp$maxchange > maxchange){
          maxchange = list_temp$maxchange
        }
      }
      b = new_b
      if(maxchange < tol){
        converged=TRUE
        break
      }
    }
    output_b[,l] = b
    if(converged == FALSE){
      print(paste("Non convergence for", lambda[l]))
      print(b)
    } else{
      new_init_b <- b
    }
    iter_list[l] = total_iter
  }
  return(list(beta = output_b, iterations = iter_list))
}

my_gd <- function(b, X, r, j_cols, n, penalty_fun, lam, gamma){
  a <- b
  Xj <- X[,j_cols,drop=FALSE]
  z <- 1/n * t(Xj) %*% r + a[j_cols]
  z_norm <- sqrt(sum(z**2))
  pen <- penalty_fun(z_norm,lam,gamma)
  maxchange <- 0
  if((pen != 0) || any(a[j_cols]!=0)){
    b[j_cols] = pen * z / z_norm
    shift <- b[j_cols] - a[j_cols]
    r = r - Xj %*% as.matrix(shift,ncol=1)
    maxchange = max(abs(shift))
  }
  return(list(b=b, r=r, maxchange=maxchange))
}

lasso_pen <- function(z,lam,gamma){
  return(ifelse(z > lam, z - lam, 0))
}

scad_pen <- function(z,lam,gamma){
  if(z < lam){
    return(0)
  } else if(z < 2*lam){
    return(z - lam)
  } else if(z < gamma*lam){
    return( (z - (gamma*lam/(gamma - 1))) / (1 - 1/(gamma-1)) )
  } else{
    return(z)
  }
}

mcp_pen <- function(z,lam,gamma){
  if(z < lam){
    return(0)
  } else if(z < gamma*lam){
    return((z-lam)/(1 - 1/gamma))
  } else{
    return(z)
  }
}

orthogonalize <- function(X, group) {
  n <- nrow(X)
  J <- max(group)
  T <- vector("list", J)
  XX <- matrix(0, nrow=nrow(X), ncol=ncol(X))
  XX[,which(group==0)] <- X[,which(group==0)]
  for (j in seq_along(numeric(J))) {
    ind <- which(group==j)
    if (length(ind)==0) next
    SVD <- svd(X[, ind, drop=FALSE], nu=0)
    r <- which(SVD$d > 1e-10)
    T[[j]] <- sweep(SVD$v[,r,drop=FALSE], 2, sqrt(n)/SVD$d[r], "*")
    XX[,ind[r]] <- X[,ind]%*%T[[j]]
  }
  nz <- !apply(XX==0,2,all)
  XX <- XX[, nz, drop=FALSE]
  attr(XX, "T") <- T
  attr(XX, "group") <- group[nz]
  XX
}


unorthogonalize <- function(b, XX, group, intercept=TRUE) {
  require(Matrix)
  ind <- !sapply(attr(XX, "T"), is.null)
  T <- bdiag(attr(XX, "T")[ind])
  if (intercept) {
    ind0 <- c(1, 1+which(group==0))
    val <- Matrix::as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else if (sum(group==0)) {
    ind0 <- which(group==0)
    val <- Matrix::as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
  } else {
    val <- as.matrix(T %*% b)
  }
}


standardize <- function(X) {
  # Get dimensions
  n <- nrow(X)
  p <- ncol(X)
  
  # Initialize result matrices and vectors
  XX <- matrix(0, n, p)
  c <- numeric(p)
  s <- numeric(p)
  
  # Convert to matrices
  X <- as.matrix(X)
  XX <- as.matrix(XX)
  c <- as.vector(c)
  s <- as.vector(s)
  
  for (j in seq_len(p)) {
    # Center
    c[j] <- mean(X[, j])
    XX[, j] <- X[, j] - c[j]
    
    # Scale
    s[j] <- sqrt(mean(XX[, j]^2))
    XX[, j] <- XX[, j] / s[j]
  }
  
  # Return list
  res <- list(XX = XX, c = c, s = s)
  return(res)
}


maxgrad <- function(X, y, K, m) {
  n <- nrow(X)
  J <- length(K) - 1
  zmax <- 0
  
  for (g in 1:J) {
    Kg <- K[g+1] - K[g]
    Z <- numeric(Kg)
    for (j in (K[g] + 1):(K[g + 1])) {
      Z[j - K[g]] <- sum(X[, j] * y)
    }
    z <- sqrt(sum(Z^2))/m[g]
    zmax <- max(zmax, z)
  }
  return(zmax)
}

loglik_new <- function(params, wS, wO, xO, yO) {
  pS = ncol(wS)
  pO = ncol(xO)
  nO = nrow(wO)
  nS = nrow(wS)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS + pO)]
  #  sigma = -log(params[pS + pO + 1])
  #  rho = atanh(params[pS + pO + 2])
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  
  res0 = -1 * (wS %*% alpha)
  res1 = wO %*% alpha
  res2 = exp(sigma)*yO - (xO %*% beta)
  res3 = cosh(rho)*res1 + sinh(rho)*res2
  
  res0 = pnorm(res0,log.p=TRUE)
  res3 = pnorm(res3,log.p=TRUE)
  
  
  out = sum(res0) + sum(res3) - sum(res2*res2) / 2 + nO * (sigma - log(2*pi)/2)
  return(out)
}



############################################################################
# grBAR selection begins from here
############################################################################


ridge_ssel <- function(W, X, s, y, group_sel, group_out,
                       nlambda=100,lambda,lambda.min=1e-3,log.lambda=TRUE,
                       eps=1e-4,max_iter=10000,
                       init_strat = "MLE",threshold_sel=0.05,threshold_out=0.05){
  require(sampleSelection)
  group <- c(0,group_sel,0,group_out,0,0)
  pS <- ncol(W) + 1
  pO <- ncol(X) + 1
  
  if(init_strat == "MLE"){
    if(is.null(colnames(W))) colnames(W) <- paste("W",1:ncol(W),sep='')
    if(is.null(colnames(X))) colnames(X) <- paste("X",1:ncol(X),sep='')
    joint_dat <- data.frame(cbind(W,X))
    s_form <- as.formula(paste("s ~",paste(colnames(W),collapse=" + ")))
    y_form <- as.formula(paste("y ~",paste(colnames(X),collapse=" + ")))
    nongroup_fit <- selection(s_form, y_form, data = joint_dat)
    est_params <- nongroup_fit$estimate
    est_params[(pS+1):(pS+pO)] <- est_params[(pS+1):(pS+pO)]/est_params[pS+pO+1]
    est_params[pS+pO+1] = -log(est_params[pS+pO+1])
    est_params[pS+pO+2] = atanh(est_params[pS+pO+2])
    init <- est_params
  }
  else{
    init <- rep(0, length(group))
  }
  
  wS <- cbind(1,W[s==0,])
  wO <- cbind(1,W[s==1,])
  xO <- cbind(1,X[s==1,])
  yO <- y[s==1]
  
  twS <- t(wS)
  twO <- t(wO)
  txO <- t(xO)
  if(missing(lambda)){
    hess <- -1*hesslik(est_params,wS,wO,xO,yO,twS,twO,txO)
    hesseig <- eigen(hess)
    pseudoX <- hesseig$vectors %*% diag((hesseig$values) ** (1/2)) %*% t(hesseig$vectors)
    pseudoY <- pseudoX %*% est_params
    
    fity <- glm(pseudoY~pseudoX[, group==0]-1 , family="gaussian")
    r <- fity$residuals
    
    lambda.max <- 0
    
    for(j in 1:max(group)){
      temp_vec <- pseudoX[group==j,]%*%r
      groupgrad <- sum(solve(hess[group==j,group==j], temp_vec)*temp_vec)/sum(group==j)
      
      if(groupgrad > lambda.max){
        lambda.max = groupgrad
      }
    }
    
    if(log.lambda==TRUE){
      lambda = exp(seq(log(lambda.max*lambda.min), log(lambda.max), length.out=nlambda))
    } else{
      lambda = seq(lambda.max*lambda.min, lambda.max, length.out=nlambda)
    }
  }
  penalty_weights <- ifelse(group == 0, 0, 1)
  new_init <- init
  iter_list <- rep(0,length(lambda))
  out_params <- matrix(0,nrow=length(init),ncol=length(lambda))
  

  for(i in 1:length(lambda)){
    
    current_lambda = lambda[i]
    
    out_list <- tryCatch({ridgefit(wS, wO, xO, yO, current_lambda, eps, max_iter,
                         new_init, penalty_weights, twS, twO, txO)}, error=function(e){return(NA)})
    
    if(is.na(out_list)){
      out_params[,i] <- rep(NA, length(init))
    }
    
    if(out_list$converged==TRUE){
      new_init <- out_list$params
    }
    
    out_params[,i] <- out_list$params
    iter_list[i] <- out_list$iter
  }
  
  loss <- rep(NA,ncol(out_params))
  bic <- rep(NA,ncol(out_params))
  df <- rep(NA,ncol(out_params))
  
  valid_idx <- which(!is.na(out_params[1,]))
  
  # In the case all of the iterations ran into errors
  if(length(valid_idx)==0){
    out <- list(params=out_params,
                lambda=lambda,
                df=df,
                loss=loss,
                bic=bic,
                n=nrow(X),
                iter=iter_list,
                opt_params=rep(NA,length(init)),
                opt_lambda=NA,
                opt_bic=NA)
    return(out)
  }
  
  valid_params <- out_params[,valid_idx]
  
  
  sel_idx <- (1:pS)[group[1:pS]!=0]
  out_idx <- ((pS+1):(pS+pO))[group[((pS+1):(pS+pO))]!=0]
  
  out_params[sel_idx,] = ifelse(abs(out_params[sel_idx,])>threshold_sel,out_params[sel_idx,],0)
  out_params[out_idx,] = ifelse(abs(out_params[out_idx,]*exp(-out_params[pS+pO+1,]))>threshold_out,out_params[out_idx,],0)
  
  loss[valid_idx] <- apply(valid_params,2,loglik_new,
                           wS=cbind(1,W[s==0,]),
                           wO=cbind(1,W[s==1,]),
                           xO=cbind(1,X[s==1,]),
                           yO=y[s==1])
  
  df[valid_idx] <- apply(valid_params,2,function(x){return(sum(x!=0))}) - 4
  bic[valid_idx] <- df[valid_idx]*log(nrow(X)) - 2*loss[valid_idx]
  valid_params[pS+pO+2,] = tanh(valid_params[pS+pO+2,])
  valid_params[pS+pO+1,] = exp(-valid_params[pS+pO+1,])
  valid_params[(pS+1):(pS+pO),] = valid_params[(pS+1):(pS+pO),]*valid_params[pS+pO+1,]
  
  out_params[,valid_idx] <- valid_params
  
  best_idx <- which(bic == min(bic,na.rm=TRUE))
  opt_params <- out_params[,best_idx]
  opt_ridge_lambda <- ridge_lambda_extend[best_idx]
  opt_lambda <- lambda_extend[best_idx]
  opt_bic <- min(bic)
  
  out <- list(params=out_params,
              lambda=lambda,
              df=df,
              loss=loss,
              bic=bic,
              n=nrow(X),
              iter=iter_list,
              opt_params=opt_params,
              opt_lambda=opt_lambda,
              opt_bic=opt_bic)
  return(out)
}

##########################


##########################################################################
# Function for group variable selection in FIML estimator of Heckman model.  
# Group Broken Adaptive Ridge (BAR) Regression method has been implemented. 
# The code can also be used for individual penalization where each variable
# is taken to belong to its own group. 
###########################################################################


#' Group variable selection in nonignorable missing data using Broken Adaptive Ridge penalty.
#'
#' @param W A matrix of covariates in selection equation (intercept is not included).
#' @param X A matrix of covariates in outcome equation (intercept is not included).
#' @param s Binary outcome for the selection equation (0/1 and false/true allowed).
#' @param y Continuous outcome for the outcome equation (0 for NA).
#' @param group_sel A vector describing the grouping of the coefficients in the selection equation. 
#' It is best if group is a vector of consecutive integers. If there are coefficients to be included 
#' in the model without being penalized, assign them to group 0.
#' @param group_out A vector describing the grouping of the coefficients in the outcome equation. 
#' It is best if group is a vector of consecutive integers. If there are coefficients to be included 
#' in the model without being penalized, assign them to group 0.
#' @param ridge_lambda A user supplied sequence of lambda values for the initial 
#' ridge estimate. Typically, this is left unspecified, and the function automatically 
#' computes a grid of lambda values.
#' @param lambda A user supplied sequence of lambda values for 
#' the BAR iteration. Typically, this is left unspecified, and the function automatically 
#' computes a grid of lambda values.

#' @param nridge_lambda The number of lambda values for the ridge. Default is 10.
#' @param nlambda The number of lambda values for the BAR iteration. Default is 50
#' which leads to 500 combinations in the grid to test.

#' @param ridge_lambda.min The smallest value of ridge_lambda. The default is 0.001.
#' @param lambda.min The smallest value for lambda in BAR iteration. The default is 0.001.
#' @param log.ridge_lambda When TRUE compute the grid values of ridge_lambda on 
#' log scale (default) or linear scale otherwise.
#' @param log.lambda When TRUE compute the grid values of lambda on 
#' log scale (default) or linear scale otherwise.
#' @param ridge_eps The value fo tolerance used in intial ridge estimation.
#' @param eps The tolerance between BAR iteration.
#' @param inner_eps The value of eps for any optimization within each BAR iteration.
#' @param max_iter Maximum number of iterations (total across entire path). Default is 10000.
#' @param init_strat Use MLE to initialize the group descent algorithm otherwise use zeros.
#' @param del Use to avoid numerical overflow. Default is 0.01. See the paper for details.

#' @param del threshold_sel The value for which to set parameter values with magnitude
#' less than it to 0. Default is 0.05.
#' @param threshold_out Same as threshold_sel.
#' @param method Three methods are implemented. Method 2 is the fastest (see details
#' in accompanying paper).

#' @return class grHeckSelect_bar containing optimal penalized coefficients, lambda values,
#' bic values etc.
#' @export
#'
#' @examples
#' #grHeckSelect_bar(W,X,s,y,group_sel, group_out, method=2)

#'


grHeckSelect_bar <- function(W, X, s, y, group_sel, group_out,
                     ridge_lambda,lambda,
                     nridge_lambda=10,nlambda=50,
                     ridge_lambda.min=1e-3,lambda.min=1e-3,
                     log.ridge_lambda=TRUE,log.lambda=TRUE,
                     ridge_eps=1e-4,eps=1e-4,inner_eps=1e-4,max_iter=10000,
                     init_strat = "MLE", del=1e-2, threshold_sel = 0.05,
                     threshold_out=0.05, method=1){
  
  require(sampleSelection)
  group <- c(0,group_sel,0,group_out,0,0)
  pS <- ncol(W) + 1
  pO <- ncol(X) + 1
  
  if(init_strat == "MLE"){
    if(is.null(colnames(W))) colnames(W) <- paste("W",1:ncol(W),sep='')
    if(is.null(colnames(X))) colnames(X) <- paste("X",1:ncol(X),sep='')
    joint_dat <- data.frame(cbind(W,X))
    s_form <- as.formula(paste("s ~",paste(colnames(W),collapse=" + ")))
    y_form <- as.formula(paste("y ~",paste(colnames(X),collapse=" + ")))
    nongroup_fit <- selection(s_form, y_form, data = joint_dat)
    est_params <- nongroup_fit$estimate
    est_params[(pS+1):(pS+pO)] <- est_params[(pS+1):(pS+pO)]/est_params[pS+pO+1]
    est_params[pS+pO+1] = -log(est_params[pS+pO+1])
    est_params[pS+pO+2] = atanh(est_params[pS+pO+2])
    init <- est_params
  }
  else{
    init <- rep(0, length(group))
  }
  
  wS <- cbind(1,W[s==0,,drop=FALSE])
  wO <- cbind(1,W[s==1,,drop=FALSE])
  xO <- cbind(1,X[s==1,,drop=FALSE])
  yO <- y[s==1]
  twS <- t(wS)
  twO <- t(wO)
  txO <- t(xO)
  
  if(missing(ridge_lambda)){
    if(log.ridge_lambda==TRUE){
      ridge_lambda = exp(seq(log(10*ridge_lambda.min), log(10), length.out=nridge_lambda))
    } else{
      ridge_lambda = seq(10*ridge_lambda.min, 10, length.out=nridge_lambda)
    }
  }
  
  if(missing(lambda)){
    hess <- -1*hesslik(est_params,wS,wO,xO,yO,twS,twO,txO)
    hesseig <- eigen(hess)
    pseudoX <- hesseig$vectors %*% diag((hesseig$values) ** (1/2)) %*% t(hesseig$vectors)
    pseudoY <- pseudoX %*% est_params
    
    fity <- glm(pseudoY~pseudoX[, group==0]-1 , family="gaussian")
    r <- fity$residuals
    
    lambda.max <- 0
    
    for(j in 1:max(group)){
      temp_vec <- pseudoX[group==j,]%*%r
      groupgrad <- sum(solve(hess[group==j,group==j], temp_vec)*temp_vec)/sum(group==j)
      
      if(groupgrad > lambda.max){
        lambda.max = groupgrad
      }
    }
    
    if(log.lambda==TRUE){
      lambda = exp(seq(log(lambda.max*lambda.min), log(lambda.max), length.out=nlambda))
    } else{
      lambda = seq(lambda.max*lambda.min, lambda.max, length.out=nlambda)
    }
  }
  
  llam <- length(lambda)
  lpsi <- length(ridge_lambda)
  ridge_lambda_extend <- rep(ridge_lambda, each=llam)
  lambda_extend <- rep(lambda, lpsi)
  out_params <- matrix(0,nrow=length(init),ncol=lpsi*llam)
  iter_list <- rep(0,lpsi*llam)
  
  
  for(i in 1:(lpsi*llam)){
    out_list <- barfit(wS, wO, xO, yO, group, ridge_lambda_extend[i], lambda_extend[i],
                         ridge_eps, eps, inner_eps, max_iter, del, init, method, twS, twO, txO)
      
    out_params[,i] <- out_list$params
    iter_list[i] <- out_list$iter
  }
  
  loss <- rep(NA,ncol(out_params))
  bic <- rep(NA,ncol(out_params))
  df <- rep(NA,ncol(out_params))
  
  valid_idx <- which(!is.na(out_params[1,]))
  
  # In the case all of the iterations ran into errors
  if(length(valid_idx)==0){
    out <- list(params=out_params,
                lambda=lambda_extend,
                ridge_lambda=ridge_lambda_extend,
                df=df,
                loss=loss,
                bic=bic,
                n=nrow(X),
                iter=iter_list,
                opt_params=rep(NA,pS+pO+2),
                opt_lambda=NA,
                opt_ridge_lambda=NA,
                opt_bic=NA,
                nvalid=0)
    return(out)
  }
  
  valid_params <- out_params[,valid_idx]
  
  for(j in 1:max(group)){
    group_idx <- which(group==j)
    if(max(group_idx) < pS+1){
      group_threshold <- threshold_sel
    } else{
      group_threshold <- threshold_out
    }
    if(length(group_idx)>1){
      group_norm <- apply(valid_params[group_idx,],2,function(x){return(sqrt(sum(x**2)))})
    }
    else{
      group_norm <- abs(valid_params[group_idx,])
    }
    
    zero_set <- which(group_norm < group_threshold)
    valid_params[group_idx,zero_set] = 0
  }
  # valid_params[sel_idx,] = ifelse(abs(valid_params[sel_idx,])>threshold_sel,valid_params[sel_idx,],0)
  # valid_params[out_idx,] = ifelse(abs(valid_params[out_idx,]*exp(-valid_params[pS+pO+1,]))>threshold_out,
  #                               valid_params[out_idx,],0)
  

  
  
  loss[valid_idx] <- apply(valid_params,2,loglik_new,
                wS=cbind(1,W[s==0,]),
                wO=cbind(1,W[s==1,]),
                xO=cbind(1,X[s==1,]),
                yO=y[s==1])
  
  df[valid_idx] <- apply(valid_params,2,function(x){return(sum(x!=0))}) - 4
  bic[valid_idx] <- df[valid_idx]*log(nrow(X)) - 2*loss[valid_idx]
  valid_params[pS+pO+2,] = tanh(valid_params[pS+pO+2,])
  valid_params[pS+pO+1,] = exp(-valid_params[pS+pO+1,])
  valid_params[(pS+1):(pS+pO),] = valid_params[(pS+1):(pS+pO),]*valid_params[pS+pO+1,]
  
  out_params[,valid_idx] <- valid_params
  
  opt_bic <- min(bic,na.rm=TRUE)
  best_idx <- which(bic == opt_bic)
  opt_params <- out_params[,best_idx]
  opt_ridge_lambda <- ridge_lambda_extend[best_idx]
  opt_lambda <- lambda_extend[best_idx]
  
  rename <- c("(Intercept)", 
            paste0("S:", colnames(W)), 
            "(Intercept)", 
            paste0("O:", colnames(X)), 
            "sigma", "rho")
   names(opt_params) <- rename

  out <- list(params=out_params,
              lambda=lambda_extend,
              ridge_lambda=ridge_lambda_extend,
              df=df,
              loss=loss,
              bic=bic,
              n=nrow(X),
              iter=iter_list,
              opt_params=opt_params,
              opt_lambda=opt_lambda,
              opt_ridge_lambda=opt_ridge_lambda,
              opt_bic=opt_bic,
              nvalid= length(valid_idx))
           class(out) <- "grHeckSelect_bar"
  
  return(out)
}



#########################################################################
# Coefficients from object of class grHeckSelect_bar
#########################################################################

#' Title Coefficients from objects returned by grHeckSelect_bar function
#'
#' @param object a fitted object of class inheriting from "grHeckSelect_bar"
#' @param ...
#'
#' @return vector of coefficients for the selection and outcome models
#' @export
#'
#' @examples
#' #coef(grHeckSelect_bar(W,X,s,y,group_sel, group_out, method=2))
#'


coef.grHeckSelect_bar<- function(object,...){
 x <- object$opt_params
return(x)
}




barfit <- function(wS, wO, xO, yO, group, ridge_lambda=1, lambda=1,
                   ridge_eps=1e-4, eps=1e-4, inner_eps=1e-4, max_iter=10000, del=1e-2,
                   init,method=1, twS, twO, txO){
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  penalty_weights <- ifelse(group == 0, 0, 1)
  out_list <- ridgefit(wS, wO, xO, yO, ridge_lambda, ridge_eps, max_iter,
                       init, penalty_weights, twS, twO, txO)
  
  ridge_init <- out_list$params
  
  penalty_weights <- rep(0, length(group))
  iter <- 0
  old_params <- ridge_init
  while(iter < max_iter){
    
    iter <- iter+1
    for(j in 1:max(group)){
      penalty_weights[group==j] = 1/(sum(old_params[which(group==j)]**2)+del**2)
    }
    new_params <- tryCatch({
      if(method==1){
        barfit_1(wS, wO, xO, yO, lambda, inner_eps, max_iter,
                            old_params, penalty_weights, twS, twO, txO)
      } else if(method==2){
        barfit_2(wS, wO, xO, yO, lambda, max_iter,
                             old_params, penalty_weights, twS, twO, txO)
      } else if(method==3){
        barfit_3(wS, wO, xO, yO, lambda, inner_eps, max_iter,
                             old_params, penalty_weights, twS, twO, txO)
      }
    }, error=function(e){return(NA)})
    
    if(any(is.na(new_params))){
      return(list(params=rep(NA,length(ridge_init)),iter=NA))
    }
    
    if(max(abs(old_params - new_params)) < eps){
      break
    }
    
    old_params <- new_params
    
  }
  return(list(params=new_params, iter=iter))
}

barfit_1 <- function(wS, wO, xO, yO, lambda, eps, max_iter,
                     old_params, penalty_weights, twS, twO, txO){
  out_list <- ridgefit(wS, wO, xO, yO, lambda, eps, max_iter,
                       old_params, penalty_weights, twS, twO, txO)
  return(out_list$params)
}

barfit_2 <- function(wS, wO, xO, yO, lambda, max_iter,
                     old_params, penalty_weights, twS, twO, txO){
  D <- diag(penalty_weights)
  H <- -1*hesslik(old_params, wS, wO, xO, yO, twS, twO, txO)
  g <- -1*gradlik(old_params, wS, wO, xO, yO, twS, twO, txO)
  new_params <- solve(H + lambda*D, H%*%old_params - g)
  return(new_params)
}

barfit_3 <- function(wS, wO, xO, yO, lambda, eps, max_iter,
                     old_params, penalty_weights, twS, twO, txO){
  pS <- ncol(wO)
  pO <- ncol(xO)
  idx <- c(1,pS+1,pS+pO+1,pS+pO+2)
  new_params <- old_params
  
  D <- diag(penalty_weights)[-idx,-idx]
  H <- -1*hesslik_3(old_params, wS, wO, xO, yO, twS, twO, txO)[-idx,-idx]
  g <- -1*gradlik_3(old_params, wS, wO, xO, yO, twS, twO, txO)[-idx]
  
  new_params[-idx] <- solve(H + lambda*D, H%*%(old_params[-idx]) - g)
  out_list <- ridgefit_2(wS, wO, xO, yO, eps, max_iter,
                         new_params, twS, twO, txO, idx)
  
  return(out_list$params)
}

# WIP of another approach
# barfit_4 <- function(wS, wO, xO, yO, group, lambda, eps, max_iter,
#                       old_params, penalty_weights, twS, twO, txO){
#   pS <- ncol(wO)
#   pO <- ncol(xO)
#   idx <- c(1,pS+1,pS+pO+1,pS+pO+2)
#   hess <- -1*hesslik(old_params,wS,wO,xO,yO,twS,twO,txO)
#   hesseig <- eigen(hess)
#   pseudoX <- hesseig$vectors %*% diag((hesseig$values) ** (1/2)) %*% t(hesseig$vectors)
#   pseudoY <- pseudoX %*% old_params
#   
#   scale <- apply(pseudoX, 2, function(x){sqrt(sum(x**2)/nrow(pseudoX))})
#   XX <- pseudoX %*% diag(1/scale)
#   
#   g_order <- order(group)
#   g_order_inv <- match(1:length(group),g_order)
#   
#   XX_ord <- XX[,g_order]
#   
#   g <- group[g_order]
#   n <- nrow(XX_ord)
#   
#   XX_orth <- orthogonalize(XX_ord,g)
#   order_params <- old_params*scale
#   order_params <- order_params[g_order]
#   order_params <- orthogonalize_init(order_params, XX_orth, g, intercept=FALSE)
#   new_params <- order_params
#   print(order_params)
#   for(j in 1:max(g)){
#     Xout <- XX_orth[,-c(g==j),drop=FALSE]
#     Xj <- XX_orth[,g==j,drop=FALSE]
#     b <- as.vector(t(Xj)%*%(pseudoY - Xout%*%order_params[-c(g==j)]))
#     c <- sum(b**2)
#     if(c < 2*n*lambda){
#       new_params[g==j] = 0
#     } else{
#       w = (c - 2*n*lambda + sqrt(c*(c - 4*n*lambda)))/2*(n**2)
#       new_params[g==j] = b * (w/(w*n + lambda))
#     }
#   }
#   print(new_params)
#   new_params <- unorthogonalize2(new_params, XX_orth, g, intercept=FALSE)[[5]]
#   print(new_params)
#   new_params <- new_params[g_order_inv]
#   print(new_params)
#   new_params <- new_params/scale
#   
#   print(new_params)
#   
#   out_list <- ridgefit_2(wS, wO, xO, yO, eps, max_iter,
#                          new_params, twS, twO, txO, idx)
#   return(out_list$params)
# }

##################

ridgefit <- function(wS, wO, xO, yO, lambda, eps, max_iter,
                     init, penalty_weights, twS, twO, txO){
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  iter <- 0
  current_pen <- lambda * penalty_weights
  converged <- FALSE
  params <- init
  old_loglik <- -1*loglik_new(params, wS, wO, xO, yO) + 0.5*sum(current_pen * params**2)
  while(iter < max_iter){
    test_change <- 0
    iter <- iter + 1
    old_params <- params
    grad <- -1*gradlik(old_params, wS, wO, xO, yO, twS, twO, txO) + current_pen*params
      
    if(sqrt(sum(grad**2)) < eps){
      converged = TRUE
      break
    }
      
    hess <- -1*hesslik(old_params, wS, wO, xO, yO, twS, twO, txO) + diag(current_pen)
      
    change <- solve(hess,grad)
    params <- old_params - change
    new_loglik <- -1*loglik_new(params, wS, wO, xO, yO) + 0.5*sum(current_pen * params**2)
    while(old_loglik < new_loglik){
      test_change = test_change+1
      change <- change/2
      params <- old_params - change
      new_loglik <- -1*loglik_new(params, wS, wO, xO, yO) + 0.5*sum(current_pen * params**2)
    }
      
    if(abs(old_loglik - new_loglik) < eps){
      converged = TRUE
      break
    }
    
    old_loglik <- new_loglik
  }
  
  if(converged==FALSE){
    print("nonconvergence")
  }

  return(list(params=params, iter=iter, converged=converged))
}

ridgefit_2 <- function(wS, wO, xO, yO, eps, max_iter,
                     init, twS, twO, txO, idx){
  
  p <- max(idx)
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  iter <- 0
  converged <- FALSE
  params <- init
  old_loglik <- -1*loglik_new(params, wS, wO, xO, yO)
  while(iter < (max_iter-1)){
    test_change <- 0
    iter <- iter + 1
    old_params <- params
    grad <- -1*gradlik_2(old_params, wS, wO, xO, yO, twS, twO, txO)[idx]
    if(sqrt(sum(grad**2)) < eps){
      converged = TRUE
      break
    }
    
    hess <- -1*hesslik_2(old_params, wS, wO, xO, yO, twS, twO, txO)[idx,idx]
    change <- rep(0,p)
    change[idx] <- solve(hess,grad)
    params <- old_params - change
    new_loglik <- -1*loglik_new(params, wS, wO, xO, yO)
    while(old_loglik < new_loglik){
      test_change = test_change+1
      change <- change/2
      params <- old_params - change
      new_loglik <- -1*loglik_new(params, wS, wO, xO, yO)
    }
    
    if(abs(old_loglik - new_loglik) < eps){
      converged = TRUE
      break
    }
    
    old_loglik <- new_loglik
  }
  
  if(converged==FALSE){
    print("nonconvergence")
  }
  
  return(list(params=params, iter=iter, converged=converged))
}

hesslik <- function(params,
                    wS,wO,xO,yO,
                    twS,twO,txO) {
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  pS = ncol(wS)
  pO = ncol(xO)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  nOO = nrow(wO)
  nS = nrow(wS)
  
  sinhrho = sinh(rho)
  coshrho = cosh(rho)
  expsigma = exp(sigma)
  
  res0 = as.vector(-1 * wS %*% alpha)
  res1 = as.vector(wO %*% alpha)
  res2 = as.vector(expsigma * yO - xO %*% beta)
  res3 = as.vector(coshrho * res1 + sinhrho * res2)
  
  res0_imr = res0
  res0_imr = exp(dnorm(res0_imr, log=TRUE) - pnorm(res0_imr, log.p=TRUE))
  res0_imr = -1 * res0_imr * (res0 + res0_imr)
  
  
  v2 = as.vector(res3)
  v2 = exp(dnorm(v2, log=TRUE) - pnorm(v2, log.p=TRUE))
  
  v3 = -1 * v2 * (v2 + res3)
  
  hess = matrix(0, nrow = pS + pO + 2, ncol = pS + pO + 2)
  
  v4 = v3 * sinhrho**2
  v5 = sinhrho * res1 + coshrho * res2
  v6 = v3 * v5
  m1 = wO * v3
  
  hess[1:pS,1:pS] = coshrho**2 * (twO %*% m1) + twS %*% (wS * res0_imr)
  hess[(pS+1):(pS+pO),1:pS] = -1 * sinh(2 * rho) / 2 * (txO %*% m1)
  hess[(pS+1):(pS+pO),(pS+1):(pS+pO)] = txO %*% (xO * (v4 - 1))
  hess[pS+pO+1,1:pS] = sinh(2 * rho) / 2 * expsigma * (twO %*% (v3 * yO))
  hess[pS+pO+1,(pS+1):(pS+pO)] = expsigma * txO %*% (yO * (1 - v4))
  hess[pS+pO+1, pS+pO+1] = expsigma * sum(yO * (sinhrho * v2 + xO %*% beta) + expsigma * (yO**2) * (v4 - 2))
  hess[pS+pO+2,1:pS] = twO %*% (coshrho * v6 + sinhrho * v2)
  hess[pS+pO+2,(pS+1):(pS+pO)] = -1 * txO %*% (sinhrho * v6 + coshrho * v2)
  hess[pS+pO+2, pS+pO+1] = expsigma * sum(yO * (sinhrho*v6 + coshrho*v2))
  hess[pS+pO+2, pS+pO+2] = sum(v3 * v5**2 + res3 * v2)
  upper <- upper.tri(hess)
  hess[upper] <- t(hess)[upper]
  return(hess)
}


gradlik <- function(params,wS,wO,xO,yO,twS,twO,txO){
  pS = ncol(wS)
  pO = ncol(xO)
  nO = nrow(wO)
  nS = nrow(wS)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  expsigma = exp(sigma)
  coshrho = cosh(rho)
  sinhrho = sinh(rho)
  
  res0 = -1 * wS %*% alpha
  res1 = wO %*% alpha
  res2 = expsigma*yO - xO %*% beta
  res3 = coshrho * res1 + sinhrho * res2
  
  v4 = sinhrho * res1 + coshrho * res2
  
  
  res0 = exp(dnorm(res0, log=TRUE) - pnorm(res0, log.p=TRUE))
  
  res3 = exp(dnorm(res3, log=TRUE) - pnorm(res3, log.p=TRUE))
  
  v5 = res2 - sinhrho * res3
  
  out <- rep(0,pS+pO+2)
  
  out[1:pS] = coshrho * (twO %*% res3) - twS %*% res0;
  out[(pS+1):(pS+pO)] = txO %*% v5
  out[pS+pO+1] = nO - expsigma * sum(yO * v5)
  out[pS+pO+2] = sum(res3 * v4)
  
  
  return(out)
}


#####################
# Helper functions specifically for method=3

gradlik_2 <- function(params,wS,wO,xO,yO,twS,twO,txO){
  pS = ncol(wS)
  pO = ncol(xO)
  nO = nrow(wO)
  nS = nrow(wS)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  expsigma = exp(sigma)
  coshrho = cosh(rho)
  sinhrho = sinh(rho)
  w1S <- wS[,1]
  w1O <- wO[,1]
  x1O <- xO[,1]
  
  
  res0 = -1 * wS %*% alpha
  res1 = wO %*% alpha
  res2 = expsigma*yO - xO %*% beta
  res3 = coshrho * res1 + sinhrho * res2
  
  v4 = sinhrho * res1 + coshrho * res2
  
  
  res0 = exp(dnorm(res0, log=TRUE) - pnorm(res0, log.p=TRUE))
  
  res3 = exp(dnorm(res3, log=TRUE) - pnorm(res3, log.p=TRUE))
  
  v5 = res2 - sinhrho * res3
  
  out <- rep(0,pS+pO+2)
  
  out[1] <- coshrho * sum(w1O * res3) - sum(w1S * res0)
  out[pS+1] <- sum(x1O * v5)
  out[pS+pO+1] = nO - expsigma * sum(yO * v5)
  out[pS+pO+2] = sum(res3 * v4)
  
  
  return(out)
}

hesslik_2 <- function(params,
                      wS,wO,xO,yO,
                      twS,twO,txO) {
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  pS = ncol(wS)
  pO = ncol(xO)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  nOO = nrow(wO)
  nS = nrow(wS)
  
  w1S <- wS[,1]
  w1O <- wO[,1]
  x1O <- xO[,1]
  
  sinhrho = sinh(rho)
  coshrho = cosh(rho)
  expsigma = exp(sigma)
  
  res0 = as.vector(-1 * wS %*% alpha)
  res1 = as.vector(wO %*% alpha)
  res2 = as.vector(expsigma * yO - xO %*% beta)
  res3 = as.vector(coshrho * res1 + sinhrho * res2)
  
  res0_imr = res0
  res0_imr = exp(dnorm(res0_imr, log=TRUE) - pnorm(res0_imr, log.p=TRUE))
  res0_imr = -1 * res0_imr * (res0 + res0_imr)
  
  
  v2 = as.vector(res3)
  v2 = exp(dnorm(v2, log=TRUE) - pnorm(v2, log.p=TRUE))
  
  v3 = -1 * v2 * (v2 + res3)
  
  hess = matrix(0, nrow = pS + pO + 2, ncol = pS + pO + 2)
  
  v4 = v3 * sinhrho**2
  v5 = sinhrho * res1 + coshrho * res2
  v6 = v3 * v5
  
  twOwO <- twO %*% wO
  
  hess[1,1] = coshrho**2 * sum(w1O**2 * v3) + sum(w1S**2 * res0_imr)
  hess[pS+1,1] = -1/2 * sinh(2*rho) * sum(x1O * w1O * v3)
  hess[1,pS+1] = hess[pS+1,1]
  hess[pS+1,pS+1] = sum(x1O**2 * (v4 - 1))
  hess[pS+pO+1,1] = 1/2 * sinh(2*rho) * expsigma * sum(w1O * v3 * yO)
  hess[1,pS+pO+1] = hess[pS+pO+1,1]
  hess[pS+pO+1,pS+1] = expsigma * sum(x1O * yO * (1 - v4))
  hess[pS+1,pS+pO+1] = hess[pS+pO+1,pS+1]
  hess[pS+pO+1, pS+pO+1] = expsigma * sum(yO * (sinhrho * v2 + xO %*% beta) + expsigma * (yO**2) * (v4 - 2))
  hess[pS+pO+2,1] = sum(w1O * (coshrho * v6 + sinhrho * v2))
  hess[1,pS+pO+2] = hess[pS+pO+2,1]
  hess[pS+pO+2,pS+1] = -1 * sum(x1O * (sinhrho * v6 + coshrho * v2))
  hess[pS+1,pS+pO+2] = hess[pS+pO+2,pS+1]
  hess[pS+pO+2, pS+pO+1] = expsigma * sum(yO * (sinhrho*v6 + coshrho*v2))
  hess[pS+pO+1,pS+pO+2] = hess[pS+pO+2,pS+pO+1]
  hess[pS+pO+2, pS+pO+2] = sum(v3 * v5**2 + res3 * v2)
  return(hess)
}


gradlik_3 <- function(params,wS,wO,xO,yO,twS,twO,txO){
  pS = ncol(wS)
  pO = ncol(xO)
  nO = nrow(wO)
  nS = nrow(wS)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  expsigma = exp(sigma)
  coshrho = cosh(rho)
  sinhrho = sinh(rho)
  
  res0 = -1 * wS %*% alpha
  res1 = wO %*% alpha
  res2 = expsigma*yO - xO %*% beta
  res3 = coshrho * res1 + sinhrho * res2
  
  
  
  res0 = exp(dnorm(res0, log=TRUE) - pnorm(res0, log.p=TRUE))
  
  res3 = exp(dnorm(res3, log=TRUE) - pnorm(res3, log.p=TRUE))
  
  v5 = res2 - sinhrho * res3
  
  out <- rep(0,pS+pO+2)
  
  out[1:pS] = coshrho * (twO %*% res3) - twS %*% res0;
  out[(pS+1):(pS+pO)] = txO %*% v5
  
  
  return(out)
}

hesslik_3 <- function(params,
                      wS,wO,xO,yO,
                      twS,twO,txO) {
  
  if(missing(twS)){
    twS <- t(wS)
  }
  if(missing(twO)){
    twO <- t(wO)
  }
  if(missing(txO)){
    txO <- t(xO)
  }
  
  pS = ncol(wS)
  pO = ncol(xO)
  alpha = params[1:pS]
  beta = params[(pS+1):(pS+pO)]
  sigma = params[pS+pO+1]
  rho = params[pS+pO+2]
  nOO = nrow(wO)
  nS = nrow(wS)
  
  sinhrho = sinh(rho)
  coshrho = cosh(rho)
  expsigma = exp(sigma)
  
  res0 = as.vector(-1 * wS %*% alpha)
  res1 = as.vector(wO %*% alpha)
  res2 = as.vector(expsigma * yO - xO %*% beta)
  res3 = as.vector(coshrho * res1 + sinhrho * res2)
  
  res0_imr = res0
  res0_imr = exp(dnorm(res0_imr, log=TRUE) - pnorm(res0_imr, log.p=TRUE))
  res0_imr = -1 * res0_imr * (res0 + res0_imr)
  
  
  v2 = as.vector(res3)
  v2 = exp(dnorm(v2, log=TRUE) - pnorm(v2, log.p=TRUE))
  
  v3 = -1 * v2 * (v2 + res3)
  
  hess = matrix(0, nrow = pS + pO + 2, ncol = pS + pO + 2)
  
  v4 = v3 * sinhrho**2
  m1 = wO * v3
  
  hess[1:pS,1:pS] = coshrho**2 * (twO %*% m1) + twS %*% (wS * res0_imr)
  hess[(pS+1):(pS+pO),1:pS] = -1 * sinh(2 * rho) / 2 * (txO %*% m1)
  hess[(pS+1):(pS+pO),(pS+1):(pS+pO)] = txO %*% (xO * (v4 - 1))
  upper <- upper.tri(hess)
  hess[upper] <- t(hess)[upper]
  return(hess)
}

########################################

# more WIP, ignore for now
# unorthogonalize2 <- function(b, XX, group, intercept=TRUE) {
#   require(Matrix)
#   ind <- !sapply(attr(XX, "T"), is.null)
#   T <- bdiag(attr(XX, "T")[ind])
#   if (intercept) {
#     ind0 <- c(1, 1+which(group==0))
#     val <- Matrix::as.matrix(rbind(b[ind0, , drop=FALSE], T %*% b[-ind0, , drop=FALSE]))
#   } else if (sum(group==0)) {
#     ind0 <- which(group==0)
#     val <- c(b[ind0], T %*% b[-ind0])
#   } else {
#     val <- as.matrix(T %*% b)
#   }
# }





