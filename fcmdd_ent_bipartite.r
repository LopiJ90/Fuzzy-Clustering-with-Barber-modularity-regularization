fcmdd_ent_bipartite <- function(X,Y,B, k, ent, ent1,gamma,
                                RS = 100, max_iter = 1e+3, 
                                sort_medoids = FALSE, init = FALSE,
                                progressbar = FALSE) {
  dx <- as.matrix(X)^2
  dy <- as.matrix(Y)^2 
  n<- nrow(dx)
  m <- nrow(dy)
  values <- vector(length(RS), mode = "numeric")
  n_iter <- vector(length(RS), mode = "numeric")
  conv <- 1e-09
  func_opty <- 10 ^ 10 * sum(dy)
  func_optx <- 10 ^ 10 * sum(dx)
  func_opt=func_opty+func_optx
  # RS iteration to reduce the risk of local optimum
  if(progressbar) pb <- txtProgressBar(min = 0, max = RS, style = 3)
  for (rs in 1:RS) {
    if(progressbar) 
      Sys.sleep(0.1)
    set.seed(rs)
    # initialization of the membership degree matrix
    if (init) {
      U <- start_u(n = n, k = k)
    } else {
      U <- matrix(runif(n * k, 0, 1), nrow = n, ncol = k)
      U <- U / apply(U, 1, sum)
    }
    if (init) {
      V <- start_u(n = m, k = k)
    } else {
      V <- matrix(runif(m * k, 0, 1), nrow = m, ncol = k)
      V <- V / apply(V, 1, sum)
    }
    
    # start of the iterative algorithm
    Uold <- U + 1
    Vold <- V + 1
    iter <- 0
    
    medoidsx= sample(1:n, k, replace=F)
    medoidsy= sample(1:m, k, replace=F)
    
    while ((sum(abs(U - Uold)) > conv | sum(abs(V-Vold))>conv) && (iter < max_iter)) {
      iter <- iter + 1
      Uold <- U  
      Vold <- V
      # update medoids
      
      medoidsx <- medoids_update(X = X, U = U, distance = dx, 
                                 k = k, n = n,start=medoidsx)
      d_icx <- dx[, medoidsx, drop = FALSE]
      # update medoid
      medoidsy <- medoids_update(X = Y, U = V, distance = dy, 
                                 k = k, n = m,start=medoidsy)
      d_icy <- dy[, medoidsy, drop = FALSE]
      # update membership degrees
      U <- membership_updateU(U = U, d_ic = d_icx, medoids = medoidsx, B,
                              n = n, k = k,ent=ent, gamma, V)
      V<- membership_updateV(V=V, d_ic=d_icy, medoids = medoidsy, t(B),
                             m = m, k = k,ent=ent1, gamma, U)
    }
    for(i in 1:k){
      U[medoidsx[i],i]=1
      U[medoidsx[i],-i]=0
      V[medoidsy[i],i]=1
      V[medoidsy[i],-i]=0
    }
    logcx<-log(U) 
    logcx[which(logcx==-Inf)]<-0
    logcy<-log(V) 
    logcy[which(logcy==-Inf)]<-0
    func <- sum(U* d_icx)+sum(V* d_icy)+ent*sum(U*logcx)+ent1*sum(V*logcy)-gamma*sum(U* B%*%V)
    values[rs] <- func
    n_iter[rs] <- iter
    if (func < func_opt) {
      U_opt <- U
      medoids_opty <- medoidsy
      V_opt <- V
      medoids_optx <- medoidsx
      func_opt <- func
      d_ic_opty <- d_icy
      d_ic_optx <- d_icx
    }
    # update progress bar
    if(progressbar) setTxtProgressBar(pb, rs)
  }
  if(progressbar) close(pb)
  names_medoidsx <- medoids_optx
  names_medoidsy <- medoids_opty   
  names(values) <- paste("Random start", 1:RS, sep = " ")
  names(n_iter) <- names(values)
  membx <- apply(U_opt, 1, max)
  hardx <- apply(U_opt, 1, which.max)
  memby <- apply(V_opt, 1, max)
  hardy <- apply(V_opt, 1, which.max)
  clusteringx <- data.frame(memb = membx, hard = hardx)
  clusteringy <- data.frame(memb = memby, hard = hardy)
  # output
  results <- list()
  results$U <- U_opt
  results$V <- V_opt
  results$clusteringx <- clusteringx
  results$clusteringy <- clusteringy
  results$medoidsx <- medoids_optx
  results$medoidsy <- medoids_opty
  results$values <- values
  results$func_obj <- min(values)
  results$n_iter <- n_iter
  results$k <- k
  results$d_icx <- d_ic_optx
  results$d_icy <- d_ic_opty
  return(results)
}

# internal functions for 'fcmdd' function ---------------------------------

# medoids update ----------------------------------------------------------

medoids_update <- function(X, U, k, n, distance,start) {
  medoids <- start
  for (c in 1:k) {
    min_dist_i <- 10^5 * sum(X ^ 2, na.rm = TRUE)
    ind=apply(U,1, which.max)
    for (i in 1:n) {
      if(ind[i]==c){
        dist_i <- sum(U[, c]* distance[i, ])
        if ((dist_i < min_dist_i) 
            &            match(i, medoids[-c], nomatch = 0) == 0
        ){
          min_dist_i <- dist_i
          medoids[c] <- i
        }
      }
    }
  }
  return(medoids)
}

membership_updateU <- function(U, d_ic, medoids, n, k,ent,  B, gamma,V) {
  for (i in 1:n) {
    for (c in 1:k) {
      U[i, c] <-(exp(-(((ent/log(k))^-1)*d_ic[i, c])+(gamma*((ent/log(k))^-1))*sum(B[i,]*V[,c])))
    }   
    U[i,]=U[i,]/sum(U[i,])
  }
  return(U)
}

membership_updateV <- function(V, d_ic, medoids, m, k,ent,  B, gamma,U) {
  for (i in 1:m) {
    for (c in 1:k) {
      V[i, c] <-(exp(-(((ent/log(k))^-1)*d_ic[i, c])+(gamma*((ent/log(k))^-1))*sum(B[i,]*U[,c])))
    }
    
    V[i,]=V[i,]/sum(V[i,])
  }
  return(V)
}

# initial U ---------------------------------------------------------------

start_u <- function(n, k) {
  U_u <- matrix(1 / k, nrow = n, ncol = k)
  U_r <- matrix(0, nrow = n, ncol = k)
  pos_ones <- sample(1:k, n, replace = TRUE)
  U_r[cbind(1:n, pos_ones)] <- 1
  U <- (1 - sqrt(2) / 2) * U_u + sqrt(2) / 2 * U_r
  return(U)
}


xie_beni_fcmdd21 <- function(fit, B) {
  c=fit$k
  n <- nrow(fit$U)
  m<- nrow(fit$V)
  cx=(n-c)
  cy=(m-c)
  a=sum(fit$U* fit$d_icx)/cx+sum(fit$V*fit$d_icy)/cy
  d_icx <- fit$d_icx
  d_icy <- fit$d_icy
  medoidsx <- unname(fit$medoidsx)
  medoidsy <- unname(fit$medoidsy)
  dist_medoidsx <- as.dist(d_icx[medoidsx,])
  dist_medoidsy <- as.dist(d_icy[medoidsy,])
  bb=min(dist_medoidsx)+min(dist_medoidsy)/c
  J <- sum(fit$U* B%*%fit$V)/c
  index=(bb+J)/a
  return(index)
}
