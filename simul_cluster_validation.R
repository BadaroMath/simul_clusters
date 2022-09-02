library("mclust")
library("cluster")
library("MixSim")
library("sClust")
library("fpc")
library("dbscan")
library("ppclust")
library("fclust")


###################### KMEANS #################################################
###############################################################################


kmeans_simul <- function(A, K){
#
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#        iter {numeric} -- number of iterations
#        nstart {numeric} -- start point
#    Returns:
#        {list} -- Rand index, correct classification proportion and information variance
    start_time <- Sys.time()
    id <- kmeans(A$X, K)$cluster
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)
    end_time <- Sys.time()
    return(c(rand, prop, varinf, end_time - start_time))
}


###################### Model Based ############################################
###############################################################################
mclust_simul <- function(A, K, model, n){
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#        model {string} -- mcluster model
#    Returns:
#        {list} -- Rand index, correct classification
#        proportion and information variance
    start_time <- Sys.time()
    id <- try(Mclust(A$X, G = K, model = "VVV", verbose = FALSE)$class)

    if (inherits(id, what = "try-error") || (is.null(id))){
        id <- try(Mclust(A$X, G = K, verbose = FALSE)$class)
        if (inherits(id, what = "try-error") || (is.null(id))){
            id <- rep(1, n)
        }
    }
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)
    end_time <- Sys.time()
    return(c(rand, prop, varinf, end_time - start_time))
}

###################### Espectral ##############################################
###############################################################################
spectral <- function(X, nn = 10, n_eig = 2){
#    Arguments:
#        X {matrix} -- matrix of data points
#        nn {numeric} -- the k nearest neighbors to consider
#        n_eig {numeric} -- number of eignenvectors to keep
#    Returns:
#        {list} -- the eigenvectors of the n_eig smallest eigenvaluesr
#    Author:
#        Nura Kawa (February 24, 2018)
#    Reference:
#        Chapter 14.5.3 of the second edition of
#        Elements of Statistical Learning: Data Mining, Inference, and
#        Prediction by Trevor Hastie, Robert Tibshirani, and Jerome Friedman.

    mutual_knn_graph <- function(X, nn = 10)
    {
        D <- as.matrix(dist(X)) # matrix of euclidean distances
                                    # between data points in X
        # intialize the knn matrix
        knn_mat <- matrix(0,
                            nrow = nrow(X),
                            ncol = nrow(X))
        # find the 10 nearest neighbors for each point
        for (i in 1:nrow(X)) { # nolint
            neighbor_index <- order(D[i,])[2:(nn + 1)]
            knn_mat[i,][neighbor_index] <- 1
        }
        # Now we note that i,j are neighbors iff K[i,j] = 1 or K[j,i] = 1
        knn_mat <- knn_mat + t(knn_mat) # find mutual knn
        knn_mat[ knn_mat == 2 ] = 1
        return(knn_mat)
    }
    graph_laplacian <- function(W, normalized = TRUE)
    {
        stopifnot(nrow(W) == ncol(W))
        g <- colSums(W) # degrees of vertices
        n <- nrow(W)
        if(normalized)
        {
            D_half <- diag(1 / sqrt(g) )
            return( diag(n) - D_half %*% W %*% D_half )
        }
        else
        {
            return( diag(g) - W)
        }
    }
    W <- mutual_knn_graph(X) # 1. matrix of similarities
    L <- graph_laplacian(W) # 2. compute graph laplacian
    ei <- eigen(L, symmetric = TRUE) # 3. Compute the eigenvectors
                                    # and values of L
    n <- nrow(L)
    return(ei$vectors[,(n - n_eig):(n - 1)])
}

spectral_simul <- function(A, K, iter, nstart){
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#        iter {numeric} -- number of iterations
#        nstart {numeric} -- start point
#    Returns:
#        {list} -- Rand index, correct classification 
#        proportion and information variance
    start_time <- Sys.time()
    X_spec <- spectral(A$X)
    id <- kmeans(X_spec, K, iter.max = iter, nstart)$cluster
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)
    end_time <- Sys.time()
    return(c(rand, prop, varinf, end_time - start_time))
}

###################### DBSCAN #################################################
###############################################################################
dbscan <- function(A, minpts, eps, n){
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#        iter {numeric} -- number of iterations
#        nstart {numeric} -- start point
#    Returns:
#        {list} -- Rand index, correct classification 
#        proportion and information variance
    id <- dbscan::dbscan(A$X, eps, minpts)$cluster
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- ifelse(round(ClassProp(A$id, id), 3) > 1, round(ClassProp(A$id, id), 3)/n, round(ClassProp(A$id, id), 3))
    varinf <- round(VarInf(A$id, id), 3)

    return(c(rand, prop, varinf))
}

dbscan_simul <- function(A, n){
#    Arguments:
#        A {list} -- generated cluster
#    Returns:
#        {list} -- Rand index, correct classification 
#        proportion and information variance
    start_time <- Sys.time()
    rand_dbscan <- matrix(NA, nrow = 50*length(seq(0, 2, 0.01)), ncol = 3)
    iter <-0
    for (minpts in seq(1, 100, 5)){
      for (eps in seq(0, 0.5, .1)){
        rand_dbscan[iter, 1] <- eps
        rand_dbscan[iter, 2] <- minpts
        rand_dbscan[iter, 3] <- dbscan(A,minpts, eps, n)[1]
        iter <- iter + 1
      }
    }
    eps_opt <- rand_dbscan[
        which(rand_dbscan[,3] == max(na.omit(rand_dbscan[, 3])))[1], 1
        ]
    minpts_opt <- rand_dbscan[
        which(rand_dbscan[,3] == max(na.omit(rand_dbscan[, 3])))[1], 2
        ]
    result = dbscan(A, minpts_opt, eps_opt)
    end_time <- Sys.time()
    
    return(c(result, end_time - start_time))
}

###################### Fuzzy C-Means ##########################################
###############################################################################
fcm_simul <- function(A, K){
#
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#    Returns:
#        {numeric} -- Fuzzy Rand index, correct classification (non fuzzy labels)
#        proportion and information variance (non fuzzy labels)
    start_time <- Sys.time()   
    fc <- fcm(A$X, K)
    rand <- round(RI.F(A$id, fc$u, "minimum"), 3)
    prop <- round(ClassProp(A$id, fc$cluster), 3)
    varinf <- round(VarInf(A$id, fc$cluster), 3)
    end_time <- Sys.time()
    
    return(c(rand, prop, varinf, end_time - start_time))
}


###################### Simulation ##########################################
###############################################################################


operator <- function(BarOmega, K, p, n, file){
  Q <- MixSim(BarOmega = BarOmega, K = K, p = p, resN = 100000)
  A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
  print("Kmeans:")
  cat(toString(
    c(n, K, p, BarOmega, kmeans_simul(A, K), "K-Means")
  ),
  file=paste("clusters_index_", file, ".txt"), append = TRUE, sep="\n")
  print("done:")
  print("mclust:")
  cat(toString(
    c(n, K, p, BarOmega, mclust_simul(A, K, "VVV", n), "Model Based")
  ),
  file=paste("clusters_index_", file, ".txt"), append = TRUE, sep="\n")
  print("done:")
  print("spectral:")
  cat(toString(
    c(n, K, p, BarOmega, spectral_simul(A, K, 100, 10), "Spectral")
  ),
  file=paste("clusters_index_", file, ".txt"), append = TRUE, sep="\n")
  print("done:")
  print("dbscan:")
  cat(toString(
    c(n, K, p, BarOmega, dbscan_simul(A, n), "DBSCAN")
  ),
  file=paste("clusters_index_", file, ".txt"), append = TRUE, sep="\n")
  print("done:")
  print("fcm:")
  cat(toString(
    c(n, K, p, BarOmega, fcm_simul(A, K), "FCM")
  ),
  file=paste("clusters_index_", file, ".txt"), append = TRUE, sep="\n")
  print("done:")
  print(c(n, K, p, BarOmega))  
}



# 1000 <= n <= 500000, step = 100  -- Número de observações
# 5 <= K <= 10000, step = 5 -- Número de componente
# 1 <= p <= 10, step = 1 -- Número de dimensões
# 0 <= BarOmega <= 0.8, step = 0.01 -- Sobreposição média.

################ General (individual)  ##############


for (BarOmega in seq(0, 0.6, 0.01)){
  operator(BarOmega, 3, 5, 1000, "BarOmega")
}

for (K in seq(2, 12, 1)){
  operator(0, K, 5, 1000, "Clusters")
}

for (n in seq(100, 10000, 100)){
  operator(0, 3, 5, n, "N")
}

for (p in seq(60, 200, 10)){
  operator(0, 3, p, 1000, "Components")
}




################ General (individual) BarOmega = 0.05  ##############

for (K in seq(2, 13, 1)){
  operator(0.05, K, 5, 1000, "Clusters_Omega5")
}

for (n in seq(1000, 500000, 1000)){
  operator(0.05, 3, 5, n, "N_Omega5")
}

for (p in seq(0, 100, 1)){
  operator(0.05, 3, p, 1000, "Components_Omega5")
}

################ General (individual) BarOmega = 0.1  ##############

for (K in seq(2, 30, 1)){
  operator(0.1, K, 5, 1000, "Clusters_Omega10")
}

for (n in seq(1000, 500000, 1000)){
  operator(0.1, 3, 5, n, "N_Omega10")
}

for (p in seq(0, 100, 1)){
  operator(0.1, 3, p, 1000, "Components_Omega10")
}


################ General (individual) BarOmega = 0.15  ##############

for (K in seq(2, 30, 1)){
  operator(0.15, K, 5, 1000, "Clusters_Omega15")
}

for (n in seq(1000, 10000, 1000)){
  operator(0.15, 3, 5, n, "N_Omega15")
}

for (p in c(3, 5, 8, 10, 15, 30, 50)){
  operator(0.15, 3, p, 1000, "Components_Omega15")
}


################ General (individual) BarOmega = 0.20  ##############

for (K in seq(2, 30, 1)){
  operator(0.20, K, 5, 1000, "Clusters_Omega20")
}

for (n in seq(1000, 500000, 1000)){
  operator(0.20, 3, 5, n, "N_Omega20")
}

for (p in c(3, 5, 8, 10, 15, 30, 50)){
  operator(0.20, 3, p, 1000, "Components_Omega20")
}



################ Cluster & Components  ##############


for (K in c(3, 5, 7, 10)){
  for (p in c(3, 5, 8, 10, 15, 30, 50)){
    operator(0, K, p, 1000, "Clusters_Comp")
  }
}



################ Obs & Components  ##############


for (n in c(50, 100, 500, 1000, 3000, 5000)){
  for (p in c(3, 5, 8, 10, 15, 30, 50)){
    operator(0, K, p, 1000, "N_Comp")
  }
}
