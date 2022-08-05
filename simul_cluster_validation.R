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
kmeans_simul <- function(A, K, iter, nstart){
#
#    Arguments:
#        A {list} -- generated cluster
#        k {numeric} -- number of clusters
#        iter {numeric} -- number of iterations
#        nstart {numeric} -- start point
#    Returns:
#        {list} -- Rand index, correct classification proportion and information variance
#   
    id <- kmeans(A$X, K, iter.max = 1000, nstart = 10)$cluster
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)

    return(c(rand, prop, varinf))
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
    return(c(rand, prop, varinf))
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

    X_spec <- spectral(A$X)
    id <- kmeans(X_spec, K, iter.max = iter, nstart = nstart)$cluster
    rand <- round(RandIndex(A$id, id)$AR, 3)
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)

    return(c(rand, prop, varinf))
}

###################### DBSCAN #################################################
###############################################################################
dbscan <- function(A, minpts, eps){
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
    prop <- round(ClassProp(A$id, id), 3)
    varinf <- round(VarInf(A$id, id), 3)

    return(c(rand, prop, varinf))
}

dbscan_simul <- function(A, n){
#    Arguments:
#        A {list} -- generated cluster
#    Returns:
#        {list} -- Rand index, correct classification 
#        proportion and information variance
    rand_dbscan <- matrix(NA, nrow = 50*length(seq(0, 2, 0.01)), ncol = 3)
    iter <-0
    for (eps in seq(0, 2, 0.01)){
        for (minpts in seq(1, n, 2)){
            rand_dbscan[iter, 1] <- eps
            rand_dbscan[iter, 2] <- minpts
            rand_dbscan[iter, 3] <- dbscan(A,minpts, eps)[1]
            iter <- iter + 1
        }
    }
    eps_opt <- rand_dbscan[which(rand_dbscan[,3] == max(na.omit(rand_dbscan[, 3])))[1], 1]
    minpts_opt <- rand_dbscan[which(rand_dbscan[,3] == max(na.omit(rand_dbscan[, 3])))[1], 2]
    return(dbscan(A, minpts_opt, eps_opt))
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
#   
    fc <- fcm(A$X, K)
    rand <- round(RI.F(A$id, fc$u, "minimum"), 3)
    prop <- round(ClassProp(A$id, fc$cluster), 3)
    varinf <- round(VarInf(A$id, fc$cluster), 3)

    return(c(rand, prop, varinf))
}

###################### Simulation ##########################################
###############################################################################

# 1000 <= n <= 500000  -- number of points
# 3 <= K <= 50 -- number of components
# 1 <= p <= 10 -- number of dimensions
# 0 <= BarOmega <= 0.8 -- value of desired average overlap.

set.seed(123)
start_time <- Sys.time()
for (n in seq(100,500000, 100)){
    for (K in seq(5,10000, 5)){
        for (p in 1:10){
            for (BarOmega in seq(0, 0.90, 0.02)){
                tryCatch({
                    Q <- MixSim(BarOmega = BarOmega, K = K, p = p)
                    A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)
                    cat(toString(
                        c(n, K, p, BarOmega, kmeans_simul(A, K, 100, 10))
                    ),
                        file="kmeans_index.txt", append = TRUE, sep="\n")
                    cat(toString(
                        c(n, K, p, BarOmega, mclust_simul(A, K, "VVV", n))
                    ),
                        file="mclust_index.txt", append = TRUE, sep="\n")
                    cat(toString(
                        c(n, K, p, BarOmega, spectral_simul(A, K, 100, 10))
                    ),
                        file="spec_index.txt", append = TRUE, sep="\n")
                    cat(toString(
                        c(n, K, p, BarOmega, dbscan_simul(A, n))
                    ),
                        file="dbscan_index.txt", append = TRUE, sep="\n")
                    cat(toString(
                        c(n, K, p, BarOmega, fcm_simul(A, K))
                    ),
                        file="fcm_index.txt", append = TRUE, sep="\n")
                    print(c(n, K, p, BarOmega))
                },
                error = function(e){
                    str(e)
                }
                )
            }
        }
    }
}
duration <-  start_time - Sys.time()
print(duration)
