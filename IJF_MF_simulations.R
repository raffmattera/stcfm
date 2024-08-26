#######################################################
#
#      Spatio-temporal clustered factor model
# 
#         Mattera, R. and Franses, P.H.
#
#                 Simulations
#
########################################################

rm(list=ls())

library(rio)
library(lmtest)
library(sandwich)
library(forecast)
library(ClustGeo)
library(TSclust)
library(POET)
library(simstudy)
library(splm)
library(ARDL)
library(sf)
library(maptools)
library(ggplot2)

######### A. Main functions

NGfactors <- function (Y, rmax=10) {
  p = nrow(Y)
  n = ncol(Y)
  Y <- Y - t(t(apply(t(Y), 2, mean))) %*% matrix(1, 1, n)
  c = seq(0.05, 5, length = 100)
  re = 20
  IC = array(0, c(2, re, rmax, 100))
  gT1HL = array(1, c(re))
  gT2HL = array(1, c(re))
  pi = array(1, c(re))
  ni = array(1, c(re))
  for (i in 1:re) {
    pi[i] = min(i * floor(p/re) + min(p, 5), p)
    ni[i] = min(i * floor(n/re) + min(n, 5), n)
    if (i == re) {
      pi[i] = p
      ni[i] = n
    }
    Yi = Y[1:pi[i], 1:ni[i]]
    frob = array(0, c(rmax))
    penal = array(0, c(rmax))
    for (k in 1:min(pi[i], ni[i], rmax)) {
      V <- eigen(t(Yi) %*% Yi)$vectors
      V = as.matrix(V)
      Dd <- eigen(t(Yi) %*% Yi)$values
      Dd = as.vector(Dd)
      F <- V[, 1:k]
      LamPCA = Yi %*% F/ni[i]
      uhat = Yi - LamPCA %*% t(F)
      frob[k] = sum(diag(uhat %*% t(uhat)))/(pi[i] * ni[i])
      gT1HL[i] = log((pi[i] * ni[i])/(pi[i] + ni[i])) * 
        (pi[i] + ni[i])/(pi[i] * ni[i])
      gT2HL[i] = log(min(pi[i], ni[i])) * (pi[i] + ni[i])/(pi[i] * 
                                                             ni[i])
      for (l in 1:100) {
        IC[1, i, k, l] = log(frob[k]) + c[l] * k * gT1HL[i]
        IC[2, i, k, l] = log(frob[k]) + c[l] * k * gT2HL[i]
      }
    }
  }
  rhat = array(0, c(2, re, 100))
  for (i in 1:re) {
    for (l in 1:100) {
      m = min(pi[i], ni[i], rmax)
      temp1 = which.min(IC[1, i, 1:m, l])
      rhat[1, i, l] = temp1
      temp2 = which.min(IC[2, i, 1:m, l])
      rhat[2, i, l] = temp2
    }
  }
  Sc1 = array(0, c(100))
  Sc2 = array(0, c(100))
  for (l in 1:100) {
    Sc1[l] = sd(rhat[1, , l])
    Sc2[l] = sd(rhat[2, , l])
  }
  c1vec = which(Sc1 == 0)
  ctemp1 = c1vec[1]
  c1 = c[ctemp1]
  K1HL = rhat[1, 1, ctemp1]
  c2vec = which(Sc2 == 0)
  ctemp2 = c2vec[1]
  c2 = c[ctemp2]
  K2HL = rhat[2, 1, ctemp2]
  c = 1
  IC = array(0, c(2, rmax))
  frob = array(0, c(rmax))
  penal = array(0, c(rmax))
  for (k in 1:rmax) {
    V <- eigen(t(Y) %*% Y)$vectors
    V = as.matrix(V)
    Dd <- eigen(t(Y) %*% Y)$values
    Dd = as.vector(Dd)
    F <- V[, 1:k]
    LamPCA = Y %*% F/n
    uhat = Y - LamPCA %*% t(F)
    frob[k] = sum(diag(uhat %*% t(uhat)))/(p * n)
    gT1BN = log((p * n)/(p + n)) * (p + n)/(p * n)
    gT2BN = log(min(p, n)) * (p + n)/(p * n)
    IC[1, k] = log(frob[k]) + k * gT1BN
    IC[2, k] = log(frob[k]) + k * gT2BN
  }
  K1BN = which.min(IC[1, ])
  K2BN = which.min(IC[2, ])
  result <- list(K1HL = K1HL, K2HL = K2HL, K1BN = K1BN, K2BN = K2BN, 
                 IC = IC)
  return(result)
}

STfacm <- function(data, coords, typeD="COR", Wmat=NULL, maxitr=100, G=NULL, Gmax=NULL, tol=0.0001){
  
  N <- ncol(data)
  
  # Step 1: Estimate K global factors
  
  K <- NGfactors(data)$K1BN
  
  # Step 2: Retrive global factors (GF)
  
  P <- sqrt(N)*eigen(cov(data))$vector[,1:K]
  GF <- (data%*%P)/N
  
  # Step 3: Residuals (ResTS)
  
  E <- data-GF%*%t(P)
  
  # Step 4: Spatio-temporal clustering of the residuals
  
  # 4.1: TS Distance:
  
  #DGmat <- DTWfeatures(ResTS)
  
  #D0<-matrix(NA, N, N) # Time evolution
  #for (i in 1:ncol(D0)) {
  #  for (j in 1:ncol(D0)) {
  #    D0[i,j] <- DTWmyf(ResTS[,i],ResTS[,j], DGmat[[i+1]],DGmat[[j+1]])$dval
  #  }
  #}
  D0 <- diss(t(E), typeD)
  D0<-D0/max(D0)
  
  # 4.2 SP Distance:
  
  if (is.null(Wmat)) {
    D1 <- as.matrix(dist(coords))
    D1 <- D1/max(D1)
  } else {
    Adj <- Wmat
    diag(Adj) <- 1
    D1 <- 1-Adj
  }
  
  # 4.3: Spatio-temporal:
  
  res.hc<-GeoClustf(D0,D1)
  plot(sort(res.hc$clust$height,decreasing = T), type="b", main="Screeplot - spatio-temporal clustering", ylab="Height",xlab="Clusters")
  
  if (is.null(Gmax)) {Gmax <- N-1}
  if (is.null(G)) {nclust <- elbow_finder(sort(res.hc$clust$height,decreasing = T),seq(1:Gmax))[2]} else {nclust <- G}
  plot(res.hc$clust)
  rect.hclust(res.hc$clust,k=nclust)
  clustering <- cutree(res.hc$clust, k=nclust)
  
  # Step 4: Estimating c-th local factors (c-LF)
  
  #dataclust <- list()
  cLF <- list()
  Pc <- list()
  resC <- list()
  
  for (c in 1:nclust) {
    if (table(clustering)[c] > K) {
      dataclust <- E[,clustering==c]
      rmax <- ifelse(min(table(clustering))==1, min(table(clustering)[-which(table(clustering)==1)]) , min(table(clustering)) )
      Kc <- NGfactors(dataclust, rmax=rmax)$K1BN
      Pc[[c]] <- sqrt(table(clustering)[c])*eigen(cov(dataclust))$vector[,1:Kc]
      cLF[[c]] <-  (dataclust%*%Pc[[c]])/table(clustering)[c] # if (covTestR::Nagao1973(dataclust)$p.value < 0.05)  { } else {}
    } else {
      dataclust <- E[,clustering==c]
      Pc[[c]] <- 0
      cLF[[c]] <- 0 # "No Local Factor"
    }
    
    #clf <- cLF[[c]]
    
    if (length(cLF[[c]])==1) {
      
      resC[[c]] <- E[,clustering==c]
      
    } else {
      
      resC[[c]] <- dataclust-cLF[[c]]%*%t(Pc[[c]])
      
    }
    
  }
  
  valvec <- numeric(maxitr)
  
  val <- (N*nrow(data))^(-1) * sum(unlist((resC))^2)
  valvec[1] <- val
  
  for (z in 1:maxitr) {
    
    cval <- val
    
    R <- matrix(NA,nrow(data),ncol(data))
    
    for (c in 1:nclust) {
      
      R[,which(clustering==c)] <- data[,which(clustering==c)] - cLF[[c]]%*%t(Pc[[c]])
      
    }
    
    K <- NGfactors(R)$K1BN
    
    P <- sqrt(N)*eigen(cov(R))$vector[,1:K]
    GF <- (R%*%P)/N
    
    # Step 3: Residuals (E-tilde)
    
    for (c in 1:nclust) {
      
      if (table(clustering)[c] > K) {
        E[,which(clustering==c)] <- R[,which(clustering==c)] - GF%*%t(P)[,which(clustering==c)]  + cLF[[c]]%*%t(Pc[[c]])
      } else {
        E[,which(clustering==c)] <- R[,which(clustering==c)] - GF%*%t(P)[,which(clustering==c)]
      }
    }
    
    #DGmat <- DTWfeatures(ResTS)
    
    D0<-matrix(NA, N, N) # Time evolution
    #for (i in 1:ncol(D0)) {
    #  for (j in 1:ncol(D0)) {
    #    D0[i,j] <- DTWmyf(ResTS[,i],ResTS[,j], DGmat[[i+1]],DGmat[[j+1]])$dval
    #  }
    #}
    D0 <- diss(t(E), typeD)
    D0<-D0/max(D0)
    res.hc<-GeoClustf(D0,D1)
    plot(sort(res.hc$clust$height,decreasing = T), type="b", main="Screeplot - spatio-temporal clustering", ylab="Height",xlab="Clusters")
    if (is.null(G)) {nclust <- elbow_finder(sort(res.hc$clust$height,decreasing = T),seq(1:Gmax))[2]} else {nclust <- G}
    plot(res.hc$clust)
    rect.hclust(res.hc$clust,k=nclust)
    clustering <- cutree(res.hc$clust, k=nclust)
    
    cLF <- list()
    resC <- list()
    Pc <- list()
    
    for (c in 1:nclust) {
      if (table(clustering)[c] > K) {
        dataclust <- E[,clustering==c]
        rmax <- ifelse(min(table(clustering))==1, min(table(clustering)[-which(table(clustering)==1)]) , min(table(clustering)) )
        Kc <- NGfactors(dataclust, rmax=rmax)$K1BN
        Pc[[c]] <- sqrt(table(clustering)[c])*eigen(cov(dataclust))$vector[,1:Kc]
        cLF[[c]] <-  (dataclust%*%Pc[[c]])/table(clustering)[c] # if (covTestR::Nagao1973(dataclust)$p.value < 0.05)  { } else { cLF[[c]] <- 0 }
      } else {
        dataclust <- E[,clustering==c]
        Pc[[c]] <- 0
        cLF[[c]] <- 0 # "No Local Factor"
      }
      
      #clf <- cLF[[c]]
      
      if (length(cLF[[c]])==1) {
        
        resC[[c]] <- E[,clustering==c]
        
      } else {
        
        resC[[c]] <- dataclust-cLF[[c]]%*%t(Pc[[c]])
        
      }
      
    }
    
    val <- (N*nrow(data))^(-1) * sum(unlist((resC))^2)
    valvec[z+1] <- val
    
    dd <- (cval-val)/val
    
    if( abs(dd) < tol ){ break }
    cat("*")
    
  }
  
  
  # Return final results
  
  STFM_results <- list("Global Factors"=GF,"Local Factors"=cLF,"Residuals"=E,"nclust"=nclust,"Global Loadings"=P,
                       "Local Loadings"=Pc,"Weights"=res.hc$Weights,"clustering"=clustering,"itr"=z,"obj"=val)
  
  return(STFM_results)
  
}

AndoYclust <- function (Y, NGfactor, NLfactors, Maxit = 100, tol = 0.0001) {
  
  Ngroups <- length(NLfactors)
  AY <- Y
  P <- nrow(AY)
  N <- ncol(AY)
  Z <- AY
  VEC <- eigen(cov(Z))$vectors
  Ltemp <- sqrt(N) * (VEC)[, 1:(sum(NLfactors))]
  Ftemp <- (Z%*%Ltemp)/N
  Km <- kmeans(Ltemp, Ngroups)
  LAB <- Km$cluster
  PP <- hist(LAB, br = 0:Ngroups, plot = FALSE)$counts
  FS <- matrix(0, nrow = P, ncol = sum(NLfactors))
  LS <- matrix(0, nrow = N, ncol = max(NLfactors))
  PredL <- matrix(0, nrow = P, ncol = N)
  for (i in 1:Ngroups) {
    index <- subset(1:N, LAB == i)
    Z <- Y[, index]
    VEC <- eigen(cov(Z))$vectors
    Ltemp <- sqrt(length(index)) * (VEC)[, 1:NLfactors[i]]
    Ftemp <- (Z%*%Ltemp)/length(index)
    LS[index, 1:NLfactors[i]] <- Ltemp
    if (i == 1) {
      FS[, 1:NLfactors[1]] <- Ftemp
    }
    if (i != 1) {
      FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Ftemp
    }
    PredL[, index] <- Ftemp %*% t(Ltemp)
  }
  Z <- AY - PredL
  VEC <- eigen(cov(Z))$vectors
  Ltemp <- (sqrt(N) * (VEC))[, 1:NGfactor]
  Ftemp <- (Z%*%Ltemp)/N
  FG <- Ftemp
  LG <- Ltemp
  PredG <- FG %*% t(LG)
  B <- (N*P)^(-1) *sum((AY - PredL - PredG)^2)
  for (ITE in 1:Maxit) {
    B.old <- B
    Y <- AY - PredL - PredG
    
    for (j in 1:N) {
      Er <- rep(10^7, len = Ngroups)
      for (i in 1:Ngroups) {
        if (NLfactors[i] != 0) {
          if (i == 1) {
            Ftemp <- FS[, 1:NLfactors[1]]
          }
          if (i != 1) {
            Ftemp <- FS[, (sum(NLfactors[1:(i - 1)]) + 
                             1):(sum(NLfactors[1:i]))]
          }
          Ltemp <- solve(t(Ftemp) %*% Ftemp) %*% t(Ftemp) %*% 
            Y[, j]
          Er[i] <- sum((Y[, j] - Ftemp %*% Ltemp)^2)
        }
        if (NLfactors[i] == 0) {
          Er[i] <- sum((Y[, j])^2)
        }
        LAB[j] <- subset(1:length(Er), Er == min(Er))
      }
    }
    
    for (i in 1:Ngroups) {
      if (NLfactors[i] != 0) {
        index <- subset(1:N, LAB == i)
        Z <- Y[, index]
        VEC <- eigen(cov(Z))$vectors
        Ltemp <- sqrt(length(index)) * (VEC)[, 1:NLfactors[i]]
        Ftemp <- (Z%*%Ltemp)/length(index)
        LS[index, 1:NLfactors[i]] <- Ltemp
        if (i == 1) {
          FS[, 1:NLfactors[1]] <- Ftemp
        }
        if (i != 1) {
          FS[, (sum(NLfactors[1:(i - 1)]) + 1):(sum(NLfactors[1:i]))] <- Ftemp
        }
        PredL[, index] <- Ftemp %*% t(Ltemp)
      }
    }
    
    Z <- AY - PredL
    VEC <- eigen(cov(Z))$vectors
    Ltemp <- (sqrt(N) * (VEC))[, 1:NGfactor]
    Ftemp <- (Z%*%Ltemp)/N
    FG <- Ftemp
    LG <- Ltemp
    PredG <- FG %*% t(LG)
    B <- (N*P)^(-1) *sum((AY - PredL - PredG)^2)
    if (abs((B.old - B)/B) < tol) {
      break
    }
  }
  Er <- AY - PredL - PredG
  
  results0 <- list("Factors"=FG,"Loadings"=LG,"nclust"=Ngroups,"GroupFactors" = FS,"GroupLoad"=LS,
                   "clustering"=LAB,"itr"=ITE,"obj"=sum(Er^2))
  
  return(results0)
}

GeoClustf<-function(D0,D1){
  nn <- ncol(D1)-1
  res.hc0<-hclustgeo(as.dist(D0),as.dist(D1), alpha=0, scale = FALSE)
  kk <- elbow_finder(seq(1:nn),sort(res.hc0$height,decreasing = T))[1]
  #kk <- elbow_finder2(seq(1:14),sort(res.hc0$height,decreasing = T))
  
  cr <- choicealpha(as.dist(D0),as.dist(D1),seq(0,1,0.1),kk,graph=TRUE, scale = FALSE)
  a <- cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  res.hc0<-hclustgeo(as.dist(D0),as.dist(D1), alpha=a, scale = FALSE)
  res.hc0$weight<-a
  res.hc0$dist.mat<-(1-a)*D0+a*D1
  
  return(list("clust"=res.hc0,"Weights"=a))
}

elbow_finder <- function(x_values, y_values) {
  # Max values to create line
  max_x_x <- max(x_values)
  max_x_y <- y_values[which.max(x_values)]
  max_y_y <- max(y_values)
  max_y_x <- x_values[which.max(y_values)]
  max_df <- data.frame(x = c(max_y_x, max_x_x), y = c(max_y_y, max_x_y))
  
  # Creating straight line between the max values
  fit <- lm(max_df$y ~ max_df$x)
  
  # Distance from point to line
  distances <- c()
  for(i in 1:length(x_values)) {
    distances <- c(distances, abs(coef(fit)[2]*x_values[i] - y_values[i] + coef(fit)[1]) / sqrt(coef(fit)[2]^2 + 1^2))
  }
  
  # Max distance point
  x_max_dist <- x_values[which.max(distances)]
  y_max_dist <- y_values[which.max(distances)]
  
  return(c(x_max_dist, y_max_dist))
}

fe_est <- function(model,varsid){
  
  df <- model$Set_Vars
  beta <- coef(model)[1,]
  id <- unique(df$category)
  fixef <- NULL
  
  for(i in 1:length(id)){
    df0 <- df[df$category==id[i],]
    fixef[i] <- mean(df0$price,na.rm=T)-beta%*%colMeans(as.matrix(df0[,varsid]),na.rm=T)
  }
  
  return(fixef)
  
}

######### B. Simulation study

# We simulate given the geography of the empirical application. 
# We deal with the US territory, so N=49.

# Scenario with G=3

ushp <- readShapePoly("US_State_Boundaries", IDvar="STATE_ABBR")
liststates <- ushp$STATE_ABBR[-c(1,12,40,48)]
liststates <- sort(liststates)
coords <- coordinates(ushp)
coords <- coords[-c(1,12,40,48),]
Dij <- as.matrix(dist(coords))
Dij <- Dij/max(Dij)
df2 <- usmap::statepop
df2 <- df2[df2$abbr %in% liststates,]
df2 <- df2[order(df2$abbr),]
df2$pop_2015 <- rep("Z",length(liststates))
colnames(df2)[4]<-"Cluster"

df2$Cluster[which(coords[, 1] >= -98)] <- 1
df2$Cluster[which(coords[, 1] < -98)] <- 2
df2$Cluster[which(coords[, 1] > -100 & coords[, 2] <= 39)] <- 3

usmap::plot_usmap(include =liststates, data = df2, values  = "Cluster") + 
  labs(title = "Simulated scenario with K=3") + 
  theme(legend.position = "right")

# DGPI

G3_dgp1 <- list()
temps <- c(100, 200, 400)

for (l in 1:length(temps)) {
  
  lengT <- temps[l]
  
  ARImat <- matrix(NA, 1000, 2)
  Errormat <- matrix(NA, 1000, 2)
  Galpha <- matrix(NA, 1000, 2)
  
  for (m in 1:1000) try({
    
    Y1 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[1])
    Y2 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[2])
    Y3 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[3])
    colnames(Y1) <- df2$abbr[df2$Cluster==1]
    colnames(Y2) <- df2$abbr[df2$Cluster==2]
    colnames(Y3) <- df2$abbr[df2$Cluster==3]
    
    # To all the series we extract K=1 global factor
    
    set.seed(m)
    GF <- runif(lengT, 0, 1)
    set.seed(m)
    LG <- runif(ncol(Y1)+ncol(Y2)+ncol(Y3), -2, 2)
    
    Y1 <- GF%*%t(LG[1:ncol(Y1)])
    Y2 <- GF%*%t(LG[(ncol(Y1)+1):(ncol(Y1)+ncol(Y2))])
    Y3 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+1):length(LG)])
    
    set.seed(m+1)
    CF1 <- rnorm(lengT, 0, 1)
    set.seed(100+m+1)
    L1 <- rnorm(ncol(Y1), 0, 1)
    Y1 <- Y1+CF1%*%t(L1)
    
    set.seed(m+2)
    CF2 <- rnorm(lengT, 0, 1)
    L2 <- rnorm(ncol(Y2), 0, 1)
    set.seed(100+m+2)
    Y2 <- Y2+CF2%*%t(L2)
    
    set.seed(m+3)
    CF3 <- rnorm(lengT, 0, 1)
    L3 <- rnorm(ncol(Y3), 0, 1)
    set.seed(100+m+2)
    Y3 <- Y3+CF3%*%t(L3)
    
    set.seed(m+1)
    Y1 <- Y1+genCorGen(n = lengT, nvars = ncol(Y1), corMatrix = genCorMat(nvars = ncol(Y1), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y1)), params2=rep(1,ncol(Y1)),
                       dist = "normal", wide = TRUE)[,-1]
    
    set.seed(m+2)
    Y2 <- Y2+genCorGen(n = lengT, nvars = ncol(Y2), corMatrix = genCorMat(nvars = ncol(Y2), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y2)), params2=rep(1,ncol(Y2)),
                       dist = "normal", wide = TRUE)[,-1]
    
    set.seed(m+3)
    Y3 <- Y3+genCorGen(n = lengT, nvars = ncol(Y3), corMatrix = genCorMat(nvars = ncol(Y3), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y3)), params2=rep(1,ncol(Y3)),
                       dist = "normal", wide = TRUE)[,-1]
    
    Y <- as.matrix(cbind(Y1,Y2,Y3))
    colnames(Y) <- c(df2$abbr[df2$Cluster==1],df2$abbr[df2$Cluster==2],
                     df2$abbr[df2$Cluster==3])
    
    # Extract partition and estimated factors with Ando and Bai approach
    
    resAB <- AndoYclust((Y), 1, c(1,1,1), tol=0.1)
    names(resAB$clustering) <- colnames(Y)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resAB$Factors
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:3) try({
      Fg[,resAB$clustering==g] <- resAB$GroupFactors[,g]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YAB1 <- mu1 + beta1%*%t(var_t)
    
    # Extract partition and estimated factors with our approach
    
    resST <- STfacm((Y), coords, G=3, tol=0.1)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resST$`Global Factors`
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:3) try({
      Fg[,resAB$clustering==g] <- resST$`Local Factors`[[g]]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YST1 <- mu1 + beta1%*%t(var_t)
    
    # Would we choose correct G and alpha?
    
    a <- STfacm(Y, coords, tol=0.1, Gmax=10)
    Galpha[m,1] <- a$nclust
    Galpha[m,2] <- a$Weights 
    
    # Evaluate clustering quality with RI
    
    ARImat[m,1] <- fossil::rand.index(as.numeric(df2$Cluster), resAB$clustering)
    ARImat[m,2] <- fossil::rand.index(as.numeric(df2$Cluster), resST$clustering)
    
    # Evaluate forecasting accuracy
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- as.matrix(GF)
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Fg[,df2$Cluster==1] <- CF1
    Fg[,df2$Cluster==2] <- CF2
    Fg[,df2$Cluster==3] <- CF3
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    Yfor1 <- mu1 + beta1%*%t(var_t)
    
    Errormat[m,1] <- sqrt(mean((YAB1-Yfor1)^2))
    Errormat[m,2] <- sqrt(mean((YST1-Yfor1)^2))
    
    print(paste0("T=",temps[l]," m=",m)) 
    
  })
  
  G3_dgp1[[l]] <- list(ARImat,Errormat,Galpha)
  
}

rio::export(G3_dgp1,"G3_dgp1.RDS")

# DGPII

G3_dgp2 <- list()
temps <- c(100, 200, 400)

for (l in 1:length(temps)) {
  
  lengT <- temps[l]
  
  ARImat <- matrix(NA, 1000, 2)
  Errormat <- matrix(NA, 1000, 2)
  Galpha <- matrix(NA, 1000, 2)
  
  for (m in 1:1000) try({
    
    Y1 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[1])
    Y2 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[2])
    Y3 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[3])
    colnames(Y1) <- df2$abbr[df2$Cluster==1]
    colnames(Y2) <- df2$abbr[df2$Cluster==2]
    colnames(Y3) <- df2$abbr[df2$Cluster==3]
    
    # To all the series we extract K=1 global factor
    
    set.seed(m)
    GF <- runif(lengT, 0, 1)
    set.seed(m)
    LG <- runif(ncol(Y1)+ncol(Y2)+ncol(Y3), -2, 2)
    
    Y1 <- GF%*%t(LG[1:ncol(Y1)])
    Y2 <- GF%*%t(LG[(ncol(Y1)+1):(ncol(Y1)+ncol(Y2))])
    Y3 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+1):length(LG)])
    
    set.seed(m+1)
    CF1 <- rnorm(lengT, 0, 1)
    set.seed(100+m+1)
    L1 <- rnorm(ncol(Y1), 0, 1)
    Y1 <- Y1+CF1%*%t(L1)
    
    set.seed(m+2)
    CF2 <- rnorm(lengT, 0, 1)
    L2 <- rnorm(ncol(Y2), 0, 1)
    set.seed(100+m+2)
    Y2 <- Y2+CF2%*%t(L2)
    
    set.seed(m+3)
    CF3 <- rnorm(lengT, 0, 1)
    L3 <- rnorm(ncol(Y3), 0, 1)
    set.seed(100+m+2)
    Y3 <- Y3+CF3%*%t(L3)
    
    Y <- as.matrix(cbind(Y1,Y2,Y3))
    
    C <- matrix(0.05, nrow=ncol(Y), ncol=ncol(Y))
    diag(C) <- 1
    C <- C^(Dij)
    set.seed(m+4)
    Y <- Y+genCorGen(n = lengT, nvars = ncol(Y), corMatrix = C, params1 = rep(0,ncol(Y)), params2=rep(1,ncol(Y)),
                     dist = "normal", wide = TRUE)[,-1]
    
    Y <- as.matrix(Y)
    colnames(Y) <- c(df2$abbr[df2$Cluster==1],df2$abbr[df2$Cluster==2],
                     df2$abbr[df2$Cluster==3])
    
    
    # Extract partition and estimated factors with Ando and Bai approach
    
    resAB <- AndoYclust((Y), 1, c(1,1,1), tol=0.1)
    names(resAB$clustering) <- colnames(Y)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resAB$Factors
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:3) try({
      Fg[,resAB$clustering==g] <- resAB$GroupFactors[,g]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YAB1 <- mu1 + beta1%*%t(var_t)
    
    # Extract partition and estimated factors with our approach
    
    resST <- STfacm((Y), coords, G=3, tol=0.1)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resST$`Global Factors`
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:3) try({
      Fg[,resAB$clustering==g] <- resST$`Local Factors`[[g]]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YST1 <- mu1 + beta1%*%t(var_t)
    
    # Would we choose correct G and alpha?
    
    a <- STfacm(Y, coords, tol=0.1, Gmax=10)
    Galpha[m,1] <- a$nclust
    Galpha[m,2] <- a$Weights 
    
    # Evaluate clustering quality with RI
    
    ARImat[m,1] <- fossil::rand.index(as.numeric(df2$Cluster), resAB$clustering)
    ARImat[m,2] <- fossil::rand.index(as.numeric(df2$Cluster), resST$clustering)
    
    # Evaluate forecasting accuracy
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- as.matrix(GF)
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Fg[,df2$Cluster==1] <- CF1
    Fg[,df2$Cluster==2] <- CF2
    Fg[,df2$Cluster==3] <- CF3
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    Yfor1 <- mu1 + beta1%*%t(var_t)
    
    Errormat[m,1] <- sqrt(mean((YAB1-Yfor1)^2))
    Errormat[m,2] <- sqrt(mean((YST1-Yfor1)^2))
    
    print(paste0("T=",temps[l]," m=",m))
  })
  
  G3_dgp2[[l]] <- list(ARImat,Errormat,Galpha)
  
}

rio::export(G3_dgp2,"G3_dgp2.RDS")

# Scenario with G=4

ushp <- readShapePoly("US_State_Boundaries", IDvar="STATE_ABBR")
liststates <- ushp$STATE_ABBR[-c(1,12,40,48)]
liststates <- sort(liststates)
coords <- coordinates(ushp)
coords <- coords[-c(1,12,40,48),]
Dij <- as.matrix(dist(coords))
Dij <- Dij/max(Dij)
df2 <- usmap::statepop
df2 <- df2[df2$abbr %in% liststates,]
df2 <- df2[order(df2$abbr),]
df2$pop_2015 <- rep("Z",length(liststates))
colnames(df2)[4]<-"Cluster"

df2$Cluster[which(coords[, 1] >= -102)] <- 1
df2$Cluster[which(coords[, 1] < -102)] <- 2
df2$Cluster[which(coords[, 1] > -92 & coords[, 2] >= 39)] <- 3
df2$Cluster[which(coords[, 1] > -89 & coords[, 2] <= 39)] <- 4

usmap::plot_usmap(include =liststates, data = df2, values  = "Cluster") + 
  labs(title = "Simulated scenario with K=4") + 
  theme(legend.position = "right")

# DGPI

G4_dgp1 <- list()
temps <- c(100, 200, 400)

for (l in 1:length(temps)) {
  
  lengT <- temps[l]
  
  ARImat <- matrix(NA, 1000, 2)
  Errormat <- matrix(NA, 1000, 2)
  Galpha <- matrix(NA, 1000, 2)
  
  for (m in 1:1000) try({
    
    Y1 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[1])
    Y2 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[2])
    Y3 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[3])
    Y4 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[4])
    
    # To all the series we extract K=1 global factor
    
    set.seed(m)
    GF <- runif(lengT, 0, 1)
    set.seed(m)
    LG <- runif(ncol(Y1)+ncol(Y2)+ncol(Y3)+ncol(Y4), -2, 2)
    
    Y1 <- GF%*%t(LG[1:ncol(Y1)])
    Y2 <- GF%*%t(LG[(ncol(Y1)+1):(ncol(Y1)+ncol(Y2))])
    Y3 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+1):(ncol(Y1)+ncol(Y2)+ncol(Y3))])
    Y4 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+ncol(Y3)+1):length(LG)])
    
    set.seed(m+1)
    CF1 <- rnorm(lengT, 0, 1)
    set.seed(100+m+1)
    L1 <- rnorm(ncol(Y1), 0, 1)
    Y1 <- Y1+CF1%*%t(L1)
    
    set.seed(m+2)
    CF2 <- rnorm(lengT, 0, 1)
    set.seed(100+m+2)
    L2 <- rnorm(ncol(Y2), 0, 1)
    Y2 <- Y2+CF2%*%t(L2)
    
    set.seed(m+3)
    CF3 <- rnorm(lengT, 0, 1)
    set.seed(100+m+3)
    L3 <- rnorm(ncol(Y3), 0, 1)
    Y3 <- Y3+CF3%*%t(L3)
    
    set.seed(m+4)
    CF4 <- rnorm(lengT, 0, 1)
    set.seed(100+m+4)
    L4 <- rnorm(ncol(Y4), 0, 1)
    Y4 <- Y4+CF4%*%t(L4)
    
    set.seed(m+1)
    Y1 <- Y1+genCorGen(n = lengT, nvars = ncol(Y1), corMatrix = genCorMat(nvars = ncol(Y1), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y1)), params2=rep(1,ncol(Y1)),
                       dist = "normal", wide = TRUE)[,-1]
    
    set.seed(m+2)
    Y2 <- Y2+genCorGen(n = lengT, nvars = ncol(Y2), corMatrix = genCorMat(nvars = ncol(Y2), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y2)), params2=rep(1,ncol(Y2)),
                       dist = "normal", wide = TRUE)[,-1]
    
    set.seed(m+3)
    Y3 <- Y3+genCorGen(n = lengT, nvars = ncol(Y3), corMatrix = genCorMat(nvars = ncol(Y3), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y3)), params2=rep(1,ncol(Y3)),
                       dist = "normal", wide = TRUE)[,-1]
    
    set.seed(m+4)
    Y4 <- Y4+genCorGen(n = lengT, nvars = ncol(Y4), corMatrix = genCorMat(nvars = ncol(Y4), rho = 0.3, corstr = "cs"), params1 = rep(0,ncol(Y4)), params2=rep(1,ncol(Y4)),
                       dist = "normal", wide = TRUE)[,-1]
    
    Y <- as.matrix(cbind(Y1,Y2,Y3,Y4))
    colnames(Y) <- c(df2$abbr[df2$Cluster==1],df2$abbr[df2$Cluster==2],
                     df2$abbr[df2$Cluster==3],df2$abbr[df2$Cluster==4])
    
    # Extract partition and estimated factors with Ando and Bai approach
    
    resAB <- AndoYclust((Y), 1, c(1,1,1,1), tol=0.1)
    names(resAB$clustering) <- colnames(Y)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resAB$Factors
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:4) try({
      Fg[,resAB$clustering==g] <- resAB$GroupFactors[,g]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YAB1 <- mu1 + beta1%*%t(var_t)
    
    # Extract partition and estimated factors with our approach
    
    resST <- STfacm((Y), coords, G=4, tol=0.1)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resST$`Global Factors`
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:4) try({
      Fg[,resAB$clustering==g] <- resST$`Local Factors`[[g]]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YST1 <- mu1 + beta1%*%t(var_t)
    
    # Would we choose correct G and alpha?
    
    a <- STfacm(Y, coords, tol=0.1, Gmax=30)
    Galpha[m,1] <- a$nclust
    Galpha[m,2] <- a$Weights 
    
    # Evaluate clustering quality with RI
    
    ARImat[m,1] <- fossil::rand.index(as.numeric(df2$Cluster), resAB$clustering)
    ARImat[m,2] <- fossil::rand.index(as.numeric(df2$Cluster), resST$clustering)
    
    # Evaluate forecasting accuracy
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- as.matrix(GF)
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Fg[,df2$Cluster==1] <- CF1
    Fg[,df2$Cluster==2] <- CF2
    Fg[,df2$Cluster==3] <- CF3
    Fg[,df2$Cluster==4] <- CF4
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    Yfor1 <- mu1 + beta1%*%t(var_t)
    colnames(Yfor1) <- colnames(Y)
    
    Errormat[m,1] <- sqrt(mean((YAB1-Yfor1)^2))
    Errormat[m,2] <- sqrt(mean((YST1-Yfor1)^2))
    
    print(paste0("T=",temps[l]," m=",m))
  })
  
  G4_dgp1[[l]] <- list(ARImat,Errormat,Galpha)
  
}

rio::export(G4_dgp1,"G4_dgp1.RDS")

# DGPII

G4_dgp2 <- list()
temps <- c(100, 200, 400)

for (l in 1:length(temps)) {
  
  lengT <- temps[l]
  
  ARImat <- matrix(NA, 1000, 2)
  Errormat <- matrix(NA, 1000, 2)
  Galpha <- matrix(NA, 1000, 2)
  
  for (m in 1:1000) try({
    
    Y1 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[1])
    Y2 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[2])
    Y3 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[3])
    Y4 <- matrix(NA, nrow=lengT, ncol=table(df2$Cluster)[4])
    colnames(Y1) <- df2$abbr[df2$Cluster==1]
    colnames(Y2) <- df2$abbr[df2$Cluster==2]
    colnames(Y3) <- df2$abbr[df2$Cluster==3]
    colnames(Y4) <- df2$abbr[df2$Cluster==4]
    
    # To all the series we extract K=1 global factor
    
    set.seed(m)
    GF <- runif(lengT, 0, 1)
    set.seed(m)
    LG <- runif(ncol(Y1)+ncol(Y2)+ncol(Y3)+ncol(Y4), -2, 2)
    
    Y1 <- GF%*%t(LG[1:ncol(Y1)])
    Y2 <- GF%*%t(LG[(ncol(Y1)+1):(ncol(Y1)+ncol(Y2))])
    Y3 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+1):(ncol(Y1)+ncol(Y2))+ncol(Y3)])
    Y4 <- GF%*%t(LG[(ncol(Y1)+ncol(Y2)+ncol(Y3)+1):length(LG)])
    
    set.seed(m+1)
    CF1 <- rnorm(lengT, 0, 1)
    set.seed(100+m+1)
    L1 <- rnorm(ncol(Y1), 0, 1)
    Y1 <- Y1+CF1%*%t(L1)
    
    set.seed(m+2)
    CF2 <- rnorm(lengT, 0, 1)
    L2 <- rnorm(ncol(Y2), 0, 1)
    set.seed(100+m+2)
    Y2 <- Y2+CF2%*%t(L2)
    
    set.seed(m+3)
    CF3 <- rnorm(lengT, 0, 1)
    L3 <- rnorm(ncol(Y3), 0, 1)
    set.seed(100+m+2)
    Y3 <- Y3+CF3%*%t(L3)
    
    set.seed(m+4)
    CF4 <- rnorm(lengT, 0, 1)
    L4 <- rnorm(ncol(Y4), 0, 1)
    Y4 <- Y4+CF4%*%t(L4)
    
    Y <- cbind(Y1,Y2,Y3,Y4)
    
    C <- matrix(0.1, nrow=ncol(Y), ncol=ncol(Y))
    diag(C) <- 1
    C <- C^(Dij)
    set.seed(m+4)
    Y <- Y+genCorGen(n = lengT, nvars = ncol(Y), corMatrix = C, params1 = rep(0,ncol(Y)), params2=rep(1,ncol(Y)),
                     dist = "normal", wide = TRUE)[,-1]
    
    Y <- as.matrix(Y)
    colnames(Y) <- c(df2$abbr[df2$Cluster==1],df2$abbr[df2$Cluster==2],
                     df2$abbr[df2$Cluster==3],df2$abbr[df2$Cluster==4])
    
    # Extract partition and estimated factors with Ando and Bai approach
    
    resAB <- AndoYclust((Y), 1, c(1,1,1,1), tol=0.1)
    names(resAB$clustering) <- colnames(Y)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resAB$Factors
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:4) try({
      Fg[,resAB$clustering==g] <- resAB$GroupFactors[,g]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YAB1 <- mu1 + beta1%*%t(var_t)
    
    # Extract partition and estimated factors with our approach
    
    resST <- STfacm((Y), coords, G=4, tol=0.1)
    
    # Panel VAR forecast:
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- resST$`Global Factors`
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    for (g in 1:4) try({
      Fg[,resAB$clustering==g] <- resST$`Local Factors`[[g]]
    },silent=T)
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    YST1 <- mu1 + beta1%*%t(var_t)
    
    # Would we choose correct G and alpha?
    
    a <- STfacm(Y, coords, tol=0.1, Gmax=30)
    Galpha[m,1] <- a$nclust
    Galpha[m,2] <- a$Weights 
    
    # Evaluate clustering quality with RI
    
    ARImat[m,1] <- fossil::rand.index(as.numeric(df2$Cluster), resAB$clustering)
    ARImat[m,2] <- fossil::rand.index(as.numeric(df2$Cluster), resST$clustering)
    
    # Evaluate forecasting accuracy
    
    id <- Y
    for (j in 1:ncol(id)) {
      id[,j] <- rep(colnames(Y)[j],nrow(id))
    }
    timeid <- Y
    for (j in 1:ncol(id)) {
      timeid[,j] <- 1:nrow(id)
    }
    
    F0 <- as.matrix(GF)
    Fg <- matrix(NA, nrow=nrow(Y), ncol=ncol(Y))
    Fg[,df2$Cluster==1] <- CF1
    Fg[,df2$Cluster==2] <- CF2
    Fg[,df2$Cluster==3] <- CF3
    Fg[,df2$Cluster==4] <- CF4
    Fg[is.na(Fg)] <- 0
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(Y),
                      "F0"=rep(F0,ncol(Y)),
                      "Fg"=as.vector(Fg))
    
    model1 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                  lags=1,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    mu1 <- fe_est(model1,c(6:8))
    beta1 <- coef(model1)[1,]
    var_t <- cbind(Y[nrow(Y),],rep(F0[nrow(F0),],ncol(Y)),Fg[nrow(Fg),])
    Yfor1 <- mu1 + beta1%*%t(var_t)
    colnames(Yfor1) <- colnames(Y)
    
    Errormat[m,1] <- sqrt(mean((YAB1-Yfor1)^2))
    Errormat[m,2] <- sqrt(mean((YST1-Yfor1)^2))
    
    print(paste0("T=",temps[l]," m=",m))
  })
  
  G4_dgp2[[l]] <- list(ARImat,Errormat,Galpha)
  
}

rio::export(G4_dgp2,"G4_dgp2.RDS")

################ C. Results Tables RI and Loss

# Tables, mean and std: DGPI

MeanRI <- matrix(NA, nrow=3, ncol=4)
MeanLoss <- matrix(NA, nrow=3, ncol=4)
SdRI <- matrix(NA, nrow=3, ncol=4)
SdLoss <- matrix(NA, nrow=3, ncol=4)

colnames(MeanRI) <- c("Ando and Bai", "Our method", "Ando and Bai", "Our method")
rownames(MeanRI) <- c("T=100", "T=200", "T=400")
colnames(MeanLoss) <- colnames(MeanRI)
rownames(MeanLoss) <- rownames(MeanRI)
colnames(SdRI) <- colnames(MeanRI)
rownames(SdRI) <- rownames(MeanRI)
colnames(SdLoss) <- colnames(MeanRI)
rownames(SdLoss) <- rownames(MeanRI)

for (i in 1:3) {
  MeanRI[i,1:2] <- round(apply(na.omit(G3_dgp1[[i]][[1]]), 2, mean),3)
  MeanRI[i,3:4] <- round(apply(na.omit(G4_dgp1[[i]][[1]]), 2, mean),3)
  MeanLoss[i,1:2] <- round(apply(na.omit(G3_dgp1[[i]][[2]]), 2, mean)*100,3)
  MeanLoss[i,3:4] <- round(apply(na.omit(G4_dgp1[[i]][[2]]), 2, mean)*100,3)
  
  SdRI[i,1:2] <- round(apply(na.omit(G3_dgp1[[i]][[1]]), 2, sd),3)
  SdRI[i,3:4] <- round(apply(na.omit(G4_dgp1[[i]][[1]]), 2, sd),3)
  SdLoss[i,1:2] <- round(apply(na.omit(G3_dgp1[[i]][[2]]), 2, sd)*100,3)
  SdLoss[i,3:4] <- round(apply(na.omit(G4_dgp1[[i]][[2]]), 2, sd)*100,3)
  
}

MeanRI
MeanLoss
SdRI
SdLoss

# Tables, mean and std: DGPII

# Subscript previous DGPI table

for (i in 1:3) {
  MeanRI[i,1:2] <- round(apply(na.omit(G3_dgp2[[i]][[1]]), 2, mean),3)
  MeanRI[i,3:4] <- round(apply(na.omit(G4_dgp2[[i]][[1]]), 2, mean),3)
  MeanLoss[i,1:2] <- round(apply(na.omit(G3_dgp2[[i]][[2]]), 2, mean)*100,3)
  MeanLoss[i,3:4] <- round(apply(na.omit(G4_dgp2[[i]][[2]]), 2, mean)*100,3)
  
  SdRI[i,1:2] <- round(apply(na.omit(G3_dgp2[[i]][[1]]), 2, sd),3)
  SdRI[i,3:4] <- round(apply(na.omit(G4_dgp2[[i]][[1]]), 2, sd),3)
  SdLoss[i,1:2] <- round(apply(na.omit(G3_dgp2[[i]][[2]]), 2, sd)*100,3)
  SdLoss[i,3:4] <- round(apply(na.omit(G4_dgp2[[i]][[2]]), 2, sd)*100,3)
  
}

MeanRI
MeanLoss
SdRI
SdLoss

## Plots

tlength <- 1 # select 1 for T=100, 2 is T=200, 3 is T=400

K3.1 <- G3_dgp1[[tlength]][[3]][,1]
K3.2 <- G3_dgp2[[tlength]][[3]][,1]
K4.1 <- G4_dgp1[[tlength]][[3]][,1]
K4.2 <- G4_dgp2[[tlength]][[3]][,1]

A3.1 <- G3_dgp1[[tlength]][[3]][,2]
A3.2 <- G3_dgp2[[tlength]][[3]][,2]
A4.1 <- G4_dgp1[[tlength]][[3]][,2]
A4.2 <- G4_dgp2[[tlength]][[3]][,2]

boxplot(cbind(K3.1,K3.2,
              K4.1,K4.2),main="Number of clusters", ylab="Selected clusters",xlab="Simulation schemes")

boxplot(cbind(A3.1,A3.2,
              A4.1,A4.2),main="Selected mixing parameter", ylab="Mixing parmeter",xlab="Simulation schemes")
