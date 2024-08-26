#######################################################
#
#      Spatio-temporal clustered factor model
# 
#         Mattera, R. and Franses, P.H.
#
#              empirical results 
#                 (Section 5)
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
library(dynlm)
library(splm)
library(ARDL)
library(dplyr)
library(tidyr)
library(panelvar)
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

STfacm <- function(data, coords, typeD="COR", Wmat=NULL, maxitr=100, G=NULL, tol=0.0001){
  
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
  
  if (is.null(G)) {nclust <- elbow_finder(sort(res.hc$clust$height,decreasing = T),seq(1:(N-1)))[2]} else {nclust <- G}
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
    

    D0<-matrix(NA, N, N) # Time evolution
    D0 <- diss(t(E), typeD)
    D0<-D0/max(D0)
    res.hc<-GeoClustf(D0,D1)
    plot(sort(res.hc$clust$height,decreasing = T), type="b", main="Screeplot - spatio-temporal clustering", ylab="Height",xlab="Clusters")
    if (is.null(G)) {nclust <- elbow_finder(sort(res.hc$clust$height,decreasing = T),seq(1:(N-1)))[2]} else {nclust <- G}
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

  cr <- choicealpha(as.dist(D0),as.dist(D1),seq(0,1,0.1),kk,graph=TRUE, scale = FALSE)
  a <- cr$range.alpha[which.min(abs(cr$Qnorm[,1]-cr$Qnorm[,2]))]
  res.hc0<-hclustgeo(as.dist(D0),as.dist(D1), alpha=a, scale = FALSE)
  res.hc0$weight<-a
  res.hc0$dist.mat<-(1-a)*D0+a*D1
  
  return(list("clust"=res.hc0,"Weights"=a))
}

elbow_finder <- function(x_values, y_values) {
  # automatic approach for elbow detection
  
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

BICval <- function(model){
  
  # Han, Philips and Sul (2017)
  # Lag Length Selection in Panel Autoregression
  
  # IC1:l,i9
  
  n <- model$nof_groups
  tk <- model$nof_observations/model$nof_groups
  k <- model$lags
  p <- length(model$dependent_vars)
  cnt <- log(sqrt(n)*tk)/(sqrt(n)*tk)
  sigma <- sum(residuals(model)^2)/(p*n*tk)
  
  val <- sigma + k*cnt
  
  return(val)
}

fe_est <- function(model,varsid){
  
  # pvar approach to price equation (applicable function to fixed effects pvar)
  
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

######## A. Application to house prices

Data <- rio::import("PriceUS.csv")
colnames(Data) <- c("State", "Year", "Quarter", "Index")
Data <- Data[order(Data$State),]
Prices <- matrix(Data$Index,nrow=193,ncol=length(unique(Data$State)))
colnames(Prices) <- unique(Data$State)

# Plot Figure 4

PricePlot <- apply(Prices[,-c(1,12)], 2, function(x){(x/x[1])*100})
Dates <- seq(as.Date("1975-01-01"),as.Date("2023-01-01"), by = "3 month")
plot(Dates, PricePlot[,1], type="l", ylim=c(min(PricePlot),max(PricePlot)), ylab="Price",xlab="Year", main="House prices in the U.S. (1975=100)")
for (i in 2:ncol(PricePlot[,-1])) {
  lines(Dates,PricePlot[,i], col=colours()[50+i])
}

LatLong <- rio::import("LatLongUS.xlsx")
LatLong <- LatLong[-40,]
data(usaww)
Wmat <- usaww/usaww
Wmat[is.nan(Wmat)] <- 0

dPrice <- apply(Prices, 2, function(x){(x-Hmisc::Lag(x))/Hmisc::Lag(x)})
dPrice <- dPrice[-1,-c(1,12)]

# Full sample clustering

FullSample <- STfacm(as.matrix(dPrice), LatLong[-c(1,12),c(2,3)], typeD="COR")

# Plot Figure 6 (Figure 5 has been obtained from an external website given the cluster labels)

plot(Dates[-1],FullSample$`Global Factors`, type="l", lwd=3, xlab="Years", ylim=c(0.2,-0.2), ylab="Value",main="Estimated factors")
lines(Dates[-1],FullSample$`Local Factors`[[1]], col="red")
lines(Dates[-1],FullSample$`Local Factors`[[3]], col="blue")
lines(Dates[-1],FullSample$`Local Factors`[[2]], col="darkgreen")
lines(Dates[-1],FullSample$`Local Factors`[[4]], col="gold")
legend("topright", c("Global", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
       col=c("black","red","darkgreen", "blue", "gold"), lty=1, lwd=c(3,1,1,1,1), cex=0.5)

# Construct panel dataset

id <- dPrice
for (i in 1:ncol(id)) {
  id[,i] <- rep(colnames(dPrice)[i],nrow(id))
}
timeid <- dPrice
for (i in 1:ncol(id)) {
  timeid[,i] <- 1:nrow(id)
}
F0 <- FullSample$`Global Factors`
Fg <- dPrice
Fg[,FullSample$clustering==1] <- FullSample$`Local Factors`[[1]]
Fg[,FullSample$clustering==2] <- FullSample$`Local Factors`[[2]]
Fg[,FullSample$clustering==3] <- FullSample$`Local Factors`[[3]]
Fg[,FullSample$clustering==4] <- FullSample$`Local Factors`[[4]]

# Load observable variables z

CPI <- rio::import("CPI_growth.csv")
CPI <- as.matrix(CPI)
RealInc <- rio::import("RealInc_growth.csv")
RealInc <- as.matrix(RealInc[,order(colnames(RealInc))])
ltIR <- rio::import("LongTerm_ir.csv")
ltIR <- as.matrix(ltIR)

# Choose the PVAR order

models_est <- list()

for (p in 1:5) {
  
  df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                    "price"=as.vector(dPrice), "cpi"=as.vector(CPI),
                    "ri"=as.vector(RealInc),"itr"=as.vector(ltIR),
                    "F0"=rep(F0,ncol(dPrice)),
                    "Fg"=as.vector(Fg))
  
  model0 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr","F0","Fg"),
                                lags=p,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  models_est[[p]] <- model0
}

BICvals <- unlist(lapply(models_est, function(x){BICval(x)}))
psel <- which.min(BICvals)

# Estimate PVAR(psel) model and its restricted versions

df <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                 "price"=as.vector(dPrice), "cpi"=as.vector(CPI),
                 "ri"=as.vector(RealInc),"itr"=as.vector(ltIR),
                 "F0"=rep(F0,ncol(dPrice)),
                 "Fg"=as.vector(Fg))

model1 <- models_est[[1]]

model2 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr","F0"),
                              lags=psel,
                              transformation = "demean",
                              data=df0,
                              panel_identifier = c("id","time"))

model3 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr"),
                              lags=psel,
                              transformation = "demean",
                              data=df0,
                              panel_identifier = c("id","time"))

model4 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                              lags=psel,
                              transformation = "demean",
                              data=df0,
                              panel_identifier = c("id","time"))

model5 <- panelvar::pvarfeols(dependent_vars = c("price","F0"),
                              lags=psel,
                              transformation = "demean",
                              data=df0,
                              panel_identifier = c("id","time"))

model6 <- panelvar::pvarfeols(dependent_vars = c("price"),
                              lags=psel,
                              transformation = "demean",
                              data=df0,
                              panel_identifier = c("id","time"))

# Evaluating in-sample fitting

# Fixed-effects

mu1 <- fe_est(model1,c(9:14))
mu2 <- fe_est(model2,c(8:12))
mu3 <- fe_est(model3,c(7:10))
mu4 <- fe_est(model4,c(6:8))
mu5 <- fe_est(model5,c(5,6))
mu6 <- fe_est(model6,c(4))

# Coefficients

beta1 <- coef(model1)[1,]
beta2 <- coef(model2)[1,]
beta3 <- coef(model3)[1,]
beta4 <- coef(model4)[1,]
beta5 <- coef(model5)[1,]
beta6 <- coef(model6)[1,]

# In-sample fit

IS1 <- matrix(rep(mu1, each=nrow(dPrice))+beta1%*%t(model1$Set_Vars_with_NAs[,c(9:14)]), nrow=nrow(dPrice), ncol=ncol(dPrice))
IS2 <- matrix(rep(mu2, each=nrow(dPrice))+beta2%*%t(model2$Set_Vars_with_NAs[,c(8:12)]), nrow=nrow(dPrice), ncol=ncol(dPrice))
IS3 <- matrix(rep(mu3, each=nrow(dPrice))+beta3%*%t(model3$Set_Vars_with_NAs[,c(7:10)]), nrow=nrow(dPrice), ncol=ncol(dPrice))
IS4 <- matrix(rep(mu4, each=nrow(dPrice))+beta4%*%t(as.matrix(model4$Set_Vars_with_NAs[,c(6:8)])), nrow=nrow(dPrice), ncol=ncol(dPrice))
IS5 <- matrix(rep(mu5, each=nrow(dPrice))+beta5%*%t(as.matrix(model5$Set_Vars_with_NAs[,c(5,6)])), nrow=nrow(dPrice), ncol=ncol(dPrice))
IS6 <- matrix(rep(mu6, each=nrow(dPrice))+beta6%*%t(as.matrix(model6$Set_Vars_with_NAs[,c(4)])), nrow=nrow(dPrice), ncol=ncol(dPrice))

# In-sample errors

E1 <- matrix(dPrice-IS1, nrow=nrow(dPrice), ncol=ncol(dPrice))
E2 <- matrix(dPrice-IS2, nrow=nrow(dPrice), ncol=ncol(dPrice))
E3 <- matrix(dPrice-IS3, nrow=nrow(dPrice), ncol=ncol(dPrice))
E4 <- matrix(dPrice-IS4, nrow=nrow(dPrice), ncol=ncol(dPrice))
E5 <- matrix(dPrice-IS5, nrow=nrow(dPrice), ncol=ncol(dPrice))
E6 <- matrix(dPrice-IS6, nrow=nrow(dPrice), ncol=ncol(dPrice))

# Accuracy in-sample

RMSE1 <- apply(E1, 2, function(x){sqrt(mean(x^2, na.rm=T))})
RMSE2 <- apply(E2, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
RMSE3 <- apply(E3, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
RMSE4 <- apply(E4, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
RMSE5 <- apply(E5, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
RMSE6 <- apply(E6, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
MAE1 <- apply(E1, 2, function(x){(mean(abs(x), na.rm=T))}) 
MAE2 <- apply(E2, 2, function(x){(mean(abs(x), na.rm=T))}) 
MAE3 <- apply(E3, 2, function(x){(mean(abs(x), na.rm=T))}) 
MAE4 <- apply(E4, 2, function(x){(mean(abs(x), na.rm=T))}) 
MAE5 <- apply(E5, 2, function(x){(mean(abs(x), na.rm=T))}) 
MAE6 <- apply(E6, 2, function(x){(mean(abs(x), na.rm=T))}) 

library(tidyr)

# Figure 7a

RMSE <- data.frame(RMSE5/RMSE4,RMSE6/RMSE4)
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",2:3,"/m1")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("In-sample fit: relative accuracy (RMSE)")

# Fgiure 7b

MAE <- data.frame(MAE5/MAE4,MAE6/MAE4)
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",2:3,"/m1")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("In-sample fit: relative accuracy (MAE)")

######### Rolling window procedure

oos <- 32 # 8 years of oos

rmse_m1 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
rmse_m2 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
rmse_m3 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
rmse_m4 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
rmse_m5 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
rmse_m6 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))

mae_m1 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
mae_m2 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
mae_m3 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
mae_m4 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
mae_m5 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
mae_m6 <- matrix(NA, nrow=oos, ncol=ncol(dPrice))

m1for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
m2for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
m3for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
m4for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
m5for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))
m6for <- matrix(NA, nrow=oos, ncol=ncol(dPrice))

library(foreach)
library(doParallel)
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

clust_res_oos <- foreach(i = 1:oos, .packages = c("TSclust","ClustGeo")) %dopar% {
  resi <- STfacm(as.matrix(dPrice[i:(159+i),]), LatLong[-c(1,12),c(2,3)], typeD="COR")
  return(resi)
}
stopCluster(cl)
selll <- NULL

for (i in 1:oos){
  
  clust0 <- clust_res_oos[[i]]
  G <- length(table(clust0$clustering))
  
  # Construct panel dataset
  
  id <- dPrice[i:(159+i),]
  for (j in 1:ncol(id)) {
    id[,j] <- rep(colnames(dPrice)[j],nrow(id))
  }
  timeid <- dPrice[i:(159+i),]
  for (j in 1:ncol(id)) {
    timeid[,j] <- 1:nrow(id)
  }
  F0 <- clust0$`Global Factors`
  Fg <- matrix(NA, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  for (g in 1:G) try({
    Fg[,clust0$clustering==g] <- clust0$`Local Factors`[[g]]
  },silent=T)
  Fg[is.na(Fg)] <- 0
  
  
  # Choose the PVAR order
  
  models_est <- list()
  
  for (p in 1:5) {
    
    df0 <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                      "price"=as.vector(dPrice[i:(159+i),]), "cpi"=as.vector(CPI[i:(159+i),]),
                      "ri"=as.vector(RealInc[i:(159+i),]),"itr"=as.vector(ltIR[i:(159+i),]),
                      "F0"=rep(F0,ncol(dPrice)),
                      "Fg"=as.vector(Fg))
    
    
    model0 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr","F0","Fg"),
                                  lags=p,
                                  transformation = "demean",
                                  data=df0,
                                  panel_identifier = c("id","time"))
    models_est[[p]] <- model0
  }
  
  BICvals <- unlist(lapply(models_est, function(x){BICval(x)}))
  psel <- which.min(BICvals)
  selll[i] <- psel
  
  # Estimate PVAR(psel) model and its restricted versions
  
  df <- data.frame("id"=as.vector(id), "time"=as.vector(timeid),
                   "price"=as.vector(dPrice[i:(159+i),]), "cpi"=as.vector(CPI[i:(159+i),]),
                   "ri"=as.vector(RealInc[i:(159+i),]),"itr"=as.vector(ltIR[i:(159+i),]),
                   "F0"=rep(F0,ncol(dPrice)),
                   "Fg"=as.vector(Fg))
  
  model1 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr","F0","Fg"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  model2 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr","F0"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  model3 <- panelvar::pvarfeols(dependent_vars = c("price","cpi","ri","itr"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  model4 <- panelvar::pvarfeols(dependent_vars = c("price","F0","Fg"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  model5 <- panelvar::pvarfeols(dependent_vars = c("price","F0"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  model6 <- panelvar::pvarfeols(dependent_vars = c("price"),
                                lags=psel,
                                transformation = "demean",
                                data=df0,
                                panel_identifier = c("id","time"))
  
  
  # Evaluating in-sample fitting
  
  mu1 <- fe_est(model1,c(9:14))
  mu2 <- fe_est(model2,c(8:12))
  mu3 <- fe_est(model3,c(7:10))
  mu4 <- fe_est(model4,c(6:8))
  mu5 <- fe_est(model5,c(5,6))
  mu6 <- fe_est(model6,c(4))
  
  beta1 <- coef(model1)[1,]
  beta2 <- coef(model2)[1,]
  beta3 <- coef(model3)[1,]
  beta4 <- coef(model4)[1,]
  beta5 <- coef(model5)[1,]
  beta6 <- coef(model6)[1,]
  
  IS1 <- matrix(rep(mu1, each=nrow(dPrice[i:(159+i),]))+beta1%*%t(model1$Set_Vars_with_NAs[,c(9:14)]), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  IS2 <- matrix(rep(mu2, each=nrow(dPrice[i:(159+i),]))+beta2%*%t(model2$Set_Vars_with_NAs[,c(8:12)]), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  IS3 <- matrix(rep(mu3, each=nrow(dPrice[i:(159+i),]))+beta3%*%t(model3$Set_Vars_with_NAs[,c(7:10)]), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  IS4 <- matrix(rep(mu4, each=nrow(dPrice[i:(159+i),]))+beta4%*%t(as.matrix(model4$Set_Vars_with_NAs[,c(6:8)])), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  IS5 <- matrix(rep(mu5, each=nrow(dPrice[i:(159+i),]))+beta5%*%t(as.matrix(model5$Set_Vars_with_NAs[,c(5,6)])), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  IS6 <- matrix(rep(mu6, each=nrow(dPrice[i:(159+i),]))+beta6%*%t(as.matrix(model6$Set_Vars_with_NAs[,c(4)])), nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E1 <- matrix(dPrice[i:(159+i),]-IS1, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E2 <- matrix(dPrice[i:(159+i),]-IS2, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E3 <- matrix(dPrice[i:(159+i),]-IS3, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E4 <- matrix(dPrice[i:(159+i),]-IS4, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E5 <- matrix(dPrice[i:(159+i),]-IS5, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  E6 <- matrix(dPrice[i:(159+i),]-IS6, nrow=nrow(dPrice[i:(159+i),]), ncol=ncol(dPrice))
  
  rmse_m1[i,] <- apply(E1, 2, function(x){sqrt(mean(x^2, na.rm=T))})
  rmse_m2[i,] <- apply(E2, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
  rmse_m3[i,] <- apply(E3, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
  rmse_m4[i,] <- apply(E4, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
  rmse_m5[i,] <- apply(E5, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
  rmse_m6[i,] <- apply(E6, 2, function(x){sqrt(mean(x^2, na.rm=T))}) 
  mae_m1[i,] <- apply(E1, 2, function(x){(mean(abs(x), na.rm=T))}) 
  mae_m2[i,] <- apply(E2, 2, function(x){(mean(abs(x), na.rm=T))}) 
  mae_m3[i,] <- apply(E3, 2, function(x){(mean(abs(x), na.rm=T))}) 
  mae_m4[i,] <- apply(E4, 2, function(x){(mean(abs(x), na.rm=T))}) 
  mae_m5[i,] <- apply(E5, 2, function(x){(mean(abs(x), na.rm=T))}) 
  mae_m6[i,] <- apply(E6, 2, function(x){(mean(abs(x), na.rm=T))}) 
  
  # OOS forecasts h=1
  
  var_t <- cbind(dPrice[159+i,],CPI[159+i,],RealInc[159+i,],rep(ltIR[159+i],ncol(dPrice)),rep(F0[nrow(F0),],ncol(dPrice)),Fg[nrow(Fg),])
  m1for[i,] <- mu1 + beta1%*%t(var_t)
  m2for[i,] <- mu2 + beta2%*%t(var_t[,-6])
  m3for[i,] <- mu3 + beta3%*%t(var_t[,-c(5,6)])
  m4for[i,] <- mu4 + beta4%*%t(var_t[,-c(2,3,4)])
  m5for[i,] <- mu5 + beta5%*%t(var_t[,-c(2,3,4,6)])
  m6for[i,] <- mu6 + beta6%*%t(var_t[,1])
  
  cat(paste0((i/oos)*100),"%")
}

# Figure 7c

RMSE <- data.frame(colMeans(rmse_m5/rmse_m4),
                   colMeans(rmse_m6/rmse_m4))
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",2:3,"/m1")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (RMSE)")

# Figure 7d

MAE <- data.frame(colMeans(mae_m5/mae_m4),
                  colMeans(mae_m6/mae_m4))
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",2:3,"/m1")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (MAE)")

# Out of sample errors

osE1 <- dPrice[161:192,]-m1for
osE2 <- dPrice[161:192,]-m2for
osE3 <- dPrice[161:192,]-m3for
osE4 <- dPrice[161:192,]-m4for
osE5 <- dPrice[161:192,]-m5for
osE6 <- dPrice[161:192,]-m6for

# Figure 8a

RMSE_oos <- rbind((apply(osE4, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE5, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE6, 2, function(x){sqrt(mean(x^2))})))
RMSE <- apply(RMSE_oos, 1, function(x){x/RMSE_oos[1,]})
RMSE <- data.frame(RMSE[,-1])
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",2:3,"/m1")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (RMSE)")

# Figure 8b

MAE_oos <- rbind((apply(osE4, 2, function(x){(mean(abs(x)))})),
                 (apply(osE5, 2, function(x){(mean(abs(x)))})),
                 (apply(osE6, 2, function(x){(mean(abs(x)))})))
MAE <- apply(MAE_oos, 1, function(x){x/MAE_oos[1,]})
MAE <- data.frame(MAE[,-1])
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",2:3,"/m1")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (MAE)")

# Table 3

xtable::xtable(cbind(RMSE,MAE), digits=4)

# Differences across clusters

RMSEtab <- cbind( (apply(osE4, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE5, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE6, 2, function(x){sqrt(mean(x^2))})))
rownames(RMSEtab) <- colnames(dPrice)
colnames(RMSEtab) <- paste0("m",1:3)
xtable::xtable(RMSEtab, digits=4)
round( colMeans(RMSEtab), 4)

MAEtab <- cbind((apply(osE4, 2, function(x){(mean(abs(x)))})),
                (apply(osE5, 2, function(x){(mean(abs(x)))})),
                (apply(osE6, 2, function(x){(mean(abs(x)))})))
rownames(MAEtab) <- colnames(dPrice)
colnames(MAEtab) <- paste0("m",1:3)
xtable::xtable(MAEtab, digits=4)
round( colMeans(MAEtab), 4)

colMeans(RMSEtab[FullSample$clustering==1,])
colMeans(RMSEtab[FullSample$clustering==2,])
colMeans(RMSEtab[FullSample$clustering==3,])
colMeans(RMSEtab[FullSample$clustering==4,])

colMeans(MAEtab[FullSample$clustering==1,])
colMeans(MAEtab[FullSample$clustering==2,])
colMeans(MAEtab[FullSample$clustering==3,])
colMeans(MAEtab[FullSample$clustering==4,])

###### OOS testing

# CW test:

CWtest1 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  fhat <- osE5[,j]^2-(osE4[,j]^2-(m5for[,j]-m4for[,j])^2)
  cwstat <- sqrt(length(fhat))*(mean(fhat)/sd((fhat-mean(fhat))))
  CWtest1[j,1] <- cwstat
  CWtest1[j,2] <- 1-pt(cwstat, length(fhat)-1)
}
sum(CWtest1[,1] > 1.282)/ncol(dPrice)

CWtest2 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  fhat <- osE6[,j]^2-(osE4[,j]^2-(m6for[,j]-m4for[,j])^2)
  cwstat <- sqrt(length(fhat))*(mean(fhat)/sd((fhat-mean(fhat))))
  CWtest2[j,1] <- cwstat
  CWtest2[j,2] <- 1-pt(cwstat, length(fhat)-1)
}
sum(CWtest2[,1] > 1.282)/ncol(dPrice)

CWtests <- cbind(CWtest1[,2],CWtest2[,2])
rownames(CWtests) <- colnames(dPrice)
xtable::xtable(CWtests[,1:2])

# DM test:

DMtest1 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  ediff <- osE5[,j]^2-osE4[,j]^2
  reg0 <- lm(ediff ~ 1)
  dmstat <- coeftest(reg0, vcov=NeweyWest(reg0, prewhite=F, lag=(length(ediff)-2)))[3]
  DMtest1[j,1] <- dmstat
  DMtest1[j,2] <- 1-pt(dmstat, df=length(ediff)-1,lower.tail = T)
}
sum(DMtest1[,1] >= 1.28)/ncol(dPrice)
sum(DMtest1[,1] > 0)/ncol(dPrice)

DMtest2 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  ediff <- osE6[,j]^2-osE4[,j]^2
  reg0 <- lm(ediff ~ 1)
  dmstat <- coeftest(reg0, vcov=NeweyWest(reg0, prewhite=F, lag=(length(ediff)-2)))[3]
  DMtest2[j,1] <- dmstat
  DMtest2[j,2] <- 1-pt(dmstat, df=length(ediff)-1,lower.tail = T)
}
sum(DMtest2[,1] >= 1.28)/ncol(dPrice)

DMtests <- cbind(DMtest1[,2],DMtest2[,2])
rownames(DMtests) <- colnames(dPrice)
xtable::xtable(DMtests)

# Table 4

xtable::xtable(cbind(CWtests, DMtests))

############### Results with observable factors

# Figure 9a

RMSE <- data.frame(colMeans(rmse_m2/rmse_m1),
                   colMeans(rmse_m3/rmse_m1))
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",5:6,"/m4")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (RMSE)")

# Figure 9b

MAE <- data.frame(colMeans(mae_m2/mae_m1),
                  colMeans(mae_m3/mae_m1))
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",5:6,"/m4")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (MAE)")

# Figure 9c

RMSE_oos <- rbind((apply(osE1, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE2, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE3, 2, function(x){sqrt(mean(x^2))})))
RMSE <- apply(RMSE_oos, 1, function(x){x/RMSE_oos[1,]})
RMSE <- data.frame(RMSE[,-1])
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",5:6,"/m4")

df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (RMSE)")

# Figure 9d
MAE_oos <- rbind((apply(osE1, 2, function(x){(mean(abs(x)))})),
                 (apply(osE2, 2, function(x){(mean(abs(x)))})),
                 (apply(osE3, 2, function(x){(mean(abs(x)))})))
MAE <- apply(MAE_oos, 1, function(x){x/MAE_oos[1,]})
MAE <- data.frame(MAE[,-1])
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",5:6,"/m4")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (MAE)")

# Table 5

xtable::xtable(cbind(RMSE,MAE), digits=4)

# Testing

CWtest1 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  fhat <- osE2[,j]^2-(osE1[,j]^2-(m2for[,j]-m1for[,j])^2)
  cwstat <- sqrt(length(fhat))*(mean(fhat)/sd((fhat-mean(fhat))))
  CWtest1[j,1] <- cwstat
  CWtest1[j,2] <- 1-pt(cwstat, length(fhat)-1)
}
sum(CWtest1[,1] > 1.282)/ncol(dPrice)

CWtest2 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  fhat <- osE3[,j]^2-(osE1[,j]^2-(m3for[,j]-m1for[,j])^2)
  cwstat <- sqrt(length(fhat))*(mean(fhat)/sd((fhat-mean(fhat))))
  CWtest2[j,1] <- cwstat
  CWtest2[j,2] <- 1-pt(cwstat, length(fhat)-1)
}
sum(CWtest2[,1] > 1.282)/ncol(dPrice)

CWtests <- cbind(CWtest1[,2],CWtest2[,2])
rownames(CWtests) <- colnames(dPrice)
xtable::xtable(CWtests[,1:2])

# DM test:

DMtest1 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  ediff <- osE2[,j]^2-osE1[,j]^2
  reg0 <- lm(ediff ~ 1)
  dmstat <- coeftest(reg0, vcov=NeweyWest(reg0, prewhite=F, lag=(length(ediff)-2)))[3]
  DMtest1[j,1] <- dmstat
  DMtest1[j,2] <- 1-pt(dmstat, df=length(ediff)-1,lower.tail = T)
}
sum(DMtest1[,1] >= 1.28)/ncol(dPrice)
sum(DMtest1[,1] > 0)/ncol(dPrice)

DMtest2 <- matrix(NA, ncol(dPrice), 2)
for (j in 1:ncol(dPrice)) {
  ediff <- osE3[,j]^2-osE1[,j]^2
  reg0 <- lm(ediff ~ 1)
  dmstat <- coeftest(reg0, vcov=NeweyWest(reg0, prewhite=F, lag=(length(ediff)-2)))[3]
  DMtest2[j,1] <- dmstat
  DMtest2[j,2] <- 1-pt(dmstat, df=length(ediff)-1,lower.tail = T)
}
sum(DMtest2[,1] >= 1.28)/ncol(dPrice)

DMtests <- cbind(DMtest1[,2],DMtest2[,2])
rownames(DMtests) <- colnames(dPrice)
xtable::xtable(DMtests)

# Table 6

xtable::xtable(cbind(CWtests, DMtests))


# Are observables alone better than cluster? Figure 10

# Figure 10a

RMSE <- data.frame(colMeans(rmse_m1/rmse_m4),
                   colMeans(rmse_m2/rmse_m4),
                   colMeans(rmse_m3/rmse_m4))
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",4:6,"/m1")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (RMSE)")

# Figure 10b

MAE <- data.frame(colMeans(mae_m1/mae_m4),
                  colMeans(mae_m2/mae_m4),
                  colMeans(mae_m3/mae_m4))
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",4:6,"/m1")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling in-sample fit: relative accuracy (MAE)")

# Figure 10c

RMSE_oos <- rbind((apply(osE4, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE1, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE2, 2, function(x){sqrt(mean(x^2))})),
                  (apply(osE3, 2, function(x){sqrt(mean(x^2))})))
RMSE <- apply(RMSE_oos, 1, function(x){x/RMSE_oos[1,]})
RMSE <- data.frame(RMSE[,-1])
rownames(RMSE) <- colnames(dPrice)
colnames(RMSE) <- paste0("m",4:6,"/m1")
df_long <- RMSE %>% 
  gather(key = "Model", value = "RMSE")
ggplot(df_long, aes(x = "", y = RMSE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "RMSE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (RMSE)")

# Figure 10d

MAE_oos <- rbind((apply(osE4, 2, function(x){(mean(abs(x)))})),
                 (apply(osE1, 2, function(x){(mean(abs(x)))})),
                 (apply(osE2, 2, function(x){(mean(abs(x)))})),
                 (apply(osE3, 2, function(x){(mean(abs(x)))})))
MAE <- apply(MAE_oos, 1, function(x){x/MAE_oos[1,]})
MAE <- data.frame(MAE[,-1])
rownames(MAE) <- colnames(dPrice)
colnames(MAE) <- paste0("m",4:6,"/m1")
df_long <- MAE %>% 
  gather(key = "Model", value = "MAE")
ggplot(df_long, aes(x = "", y = MAE, fill = Model)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(x = "", y = "MAE", fill = "Model")+
  ggtitle("Rolling out-of-sample: relative accuracy (MAE)")
