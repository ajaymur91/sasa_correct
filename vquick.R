#!/usr/bin/env Rscript
# Load libraries
options(warn=-1)
suppressMessages(library(igraph))
suppressMessages(library(dplyr))

################ Read Command line inputs ###########
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)<3) {
  stop("At least 3 arguments must be supplied (input file).n", call.=FALSE)
} 
############### Inputs #########################
Ns <- as.integer(args[1])   # Generate Ns configurations
Np <- as.integer(args[2])   # Np particles in Box
Si <- as.integer(args[3])   # Required number of successful insertions 

Ld <- 0.35 # Connected if distance less than Ld
Lp <- 0.0210870
Ln <- 0.0169807

# Box dimensions
################################################
Lx <- read.table('./L',header=FALSE)$V1
Ly <- read.table('./L',header=FALSE)$V2
Lz <- read.table('./L',header=FALSE)$V3
################################################

# Read positions of parent
GX <- read.table('./X',header=FALSE)
GY <- read.table('./Y',header=FALSE)
GZ <- read.table('./Z',header=FALSE)

################################################
# Adjacency correction factor (Acf)
# AdjM_true <- Argon-like-AdjM ** Acf
# i.e. Cation-cation connections have to be neglected
# i.e. anion-anion connections have to be neglected
dim <- Np
Acf1 <- matrix(rep(1,dim*dim), nrow = dim, ncol = dim)
for(i in seq(1,dim,1)){
  for(j in seq(1,dim,1)){
    if((i+j) %% 2==0) { Acf1[i,j]=0 }
  }
}
################################################

# Function to find if cluster is connected
iscon1 <- function(X,Y,Z,Ld) {
  Dim <- length(X)
  # Calc adj matrix (Notice Acf factor)
  AdjM <- 1*(as.matrix(dist(cbind(X,Y,Z),method = "euclidian"))<Ld) * Acf1
  # Convert to graph
  G <- graph_from_adjacency_matrix(AdjM,mode = "undirected",diag = FALSE,add.colnames = NA)
  # Check connected
  return(1*(is_connected(G)))
}

mratio1 <- function(X,Y,Z) {
  Dim <- length(X)
  # Calc adj matrix (Notice Acf factor)
  AdjM <- 1*(as.matrix(dist(cbind(X,Y,Z),method = "euclidian"))<Ld) * Acf
  # Convert to graph
  G <- graph_from_adjacency_matrix(AdjM,mode = "undirected",diag = FALSE,add.colnames = NULL)
  # Calc mratio
  CV <- seq(Dim,0.5,-2)
  L <- length(CV)
  M <- rep(0,L)
  n=1
  for(l in CV){
  M[n] <- is_connected(G-vertex(l))
  n=n+1
  }
  return(L/sum(M))
}


# Save Connectivity for each configuration in a vector: C

################################################
################################################
k=0
SUM=0
while(SUM < Si)
{
  C <- rep(0,Ns)
  X <- matrix(round(runif(Ns*1,min=0,max = Lx),3), ncol = 1)
  Y <- matrix(round(runif(Ns*1,min=0,max = Ly),3), ncol = 1)
  Z <- matrix(round(runif(Ns*1,min=0,max = Lz),3), ncol = 1)
  ################################################
  
  # Combine fixed particle with configurations to create cluster of size (N+1)
  X2 <- cbind(X,rev(GX))
  Y2 <- cbind(Y,rev(GY))
  Z2 <- cbind(Z,rev(GZ))
  
  X <- X2
  Y <- Y2
  Z <- Z2
  rm(X2,Y2,Z2)
  
  #I th configuration
  for(i in seq(1,Ns,1)){
    C[i] <-  iscon1(c(X[i,]),c(Y[i,]),c(Z[i,]),Ld)
  }
  IND <- which(C==1)
  if(exists("SX"))
  {
    SIND <- rbind(SIND,C)
    SX <- rbind(SX,X[IND,])
    SY <- rbind(SY,Y[IND,])
    SZ <- rbind(SZ,Z[IND,])
  } else {
    SIND <- C
    SX <- X[IND,]
    SY <- Y[IND,]
    SZ <- Z[IND,]
  }
  k=k+1
  SUM <- SUM + sum(C)
}

## Mratio
#M1 <- rep(0,Si)
#for(j in seq(1,Si,1)){
# M1[j] <-  mratio1(c(SX[j,]),c(SY[j,]),c(SZ[j,]))
#}
write.table(SX[1:Si,],file = paste0('XC'),row.names=F,col.names=F)
write.table(SY[1:Si,],file = paste0('YC'),row.names=F,col.names=F)
write.table(SZ[1:Si,],file = paste0('ZC'),row.names=F,col.names=F)
#write.table(M1,file = paste0('M1'),row.names=F,col.names=F)
################################################
################################################
# Adjacency correction factor (Acf)
# AdjM_true <- Argon-like-AdjM ** Acf
# i.e. Cation-cation connections have to be neglected
# i.e. anion-anion connections have to be neglected
dim <- Np + 1
Acf <- matrix(rep(1,dim*dim), nrow = dim, ncol = dim)
for(i in seq(1,dim,1)){
  for(j in seq(1,dim,1)){
    if((i+j) %% 2==0) { Acf[i,j]=0 }
  }
}
################################################

# Function to find if cluster is connected
iscon2 <- function(X,Y,Z,Ld) {
  Dim <- length(X)
  # Calc adj matrix (Notice Acf factor)
  AdjM <- 1*(as.matrix(dist(cbind(X,Y,Z),method = "euclidian"))<Ld) * Acf
  # Convert to graph
  G <- graph_from_adjacency_matrix(AdjM,mode = "undirected",diag = FALSE)
  # Check connected
  return(1*(is_connected(G)))
}


mratio2 <- function(X,Y,Z) {
  Dim <- length(X)
  # Calc adj matrix (Notice Acf factor)
  AdjM <- 1*(as.matrix(dist(cbind(X,Y,Z),method = "euclidian"))<Ld) * Acf
  # Convert to graph
  G <- graph_from_adjacency_matrix(AdjM,mode = "undirected",diag = FALSE,add.colnames = NULL)
  # Calc mratio
  CV1 <- seq(Dim,0.5,-2)
  CV2 <- seq(Dim-1,0.5,-2)
  L1 <- length(CV1)
  L2 <- length(CV2)
  M <- rep(0,L1*L2)
  n=1
  for(l in CV1){
    for(o in CV2){
      M[n] <- is_connected(G-vertex(l,o))
      n=n+1
    }
  }
  return(L1*L2/sum(M))
}
################################################
################################################
Ns <- 15
Ti <- 1
for(j in seq(1,Si,1)){
k=0
SUM=0
while(SUM < Ti)
{
  C <- rep(0,Ns)
  X <- matrix(round(runif(Ns*1,min=0,max = Lx),3), ncol = 1)
  Y <- matrix(round(runif(Ns*1,min=0,max = Ly),3), ncol = 1)
  Z <- matrix(round(runif(Ns*1,min=0,max = Lz),3), ncol = 1)
  ################################################
  
  # Combine fixed particle with configurations to create cluster of size (N+1)
  X2 <- cbind(X,rev(SX[j,]))
  Y2 <- cbind(Y,rev(SY[j,]))
  Z2 <- cbind(Z,rev(SZ[j,]))
  
  X <- X2
  Y <- Y2
  Z <- Z2
  rm(X2,Y2,Z2)
  
  #I th configuration
  for(i in seq(1,Ns,1)){
    C[i] <-  iscon2(c(X[i,]),c(Y[i,]),c(Z[i,]),Ld)
  }
  IND <- which(C==1)
  SUMC <- sum(C)
  if(SUMC>0)
  {
  if(exists("TX"))
  {
    TIND <- rbind(TIND,C)
    TX <- rbind(TX,X[IND,][1,])
    TY <- rbind(TY,Y[IND,][1,])
    TZ <- rbind(TZ,Z[IND,][1,])
  } else {
    TIND <- C
    TX <- X[IND,][1,]
    TY <- Y[IND,][1,]
    TZ <- Z[IND,][1,]
  }
  }
  k=k+1
  SUM <- SUM + sum(C)
}
}
write.table(TX,file = paste0('TX'),row.names=F,col.names=F)
write.table(TY,file = paste0('TY'),row.names=F,col.names=F)
write.table(TZ,file = paste0('TZ'),row.names=F,col.names=F)

  Ltx <- length(TX[,1])
  Lty <- length(TX[1,])
  Cl <- rep(0,Ltx)
for(i in seq(1,Ltx,1)){
   Cl[i] <-  iscon1(c(TX[i,1:(Lty-1)]),c(TY[i,1:(Lty-1)]),c(TZ[i,1:(Lty-1)]),Ld)
 }
#M12 <- rep(0,Si)
#for(j in seq(1,Si,1)){
#  M12[j] <- mratio1(c(TX[j,]),c(TY[j,]),c(TZ[j,]))
#}
#
#M2 <- rep(0,Si)
#for(j in seq(1,Si,1)){
#  M2[j] <- mratio2(c(TX[j,]),c(TY[j,]),c(TZ[j,]))
#}
write.table(Cl,file = paste0('Cl.txt'),row.names=F,col.names=F)

#write.table(M12,file = paste0('M12'),row.names=F,col.names=F)
#write.table(M2,file = paste0('M2'),row.names=F,col.names=F)
######################################################################

# Bootstrapping for error analysis
#Nb <- 100
#Veff <- rep(0,Nb)
# Also Resample are with replacement
#frac <- 1
#for(i in seq(1,Nb,1)){
#  # Find length of sample for 10% of the data
#  S <- sample_frac(data.frame(SIND),frac,replace = TRUE)
#  Veff[i] <- (Lx*Ly*Lz)^1 * (sum(S) / (Ns*k))  
#}
#StatsV <- cbind(mean(Veff),sd(Veff),min(Veff),max(Veff))
#cat(mean(Veff))
#write.table(mean(Veff),file = paste0('volume.txt'),row.names = F,col.names = F)
