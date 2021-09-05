#!/usr/bin/env Rscript
v1 <- read.table('child_0/volume_Cl.txt')$V1
v2 <- read.table('child_1/volume_Na.txt')$V1

L <- length(read.table('child_0/volume.txt')$V1)
v <- c(read.table('child_0/volume.txt')$V1,read.table('child_1/volume.txt')$V1)
Cl <- c(read.table('child_0/Cl.txt')$V1,read.table('child_1/Cl.txt')$V1)+1
V1 <- v1+v2

#V1 <- c(rep(v1,L),rep(v2,L))

#M1 <- c(read.table('child_0/M1')$V1,read.table('child_1/M1')$V1)
#M2 <- c(read.table('child_0/M2')$V1,read.table('child_1/M2')$V1)
#M12 <- c(read.table('child_0/M12')$V1,read.table('child_1/M12')$V1)
#VE <- (exp(-(-log(V1) -log(v*M12) + log(M2))))
#VE2 <- VE[is.finite(VE)]

print(mean( V1*(v/Cl)  )) 
