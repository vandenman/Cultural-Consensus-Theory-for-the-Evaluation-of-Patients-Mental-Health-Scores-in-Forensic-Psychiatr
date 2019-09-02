library(rjags)

model <- function() {
  
  #Data
  for (i in 1:n) {
    for (k in 1:m) {
      tau[i,k] <- E[i] / lam[k]
      pX[i,k,1] <- pnorm((a[i]*g[1]) + b[i],Truth[k],tau[i,k])
      for (c in 2:(C-1)) {
        pX[i,k,c] <- pnorm((a[i]*g[c]) + b[i],Truth[k],tau[i,k]) - sum(pX[i,k,1:(c-1)])
      }
      pX[i,k,C] <- (1 - sum(pX[i,k,1:(C-1)]))
      X[i,k] ~ dcat(pX[i,k,1:C])
    }
  }
  
  #Parameters
  for (k in 1:m) {
    Truth[k] ~ dnorm(Tmu,Ttau)
    ilogitT[k] <- ilogit(Truth[k])
    lam[k] ~ dgamma(lamtau,lamtau)
  }
  for (c in 1:(C-2)) {
    dg[c] ~ dnorm(0,.1)
  }
  dg2[1:(C-2)] <- dg[1:(C-2)]
  dg2[C-1] <- -sum(dg[1:(C-2)])
  g <- sort(dg2)
  
  for(c in 1:(C-1)) {
    ilogitg[c] <- ilogit(g[c])
  }
  for (i in 1:n) {
    E[i] ~ dgamma(pow(Emu,2)*Etau,Emu*Etau)
    a[i] ~ dgamma(atau,atau)
    b[i] ~ dnorm(bmu,btau)
  }
  #Hyperparameters
  Tmu ~ dnorm(0,.1)
  Ttau ~ dgamma(1,.1)
  bmu <- 0
  btau ~ dgamma(1,.1)
  
  amu <- 1
  atau ~ dgamma(1,.1)
  Emu ~ dgamma(4,4)
  Etau ~ dgamma(4,4)
  lammu <- 1
  lamtau ~ dgamma(4,4)
  
}

nr = 50
ni = 20
nc = 5

test <- genData(nr, ni, nc)
datList <- list(
	# data
	X = test$x,
	n = nr,
	m = ni,
	C = nc
)

rr <- jags(data = datList, parameters.to.save = c("E", "a", "b", "T"), model.file = model)
