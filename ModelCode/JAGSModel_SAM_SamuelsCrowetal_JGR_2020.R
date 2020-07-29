#-------------------------------------------------------------------#
# this branch aims at studying the original code (w/o modification) 
#-------------------------------------------------------------------#

model{
  ##Placeholder for initial part of data set. 
  for(i in 1:(Nstart-1)){
    
    mu[i] <- 0
    
    
    Y.rep[i] <- 0
    
    for(j in 1:7){
      X[j,i] <- 0
    }
  }
  
  # Priors for regression parameters:
  
  # Overall intercept:
  beta0 ~ dnorm(0,0.00001)
  wP[1] <- 0    # precipitation in first time step i assumed to have no effect – moisture effects incorporated via soil water of current week
  
  
  # Main effects:
  for(j in 1:Nparms){
    beta1[j] ~ dnorm(0,0.00001)
  }
  
  #non-linear (square) effects
  for(j in 1:2){
    beta1a[j] ~ dnorm(0,0.00001)
  }
  
  # Two-way interaction effects:
  for(j in 1:10){
    
    beta2[j,1] ~ dnorm(0,0.00001)
  }
  
  
  #Priors for importance weights for each covariate (imposes Dirichlet(1) priors):
  for(j in 1:Nlag){
  dV[j]    ~ dgamma(1,1)
  dT[j]    ~ dgamma(1,1)
  dSs[j]   ~ dgamma(1,1)
  dSd[j]   ~ dgamma(1,1)
  dVrng[j] ~ dgamma(1,1)
  
  wV[j]    <- dV[j]/sum(dV[])
  wT[j]    <- dT[j]/sum(dT[])
  wSs[j]   <- dSs[j]/sum(dSs[])
  wSd[j]   <- dSd[j]/sum(dSd[])
  wVrng[j] <- dVrng[j]/sum(dVrng[])
  }

#Priors for importance weights for precipitation:

for(j in 1:(NlagP-1)){
  dP[j] ~ dgamma(1,1)
  
  wP[j+1] <- dP[j]/sum(dP[])      
  
}

#rearrange precip weights into weights at the weekly and monthly scales.

for(j in 1:4){
  
  wP.weekly[j] <- wP[j]
}

for(j in 1:6){
  wP.monthly[j] <- equals(j,1)*sum(wP[1:4]) + (1-equals(j,1))*wP[j+3]
  
}

for(i in Nstart:Nend){
  
  Y[i] ~ dnorm(mu[i],tau)     # Likelihood of observed vector of Y values (defined in script file)
  Y.rep[i] ~ dnorm(mu[i],tau) # Replicated data for evaluating model fit
  
  mu[i] <- beta0 + main.effects[i] + squared.terms[i] + interactions[i] ## regression model
  
  # Define parts involving main effects, quadratic effects, and 2-way interactions
  
  main.effects[i]  <- sum(X.effect[,i])
  interactions[i]  <- sum(sum.XX.int[,i])
  squared.terms[i] <- sum(X2.effect[,i])
  
  ## Calculate net sensitivities (derivative)
  
  dYdVPD[i] <- beta1[1] + 2*beta1a[1]*VPDant[i] + beta2[1,1]*TAant[i] + beta2[2,1]*PPTant[i] + beta2[3,1]*Sshall_ant[i] + beta2[4,1]*Sdeep_ant[i]
  dYdT[i]   <- beta1[2] + 2*beta1a[2]*TAant[i] + beta2[1,1]*VPDant[i] + beta2[5,1]*PPTant[i] + beta2[6,1]*Sshall_ant[i] + beta2[7,1]*Sdeep_ant[i]
  dYdP[i]   <- beta1[3] + beta2[2,1]*VPDant[i] + beta2[5,1]*TAant[i] + beta2[8,1]*Sshall_ant[i] + beta2[9,1]*Sdeep_ant[i]
  dYdSs[i]  <- beta1[4] + beta2[3,1]*VPDant[i] + beta2[6,1]*TAant[i] + beta2[8,1]*PPTant[i] + beta2[10,1]*Sdeep_ant[i]
  dYdSd[i]  <- beta1[5] + beta2[4,1]*VPDant[i] + beta2[7,1]*TAant[i] + beta2[9,1]*PPTant[i] + beta2[10,1]*Sshall_ant[i]
  
  dYdX[i,1] <- dYdVPD[i]
  dYdX[i,2] <- dYdT[i]
  dYdX[i,3] <- dYdP[i]
  dYdX[i,4] <- dYdSs[i]
  dYdX[i,5] <- dYdSd[i]
  
  # creating antecedent covariates; This matrix of values will end up being used to calculate the parts involving main effects, interactions, and squared terms in the regression model:
  X[1,i] <- VPDant[i]       ## Also included as squared term
  X[2,i] <- TAant[i]        ## Also included as squared term
  X[3,i] <- PPTant[i]
  X[4,i] <- Sshall_ant[i]
  X[5,i] <- Sdeep_ant[i]
  X[6,i] <- PAR[Yday[i]]    ## not included in interactions
  X[7,i] <- Vrngant[i]      ## not included in interactions
  
  # Computed antecedent values. PAR is assumed to instantaneously affect ecosystem fluxes, so no antecedent term calculated
  
  VPDant[i]     <- sum(VPDtemp[i,])
  TAant[i]      <- sum(TAtemp[i,])
  PPTant[i]     <- sum(PPTtemp[i,])
  Sshall_ant[i] <- sum(Sshalltemp[i,])
  Sdeep_ant[i]  <- sum(Sdeeptemp[i,])
  Vrngant[i]    <- sum(Vrngtemp[i,])
  
  ## Intermediate weighted values of covariates with influence over flux over the past few days (or months for ppt)
  
  for(j in 1:Nlag){
    VPDtemp[i,j]    <- wV[j]*VPD_F[Yday[i]-j+1]
    TAtemp[i,j]     <- wT[j]*TA_F[Yday[i]-j+1]
    Sshalltemp[i,j] <- wSs[j]*Sshall[Yday[i]-j+1]
    Sdeeptemp[i,j]  <- wSd[j]*Sdeep[Yday[i]-j+1]
    Vrngtemp[i,j]    <- wVrng[j]*V_rng[Yday[i]-j+1]
  }
  
  for(j in 1:NlagP){
    PPTtemp[i,j] <- wP[j]*P_Ftemp[i,j]
    P_Ftemp[i,j] <- sum(P_F[(Yday[i]-P1[j]):(Yday[i]-P2[j])])
  }
  
  for(j in 1:Nparms){
    X.effect[j,i]<- beta1[j]*X[j,i]
  }
  
  # Individual 2-way interaction terms:
  for(j in 1:10){
    
    XX.int[j,i] <- beta2[j,1]*X[ID1[j],i]*X[ID2[j],i]
    sum.XX.int[j,i] <- sum(XX.int[j,i])
  }
  # Squared terms:
  for(j in 1:2){
    X2.effect[j,i] <- beta1a[j]*pow(X[j,i],2)
  }
}  

#Priors for standard deviations in likelihood models 

tau     <- pow(sig,-2)
sig     ~ dunif(0,1000)

}
