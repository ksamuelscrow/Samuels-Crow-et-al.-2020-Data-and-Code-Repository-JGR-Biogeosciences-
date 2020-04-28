### This script file will allow the user to run the SAM model to evaluate controls on NEE or ET for 
### either the US-Vcp or US-Mpj site. The user should update the file to reflect the working directory,
### input data files, and flux of interest. Items that the user should change are indicated by “XXXX”.

#Part 1: Load data, construct data list for JAGS model, generate initial values

setwd(“XXXX”) # set to your working directory
library(rjags) # install rjags if necessary

dataIN = read.csv(“XXXX.csv”) #Choose the covariate file of interest. This should be a file that has all covariates for all time periods to ensure that 6-months of precipitation is supplied to the model – not just growing season precipitation
YIN    = read.csv(“XXXX.csv”) #Choose the response variable file of interest. This file includes growing-season response (Y) variables of interest along with indices that link the Y variables back to covariates of interest.

# The covariate data file (for dataIN) should include (1) the serial date, and (2) the following covariates (as columns in the data file):
# (a) precipitation
# (b) air temperature
# (c) VPD (average)
# (d) VPD (min)
# (e) VPD (max)
# (f) PAR
# (g) shallow soil moisture
# (h) deep soil moisture

# The response variable data file (for YIN) should include (as columns):
# (a) ET
# (b) NEE
# (c) index to link the growing season Y variables back to the appropriate row in the covariate data set

# Nparms is the number of driving variables included to calculate main effects. 

# Below, define the range of VPD in any given day (β¬6¬ in Tables S1a, S1b, S4) and assign column for Y values

V1 = dataIN[,”XX1”] - dataIN[,”XX2”] ## Calculate the difference between minimum and maximum daily VPD (XX1 = column number for VPD_max, XX2 = column number for VPDmin)
Y  = YIN[,”XX3”] # Change "XX3" to the column for the response variable of interest (ET or NEE)

Nj <- 10

jIND <- read.csv("InteractionIND.csv") ## This file (see E.1, below) provides indices to calculate interactions between covariates

Nstart = “XXXX”  ## Choose the starting index. This is an index for a row in the Y data file. The value in the indexed row should be greater than 1 to accommodate calculation of antecedent values.
Nend   = “XXXX”  ## Choose the ending index. 
Yday = YIN[,“XXXX”] ## Choose column in YIN that provides indices linking response variables with covariates

# Data read into model -- covariates are centered and standardized in the data list (IMPORTANT: replace name of covariate in quotes with the appropriate column number in the dataIN file):

data = list(Nstart = Nstart, Nend = Nend, Nlag = 7, NlagP = 9, Nparms = 7, Yday = Yday, ID1 = jIND[,2], ID2 = jIND[,3], Y = Y, VPD_F = (dataIN[,”VPDavg”]-mean(dataIN[,”VPDavg”], na.rm = TRUE))/sd(dataIN[,”VPDavg”], na.rm = TRUE),
            TA_F = (dataIN[,”Tair”]-mean(dataIN[,”Tair”], na.rm = TRUE))/sd(dataIN[,”Tair”], na.rm = TRUE), P_F = (dataIN[,”precip”]-mean(dataIN[,”precip”], na.rm = TRUE))/sd(dataIN[,”precip”], na.rm = TRUE), PAR = (dataIN[,”PAR”]-mean(dataIN[,”PAR”], na.rm = TRUE))/sd(dataIN[,”PAR”], na.rm = TRUE), Sshall = (dataIN[,”Sshall”]-mean(dataIN[,”Sshall”], na.rm = TRUE))/sd(dataIN[,”Sshall”], na.rm = TRUE), Sdeep  = (dataIN[,”Sdeep”]-mean(dataIN[,”Sdeep”], na.rm = TRUE))/sd(dataIN[,”Sdeep”], na.rm = TRUE), 
            V_rng = (V1-mean(V1))/sd(V1),
            P1 = c(6, 13, 20, 27, 55, 83, 111, 139, 167),
            P2 = c(0, 7, 14, 21, 28, 56, 84, 112, 140))

#Initial values are estimated using a linear model. As in the data list (above), covariates are centered and standardized. Replace name of covariate in quotes with the appropriate column number in the dataIN file.

X1  = cbind((dataIN[Yday[Nstart:Nend],”VPDavg”]-mean(dataIN[Yday[Nstart:Nend], ”VPDavg”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”VPDavg”], na.rm = TRUE), (dataIN[Yday[Nstart:Nend], ”Tair”]-mean(dataIN[Yday[Nstart:Nend], ”Tair”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”Tair”], na.rm = TRUE), (dataIN[Yday[Nstart:Nend], ”precip”]-mean(dataIN[Yday[Nstart:Nend], ”precip”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”precip”], na.rm = TRUE), (dataIN[Yday[Nstart:Nend], ”PAR”]-mean(dataIN[Yday[Nstart:Nend], ”PAR”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”PAR”], na.rm = TRUE), (dataIN[Yday[Nstart:Nend], ”Sshall”]-mean(dataIN[Yday[Nstart:Nend], ”Sshall”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”Sshall”], na.rm = TRUE), (dataIN[Yday[Nstart:Nend], ”Sdeep”]-mean(dataIN[Yday[Nstart:Nend], ”Sdeep”], na.rm = TRUE))/sd(dataIN[Yday[Nstart:Nend], ”Sdeep”], na.rm = TRUE), (V1[Yday[Nstart:Nend]]-mean(V1[Yday[Nstart:Nend]], na.rm = TRUE))/sd(V1[Yday[Nstart:Nend]], na.rm = TRUE))
X1a = cbind(X1[,1]^2, X1[,2]^2) #non-linear (square) terms calculated for VPD and Tair
X2  = cbind(X1[,1]*X1[,2], X1[,1]*X1[,3], X1[,1]*X1[,4], X1[,1]*X1[,5], X1[,2]*X1[,3], X1[,2]*X1[,4], X1[,2]*X1[,5], X1[,3]*X1[,4], X1[,3]*X1[,5], X1[,4]*X1[,5]) #Interactions incorporated into linear model used to estimate initial values

fit <- lm(Y[Nstart:Nend] ~ X1[,1] + X1[,2] + X1[,3] + X1[,4] + X1[,5] + X1[,6] + X1[,7] + X1a[,1] + X1a[,2] + X2[,1] + X2[,2]  + X2[,3]  + X2[,4]  + X2[,5]  + X2[,6]  + X2[,7]  + X2[,8]  + X2[,9]  + X2[,10])

beta0  = fit$coefficients[1]
beta1  = fit$coefficients[2:8]
beta1a = fit$coefficients[9:10]

beta2 <- matrix(data = 0, nrow = 10, ncol = 1)

beta2[1:4,1]   = as.numeric(fit$coefficients[11:14])
beta2[5:7,1]   = as.numeric(fit$coefficients[15:17])
beta2[8:9,1]   = as.numeric(fit$coefficients[18:19])
beta2[10,1]    = as.numeric(fit$coefficients[20])

inits = list(list(beta0 = beta0, beta1 = beta1, beta1a = beta1a, beta2 = beta2, sig = 1, sigV = 1, sigT = 1, sigP = 1, sigPAR = 1, sigSs = 1, sigSd = 1, sigVrng = 1, muV=-1, muT=0.1, muP=0.3, muPAR = 0.5, muSs = -0.1, muSd = 0.1, muVrng = -1),
             list(beta0 = beta0/10, beta1 = beta1/10, beta1a = beta1a/10, beta2 = beta2/10, sig = 1/2, sigV = 1/2, sigT = 1/2, sigP = 1/2, sigPAR = 1/2, sigSs = 1/2, sigSd = 1/2, sigVrng = 1/2, muV=-1/2, muT=0.01, muP=0.03, muPAR = 0.05, muSs = -0.01, muSd = 0.01, muVrng = -0.1),
             list(beta0 = beta0*10, beta1 = beta1*10, beta1a = beta1a*10, beta2 = beta2*10, sig = 2, sigV = 2, sigT = 2, sigP = 2, sigPAR = 2, sigSs = 2, sigSd = 2, sigVrng = 2, muV=-2, muT=1, muP=3, muPAR = 5, muSs = -1, muSd = 1, muVrng = -10))

#####################################################################

#Part 2: Initialize JAGS Model

n.adapt = 100 # adjust this number (and n.iter) as appropriate 
n.iter = 1000
n.chains = 3

jm1.b=jags.model("JAGS_Model_FluxData.R",
                 data=data,
                 n.chains=n.chains,
                 n.adapt=n.adapt,
                 inits = inits)

#####################################################################

#Part 3: Run JAGS Model

load.module("dic")

### Choose the parameters to monitor. For this analysis, we allowed the model to converge while monitoring variables included in "zc1"
### below before monitoring dYdX (zc1dYdX), Xant (zc1X), or Y (zc1Y)

n.iter = 40000

zc1 = coda.samples(jm1.b,variable.names=c("deviance","beta0","beta1","beta1a","beta2", "wT","wV","wP","wSs","wSd","wP.weekly","wP.monthly",'wVrng', "sig"),
                   n.iter=n.iter,thin = 40, n.adapt=1000)

### Evaluate convergence before monitoring the variables below

n.iter = 4000

zc1dYdX = coda.samples(jm1.b,variable.names = c("dYdVPD","dYdT","dYdSs","dYdSd","dYdP"),n.iter)

zc1X = coda.samples(jm1.b,variable.names = c("VPDant","TAant","PPTant","Sshall_ant","Sdeep_ant","PAR","Vrngant"),n.iter)

n.iter = 5000

zc1Y = coda.samples(jm1.b,variable.names=c("Y.rep"),n.iter=n.iter)
