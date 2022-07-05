
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
##***   Example of using predictJAGS
##*  Jays data from the Swiss MHB
#   example extracted from the session 7.9.7 in the chapter 7
#   of the AHM book vol.1 (Kéry and Royle, 2016)
#   code source (https://github.com/mikemeredith/AHM_code/blob/master/AHM1_ch07/AHM1_07.09.R)
###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###
#* Ismael V. Brack jul/2022


library(unmarked)
library(jagsUI)
source("func_predictJAGS.R")



# AHM1 book code... get Jays data and fit model ---------------------------------
  
data(jay)

crPiFun <- function(p) {
  p1 <- p[,1] # Extract the columns of the p matrix, one for
  p2 <- p[,2] #   each of J = 3 sample occasions
  p3 <- p[,3]
  cbind(      # define multinomial cell probabilities:
    "100" = p1 * (1-p2) * (1-p3),
    "010" = (1-p1) * p2 * (1-p3),
    "001" = (1-p1) * (1-p2) * p3,
    "110" = p1 * p2 * (1-p3),
    "101" = p1 * (1-p2) * p3,
    "011" = (1-p1) * p2 * p3,
    "111" = p1 * p2 * p3,
    "10x" = p1*(1-p2),
    "01x" = (1-p1)*p2,
    "11x" = p1*p2)
}


o2y <- matrix(1, 3, 10)

# Grab the data objects
covinfo <- jay$covinfo
gridinfo <- jay$gridinfo
sitecovs <- jay$sitecovs
caphist <- jay$caphist

# Get observation covariates to use in model
# Day of year, sample intensity and survey duration.
# Standardize them.
day <- scale(covinfo[,c("date1", "date2", "date3")])
dur <- as.matrix(covinfo[,c("dur1", "dur2", "dur3")])
dur[is.na(dur)] <- mean(dur, na.rm=TRUE) # Pad the 6 missing values
intensity <- dur / sitecovs[,"length"] # Sample rate = duration/length
dur <- scale(dur)
intensity <- scale(intensity)


reps <- apply(!is.na(day), 1, sum)
day[reps==2,3] <- 0
dur[reps==2,3] <- 0
intensity[reps==2, 3] <- 0

# Store the observation covariates in a list
obscovs <- list(intensity = intensity, dur = dur, day = day)

# Standardize site covariates
sitecovs[,"elev"] <- scale(sitecovs[,"elev"])
sitecovs[,"forest"] <- scale(sitecovs[,"forest"])
# NOTE: length is NOT standardized, but over-written with its reciprocal
sitecovs[,"iLength"] <- 1 / sitecovs[,"length"]

# Create unmarkedFrame (need crPiFun above and unmarked loaded)
caphist <- as.matrix(caphist)
mhb.umf <- unmarkedFrameMPois(y=caphist, siteCovs=as.data.frame(sitecovs),
                              obsCovs=obscovs, obsToY=o2y, piFun="crPiFun")

# 7.9.7 Bayesian Analysis of the MHB Data
# ------------------------------------------------------------------------
# Extract data and do data augmentation up to M = 400
y <- as.matrix(getY(mhb.umf))

# Now  we have to stretch out the encounter frequencies into individuals...
# There were 439 unique individuals observed during the survey
eh <- unlist(dimnames(y)[2])
ehid <- col(y)[y>0 & !is.na(y)]  # Column ids, index to encounter history
eh <- eh[ehid]
siteid <- row(y)[y>0 & !is.na(y)] # Site ids
y <- y[y > 0 & !is.na(y)]   # Positive counts
eh <- rep(eh, y)
siteid <- rep(siteid, y)

eh.mat <- matrix(NA,nrow=length(eh),ncol=3)
for(i in 1:length(eh)){
  eh.mat[i,]<- as.numeric(unlist(strsplit(eh[i],split="")))
}


# Define some things and do the data augmentation
nsites = nrow(sitecovs)
nind <- nrow(eh.mat)
M <- 800
y <- rbind(eh.mat, matrix(0, nrow=(M-nind), ncol=3))

# Augment site ID
site <- c(siteid, rep(NA, M-nind))
sitecovs <- siteCovs(mhb.umf)

obscovs <- obsCovs(mhb.umf)
Intensity <- matrix(obscovs[,"intensity"], nrow=nsites, ncol=3,byrow=TRUE)

# Bundle data for BUGS
# go to 'RELATIONSHIPS PREDICTION' and import the JAGS output if you want to skip this part
data <- list(y = y, J = 3, M = M , nsites=nsites, X = as.matrix(sitecovs),
             Intensity= Intensity, group=site)
str(data)

# Specify model in BUGS language
cat("
model {
  # Prior distributions
  p0 ~ dunif(0,1)
  alpha0 <- log(p0 / (1-p0))
  alpha1 ~ dnorm(0, 0.01)
  beta0 ~ dnorm(0,0.01)
  beta1 ~ dnorm(0,0.01)
  beta2 ~ dnorm(0,0.01)
  beta3 ~ dnorm(0,0.01)
  psi <- sum(lambda[]) / M   # psi is a derived parameter
  # Model for abundance: lambda depends on Elev, Length, Forest
  for(s in 1:nsites){
    log(lambda[s]) <- beta0 + beta1 * X[s,1] + beta2*X[s,4] + beta3*X[s,3] # modified from the original code to use iLength
    probs[s] <- lambda[s] / sum(lambda[])
  }
  # Model for individual encounter histories
  for(i in 1:M){
    group[i] ~ dcat(probs[])  # Group == site membership
    z[i] ~ dbern(psi)         # Data augmentation variables
    # Observation model: p depends on Intensity
    for(j in 1:J){
      logit(p[i,j]) <- alpha0 + alpha1 * Intensity[group[i],j]
      pz[i,j] <- p[i,j] * z[i]
      y[i,j] ~ dbern(pz[i,j])
    }
  }
}
",fill=TRUE,file="model.txt")

# Parameters monitored
parameters <- c("p0", "alpha0", "alpha1", "beta0", "beta1", "beta2",
                "beta3", "psi")

# Initial values
zst <- c(rep(1,M-100), rep(0,100)) ## see errata
inits <- function(){
  list (p0 = runif(1), alpha1 = runif(1), beta0=runif(1),
        beta1=rnorm(1), z= zst ) }

# MCMC settings
ni <- 11000   ;   nb <- 1000   ;   nt <- 4   ;   nc <- 3

# Call JAGS from R and summarize marginal posteriors

out <- jags(data, inits, parameters, "model.txt",
            # n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni)
            n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni, parallel=TRUE)  # ~~~ for testing

print(out, digits = 2)

saveRDS(out, "outJAGS.rds")

# code from AHM1 book up to here ---******************************

# RELATIONSHIPS PREDICTION ----------------------------------------------
library(ggplot2)

out <- readRDS("outJAGS.rds")

##* Effect of Intensity (sample rate) on detection probability

# Create an object containing Intensity values to predict detection
newInten <- data.frame(Intensity=seq(min(Intensity),max(Intensity),,100))

# Predicting...
pred.Inten<- predictJAGS(fm=out,params=c("alpha0","alpha1"),
                           newdata=newInten,link="logit",quants=c(0.025,0.5,0.975),appendData=T)

names(pred.Inten)[2:4] <- c("lower","median","upper") # rename quantiles

# Figure!
ggplot(pred.Inten, aes(x=Intensity,y=mean,ymax=upper,ymin=lower)) +
  geom_line(size=1.4) + geom_ribbon(alpha=0.4, fill="gray40") +
  theme_classic()

##* Effect of elevation on local abundance
# Create an object containing elevation values to predict detection (keeping the mean value for forest and length)
newElev <- data.frame(elev=seq(min(sitecovs$elev),max(sitecovs$elev),,100),
                       length=mean(sitecovs$iLength),
                       forest=mean(sitecovs$forest))

# Predicting...
pred.elev<- predictJAGS(fm=out,params=c("beta0","beta1","beta2","beta3"),
                         newdata=newElev,link="log",quants=c(0.025,0.5,0.975),appendData=T)

names(pred.elev)[2:4] <- c("lower","median","upper") # rename quantiles

# Figure!
ggplot(pred.elev, aes(x=elev,y=mean,ymax=upper,ymin=lower)) +
  geom_line(size=1.4) + geom_ribbon(alpha=0.4, fill="gray40") +
  theme_classic()

##* Effect of forest on local abundance
# Create an object containing elevation values to predict detection (keeping the mean value for forest and length)
newForest <- data.frame(elev=mean(sitecovs$elev),
                      length=mean(sitecovs$iLength),
                      forest=seq(min(sitecovs$forest),max(sitecovs$forest),,100))

# Predicting...
pred.forest<- predictJAGS(fm=out,params=c("beta0","beta1","beta2","beta3"),
                        newdata=newForest,link="log",quants=c(0.025,0.5,0.975),appendData=T)

names(pred.forest)[2:4] <- c("lower","median","upper") # rename quantiles

# Figure!
ggplot(pred.forest, aes(x=forest,y=mean,ymax=upper,ymin=lower)) +
  geom_line(size=1.4,col="darkgreen") + geom_ribbon(alpha=0.4, fill="darkgreen") +
  theme_classic()
