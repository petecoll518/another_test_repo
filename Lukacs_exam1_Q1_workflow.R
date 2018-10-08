##### Question 1 (whales)

# Load packages
library(R2jags)
library(mcmcplots)
library(knitr)

# Define a function to generate the initial values and return them as a list
mde_ch_init <- function( ch, f ){
  for( i in 1:dim(ch)[1] ){
    ch[i,1:f[i]] <- NA
  }
  return(ch)
}

# a function to identify the time period of 
# the first capture from an encounter history matrix
mde_get_first <- function(x){
  
  return( min ( which( x != 0 ) ) )
}

# read in Whale data
load("whales.RData")

# Make encounter history matrix
EH <- whale$EH

nyears <- ncol( EH )

fall <- apply( EH, 1, mde_get_first ) # comes from mde_get_first.r file

# remove animals caught only in the last year 
length(which(fall == nyears)) # No indinviduals caught in last year only
f <- fall

nAnimal <- nrow( EH )

# Define latent state
z <- mde_ch_init( EH, f )   
z <- ifelse( !is.na(z), 1, z )

####################################################################################

# Hypothesis 1: phi and p are constant across time and sex.
# phi(.)p(.)

# Define a list of data to be passed to JAGS (in JAGS terminology)
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nyears )

### inital values

cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi", "b0.p", "mean.phi", "mean.p" )

# set up for MCMC run
ni <- 10000
nt <- 4
nb <- 5000
nc <- 3

# run the MCMC chain in JAGS. 
cjs.result.H1 <- jags( cjs.data, 
                       cjs.inits,
                       cjs.parms,
                       here::here("jags", "whales.phidot.pdot.txt"), 
                       n.iter=ni, 
                       n.burnin=nb,
                       n.thin=nt
)
cjs.result.H1
mcmcplot(cjs.result.H1)
H1.results <- (cjs.result.H1$BUGSoutput$summary)
H1.table <- kable(H1.results)

##############################################################################################

#### Hypothesis 2

# phi is lower in females than males, constant across time
# p is lower in females than males (due to accumulation of injuries)

# Define a sex vector
male <- whale$male

# Define a list of data to be passed to JAGS (in JAGS terminology)
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nyears, male = male )

### inital values

cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b1.phi = runif(1, -1, 1),
    b0.p = runif(1, -3, 3),
    b1.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi", "b1.phi", "b0.p", "b1.p",
               "mean.phi", "mean.p" )

# set up for MCMC run
ni <- 20000
nt <- 4
nb <- 10000
nc <- 3

# run the MCMC chain in JAGS. 
cjs.result.H2 <- jags( cjs.data, 
                       cjs.inits,
                       cjs.parms,
                       here::here("jags", "whales.phisex.psex.txt"), 
                       n.iter=ni, 
                       n.burnin=nb,
                       n.thin=nt
)
cjs.result.H2
mcmcplot(cjs.result.H2)

#################################################################################################

#### Hypothesis 3

# phi and p decrease with time due to increasing ship traffic and increasing
# misidentifications due to injury


# Define a list of data to be passed to JAGS (in JAGS terminology)
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nyears )

### inital values

cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b1.phi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3),
    b1.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi", "b1.phi", "b0.p", "b1.p",
               "mean.phi", "mean.p" )

# set up for MCMC run
ni <- 20000
nt <- 4
nb <- 10000
nc <- 3

# run the MCMC chain in JAGS. 
cjs.result.H3 <- jags( cjs.data, 
                       cjs.inits,
                       cjs.parms,
                       here::here("jags", "whales.phitime.ptime.txt"), 
                       n.iter=ni, 
                       n.burnin=nb,
                       n.thin=nt
)
cjs.result.H3
mcmcplot(cjs.result.H3)

########################################################################################

#### Hypothesis 4

# phi starts the same but then decreases with time for females only
# p is constant

# Define a sex vector
female <- whale$female

# Define a list of data to be passed to JAGS (in JAGS terminology)
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nyears, female = female )

### inital values

cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b1.phi = runif(1, -3, 3),
    b2.phi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3),
    b1.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi", "b1.phi", "b2.phi", "b0.p", 
               "mean.phi", "mean.p" )

# set up for MCMC run
ni <- 10000
nt <- 4
nb <- 5000
nc <- 3

# run the MCMC chain in JAGS. 
cjs.result.H4 <- jags( cjs.data, 
                       cjs.inits,
                       cjs.parms,
                       here::here("jags", "whales.phisextime.pdot.txt"), 
                       n.iter=ni, 
                       n.burnin=nb,
                       n.thin=nt
)
cjs.result.H4
mcmcplot(cjs.result.H4)

############################################################################################

# Hypothesis 5: phi is constant and p changes with time
# phi(.)p(time)

# Define a list of data to be passed to JAGS (in JAGS terminology)
cjs.data <- list( y=EH, f=f, nind=nAnimal, nocc=nyears )

### inital values

cjs.inits <- function(){
  list(
    z = z,
    b0.phi = runif(1, -3, 3),
    b0.p = runif(1, -3, 3)
  )
}

# set parameters to track in JAGS
cjs.parms <- c("b0.phi", "b0.p", "b1.p", "mean.phi", "mean.p" )

# set up for MCMC run
ni <- 10000
nt <- 4
nb <- 5000
nc <- 3

# run the MCMC chain in JAGS. 
cjs.result.H5 <- jags( cjs.data, 
                       cjs.inits,
                       cjs.parms,
                       here::here("jags", "whales.phidot.ptime.txt"), 
                       n.iter=ni, 
                       n.burnin=nb,
                       n.thin=nt
)
cjs.result.H5
mcmcplot(cjs.result.H5)