# From Skarpaas, O. and Shea, K. 2007. Dispersal Patterns, Dispersal Mechanisms, and Invasion Wave Speeds for Invasive Thistles. - The American Naturalist 170: 421–430.

#Eq. 6 and 7:
nu <- function(H, U, Fm) {H*U/Fm}
lambda <- function(H, sigma) {(H/sigma)^2}

# Eq. A1: (note that the lower bound of integration was changed to match the conifer dispersal functions.)
integral <- function(z, ustar, K=0.4, d, z0) {
  (ustar/K)*log((z-d)/z0)
}
U <- function(H, ustar, K=0.4, d, z0) {
  (1/H)*integrate(integral, 
                  lower=z0+d, upper=H, # lower limit set to z0+dhat to avoid log(0) and log(negative)
                  ustar=ustar, d=d, z0=z0)$value 
}

# Eq. A2:
ustar <- function(K=0.4, Um, z, d, z0) {
  K*Um/log((z-d)/z0)
}

# Eq. A4:
sigma <- function(Aw=1.3, K=0.4, z, d, ustar, C0=3.125, U) {
  (2*Aw^2)*sqrt((K*(z-d)*ustar)/(C0*U))
}
