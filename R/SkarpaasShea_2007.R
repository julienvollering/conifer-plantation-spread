# From Skarpaas, O. and Shea, K. 2007. Dispersal Patterns, Dispersal Mechanisms, and Invasion Wave Speeds for Invasive Thistles. - The American Naturalist 170: 421–430.

#Eq. 6:
# H seed release height
# U mean horizontal wind velocity between H and the ground
# Fm seed terminal velocity
calc_nu <- function(H, U, Fm) {H*U/Fm}

#Eq 7:
# H seed release height
# sigma turbulent flow parameter
calc_lambda <- function(H, sigma) {(H/sigma)^2}

# Eq. A1: (lower bound of integration changed to avoid undefined values)
# H seed release height
# ustar friction velocity
# d surface roughness parameter d
# z0 surface roughness parameter z0
calc_U <- function(H, ustar, d, z0) {
  (1/H)*integrate(function(z, ustar, d, z0) {(ustar/0.4)*log((z-d)/z0)}, 
                  lower=z0+d, upper=H, # lower limit set to z0+d to avoid log(<1)
                  ustar=ustar, d=d, z0=z0)$value 
}

# Eq. A2:
# Um mean wind speed at measurement height
# zm wind speed measurement height
# d surface roughness parameter d
# z0 surface roughness parameter z0
calc_ustar <- function(Um, zm, d, z0) {
  0.4*Um/log((zm-d)/z0)
}

# Eq. A4:
# z height above ground
# d surface roughness parameter d
# ustar friction velocity
# U mean horizontal wind velocity between H and the ground
calc_sigma <- function(z, d, ustar, U) {
  (2*1.3^2)*sqrt((0.4*(z-d)*ustar)/(3.125*U))
}

# h vegetation height
calc_d <- function(h) {
  0.7*h + 0.01 # Avoids NaN in calc_U integration
}

# h vegetation height
calc_z0 <- function(h) {
  0.1*h + 0.01 # Avoids NaN in calc_U integration
}

parameterize_WALD <- function(H, Fm, Um, zm, h) {
  if (h > zm) { zm <- h } # Avoids unrealistic estimates of friction velocity (ustar)
  d <- calc_d(h)
  z0 <- calc_z0(h)
  ustar <- calc_ustar(Um, zm, d, z0)
  if (h >= 0.8*H) { # Avoids wind speeds below friction velocity (ustar), else they decrease to zero
    U <- ustar 
  } else {
    U <- calc_U(H, ustar, d, z0)
  }
  sigma <- calc_sigma(z = max(H, h), d, ustar, U) # Avoids sqrt(<0), z = turbulent flow height above ground
  nu <- calc_nu(H, U, Fm)
  lambda <- calc_lambda(H, sigma)
  function(x) {statmod::dinvgauss(x, mean=nu, shape=lambda)}
}
