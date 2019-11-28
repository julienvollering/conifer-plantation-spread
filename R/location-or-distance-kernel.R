# Test params
r <- 10
a <- c(1:10)
b <- c(1:10)

# From Nathan et al. 2012
# Wald dispersal *location* kernel
kLr <- function(r, a, b) {
  (sqrt(b) / sqrt(8*pi^3*r^5)) * exp(-(b*(r - a)^2) / (2*a^2*r))
}

kLr(r, a, b)

# Dispersal location kernal is inverse gaussion distribution divided by 2*pi*r (eq. 15.3 in Nathan 2012)
statmod::dinvgauss(r, mean=a, shape=b)/(2*pi*r)
