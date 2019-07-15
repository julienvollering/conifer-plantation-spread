angle <- function(dx,dy)
{
  a <- atan(abs(dx/dy))*180/pi
  a[dx>0 & dy<0] <- 180-a[dx>0 & dy<0]
  a[dx<0 & dy<0] <- 180+a[dx<0 & dy<0]
  a[dx<0 & dy>0] <- 360-a[dx<0 & dy>0]
  a[dx==0 & dy<0] <- 180
  a[dx<0 & dy==0] <- 270
  return(a)
}

matimage <- function(matrix) {
  require(tidyverse)
  matrix %>%
    apply(2, rev) %>% 
    t() %>% 
    image()
}
