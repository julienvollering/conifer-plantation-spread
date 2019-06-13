angle <- function(dx,dy) #JV: problem with this function: gives 0 deg for N, and 90 deg for E --- should be 180 and 270 deg
{
  a <- atan(abs(dx/dy))*180/pi
  a[dx>0 & dy<0] <- 180-a[dx>0 & dy<0]
  a[dx<0 & dy<0] <- 180+a[dx<0 & dy<0]
  a[dx<0 & dy>0] <- 360-a[dx<0 & dy>0]
  a[dx==0 & dy<0] <- 180
  a[dx<0 & dy==0] <- 270
  return(a)
}