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

assign_sectors <- function(to, around, between) {
  # compass bearing FROM "around" TO "to (matches wind direction DD)
  dX <- st_coordinates(around)[,1] - st_coordinates(to)[,1]
  dY <- st_coordinates(around)[,2] - st_coordinates(to)[,2]
  sectors <- cut(angle(dX, dY), between, include.lowest = TRUE)
  sectors[is.na(sectors)] <- levels(sectors)[1]
  return(sectors)
}

sum_third_dim <- function(array, weights = NULL) {
  if (!is.null(weights)) {
    for (i in 1:dim(array)[3]) {
      array[, , i] <- array[, , i] * weights[i]
    }
  }
  rowSums(array, na.rm = TRUE, dims = 2)
}

sum_second_dim <- function(array, weights = NULL) {
  if (!is.null(weights)) {
    for (i in 1:dim(array)[2]) {
      array[, i] <- array[, i] * weights[i]
    }
  }
  rowSums(array, na.rm = TRUE)
}
