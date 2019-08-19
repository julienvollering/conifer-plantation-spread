library(tidyverse)
library(units)
library(sf)
library(here)

calc_density_by_field <- function(poly, pts, field) {
  
  field.enq <- enquo(field)
  poly <- poly %>%
    group_by(!!field.enq) %>%
    summarize()
  poly <- mutate(poly, area = st_area(poly))
  units(poly$area) <- with(ud_units, daa)

  pts <- st_transform(pts, st_crs(poly))
  browser()
  suppressWarnings({
    if (nrow(st_intersection(pts, poly)) > 0) {
      tally <- pts %>%
        st_intersection(poly) %>%
        group_by(!!field.enq) %>%
        summarize(n = n()) %>%
        st_set_geometry(NULL)
    } else {
      tally <- poly %>% 
        group_by(!!field.enq) %>%
        summarize(n = 0) %>%
        st_set_geometry(NULL)
    }
  })

  density <- left_join(poly, tally, by = quo_name(field.enq)) %>%
    st_set_geometry(NULL) %>%
    mutate(n = replace_na(n, 0)) %>%
    mutate(density = round(n/area, digits = 1)) %>%
    arrange(desc(density), desc(n), area)
  
  return(density)
}


make_grid <- function(poly, gridres) {
  
  datageom <- st_geometry(poly)
  bb <- st_bbox(datageom)
  grd <- st_sf(st_make_grid(datageom, 
                            cellsize = gridres, 
                            offset = floor(bb[c("xmin", "ymin")])))
  rename(grd, geometry = 1) 
}


# Note: grid must not contain column "ID"
add_dummies_to_grid <- function(grid, poly, field) {
  
  field.enq <- enquo(field)
  grid <- add_column(grid, ID = 1:nrow(grid))
  poly <- lwgeom::st_make_valid(poly) %>% 
    group_by(!!field.enq) %>% 
    summarize()
  
  suppressWarnings({
    int <- st_intersection(poly, grid) %>% 
      mutate(area = st_area(geometry))
  })
  
  groups <- pull(int, !!field.enq) %>% unique
  for (i in seq_along(groups)) {
    type <- groups[i]
    summ <- int %>% 
      st_set_geometry(NULL) %>% 
      filter(!!field.enq == type) %>%
      group_by(ID) %>% 
      summarize(!!type := sum(area))
    grid <- left_join(grid, summ, by = "ID")
  }
  
  grid <- grid %>% 
    select(-ID) %>% 
    mutate_at(groups, `/`, st_area(grid)) %>% 
    mutate_at(groups, as.numeric)
  
  rowS <- select(grid, groups) %>% 
    st_set_geometry(NULL) %>% 
    rowSums(na.rm = TRUE)
  fact <- ifelse(rowS > 1, 1/rowS, 1)
  grid <- mutate_at(grid, groups, function(x) {x * fact})
  
  return(grid)
}


# From Bullock et al. 2017
exp_power <- function(d, a = 2.7825, b = 0.8346) {
  b*(2*pi*a^2*gamma(2/b))^-1*exp(-(d^b/a^b))
}


add_seed_shadow <- function(grid, sourcepts) {
  suppressWarnings({
    centr <- st_centroid(grid) 
  })
  if (st_crs(sourcepts) != st_crs(centr)) {
    sourcepts <- st_transform(sourcepts, st_crs(centr))
  }
  src.array <- matrix(nrow = nrow(grid), ncol = nrow(sourcepts))
  for (i in 1:nrow(sourcepts)) {
    src.array[, i] <- st_distance(sourcepts[i,], centr) %>% 
      as.numeric() %>% 
      exp_power()
  }
  
  add_column(grid, seeds = rowSums(src.array, na.rm = TRUE), .before = "geometry")
}


add_relative_elevation <- function(grid, sourcepoly, dtm) {
  
  sourcepoly <- as(sourcepoly, "Spatial") %>% 
    sp::spTransform(CRSobj = raster::projection(dtm))
  elevref <- raster::extract(dtm, sourcepoly, fun = max) %>% 
    as.numeric() %>% 
    max()
  
  grd.proj <- grid %>% 
    as("Spatial") %>%   
    sp::spTransform(CRSobj = raster::projection(dtm))
  elev <- raster::extract(dtm, grd.proj, fun = mean) %>% 
    as.numeric
  grid %>% 
    add_column(relelev = elev, .before = "geometry") %>% 
    mutate(relelev = relelev - elevref)
}


# Note: grid must not contain column "ID"
add_wildlings <- function(grid, pts) {
  
  grid <- add_column(grid, ID = 1:nrow(grid))
  
  pts <- st_transform(pts, st_crs(grid))
  suppressWarnings({
    tally <- grid %>% 
      st_intersection(pts) %>% 
      st_set_geometry(NULL) %>% 
      group_by(ID) %>% 
      summarize(wildlings = n())
  })  
  
  grid %>% 
    left_join(tally, by = "ID") %>% 
    select(-ID) %>% 
    select(-geometry, geometry)
}