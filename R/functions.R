library(tidyverse)
library(units)
library(sf)
library(here)

interpolate_height <- function(locality, localities) {
  sp <- pull(locality, species)
  age <- pull(locality, age.at.registration)
  localities <- filter(localities, species == sp) %>% 
    select(height.source, age.at.registration) %>% 
    na.omit() %>% 
    mutate(weight = ifelse(age == age.at.registration, 2,
                           1/abs(age.at.registration - age)))
  weighted.mean(localities$height.source, localities$weight) %>% round()
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
  poly <- sf::st_make_valid(poly) %>%
    mutate(!!field.enq := as.character(!!field.enq)) %>% 
    group_by(!!field.enq) %>% 
    summarize(.groups = "drop") # geometry union
  suppressWarnings({
    int <- st_intersection(poly, grid) %>% 
      mutate(area = st_area(geometry))
  })
  
  groups <- pull(int, !!field.enq) %>% unique
  for (i in seq_along(groups)) {
    name <- groups[i]
    summ <- int %>% 
      st_set_geometry(NULL) %>% 
      filter(!!field.enq == name) %>%
      group_by(ID) %>% 
      summarize(!!name := sum(area), .groups = "drop_last")
    grid <- left_join(grid, summ, by = "ID")
  }
  
  grid <- grid %>% 
    select(-ID) %>% 
    mutate_at(all_of(groups), `/`, st_area(grid)) %>% 
    mutate_at(all_of(groups), as.numeric)
  
  rowS <- select(grid, all_of(groups)) %>% 
    st_set_geometry(NULL) %>% 
    rowSums(na.rm = TRUE)
  fact <- ifelse(rowS > 1, 1/rowS, 1)
  grid <- grid %>% 
    mutate_at(all_of(groups), function(x) {x * fact}) %>% 
    mutate_at(all_of(groups), ~ replace_na(., replace = 0)) %>% 
    mutate_at(all_of(groups), function(x) {attributes(x) <- NULL; x})
  
  return(grid)
}


# From Bullock et al. 2017
ExP <- function(d, a = 2.7825, b = 0.8346) {
  b*(2*pi*a^2*gamma(2/b))^-1*exp(-(d^b/a^b))
}


add_seeds_ExP <- function(grid, sourcepts) {
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
      ExP()
  }
  
  add_column(grid, seeds.ExP = rowSums(src.array, na.rm = TRUE), .before = "geometry")
}


add_seeds_WALD <- function(grid, raster) {
  grd.proj <- sf::st_transform(grid, crs = sf::st_crs(raster)) 
  WALD <- exactextractr::exact_extract(raster, grd.proj, fun = 'mean', progress = FALSE)
  add_column(grid, seeds.WALD = WALD, .before = "geometry")
}


add_relative_elevation <- function(grid, sourcepoly, dtm) {
  
  sourcepoly.proj <- sf::st_transform(sourcepoly, crs = sf::st_crs(dtm))
  elevref <- exactextractr::exact_extract(dtm, sourcepoly.proj, fun = 'max', 
                                          progress = FALSE) %>% 
    max()
  
  grd.proj <- sf::st_transform(grid, crs = sf::st_crs(dtm)) 
  elev <- exactextractr::exact_extract(dtm, grd.proj, fun = 'mean', progress = FALSE)
  add_column(grid, relelev = elev, .before = "geometry") %>% 
    mutate(relelev = relelev - elevref)
}


# Note: grid must not contain column "ID"
add_wildlings <- function(grid, pts) {
  if(nrow(pts) == 0) {
    out <- grid %>% 
      mutate(wildlings = 0, .before = locality)
  } else {
  idx <- st_geometry(grid) %>% 
    st_intersection(st_geometry(pts)) %>% 
    attr(., "idx")
  tally <- table(idx[,1]) %>% 
    as.data.frame(responseName = "wildlings", stringsAsFactors = FALSE) %>% 
    mutate(Var1 = as.integer(Var1))
  out <- grid %>% 
    mutate(Var1 = 1:nrow(grid), .before = 1) %>% 
    mutate(wildlings = 0, .before = locality) %>% 
    as_tibble() %>% 
    rows_update(tally, by = "Var1") %>%
    select(-Var1) %>% 
    st_as_sf()
  }
  return(out)
}

remove_zero_cols <- function(sf, na.rm = TRUE) {
  df <- st_drop_geometry(sf)
  nonnum <- select(df, !where(is.numeric)) %>% 
    names()
  num <- select(df, where(is.numeric)) %>% 
    colSums(na.rm = na.rm) %>% 
    enframe() %>% 
    filter(value != 0) %>% 
    pull(name)
  return(select(sf, all_of(c(nonnum, num))))
}

identify_completeseparation <- function(dat, typenames) {
  pivot_longer(dat, one_of(typenames), names_to = "nin", values_to = "are") %>%
    group_by(nin) %>% 
    summarize(n = sum(wildlings*are), .groups = "drop_last") %>% 
    filter(n == 0) %>% 
    pull(nin)
}

paste_modelformula <- function(dat, typenames, contrasttype) {
  cs <- identify_completeseparation(dat, typenames)
  mtypes <- typenames[!(typenames %in% c(cs, contrasttype))]
  fstring <- paste("wildlings ~ age + bio01 + bio19 + relelev +",
                   paste(mtypes, collapse = " + "),
                   "+ (1 | locality)")
  return(formula(fstring))
}

fit_genpoisWALDmodel <- function(dat, typenames, contrasttype) {
  formula <- paste_modelformula(dat, typenames, contrasttype)
  mod <- glmmTMB(formula, dat, offset = seeds.WALD, family = "genpois", 
                 ziformula = ~ age + (1 | locality))
  return(mod)
}

calculate_wildlings_per_nin <- function(dat, typenames, contrasttype) {
  dens <- pivot_longer(dat, one_of(typenames), names_to = "nin", values_to = "are") %>%
    group_by(nin) %>% 
    summarize(daa = sum(are)/10, # units from are to decare
              n = sum(wildlings*are), .groups = "drop_last") %>% 
    mutate(density = n/daa)
  scaler <- filter(dens, nin == contrasttype) %>% 
    pull(density)
  dens <- mutate(dens, fit = density/scaler)
  return(dens)
}
