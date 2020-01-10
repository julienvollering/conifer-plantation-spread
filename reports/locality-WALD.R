# setup #### 
speciesi <- "Larix" # adjust for RStudio jobs partitioning
localitiesi <- "voren" # adjust for RStudio jobs partitioning
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(here)
library(sf)
library(units)
source(here("R","locality-functions.R"))
source(here("R","dispersal-functions.R"))
source(here("R","SkarpaasShea_2007.R"))

# species dispersal traits ####
traits <- tribble(
  ~species, ~terminal.velocity, ~dispersal.season,
  "Larix", 1.0, c(12,1:5), # After Sandvik 2012 (for L. decidua), and Sullivan 1994 (for L. decidua)
  "P.abies", 0.58, c(11:12,1:5), # After Sandvik 2012 (but see Kaliniewicz et al. 2018), and Sullivan 1994
  "P.contorta", 0.82, c(9:12), # After Sandvik 2012
  "P.sitchensis-lutz", 0.94, c(10:12,1:2) # After Sandvik 2012 (but see Kaliniewicz et al. 2018), and Harris 1990
)

# localities table ####
loc.tab <- read_csv(here("data","localities.csv"))
loc.tab[is.na(loc.tab$age.at.registration),]
loc.tab[loc.tab$locality=="Håkøya", "age.at.registration"] <- loc.tab %>% 
  filter(species == "P.sitchensis-lutz", county %in% c("Troms","Nordland")) %>% 
  pull(age.at.registration) %>% mean(na.rm = TRUE) %>% round()
loc.tab[is.na(loc.tab$height.source),]

for (i in which(is.na(loc.tab$height.source))) {
  loc.tab[i, "height.source"] <- interpolate_height(loc.tab[i,], loc.tab)
}

# wind data ####
noraxy <- st_read(here("data","raw","NORA10","NORAsites.shp"))
nora <- read_csv(here("data","raw","NORA10","NORAwind.csv"))

eklimaxy <- st_read(here("data","raw","eKlima","EKLIMAsites.shp"))
eklima <- read_csv(here("data","raw","eKlima","EKLIMAwind.csv"))

length(intersect(nora$SITE, eklima$St.no))

sector.breaks <- seq(0, 360, by=20)

# vegetation height ####
grouping <- tribble(
  ~code, ~group, ~vegheight,
  "T1", "T1", 0,
  "T2", "T2", 0.5,
  "T3", "T3", 0.5,
  "T4", "T4", 10,
  "T6", "T6", 0,
  "T12", "T12", 0.5,
  "T13", "T13", 0,
  "T16", "T16", 0.5,
  "T17", "T17", 0,
  "T18", "T18", 0,
  "T21", "T21", 0,
  "T24", "T24", 0.5,
  "T27", "T27", 0,
  "T29", "T29", 0,
  "T30", "T30", 10,
  "T31", "T31", 0.5,
  "T32", "T32", 0.5,
  "T33", "T33", 0.5,
  "T34", "T34", 0.5,
  "T35", "disturbed", 0,
  "T36", "disturbed", 0.5,
  "T37", "disturbed", 0,
  "T38", "plantation forest", 10,
  "T39", "disturbed", 0,
  "T40", "disturbed", 0,
  "T41", "cultivated", 0.5,
  "T42", "disturbed", 0,
  "T43", "disturbed", 0,
  "T44", "cultivated", 0.5,
  "T45", "cultivated", 0.5,
  "T35.T37.T39.T43", "disturbed", 0,
  "T35.T37", "disturbed", 0,
  "V1", "fen", 0,
  "V2", "swamp", 10,
  "V3", "bog", 0,
  "V4", "spring", 0,
  "V8", "swamp", 10,
  "V9", "fen", 0,
  "V10", "fen", 0,
  "V11", "disturbed", 0,
  "V12", "disturbed", 0,
  "V13", "disturbed", 0) 

vegheights <- grouping %>% 
  group_by(group) %>% 
  summarize(vegheight = mean(vegheight))

# WALD seed shadows ####
species <- list.files(here("data","qc"))

for (i in speciesi) { 
  message(i)
  terminal.velocity <- filter(traits, species == i) %>% 
    pull(terminal.velocity)
  dispersal.season <- filter(traits, species == i) %>% 
    pull(dispersal.season) %>% 
    unlist()
  
  localities <- list.files(here("data","qc",i))
  
  for (j in localitiesi) { 
    message(j)
    nin <- st_read(here("data","qc",i,j,"nin.shp"), quiet = TRUE) %>% 
      as_tibble %>% st_as_sf() %>% 
      mutate(code = as.character(code), group = NULL) %>% 
      left_join(grouping, by = "code") %>% 
      filter(group != "plantation forest")
    grd <- make_grid(nin, gridres = 10)
    grd <- add_dummies_to_grid(grd, nin, field = group) %>% 
      as_tibble() %>% st_as_sf()
    grd <- grd %>% 
      mutate(rsum = rowSums(st_set_geometry(grd, NULL), na.rm = TRUE)) %>% 
      filter(rsum >= 0.5) %>% 
      select(-rsum) %>% 
      mutate_at(vars(-geometry), round, digits = 3) 
    vegheights.j <- vegheights %>% 
      filter(group %in% colnames(grd)) %>% 
      arrange(match(group, colnames(grd))) %>% 
      pull(vegheight)
    grd <- grd %>%
      mutate(normalizer = st_set_geometry(., NULL) %>% rowSums(., na.rm = TRUE))
    grd <- grd %>% 
      mutate(vegheight = grd %>% 
               st_set_geometry(NULL) %>%
               select(-normalizer) %>% 
               replace(is.na(.), 0) %>% 
               as.matrix() %*% vegheights.j %>% 
               as.vector())
    grd <- grd %>%
      mutate(vegheight = vegheight / normalizer) %>% 
      select(-normalizer)
    grd.pts <- st_centroid(st_geometry(grd))
    
    source <- st_read(here("data","qc",i,j,"source_polygon.shp"), quiet = TRUE) 
    source <- mutate(source, nsources = st_area(source) %>% 
                       `/`(100) %>% # 1 source per 100m2 (gridres 10x10m)
                       drop_units() %>% 
                       round())
    source.nonz <- filter(source, nsources > 0)
    set.seed(42)
    sources <- st_sample(source.nonz, source.nonz$nsources, type = "hexagonal") %>% 
      st_sf() %>% 
      st_transform(st_crs(grd))
    height <- loc.tab %>% 
      mutate(locality = tolower(locality)) %>% 
      filter(species == i, locality == j) %>% 
      pull(height.source)
    
    grdxy <- grd[1,] %>%  st_transform(st_crs(eklimaxy)) # EPSG 25833
    eklima.near <- eklimaxy[st_nearest_feature(grdxy, eklimaxy), ]
    nora.near <- noraxy[st_nearest_feature(grdxy, noraxy), ]
    if (st_distance(grdxy, eklima.near, by_element = TRUE) > set_units(10e3, "m") & 
        st_distance(grdxy, nora.near, by_element = TRUE) < set_units(10e3, "m")) {
      windobs <- filter(nora, SITE == nora.near$SITE)
      st_distance(grdxy, nora.near, by_element = TRUE) %>% 
        `/`(1e3) %>% 
        round() %>% 
        paste("km away, NORA") %>% 
        message()
    } else {
      windobs <- filter(eklima, St.no == eklima.near$Stnr)
      st_distance(grdxy, eklima.near, by_element = TRUE) %>% 
        `/`(1e3) %>% 
        round() %>% 
        paste("km away, eKlima") %>% 
        message()
    }
    windobs <- windobs %>% 
      select(Mnth, DD, FF) %>% 
      filter(Mnth %in% dispersal.season) %>% 
      mutate(sector = cut(DD, sector.breaks, include.lowest = TRUE))
    # Assume probability of abscission is UNRELATED to wind speed (constant)
    set.seed(42)
    sectors360 <- tibble(sector = unique(cut(0:360, sector.breaks, include.lowest = TRUE))) %>% 
      add_column(weight = sample_frequency(among = windobs, of = quo(sector), weightedby = NULL))
    set.seed(42)
    winds <- lapply(sectors360$sector, function(x) {
      filter(windobs, sector == x) %>% 
        sample_n(100, replace = TRUE, weight = NULL) %>% 
        pull(FF)
    })
    names(winds) <- sectors360$sector
    
    pb <- txtProgressBar(0, nrow(sources), style = 3)
    
    src.arr <- array(dim = c(nrow(grd), nrow(sources)))
    for (k in seq_along(sources$geometry)) { # k <- 1
      grd$sector <- assign_sectors(to = grd.pts, around = sources[k, ], 
                                   between = sector.breaks)
      distances.k <- as.numeric(st_distance(grd.pts, sources[k,]))
      sectors.k <- sectors360$sector[sectors360$sector %in% grd$sector]
      
      sec.arr <- array(dim = c(nrow(grd), length(sectors.k)))
      for (l in seq_along(sectors.k)) { # l <- 1
        vegheight <- grd %>% 
          filter(sector == sectors.k[l]) %>% 
          pull(vegheight) %>% 
          mean()
        distances.l <- distances.k
        distances.l[grd$sector != sectors.k[l]] <- NA
        
        wind.arr <- array(data = NA, dim = c(nrow(grd), 100))
        for (m in 1:100) {
          WALD <- parameterize_WALD(H = height, 
                                    Fm = terminal.velocity,
                                    Um = winds[[l]][m], 
                                    zm = 10, # 10 m above ground 
                                    h = vegheight)
          
          # Divisor below accounts for the fact that the WALD formulation used here (inverse gaussian dist) is the "dispersal DISTANCE kernel" (eq. 15.3 Nathan et al. 2012). 
          # 2πr reduces the probability density by spreading it across a 2D ring
          # nrow(sectors360) further reduces the probability density because we take only a fraction of the 2D ring (an arc)
          # Essentially, far cells are reduced relative to close cells, because the distance kernel is spread over a longer arc at distance.
          wind.arr[, m] <- WALD(distances.l) / (2*pi*distances.l*nrow(sectors360))
        }
        sec.arr[, l] <- sum_second_dim(wind.arr)
      }
      src.arr[, k] <- sum_second_dim(sec.arr, 
                                     sectors360$weight[which(sectors360$sector %in% sectors.k)])
      setTxtProgressBar(pb, k)
    }      
    close(pb)
    
    grd <- mutate(grd, wald = sum_second_dim(src.arr))
    ras <- grd %>% 
      select(wald) %>% 
      as("Spatial") %>% 
      raster::rasterize(raster::raster(grd, resolution = 10), "wald")
    raster::writeRaster(ras, here("data","qc",i,j,"wald.tif"), overwrite = TRUE)
  }
}
