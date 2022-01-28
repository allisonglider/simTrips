require(track2KBA)
library(crawl)
library(momentuHMM)
library(rgdal)
library(sp)
library(dplyr)
library(ggplot2)
theme_set(theme_light())

#1897273090

dataset <- move2KBA(movebankID = 	1895716931,
                    user = rstudioapi::askForPassword(prompt = 'Movebank username:'),
                    password = rstudioapi::askForPassword(prompt = 'Movebank password:'),
                    filename = NULL)

tracks <- dataset$data
colony <- dataset$site

trips <- tripSplit(
  dataGroup  = tracks, # data formatted using formatFields()
  colony     = colony, # data on colony location - can be extracted from movebank data using move2KBA()
  innerBuff  = 5,      # (km) minimum distance from the colony to be in a trip
  returnBuff = 10,     # (km) outer buffer to capture incomplete return trips
  duration   = 6,      # (hrs) minimum trip duration
  rmNonTrip  = TRUE    # T/F removes times when not in trips
)

trips <- subset(trips, trips$Returns == 'Yes')

sumTrips <- tripSummary(trips = trips, colony = colony)
sumTrips

# -----

data <- data.frame(
  ID = trips$tripID,
  time = as.POSIXct(trips$DateTime),
  lon = trips$Longitude,
  lat =trips$Latitude
)
rawData <- data

utmcoord <- sp::spTransform(trips, CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))

rawData$x <- attr(utmcoord,"coords")[,1]
rawData$y <- attr(utmcoord,"coords")[,2]

crwOut <- crawlWrap(obsData=rawData, timeStep="20 min",
                    theta=c(4, 0), fixPar=c(NA,NA))

pred <- data.frame(crwOut$crwPredict)
predcoord <- sp::SpatialPoints(pred[pred$locType == 'p',c('mu.x','mu.y')],
                               proj4string=CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))
mapview::mapview(predcoord)

center <- as.matrix(colony, by.row = F)
center <- sp::SpatialPoints(center,
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
center <- sp::spTransform(center,CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))
center <- as.matrix(data.frame(coordinates(center)))

data <- prepData(data=crwOut, type="UTM", centers=center)

head(data)

data <- data %>% 
  group_by(ID) %>% 
  mutate(
    dt = as.numeric(difftime(time, lag(time, 1), units = 'hours')),
    dt = ifelse(is.na(dt), 0, dt),
    cum_time = cumsum(dt),
    d_dist = center1.dist - lag(center1.dist, 1)
  )

ggplot(data, aes(y = d_dist, x = cum_time)) +
  geom_point()

center <- as.matrix(colony, by.row = F)
col_loc <- sp::SpatialPoints(center,
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
col_loc <- sp::spTransform(col_loc,CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))


samps <- 1000
steps <- data[!is.na(data$angle),c('step', 'angle')]
dt <- median(data$dt)

k <- 0
if (k == 0) {
  out <- data.frame(time = 0, x = coordinates(col_loc)[,1], y = coordinates(col_loc)[,2],
                    dist = 0, d_dist = 0)
  k <- k + 1
}

obs_dist <- na.omit(data$d_dist[data$ID == unique(data$ID)[3]])

while (k < 200) {
  
  dd_dist <- density(obs_dist[k:(k + 15)], from = -50000, to =50000)
  dd_dist$y <- dd_dist$y/sum(dd_dist$y)
  dist_fun <- approxfun(x = dd_dist$x, y = dd_dist$y)
  
  # if (k > length(obs_dist) - 10) {
  #   mean_dd <- -3000
  #   sd_dd <- 1000 
  # } else {
  #   mean_dd <- 0
  #   sd_dd <- 3000
  #   # mean_dd <- mean(obs_dist[k:(k + 5)], na.rm = T)
  #   # sd_dd <- sd(obs_dist[k:(k + 5)], na.rm = T)
  # }
  # 
  # if (k < 20) mean_dd <- 3000
  # if (k < 20) sd_dd <- 1000
  
  part_idx <- sample(1:nrow(steps), samps)
  tt <- steps[part_idx,]
  tt$x <- out$x[k] + (tt$step * cos(tt$angle * 180/pi))
  tt$y <- out$y[k] + (tt$step * sin(tt$angle * 180/pi))
  pp <- SpatialPoints(tt[, c('x','y')], proj4string = CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))
  tt$dist <- raster::pointDistance(col_loc, pp)
  tt$d_dist <- tt$dist - out$dist[k]
  
  
  # p_grid <- seq(from = -20000, to = 20000)
  # liklihood <- dnorm(tt$y,mean = mean_dd, sd_dd)
  # unstd.posterior <- liklihood
  posterior <- dist_fun(tt$d_dist)
  samp_idx <- sample(1:samps, 1, prob = posterior)
  tt[samp_idx,]

  temp <- data.frame(time = out$time[k] + dt, x = tt$x[samp_idx], y = tt$y[samp_idx],
                     dist = tt$dist[samp_idx], d_dist = tt$d_dist[samp_idx])
  out <-rbind(out, temp)
  
  k <- k + 1
}

ss <- SpatialPoints(out[, c('x','y')], proj4string = CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))

mapview::mapview(ss) + mapview::mapview(col_loc, col.regions = 'red') 

ggplot() +
  geom_sf(data = as(col_loc, 'sf')) +
  geom_path(data = out, aes(x = x, y = y))

ggplot(out, aes(x = time, y = dist)) +
  geom_line()

# ------
library(adehabitatLT)

bath <- raster::raster('E:/ECCC_OPP/bathymetry/bathymetry.tiff')

origin_coord <- data.frame(x = -131.21, y= 52.349)
origin_col <- as.matrix(origin_coord, by.row = F)
origin_col <- sp::SpatialPoints(origin_col, proj4string=CRS("+proj=longlat +datum=WGS84"))
origin_col <- sp::spTransform(origin_col,CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))

sim_coord <- data.frame(x = -133.179, y= 53.936)
sim_col <- as.matrix(sim_coord, by.row = F)
sim_col <- sp::SpatialPoints(sim_col, proj4string=CRS("+proj=longlat +datum=WGS84"))
sim_col <- sp::spTransform(sim_col,CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))

water <- bath
water <- raster::projectRaster(water, crs = origin_col)
raster::values(water) <- ifelse(raster::values(water) <0, 1, 0)

dist_origin <- raster::raster(res = 500, ext = raster::extent(water), crs = water)
dist_origin <- raster::distanceFromPoints(dist_origin, origin_col)

dist_sim <- raster::raster(res = 500, ext = raster::extent(water), crs = water)
dist_sim <- raster::distanceFromPoints(dist_sim, sim_col)

new_water <- raster::resample(x = water, y = dist_sim, method="bilinear")
raster::values(new_water) <- ifelse(raster::values(new_water) > 0, 1, 0)
plot(new_water)

plot(new_water)
water <- new_water
# ------

myCS = suppressWarnings(gdistance::transition(new_water, mean, directions=16))
myCS = suppressWarnings(gdistance::geoCorrection(myCS))

dist_origin <- suppressWarnings(gdistance::accCost(myCS, fromCoords = origin_col))
raster::values(dist_origin)[is.infinite(raster::values(dist_origin))] <- NA

dist_sim <- suppressWarnings(gdistance::accCost(myCS, fromCoords = sim_col))
raster::values(dist_sim)[is.infinite(raster::values(dist_sim))] <- NA

water <- dist_sim
raster::values(water) <- ifelse(is.na(raster::values(water)), 0, 1)

# -----

lt <- as.ltraj(xy = data[,c('x','y')], date = data$time, id = factor(data$ID))
lt_dat <- ld(lt) %>% 
  dplyr::select(id,date,x,y, dist,abs.angle,rel.angle,dt) %>% 
  mutate(
    #col.dist = sqrt((x - sp::coordinates(col_loc)[,1])^2 + (y - sp::coordinates(col_loc)[,2])^2),
    col.dist = raster::extract(dist_origin, data.frame(x,y)),
    d.dist = col.dist - lag(col.dist)
  )

samps <- 1000
dt <- median(lt_dat$dt, na.rm = T)

start_loc <- sample(which(is.na(lt_dat$dt)), 1)
sample_trip <- sample(unique(lt_dat$id), 1)
sample_trip <- subset(lt_dat, lt_dat$id == sample_trip) 
sample_trip <- subset(sample_trip, !is.na(sample_trip$abs.angle)) 

sim_id <- 'sim_1'
sim_data <- data.frame(
  id = sim_id,
  date = seq(from = 0, by = dt, length.out = nrow(sample_trip) - 10),
  x = NA, y = NA, dist = NA, abs.angle = NA, rel.angle = NA, 
  dt = dt, col.dist = NA, d.dist = NA
)

k <- 1
if (k == 1) {
  start_dist <- rnorm(samps, mean = 5000, 500)
  start_angle <- runif(samps, min = -3, max = 3)
  start_x <- coordinates(sim_col )[,1] + (start_dist * cos(start_angle))
  start_y <- coordinates(sim_col )[,2] + (start_dist * sin(start_angle))
  start_water <- raster::extract(water, data.frame(start_x, start_y))
  start_col.dist <- raster::extract(dist_sim, data.frame(start_x, start_y))
  start_col.dist[is.na(start_col.dist)] <- Inf
  
  samp_idx <- sample(1:samps, 1, prob = start_water)
  sim_data$x[k] <- start_x[samp_idx]
  sim_data$y[k] <- start_y[samp_idx]
  sim_data$dist[k] <- start_dist[samp_idx]
  sim_data$abs.angle[k] <- start_angle[samp_idx]
  sim_data$col.dist[k] <- start_col.dist[samp_idx]
  
  # sim_data$x[k] <- sample_trip$x[k]
  # sim_data$y[k] <- sample_trip$y[k]
  # sim_data$abs.angle[k] <- sample_trip$abs.angle[k]
  # sim_data$col.dist[k] <- sample_trip$col.dist[k]
  k <- k + 1
}

steps <- lt_dat %>% filter(!is.na(rel.angle))
trip_pos <- round(1:nrow(sample_trip)/nrow(sample_trip) * 10)/10

while (k <= nrow(sim_data)) {
  
  pts_idx <- sample(1:nrow(steps), samps, replace = T)
  samp_dist <- steps$dist[pts_idx]
  samp_rel.angle <- steps$rel.angle[pts_idx]
  samp_abs.angle <- samp_rel.angle + sim_data$abs.angle[k - 1]
  samp_x <- sim_data$x[k - 1] + (samp_dist * cos(samp_abs.angle))
  samp_y <- sim_data$y[k - 1] + (samp_dist * sin(samp_abs.angle))
  #samp_col.dist <- sqrt((samp_x - sp::coordinates(col_loc)[,1])^2 + (samp_y - sp::coordinates(col_loc)[,2])^2)
  samp_col.dist <- raster::extract(dist_sim , data.frame(samp_x, samp_y))
  samp_col.dist[is.na(samp_col.dist)] <- Inf
  samp_d.dist <- samp_col.dist - sim_data$col.dist[k - 1]
    
  if (trip_pos[k] < 0.9) {
    dd_dist <- density(sample_trip$col.dist[trip_pos == trip_pos[k]], from = 0, to =max(steps$col.dist))
    dd_dist$y <- dd_dist$y/sum(dd_dist$y)
    dist_fun <- approxfun(x = dd_dist$x, y = dd_dist$y, rule = 2)
    posterior <- dist_fun(samp_col.dist)
  } else {
    dd_dist <- density(sample_trip$d.dist[trip_pos == trip_pos[k]], from = -20000, to =0)
    dd_dist$y <- dd_dist$y/sum(dd_dist$y)
    dist_fun <- approxfun(x = dd_dist$x, y = dd_dist$y, yright = 0, yleft = 0)
    posterior <- dist_fun(samp_d.dist)
    
    p_grid <- seq(from = -20000, to = 20000)
    liklihood <- dnorm(samp_d.dist, mean = -5000, sd = 1000)
    posterior <- liklihood/sum(liklihood)
  }
  
  water_probs <- raster::extract(water, data.frame(samp_x, samp_y))
  posterior <- posterior * water_probs
  
  samp_idx <- sample(1:samps, 1, prob = posterior)
    
  sim_data$x[k] <- samp_x[samp_idx]
  sim_data$y[k] <- samp_y[samp_idx]
  sim_data$dist[k] <- samp_dist[samp_idx]
  sim_data$abs.angle[k] <- samp_abs.angle[samp_idx]
  sim_data$rel.angle[k] <- samp_rel.angle[samp_idx]
  sim_data$col.dist[k] <- samp_col.dist[samp_idx]
  sim_data$d.dist[k] <- samp_d.dist[samp_idx]
  
  k <- k + 1
}

ss <- SpatialPoints(sim_data[1:(k-1), c('x','y')], proj4string = CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))
oo <- SpatialPoints(sample_trip[1:(k-1), c('x','y')], proj4string = CRS("+proj=utm +zone=8 +ellps=WGS84 +units=m +no_defs "))

mapview::mapview(oo, col.regions = 'black') + 
  mapview::mapview(ss) + 
  mapview::mapview(origin_col , col.regions = 'red')  + 
  mapview::mapview(sim_col , col.regions = 'red') 


ggplot(sim_data, aes(x = date, y = col.dist)) +
  geom_line() 

ggplot(sample_trip, aes(x = date, y = col.dist)) +
  geom_line() 
