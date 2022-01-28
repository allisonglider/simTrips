#' Creates raster stack of coastline, overwater distance, and bathymetry
#'
#' @param obs_center Coordinates (decimal degrees) of central place, from which
#' over water distances will be calculated
#' @param max_range Maximum distance (m) from the central place the species is
#' expected to travel.
#' @param bathymetry raster where depth values should be meters below 0,
#' used to define coastline
#' @param analysis_crs a proj4string for the projection of raster outputs, should have units in m

create_raster_inputs <- function(obs_center,
                                 max_range,
                                 bathymetry,
                                 analysis_crs,
                                 raster_res = 500) {


  ctrs <- data.frame(x = obs_center[1], y = obs_center[2])

  origin_col <- sp::SpatialPoints(as.matrix(ctrs, by.row = F), proj4string=raster::crs('+proj=longlat'))
  origin_col <- sp::spTransform(origin_col, CRSobj = sp::CRS(analysis_crs))
  origin_buff <- raster::buffer(origin_col, width = max_range * 2)

  bathymetry <- bathymetry
  bathymetry <- raster::projectRaster(bathymetry, crs = origin_col)

  water <- raster::crop(bathymetry, origin_buff)

  dist_origin <- raster::raster(res = raster_res, ext = raster::extent(water), crs = water)

  water <- raster::resample(x = water, y = dist_origin, method="bilinear")
  raster::values(water) <- ifelse(raster::values(water) < 0, 1, 0)
  dist_origin <- raster::distanceFromPoints(dist_origin, origin_col)

  # Buffer raster_res*4 around origin and set water ==1
  col_area <- raster::cellFromPolygon(water, raster::buffer(origin_col, raster_res*4))[[1]]
  raster::values(water)[col_area] <- 1

  myCS <- suppressWarnings(gdistance::transition(water, mean, directions=16))
  myCS <- suppressWarnings(gdistance::geoCorrection(myCS))

  dist_origin <- suppressWarnings(gdistance::accCost(myCS, fromCoords = origin_col))
  #raster::values(dist_origin)[is.infinite(raster::values(dist_origin))] <- NA

  bath <- raster::resample(bathymetry, water)
  raster::values(bath)[raster::values(bath) > 0] <- NA
  s <- raster::stack(water, dist_origin, bath)
  names(s) <- c('water','dist','bathymetry')

  raster::plot(dist_origin, main = "Over water distance from original location (m)")
  points(origin_col, col = 'black', pch = 19)
  par(mfrow = c(1,1))

  s
  #' @export
}
