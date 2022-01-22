

create_raster_inputs <- function(obs_center, sim_center = NULL, max_range,
                                 bathymetry,
                                 analysis_crs,
                                 raster_res = 500) {

  if (is.null(sim_center)) sim_center <- obs_center

  ctrs <- data.frame(
    x = c(obs_center[1], sim_center[1]),
    y = c(obs_center[2], sim_center[2])
  )

  out <- vector(mode = 'list', length = 2)
  names(out) <- c('origin', 'sim')

  for (i in 1:2) {

    if (i == 2 & identical(sim_center, obs_center)) {
      s <- raster::stack(water, dist_origin)
      names(s) <- c('water','dist')
      out[[i]] <- s

      raster::plot(dist_origin, main = "Over water distance from simulation location (m)")
      points(origin_col, col = 'black', pch = 19)
      par(mfrow = c(1,1))

    } else {

      origin_col <- sp::SpatialPoints(as.matrix(ctrs[i,], by.row = F), proj4string=raster::crs('+proj=longlat'))
      origin_col <- sp::spTransform(origin_col, CRSobj = sp::CRS(analysis_crs))
      origin_buff <- raster::buffer(origin_col, width = max_range * 2)

      bathymetry <- bathymetry
      bathymetry <- raster::projectRaster(bathymetry, crs = origin_col)

      water <- raster::crop(bathymetry, origin_buff)

      dist_origin <- raster::raster(res = raster_res, ext = raster::extent(water), crs = water)
      dist_origin <- raster::distanceFromPoints(dist_origin, origin_col)

      water <- raster::resample(x = water, y = dist_origin, method="bilinear")
      raster::values(water) <- ifelse(raster::values(water) < 0, 1, 0)

      myCS <- suppressWarnings(gdistance::transition(water, mean, directions=16))
      myCS <- suppressWarnings(gdistance::geoCorrection(myCS))

      dist_origin <- suppressWarnings(gdistance::accCost(myCS, fromCoords = origin_col))
      #raster::values(dist_origin)[is.infinite(raster::values(dist_origin))] <- NA

      s <- raster::stack(water, dist_origin)
      names(s) <- c('water','dist')
      out[[i]] <- s

      raster::plot(dist_origin, main = "Over water distance from original location (m)")
      points(origin_col, col = 'black', pch = 19)
      par(mfrow = c(1,1))


    }
  }

  out
  #' @export
}
