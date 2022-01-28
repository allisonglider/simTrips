
#' Uses a particle filter to simulate trips based on observed movement rates,
#' conditioned on observed PDF of depth and distance from colony
#'

sim_trips <- function(data,
                      origin_col,
                      sim_col,
                      dist_origin,
                      dist_sim,
                      water_origin,
                      water_sim,
                      bath_origin,
                      bath_sim,
                      n_trips = 50,
                      samps = 1000) {

  lt <- adehabitatLT::as.ltraj(xy = data[,c('x','y')], date = data$time, id = factor(data$ID))
  lt_dat <- adehabitatLT::ld(lt) %>%
    dplyr::select(id,date,x,y, dist,abs.angle,rel.angle,dt) %>%
    dplyr::group_by(id) %>%
    dplyr::mutate(
      col.dist = raster::extract(dist_origin, data.frame(x,y)),
      bath = raster::extract(bath_origin, data.frame(x,y)),
      d.dist = col.dist - lag(col.dist),
      n = dplyr::n()
    )

  dt <- median(lt_dat$dt, na.rm = T)

  out <- data.frame()

  steps <- lt_dat %>%
    dplyr::filter(!is.na(dt)) %>%
    dplyr::mutate(rel.angle = ifelse(is.na(rel.angle), 0, rel.angle))

  for (i in 1:n_trips) {

    start_loc <- sample(which(is.na(lt_dat$dt)), 1)
    sample_trip <- sample(unique(lt_dat$id), 1)
    sample_trip <- subset(lt_dat, lt_dat$id == sample_trip)
    sample_trip <- subset(sample_trip, !is.na(sample_trip$abs.angle))
    trip_pos <- round(1:nrow(sample_trip)/nrow(sample_trip) * 10)/10

    sim_id <- i
    sim_data <- data.frame(
      id = sim_id,
      date = seq(from = 0, by = dt, length.out = nrow(sample_trip)),
      x = NA, y = NA, dist = NA, abs.angle = NA, rel.angle = NA,
      dt = dt, col.dist = NA, bath = NA, d.dist = NA
    )

    k <- 1
    if (k == 1) {
      start_dist <- rnorm(samps, mean = 5000, 500)
      start_angle <- runif(samps, min = -3, max = 3)
      start_x <- sp::coordinates(sim_col)[,1] + (start_dist * cos(start_angle))
      start_y <- sp::coordinates(sim_col)[,2] + (start_dist * sin(start_angle))
      start_water <- raster::extract(water_sim, data.frame(start_x, start_y))
      start_col.dist <- raster::extract(dist_sim, data.frame(start_x, start_y))
      start_col.dist[is.na(start_col.dist)] <- Inf
      start_bath <- raster::extract(bath_sim, data.frame(start_x, start_y))
      start_bath[is.na(start_bath)] <- Inf

      samp_idx <- sample(1:samps, 1, prob = start_water)
      sim_data$x[k] <- start_x[samp_idx]
      sim_data$y[k] <- start_y[samp_idx]
      sim_data$dist[k] <- start_dist[samp_idx]
      sim_data$abs.angle[k] <- start_angle[samp_idx]
      sim_data$col.dist[k] <- start_col.dist[samp_idx]
      sim_data$bath[k] <- start_bath[samp_idx]

      k <- k + 1
    }


    while (k <= nrow(sim_data)) {

      pts_idx <- sample(1:nrow(steps), samps, replace = T)
      samp_dist <- steps$dist[pts_idx]
      samp_rel.angle <- steps$rel.angle[pts_idx]
      samp_abs.angle <- samp_rel.angle + sim_data$abs.angle[k - 1]
      samp_x <- sim_data$x[k - 1] + (samp_dist * cos(samp_abs.angle))
      samp_y <- sim_data$y[k - 1] + (samp_dist * sin(samp_abs.angle))
      samp_col.dist <- raster::extract(dist_sim , data.frame(samp_x, samp_y))
      samp_col.dist[is.na(samp_col.dist)] <- Inf
      samp_d.dist <- samp_col.dist - sim_data$col.dist[k - 1]
      samp_bath <- raster::extract(bath_sim, data.frame(samp_x, samp_y))
      samp_bath[is.na(samp_bath)] <- Inf

      if (trip_pos[k] < 0.9) {
        # dd_dist <- density(sample_trip$col.dist[trip_pos == trip_pos[k]], from = 0, to =max(steps$col.dist) * 1.25)
        dd_dist <- density(steps$col.dist, from = 0, to =max(steps$col.dist) * 1.25, adjust = 1)

        dist_fun <- approxfun(x = dd_dist$x, y = dd_dist$y, rule = 2)
        posterior <- dist_fun(samp_col.dist)/sum(dist_fun(samp_col.dist))

        bath_dens <- density(steps$bath,
                             from = min(steps$bath, na.rm = T) * 1.25,
                             to =0, adjust = 1)
        bath_fun <- approxfun(x = bath_dens$x, y = bath_dens$y, rule = 2)

        posterior <- (dist_fun(samp_col.dist)/sum(dist_fun(samp_col.dist)) * bath_fun(samp_bath)/sum(bath_fun(samp_bath)))
      } else {
        p_grid <- seq(from = -20000, to = 20000)
        liklihood <- dnorm(samp_d.dist, mean = -5000, sd = 1000)
        posterior <- liklihood/sum(liklihood)
      }

      water_probs <- raster::extract(water_sim, data.frame(samp_x, samp_y))
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
    out <- rbind(out, sim_data)

  }

  out <- subset(out, out$col.dist > 5000)

  return(out)
  #' @export
}
