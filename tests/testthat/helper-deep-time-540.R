make_test_timecube_540 <- function(times = c(540, 0), lon = c(-10, 10), lat = c(-5, 5)) {
  e <- new.env(parent = emptyenv())
  e$age_ma <- times
  e$lon <- lon
  e$lat <- lat
  e$variable_names <- c("bio1", "bio12", "elev", "land_mask_dem", "land_area_km2")
  e$product <- list(grid = "test")
  grid <- expand.grid(lon = lon, lat = lat)
  e$get_slice <- function(age, variables = NULL, as_data_table = TRUE) {
    if (is.null(variables)) variables <- e$variable_names
    d <- data.frame(time_ma = age, lon = grid$lon, lat = grid$lat)
    if ("bio1" %in% variables) d$bio1 <- 15 + grid$lat / 5 - age / 540
    if ("bio12" %in% variables) d$bio12 <- 900 + grid$lon * 2 + age / 10
    if ("elev" %in% variables) d$elev <- 100 + grid$lat * grid$lon
    if ("land_mask_dem" %in% variables) d$land_mask_dem <- c(1, 1, 0, 1)[seq_len(nrow(grid))]
    if ("land_area_km2" %in% variables) d$land_area_km2 <- 100
    d
  }
  e$get_metadata <- function(variable = NULL) {
    data.frame(variable = e$variable_names, unit = "test", stringsAsFactors = FALSE)
  }
  e
}
