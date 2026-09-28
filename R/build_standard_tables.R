#' Build a species-region history table from standardized events
#'
#' @param raw_events Standardized event table from `extract_biogeobears_events()`.
#' @param tip_ranges Optional present-day species-region ranges.
#' @param species Character vector of species to retain.
#' @param regions Character vector of valid regions.
#' @return A species-region history data.frame.
#' @export
build_species_region_history <- function(raw_events, tip_ranges = NULL,
                                         species = NULL, regions = NULL) {
  if (is.null(raw_events) || nrow(raw_events) == 0) return(NULL)
  d <- raw_events
  d$event_prob <- .std_prob(d$event_prob, 1)

  disp <- d[d$event_type %in% c("dispersal", "founder") & !is.na(d$to_region), , drop = FALSE]
  if (!is.null(species)) disp <- disp[is.na(disp$species) | disp$species %in% species, , drop = FALSE]
  if (!is.null(regions)) disp <- disp[disp$to_region %in% regions, , drop = FALSE]

  # If species is missing from BSM rows, this table cannot be species resolved.
  disp <- disp[!is.na(disp$species) & disp$species != "", , drop = FALSE]

  rows <- list()
  if (nrow(disp) > 0) {
    key <- paste(disp$species, disp$to_region, sep = "\r")
    for (k in unique(key)) {
      tab <- disp[key == k, , drop = FALSE]
      sp <- tab$species[1]
      reg <- tab$to_region[1]
      w <- tab$event_prob
      rows[[length(rows) + 1]] <- data.frame(
        species = sp,
        region = reg,
        entry_time = .weighted_mean_safe(tab$time, w),
        entry_time_dispersion = .weighted_sd_safe(tab$time, w),
        source_region = .weighted_mode(tab$from_region, w),
        entry_prob = min(1, sum(w, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
  }

  out <- .safe_bind_rows(rows)
  if (is.null(out)) out <- data.frame(species = character(), region = character())

  # Add current ranges as retained information when available.
  tr <- .standardize_tip_ranges(tip_ranges)
  if (!is.null(tr) && nrow(tr) > 0) {
    tr <- tr[tr$presence > 0, , drop = FALSE]
    tr$retained <- 1
    tr$entry_time <- NA_real_
    tr$entry_time_dispersion <- NA_real_
    tr$source_region <- NA_character_
    tr$entry_prob <- tr$presence
    tr <- tr[c("species", "region", "entry_time", "entry_time_dispersion", "source_region", "entry_prob", "retained")]
    out$retained <- NA_real_
    out <- .safe_bind_rows(list(out, tr))
    out <- .collapse_species_region_history(out)
  } else {
    out$retained <- NA_real_
  }

  # Add loss time from loss events.
  loss <- d[d$event_type == "loss" & !is.na(d$region) & !is.na(d$species), , drop = FALSE]
  if (nrow(loss) > 0 && nrow(out) > 0) {
    loss_key <- paste(loss$species, loss$region, sep = "\r")
    loss_rows <- lapply(unique(loss_key), function(k) {
      tab <- loss[loss_key == k, , drop = FALSE]
      data.frame(
        species = tab$species[1],
        region = tab$region[1],
        loss_time = .weighted_mean_safe(tab$time, tab$event_prob),
        stringsAsFactors = FALSE
      )
    })
    loss_df <- .safe_bind_rows(loss_rows)
    out <- merge(out, loss_df, by = c("species", "region"), all.x = TRUE, sort = FALSE)
  } else {
    out$loss_time <- NA_real_
  }

  out <- out[order(out$species, out$region), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.collapse_species_region_history <- function(x) {
  if (is.null(x) || nrow(x) == 0) return(x)
  key <- paste(x$species, x$region, sep = "\r")
  rows <- lapply(unique(key), function(k) {
    tab <- x[key == k, , drop = FALSE]
    w <- .std_prob(tab$entry_prob, 1)
    data.frame(
      species = tab$species[1],
      region = tab$region[1],
      entry_time = .weighted_mean_safe(tab$entry_time, w),
      entry_time_dispersion = .weighted_sd_safe(tab$entry_time, w),
      source_region = .weighted_mode(tab$source_region, w),
      entry_prob = min(1, sum(w, na.rm = TRUE)),
      retained = if (all(is.na(tab$retained))) NA_real_ else max(tab$retained, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  .safe_bind_rows(rows)
}

.standardize_tip_ranges <- function(tip_ranges) {
  if (is.null(tip_ranges)) return(NULL)
  if (is.matrix(tip_ranges)) {
    if (is.null(rownames(tip_ranges)) || is.null(colnames(tip_ranges))) {
      stop("tip_ranges matrix must have species rownames and region colnames.", call. = FALSE)
    }
    idx <- which(tip_ranges > 0, arr.ind = TRUE)
    if (nrow(idx) == 0) return(data.frame(species = character(), region = character(), presence = numeric()))
    return(data.frame(
      species = rownames(tip_ranges)[idx[, 1]],
      region = colnames(tip_ranges)[idx[, 2]],
      presence = as.numeric(tip_ranges[idx]),
      stringsAsFactors = FALSE
    ))
  }
  x <- as.data.frame(tip_ranges, stringsAsFactors = FALSE)
  if (all(c("species", "region") %in% names(x))) {
    if (!"presence" %in% names(x)) {
      prob_col <- intersect(c("prob", "range_prob", "value", "present"), names(x))
      x$presence <- if (length(prob_col) > 0) as.numeric(x[[prob_col[1]]]) else 1
    }
    return(x[c("species", "region", "presence")])
  }
  stop("tip_ranges must be a species x region matrix or a data.frame with species and region columns.", call. = FALSE)
}

#' Build a speciation event table from standardized events
#'
#' @param raw_events Standardized event table.
#' @param phy Optional tree. Not used in v0.1.
#' @return A speciation event data.frame or NULL.
#' @export
build_speciation_events <- function(raw_events, phy = NULL) {
  if (is.null(raw_events) || nrow(raw_events) == 0) return(NULL)
  d <- raw_events[raw_events$event_type %in% c("speciation", "founder"), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  desc <- d$descendant_species
  desc[is.na(desc) | desc == ""] <- d$species[is.na(desc) | desc == ""]
  out <- data.frame(
    event_id = d$event_id,
    region = ifelse(!is.na(d$region), d$region, d$to_region),
    time = d$time,
    descendant_species = desc,
    event_prob = .std_prob(d$event_prob, 1),
    event_type = d$event_type,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$descendant_species) & out$descendant_species != "" & !is.na(out$region), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build a history event table from standardized events
#'
#' @param raw_events Standardized event table.
#' @return A dispersal/loss event data.frame or NULL.
#' @export
build_history_events <- function(raw_events) {
  if (is.null(raw_events) || nrow(raw_events) == 0) return(NULL)
  d <- raw_events[raw_events$event_type %in% c("dispersal", "loss", "founder"), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  d$event_type[d$event_type == "founder"] <- "dispersal"
  out <- data.frame(
    event_id = d$event_id,
    event_type = d$event_type,
    species = d$species,
    from_region = d$from_region,
    to_region = d$to_region,
    region = d$region,
    time = d$time,
    event_prob = .std_prob(d$event_prob, 1),
    stringsAsFactors = FALSE
  )
  # Keep usable rows only.
  ok_disp <- out$event_type == "dispersal" & !is.na(out$from_region) & !is.na(out$to_region)
  ok_loss <- out$event_type == "loss" & !is.na(out$region)
  out <- out[ok_disp | ok_loss, , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Build or standardize species-level trait history
#'
#' @param traits Optional species trait table.
#' @param ancestral_traits Optional precomputed ancestral trait summaries.
#' @param trait_history Optional already prepared species-level trait history table.
#' @param species Species vector for alignment.
#' @return A trait history data.frame or NULL.
#' @export
build_trait_history <- function(traits = NULL, ancestral_traits = NULL,
                                trait_history = NULL, species = NULL) {
  if (!is.null(trait_history)) {
    x <- as.data.frame(trait_history, stringsAsFactors = FALSE)
    if ("species" %in% names(x)) rownames(x) <- as.character(x$species)
    if (!"species" %in% names(x)) x$species <- rownames(x)
    if (!is.null(species)) x <- x[x$species %in% species, , drop = FALSE]
    rownames(x) <- x$species
    return(x)
  }

  # v0.1 does not infer ancestral traits. If ancestral_traits is already a
  # species-level summary, standardize it.
  if (!is.null(ancestral_traits)) {
    x <- as.data.frame(ancestral_traits, stringsAsFactors = FALSE)
    if ("species" %in% names(x)) rownames(x) <- as.character(x$species)
    if (is.null(rownames(x))) stop("ancestral_traits must have rownames or a species column.", call. = FALSE)
    if (!"species" %in% names(x)) x$species <- rownames(x)
    if (!is.null(species)) x <- x[x$species %in% species, , drop = FALSE]
    rownames(x) <- x$species
    return(x)
  }

  if (!is.null(traits)) {
    x <- .align_traits(traits, species, "traits")
    if (is.null(x) || nrow(x) == 0) return(NULL)
    x$species <- rownames(x)
    rownames(x) <- x$species
    return(x)
  }
  NULL
}
