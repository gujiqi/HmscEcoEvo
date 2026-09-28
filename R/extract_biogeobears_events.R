#' Extract standardized event rows from BioGeoBEARS/BSM-like input
#'
#' This is a conservative extractor. It accepts a data.frame or a list of
#' data.frames with event-like columns and maps them to HmscEcoEvo's internal
#' event format. It does not run BioGeoBEARS and does not infer events.
#'
#' @param biogeobears Optional BioGeoBEARS result object. Event-like data.frames
#'   are searched recursively.
#' @param bsm_events Optional BSM events as data.frame or list of data.frames.
#' @return A standardized event data.frame.
#' @export
extract_biogeobears_events <- function(biogeobears = NULL, bsm_events = NULL) {
  candidates <- list()

  if (!is.null(bsm_events)) {
    if (is.data.frame(bsm_events)) {
      candidates <- c(candidates, list(bsm_events))
    } else if (is.list(bsm_events)) {
      dfs <- .find_data_frames(bsm_events)
      candidates <- c(candidates, unname(dfs))
    } else {
      stop("bsm_events must be a data.frame or a list containing data.frames.", call. = FALSE)
    }
  }

  if (!is.null(biogeobears)) {
    dfs <- .find_data_frames(biogeobears)
    if (length(dfs) > 0) {
      scores <- vapply(dfs, .event_like_score, numeric(1))
      dfs <- dfs[scores >= 3]
      candidates <- c(candidates, unname(dfs))
    }
  }

  if (length(candidates) == 0) {
    stop("No event-like data.frames found. Provide bsm_events as a data.frame, or provide manual history tables.", call. = FALSE)
  }

  std <- lapply(seq_along(candidates), function(i) .standardize_event_like_df(candidates[[i]], map_id = i))
  out <- .safe_bind_rows(std)
  if (is.null(out) || nrow(out) == 0) {
    stop("Event-like tables were found, but none could be standardized. Please provide manual standard tables.", call. = FALSE)
  }
  rownames(out) <- NULL
  out
}

.standardize_event_like_df <- function(d, map_id = NA_integer_) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)

  d <- .rename_first_existing(d, "event_type", c("type", "event", "eventName", "event_name", "etype"))
  d <- .rename_first_existing(d, "time", c("age", "event_time", "abs_time", "t", "ma"))
  d <- .rename_first_existing(d, "from_region", c("from", "source", "source_region", "fromArea", "from_area", "src"))
  d <- .rename_first_existing(d, "to_region", c("to", "target", "dest", "destination", "toArea", "to_area"))
  d <- .rename_first_existing(d, "region", c("area", "lost_region", "range", "state"))
  d <- .rename_first_existing(d, "species", c("tip", "taxon", "taxa", "terminal"))
  d <- .rename_first_existing(d, "descendant_species", c("descendant", "descendants", "desc", "tips"))
  d <- .rename_first_existing(d, "event_prob", c("prob", "probability", "weight", "posterior", "freq", "frequency"))
  d <- .rename_first_existing(d, "node", c("node_id", "nodenum", "node_number"))
  d <- .rename_first_existing(d, "branch", c("branch_id", "edge", "edge_id"))
  d <- .rename_first_existing(d, "event_id", c("id", "eventid"))

  if (!"event_type" %in% names(d)) {
    # If no event type exists, this is not usable.
    return(NULL)
  }

  d$event_type <- .map_event_type(d$event_type)
  d <- .add_if_missing(d, "event_id", paste0("event_", map_id, "_", seq_len(nrow(d))))
  d <- .add_if_missing(d, "species", NA_character_)
  d <- .add_if_missing(d, "descendant_species", NA_character_)
  d <- .add_if_missing(d, "node", NA)
  d <- .add_if_missing(d, "branch", NA)
  d <- .add_if_missing(d, "from_region", NA_character_)
  d <- .add_if_missing(d, "to_region", NA_character_)
  d <- .add_if_missing(d, "region", NA_character_)
  d <- .add_if_missing(d, "time", NA_real_)
  d <- .add_if_missing(d, "event_prob", 1)

  d$event_prob <- .std_prob(d$event_prob, 1)
  d$time <- suppressWarnings(as.numeric(d$time))
  d$map_id <- map_id

  keep <- c("event_id", "event_type", "species", "descendant_species", "node", "branch",
            "from_region", "to_region", "region", "time", "event_prob", "map_id")
  d[keep]
}
