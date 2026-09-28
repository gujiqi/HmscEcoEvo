# Derived movement-network results; plate carriage is never biological flow.

#' Aggregate joint posterior movement mass into patch corridors
#'
#' This function requires posterior *edge* mass, not marginal location maps.
#' A caller must obtain the source-to-target joint distribution from a
#' transition model and its forward/backward messages, or from weighted path
#' draws. It cannot be reconstructed by multiplying two location marginals.
#' Carriage and same-cell retention are excluded from active-dispersal
#' summaries. Patch IDs are time-specific, so a source and target patch with
#' different IDs are not automatically different historical patch networks.
#'
#' @param edge_flow Table with `from_cell_id`, `to_cell_id`, `time_from_ma`,
#'   `time_to_ma`, `posterior_edge_mass`, and `edge_kind`.
#' @param from_membership,to_membership Patch-membership tables for the two
#'   time endpoints. Each cell may belong to at most one patch per endpoint.
#' @return A `corridors` table of directed flow among time-specific patch
#'   endpoints, plus mass on active edges with no supported patch at one or
#'   both endpoints. No demographic source-sink or persistent-network label
#'   is inferred.
#' @export
hee_result_corridor_flow <- function(edge_flow, from_membership,
                                     to_membership) {
  f <- as.data.frame(edge_flow); a <- as.data.frame(from_membership)
  b <- as.data.frame(to_membership)
  need <- c("from_cell_id", "to_cell_id", "time_from_ma", "time_to_ma",
            "posterior_edge_mass", "edge_kind")
  if (!all(need %in% names(f)) ||
      !all(c("cell_id", "patch_id") %in% names(a)) ||
      !all(c("cell_id", "patch_id") %in% names(b)) ||
      anyDuplicated(a$cell_id) || anyDuplicated(b$cell_id) ||
      any(!is.finite(f$posterior_edge_mass) | f$posterior_edge_mass < 0) ||
      any(!is.finite(f$time_from_ma) | !is.finite(f$time_to_ma) |
          f$time_from_ma <= f$time_to_ma)) {
    stop("Joint edge-flow input is missing or invalid; marginal maps are insufficient.",
         call. = FALSE)
  }
  empty <- data.frame(from_patch_id = character(), to_patch_id = character(),
    time_from_ma = numeric(), time_to_ma = numeric(),
    posterior_movement_mass = numeric(), stringsAsFactors = FALSE)
  if (!nrow(f)) return(list(corridors = empty, unassigned_edge_mass = 0))
  f <- f[f$edge_kind == "active_dispersal" &
           f$from_cell_id != f$to_cell_id & f$posterior_edge_mass > 0, , drop = FALSE]
  if (!nrow(f)) return(list(corridors = empty, unassigned_edge_mass = 0))
  total_active_mass <- sum(f$posterior_edge_mass)
  x <- merge(f, setNames(a[, c("cell_id", "patch_id")],
                          c("from_cell_id", "from_patch_id")), by = "from_cell_id")
  x <- merge(x, setNames(b[, c("cell_id", "patch_id")],
                          c("to_cell_id", "to_patch_id")), by = "to_cell_id")
  if (!nrow(x)) return(list(corridors = empty,
                            unassigned_edge_mass = total_active_mass))
  out <- stats::aggregate(x$posterior_edge_mass,
    x[, c("from_patch_id", "to_patch_id", "time_from_ma", "time_to_ma")], sum)
  names(out)[[5L]] <- "posterior_movement_mass"
  list(corridors = out[, names(empty), drop = FALSE],
       unassigned_edge_mass = total_active_mass - sum(x$posterior_edge_mass))
}

#' Classify patch dispersal exchange without demographic source-sink claims
#'
#' Incoming and outgoing numbers are expected active-movement mass from joint
#' edge-flow inference. A patch can be an exporter, recipient, exchange hub,
#' or weakly connected patch relative to a declared absolute threshold.
#' These labels are *not* local reproduction-based source/sink categories.
#'
#' @param corridors Output `$corridors` from [hee_result_corridor_flow()] or
#'   a table with the same columns.
#' @param minimum_flow Minimum expected movement mass for a high-flow role.
#' @return Patch-level incoming/outgoing flow and role table.
#' @export
hee_result_donor_recipient <- function(corridors, minimum_flow = 0.05) {
  x <- as.data.frame(corridors)
  need <- c("from_patch_id", "to_patch_id", "posterior_movement_mass")
  if (!all(need %in% names(x)) || length(minimum_flow) != 1L ||
      !is.finite(minimum_flow) || minimum_flow < 0 ||
      any(!is.finite(x$posterior_movement_mass) |
          x$posterior_movement_mass < 0)) {
    stop("Invalid posterior corridor flow or minimum_flow.", call. = FALSE)
  }
  empty <- data.frame(patch_id = character(), outgoing_mass = numeric(),
    incoming_mass = numeric(), net_export_mass = numeric(), role = character())
  if (!nrow(x)) return(empty)
  ids <- sort(unique(c(as.character(x$from_patch_id),
                       as.character(x$to_patch_id))))
  outgoing <- stats::setNames(rep(0, length(ids)), ids)
  incoming <- outgoing
  out_sum <- rowsum(x$posterior_movement_mass, x$from_patch_id)
  in_sum <- rowsum(x$posterior_movement_mass, x$to_patch_id)
  outgoing[rownames(out_sum)] <- out_sum[, 1L]
  incoming[rownames(in_sum)] <- in_sum[, 1L]
  high_out <- outgoing >= minimum_flow & outgoing > 0
  high_in <- incoming >= minimum_flow & incoming > 0
  role <- ifelse(high_out & high_in, "exchange_hub",
          ifelse(high_out, "dispersal_exporter",
          ifelse(high_in, "dispersal_recipient", "weak_exchange")))
  data.frame(patch_id = ids, outgoing_mass = unname(outgoing),
    incoming_mass = unname(incoming),
    net_export_mass = unname(outgoing - incoming), role = unname(role),
    stringsAsFactors = FALSE)
}
