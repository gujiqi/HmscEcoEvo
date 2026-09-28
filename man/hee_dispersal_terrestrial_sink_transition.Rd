\name{hee_dispersal_terrestrial_sink_transition}
\alias{hee_dispersal_terrestrial_sink_transition}
\title{Enforce terrestrial states after diffusion and plate carriage}
\usage{
hee_dispersal_terrestrial_sink_transition(
  transition, source_land, target_land,
  sink_state = "__marine_or_lost_land__"
)
}
\arguments{
\item{transition}{Named target-by-source column-stochastic matrix or location-transition object.}
\item{source_land}{Logical availability at the older endpoint, named by state or in matrix order.}
\item{target_land}{Logical availability at the younger endpoint, named by state or in matrix order.}
\item{sink_state}{Name for the absorbing unavailable-geography state.}
}
\value{A list containing a sparse source-normalised \code{transition}, the \code{sink_state}, and diagnostics.}
\description{
Marine-target mass enters an explicit absorbing sink rather than being
renormalised onto land. Direct land-to-land crossings remain possible if
they were already specified by the active dispersal transition. The sink
describes unavailable geography, not inferred biological extinction.
}
