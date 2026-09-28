\name{hee_dispersal_rare_land_jump_transition}
\alias{hee_dispersal_rare_land_jump_transition}
\title{Add a rare land-to-land jump tail to a location transition}
\usage{
hee_dispersal_rare_land_jump_transition(
  transition, land_cells, hub_cells, jump_rate_per_myr, delta_t_myr
)
}
\arguments{
\item{transition}{Named target-by-source column-stochastic matrix or location-transition object.}
\item{land_cells}{Logical terrestrial availability at the movement stage.}
\item{hub_cells}{Declared terrestrial landing-state IDs.}
\item{jump_rate_per_myr}{Non-negative rare-jump rate per Myr.}
\item{delta_t_myr}{Positive movement-stage duration in Myr.}
}
\value{A sparse source-normalised transition matrix on the original location states.}
\description{
Mixes a low-rate direct land-to-land jump into an existing active movement
transition. It does not create a marine ancestral location. The rate is a
documented sensitivity parameter unless calibrated with independent data.
}
