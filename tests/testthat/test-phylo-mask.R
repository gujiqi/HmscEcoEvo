test_that("phylogenetic time mask blocks species before origin", {
  mask <- hee_phylo_time_mask(times = c(50, 2, 0),
                              species_origin = c(sp1 = 8, sp2 = 60),
                              level = "species")
  expect_equal(mask["sp1", "50Ma"], 0)
  expect_equal(mask["sp1", "2Ma"], 1)
  expect_equal(mask["sp2", "50Ma"], 1)
})

test_that("phylogenetic time mask names unnamed origin vectors from species", {
  mask <- hee_phylo_time_mask(times = c(10, 0),
                              species_origin = c(8, 3),
                              species = c("sp1", "sp2"))
  expect_equal(mask["sp1", "10Ma"], 0)
  expect_equal(mask["sp1", "0Ma"], 1)
  expect_equal(mask["sp2", "10Ma"], 0)
  expect_equal(mask["sp2", "0Ma"], 1)
  expect_error(
    hee_phylo_time_mask(times = c(10, 0), species_origin = c(8, 3)),
    "Unnamed species_origin"
  )
})

test_that("species-level phylogenetic masks require explicit origin ages by default", {
  expect_error(
    hee_phylo_time_mask(times = c(540, 0), species = "sp1"),
    "species_origin"
  )

  mask <- hee_phylo_time_mask(times = c(540, 0), species = "sp1",
                              missing_origin = "oldest")
  expect_equal(mask["sp1", "540Ma"], 1)
  expect_equal(mask["sp1", "0Ma"], 1)
})
