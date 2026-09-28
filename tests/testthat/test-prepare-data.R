test_that("prepare-data checks align HMSC four matrices", {
  comm <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 3,
                 dimnames = list(paste0("site", 1:3), c("sp1", "sp2")))
  env <- data.frame(bio1 = c(0.1, 0.2, 0.3), bio12 = c(10, 20, 30),
                    row.names = rownames(comm))
  traits <- data.frame(body = c(1, 2), seed = c(0.5, 0.8),
                       row.names = colnames(comm))
  dat <- hee_prepare_data(comm, env_now = env, traits = traits)
  expect_s3_class(dat, "hee_data")
  expect_true(dat$diagnostics$comm$ok)
  expect_true(dat$diagnostics$env_now$ok)
  expect_true(dat$diagnostics$traits$ok)
})

test_that("prepare-data errors when site and environment rows are not aligned", {
  comm <- matrix(c(1, 0, 0, 1), nrow = 2,
                 dimnames = list(c("site1", "site2"), c("sp1", "sp2")))
  env <- data.frame(bio1 = c(0.1, 0.2), row.names = c("site1", "siteX"))
  expect_error(hee_prepare_data(comm, env_now = env), "Missing site")
})
