test_that("terminal soft likelihood favours its supported state", {
  out <- hee_terminal_soft_likelihood(matrix(c(0.1, 0.9), ncol = 1))
  expect_gt(out$state0_likelihood[1, 1], out$state1_likelihood[1, 1])
  expect_gt(out$state1_likelihood[2, 1], out$state0_likelihood[2, 1])
})

test_that("binary CTMC message is propagated backward correctly", {
  back <- hee_ctmc_backward_message(
    transition_p01 = matrix(0.2, 1, 1),
    transition_p11 = matrix(0.8, 1, 1),
    next_state0_likelihood = matrix(0.9, 1, 1),
    next_state1_likelihood = matrix(0.1, 1, 1)
  )
  # A later observation favouring absence makes an initially absent state
  # more plausible than an initially occupied state after this transition.
  expect_gt(back$state0_likelihood[1, 1], back$state1_likelihood[1, 1])
  q <- hee_ctmc_apply_message(matrix(0.5, 1, 1),
                                back$state0_likelihood,
                                back$state1_likelihood)
  expect_lt(q[1, 1], 0.5)
})

test_that("terminal smoothing validates incompatible matrices", {
  expect_error(
    hee_ctmc_backward_message(matrix(0.2, 1, 1), matrix(0.8, 1, 1),
                              matrix(0.9, 2, 1), matrix(0.1, 1, 1)),
    "identical dimensions"
  )
})

test_that("CTMC interval composition agrees with two successive transitions", {
  out <- hee_ctmc_compose_transition(
    first_p01 = matrix(0.2, 1, 1), first_p11 = matrix(0.8, 1, 1),
    second_p01 = matrix(0.1, 1, 1), second_p11 = matrix(0.9, 1, 1)
  )
  # Starting absent: 0 -> 0 -> 1 or 0 -> 1 -> 1.
  expect_equal(out$transition_p01[1, 1], 0.8 * 0.1 + 0.2 * 0.9)
  # Starting present: 1 -> 0 -> 1 or 1 -> 1 -> 1.
  expect_equal(out$transition_p11[1, 1], 0.2 * 0.1 + 0.8 * 0.9)
})
