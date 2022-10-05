test_that("getDist works as expected", {
  expect_equal(round(getDist(simdat, simsurvdat)[[1]][5,1],5), 0.11824)
})
