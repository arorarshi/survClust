test_that("combineDist works as expected", {
    expect_equal(round(combineDist(getDist(simdat, simsurvdat))[5,1],5), 0.11824)
})
