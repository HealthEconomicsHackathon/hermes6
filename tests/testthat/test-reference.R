test_that("output is the same as expected", {
  reference <- readRDS(file.path("testdata", "reference.RDS"))
  markov <- markov()
  expect_equal(markov, reference)
})
