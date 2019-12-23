library(BioNet)

nDraws <- 1E4

fakeShapeParam <- 0.4
fakeUniformFraction <- 0.6

randomBUdraws <- rbetauniform(nDraws, 0.4, 0.6)  

test_that("Inspect correctness of fit against known values", {
  
  refittedBetaUniformModel <- fitBetaUniformMixtureDistribution(randomBUdraws, nStarts = 5)
  
  expect_s4_class(refittedBetaUniformModel, "BetaUniformModel")
  
  expect_equivalent( as.numeric(refittedBetaUniformModel@a), fakeShapeParam, tolerance = 0.05)
  expect_equivalent( as.numeric(refittedBetaUniformModel@lambda), fakeUniformFraction, tolerance = 0.05)
})


test_that("Contrast against another implementation from BioNet", {
  
  refittedBetaUniformModel <- fitBetaUniformMixtureDistribution(randomBUdraws)
  
  suppressWarnings(refittedBUM <- fitBumModel(randomBUdraws, plot = FALSE))
  
  expect_equivalent( refittedBetaUniformModel, convertFromBUM(refittedBUM), tolerance = 0.05)
})

