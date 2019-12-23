data("siTAZdiffex_HFL1_diffexDT")

siTAZ_betaUnifModel <- fitBetaUniformMixtureDistribution(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue, nStarts = 5)

test_that("Proportions of object should add to 1", {
  
  FP <- falsePositiveFraction(siTAZ_betaUnifModel, 0.05)
  TP <- truePositiveFraction(siTAZ_betaUnifModel, 0.05)
  FN <- falseNegativeFraction(siTAZ_betaUnifModel, 0.05)
  TN <- trueNegativeFraction(siTAZ_betaUnifModel, 0.05)
  
  expect_equal(sum(c(FP,TP,FN,TN)), 1 , tolerence = 1E-6)
  
  expect_lt(noiseFractionUpperBound(siTAZ_betaUnifModel), 1)
})

test_that("Inspect plots produced", {
  
  metaDEGthPlotList <- plot(siTAZ_betaUnifModel)
  
  expect_true(is.ggplot(metaDEGthPlotList$pValHist))
  expect_true(is.ggplot(metaDEGthPlotList$qqPlot))
  expect_true(is.ggplot(metaDEGthPlotList$LLsurface))
})