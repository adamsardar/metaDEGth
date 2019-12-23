library(BioNet)
data("siTAZdiffex_HFL1_diffexDT")

siTAZ_betaUnifModel <- fitBetaUniformMixtureDistribution(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue, nStarts = 5)

test_that("Ensure equivilence with BioNet implementation", {
  
  metaDEGth_betaUniformScores <- betaUniformScore(siTAZ_betaUnifModel, FDR = 0.05)
  
  BioNet_betaUniformScores <- scoreFunction(list(a = siTAZ_betaUnifModel@a,
            lambda = siTAZ_betaUnifModel@lambda,
            pvalues = siTAZ_betaUnifModel@pvalues@.Data), fdr = 0.05)
  
  expect_equal(metaDEGth_betaUniformScores, BioNet_betaUniformScores, label = 'Log-ratios of likelihoods should be equivilent')
})

test_that("Ensure that smallest P-value is given the best score (most likely from signal)",{
  
  bestScore <- betaUniformScore(sort(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue)[1], siTAZ_betaUnifModel, FDR = 0.05)
  
  lowerScores <- betaUniformScore( sort(siTAZdiffex_HFL1_diffexDT$siTAZ1_pValue)[-1], siTAZ_betaUnifModel, FDR = 0.05)
  
  expect_true(all(lowerScores < rep(bestScore, times = length(lowerScores))))
})

