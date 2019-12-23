test_that("Contrast against a chi-square hand-rolled version",{
  
  nDraws = 1000
  randomUnifDraws <- runif(nDraws)  
  
  fishersMethodByChiSq <- pchisq(sum(-2*log(randomUnifDraws)), df = 2*nDraws, lower.tail = FALSE)
  
  #S4 class Pvalue allows for sanity checks on return
  expect_s4_class( fishersPvalueSumTest( new("Pvalues", randomUnifDraws)), 
                   "Pvalue")
  
  expect_identical(fishersPvalueSumTest( new("Pvalues", randomUnifDraws)),
                    new("Pvalue", fishersMethodByChiSq),
                   info = "Inspect Pvalues dispatch")

  expect_equivalent(fishersPvalueSumTest( sample(c(randomUnifDraws, NA, NA, NA, NA)) ),
                    new("Pvalue", fishersMethodByChiSq),
                    info = "Inspect numeric dispatch with NA values injected in")
})

