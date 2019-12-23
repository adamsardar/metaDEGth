nDraws <- 1000

fakeShapeParam <- 0.4
fakeUniformFraction <- 0.6

randomBUdraws <- rbetauniform(nDraws, 0.4, 0.6)  

fakeBUM_S3obj <- list(a = mkScalar(fakeShapeParam), 
                      lambda = fakeUniformFraction, 
                      pvalues = randomBUdraws,
                      negLL = -sum(dbetauniform(randomBUdraws, 0.4, 0.6, log = TRUE)))
class(fakeBUM_S3obj) <- "bum"


test_that("Contrast Fisher and BetaUniformSum tests",{

  expect_error(rbetauniform(nDraws, -0.4, 0.6))
  expect_error(rbetauniform(nDraws, 0.4, 1.2))
  
  expect_s4_class(convertFromBUM(fakeBUM_S3obj), "BetaUniformModel")
  
  fakeBetaUniformObj <- convertFromBUM(fakeBUM_S3obj)
  
  #S4 class Pvalue allows for sanity checks on return
  expect_s4_class( betaUniformPvalueSumTest( new("Pvalues", sample(c(sample( seq(0.01,0.1,0.01), 10), NA, NA, NA, NA))), 
                                             fakeBetaUniformObj), 
                   "Pvalue")
  
  expect_warning(betaUniformPvalueSumTest( new("Pvalues", c(1E-320, sample(c(sample( seq(0.01,0.1,0.01), 10), NA, NA, NA, NA)))), fakeBetaUniformObj), regexp = "machine precision")
  
  for(i in 1:10){

    fakePvalSet <- sample( seq(0.01,0.1,0.01), 10)
    
    expect_lt(fishersPvalueSumTest( new("Pvalues", fakeSignificantSet)),
              betaUniformPvalueSumTest( new("Pvalues", fakeSignificantSet), fakeBetaUniformObj),
              label = "The beta-uniform parameterisation P-value should always be larger than the Fisher")    
  }
  
})



test_that("Inspect transition rate matrix", {
  
  nVal <- 10

  # The phase-type distribution is defined by its transition matrix
  # https://en.wikipedia.org/wiki/Transition_rate_matrix  
  transitionMatrix <- constructHyperexponentialSumTransitionMatrix(nVal, fakeUniformFraction, fakeShapeParam)
  
  # A 2m x 2m matrix
  expect_s4_class(transitionMatrix,
                  "Matrix")
  
  expect_equal(nrow(transitionMatrix), ncol(transitionMatrix))
  expect_equal(nrow(transitionMatrix), 2*nVal)

  #Row sums of zero except last two
  expect_equal(rowSums(transitionMatrix[-(2*nVal-1):(-2*nVal),]),
               rep(0, times = 2*(nVal-1)),
               label = "Row sums of a laplacian are zero")
  
  expect_true(all(diag(transitionMatrix) <= 0), label = "Diagonal entries are all negative semi-definite")

  transitionMatrixCopy <- transitionMatrix
  diag(transitionMatrixCopy) <- 0
  expect_true(all(transitionMatrixCopy >= 0 ), label = "All off-diagonal entries are positive semi-definite")
})