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


test_that("Inspect transition rate matrix - should have known properties", {
  
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


# An older but much slower version of transition matrix generator. Used to test against
constructHyperexponentialSumTransitionMatrixReference <- function(nValues, uniformProportion, fittedBetaShape){
  
  recurrentStateMatrixTile <- Matrix(c(-1,0,0,-fittedBetaShape), ncol = 2, nrow = 2)
  transitionStateMatrixTile <- Matrix(c(uniformProportion, uniformProportion*fittedBetaShape,
                                        (1-uniformProportion),(1-uniformProportion)*fittedBetaShape), nrow = 2, ncol = 2)
  
  hyperexponentialSumTransitionMatrix <- kronecker(Matrix::diag(nrow = nValues, ncol = nValues), recurrentStateMatrixTile)
  
  # If the number of values is 1, then we have a hyperexponeitial distribution and this transition matrix is suitable
  
  # For nValues = 2, we can just insert the tile describing transition from state 1 to 2 only
  if(nValues == 2){  hyperexponentialSumTransitionMatrix[1:2,3:4] <- transitionStateMatrixTile }
  
  # for the general case, 
  if(nValues > 2){
    
    hyperexponentialSumTransitionMatrix <- hyperexponentialSumTransitionMatrix +
      # Create an off diagonal matrix and add the transition tiles in
      kronecker({M <- Matrix(0, nrow = nValues, ncol = nValues); diag(M[,-1]) <- 1; M}, transitionStateMatrixTile) }
  
  return(hyperexponentialSumTransitionMatrix)
}

test_that("Contrast transition rate matrix against an alternative implemetnation", {
  
  nVal <- 1
  
  # 50 runs because this transition matrix is essential to be correct!
  for(i in 1:50){
  
    fakeShape <- runif(1, 0.1, 0.8)
    fakeMixture <- runif(1, 0.1, 0.8)
    
    mod <- constructHyperexponentialSumTransitionMatrix(nVal, fakeMixture, fakeShape)
    ref <- constructHyperexponentialSumTransitionMatrixReference(nVal, fakeMixture, fakeShape)
    
    expect_equivalent(mod, ref)
    
    #No zero counts
    # nVal =1 and nVal=2 are odd edge cases that *must* be checked
    if(i == 1){nVal <- 2}else{nVal <- rpois(1,10) + 1}
  }
  
})


test_that("Compare accuracy against known phase-type distribution P-value (actuar package)",{
  
  nVal <- 5
  
  fakeBetaUniformObj <- convertFromBUM(fakeBUM_S3obj)
  
  fakePvalSet <- sample( seq(0.01,0.1,0.01), nVal)
  
  # The phase-type distribution is defined by its transition matrix
  # https://en.wikipedia.org/wiki/Transition_rate_matrix  
  transitionMatrix <- constructHyperexponentialSumTransitionMatrix(nVal, fakeUniformFraction, fakeShapeParam)
  
  expect_equal(
    as.numeric(betaUniformPvalueSumTest( new("Pvalues", fakePvalSet), fakeBetaUniformObj)),
  
    actuar::pphtype(sum(-log(fakePvalSet)),
          prob = c(c(fakeUniformFraction,1-fakeUniformFraction),rep(0, times = 2*(nVal-1))),
          rates = as.matrix(transitionMatrix),
          lower.tail = FALSE),
    
    tolerance = 1E-5)
  
})

