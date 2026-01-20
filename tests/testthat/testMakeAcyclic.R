context("Make Acyclic")

test_that("make_acyclic can reverse an edge when forced", {
  
  sample_A <- matrix(c(0, 0, 0,
                       1, 0, 0,
                       1, 1, 0),
                     byrow = TRUE,
                     nrow = 3)
  
  coord <- coordinates(adjMatrix = sample_A)
  
  nEdges <- dim(coord)[2]
  nNodes <- ncol(sample_A)
  
  # Create a cycle 1 -> 2 -> 3 -> 1
  currentES <- c(0,1,0)
  
  # The probability of edge absence is 0
  # This should force only a reverse of an edge in the cycle
  prior <- c(0.5, 0.5, 0.0)
  
  # Verify the correct adj matrix
  #       [,1] [,2] [,3]
  # [1,]    0    1    0
  # [2,]    0    0    1
  # [3,]    1    0    0
  build_adj(coord, currentES, nNodes)
  
  edgeType <- detType(coord = coord,
                      nEdges = nEdges,
                      nCPh = 0,
                      nGV = 0,
                      nNodes = nNodes,
                      pmr = 0)
  
  res <- cycleRmvr(coord, currentES, nNodes, prior, edgeType, pmr, nCPh)
  
  # Should not see any edge absence 
  expect_false(any(res == 2))
})

test_that("make_acyclic removes multiple cycles", {
  
  sample_A <- matrix(c(0, 0, 0, 0,
                       1, 0, 0, 0,
                       1, 1, 0, 0,
                       1, 1, 1, 0),
                     byrow = TRUE,
                     nrow = 4)
  
  coord <- coordinates(adjMatrix = sample_A)
  
  nEdges <- dim(coord)[2]
  nNodes <- ncol(sample_A)
  
  # 2 Cycles
  # 2 -> 4 -> 3 -> 2
  # 1 -> 2 -> 4 -> 1
  currentES <- c(0, 2, 1, 1, 0, 1)
  
  # Approximately equal distribution of the prior for testing
  prior <- c(0.3, 0.3, 0.4)
  
  # Verify the correct adj matrix
  #       [,1] [,2] [,3] [,4]
  # [1,]    0    1    0    0
  # [2,]    0    0    0    1
  # [3,]    0    1    0    0
  # [4,]    1    0    1    0
  build_adj(coord, currentES, nNodes)
  
  edgeType <- detType(coord = coord,
                      nEdges = nEdges,
                      nCPh = 0,
                      nGV = 0,
                      nNodes = nNodes,
                      pmr = 0)
  
  res <- cycleRmvr(coord, currentES, nNodes, prior, edgeType, pmr, nCPh)
  
  # dag_constraint requires an adj matrix
  adj_res <- build_adj(coord, res, nNodes)
  
  # If there are no cycles dag_constraint will return a 0
  expect_equal(dag_constraint(adj_res), 0)
})