## Pooled variance function
# varVec <- c(1, 2, 3)
# nVec <- c(5, 6, 7)

pooledVar <- function(varVec, nVec) {
  numVec <- rep(NA, length(varVec))
  for (i in seq_along(varVec)) {
    numVec[i] <- (nVec[i] - 1) * varVec[i]
  }
  denVec <- sapply(nVec, function(x) {x - 1})
  sum(numVec) / sum(denVec)
}
