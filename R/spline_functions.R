

constructBS <- function(x, knots, degree, ref_value){
  L <- length(knots)
  B <- splines::bs(x, knots = knots[-c(1,L)], Boundary.knots = knots[c(1,L)], intercept = T, degree = degree)

  ref_pos <- which(knots == ref_value)
  to_rm <- ref_pos
  if(degree >= 2) to_rm <- c(to_rm, ref_pos + 1)
  if(degree >= 3) to_rm <- c(ref_pos - 1, to_rm)
  B <- B[, -to_rm]
  if(degree > 1) B <- cbind(poly(x - ref_value, degree=degree-1, raw=T), B)
  B
}
