x <- seq(0,100,.1)
ref_value <- 50
knots <- 20
a <- 2
c <- .001
poly_degree <- 1

c_new <- c + ref_value
lambda <- (a - 1) / a
x_new <- ((x + c)^lambda - c_new^lambda)/(lambda * c_new^(lambda - 1))


# design matrix
x_pos <- pmax(x - ref_value, 0)
x_neg <- pmax(ref_value - x, 0)

B_pos <- Compute_Design(x = x_pos, k = knots, region = range(x_pos))
B_neg <- Compute_Design(x = x_neg, k = knots, region = range(x_neg))
matplot(sort(x), cbind(B_neg,B_pos)[order(x),], type = "l")

# interpolation
X_int1 <- Matrix::sparse.model.matrix(~ 0 + x_new)
X_int2 <- stats::poly(x_new, degree = poly_degree, raw = TRUE)


# precision matrix
P_pos <- Compute_Prec(a = a, c = c_new, k = knots, region = range(x_pos))
P_neg <- Compute_Prec_rev(a = a, c = c_new, k = knots, region = range(x_neg))
P <- Matrix::bdiag(P_pos, P_neg)

Matrix::image(P)
Matrix::diag(P); plot(Matrix::diag(P))
