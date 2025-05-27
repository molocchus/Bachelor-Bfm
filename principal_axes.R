m1 = matrix(0.3,nrow = 6, ncol=6)
m1[1:2,1:2] = m1[1:2,1:2] * 3
m1[3:4,3:4] = m1[3:4,3:4] * 2
m1[5:6,5:6] = m1[5:6,5:6] * 4
m1[1:4,1:4] = m1[1:4,1:4] - 0.5
diag(m1) = 1
m1

R = Holzinger

library(psych)


pa <- function(R,n = 3, i = 1) {
  for( k in 1:i){
    SVD = svd(R)
    A =  diag(sqrt(SVD$d[1:n])) %*% t(SVD$v[,1:n])
    res1 = R - t(A) %*% A
    U2 = diag(res1)
    res = res1
    diag(res1) = 0
    err = sum(abs(res1))
    R = R - diag(U2)
  }
  return(list(t(A) %*% A, err, U2, A))
}

X = pa(R,n = 3, i = 10)




FA = fa(R, 4, rotate = "none")
L = FA$loadings

sl = SL(L)



A = sl %*% t(sl)
A = L %*% t(L)

diag(A) = 0
sum(abs(R - A)) - 14


data = round(sl %*% t(sl),2)
latex_code <- "\begin{tabular}"
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
      latex_code <- paste0(latex_code, data[i, j], " ")
    if (j < ncol(data)) {
      latex_code <- paste0(latex_code, "& ")
    }
  }
  latex_code <- paste0(latex_code, "\\")
}
latex_code <- paste0(latex_code, "\\end{tabular}")
latex_code

