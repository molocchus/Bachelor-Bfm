m1 = matrix(0.1,nrow = 6, ncol=6)
m1[1:2,1:2] = m1[1:2,1:2] * 3
m1[3:4,3:4] = m1[1:2,1:2] * 2
m1[5:6,5:6] = m1[1:2,1:2] * 4

m1 = m1 + matrix(runif(36), nrow = 6, ncol = 6)/50

diag(m1) = c(rep(1,6))
m1

library(psych)


#######Pierwsze rozwiązanie SL
SL = schmid(m1, 3)
yee = round(SL$sl, 2)[1:6,1:4]
yee


data = yee
latex_code <- "\begin{tabular}"
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if (abs(data[i, j]) > 0.3) {
      latex_code <- paste0(latex_code, " $\\mathbf{", data[i, j], "}$ ")
    } else {
      latex_code <- paste0(latex_code, data[i, j], " ")
    }
    if (j < ncol(data)) {
      latex_code <- paste0(latex_code, "& ")
    }
  }
  latex_code <- paste0(latex_code, "\\")
}
latex_code <- paste0(latex_code, "\\end{tabular}")
latex_code



#######Normalizacja danych

normalized = yee[,2:4]

for( i in 1:nrow(normalized)){
  for(j in 1:ncol(normalized)){
    normalized[i,j] = normalized[i,j] / sqrt(rowSums(yee[,2:4]^2))[i]
  }
}


data = round(normalized^2,2)

latex_code <- "\begin{tabular}"
for (i in 1:nrow(data)) {
  for (j in 1:ncol(data)) {
    if (abs(data[i, j]) > 0.3) {
      latex_code <- paste0(latex_code, " $\\mathbf{", data[i, j], "}$ ")
    } else {
      latex_code <- paste0(latex_code, data[i, j], " ")
    }
    if (j < ncol(data)) {
      latex_code <- paste0(latex_code, "& ")
    }
  }
  latex_code <- paste0(latex_code, "\\")
}
latex_code <- paste0(latex_code, "\\end{tabular}")
latex_code



data = round(normalized^2,2)
column_indices <- apply(data, 2, order, decreasing=TRUE)

latex_code <- "\\begin{tabular}\n"
for (i in 1:nrow(data)) {
  print(i)
  for (j in 1:ncol(data)) {
    sorted_row_idx <- column_indices[i, j]
    sorted_value <- data[sorted_row_idx, j]
    if (abs(sorted_value) > 0.3) {
      latex_code <- paste0(latex_code, sorted_row_idx, " & $ \\mathbf{", sorted_value, "}$  ")
    } else {
      latex_code <- paste0(latex_code, sorted_row_idx, " & ", sorted_value)
    }
    if (j < ncol(data)) {
      latex_code <- paste0(latex_code, "& ")
    }
  }
  latex_code <- paste0(latex_code, "\\\\\n")
}
latex_code <- paste0(latex_code, "\\end{tabular}")

print(latex_code)


sorted_df <- data.frame(apply(round(normalized^2, 2), 2, function(x) x[order(-x)]))

data = sorted_df
data[1:5,] = data[1:5,] - data[2:6,]
data


latex_code <- "\\begin{tabular}\n"
for (i in 1:nrow(data)) {
  print(i)
  for (j in 1:ncol(data)) {
    sorted_row_idx <- column_indices[i, j]
    sorted_value <- data[i, j]
    if (abs(sorted_value) > 0.3) {
      latex_code <- paste0(latex_code, sorted_row_idx, " & $ \\mathbf{", sorted_value, "}$  ")
    } else {
      latex_code <- paste0(latex_code, sorted_row_idx, " & ", sorted_value)
    }
    if (j < ncol(data)) {
      latex_code <- paste0(latex_code, "& ")
    }
  }
  latex_code <- paste0(latex_code, "\\\\\n")
}
latex_code <- paste0(latex_code, "\\end{tabular}")

print(latex_code)








