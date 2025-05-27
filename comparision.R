omega_h = function(rot, R){
  sum(rot[,1])^2/sum(R)
}

ECV = function(rot){
  sq = colSums(rot^2)
  return(sq[1]/sum(sq))
}

omega_t = function(rot, R){
  sum(colSums(rot)^2)/sum(R)
}

SL = function(L){
  lp_rotated = perform_rotation(L, method = "geominQ")
  lp1 <- fa(lp_rotated$Phi, nfactors=1, fm = "pa")
  lpSL1          <- lp_rotated$loadings %*% lp1$loadings
  psl1           <- matrix (0, dim(lp1$loadings), dim(lp1$loadings))
  diag(psl1)     <- sqrt(1- lp1$loadings^2)
  
  if (any(is.na(diag(psl1)))== TRUE) {
    stop(return(list(SL_loadings = NA,loadings = NA)))
  }
  
  lpsl2          <- lp_rotated$loadings %*% psl1
  SL_loadings    <- cbind (lpSL1,lpsl2)
}


library(psych)
library(GPArotation)

#dane z porównania:
data(Thurstone) #T 9x9 N= 213 prawie to
data("Holzinger") # HS 14x14 N = 355 to chyba to
data(Thurstone.33) #B/T 9x9 N=4175 problem z SL i SLiD
data("Bechtoldt") #T/B 17x17 N=213 to chyba to albo Bechtoldt.1
data("Harman74.cor") #H/H 24x24 N=145 ok
data("Reise") #R/M/H 16x16 N=35000 ok
data("Schmid") #Chen 18x18 N=403 jest idealnie tylko trzeba zamienić na macierz korelcji


label = "Reise"
data = Reise
R = cov2cor(data)

scree(R)
N = 5

Omega_h = label
Omega_t = label
ECV_ = label


lp = fa(R, nfactors = N + 1, rotate = "none")
L = lp$loadings


rot = perform_rotation(L, method = 'bigeominT')
Omega_h = paste(Omega_h, " & ", round(omega_h(rot$loadings, R),3))
Omega_t = paste(Omega_t, " & ", round(omega_t(rot$loadings, R),3))
ECV_ = paste(ECV_, " & ", round(ECV(rot$loadings),3))
print(paste('bigeominT omega_h: ', round(omega_h(rot$loadings, R),2)))
print(paste('bigeominT omega_t: ', round(omega_t(rot$loadings, R),2)))
print(paste('bigeominT ECV: ', round(ECV(rot$loadings),2)))



rot = perform_rotation(L, method = 'bifactorT')
Omega_h = paste(Omega_h, " & ", round(omega_h(rot$loadings, R),3))
Omega_t = paste(Omega_t, " & ", round(omega_t(rot$loadings, R),3))
ECV_ = paste(ECV_, " & ", round(ECV(rot$loadings),3))
print(paste('bifactorT omega_h: ', round(omega_h(rot$loadings, R),2)))
print(paste('bifactorT omega_t: ', round(omega_t(rot$loadings, R),2)))
print(paste('bifactorT ECV: ', round(ECV(rot$loadings),2)))

rot = perform_rotation(L, method = 'bifactorQ')
Omega_h = paste(Omega_h, " & ", round(omega_h(rot$loadings, R),3))
Omega_t = paste(Omega_t, " & ", round(omega_t(rot$loadings, R),3))
ECV_ = paste(ECV_, " & ", round(ECV(rot$loadings),3))
print(paste('bifactorQ omega_h: ', round(omega_h(rot$loadings, R),2)))
print(paste('bifactorQ omega_t: ', round(omega_t(rot$loadings, R),2)))
print(paste('bifactorQ ECV: ', round(ECV(rot$loadings),2)))

lp = fa(R, nfactors = N, rotate = "none")
L = lp$loadings

rot = SL(L)
Omega_h = paste(Omega_h, " & ", round(omega_h(rot, R),3))
Omega_t = paste(Omega_t, " & ", round(omega_t(rot, R),3))
ECV_ = paste(ECV_, " & ", round(ECV(rot),3))
print(paste('SL omega_h: ', round(omega_h(rot, R),2)))
print(paste('SL omega_t: ', round(omega_t(rot, R),2)))
print(paste('SL ECV: ', round(ECV(rot),2)))

rot = SLi(matrix = R, specific_factors = N)
Omega_h = paste(Omega_h, " & ", round(omega_h(rot$loadings, R),3))
Omega_t = paste(Omega_t, " & ", round(omega_t(rot$loadings, R),3))
ECV_ = paste(ECV_, " & ", round(ECV(rot$loadings),3))
print(paste('SLiD omega_h: ', round(omega_h(rot$loadings, R),2)))
print(paste('SLiD omega_t: ', round(omega_t(rot$loadings, R),2)))
print(paste('SLiD ECV: ', round(ECV(rot$loadings),2)))

Omega_h
Omega_t
ECV_



mean(c(0.865,     0.795,     0.809,        0.901,        0.881,   0.874,  0.854))
mean(c(0.838,     0.812,     0.765,        0.895,        0.881,   0.897,  0.874))
mean(c(0.828,      0.82,     0.758,        0.902,        0.887,   0.895,  0.865))
mean(c(0.715,     0.717,     0.628,                   0.664,   0.825,  0.781))
mean(c(0.715,     0.774,     0.713,                  0.835,    0.86,  0.841))


mean(c(0.942,     0.947,     0.903,        0.934,        0.942,   0.931,  0.958))
mean(c(0.942,     0.947,     0.903,        0.934,        0.942,   0.931,  0.958))
mean(c(0.909,     0.944,     0.908,        0.934,        0.937,    0.94,  0.947))
mean(c( 0.93,     0.949,       0.9,           0.937,    0.93,  0.954))
mean(c(0.941,     0.946,     0.903,            0.942,   0.931,  0.957))


mean(c(0.658,     0.456,     0.58,         0.772,        0.556,   0.721,  0.641))
mean(c(0.631,     0.472,     0.567,        0.77,         0.553,   0.741,  0.654))
mean(c(0.63,      0.486,     0.555,        0.775,        0.56,    0.737,  0.652))
mean(c(0.581,     0.417,     0.456,                0.451,   0.696,  0.602))
mean(c(0.548,     0.446,     0.529,                0.524,   0.709,  0.632))

