# setwd("C:/Users/grzeg/Desktop/")
# 
#HEXACO = read.table("data.csv", sep = '\t', header = T)
# 
# length(unique(HEXACO$country))
# 
# HEXACO = HEXACO[,1:(ncol(HEXACO)-4)]
# 
# 
# 
# Honesty_Humility = HEXACO[,1:40]
# Emotionality = HEXACO[,41:80]
# eXtraversion = HEXACO[,81:120]
# Agreeableness = HEXACO[,121:160]
# Conscientiousness = HEXACO[,161:200]
# Openness = HEXACO[,201:240]
# 


#HEXACO = HEXACO[,(1:120)*2]

Honesty_Humility = HEXACO[,1:20]
Emotionality = HEXACO[,21:40]
eXtraversion = HEXACO[,41:60]
Agreeableness = HEXACO[,61:80]
Conscientiousness = HEXACO[,81:100]
Openness = HEXACO[,101:120]





label = "Openness"
data = Openness
R = cor(data)


N = 4

Omega_h = label
Omega_t = label
ECV_ = label

l1 = fa(R, nfactors = 1, rotate = "none")
keys = ifelse(l1$loadings >= 0, 1, -1)

for( i in 1:ncol(R)){
  for(j in 1:ncol(R)){
    R[i,j] = R[i,j] * keys[i] * keys[j]
  }
}


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

