library(sinatra)
library(varbvs)

WT_data <- data.matrix(read.table("D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_WT_65_213_4_15_01_4_200.txt",row.names = NULL))
L67S_data <- data.matrix(read.table("D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_R164S_65_213_4_15_01_4_200.txt",row.names = NULL))

dimnames(WT_data) <- NULL
dimnames(L67S_data) <- NULL
X <- rbind(WT_data,L67S_data)

Xmean=apply(X, 2, mean)
Xsd=apply(X, 2, sd)
X[is.na(X)] <- 0

write.table(X[1:101,], file = "D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_WT_65_213_4_15_01_4_200_norm_R164S.txt", sep = " ",row.names = FALSE, col.names = FALSE)
write.table(X[102:202,], file = "D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_R164S_65_213_4_15_01_4_200_norm_WT.txt", sep = " ",row.names = FALSE, col.names = FALSE)

y <- c(rep(0,101),rep(1,101))
fit <- varbvs(X=X,Z=NULL,y=y)
print(summary(fit))
output <- fit[["pip"]]

write.table(output, file = "D:/share/beta_lactamase_test/gabriel_data/pip_WT_65_213_R164S_ect_nofaces_200.txt", sep = " ",row.names = FALSE, col.names = FALSE)




####
#### Control
####

WT_data <- data.matrix(read.table("D:/share/ubq/PH/sect_L67S_3_ph_noh_unit_15_01_4_25.txt",row.names = NULL))
dimnames(WT_data) <- NULL
X <- WT_data
Xmean=apply(X, 2, mean)
Xsd=apply(X, 2, sd)
X=t((t(X)-Xmean)/Xsd)
X[is.na(X)] <- 0
#X <- rbind(WT_data,L67S_data)
y <- c(rep(0,51),rep(1,50))
fit <- varbvs(X=X,Z=NULL,y=y)
print(summary(fit))

output <- fit[["pip"]]
write.table(output, file = "D:/share/ubq/PH/pip_noh_norm_L67S_3_control.txt", sep = " ",row.names = FALSE, col.names = FALSE)
