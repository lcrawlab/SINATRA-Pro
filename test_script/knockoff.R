library(knockoff)

WT_data <- data.matrix(read.table("D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_WT_65_213_4_15_01_4_100.txt",row.names = NULL))
L67S_data <- data.matrix(read.table("D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_R164S_65_213_4_15_01_4_100.txt",row.names = NULL))

dimnames(WT_data) <- NULL
dimnames(L67S_data) <- NULL
X <- rbind(WT_data,L67S_data)

Xmean=apply(X, 2, mean)
Xsd=apply(X, 2, sd)
X[is.na(X)] <- 0

write.table(X[1:101,], file = "D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_WT_65_213_4_15_01_4_100_norm_R164S.txt", sep = " ",row.names = FALSE, col.names = FALSE)
write.table(X[102:202,], file = "D:/share/beta_lactamase_test/gabriel_data/ect_nofaces_R164S_65_213_4_15_01_4_100_norm_WT.txt", sep = " ",row.names = FALSE, col.names = FALSE)

y <- c(rep(0,101),rep(1,101))

result = knockoff.filter(X, y)
