library('SNFtool')

data1 = read.delim('matrices/gene_exp_matrix.tsv',sep = '\t')
data2 = read.delim('matrices/gene_cn_matrix.tsv.gz',sep = '\t')
data3 = read.delim('matrices/gene_fusion_matrix.tsv',sep = '\t')

data.labels = read.delim('annotations/fhr-annotations.tsv')

data = Reduce(function(...) merge(...,by="PUBLIC_ID"), list(data.labels, data1, data2))

truelabel = c("-1"="NA","0"="SR","1"="GHR","2"="SR")[as.character(data$risk)]

matrix1 = standardNormalization(as.matrix(data[, 8:18632]))
matrix2 = standardNormalization(as.matrix(data[,18634:37420]))
matrix3 = standardNormalization(as.matrix(data[,37422:]))

Dist1 = (dist2(matrix1, matrix1))^(1/2)
Dist2 = (dist2(matrix2, matrix2))^(1/2)

W1 = affinityMatrix(Dist1, K=20, sigma=0.5)
W2 = affinityMatrix(Dist2, K=20, sigma=0.5)

displayClusters(W1, truelabel)
displayClusters(W2, truelabel)

W = SNF(list(W1,W2),K=20,t=20)

displayClusters(W, truelabel)

