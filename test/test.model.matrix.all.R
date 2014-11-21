library(RUnit)
library(RUVSeq)
library(zebrafishRNASeq)
library(RColorBrewer)
## load the actual script
data(zfGenes)

## setUp 'mock objects'


filter <- apply(zfGenes, 1, function(x) length(x[x>5]))
filtered <- zfGenes[filter, ]
genes <- rownames(filtered)[grep("^ENS", rownames(filtered))]
spikes <- rownames(filtered)[grep("^ERCC", rownames(filtered))]

x <- as.factor(rep(c("Ctl", "Trt"), each=3))
set <-newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x,
  row.names=colnames(filtered)))

design <- model.matrix(~x, data=pData(set))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

mock_set4 <- RUVr(set, genes, k=3, res)

# create model for each of the kmers
k_1 <- model.matrix(~0 + x + W_1, data=pData(mock_set4))
k_2 <- model.matrix(~0 + x + W_1 + W_2, data=pData(mock_set4))
k_3 <- model.matrix(~0 + x + W_1 + W_2 + W_3, data=pData(mock_set4))

## initialize our class
testset <- new("multipleExpression", set=mock_set4, classification=x)
testset <- model.matrix.all(testset)

#
test.model.matrix.all <- function() {
  checkEquals(testset@model[[1]], k_1);
  checkEquals(testset@model[[2]], k_2);
  checkEquals(testset@model[[3]], k_3);
}
