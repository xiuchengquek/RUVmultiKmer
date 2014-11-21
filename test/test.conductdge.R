
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
  testset <- new("multipleExpression", kmer=3, set=mock_set4, classification=x)
  testset <- model.matrix.all(testset)

  #
  counts <- counts(mock_set4);
  groups <- x;
  y <- DGEList(counts=counts, group=groups);
  y <- calcNormFactors(y, method="upperquartile")
  y1 <- estimateGLMCommonDisp(y, k_1)
  y1 <- estimateGLMTagwiseDisp(y1, k_1)
  fit1 <- glmFit(y1, k_1)

  y <- DGEList(counts=counts, group=groups);
  y <- calcNormFactors(y, method="upperquartile")
  y2 <- estimateGLMCommonDisp(y, k_2)
  y2 <- estimateGLMTagwiseDisp(y2, k_2)
  fit2 <- glmFit(y2, k_2)

  y <- DGEList(counts=counts, group=groups);
  y <- calcNormFactors(y, method="upperquartile")
  y3 <- estimateGLMCommonDisp(y, k_3)
  y3 <- estimateGLMTagwiseDisp(y3, k_3)
  fit3 <- glmFit(y3, k_3)


  testfit <- conductDGE(testset)

  test.conductDGE <- function(){
    checkEquals(fit1 , testfit@fits[[1]])
    checkEquals(fit2 , testfit@fits[[2]])
    checkEquals(fit3 , testfit@fits[[3]])

  }
