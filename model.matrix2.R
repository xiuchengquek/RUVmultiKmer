

## multipleExpression is a class object with the following slot
setClass("multipleExpression", representation(
                                  kmer="numeric",
                                  set="SeqExpressionSet",
                                  model="list",
                                  fits="list",
                                  classification='numeric'
                                ));

setGeneric('model.matrix.all', function(set) standardGeneric("model.matrix.all"));

setMethod('model.matrix.all', c(set="multipleExpression"),

  function(set){

    ##function of constructing fitting
    constructFitting <- function (df, coefficient){
      ## construct basic design matrix without including the cofficient
      baseModel <- model.matrix(~0 + x, data=df);

      ## add in the coefficient
      reBasedModel <- cbind(baseModel, df[2:coefficient]);
      assignAttr <- attr(baseModel, 'assign');
      attr(reBasedModel, 'assign') <- c(assignAttr, 2:coefficient);
      return(reBasedModel)
    }

    x <- pData(set);
    columns <- ncol(x);
    coefficients <- 2:columns;

    ## return a list of fits
    designMatrices = lapply(coefficients, constructfitting);

    slot(set, 'model', check=TRUE) = designMatrices;
    return(set)
})

setGeneric('conductDGE', function(set) standardGeneric('conductDGE"));

setMethod('conductDGE', c(set='multipleExpression'),

  function(set) {

    singleDGE <- function(design, set){
      counts <- counts(set@set);
      groups <- set@classification;
      y <- DGEList(counts=counts, group=groups);
      y <- calcNormFactors(y, method="upperquartile")
      y <- estimateGLMCommonDisp(y, design)
      y <- estimateGLMTagwiseDisp(y, design)
      fit <- glmFit(y, design)
      return(fit)
    }

    if (lenslot('multipleExpression', 'model') == 0){
      warn('model slot not defined')
    }
    else {
      fits = lapply(set@model, singleDGE, set=set@set);
      slot(set, 'fits', check=TRUE) <- fits;
      return(set)
    }
  }

)
