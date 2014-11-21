setGeneric('conductDGE', function(set) standardGeneric('conductDGE'));

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

    if (length(slot(set, 'model')) == 0){
      warning('model slot not defined')
    }
    else {
      fits <- lapply(set@model, singleDGE, set=set);
      slot(set, 'fits', check=TRUE) <- fits;
      return(set)
    }
  }

)
