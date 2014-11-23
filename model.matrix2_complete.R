

## multipleExpression is a class object with the following slot
setClass("multipleExpression", representation(
                                  set="SeqExpressionSet",
                                  model="list",
                                  fits="list",
                                  classification='factor'
                                ));

setGeneric('model.matrix.all', function(set) standardGeneric("model.matrix.all"));

setMethod('model.matrix.all', c(set="multipleExpression"),

  function(set){

    constructFitting <- function (coefficient, df){
      ## construct basic design matrix without including the cofficient
      baseModel <- model.matrix(~0 + x, data=df);

      ## add in the coefficient
      message(df)
      coefficient_values <- as.matrix(df[2:coefficient])
      reBasedModel <- cbind(baseModel,  coefficient_values)

      assignAttr <- attr(baseModel, 'assign');
      contrastAttr <- attr(baseModel, 'contrasts')

      attr(reBasedModel, 'assign') <- c(assignAttr, 2:coefficient);
      attr(reBasedModel, 'contrasts') <- contrastAttr
      return(reBasedModel)
    }

    ##function of constructing fittingncol


    x <- pData(set@set);
    columns <- ncol(x);
    coefficients <- 2:columns;

    ## return a list of fits
    designMatrices = lapply(coefficients, constructFitting, df = x);

    slot(set, 'model', check=TRUE) = designMatrices;
    return(set)
})
