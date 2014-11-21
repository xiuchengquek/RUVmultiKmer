#Wrapper for RUVSeq Assessment of multiple unwanted Variation


###Introduction


Wrapper for model matrix and edgeR to allow users to build multiple design model
and conduct glm fitting based on the number of unwanted variation declared in
the normalization step in RUVSeq </br>


###Usage

```R

  source('model.matrix2_complete.R')
  source('conductMultipleDGE.R')

  # run your usually RUVSeq steps...


  RuvRset <- RUVr(set, genes, k=3, res)

  ## initialize your new "multipleExpression" object
  multipleSet <- new("multipleExpression", kmer=3, set=RuvRset, classification=x)
  multipleSet <- model.matrix.all(multipleSet)

  ##conduct DGE, return a list DGEGLM object in the slot "fits"
  multipleFit <- conductDGE(multipleSet)

  ## To access the fits list
  multiplefit@fits

```


### Implementation

Create a S4 object with the signature "multipleExpression". </br>
The s4 object will consist of the following slots
- kmer : the number of unwanted factors
- set : the SeqExpression Set
- model : list of design matrix
- fits : list of fitted models
- Classification : factors of classification

There are 2 S4 methods :
  - S4 Generic : "model.matrix.all"
  - S4 Generic : "conductDGE"


S4 Generic : "model.matrix.all"
  - Take a S4 Object with Signature "multipleExpression" and fill the slot "model"
  with a list of design matrix. The number of design matrix depends on the k
  declared in the normalization step

S4 Generic : "conductDGE"
  - Take a S4 Object with Signature "multipleExpression" , look at the slot "model"
    and conduct model fitting. Return a list of fitted model in the slot "fits"


  ```R

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


  ```


### Tests

Test have been down with RUnit package. Uses the zerbrafish package as the mock.

The test bascically just test that the output from a for loop and the S4 methods here
will produce the same results.


To run the test


```bash
   Rscript run_test.R
```
