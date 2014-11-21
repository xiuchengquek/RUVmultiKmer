library(RUnit)
library(RUVSeq)
library(zebrafishRNASeq)
library(RColorBrewer)
## load the actual script
data(zfGenes)

source('model.matrix2_complete.R')
source('conductMultipleDGE.R')


testsuite.testSuiteMutliple <- defineTestSuite("testSuiteMutliple",
                  dir=file.path('test'), testFileRegexp = '*.R')

test.result <- runTestSuite(testsuite.testSuiteMutliple)
printTextProtocol(test.result)
