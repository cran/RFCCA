## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message=FALSE-----------------------------------------------------------
library(RFCCA)
data(data)

set.seed(2345)
smp <- sample(dim(data$X)[1], round(dim(data$X)[1]*0.7))
train.X <- data$X[smp,]
train.Y <- data$Y[smp,]
train.Z <- data$Z[smp,]
test.Z <- data$Z[-smp,]

## -----------------------------------------------------------------------------
rf.obj <- rfcca(X = train.X, Y = train.Y, Z = train.Z, ntree = 100)

## -----------------------------------------------------------------------------
test.obj <- global.significance(X = train.X, Y = train.Y, Z = train.Z, 
                                ntree = 100, nperm = 10)
test.obj$pvalue

## -----------------------------------------------------------------------------
pred.oob <- rf.obj$predicted.oob
head(pred.oob)

## -----------------------------------------------------------------------------
vimp.obj <- vimp(rf.obj)
vimp <- vimp.obj$importance
vimp

## ---- fig.show='hold', fig.width=6, fig.height=4, fig.align='center'----------
plot.vimp(vimp.obj)

## -----------------------------------------------------------------------------
pred.obj <- predict(rf.obj, test.Z)
pred <- pred.obj$predicted
head(pred)

## -----------------------------------------------------------------------------
head(pred.obj$predicted.coef$coefx)
head(pred.obj$predicted.coef$coefy)

## -----------------------------------------------------------------------------
pred.obj2 <- predict(rf.obj, test.Z, finalcca = "scca")
pred2 <- pred.obj2$predicted
head(pred2)

## -----------------------------------------------------------------------------
head(pred.obj2$predicted.coef$coefx)
head(pred.obj2$predicted.coef$coefy)

## ---- fig.show='hold', fig.width=5, fig.height=5, fig.align='center'----------
plot(pred, pred2, xlab="Classical CCA", ylab="Sparse CCA", 
     xlim=c(min(pred,pred2),max(pred,pred2)), ylim=c(min(pred,pred2),max(pred,pred2)),
     pch = 20) 
abline(0,1,col = "red") 

## -----------------------------------------------------------------------------
sessionInfo()

