library(spdep)                               # Spatial Analysis; also opens libraries "maptools" and "sp"
library(RColorBrewer)                        # see "http://colorbrewer2.org/"
library(classInt)                            # functions of data classification
library(car)
library(maptools)


plotColorQual <- function(var.name,shape,my.title="",
                          my.legend=deparse(substitute(var.name)),
                          addToMap=F) {
  ##
  ## Plot a qualitative colors for factor "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  if (!is.factor(var.name)) stop("plotColorQual: Not a factor.")
  
  qualVal <- as.numeric(unclass(var.name))
  qualName <- levels(var.name)
  pal.Qual <- brewer.pal(12,"Set3")
  map.col <- pal.Qual[qualVal]
  
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend, legend=qualName,
         fill=pal.Qual[1:length(qualName)],bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotColorQual

plotColorRamp <- function(var.name,shape,n.breaks=8,my.title="",
                          my.legend=deparse(substitute(var.name)),
                          addToMap=F) {
  ##
  ## Plot a color ramp variable "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  
  ## define breaks and color assignment
  q.breaks <- classIntervals(var.name, n=n.breaks, style="quantile")
  pal.YlOrRd <- brewer.pal(n.breaks, "Oranges")
  #pal.YlOrRd <- brewer.pal(n.breaks, "YlOrRd")
  map.col <- pal.YlOrRd[findInterval(var.name,q.breaks$brks,rightmost.closed=T)]
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend,legend=leglabs(round(q.breaks$brks,digits=3)),
         fill=pal.YlOrRd,bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotColorRamp 

plotBiPolar <- function(var.name,shape,
                        neg.breaks=4,pos.breaks=neg.breaks,break.value=0,
                        my.title="",my.legend=deparse(substitute(var.name)),
                        addToMap=F) {
  ##
  ## Plot bipolar map theme for variable "var.name"
  ##
  require(spdep); require(RColorBrewer); require(classInt)
  
  ## define quantile breaks and color assignment
  q.neg.breaks <- classIntervals((var.name[var.name < break.value]), n=neg.breaks, style="quantile")
  q.pos.breaks <- classIntervals((var.name[var.name > break.value]), n=pos.breaks, style="quantile")
  q.breaks <- c(q.neg.breaks$brks[-(neg.breaks+1)],break.value,q.pos.breaks$brks[-1])     # combine neg and pos over zero
  
  pal.neg <- brewer.pal(neg.breaks, "Blues")
  pal.pos <- brewer.pal(pos.breaks, "Reds")
  pal <- c(rev(pal.neg),pal.pos)                                                # combine palettes
  
  map.col <- pal[findInterval(var.name,q.breaks,rightmost.closed=T)]
  ## generate choropleth map
  plot(shape,col=map.col,border=grey(0.9),axes=T,add=addToMap)
  legend("bottomleft", title=my.legend,legend=leglabs(round(q.breaks,digits=3)),
         fill=pal,bty="n",ncol=1)
  title(my.title)
  box()
} # end:plotBiPolar

setwd("G:\\Data\\QE\\TXCntyCrime2015")
TX.shp <- readShapePoly("TXCntyCrimeShp.shp",               
                        proj4string=CRS("+proj=longlat"))

TX.bbox <- bbox(TX.shp)                                        
plot(TX.shp,axes=T,col=grey(0.9),border="white",                
     xlim=TX.bbox[1,],ylim=TX.bbox[2,])                     
plotColorRamp(TX.shp$CRIMERATE,TX.shp,my.title="crimerate",
              my.legend="CRIMERATE", addToMap=T)   

scatterplotMatrix(~CRIMERATE+NODIPLOMA+UNEMPL+OWNHOUSE+POVERTY,data = TX.shp,id.n=2)

#lOG for NODIPLOMA
e1071::skewness(bcPower(TX.shp$NODIPLOMA,1))
e1071::skewness(bcPower(TX.shp$NODIPLOMA,0.5))
e1071::skewness(bcPower(TX.shp$NODIPLOMA,0))

#No for UNEMPL
TX.shp$UNEMPL<-TX.shp$UNEMPL+0.01
e1071::skewness(bcPower(TX.shp$UNEMPL,1))
e1071::skewness(bcPower(TX.shp$UNEMPL,0.5))
e1071::skewness(bcPower(TX.shp$UNEMPL,0))

#No for OWNHOUSE
e1071::skewness(bcPower(TX.shp$OWNHOUSE,1))
e1071::skewness(bcPower(TX.shp$OWNHOUSE,0.5))
e1071::skewness(bcPower(TX.shp$OWNHOUSE,0))


#0.5 for POVERTY
TX.shp$POVERTY<-TX.shp$POVERTY+0.01
e1071::skewness(bcPower(TX.shp$POVERTY,1))
e1071::skewness(bcPower(TX.shp$POVERTY,0.5))
e1071::skewness(bcPower(TX.shp$POVERTY,0))

lm1 <- lm(bcPower(CRIMERATE,0.5)~log(NODIPLOMA)+UNEMPL+OWNHOUSE+bcPower(POVERTY,0.5)+URBRURAL,data=TX.shp)
summary(lm1)
lm2 <- lm(bcPower(CRIMERATE,0.5)~log(NODIPLOMA)+UNEMPL+OWNHOUSE+bcPower(POVERTY,0.5),data=TX.shp)
anova(lm2,lm1)
vif(lm1)

#test the normalty of residual
car::qqPlot(lm1,id.n=2)
#determin whether need quadratic function
car::residualPlots(lm1, id.n=2)

#detect influential cases and outliers
car::influenceIndexPlot(lm1,id.n=4)
TX.df <- as.data.frame(TX.shp)
TX.shp<-TX.shp[-8,]
TX.shp<-TX.shp[-78,]
lmrefine <- lm(bcPower(CRIMERATE,0.5)~log(NODIPLOMA)+I(UNEMPL^2)+OWNHOUSE+bcPower(POVERTY,0.5)+URBRURAL,data=TX.shp)
summary(lmrefine)

vif(lmrefine)
car::qqPlot(lmrefine,id.n=2)
car::residualPlots(lmrefine, id.n=2)
car::influenceIndexPlot(lmrefine,id.n=4)


#############################################weighted model################
lmHetero <- function(X, ...){UseMethod("lmHetero")}                  # Declare function as S3 object

lmHetero.default <- function(formula, hetero=~1, data, subset, na.action, 
                             contrasts = NULL, iter=FALSE, ...) {
  ################################################################################################
  ##
  ## Purpose: Calculate multiplicately weighted regression models with respect to externally
  ##          given exogenous variables
  ##
  ## Source: The ML estimation procedure for multiplicately weighted regression is given
  ##         in Greene W. H. (2000). Econometric Analysis. 4th edition. Upper Saddle River:
  ##         Prentice Hall. pp 516-521 (Note: page numbers will differ for other editions)
  ##
  ## Objective: estimate gamma values for mulitiplicative heteroscedasticity models
  ##
  ## Syntax:
  ## [1] lmHetero(y~x1+x2, hetero= ~z1+z2, data=mydata)
  ## [2] lmHetero(y~x1+x2 | z1+z2, data=mydata)
  ## Note: An expression for Z must be present.
  ## y : Dependent variable 
  ## X : Independent variable(s) with added intercept
  ## Z : Weight variable(s) with intercept added.
  ##
  ##     !!! Important user input: weighte variables Z must be enter log-transformed !!!
  ##
  ## Authors: Michael Tiefelsdorf (tiefelsdorf@utd.edu) & Yongwan Chun
  ############################################################################################### 
  
  ## Parseing the function call
  require(Formula) 
  cl <- match.call()
  if (missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  if (!missing(hetero)) {
    formula <- as.Formula(formula, hetero)
    cl$hetero <- NULL
    cl$formula <- formula(formula)
  }
  else {
    formula <- as.Formula(formula)
  }
  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in% 1:2)
  has_dot <- function(formula) inherits(try(terms(formula), silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2))
      formula <- as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  y <- model.response(mf, "numeric")
  mt <- terms(formula, data = data)
  mtX <- terms(formula, data = data, rhs = 1)
  X <- model.matrix(mtX, mf, contrasts)
  if (length(formula)[2] < 2L) {
    #stop("Set of weights variables is missing in function call!")
    #mtZ <- NULL
    #Z <- NULL
    zNames <- "(Intercept)"
    Z <- matrix(1, nrow=length(y), ncol=1)
  }
  else {
    mtZ <- delete.response(terms(formula, data = data, rhs = 2))
    Z <- model.matrix(mtZ, mf, contrasts)
    zNames <- colnames(Z)
    if (!all(exp(Z) > 0)) stop("All weight variables must be positive!")
    #Z[,-1] <- log(Z[,-1])
  }
  ## done::parsing and starting the estimation
  
  ##
  ## Calculate start values for Beta and Gamma
  ##
  nofreg <- length(y)
  Beta <- qr.solve(crossprod(X),crossprod(X,y))
  res <- y - X %*% Beta
  res <- log(res^2)
  
  Gamma <- qr.solve(crossprod(Z),crossprod(Z,res))
  Gamma[1] <- Gamma[1]+1.2704    # Harvey (1976) correction
  if (iter==TRUE) { cat("Beta:",Beta);cat("\n");cat("Gamma:",Gamma);cat("\n") }
  
  qrcz <- qr(crossprod(Z))
  
  ## Start interative estimation procedure
  MaxDiff <- Inf
  while (MaxDiff > 0.001) {
    Sigma2 <- exp(Z %*% Gamma) 
    W <- diag(1 / Sigma2[,1])
    
    ## Calculate the values at the next step
    #  newBeta <- solve(t(x) %*% W %*% x, t(x) %*% W %*% y)
    xW <- crossprod(X,W)
    newBeta <- qr.solve(xW %*% X, xW %*% y)
    res2 <- (y - X %*% newBeta)^2
    newGamma <- Gamma + qr.coef(qrcz,crossprod(Z, res2 / Sigma2 - 1))        
    MaxDiff <- max(abs(newGamma - Gamma))
    if (iter==TRUE) {cat("Beta:",newBeta);cat("\n");cat("Gamma:",newGamma);cat("\n")}
    Beta <- newBeta
    Gamma <- newGamma
  }  # end while
  
  Sigma2 <- exp(Z %*% newGamma) 
  cBeta <- qr.solve(xW %*% X)      # The global covariance matrix is block diagonal
  cGamma <- 2*qr.solve(crossprod(Z))
  logLikeH1 <- -(nofreg/2)*log(2*pi) - 0.5*sum(log(Sigma2[,1])) - 0.5* sum(res2/Sigma2[,1])
  logLikeH0 <- logLik(lm(y~X))[1]
  ##
  ## Report results. As ML-estimates the values are slightly biased
  ##
  rval <- list()
  rval$sigma2 <- exp(newGamma[1]) 
  rval$gamma <- newGamma
  rval$namesGamma <- zNames
  rval$beta <- newBeta 
  rval$weights <- 1/Sigma2[,1]
  rval$covBeta <- cBeta
  rval$covGamma <- cGamma
  rval$logLikeH1 <- logLikeH1
  rval$logLikeH0 <- logLikeH0
  #   rval$call <- cl
  #   rval$formula <- formula(formula)
  #   rval$terms <- list(regressors = mtX, hetero = mtZ, full = mt)
  #   rval$na.action <- attr(mf, "na.action")
  #   rval$levels <- .getXlevels(mt, mf)
  #   rval$contrasts <- list(regressors = attr(X, "contrasts"),
  #                        hetero = attr(Z, "contrasts"))
  class(rval) <- "lmHetero"
  return(rval)
} # end::lmHetero
summary.lmHetero <- function(object, ...){
  ## Regression coefficent info
  b <- object$beta
  se <- sqrt(diag(object$covBeta))
  z <- b/se
  table1 <- cbind(b, se, z, 2*(1-pnorm(abs(z))))
  colnames(table1) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
  rownames(table1) <- rownames(object$beta)
  
  ## Gamma weights coefficient infor
  g <- object$gamma
  seg <- sqrt(diag(object$covGamma))
  zg <- g/seg
  table2 <- cbind(g, seg, zg, 2*(1-pnorm(abs(zg))))
  colnames(table2) <- c("Gamma", "Std.Err", "z-value", "Pr(>|z|)")
  rownames(table2) <- object$namesGamma
  
  result <- list(coef=table1, gamma=table2, 
                 logLikeH1=object$logLikeH1,logLikeH0=object$logLikeH0)
  class(result) <- "summary.lmHetero" 
  result
} # end::summary.lmHetero

print.summary.lmHetero <- function(X, ...){
  cat("\n===============================================================")
  cat("\nMultiplicatively Weighted Heteroscedasticity ML-Regrssion Model")
  cat("\n===============================================================\n\n")
  cat("Regression Coefficients:\n")
  printCoefmat(X$coef)
  cat("\nGamma Coefficients:\n")
  printCoefmat(X$gamma)
  cat("\nlog-likelihood =", X$logLikeH1,"\n")
  
  ## Likelihood ratio test for heteroscedasticity
  df <- nrow(X$gamma)-1                      # with df=dim(Gamma[-Intercept])
  if (df > 0L)
  {
    LR <- 2*abs(X$logLikeH1 - X$logLikeH0)     # chi^2 distributed 
    table <- cbind(LR,df,pchisq(LR,df,lower.tail=FALSE))
    colnames(table) <- c("LR","  df"," Pr(Chi > LR)")
    rownames(table) <- ""
    LRTest <- c("LR"=round(LR,2),"df"=round(df,0),"Pr(Chi > LR)"=round(pchisq(LR,df,lower.tail=FALSE),5))
    cat("\nHeteroscedasticity likelihood ratio test:\n")
    print(table)
  }
  cat("\n")
  invisible(X)
} # end::print.summary.lmHetero


lmheter <-lmHetero(bcPower(CRIMERATE,0.5)~log(NODIPLOMA)+I(UNEMPL^2)+OWNHOUSE+bcPower(POVERTY,0.5)+URBRURAL,hetero = ~log(POPCRIME),data=TX.shp)
summary(lmheter)


medResid <- residuals(lmrefine)
hist(medResid) 
plot(TX.shp,axes=T,col=grey(0.9),border="white",
     xlim=TX.bbox[1,],ylim=TX.bbox[2,])               
plotBiPolar(medResid,TX.shp,                           
            neg.breaks=4,pos.breaks=4,break.value=0.0, 
            my.title="CRIMERATE Model Residuals",
            my.legend="Residuals",addToMap=T)


tx.centroid <- coordinates(TX.shp)
tx.link <- poly2nb(TX.shp, queen=F)                          
plot(TX.shp,axes=T,col=grey(0.9),border="white",
     xlim=TX.bbox[1,],ylim=TX.bbox[2,])                      
plot(TX.shp,col="palegreen3" ,border=grey(0.9), axes=T, add=T) 
plot(tx.link,coords=tx.centroid, pch=19, cex=0.1,            
     col="blue", add=T)
title("Augmented Spatial Links among counties")                
box()   


tx.linkW <- nb2listw(tx.link, style="W")                  
outlierSpatial <- moran.plot(as.vector(residuals(lmrefine)), tx.linkW, labels=TX.shp$NAME) 
lm.morantest(lmrefine, tx.linkW)


refine.SAR <- spautolm(lmrefine, data=TX.shp,listw=tx.linkW, family="SAR")
summary(refine.SAR)

(likeH0 <- logLik(lmrefine))
(likeH1 <- logLik(refine.SAR))
cat("chi-square value:  ", chi <- -2*(likeH0[1]-likeH1[1]))
cat("error-probability: ", pchisq(chi, df=1, lower.tail=F))

moran.mc(residuals(refine.SAR), tx.linkW, nsim=10000) 
