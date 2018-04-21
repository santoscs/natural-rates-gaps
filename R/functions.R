# Functions

ApplyOLS <- function(mY, mX) {
  ct  <- nrow(mY)
  mxx <- crossprod(mX)
  mxy <- crossprod(mX, mY)
  mb  <- t(solve(mxx, mxy))
  me  <- mY - mX %*% t(mb)
  mee <- crossprod(me)
  ms  <- mee / ct
  return(list(mB = mb, mXX = mxx, mEE = mee, mS = ms))
}

SpecifyPrior <- function(mY, iP, dLambda, dN_0) {
  cn   <- ncol(mY)
  mm_0 <- matrix(0, cn, iP * cn)
  vs2  <- NULL
  for (i in 1:cn) {
    dy.ar <- ar(mY[, i], FALSE, iP, method = "ols", demean = FALSE)
    vs2   <- c(vs2, dy.ar$var.pred)
  }
  md_0 <- (1 / dLambda ^ 2) * diag((1:iP) ^ 2) %x% diag(vs2)
  ms_0 <- (dN_0 - cn - 1) * diag(vs2)
  return(list(mM_0 = mm_0, mD_0 = md_0, dN_0 = dN_0, mS_0 = ms_0))
}

SpecifyPosterior <- function(mM_0, mD_0, dN_0, mS_0, mY, mX) {
  ct       <- nrow(mY)
  lols     <- ApplyOLS(mY, mX)
  mphi_OLS <- lols$mB
  mxx      <- lols$mXX
  mee      <- lols$mEE
  mxxInv   <- solve(mxx)
  md_1     <- mxx + mD_0
  md_1Inv  <- solve(md_1)
  mm_1     <- (mphi_OLS %*% mxx + mM_0 %*% mD_0) %*% md_1Inv
  dn_1     <- ct + dN_0
  md_0Inv  <- solve(mD_0)
  ms_1     <- (mphi_OLS - mM_0) %*% solve(mxxInv + md_0Inv) %*% t(mphi_OLS - mM_0) + mee + mS_0
  return(list(mM_1 = mm_1, mD_1 = md_1, dN_1 = dn_1, mS_1 = ms_1))
}

LogMultiGamma <- function(cN, dAlpha) {
  (cN * (cN - 1) / 4) * log(pi) + sum(lgamma(dAlpha - (1: cN - 1) / 2))
}

ComputeLogMarginalLikelihood <- function(mD_0, dN_0, mS_0, mD_1, dN_1, mS_1) {
  cn          <- nrow(mS_0)
  dlprior     <- -(dN_0 / 2) * log(det(mS_0 / 2)) - (cn / 2) * log(det(mD_0)) + LogMultiGamma(cn, dN_0 / 2)
  dlposterior <- -(dN_1 / 2) * log(det(mS_1 / 2)) - (cn / 2) * log(det(mD_1)) + LogMultiGamma(cn, dN_1 / 2)
  dlml        <- dlposterior - dlprior
  return(dlml)
}

ComputeSDRatio <- function(mD_0Inv, dN_0, mS_0, mM_1, mD_1Inv, dN_1, mS_1) {
  cn    <- nrow(mM_1)
  ip    <- ncol(mM_1) / cn
  mR    <- diag(c(t(matrix(1, 1, ip) %x% (matrix(1, cn, cn) - diag(cn)))))
  mR    <- unique(mR)[-1,]
  vmu_1 <- c(t(mM_1))
  vx_1  <- c(mR %*% vmu_1)
  vx_0  <- 0 * vx_1
  dp_0  <- 0
  dp_1  <- 0
  cR    <- 1000
  for (r in 1:cR) {
    msigma_0 <- riwish(dN_0, mS_0)
    msigma_1 <- riwish(dN_1, mS_1)
    dp_0     <- dp_0 + dmvnorm(vx_0,, mR %*% (msigma_0 %x% mD_0Inv) %*% t(mR))
    dp_1     <- dp_1 + dmvnorm(vx_1,, mR %*% (msigma_1 %x% mD_1Inv) %*% t(mR))
  }
  return(dp_1 / dp_0)
}

GetCompanionMatrix <- function(mPhi) {
  cn <- nrow(mPhi)
  ip <- ncol(mPhi) / cn
  ma <- rbind(mPhi, cbind(diag((ip - 1) * cn), matrix(0, (ip - 1) * cn, cn)))
  return(ma)
}

ComputeInitStateVariance <- function(mPhi, mSigma) {
  ma     <- GetCompanionMatrix(mPhi)
  cm     <- nrow(ma)
  cn     <- nrow(mPhi)
  mb     <- rbind(t(chol(mSigma)), matrix(0, cm - cn, cn))
  mc     <- diag(cm ^ 2) - ma %x% ma
  #	mc     <- Matrix(mc)
  mgamma <- matrix(solve(mc) %*% c(mb %*% t(mb)), cm, cm)
  return(mgamma)
}

DrawPhi <- function(mSigma, mM_1, mD_1Inv) {
  cn          <- nrow(mM_1)
  ip          <- ncol(mM_1) / cn
  ck          <- ip * cn ^ 2
  vmean       <- c(t(mM_1))
  mvariance   <- mSigma %x% mD_1Inv
  dlambda_max <- 1
  while (dlambda_max >= 1) {
    vphi        <- rmvnorm(1, vmean, mvariance)
    mphi        <- t(matrix(vphi, ip * cn, cn))
    ma          <- GetCompanionMatrix(mphi)
    ma.eigen    <- eigen(ma)
    vlambda     <- ma.eigen$values
    dlambda_max <- max(abs(vlambda))
  }
  return(mphi)
}

ComputePmove <- function(mGamma, mGamma_star, vS_0) {
  dlnf0      <- dmvnorm(vS_0,, mGamma,      log = TRUE)
  dlnf0_star <- dmvnorm(vS_0,, mGamma_star, log = TRUE)
  dalpha     <- min(exp(dlnf0_star - dlnf0), 1)
  return(dalpha)
}

lnML <- function(x) {
  dlambda    <- x[1]
  dn_0       <- x[2]
  #	dlambda    <- x
  #	dn_0       <- cN + 2
  lprior     <- SpecifyPrior(mdY, iP, dlambda, dn_0)
  mm_0       <- lprior$mM_0
  md_0       <- lprior$mD_0
  dn_0       <- lprior$dN_0
  ms_0       <- lprior$mS_0
  lposterior <- SpecifyPosterior(mm_0, md_0, dn_0, ms_0, my, mx)
  md_1       <- lposterior$mD_1
  dn_1       <- lposterior$dN_1
  ms_1       <- lposterior$mS_1
  dlml       <- ComputeLogMarginalLikelihood(md_0, dn_0, ms_0, md_1, dn_1, ms_1)
  return(-dlml)
}

#' Plota series temporais com ggplot2
#' 
#' @param x objeto ts ou mts com as series temporais
#' @param y (opicional) objeto ts ou mts com dimensao de x para ser
#' plotado junto com x no mesmo grafico 
#' @param escala Are scales shared across all facets
#'  ("fixed"), or do they vary across 
#'  rows ("free_x"), columns (the default, "free_y"), or both 
#'  rows and columns ("free")
#' @param facet as series sao plotadas em graficos diferente (facet = TRUE, the default),
#' ou no mesmo grafico (facet = FALSE)
#' @param name optional name for ts univariate
#' 
#' @return ggplot das series
#' 
#' @import ggplot2 zoo
#' 
#' @export
#' 

tsplot <- function(x, y = NULL, escala = 'free_y', facet = TRUE, name = NULL){
  nseries <- NCOL(x)
  ntime <- NROW(x)
  x <- zoo::as.zoo(x)
  df.x <- zoo::fortify.zoo(x, melt = TRUE)
  if(nseries==1 & !is.null(name)){
    df.x[,"Series"] <- rep(name, ntime)
  }
  if(!is.null(y)){
    y <- zoo::as.zoo(y)
    df.y <- zoo::fortify.zoo(y, melt = TRUE)
    if(facet){
      df <- ggplot2::fortify(cbind(df.x, Value2=df.y[,3]), index.name = "Index")
      p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = Value))
      p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Index, y = Value2),
                                  linetype=2, colour="red", size = 1/2, alpha = 1)
      p <- p + ggplot2::geom_line(size = 1/2, alpha = 1, colour="blue")  
      p <- p + ggplot2::facet_grid(Series ~ ., scales = "free_y") 
      #p <- p + ggplot2::facet_wrap(~ Series, scales = "free_y")
    }else{
      p <- ggplot(df.x, aes(x = Index, y = Value))
      p <- p + geom_line(data = df.y, aes(x = Index, y = Value, group = Series, size = Series), size = 1/2, alpha = 1, colour="blue")
      p <- p + geom_line(linetype=2, size = 1/2, alpha = 1, colour="blue")  # Drawing the "overlayer"
    }
    p <- p + ggplot2::labs(y="", x="")
    p <- p + ggplot2::theme_bw(base_size=12)
    return(p)  
  }
  if(!facet){
    p <- ggplot2::ggplot(data = df.x, ggplot2::aes(x = Index, y = Value, color=Series, linetype=Series))
    p <- p + ggplot2::geom_line(size = 3/4)  
    p <- p + ggplot2::labs(y="", x="")
    p <- p + ggplot2::theme_bw(base_size=12)
    return(p)  
  }else{
    p <-ggplot2::ggplot(df.x, ggplot2::aes(x=Index, y=Value, group_by())) +
      ggplot2::geom_line(size = 1/2, alpha = 1, colour="blue") +
      ggplot2::facet_grid(Series ~ ., scales = escala) +
      ggplot2::labs(y="", x="") +
      ggplot2::theme_bw(base_size=12)
  }
  return(p)
}


#' Intervalo de confianca 
#' 
#' Calcula o intervalo de confianca para uma estimativa pontual
#' 
#' @param x estimativa pontual
#' @param dp desvio padrao da estimativa pontual
#' @param nc nivel de confianca do intervalo, entre 0 e 1
#' 
#' @return uma lista com 
#' \item{\code{upper}}{limite superior do intervalo}
#' \item{\code{lower}}{limite inferior do intervalo}
#' 
#' @export
#'    

ic <- function(x, dp, nc = 0.95){
  if(nc >= 1 | nc <=0){
    stop("nc fora do intervalo 0 < nc < 1")
  }
  n <- dim(x)[1]
  # margem de erro
  tn <- abs(qt(((1-nc)/2), df = n))
  upper <- x + tn*dp
  lower <- x - tn*dp
  return(list(upper = upper, lower = lower))
}

#' Plot com ggplot2
#' 
#' @param x series a serem plotadas
#' @param y opcional serie ou series a serem plotadas juntas para comparacao
#' @param upper limite superior para regiao sombreada
#' @param lower limite inferior para regiao sombreada
#' 
#' @return ggplot das series
#' 
#' @import ggplot2 zoo
#' 
#' @export
#' 

bitsplot <- function(x, y=NULL, upper, lower){
  df.x <- zoo::fortify.zoo(x, melt = TRUE)
  ic.upper <- zoo::fortify.zoo(upper, melt = TRUE)
  ic.lower <- zoo::fortify.zoo(lower, melt = TRUE)
  if(is.null(y)){
    df <- ggplot2::fortify(cbind(df.x, upper=ic.upper[,3], lower=ic.lower[,3]), index.name = "Index")
  }else{
    df.y <- zoo::fortify.zoo(y, melt = TRUE)
    df <- ggplot2::fortify(cbind(df.x, Value2=df.y[,3], upper=ic.upper[,3], lower=ic.lower[,3]), index.name = "Index")
  }
  
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = Value))
  if(!is.null(y)){
    p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Index, y = Value2),
                                linetype=2, size = 1/2, color = "red")
  }
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),
                                alpha=0.3)
  p <- p + ggplot2::geom_line(size = 1/2, color = "blue")
  p <- p + ggplot2::facet_grid(Series ~ ., scales = "free_y") 
  p <- p + ggplot2::labs(y="", x="")
  p <- p + ggplot2::theme_bw()
  return(p)  
}




#' ts theme ggplot2
#'
#' ts theme set the general aspect of the plot such as 
#' the colour of the background, gridlines, the size and colour of fonts
#' 
#' @param base_size base font size
#' @param base_family base font family
#' 
#' @details theme_ts is based in the classic dark-on-light ggplot2 theme. 
#' May work better for presentations displayed with a projector
#' 
#' @examples 
#' \dontrun{ 
#' p <- ggplot(mtcars) + geom_point(aes(x = wt, y = mpg, 
#' colour=factor(gear))) + facet_wrap(~am)
#' p
#' p + theme_ts()
#' }
#' 
#' @import ggplot2
#' @export
#' 

theme_ts <- function (base_size = 12, base_family = "") 
{
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.text = element_text(size = rel(0.8)), 
          strip.text = element_text(size = rel(0.9)),
          axis.ticks = element_line(colour = "black"),
          legend.key = element_rect(colour = "grey80"), 
          panel.background = element_rect(fill = "white", colour = NA), 
          panel.border = element_rect(fill = NA, colour = "grey50"), 
          panel.grid.major = element_line(colour = "grey88", size = 0.2),
          panel.grid.minor = element_line(colour = "grey98", size = 0.5), 
          strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
}



DrawScatter <- function(x, name1, name2){
  # nlim <- round(max(max(x), abs(min(x))), digits=1) 
  # nlim <- nlim +0.05*nlim
  x <- x[,c(name1, name2)]
  p <- ggplot(as.data.frame(x), aes(x=x[,1], y=x[,2])) +
    geom_point(shape=21, size=0.5) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, colour="grey") + geom_vline(xintercept = 0, colour="grey") +
    xlab(name1) + ylab(name2)
    # xlim(c(-nlim,nlim)) + ylim(c(-nlim,nlim))
  return(p)
}

DrawPDF <- function(x, name1){
  p <- ggplot(as.data.frame(x), aes(x=x)) +
    geom_density() +
    theme_bw()  +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_hline(yintercept = 0, colour="grey") + geom_vline(xintercept = 0, colour="grey") +
    xlab(name1)
  return(p)
}


