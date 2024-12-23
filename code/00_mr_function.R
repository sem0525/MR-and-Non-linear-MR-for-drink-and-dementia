#  Perform 2 sample IV using Wald ratio.
mr_wald_ratio <- function(SNP, b_exp, b_out, se_exp, se_out)
{
  b <- b_out / b_exp
  se <- se_out / abs(b_exp)
  # sqrt((segd^2/gp^2) + (gd^2/gp^4)*segp^2 - 2*(gd/gp^3)) #full delta method with cov set to 0
  pval <- stats::pnorm(abs(b) / se, lower.tail=FALSE) * 2
  return(data.frame(SNP=SNP,b=b, se=se, pval=pval, nsnp=1))
}

#' Perform 2 sample IV using simple standard error
mr_meta_fixed_simple <- function(b_exp, b_out, se_exp, se_out)
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
  {
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  }
  b <- sum(b_exp*b_out / se_out^2) / sum(b_exp^2/se_out^2)
  se <- sqrt(1 / sum(b_exp^2/se_out^2))
  pval <- 2 * stats::pnorm(abs(b) / se, lower.tail=FALSE)
  return(data.frame(b=b, se=se, pval=pval, nsnp=length(b_exp)))
}

#' Perform 2 sample IV using fixed effects meta analysis and delta method for standard errors
mr_meta_fixed <- function(b_exp, b_out, se_exp, se_out, Cov=0)
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 1)
  {
    return(list(b=NA, se=NA, pval=NA, nsnp=NA, Q =NA, Q_df =NA, Q_pval =NA))
  }
  ratio <- b_out / b_exp
  ratio.se <- sqrt((se_out^2/b_exp^2) + (b_out^2/b_exp^4)*se_exp^2 - 2*(b_out/b_exp^3)*Cov)
  res <- metagen(ratio, ratio.se) #!!!!
  b <- res$TE.fixed
  se <- res$seTE.fixed
  pval <- res$pval.fixed
  Q_pval <- stats::pchisq(res$Q, res$df.Q, lower.tail=FALSE)
  return(list(b=b, se=se, pval=pval, nsnp=length(b_exp), Q = res$Q, Q_df = res$df.Q, Q_pval = Q_pval))
}



#' Maximum likelihood MR method
mr_two_sample_ml <- function(b_exp, b_out, se_exp, se_out)
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
  {
    return(data.frame(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
  }
  loglikelihood <- function(param) {
    return(1/2*sum((b_exp-param[1:length(b_exp)])^2/se_exp^2)+1/2*sum((b_out-param[length(b_exp)+1]*param[1:length(b_exp)])^2/se_out^2))
  }
  opt <- try(stats::optim(
    c(b_exp, sum(b_exp*b_out/se_out^2)/sum(b_exp^2/se_out^2)),
    loglikelihood,
    hessian=TRUE,
    control = list(maxit=25000)), silent=TRUE)
  if(inherits(opt, "try-error"))
  {
    message("mr_two_sample_ml failed to converge")
    return(data.frame(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
  }
  
  b <- opt$par[length(b_exp)+1]
  se <- try(sqrt(solve(opt$hessian)[length(b_exp)+1,length(b_exp)+1]))
  if(inherits(se, "try-error"))
  {
    message("mr_two_sample_ml failed to converge")
    return(data.frame(b=NA, se=NA, pval=NA, nsnp=NA, Q=NA, Q_df=NA, Q_pval=NA))
  }
  
  pval <- 2 * stats::pnorm(abs(b) / se, lower.tail=FALSE)
  
  Q <- 2 * opt$value
  Q_df <- length(b_exp) - 1
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
  
  return(data.frame(b=b, se=se, pval=pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

#' Egger's regression for Mendelian randomization
mr_egger_regression <- function(b_exp, b_out, se_exp, se_out)
{
  stopifnot(length(b_exp) == length(b_out))
  stopifnot(length(se_exp) == length(se_out))
  stopifnot(length(b_exp) == length(se_out))
  
  nulllist <- list(
    b = NA,
    se = NA,
    pval = NA,
    nsnp = NA,
    b_i = NA,
    se_i = NA,
    pval_i = NA,
    Q = NA,
    Q_df = NA,
    Q_pval = NA,
    mod = NA,
    smod = NA,
    dat = NA
  )
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
  {
    return(nulllist)
  }
  
  sign0 <- function(x)
  {
    x[x==0] <- 1
    return(sign(x))
  }
  
  to_flip <- sign0(b_exp) == -1
  b_out = b_out*sign0(b_exp)
  b_exp = abs(b_exp)
  dat <- data.frame(b_out=b_out, b_exp=b_exp, se_exp=se_exp, se_out=se_out, flipped=to_flip)
  mod <- stats::lm(b_out ~ b_exp, weights=1/se_out^2)
  smod <- summary(mod)
  if(nrow(stats::coefficients(smod)) > 1)
  {
    b <- stats::coefficients(smod)[2,1]
    se <- stats::coefficients(smod)[2,2] / min(1,smod$sigma)
    pval <- 2 * stats::pt(abs(b / se), length(b_exp) - 2, lower.tail = FALSE)
    b_i <- stats::coefficients(smod)[1,1]
    se_i <- stats::coefficients(smod)[1,2] / min(1,smod$sigma)
    pval_i <- 2 * stats::pt(abs(b_i / se_i), length(b_exp) - 2, lower.tail = FALSE)
    
    Q <- smod$sigma^2 * (length(b_exp) - 2)
    Q_df <- length(b_exp) - 2
    Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
  } else {
    warning("Collinearities in MR Egger, try LD pruning the exposure variables.")
    return(nulllist)
  }
  return(list(b = b, se = se, pval = pval, nsnp = length(b_exp), b_i = b_i, se_i = se_i, pval_i = pval_i, Q = Q, Q_df = Q_df, Q_pval = Q_pval, mod = smod, dat = dat))
}

#' Weighted median method
#'
#' Perform MR using summary statistics. Bootstraps used to calculate standard error. 
mr_weighted_median <- function(b_exp, b_out, se_exp, se_out,nboot = 1000)
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 3)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  b_iv <- b_out / b_exp
  VBj <- ((se_out)^2)/(b_exp)^2 + (b_out^2)*((se_exp^2))/(b_exp)^4
  b <- weighted_median(b_iv, 1 / VBj)
  se <- weighted_median_bootstrap(b_exp, b_out, se_exp, se_out, 1 / VBj, nboot)
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
  return(data.frame(b=b, se=se, pval=pval, Q=NA, Q_df=NA, Q_pval=NA, nsnp=length(b_exp)))
}
weighted_median <- function(b_iv, weights)
{
  betaIV.order <- b_iv[order(b_iv)]
  weights.order <- weights[order(b_iv)]
  weights.sum <- cumsum(weights.order)-0.5*weights.order
  weights.sum <- weights.sum/sum(weights.order)
  below <- max(which(weights.sum<0.5))
  b = betaIV.order[below] + (betaIV.order[below+1]-betaIV.order[below])*
    (0.5-weights.sum[below])/(weights.sum[below+1]-weights.sum[below])
  return(b)
}

weighted_median_bootstrap <- function(b_exp, b_out, se_exp, se_out, weights, nboot)
{
  med <- rep(0, nboot)
  for(i in 1:nboot){
    b_exp.boot = stats::rnorm(length(b_exp), mean=b_exp, sd=se_exp)
    b_out.boot = stats::rnorm(length(b_out), mean=b_out, sd=se_out)
    betaIV.boot = b_out.boot/b_exp.boot
    med[i] = weighted_median(betaIV.boot, weights)
  }
  return(stats::sd(med))
}


mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  se <- ivw.res$coef["b_exp","Std. Error"]/min(1,ivw.res$sigma) #sigma is the residual standard error
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(data.frame(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}

mr_ivw_mre <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  se <- ivw.res$coef["b_exp","Std. Error"]
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(data.frame(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


mr_ivw_fe <- function(b_exp, b_out, se_exp, se_out, parameters=default_parameters())
{
  if(sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2)
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  
  ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, weights = 1/se_out^2))
  b <- ivw.res$coef["b_exp","Estimate"]
  se <- ivw.res$coef["b_exp","Std. Error"]/ivw.res$sigma
  pval <- 2 * stats::pnorm(abs(b/se), lower.tail=FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail=FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(data.frame(b = b, se = se, pval = pval, nsnp=length(b_exp), Q = Q, Q_df = Q_df, Q_pval = Q_pval))
}


#' Robust adjusted profile score
#'
#' @inheritParams mr_ivw
#' @param parameters A list of parameters. Specifically, `over.dispersion` and `loss.function`.
#' `over.dispersion` is a logical concerning should the model consider overdispersion (systematic pleiotropy).
#' And `loss.function` allows using either the squared error loss (`"l2"`) or robust loss functions/scores (`"huber"` or `"tukey"`).
#' The default is `parameters=list(overdispersion = TRUE, loss.function = "tukey")`.
#'
#' @details This function calls the \code{mr.raps} package. Please refer to the documentation of that package for more detail.
#'
#' @references Qingyuan Zhao, Jingshu Wang, Jack Bowden, Dylan S. Small. Statistical inference in two-sample summary-data Mendelian randomization using robust adjusted profile score. Forthcoming.
#'
#' @return List with the following elements:
#' \describe{
#' \item{b}{MR estimate}
#' \item{se}{Standard error}
#' \item{pval}{p-value}
#' \item{nsnp}{Number of SNPs}
#' }
#'
#' @export
mr_raps <- function(b_exp, b_out, se_exp, se_out, parameters = default_parameters()) {
  
  data <- data.frame(beta.exposure = b_exp,
                     beta.outcome = b_out,
                     se.exposure = se_exp,
                     se.outcome = se_out)
  out <- suppressMessages(
    mr.raps::mr.raps(data,
                     diagnostics = FALSE,
                     over.dispersion = parameters$over.dispersion,
                     loss.function = parameters$loss.function,
                     shrinkage = parameters$shrinkage))
  list(b = out$beta.hat,
       se = out$beta.se,
       pval = stats::pnorm(- abs(out$beta.hat / out$beta.se)) * 2,
       nsnp = length(b_exp))
  
}

#' MR sign test
#'
#' Tests how often the SNP-exposure and SNP-outcome signs are concordant.
#' This is to avoid the problem of averaging over all SNPs, which can suffer bias due to outliers with strong effects; and to avoid excluding SNPs which is implicit in median and mode based estimators.
#' The effect estimate here is not to be interpreted as the effect size - it is the proportion of SNP-exposure and SNP-outcome effects that have concordant signs.
#' e.g. +1 means all have the same sign, -1 means all have opposite signs, and 0 means that there is an equal number of concordant and discordant signs.
#' Restricted to only work if there are 6 or more valid SNPs.
#'
#' @param b_exp Vector of genetic effects on exposure
#' @param b_out Vector of genetic effects on outcome
#' @param se_exp Not required
#' @param se_out Not required
#' @param parameters Not required
#'
#' @export
#' @return List with the following elements:
#' \describe{
#' \item{b}{Concordance (see description)}
#' \item{se}{NA}
#' \item{pval}{p-value}
#' \item{nsnp}{Number of SNPs (excludes NAs and effect estimates that are 0)}
#' }
mr_sign <- function(b_exp, b_out, se_exp=NULL, se_out=NULL, parameters=NULL)
{
  b_exp[b_exp == 0] <- NA
  b_out[b_out == 0] <- NA
  if(sum(!is.na(b_exp) & !is.na(b_out)) < 6)
  {
    return(list(b=NA, se=NA, pval=NA, nsnp=NA))
  }
  x <- sum(sign(b_exp) == sign(b_out), na.rm=TRUE)
  n <- sum(!is.na(b_exp) & !is.na(b_out))
  
  out <- stats::binom.test(x=x, n=n, p=0.5)
  b <- (out$estimate - 0.5) * 2
  names(b) <- NULL
  pval <- out$p.value
  return(list(b=b, se=NA, pval=pval, nsnp=n))
}

default_parameters <- function()
{
  list(
    test_dist = "z",
    nboot = 1000,
    Cov = 0,
    penk = 20,
    phi = 1,
    alpha = 0.05,
    Qthresh = 0.05,
    over.dispersion = TRUE,
    loss.function = "huber",
    shrinkage = FALSE
  )
}



mr_leaveoneout <- function(dat, parameters=default_parameters(), method=mr_ivw)
{
  if(!"samplesize.outcome" %in% names(dat))
  {
    dat$samplesize.outcome <- NA
  }
  
  stopifnot("outcome" %in% names(dat))
  stopifnot("exposure" %in% names(dat))
  stopifnot("beta.exposure" %in% names(dat))
  stopifnot("beta.outcome" %in% names(dat))
  stopifnot("se.exposure" %in% names(dat))
  stopifnot("se.outcome" %in% names(dat))
  
  
  res <- plyr::ddply(dat, c("id.exposure", "id.outcome"), function(X)
  {
    x <- subset(X, mr_keep)
    nsnp <- nrow(x)
    if(nsnp == 0)
    {
      x <- X[1,]
      d <- data.frame(
        SNP = "All",
        b = NA,
        se = NA,
        p = NA,
        samplesize = NA,
        outcome = x$outcome[1],
        exposure = x$exposure[1]
      )
      return(d)
    }
    if(nsnp > 2)
    {
      l <- lapply(1:nsnp, function(i)
      {
        with(x, method(beta.exposure[-i], beta.outcome[-i], se.exposure[-i], se.outcome[-i], parameters))
      })
      l[[nsnp+1]] <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = c(as.character(x$SNP), "All"),
        b = sapply(l, function(y) y$b),
        se = sapply(l, function(y) y$se),
        p = sapply(l, function(y) y$pval),
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$outcome[1]
      d$exposure <- x$exposure[1]
      
    } else {
      a <- with(x, method(beta.exposure, beta.outcome, se.exposure, se.outcome, parameters))
      d <- data.frame(
        SNP = "All",
        b = a$b,
        se = a$se,
        p = a$pval,
        samplesize = x$samplesize.outcome[1]
      )
      d$outcome <- x$outcome[1]
      d$exposure <- x$exposure[1]
    }
    return(d)
  })
  res <- subset(res, select=c(exposure, outcome, id.exposure, id.outcome, samplesize, SNP, b, se, p))
  return(res)
}



#' Plot results from leaveoneout analysis
#' 
#' Plot results from leaveoneout analysis.
#'
#' @param leaveoneout_results Output from [mr_leaveoneout()].
#'
#' @export
#' @return List of plots
mr_leaveoneout_plot <- function(leaveoneout_results)
{
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    # Need to have at least 3 SNPs because IVW etc methods can't be performed with fewer than 2 SNPs
    if(sum(!grepl("All", d$SNP)) < 3) {
      return(
        blank_plot("Insufficient number of SNPs")
      )
    }
    d$up <- d$b + 1.96 * d$se
    d$lo <- d$b - 1.96 * d$se
    d$tot <- 1
    d$tot[d$SNP != "All"] <- 0.01
    d$SNP <- as.character(d$SNP)
    nom <- d$SNP[d$SNP != "All"]
    nom <- nom[order(d$b)]
    d <- rbind(d, d[nrow(d),])
    d$SNP[nrow(d)-1] <- ""
    d$b[nrow(d)-1] <- NA
    d$up[nrow(d)-1] <- NA
    d$lo[nrow(d)-1] <- NA
    d$SNP <- ordered(d$SNP, levels=c("All", "", nom))
    
    ggplot2::ggplot(d, ggplot2::aes(y=SNP, x=b)) +
      ggplot2::geom_vline(xintercept=0, linetype="dotted") +
      # ggplot2::geom_errorbarh(ggplot2::aes(xmin=pmax(lo, min(d$b, na.rm=T)), xmax=pmin(up, max(d$b, na.rm=T)), size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=up, size=as.factor(tot), colour=as.factor(tot)), height=0) +
      ggplot2::geom_point(ggplot2::aes(colour=as.factor(tot))) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% "")), colour="grey") +
      ggplot2::scale_colour_manual(values=c("black", "red")) +
      ggplot2::scale_size_manual(values=c(0.3, 1)) +
      # xlim(c(min(c(0, d$b), na.rm=T), max(c(0, d$b), na.rm=T))) +
      ggplot2::theme(
        legend.position="none", 
        axis.text.y=ggplot2::element_text(size=6), 
        axis.ticks.y=ggplot2::element_line(size=0),
        axis.title.x=ggplot2::element_text(size=8)) +
      ggplot2::labs(y="", x=paste0("MR leave-one-out sensitivity analysis for\n'", d$exposure[1], "' on '", d$outcome[1], "'"))+
      ggplot2::theme_bw()
  })
  res
}


#' Create scatter plot with lines showing the causal estimate for different MR tests
#'
#' Requires dev version of ggplot2
#' 
#' @param mr_results Output from [mr()].
#' @param dat Output from [harmonise_data()].
#' @export
#' @return List of plots
mr_scatter_plot <- function(mr_results, dat)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position="top", legend.direction="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
  })
  mrres
}


blank_plot <- function(message)
{
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + 
    ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + 
    ggplot2::labs(x=NULL,y=NULL) + 
    ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}
