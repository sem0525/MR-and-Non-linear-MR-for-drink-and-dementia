

summ_data$summary
by = summ_data$summary$by
bx = summ_data$summary$bx
byse = summ_data$summary$byse
bxse = summ_data$summary$bxse
xmean = summ_data$summary$xmean
frac_poly_summ_mr <- function(by, bx, byse, bxse, xmean, method = "FE", d = "both",
                              powers = c(0, -2, -1.5, -1, -0.5, 1, 2),
                              pd = 0.05, average.exposure.associations = FALSE, ci = "model_se", nboot = 100,
                              fig = FALSE, family = "binomial", offset = 0,
                              pref_x = "x", pref_y = "y", ref = NA,
                              ci_type = "overall", breaks = NULL,
                              ylim_lower = NA, ylim_upper = NA,
                              xlim_lower = NA, xlim_upper = NA, seed=335) {
  
  #### Start ###
  frac_coef <- by
  frac_se <- byse
  xcoef_sub <- bx
  xcoef_sub_se <- bxse
  ##### Best-fitting fractional polynomial of degree 1 #####
  p <- NULL
  ML <- NULL
  j <- 1
  for (p1 in powers) {
    if (p1 == -1) {
      x1 <- xmean^p1
    } else {
      x1 <- (p1 + 1) * xmean^p1
    }
    p[j] <- p1
    cc <- try(rma(frac_coef / xcoef ~ -1 + x1,
                  vi = (frac_se / xcoef)^2,
                  method = method
    ), silent = TRUE)
    if (methods::is(cc, "try-error") == T) {
      ML[j] <- NA
    }
    if (methods::is(cc, "try-error") == F) {
      ML[j] <- summary(rma(frac_coef / xcoef ~ -1 + x1,
                           vi = (frac_se / xcoef)^2,
                           method = method
      ))$fit.stats[1, 1]
    }
    
    j <- j + 1
  }
  fits <- data.frame(p, ML)
  
  # Debugging print statements
  print("Fits data frame:")
  print(fits)
  
  fits$max <- 0
  fits$max[fits$ML == max(fits$ML, na.rm = T)] <- 1
  p_ML <- fits$p[fits$max == 1]
  
  print("p_ML:")
  print(p_ML)
  
  ##### Best-fitting fractional polynomial of degree 2 #####
  if (d == 1 | d == 2 | d == "both") {
    powers1 <- powers
    powers2 <- powers
    p1 <- NULL
    p2 <- NULL
    ML <- NULL
    j <- 1
    for (p11 in powers1) {
      if (p11 == -1) {
        x1 <- xmean^p11
      } else {
        x1 <- (p11 + 1) * xmean^p11
      }
      for (p21 in powers2) {
        if (p11 == p21) {
          if (p21 == -1) {
            x2 <- 2 * (xmean^p21) * log(xmean)
          } else {
            x2 <- ((p21 + 1) * (xmean^p21) * log(xmean) + xmean^p21)
          }
        } else {
          if (p21 == -1) {
            x2 <- xmean^p21
          } else {
            x2 <- (p21 + 1) * xmean^p21
          }
        }
        p1[j] <- p
        
        
        
  
  if (average.exposure.associations == TRUE) { xcoef <- sum(bx * (bxse^(-2))) / sum(bxse^(-2)) }
  else { xcoef <- bx }
  q <- length(by)
  ##### Best-fitting fractional polynomial of degree 1 #####
  p <- NULL
  ML <- NULL
  j <- 1
  for (p1 in powers) {
    if (p1 == -1) {
      x1 <- xmean^p1
    } else {
      x1 <- (p1 + 1) * xmean^p1
    }
    p[j] <- p1
    rma(frac_coef / xcoef ~ -1 + x1,
                  vi = (frac_se / xcoef)^2,
                  method = method
    )
 
      ML[j] <- summary(rma(frac_coef / xcoef ~ -1 + x1,
                           vi = (frac_se / xcoef)^2,
                           method = method
      ))$fit.stats[1, 1]
    }
    j <- j + 1
  }
  fits <- data.frame(p, ML)
  fits$max <- 0
  fits$max[fits$ML == max(fits$ML, na.rm = T)] <- 1
  p_ML <- fits$p[fits$max == 1]
  
  ##### Best-fitting fractional polynomial of degree 2 #####
  if (d == 1 | d == 2 | d == "both") {
    powers1 <- powers
    powers2 <- powers
    p1 <- NULL
    p2 <- NULL
    ML <- NULL
    j <- 1
    for (p11 in powers1) {
      if (p11 == -1) {
        x1 <- xmean^p11
      } else {
        x1 <- (p11 + 1) * xmean^p11
      }
      for (p21 in powers2) {
        if (p11 == p21) {
          if (p21 == -1) {
            x2 <- 2 * (xmean^p21) * log(xmean)
          } else {
            x2 <- ((p21 + 1) * (xmean^p21) * log(xmean) + xmean^p21)
          }
        } else {
          if (p21 == -1) {
            x2 <- xmean^p21
          } else {
            x2 <- (p21 + 1) * xmean^p21
          }
        }
        p1[j] <- p11
        p2[j] <- p21
        cc <- try(rma(frac_coef / xcoef ~ -1 + x1 + x2,
                      vi = (frac_se / xcoef)^2,
                      method = method
        ), silent = TRUE)
        if (methods::is(cc, "try-error") == T) {
          ML[j] <- NA
        }
        if (methods::is(cc, "try-error") == F) {
          ML[j] <- summary(rma(frac_coef / xcoef ~ -1 + x1 + x2,
                               vi = (frac_se / xcoef)^2,
                               method = method
          ))$fit.stats[1, 1]
        }
        j <- j + 1
      }
      powers2 <- powers2[-1]
    }
    fits <- data.frame(p1, p2, ML)
    fits$max <- 0
    fits$max[fits$ML == max(fits$ML, na.rm = T)] <- 1
    p1_ML <- fits$p1[fits$max == 1]
    p2_ML <- fits$p2[fits$max == 1]
  }
  
  ##### Best-fitting fractional polynomial of either degree 1 or degree 2 ###
  p_d1_d2 <- NA
  if (d == 1 | d == 2 | d == "both") {
    if (p_ML == -1) {
      x1 <- xmean^p_ML
    } else {
      x1 <- (p_ML + 1) * xmean^p_ML
    }
    best_fracp_d1 <- rma(frac_coef / xcoef ~ -1 + x1,
                         vi = (frac_se / xcoef)^2,
                         method = method
    )
    dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2, 1]
    if (p1_ML == -1) {
      x1 <- xmean^p1_ML
    } else {
      x1 <- (p1_ML + 1) * xmean^p1_ML
    }
    if (p1_ML == p2_ML) {
      if (p2_ML == -1) {
        x2 <- 2 * (xmean^p2_ML) * log(xmean)
      } else {
        x2 <- ((p2_ML + 1) * (xmean^p2_ML) * log(xmean) + xmean^p2_ML)
      }
    } else {
      if (p2_ML == -1) {
        x2 <- xmean^p2_ML
      } else {
        x2 <- (p2_ML + 1) * xmean^p2_ML
      }
    }
    best_fracp_d2 <- rma(frac_coef / xcoef ~ -1 + x1 + x2,
                         vi = (frac_se / xcoef)^2,
                         method = method
    )
    dev_best_fracp_d2 <- best_fracp_d2$fit.stats[2, 1]
    p_d1_d2 <- 1 - stats::pchisq((dev_best_fracp_d1 - dev_best_fracp_d2),
                                 df = 2
    )
    if (p_d1_d2 >= pd) {
      d1 <- 1
    } else {
      d1 <- 2
    }
    if (d == "both") {
      d <- d1
    }
  }
  
  ##### Model #####
  if (d == 1) {
    if (p_ML == -1) {
      x1 <- xmean^p_ML
    } else {
      x1 <- (p_ML + 1) * xmean^p_ML
    }
    model <- rma(frac_coef / xcoef ~ -1 + x1,
                 vi = (frac_se / xcoef)^2,
                 method = method
    )
  }
  if (d == 2) {
    if (p1_ML == -1) {
      x1 <- xmean^p1_ML
    } else {
      x1 <- (p1_ML + 1) * xmean^p1_ML
    }
    if (p1_ML == p2_ML) {
      if (p2_ML == -1) {
        x2 <- 2 * (xmean^p2_ML) * log(xmean)
      } else {
        x2 <- ((p2_ML + 1) * (xmean^p2_ML) * log(xmean) + xmean^p2_ML)
      }
    } else {
      if (p2_ML == -1) {
        x2 <- xmean^p2_ML
      } else {
        x2 <- (p2_ML + 1) * xmean^p2_ML
      }
    }
    model <- rma(frac_coef / xcoef ~ -1 + x1 + x2,
                 vi = (frac_se / xcoef)^2,
                 method = method
    )
  }
  
  ##### Bootstrap #####
  
  if (ci == "bootstrap_per" | ci == "bootstrap_se") {
    if (d == 1) {
      frac_coef_boot <- NULL
      for (i in 1:nboot) {
        frac_coef1 <- by + stats::rnorm(q, 0, byse)
        frac_se1 <- byse
        xmean1 <- xmean
        if (p_ML == -1) {
          x111 <- xmean1^p_ML
        } else {
          x111 <- (p_ML + 1) * xmean1^p_ML
        }
        mod <- rma.uni(frac_coef1 / xcoef ~ -1 + x111,
                       vi = (frac_se1 / xcoef)^2, method = method
        )
        frac_coef_boot[i] <- mod$b[1]
      }
    }
    if (d == 2) {
      frac_coef_boot <- matrix(, nrow = nboot, ncol = 2)
      for (i in 1:nboot) {
        frac_coef1 <- by + stats::rnorm(q, 0, byse)
        frac_se1 <- byse
        xmean1 <- xmean
        if (p1_ML == -1) {
          x111 <- xmean1^p1_ML
        } else {
          x111 <- (p1_ML + 1) * xmean1^p1_ML
        }
        if (p1_ML == p2_ML) {
          if (p2_ML == -1) {
            x211 <- 2 * (xmean1^p2_ML) * log(xmean1)
          } else {
            x211 <- ((p2_ML + 1) * (xmean1^p2_ML) * log(xmean1)
                     + xmean1^p2_ML)
          }
        } else {
          if (p2_ML == -1) {
            x211 <- xmean1^p2_ML
          } else {
            x211 <- (p2_ML + 1) * xmean1^p2_ML
          }
        }
        mod <- rma.uni(frac_coef1 / xcoef ~ -1 + x111 + x211,
                       vi = (frac_se1 / xcoef)^2, method = method
        )
        frac_coef_boot[i, 1] <- mod$b[1]
        frac_coef_boot[i, 2] <- mod$b[2]
      }
    }
  }
  
  ##### Fractional polynomial degree 1 test against linearity #####
  if (p_ML == -1) {
    x1 <- xmean^p_ML
  } else {
    x1 <- (p_ML + 1) * xmean^p_ML
  }
  linear <- rma(frac_coef / xcoef ~ 1, vi = (frac_se / xcoef)^2, method = method)
  dev_linear <- linear$fit.stats[2, 1]
  best_fracp_d1 <- rma(frac_coef / xcoef ~ -1 + x1,
                       vi = (frac_se / xcoef)^2,
                       method = method
  )
  dev_best_fracp_d1 <- best_fracp_d1$fit.stats[2, 1]
  p_fp <- 1 - stats::pchisq((dev_linear - dev_best_fracp_d1), df = 1)
  
  ##### Other tests #####
  p_quadratic <- rma(frac_coef / xcoef ~ xmean,
                     vi = (frac_se / xcoef)^2,
                     method = method
  )$pval[2]
  p_q <- 1 - stats::pchisq(rma(frac_coef / xcoef, vi = (frac_se / xcoef)^2)$QE,
                           df = (q - 1)
  )
  
  ##### Results #####
  beta <- as.numeric(model$b)
  if (ci == "model_se") {
    if (d == 1) {
      powers <- p_ML + 1
    }
    if (d == 2) {
      powers <- c(p1_ML, p2_ML)
      powers <- powers + 1
    }
    cov <- model$vb
    se <- model$se
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * stats::pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_se") {
    if (d == 1) {
      powers <- p_ML + 1
      cov <- stats::var(frac_coef_boot)
      se <- sqrt(cov)
    }
    if (d == 2) {
      powers <- c(p1_ML, p2_ML)
      powers <- powers + 1
      cov <- cov(frac_coef_boot)
      se <- sqrt(diag(cov))
    }
    lci <- beta - 1.96 * se
    uci <- beta + 1.96 * se
    pval <- 2 * stats::pnorm(-abs(beta / se))
  }
  if (ci == "bootstrap_per") {
    if (d == 1) {
      powers <- p_ML + 1
      se <- NA
      lci <- stats::quantile(frac_coef_boot, probs = 0.025)
      uci <- stats::quantile(frac_coef_boot, probs = 0.975)
      pval <- NA
    }
    if (d == 2) {
      powers <- c(p1_ML, p2_ML)
      powers <- powers + 1
      se <- rep(NA, 2)
      lci <- NULL
      uci <- NULL
      pval <- NULL
      lci[1] <- stats::quantile(frac_coef_boot[, 1], probs = 0.025)
      lci[2] <- stats::quantile(frac_coef_boot[, 2], probs = 0.025)
      uci[1] <- stats::quantile(frac_coef_boot[, 1], probs = 0.975)
      uci[2] <- stats::quantile(frac_coef_boot[, 2], probs = 0.975)
      pval <- rep(NA, 2)
    }
  }
  lci <- as.numeric(lci)
  uci <- as.numeric(uci)
  if (ci == "model_se") {
    nboot <- NA
  }
  
  ##### Figure #####
  if (fig == TRUE) {
    if (ci_type == "overall") {
      if (is.na(ref)) {
        ref <- mean(xmean)
      }
      plot_data <- data.frame(x = stats::runif(
        10000, min(xmean),
        max(xmean)
      ))
      plot_data_1 <- data.frame(x = ref, y = 0)
      if (d == 1) {
        if (p_ML == -1) {
          plot_data$yest <- beta * log(plot_data$x) - (beta * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) - log(ref))^2 *
                                    as.vector(cov))
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot) -
              reprow(log(ref) %*% t(frac_coef_boot),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        } else {
          plot_data$yest <- beta * plot_data$x^(p_ML + 1) -
            beta * ref^(p_ML + 1)
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p_ML + 1) -
                                     ref^(p_ML + 1))^2 * as.vector(cov))
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p_ML + 1) %*% t(frac_coef_boot) -
              reprow(ref^(p_ML + 1) %*% t(frac_coef_boot),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
      }
      if (d == 2) {
        if (p1_ML == -1 & p2_ML == -1) {
          plot_data$yest <- beta[1] * log(plot_data$x) +
            beta[2] * log(plot_data$x) * log(plot_data$x) -
            (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) -
                                     log(ref))^2 * cov[1, 1] +
                                    2 * (log(plot_data$x) -
                                           log(ref)) * (log(plot_data$x) *
                                                          log(plot_data$x) -
                                                          log(ref) * log(ref)) * cov[1, 2] +
                                    (log(plot_data$x) * log(plot_data$x) -
                                       log(ref) * log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
              log(plot_data$x) *
              log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
              reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                       log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML == -1 & p2_ML != -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * log(plot_data$x) +
            beta[2] * plot_data$x^(p2_ML + 1) - (beta[1] * log(ref) +
                                                   beta[2] * ref^(p2_ML + 1))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) -
                                     log(ref))^2 * cov[1, 1] +
                                    2 * (log(plot_data$x) -
                                           log(ref)) * (plot_data$x^(p2_ML + 1) -
                                                          ref^(p2_ML + 1)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
              plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
              reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                       ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML == -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * log(plot_data$x) -
            (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (log(plot_data$x) -
                                                                 log(ref)) * cov[1, 2] +
                                    (log(plot_data$x) -
                                       log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       log(ref) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML != -1 & p1_ML == p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * plot_data$x^(p2_ML + 1) * log(plot_data$x) -
            (beta[1] * ref^(p1_ML + 1) +
               beta[2] * ref^(p2_ML + 1) * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) *
                                                                 log(plot_data$x) -
                                                                 ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                                       ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              (plot_data$x^(p2_ML + 1) *
                 log(plot_data$x)) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML != -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * plot_data$x^(p2_ML + 1) -
            (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) -
                                                                 ref^(p2_ML + 1)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
      }
      plot_data$x <- plot_data$x + offset
      plot_data_1$x <- plot_data_1$x + offset
      ref <- ref + offset
      if (family != "binomial") {
        figure <- ggplot2::ggplot(plot_data, aes(x = x))
        figure <- figure +
          geom_hline(aes(yintercept = 0),linetype="dashed", color = "black") +
          
          geom_line(aes(y = lci), color = "#56B4E9") +
          geom_line(aes(y = uci), color = "#56B4E9") +
          geom_ribbon(aes(x,ymin = lci, ymax = uci),
                      alpha = 1,fill="#56B4E9") +
          geom_line(aes(y = yest), color = "black") +
          theme_bw() +
          labs(x = pref_x, y = pref_y) +
          theme(
            axis.title.x = element_text(vjust = 0.5, size = 20),
            axis.title.y = element_text(vjust = 0.5, size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18)
          ) +
          geom_point(aes(x = x, y = y),
                     data = plot_data_1,
                     colour = "grey",
                     size = 1
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        if (!is.null(breaks)) {
          suppressMessages(figure <- figure +
                             scale_y_continuous(breaks = breaks))
        }
      }
      if (family == "binomial") {
        plot_data$yest <- exp(plot_data$yest)
        plot_data$uci <- exp(plot_data$uci)
        plot_data$lci <- exp(plot_data$lci)
        plot_data_1$y <- exp(0)
        figure <- ggplot2::ggplot(plot_data, aes(x = x))
        figure <- figure +
          geom_hline(aes(yintercept = 1), linetype="dashed", color = "black") +
          
          geom_line(aes(y = lci), color = "#56B4E9") +
          geom_line(aes(y = uci), color = "#56B4E9") +
          geom_ribbon(aes(x,ymin = lci, ymax = uci),
                      alpha = 1,fill="#56B4E9") +
          geom_line(aes(y = yest), color = "black") +
          theme_bw() +
          labs(x = pref_x, y = pref_y) +
          theme(
            axis.title.x = element_text(vjust = 0.5, size = 16),
            axis.title.y = element_text(vjust = 0.5, size = 16),
            axis.text.x = element_text(size = 14),
            axis.text.y = element_text(size = 14)
          ) +
          geom_point(aes(x = x, y = y),
                     data = plot_data_1,
                     colour = "red",
                     size = 4
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        if (!is.null(breaks)) {
          figure <- figure + scale_y_continuous(breaks = breaks)
        }
        figure <- figure + coord_trans(y = "log")
      }
    }
    if (ci_type == "quantile") {
      if (is.na(ref)) {
        ref <- mean(xmean)
      }
      xmean_ci <- xmean
      plot_data <- data.frame(x = c(ref, xmean_ci))
      plot_data_1 <- data.frame(x = stats::runif(
        10000, min(xmean),
        max(xmean)
      ))
      if (d == 1) {
        if (p_ML == -1) {
          plot_data$yest <- beta * log(plot_data$x) -
            (beta * log(ref))
          plot_data_1$yest <- beta * log(plot_data_1$x) -
            (beta * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) - log(ref))^2 *
                                    as.vector(cov))
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot) -
              reprow(log(ref) %*% t(frac_coef_boot),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p_ML != -1) {
          plot_data$yest <- beta * plot_data$x^(p_ML + 1) -
            beta * ref^(p_ML + 1)
          plot_data_1$yest <- beta * plot_data_1$x^(p_ML + 1) -
            beta * ref^(p_ML + 1)
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p_ML + 1) -
                                     ref^(p_ML + 1))^2 * as.vector(cov))
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p_ML + 1) %*% t(frac_coef_boot) -
              reprow(ref^(p_ML + 1) %*% t(frac_coef_boot),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
      }
      if (d == 2) {
        if (p1_ML == -1 & p2_ML == -1) {
          plot_data$yest <- beta[1] * log(plot_data$x) +
            beta[2] * log(plot_data$x) * log(plot_data$x) -
            (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
          plot_data_1$yest <- beta[1] * log(plot_data_1$x) +
            beta[2] * log(plot_data_1$x) * log(plot_data_1$x) -
            (beta[1] * log(ref) + beta[2] * log(ref) * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) -
                                     log(ref))^2 * cov[1, 1] +
                                    2 * (log(plot_data$x) -
                                           log(ref)) * (log(plot_data$x) * log(plot_data$x) -
                                                          log(ref) * log(ref)) * cov[1, 2] +
                                    (log(plot_data$x) * log(plot_data$x) -
                                       log(ref) * log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
              log(plot_data$x) *
              log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
              reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                       log(ref) * log(ref) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML == -1 & p2_ML != -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * log(plot_data$x) +
            beta[2] * plot_data$x^(p2_ML + 1) -
            (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
          plot_data_1$yest <- beta[1] * log(plot_data_1$x) +
            beta[2] * plot_data_1$x^(p2_ML + 1) -
            (beta[1] * log(ref) + beta[2] * ref^(p2_ML + 1))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((log(plot_data$x) -
                                     log(ref))^2 * cov[1, 1] +
                                    2 * (log(plot_data$x) - log(ref)) *
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- log(plot_data$x) %*% t(frac_coef_boot[, 1]) +
              plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
              reprow(log(ref) %*% t(frac_coef_boot[, 1]) +
                       ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML == -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * log(plot_data$x) -
            (beta[1] * ref^(p1_ML + 1) + beta[2] * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (log(plot_data) -
                                                                 log(ref)) * cov[1, 2] +
                                    (log(plot_data$x) -
                                       log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              log(plot_data$x) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       log(ref) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML != -1 & p1_ML == p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * plot_data$x^(p2_ML + 1) * log(plot_data$x) -
            (beta[1] * ref^(p1_ML + 1) +
               beta[2] * ref^(p2_ML + 1) * log(ref))
          plot_data_1$yest <- beta[1] * plot_data_1$x^(p1_ML + 1) +
            beta[2] * plot_data_1$x^(p2_ML + 1) * log(plot_data_1$x) -
            (beta[1] * ref^(p1_ML + 1) +
               beta[2] * ref^(p2_ML + 1) * log(ref))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) *
                                                                 log(plot_data$x) -
                                                                 ref^(p2_ML + 1) * log(ref)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) * log(plot_data$x) -
                                       ref^(p2_ML + 1) * log(ref))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              (plot_data$x^(p2_ML + 1) *
                 log(plot_data$x)) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       (ref^(p2_ML + 1) * log(ref)) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
        if (p1_ML != -1 & p2_ML != -1 & p1_ML != p2_ML) {
          plot_data$yest <- beta[1] * plot_data$x^(p1_ML + 1) +
            beta[2] * plot_data$x^(p2_ML + 1) -
            (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
          plot_data_1$yest <- beta[1] * plot_data_1$x^(p1_ML + 1) +
            beta[2] * plot_data_1$x^(p2_ML + 1) -
            (beta[1] * ref^(p1_ML + 1) + beta[2] * ref^(p2_ML + 1))
          if (ci != "bootstrap_per") {
            plot_data$yse <- sqrt((plot_data$x^(p1_ML + 1) -
                                     ref^(p1_ML + 1))^2 * cov[1, 1] +
                                    2 * (plot_data$x^(p1_ML + 1) -
                                           ref^(p1_ML + 1)) * (plot_data$x^(p2_ML + 1) -
                                                                 ref^(p2_ML + 1)) * cov[1, 2] +
                                    (plot_data$x^(p2_ML + 1) -
                                       ref^(p2_ML + 1))^2 * cov[2, 2])
            plot_data$lci <- plot_data$yest - 1.96 * plot_data$yse
            plot_data$uci <- plot_data$yest + 1.96 * plot_data$yse
          } else {
            boot <- plot_data$x^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
              plot_data$x^(p2_ML + 1) %*% t(frac_coef_boot[, 2]) -
              reprow(ref^(p1_ML + 1) %*% t(frac_coef_boot[, 1]) +
                       ref^(p2_ML + 1) %*% t(frac_coef_boot[, 2]),
                     n = nrow(plot_data)
              )
            plot_data$lci <- rowQuantiles(boot, probs = 0.025)
            plot_data$uci <- rowQuantiles(boot, probs = 0.975)
          }
        }
      }
      highlight <- c("red", rep("black", (nrow(plot_data) - 1)))
      plot_data$x <- plot_data$x + offset
      plot_data_1$x <- plot_data_1$x + offset
      if (family != "binomial") {
        figure <- ggplot2::ggplot(plot_data, aes(x = x))
        figure <- figure +
          geom_hline(aes(yintercept = 0), linetype="dashed", color = "black") +
          
          geom_errorbar(
            mapping = aes(x = x, ymin = lci, ymax = uci),
            color = "#56B4E9", width = 0.025
          ) +
          geom_ribbon(aes(x,ymin = lci, ymax = uci),
                      alpha = 1,fill="#56B4E9") +
          geom_line(aes(x = x, y = yest),
                    color = "black",
                    data = plot_data_1
          ) +
          geom_point(aes(y = yest), color = highlight, size = 4) +
          theme_bw() +
          labs(x = pref_x, y = pref_y) +
          theme(
            axis.title.x = element_text(vjust = 0.5, size = 20),
            axis.title.y = element_text(vjust = 0.5, size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18)
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        if (!is.null(breaks)) {
          figure <- figure + scale_y_continuous(breaks = breaks)
        }
      }
      if (family == "binomial") {
        plot_data$yest <- exp(plot_data$yest)
        plot_data$uci <- exp(plot_data$uci)
        plot_data$lci <- exp(plot_data$lci)
        figure <- ggplot2::ggplot(plot_data, aes(x = x))
        plot_data_1$yest <- exp(plot_data_1$yest)
        figure <- figure +
          geom_hline(aes(yintercept = 1), linetype="dashed", color = "black") +
          geom_line(aes(x = x, y = yest),
                    color = "black",
                    data = plot_data_1
          ) +
          geom_errorbar(
            mapping = aes(x = x, ymin = lci, ymax = uci),
            color = "#56B4E9", width = 0.025
          ) +
          geom_ribbon(aes(x,ymin = lci, ymax = uci),
                      alpha = 1,fill="#56B4E9") +
          geom_point(aes(y = yest), color = highlight, size = 4) +
          theme_bw() +
          labs(x = pref_x, y = pref_y) +
          theme(
            axis.title.x = element_text(vjust = 0.5, size = 20),
            axis.title.y = element_text(vjust = 0.5, size = 20),
            axis.text.x = element_text(size = 18),
            axis.text.y = element_text(size = 18)
          ) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()
          )
        if (!is.null(breaks)) {
          figure <- figure + scale_y_continuous(breaks = breaks)
        }
        figure <- figure + coord_trans(y = "log")
      }
    }
    figure <- figure + theme(
      panel.border = element_blank(),
      axis.line = element_line()
    )
    if (!is.na(ylim_lower) | !is.na(ylim_upper)) {
      figure <- figure +
        scale_y_continuous(
          limits = c(ylim_lower, ylim_upper),
          breaks = breaks
        )
    }
    if (!is.na(xlim_lower) | !is.na(xlim_upper)) {
      figure <- figure + xlim(xlim_lower, xlim_upper)
    }
  }
  
  
  ##### Return #####
  model <- as.matrix(data.frame(
    q = q, xpos = "user",
    ci_type = ci, nboot = nboot
  ))
  coefficients <- as.matrix(data.frame(
    beta = beta, se = se,
    lci = lci, uci = uci, pval = pval
  ))
  rownames(coefficients) <- powers
  lace <- as.matrix(data.frame(
    beta = (frac_coef / xcoef),
    se = (abs(frac_se / xcoef)),
    lci = (frac_coef / xcoef - 1.96 * (abs(frac_se / xcoef))),
    uci = (frac_coef / xcoef + 1.96 * (abs(frac_se / xcoef))),
    pval = (2 * stats::pnorm(-abs(frac_coef / frac_se)))
  ))
  rownames(lace) <- 1:nrow(lace)
  xcoef_quant <- as.matrix(data.frame(beta = xcoef_sub, se = xcoef_sub_se))
  rownames(xcoef_quant) <- 1:nrow(xcoef_quant)
  p_tests <- as.matrix(data.frame(
    fp_d1_d2 = p_d1_d2, fp = p_fp,
    quad = p_quadratic, Q = p_q
  ))
  
  
  if (fig == TRUE) {
    results <- list(
      n = NA, model = model, powers = powers,
      coefficients = coefficients, lace = lace,
      xcoef = xcoef_quant, p_tests = p_tests, figure = figure
    )
  }
  if (fig == FALSE) {
    results <- list(
      n = NA, model = model, powers = powers,
      coefficients = coefficients, lace = lace,
      xcoef = xcoef_quant, p_tests = p_tests
    )
  }
  class(results) <- "frac_poly_mr"
  return(results)
}

    