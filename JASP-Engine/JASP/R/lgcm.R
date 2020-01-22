#
# Copyright (C) 2013-2019 University of Amsterdam
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

LatentGrowthCurve <- function(jaspResults, dataset, options, ...) {
  ready <- length(options[["variables"]]) > 2
  
  # Read dataset
  dataset <- .lgcmReadData(dataset, options)
  
  # Preprocess options
  options <- .lgcmPreprocessOptions(dataset, options)
  
  # Enrich dataset
  dataset <- .lgcmEnrichData(dataset, options)
  
  # Error checking
  errors <- .lgcmCheckErrors(dataset, options)
  
  # Create model container
  modelContainer <- .lgcmModelContainer(jaspResults, options)
  
  # output
  .lgcmFitTable(            modelContainer, dataset, options, ready )
  .lgcmParameterTables(     modelContainer, dataset, options, ready )
  .lcgmAdditionalFitTables( modelContainer, dataset, options, ready )
  .lgcmCurvePlot(           modelContainer, dataset, options, ready )
  .lgcmPathPlot(            modelContainer, dataset, options, ready )
  .lgcmSyntax(              modelContainer, dataset, options, ready )
  #.lgcmMisfitPlot(      modelContainer, dataset, options, ready )
  
}

# Preprocessing functions ----
.lgcmPreprocessOptions <- function(dataset, options) {
  # add dummy names
  if (length(options[["categorical"]]) > 0) {
    frml <- as.formula(paste("~", paste(.v(options[["categorical"]]), collapse = "+")))
    dumnames <- colnames(model.matrix(frml, dataset))[-1]
    options[["dummy_names"]] <- stringr::str_replace_all(
      string      = dumnames, 
      pattern     = .v(options[["categorical"]]), 
      replacement = options[["categorical"]]
    )
  }
  return(options)
}

.lgcmReadData <- function(dataset, options) {
  if (!is.null(dataset)) return(dataset)
  .readDataSetToEnd(columns = c(
    options[["variables"]], 
    options[["regressions"]], 
    options[["categorical"]], 
    options[["covariates"]]
  ))
}

.lgcmEnrichData <- function(dataset, options) {
  if (interactive()) browser()
  # Add dummies
  if (length(options[["categorical"]]) > 0) {
    frml <- as.formula(paste("~", paste(.v(options[["categorical"]]), collapse = "+")))
    mm <- model.matrix(frml, dataset)[,-1, drop = FALSE]
    colnames(mm) <- .v(options[["dummy_names"]])
    dataset <- cbind(dataset, mm)
    
  }
  return(dataset)
}

.lgcmCheckErrors <- function(dataset, options) {
  # some error check
  return(TRUE)
}

# Results functions ----
.lgcmComputeResults <- function(modelContainer, dataset, options) {
  lgcmResult <- try(lavaan::growth(
    model     = .lgcmOptionsToMod(options),
    data      = dataset,
    se        = ifelse(options[["se"]] == "bootstrap", "standard", options[["se"]]),
    mimic     = options[["mimic"]],
    estimator = options[["estimator"]],
    missing   = options[["missing"]]
  ))
  if (inherits(lgcmResult, "try-error")) {
    modelContainer$setError(paste(
      "Model error:", 
      .decodeVarsInMessage(names(dataset), attr(lgcmResult, "condition")$message))
    )
    return()
  }
  
  admissible <- .withWarnings(lavaan:::lav_object_post_check(lgcmResult))
  
  # if (!admissible$value) {
  #   modelContainer$setError(paste(
  #     "The model is not admissible:", 
  #     .decodeVarsInMessage(names(dataset), admissible$warnings[[1]]$message))
  #   )
  # }
  
  if (!lgcmResult@optim$converged) {
    modelContainer$setError("The model could not be estimated due to nonconvergence.")
  }
  
  if (lgcmResult@test[[1]]$df < 0) {
    modelContainer$setError("The model could not be estimated: No degrees of freedom left.")
  }
  
  
  # Bootstrapping with interruptible progress bar
  if (options[["se"]] == "bootstrap") {
    startProgressbar(options[["bootstrapNumber"]])
    
    boot_1      <- lavaan::bootstrapLavaan(lgcmResult, R = 1)
    bootres     <- matrix(0, options[["bootstrapNumber"]], length(boot_1))
    bootres[1,] <- boot_1
    for (i in 2:options[["bootstrapNumber"]]) {
      bootres[i,] <- lavaan::bootstrapLavaan(lgcmResult, 1)
      progressbarTick()
    }
    
    lgcmResult@boot       <- list(coef = bootres)
    lgcmResult@options[["se"]] <- "bootstrap"
  }
  
  # Save cfaResult as state so it's available even when opts don't change
  modelContainer[["model"]] <- createJaspState(lgcmResult)
  return(lgcmResult)
}

.lgcmOptionsToMod <- function(options, base64 = TRUE) {
  if (!base64) .v <- I
  timings <- sapply(options[["timings"]], function(t) t$timing)
  
  # Header info
  Hed <- paste0(
    "# -------------------------------------------\n",
    "# Latent Growth Curve model generated by JASP\n", 
    "# -------------------------------------------\n"
  )
  
  # Basic LGCM curve information
  Int <- if (options[["intercept"]])
    paste("I =~", paste0("1*", .v(options[["variables"]]), collapse = " + "))
  else NULL
  Lin <- if (options[["linear"]]) 
    paste("\nL =~", paste0(timings, "*", .v(options[["variables"]]), collapse = " + "))
  else NULL
  Qua <- if (options[["quadratic"]])
    paste("\nQ =~", paste0(timings^2, "*", .v(options[["variables"]]), collapse = " + "))
  else NULL
  Cub <- if (options[["cubic"]])
    paste("\nC =~", paste0(timings^3, "*", .v(options[["variables"]]), collapse = " + "))
  else NULL
  LGC <- paste0("\n# Growth curve\n", Int, Lin, Qua, Cub)
  
  curve <- c("I", "L", "Q", "C")[c(options[["intercept"]], options[["linear"]], options[["quadratic"]], options[["cubic"]])]
  
  # Covarying latents
  if (!options[["covar"]]) {
    Cov <- "\n\n# Suppress latent covariance"
    for (i in seq_along(curve))
      for (j in seq_along(curve))
        if (i < j) Cov <- paste0(Cov, "\n", curve[i], " ~~ 0*", curve[j])
  } else {
    Cov <- NULL
  }
  
  # Add regressions
  Reg <- if (length(options[["regressions"]]) > 0)
    paste0("\n\n# Regressions\n", paste(curve, collapse = " + "), " ~ ", 
           paste(.v(options[["regressions"]]), collapse = " + "))
  else NULL
  
  # Add dummy variables
  Dum <- if (length(options[["dummy_names"]]) > 0)
    paste0("\n\n# Dummy-coded categorical predictors\n", paste(curve, collapse = " + "), " ~ ", 
           paste(.v(options[["dummy_names"]]), collapse = " + "))
  
  # Add time-varying covariates
  # eww this is hard
  
  # Put everything together
  paste0(Hed, LGC, Cov, Reg, Dum)
}

# Output functions ----
.lgcmModelContainer <- function(jaspResults, options) {
  if (!is.null(jaspResults[["modelContainer"]])) {
    modelContainer <- jaspResults[["modelContainer"]]
  } else {
    modelContainer <- createJaspContainer()
    modelContainer$dependOn(c(
      "variables", "regressions", "covariates", "categorical", "timings", 
      "intercept", "linear", "quadratic", "cubic", "covar",
      "se", "bootstrapNumber"
    ))
    jaspResults[["modelContainer"]] <- modelContainer
  }
  
  return(modelContainer)
}

.lgcmFitTable <- function(modelContainer, dataset, options, ready) {
  if (!is.null(modelContainer[["maintab"]])) return()
  maintab <- createJaspTable("Chi-square Test")
  maintab$addColumnInfo(name = "mod",    title = "Model",        type = "string")
  maintab$addColumnInfo(name = "chisq",  title = "\u03a7\u00b2", type = "number", format = "dp:3")
  maintab$addColumnInfo(name = "df",     title = "df",           type = "integer")
  maintab$addColumnInfo(name = "pvalue", title = "p",            type = "number", format = "dp:3;p:.001")
  
  modelContainer[["maintab"]] <- createJaspContainer("Model fit")
  modelContainer[["maintab"]]$position <- 1
  modelContainer[["maintab"]][["chisqtab"]] <- maintab
  
  # add data to the table!
  if (!ready) return()
  lgcmResult <- .lgcmComputeResults(modelContainer, dataset, options)
  if (modelContainer$getError()) return()
  
  fm <- lavaan::fitMeasures(lgcmResult)
  maintab[["mod"]]    <- c("Baseline model", "Growth curve model")
  maintab[["chisq"]]  <- fm[c("baseline.chisq", "chisq")]
  maintab[["df"]]     <- fm[c("baseline.df", "df")]
  maintab[["pvalue"]] <- c(NA, fm["pvalue"])
}

.lgcmParameterTables <- function(modelContainer, dataset, options, ready) {
  partabs <- if (!is.null(modelContainer[["partabs"]])) {
    modelContainer[["partabs"]]
  } else {
    modelContainer[["partabs"]] <- createJaspContainer("Parameter estimates")
  }
  
  partabs$position <- 2
  
  # create tables
  
  # latent curve
  latcur <- createJaspTable("Latent curve")
  latcur$addColumnInfo("component", title = "Component",  type = "string", combine = TRUE)
  latcur$addColumnInfo("param",     title = "Parameter",  type = "string")
  latcur$addColumnInfo("est",       title = "Estimate",   type = "number", format = "dp:3")
  latcur$addColumnInfo("se" ,       title = "Std. Error", type = "number", format = "dp:3")
  latcur$addColumnInfo("zval",      title = "z-value" ,   type = "number", format = "dp:3")
  latcur$addColumnInfo("pval",      title = "p" ,         type = "number", format = "dp:3;p:.001")
  latcur$addColumnInfo("cilo",      title = "Lower" ,     type = "number", format = "dp:3", 
                       overtitle = "95% Confidence Interval")
  latcur$addColumnInfo("ciup",      title = "Upper" ,     type = "number", format = "dp:3", 
                       overtitle = "95% Confidence Interval")
  modelContainer[["partabs"]][["latcur"]] <- latcur
  
  # covariance
  if (options[["covar"]]) {
    latcov <- createJaspTable("Latent covariances")
    latcov$addColumnInfo("cov",  title = "Covariance", type = "string")
    latcov$addColumnInfo("est",  title = "Estimate",   type = "number", format = "dp:3")
    latcov$addColumnInfo("se" ,  title = "Std. Error", type = "number", format = "dp:3")
    latcov$addColumnInfo("zval", title = "z-value" ,   type = "number", format = "dp:3")
    latcov$addColumnInfo("pval", title = "p" ,         type = "number", format = "dp:3;p:.001")
    latcov$addColumnInfo("cilo", title = "Lower" ,     type = "number", format = "dp:3", 
                         overtitle = "95% Confidence Interval")
    latcov$addColumnInfo("ciup",      title = "Upper" ,     type = "number", format = "dp:3", 
                         overtitle = "95% Confidence Interval")
    modelContainer[["partabs"]][["latcov"]] <- latcov

  }
  
  # regressions
  if (length(c(options[["regressions"]], options[["categorical"]])) > 0) {
    latreg <- createJaspTable("Regressions")
    latreg$addColumnInfo("component", title = "Component",  type = "string", combine = TRUE)
    latreg$addColumnInfo("predictor", title = "Predictor",  type = "string")
    latreg$addColumnInfo("est",       title = "Estimate",   type = "number", format = "dp:3")
    latreg$addColumnInfo("se" ,       title = "Std. Error", type = "number", format = "dp:3")
    latreg$addColumnInfo("zval",      title = "z-value" ,   type = "number", format = "dp:3")
    latreg$addColumnInfo("pval",      title = "p" ,         type = "number", format = "dp:3;p:.001")
    latreg$addColumnInfo("cilo",      title = "Lower" ,     type = "number", format = "dp:3", 
                         overtitle = "95% Confidence Interval")
    latreg$addColumnInfo("ciup",      title = "Upper" ,     type = "number", format = "dp:3", 
                         overtitle = "95% Confidence Interval")
    modelContainer[["partabs"]][["latreg"]] <- latreg
  }
  
  
  if (!ready || modelContainer$getError()) return()
  
  lgcmResult <- modelContainer[["model"]][["object"]]
  pe <- lavaan::parameterestimates(lgcmResult, level = options[["ciWidth"]])
  slope_names <- c(
    "^I$" = "Intercept", 
    "^L$" = "Linear slope", 
    "^Q$"= "Quadratic slope", 
    "^C$" = "Cubic slope"
  )
  
  pe[["lhs"]] <- stringr::str_replace_all(pe[["lhs"]], slope_names)
  pe[["rhs"]] <- stringr::str_replace_all(pe[["rhs"]], slope_names)
  
  # latent curve
  pecur <- pe[pe$lhs %in% slope_names & (pe$op == "~1" | pe$rhs == pe$lhs),]
  pecur <- pecur[order(pecur$lhs, rev(pecur$op)),]
  latcur[["component"]] <- pecur[["lhs"]]
  latcur[["param"]]     <- ifelse(pecur[["op"]] == "~1", "Mean", "Variance")
  latcur[["est"]]       <- pecur[["est"]]
  latcur[["se" ]]       <- pecur[["se"]]
  latcur[["zval"]]      <- pecur[["z"]]
  latcur[["pval"]]      <- pecur[["pvalue"]]
  latcur[["cilo"]]      <- pecur[["ci.lower"]]
  latcur[["ciup"]]      <- pecur[["ci.upper"]]
  
  # covariance
  if (options[["covar"]]) {
    pecov <- pe[pe$lhs %in% slope_names & pe$op == "~~" & pe$lhs != pe$rhs,]
    latcov[["cov"]]  <- paste(pecov[["lhs"]], "\u2B64\u00A0", pecov[["rhs"]])
    latcov[["est"]]  <- pecov[["est"]]
    latcov[["se" ]]  <- pecov[["se"]]
    latcov[["zval"]] <- pecov[["z"]]
    latcov[["pval"]] <- pecov[["pvalue"]]
    latcov[["cilo"]] <- pecov[["ci.lower"]]
    latcov[["ciup"]] <- pecov[["ci.upper"]]
  }
  
  # regressions
  if (length(c(options[["regressions"]], options[["categorical"]])) > 0) {
    pereg <- pe[pe$lhs %in% slope_names & pe$op == "~",]
    latreg[["component"]] <- pereg[["lhs"]]
    latreg[["predictor"]] <- .unv(pereg[["rhs"]])
    latreg[["est"]]       <- pereg[["est"]]
    latreg[["se" ]]       <- pereg[["se"]]
    latreg[["zval"]]      <- pereg[["z"]]
    latreg[["pval"]]      <- pereg[["pvalue"]]
    latreg[["cilo"]]      <- pereg[["ci.lower"]]
    latreg[["ciup"]]      <- pereg[["ci.upper"]]
  }
}

.lcgmAdditionalFitTables <- function(modelContainer, dataset, options, ready) {
  if (!options[["outputAdditionalFitMeasures"]]) return()
  
  fitms <- createJaspContainer(gettext("Additional Fit measures"))
  fitms$dependOn("outputAdditionalFitMeasures")
  fitms$position <- 3
  modelContainer[["fitMeasures"]] <- fitms
  
  # Fit indices
  fitms[["indices"]] <- fitin <- createJaspTable(gettext("Fit indices"))
  fitin$addColumnInfo(name = "index", title = gettext("Index"), type = "string")
  fitin$addColumnInfo(name = "value", title = gettext("Value"), type = "number", format = "sf:4;dp:3")
  fitin$setExpectedSize(rows = 1, cols = 2)
  
  # information criteria
  fitms[["incrits"]] <- fitic <- createJaspTable(gettext("Information criteria"))
  fitic$addColumnInfo(name = "index", title = "",               type = "string")
  fitic$addColumnInfo(name = "value", title = gettext("Value"), type = "number", format = "sf:4;dp:3")
  fitic$setExpectedSize(rows = 1, cols = 2)
  
  # other fit measures
  fitms[["others"]] <- fitot <- createJaspTable(gettext("Other fit measures"))
  fitot$addColumnInfo(name = "index", title = gettext("Metric"), type = "string")
  fitot$addColumnInfo(name = "value", title = gettext("Value"),  type = "number", format = "sf:4;dp:3")
  fitot$setExpectedSize(rows = 1, cols = 2)
  
  if (!ready || modelContainer$getError()) return()
  
  # actually compute the fit measures
  fm <- lavaan::fitmeasures(modelContainer[["model"]][["object"]])
  
  # Fit indices
  fitin[["index"]] <- c(
    gettext("Comparative Fit Index (CFI)"),
    gettext("Tucker-Lewis Index (TLI)"),
    gettext("Bentler-Bonett Non-normed Fit Index (NNFI)"),
    gettext("Bentler-Bonett Normed Fit Index (NFI)"),
    gettext("Parsimony Normed Fit Index (PNFI)"),
    gettext("Bollen's Relative Fit Index (RFI)"),
    gettext("Bollen's Incremental Fit Index (IFI)"),
    gettext("Relative Noncentrality Index (RNI)")
  )
  fitin[["value"]] <- fm[c("cfi", "tli", "nnfi", "nfi", "pnfi", "rfi", "ifi", "rni")]
  
  # information criteria
  fitic[["index"]] <- c(
    gettext("Log-likelihood"),
    gettext("Number of free parameters"),
    gettext("Akaike (AIC)"),
    gettext("Bayesian (BIC)"),
    gettext("Sample-size adjusted Bayesian (SSABIC)")
  )
  fitic[["value"]] <- fm[c("logl", "npar", "aic", "bic", "bic2")]
  
  # other fitmeasures
  fitot[["index"]] <- c(
    gettext("Root mean square error of approximation (RMSEA)"),
    gettext("RMSEA 90% CI lower bound"),
    gettext("RMSEA 90% CI upper bound"),
    gettext("RMSEA p-value"),
    gettext("Standardized root mean square residual (SRMR)"),
    gettextf("Hoelter's critical N (%s = .05)","\u03B1"),
    gettextf("Hoelter's critical N (%s = .01)","\u03B1"),
    gettext("Goodness of fit index (GFI)"),
    gettext("McDonald fit index (MFI)"),
    gettext("Expected cross validation index (ECVI)")
  )
  fitot[["value"]] <- fm[c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue",
                           "srmr", "cn_05", "cn_01", "gfi", "mfi", "ecvi")]
  
  return()
}


.lgcmCurvePlot <- function(modelContainer, dataset, options, ready) {
  if (!options[["curveplot"]] || !is.null(modelContainer[["curveplot"]])) return()
  curveplot <- createJaspPlot(title = "Curve plot", width = 480, height = 320)
  curveplot$dependOn(c("curveplot", "plot_categorical", "plot_max_n", "colorPalette"))
  curveplot$position <- 8
  modelContainer[["curveplot"]] <- curveplot
  if (!ready || modelContainer$getError()) return()
  
  lgcmResult <- modelContainer[["model"]][["object"]]
  plt <- .lgcmComputeCurvePlot(lgcmResult, dataset, options)
  curveplot$plotObject <- plt
}

.lgcmComputeCurvePlot <- function(lgcmResult, dataset, options) {
  N   <- lgcmResult@Data@nobs[[1]]
  P   <- length(options[["variables"]])
  idx <- 1:N
  if (N > options[["plot_max_n"]]) {
    idx <- 1:options[["plot_max_n"]]
    N   <- options[["plot_max_n"]]
  }
  ctgcl <- options[["plot_categorical"]] != ""

  if (interactive()) browser()
  
  # plot the individual-level growth curves
  preds   <- lavaan::lavPredict(lgcmResult)[idx, , drop = FALSE]
  preds   <- cbind(preds, matrix(0, nrow(preds), 4 - ncol(preds)))
  timings <- sapply(options[["timings"]], function(t) t$timing)
  xrange  <- range(timings)
  xx      <- seq(xrange[1], xrange[2], length.out = 1000)
  df_wide <- data.frame(xx = xx, apply(preds, 1, function(b) b[1] + xx*b[2] + xx^2*b[3] + xx^3*b[4]))
  df_long <- tidyr::gather(df_wide, key = "Participant", value = "Val", -"xx")
  
  if (ctgcl) df_long[[options[["plot_categorical"]]]] <- rep(dataset[[.v(options[["plot_categorical"]])]][idx], each = 1000)
  
  # create raw data points data frame
  points <- data.frame(lgcmResult@Data@X[[1]])[idx, lgcmResult@Data@ov.names[[1]] %in% .v(options[["variables"]])]
  names(points) <- timings
  points[["Participant"]] <- paste0("X", 1:nrow(points))
  points_long <- tidyr::gather(points, key = "xx", value = "Val", -"Participant")
  points_long[["xx"]] <- as.numeric(points_long[["xx"]])
  
  if (ctgcl) 
    points_long[[options[["plot_categorical"]]]] <- rep(dataset[[.v(options[["plot_categorical"]])]][idx], length(timings))
  
  # points may need to be jittered
  jitwidth <- if (N > 30) diff(range(timings) / (15 * P)) else 0
  pos <- ggplot2::position_jitter(width = jitwidth)
  
  # lines may need to be transparent
  cc <- if (ctgcl) length(unique(points_long[[options[["plot_categorical"]]]])) else 1
  transparency <- min(1, (log(cc) + 1) / log(N))
  
  # create the plot
  p <- 
    ggplot2::ggplot(df_long, ggplot2::aes(x = xx, y = Val, group = Participant)) + 
    ggplot2::geom_point(data = points_long, position = pos) +
    ggplot2::geom_line(alpha = transparency) +
    ggplot2::labs(y = "Value", x = "Time") +
    JASPgraphs::themeJaspRaw()
  
  if (ctgcl) 
    return(p + 
             ggplot2::aes_(colour = as.name(options[["plot_categorical"]]), 
                             shape  = as.name(options[["plot_categorical"]])) +
             JASPgraphs::scale_JASPcolor_discrete(options[["colorPalette"]]))
  else 
    return(p)
}

.lgcmPathPlot <- function(modelContainer, dataset, options, ready) {
  if (!options$pathplot || !is.null(modelContainer[["pathplot"]])) return()
  
  modelContainer[["pathplot"]] <- createJaspPlot(title = "Model plot", height = 500, width = 640)
  modelContainer[["pathplot"]]$dependOn(c("pathplot", "plotmeans", "plotpars"))
  modelContainer[["pathplot"]]$position <- 9
  
  if (!ready || modelContainer$getError()) return()
  
  lgcmResult <- modelContainer[["model"]][["object"]]
  png() # semplot opens a device even though we specify doNotPlot, so we hack.
  pathplot <- semPlot::semPaths(
    object         = .lavToPlotObj(lgcmResult),
    DoNotPlot      = TRUE,
    ask            = FALSE,
    layout         = "tree2",
    intercepts     = options$plotmeans,
    whatLabels     = ifelse(options$plotpars, "par", "name"),
    edge.color     = "black",
    color          = list(lat = "#EAEAEA", man = "#EAEAEA", int = "#FFFFFF"),
    border.width   = 1.5,
    edge.label.cex = 0.9,
    lty            = 2,
    title          = FALSE
  )
  dev.off()
  
  modelContainer[["pathplot"]]$plotObject <- pathplot
}

.lgcmSyntax <- function(modelContainer, dataset, options, ready) {
  if (!ready || !options[["showSyntax"]]) return()
  
  modelContainer[["model_syntax"]] <- createJaspHtml(
    text     = .lgcmOptionsToMod(options, FALSE), 
    class    = "jasp-code", 
    position = 10, 
    title    = "Model Syntax"
  )
}


# Unused functions ----
.lgcmMisfitPlot <- function(modelContainer, dataset, options, ready) {
  if (!options[["misfitplot"]] || !is.null(modelContainer[["misfitplot"]])) return()
  wh <- 50 + 50*length(options[["variables"]])
  misplot <- createJaspPlot(title = "Misfit plot", width = wh, height = wh)
  misplot$dependOn("misfitplot")
  misplot$position <- 9
  modelContainer[["misfitplot"]] <- misplot
  if (!ready || modelContainer$getError()) return()
  
  lgcmResult <- modelContainer[["model"]][["object"]]
  rescor <- lavaan::residuals(lgcmResult, type = "cor")
  cc <- rescor[["cov"]]
  cc[upper.tri(cc)] <- NA
  gg <- .resCorToMisFitPlot(cc)
  misplot$plotObject <- gg
}

.lgcmPlotRibbon <- function(lgcmResult, options) {
  # plot uncertainty ribbon conditional on regressors == 0
  
  # get parameter values
  pe <- lavaan::parameterestimates(lgcmResult)
  mu <- pe[pe$lhs %in% c("I", "L", "Q", "C") & pe$rhs == "", ]
  addrow <- matrix(0, 4 - nrow(mu), ncol(mu))
  colnames(addrow) <- names(mu)
  mu <- rbind(mu, addrow)
  s2 <- pe[pe$lhs %in% c("I", "L", "Q", "C") & pe$lhs == pe$rhs,]
  s2 <- rbind(s2, addrow)
  
  # create x range
  timings <- sapply(options[["timings"]], function(t) t$timing)
  xrange  <- range(timings)
  xx      <- seq(xrange[1], xrange[2], length.out = 1000)
  
  # inner ribbon (with only parameter uncertainty, no variance)
  # growth curve for a typical person with covariates at 0
  mu_mu <- mu$est[1] + xx*mu$est[2] + xx^2*mu$est[3] + xx^3*mu$est[4]
  mu_up <- mu$ci.upper[1] + xx*mu$ci.upper[2] + xx^2*mu$ci.upper[3] + xx^3*mu$ci.upper[4]
  mu_lo <- mu$ci.lower[1] + xx*mu$ci.lower[2] + xx^2*mu$ci.lower[3] + xx^3*mu$ci.lower[4]
  
  
  # growth curve for draws as if our group were the population
  mp <- qnorm(0.975)
  s2_up <- (mu$est[1] + mp*sqrt(s2$est[1])) + 
    xx   * (mu$est[2] + mp*sqrt(s2$est[2])) + 
    xx^2 * (mu$est[3] + mp*sqrt(s2$est[3])) + 
    xx^3 * (mu$est[4] + mp*sqrt(s2$est[4]))
  s2_lo <- (mu$est[1] - mp*sqrt(s2$est[1])) + 
    xx   * (mu$est[2] - mp*sqrt(s2$est[2])) + 
    xx^2 * (mu$est[3] - mp*sqrt(s2$est[3])) + 
    xx^3 * (mu$est[4] - mp*sqrt(s2$est[4]))
  
  # growth curve for groups of people with covariates 0 and parameter uncertainty
  ms_up <- (mu$ci.upper[1] + mp*sqrt(s2$ci.upper[1])) + 
    xx   * (mu$ci.upper[2] + mp*sqrt(s2$ci.upper[2])) + 
    xx^2 * (mu$ci.upper[3] + mp*sqrt(s2$ci.upper[3])) + 
    xx^3 * (mu$ci.upper[4] + mp*sqrt(s2$ci.upper[4]))
  ms_lo <- (mu$ci.lower[1] - mp*sqrt(s2$ci.upper[1])) + 
    xx   * (mu$ci.lower[2] - mp*sqrt(s2$ci.upper[2])) + 
    xx^2 * (mu$ci.lower[3] - mp*sqrt(s2$ci.upper[3])) + 
    xx^3 * (mu$ci.lower[4] - mp*sqrt(s2$ci.upper[4]))
  
  fac <- forcats::fct_rev(forcats::as_factor(rep(c("mu", "s2", "ms"), each = 1000)))
  
  df    <- data.frame(xx = rep(xx, 3), up = c(mu_up, s2_up, ms_up), lo = c(mu_lo, s2_lo, ms_lo), which = fac)
  mu_df <- data.frame(xx = xx, y = mu_mu) 
  
  # create raw data points data frame
  points <- data.frame(lgcmResult@Data@X[[1]])[idx, lgcmResult@Data@ov.names[[1]] %in% options[["variables"]]]
  names(points) <- timings
  points[["Participant"]] <- paste0("X", 1:nrow(points))
  points_long <- tidyr::gather(points, key = "xx", value = "Val", -"Participant")
  points_long[["xx"]] <- as.numeric(points_long[["xx"]])
  
  # points may need to be jittered
  jitwidth <- if (N > 30) diff(range(timings) / (15 * P)) else 0
  pj <- ggplot2::position_jitter(width = jitwidth)
  
  ggplot2::ggplot(df, mapping = ggplot2::aes(x = xx)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymax = up, ymin = lo, fill = which)) +
    ggplot2::geom_line(data = mu_df, col = "#454545", ggplot2::aes(y = y)) +
    ggplot2::geom_point(ggplot2::aes(y = Val), data = points_long, position = pj, col = "#454545") +
    ggplot2::scale_fill_manual(
      values = c("#ABABAB", "#909090", "#777777"),
      guide  = 'legend', 
      labels = c("Curve + variance + uncertainty", 
                 "Curve + variance", 
                 "Curve + parameter uncertainty")
    ) +
    labs(fill = "", x = "Time", y = "Value") +
    theme_minimal()
}
