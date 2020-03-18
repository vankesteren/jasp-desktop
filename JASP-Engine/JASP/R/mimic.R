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

MIMIC <- function(jaspResults, dataset, options, ...) {
  jaspResults$addCitation("Rosseel, Y. (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36. URL http://www.jstatsoft.org/v48/i02/")
  
  # Read dataset
  dataset <- .mimicReadData(dataset, options)
  ready   <- .mimicCheckErrors(dataset, options)
  
  modelContainer <- .mimicModelContainer(jaspResults)
  
  # Output functions
  .mimicParTable(modelContainer, dataset, options, ready)
  .mimicRsquared(modelContainer, options, ready)
  .mimicPathPlot(modelContainer, options, ready)
  .mimicSyntax(  modelContainer, options, ready)
  
}

# Preprocessing functions ----
.mimicReadData <- function(dataset, options) {
  if (!is.null(dataset)) return(dataset)
  
  vars <- c(options$predictors, options$indicators)
  return(.readDataSetToEnd(columns = vars))
}

.mimicCheckErrors <- function(dataset, options) {
  if (length(options$indicators) < 3 || length(options$predictors) == 0) return(FALSE)
  
  # Check for missing value handling
  if (options$estimator %in% c("GLS", "WLS", "ULS", "DWLS") && options$missing == "fiml")
    .quitAnalysis(gettext("FIML only available with ML-type estimators."))
  
  # Exogenous variables can be binary or continuous
  exo <- options$predictors
  # Endogenous variables need to be scale or ordinal
  endo <- options$indicators
  
  customChecks <- list(
    checkExogenous = function() {
      admissible <- vapply(exo, function(exo_var) {
        var <- na.omit(dataset[[.v(exo_var)]])
        if (is.ordered(var)) return(FALSE)
        if ((is.character(var) || is.factor(var)) && length(unique(var)) != 2) return(FALSE)
        return(TRUE)
      }, TRUE)
      if (!all(admissible))
        gettextf("Not all exogenous variables are admissible. Inadmissible exogenous variables: %s. Only binary or continuous exogenous variables allowed.", paste(exo[!admissible], collapse = ", "))
    },
    
    checkEndogenous = function() {
      if (length(options$confounds) > 0) endo <- c(endo, options$predictor)
      admissible <- vapply(endo, function(endo_var) {
        var <- na.omit(dataset[[.v(endo_var)]])
        if (!(is.ordered(var) || is.numeric(var))) {
          return(FALSE)
        }
        return(TRUE)
      }, TRUE)
      if (!all(admissible))
        gettextf("Not all endogenous variables are admissible. Inadmissible endogenous variables: %s. Only scale or ordinal endogenous variables allowed.", paste(endo[!admissible], collapse = ", "))
    },
    
    checkCategoricalEndo = function() {
      if (length(options$confounds) > 0) endo <- c(endo, options$predictor)
      
      admissible <- vapply(endo, function(endo_var) {
        var <- na.omit(dataset[[.v(endo_var)]])
        if (is.ordered(var) && options$missing == "fiml") {
          return(FALSE)
        }
        return(TRUE)
      }, TRUE)
      
      if (!all(admissible))
        gettextf("FIML missing value handling only available when all endogenous variables are of scale type. Ordinal endogenous variables in the model: %s", paste(endo[!admissible], collapse = ", "))
    }
    
  )
  
  .hasErrors(dataset, type = c('observations', 'variance', 'infinity'), custom = customChecks, 
             all.target = c(endo, exo), observations.amount = paste('<', length(c(endo, exo))), 
             exitAnalysisIfErrors = TRUE)
  
  return(TRUE)
}

# Results functions ----

.mimicComputeResults <- function(modelContainer, dataset, options, ready) {
  mimicResult <- try(lavaan::sem(
    model           = .mimicToLavMod(options),
    data            = dataset,
    se              = ifelse(options$se == "bootstrap", "standard", options$se),
    mimic           = options$mimic,
    estimator       = options$estimator,
    missing         = options$missing,
    std.lv          = TRUE
  ))
  
  if (inherits(mimicResult, "try-error")) {
    errmsg <- gettextf("Estimation failed\nMessage:\n%s", attr(mimicResult, "condition")$message)
    modelContainer$setError(.decodeVarsInMessage(names(dataset), errmsg))
  }
  
  if (options$se == "bootstrap") {
    mimicResult <- lavBootstrap(mimicResult, options$bootstrapNumber)
  }
  
  modelContainer[["model"]] <- createJaspState(mimicResult)
  return(mimicResult)
}

.mimicToLavMod <- function(options, base64 = TRUE) {
  
  if (!base64) .v <- I
  
  header <- "
  # -----------------------------
  # MIMIC model generated by JASP
  # -----------------------------
  "
  measurement <- paste(
    "  Y =~", 
    paste(
      paste("lambda", seq_along(options[["indicators"]]), "*", sep = ""),
      .v(options[["indicators"]]), 
      collapse = " + ", sep = ""
    )
  )
  structural <- paste(
    "  Y ~ ",
    paste(
      paste("beta", seq_along(options[["predictors"]]), "*", sep = ""),
      .v(options[["predictors"]]), 
      collapse = " + ", sep = ""
    )
  )
  
  return(paste(header, measurement, structural, sep = "\n"))
}

# Output functions ----

.mimicModelContainer <- function(jaspResults, options) {
  if (!is.null(jaspResults[["modelContainer"]])) {
    modelContainer <- jaspResults[["modelContainer"]]
  } else {
    modelContainer <- createJaspContainer()
    modelContainer$dependOn(c(
      "predictors", "indicators", "includemeanstructure", 
      "bootstrapNumber", "fixManifestInterceptsToZero", "mimic", "se", "estimator", 
      "std", "missing")
    )
    jaspResults[["modelContainer"]] <- modelContainer
  }
  
  return(modelContainer)
}

.mimicParTable <- function(modelContainer, dataset, options, ready) {
  if (!is.null(modelContainer[["parest"]])) return()
  modelContainer[["parest"]] <- pecont <- createJaspContainer(gettext("Parameter estimates"))
  pecont$dependOn(options = c("ciWidth", "bootCItype"))
  pecont$position <- 0
  
  ## betas
  bettab <- createJaspTable(title = gettext("Predictor coefficients"))
  
  bettab$addColumnInfo(name = "rhs",      title = gettext("Predictor"),  type = "string")
  bettab$addColumnInfo(name = "est",      title = gettext("Estimate"),   type = "number", format = "sf:4;dp:3")
  bettab$addColumnInfo(name = "se",       title = gettext("Std. Error"), type = "number", format = "sf:4;dp:3")
  bettab$addColumnInfo(name = "z",        title = gettext("z-value"),    type = "number", format = "sf:4;dp:3")
  bettab$addColumnInfo(name = "pvalue",   title = gettext("p"),          type = "number", format = "dp:3;p:.001")
  bettab$addColumnInfo(name = "ci.lower", title = gettext("Lower"),      type = "number", format = "sf:4;dp:3",
                       overtitle = gettextf("%s%% Confidence Interval", options$ciWidth * 100))
  bettab$addColumnInfo(name = "ci.upper", title = gettext("Upper"),      type = "number", format = "sf:4;dp:3",
                       overtitle = gettextf("%s%% Confidence Interval", options$ciWidth * 100))
  
  pecont[["bet"]] <- bettab
  
  ## lambdas
  lamtab <- createJaspTable(title = gettext("Indicator coefficients"))
  
  lamtab$addColumnInfo(name = "rhs",      title = gettext("Predictor"),  type = "string")
  lamtab$addColumnInfo(name = "est",      title = gettext("Estimate"),   type = "number", format = "sf:4;dp:3")
  lamtab$addColumnInfo(name = "se",       title = gettext("Std. Error"), type = "number", format = "sf:4;dp:3")
  lamtab$addColumnInfo(name = "z",        title = gettext("z-value"),    type = "number", format = "sf:4;dp:3")
  lamtab$addColumnInfo(name = "pvalue",   title = gettext("p"),          type = "number", format = "dp:3;p:.001")
  lamtab$addColumnInfo(name = "ci.lower", title = gettext("Lower"),      type = "number", format = "sf:4;dp:3",
                       overtitle = gettextf("%s%% Confidence Interval", options$ciWidth * 100))
  lamtab$addColumnInfo(name = "ci.upper", title = gettext("Upper"),      type = "number", format = "sf:4;dp:3",
                       overtitle = gettextf("%s%% Confidence Interval", options$ciWidth * 100))
  
  pecont[["lam"]] <- lamtab
  
  if (!ready) return()
  
  # add data to the tables!
  .mimicComputeResults(modelContainer, dataset, options, ready)
  
  if (modelContainer$getError()) return()
  
  
  pe <- lavaan::parameterEstimates(modelContainer[["model"]][["object"]], 
                                   boot.ci.type = options$bootCItype, 
                                   level = options$ciWidth)
  
  pe_bet <- pe[substr(pe$label, 1, 1) == "b", ]
  bettab[["rhs"]]      <- .unv(pe_bet$rhs)
  bettab[["est"]]      <- pe_bet$est
  bettab[["se"]]       <- pe_bet$se
  bettab[["z"]]        <- pe_bet$z
  bettab[["pvalue"]]   <- pe_bet$pvalue
  bettab[["ci.lower"]] <- pe_bet$ci.lower
  bettab[["ci.upper"]] <- pe_bet$ci.upper
  
  pe_lam <- pe[substr(pe$label, 1, 1) == "l", ]
  lamtab[["rhs"]]      <- .unv(pe_lam$rhs)
  lamtab[["est"]]      <- pe_lam$est
  lamtab[["se"]]       <- pe_lam$se
  lamtab[["z"]]        <- pe_lam$z
  lamtab[["pvalue"]]   <- pe_lam$pvalue
  lamtab[["ci.lower"]] <- pe_lam$ci.lower
  lamtab[["ci.upper"]] <- pe_lam$ci.upper
}

.mimicRsquared <- function(modelContainer, options, ready) {
  if (!options$rsquared || !is.null(modelContainer[["rsquared"]])) return()
  
  tabr2 <- createJaspTable(gettext("R-Squared"))
  tabr2$addColumnInfo(name = "__var__", title = "", type = "string")
  tabr2$addColumnInfo(name = "rsq", title = "R\u00B2", type = "number", format = "sf:4;dp:3")
  tabr2$dependOn(options = "rsquared")
  tabr2$position <- 1
  
  modelContainer[["rsquared"]] <- tabr2 
  
  if (!ready || modelContainer$getError()) return()
  
  r2res              <- lavaan::inspect(modelContainer[["model"]][["object"]], "r2")
  tabr2[["__var__"]] <- .unv(names(r2res))
  tabr2[["rsq"]]     <- r2res
}

.mimicPathPlot <- function(modelContainer, options, ready) {
  if (!options$pathplot || !ready || !is.null(modelContainer[["plot"]])) return()
  
  plt <- createJaspPlot(title = gettext("Path plot"), width = 600, height = 400)
  plt$dependOn(options = c("pathplot", "plotpars", "plotlegend"))
  plt$position <- 2
  
  modelContainer[["plot"]] <- plt
  
  if (modelContainer$getError()) return()
  
  # create a qgraph object using semplot
  po <- .lavToPlotObj(modelContainer[["model"]][["object"]])
  pp <- .suppressGrDevice(semPlot::semPaths(
    object         = po,
    layout         = "tree2",
    intercepts     = FALSE,
    reorder        = FALSE,
    whatLabels     = ifelse(options$plotpars, "par", "name"),
    edge.color     = "black",
    color          = list(lat = "#EAEAEA", man = "#EAEAEA", int = "#FFFFFF"),
    title          = FALSE,
    legend         = options$plotlegend,
    legend.mode    = "names",
    legend.cex     = 0.6,
    nodeNames      = decodeColNames(po@Vars$name),
    nCharNodes     = 3,
    normalize      = TRUE
  ))
  
  plt$plotObject <- pp
}

.mimicSyntax <- function(modelContainer, options, ready) {
  if (!options$showSyntax || !ready) return()
  modelContainer[["syntax"]] <- createJaspHtml(.mimicToLavMod(options, FALSE), class = "jasp-code", title = gettext("Model syntax"))
  modelContainer[["syntax"]]$dependOn("showSyntax")
  modelContainer[["syntax"]]$position <- 3
}