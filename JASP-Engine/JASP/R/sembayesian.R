#
# Copyright (C) 2013-2018 University of Amsterdam
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

SEMBayes <- function(jaspResults, dataset, options, ...) {
  jaspResults$addCitation("Merkle EC, Rosseel Y (2018). blavaan: Bayesian Structural Equation Models via Parameter Expansion. Journal of Statistical Software, 85(4), 1-30. doi: 10.18637/jss.v085.i04")

  # Read dataset
  dataset <- .bsemReadData(dataset, options)
  ready   <- .bsemCheckErrors(dataset, options)
  
  modelContainer <- .bsemModelContainer(jaspResults)
  
  # destroy existing shiny process
  
  
  # Output functions
  .bsemParTable(  modelContainer, dataset, options, ready)
  .bsemParDiag(   modelContainer, dataset, options, ready)
  .bsemTracePlots(modelContainer, dataset, options, ready)
  .bsemPostPlots( modelContainer, dataset, options, ready)
  
  # .bsemShinyLaunch(modelContainer, dataset, options, ready)
}


# Preprocessing functions ----
.bsemReadData <- function(dataset, options) {
  if (!is.null(dataset)) return(dataset)
  return(.readDataSetToEnd(all.columns = TRUE))
}

.bsemCheckErrors <- function(dataset, options) {
  nchar(options[["model"]]) > 0
}

# Results functions ----

.bsemComputeResults <- function(modelContainer, dataset, options, ready) {
  mod <- .bsemToMod(options, dataset)
  blavaan <- blavaan::blavaan
  bsemResult <- try(blavaan::bsem(
    model    = mod, 
    data     = dataset, 
    burnin   = options[["burnin"]],
    sample   = options[["nsamples"]],
    n.chains = options[["nchains"]]
  ))
  
  if (inherits(bsemResult, "try-error")) {
    errmsg <- paste("Estimation failed\nMessage:\n", attr(bsemResult, "condition")$message)
    modelContainer$setError(.decodeVarsInMessage(names(dataset), errmsg))
  }
  
  modelContainer[["model"]] <- createJaspState(bsemResult)
  return(bsemResult)
}

.bsemToMod <- function(options, dataset) {
  mod      <- options[["model"]]
  usedVars <- .bsemUsedVars(mod, names(dataset))
  return(.bsemTranslateModel(mod, usedVars))
}

.bsemUsedVars <- function(model, vars) {
  vv <- .unv(vars)
  findpattern <- paste0("(?<=[\\s\\+\\^\\=\\~\\<\\*\\>\\:\\%\\|\\+]|^)\\Q",
                        vv,
                        "\\E(?=[\\s\\+\\^\\=\\~\\<\\*\\>\\:\\%\\|\\+]|$)")
  return(vv[vapply(findpattern,
                   function(p) stringr::str_detect(model, p),
                   FUN.VALUE = TRUE,
                   USE.NAMES = FALSE)])
}

.bsemTranslateModel <- function(model, variables) {
  if (length(variables) == 0) {
    return(model)
  }
  
  variables <- variables[order(nchar(variables), decreasing = TRUE)]
  with.s.quotes <- paste("\\b'", variables, "'\\b", sep="")
  with.d.quotes <- paste('\\b"', variables, '"\\b', sep="")
  
  new.names <- .v(variables)
  
  for (i in 1:length(variables)) {
    model <- gsub(with.d.quotes[i], new.names[i], model)
  }
  
  for (i in 1:length(variables)) {
    model <- gsub(with.s.quotes[i], new.names[i], model)
  }
  
  for (i in 1:length(variables)) {
    model <- gsub(paste0("\\b",variables[i], "\\b"), new.names[i], model)
  }
  
  return(model)
}

# Output functions ----
.bsemModelContainer <- function(jaspResults) {
  if (!is.null(jaspResults[["modelContainer"]])) {
    modelContainer <- jaspResults[["modelContainer"]]
  } else {
    modelContainer <- createJaspContainer()
    modelContainer$dependOn(c("model", "nchains", "nsamples", "burnin"))
    jaspResults[["modelContainer"]] <- modelContainer
  }
  
  return(modelContainer)
}

.bsemParTable <- function(modelContainer, dataset, options, ready) {
  if (!is.null(modelContainer[["partab"]])) return()
  
  ## parameter summary
  partab <- createJaspTable(title = "Parameter estimates")
  
  partab$addColumnInfo(name = "cat",    title = "",       type = "string", combine = TRUE)
  partab$addColumnInfo(name = "par",    title = "",       type = "string")
  partab$addColumnInfo(name = "mean",   title = "Mean",   type = "number", format = "sf:4;dp:3")
  partab$addColumnInfo(name = "median", title = "Median", type = "number", format = "sf:4;dp:3")
  partab$addColumnInfo(name = "lower",  title = "Lower",  type = "number", format = "sf:4;dp:3", 
                       overtitle = "95% HPD interval")
  partab$addColumnInfo(name = "upper",  title = "Upper",  type = "number", format = "sf:4;dp:3", 
                       overtitle = "95% HPD interval")
  
  partab$position <- 1
  
  modelContainer[["partab"]] <- partab
  
  if (!ready) return()
  
  # add data to the tables!
  bsemResult <- .bsemComputeResults(modelContainer, dataset, options, ready)
  
  if (modelContainer$getError()) return()
  
  hpd <- blavaan::blavInspect(bsemResult, "hpd")
  
  parinfo  <- .bsemParseParNames(rownames(hpd))
  parsort  <- order(parinfo[2,])
  parnames <- sapply(rownames(hpd), .decodeVarsInMessage, encodedVars = names(dataset))
  
  partab[["cat"]]   <- parinfo[1, parsort]
  partab[["par"]]   <- parnames[parsort]
  partab[["mean"]]  <- blavaan::blavInspect(bsemResult, "postmean")[parsort]
  partab[["median"]]<- blavaan::blavInspect(bsemResult, "postmedian")[parsort]
  partab[["lower"]] <- hpd[parsort, "lower"]
  partab[["upper"]] <- hpd[parsort, "upper"]
}

.bsemParseParNames <- function(parnames) {
  vapply(parnames, function(pn) {
    op <- na.omit(stringr::str_extract(pn, c("~~", "~1", "=~", "~")))[1]
    if (is.na(op)) {
      split <- c(pn, "") 
    } else {
      split <-  stringr::str_split_fixed(pn, op, 2)
    }
    if (is.na(op)) {
      category <- "Named"
      position <- 1
    } else if (op == "=~") {
      category <- "Loadings"
      position <- 2
    } else if (op == "~") {
      category <- "Regressions"
      position <- 3
    }else if (op == "~~") {
      if (split[1] == split[2]) {
        category <- "Variances"
        position <- 4
      } else {
        category <- "Covariances"
        position <- 5
      }
    } else if (op == "~1") {
      category <- "Intercepts"
      position <- 6
    }
    c(category, position)
  }, character(2))
}

.bsemParDiag <- function(modelContainer, dataset, options, ready) {
  if (!options[["pardiag"]] || !is.null(modelContainer[["diatab"]])) return()
  diatab <- createJaspTable(title = "Parameter diagnostics")
  
  diatab$addColumnInfo(name = "cat",  title = "",     type = "string", combine = TRUE)
  diatab$addColumnInfo(name = "par",  title = "",     type = "string")
  diatab$addColumnInfo(name = "psrf", title = "PSRF", type = "number", format = "sf:4;dp:3")
  diatab$addColumnInfo(name = "neff", title = "NEFF", type = "number", format = "sf:4;dp:3")
  diatab$addFootnote(message = paste("PSRF: Gelman-Rubin potential scale reduction factor,", 
                                     "NEFF: Effective sample size"))
  diatab$position <- 2
  modelContainer[["diatab"]] <- diatab
  
  if (!ready || modelContainer$getError()) return()
  
  bsemResult <- modelContainer[["model"]][["object"]]
  
  psrf     <- blavaan::blavInspect(bsemResult, "psrf")
  neff     <- blavaan::blavInspect(bsemResult, "neff")
  parinfo  <- .bsemParseParNames(names(psrf))
  parsort  <- order(parinfo[2,])
  parnames <- sapply(names(psrf), .decodeVarsInMessage, encodedVars = names(dataset))
  
  diatab[["cat"]]  <- parinfo[1, parsort]
  diatab[["par"]]  <- parnames[parsort]
  diatab[["psrf"]] <- psrf[parsort]
  diatab[["neff"]] <- neff[parsort]
}

.bsemTracePlots <- function(modelContainer, dataset, options, ready) {
  if (!options[["traceplot"]] || !is.null(modelContainer[["traceplot"]]) || !ready) return()
  bsemResult <- modelContainer[["model"]][["object"]]
  traceplot <- blavaan:::plot.blavaan(bsemResult)
  traceplot$position <- 3
  parnames <- sapply(names(blavaan::blavInspect(bsemResult, "postmean")), .decodeVarsInMessage, 
                     encodedVars = names(dataset))
  levels(traceplot[["data"]][["parameter"]]) <- parnames
  modelContainer[["traceplot"]] <- createJaspPlot(traceplot, "Trace plot", width = 640, height = 320)
  modelContainer[["traceplot"]]$dependOn("traceplot")
}

.bsemPostPlots <- function(modelContainer, dataset, options, ready) {
  if (!options[["postplot"]] || !is.null(modelContainer[["postplot"]]) || !ready) return()
  bsemResult <- modelContainer[["model"]][["object"]]
  postplot <- blavaan:::plot.blavaan(bsemResult, plot.type = "dens", fill = "dark gray")
  postplot$position <- 4
  parnames <- sapply(names(blavaan::blavInspect(bsemResult, "postmean")), .decodeVarsInMessage, 
                     encodedVars = names(dataset))
  levels(postplot[["data"]][["parameter"]]) <- parnames
  modelContainer[["postplot"]] <- createJaspPlot(postplot, "Posterior plot", width = 640, height = 320)
  modelContainer[["postplot"]]$dependOn("postplot")
}

# .bsemShinyLaunch <- function(modelContainer, dataset, options, ready) {
#   if (modelContainer$getError() || !ready) return()
#   blav <- modelContainer[["model"]][["object"]]
#   
#   sso <- shinystan::as.shinystan(blav@external$mcmcout)
#   .launchShinyApp(sso)
# }

# .launchShinyApp <- function(sso) {
#     .SHINYSTAN_OBJECT <<- sso
#     on.exit(.SHINYSTAN_OBJECT <<- NULL, add = TRUE)
#     
#     p <- processx::process$new("Rscript", c('-e','shiny::runApp(system.file("ShinyStan", package = "shinystan"), launch.browser = FALSE, port = 8080, host = "127.0.0.1", display.mode = "normal")'),
#                                cleanup = TRUE)
# }

