# 1 Copyright
#
# This file is part of PFunc. PFunc provides a set of simple tools for users
# to analyze preference functions and other function-valued traits.
#
# Copyright 2016-2022 Joseph Kilmer
#
# PFunc is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PFunc is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# 2 Author comment
#     Thanks for using my splining preference function script. You can contact
#     me at joeykilmer@gmail.com with comments or questions.
#     Visit http://github.com/joccalor/pfunc to find the latest version.

# 3 source() and library() statements
#
# Read this whole file into R with the following command:
# source("PFunc_RCode.R")

if (!"mgcv" %in% installed.packages()) {
  stop("Please install the package 'mgcv' to use this script.
        Use the command: install.packages('mgcv')")
}

if (!"package:mgcv" %in% search()) {
  library(mgcv)
}

# 4 Functions
PFunc<-function(input.data, diagnose.col = 0, diagnose.sp = -1,
                auto.sort = TRUE,
                blocklist = NULL, allowlist = NULL,
                sp.binding = TRUE, sp.assign = 0,
                max.sp = 5, min.sp = 0.05, k = -1,
                peak.within = 1,
                summary.out = "spline_summaries.csv",
                graph.out = "spline_graphs.pdf",
                points.out = FALSE,
                max.y = 0, pdf.row = 4, pdf.col = 3,
                graph.points = TRUE, graph.peak = TRUE, graph.tol = TRUE,
                graph.sp = FALSE, graph.se = FALSE,
                drop = 1/3, tol.mode = "broad", tol.floor = 0,
                n.predictions = 0, ghost = FALSE, allfromsplines = TRUE,
                forgui = FALSE) {
  # This is the primary function to call. All others below are secondary.
  # See the README file for argument definitions and examples of use.
  k <- CheckValues(input.data, k, blocklist,
                   allowlist, sp.assign, max.sp, min.sp,
                   peak.within, forgui)
  if (max.y == 0) {
    max.y <- CalculateMaxY(input.data)
  }
  orig.input.names <- names(input.data)

  input.data <- BlockOrAllowList(input.data, blocklist, allowlist)

  names(input.data)[1] <- "stimulus"
  input.stimuli <- input.data$stimulus

  if (auto.sort == TRUE) {
    input.data <- input.data[order(input.data$stimulus), ]
  }

  pred.x.vals <- GeneratePredictionLocations(input.data, n.predictions)

  if (diagnose.col > 1) {
      Diagnose(input.data, diagnose.col, diagnose.sp, peak.within, drop,
               tol.mode, max.y, pred.x.vals, allfromsplines, k,
               sp.binding, min.sp, max.sp, sp.assign, graph.points,
               graph.se, forgui, tol.floor, points.out)
  } else {
    output<-data.frame(name = names(input.data[2:length(names(input.data))]),
                       peak_pref=NA, peak_height=NA, tolerance=NA,
                       HD_strength=NA, HI_strength=NA,
                       responsiveness=NA, smoothing=NA)
    if (graph.out != FALSE) {
      pdf(file = graph.out, onefile=T)
      par(mfrow=c(pdf.row, pdf.col), mar = c(1.5, 1.1, 2, 1.1),
          oma = c(1, 1.5, 0, .5))
    }
    for (response.column in 2:ncol(input.data)) {
      response.column.name <- names(input.data)[response.column]
      orig.col.num <- which(orig.input.names == response.column.name)
      output.row <- response.column - 1

      is.flat = CheckForFlat(input.data, response.column)

      preference.function <- gam(input.data[, response.column] ~
                                 s(stimulus, k = k),
                                 data = input.data, scale = -1)
      smoothing.parameter <- preference.function$sp

      if (sp.binding == TRUE) {
        smoothing.parameter <- SPBinding(smoothing.parameter, max.sp, min.sp)
      }
      ghost.bundle <- list(NULL)
      if (is.matrix(sp.assign)) {
        if (ghost == TRUE) {
          ghost.bundle <- Ghost(input.data, response.column, stimulus, k,
                                preference.function, smoothing.parameter,
                                orig.col.num, sp.assign, input.stimuli,
                                peak.within, drop, is.flat, tol.floor)
        }
        smoothing.parameter <- SPAssign(smoothing.parameter, orig.col.num,
                                        sp.assign)
      }
      preference.function <- gam(input.data[, response.column] ~
                                 s(stimulus, k = k), data = input.data,
                                 scale = -1, sp = smoothing.parameter)

  # Predicted Points
      pred.y.vals <- predict.gam(preference.function, pred.x.vals,
                                 se.fit = FALSE)
      predicted.points <- cbind(pred.x.vals, pred.y.vals)
      names(predicted.points)[ncol(predicted.points)] <- names(input.data)[response.column]

  # Peak
      peak.bundle <- Peak(input.stimuli, preference.function,
                          peak.within, is.flat)

  # Tolerance
      tolerance.bundle <- Tolerance(drop, peak.bundle, is.flat,
                                    preference.function, tol.floor)

  # Graphing
      if (graph.out != FALSE) {
        GraphSpline(input.data, peak.bundle, tolerance.bundle,
                    response.column.name, response.column, graph.points,
                    graph.peak, graph.tol, graph.sp, tol.mode, max.y, ghost,
                    ghost.bundle, is.flat, smoothing.parameter, orig.col.num)
      }

  # Strength and Responsiveness
      str.resp.bundle <- StrResp(input.data[,response.column],
                                 predicted.points[,ncol(predicted.points)],
                                 is.flat, allfromsplines)

  # Output
      if (tol.mode == "strict") {
        final.tol <- tolerance.bundle$strict.tolerance
      }
      else {
        final.tol <- tolerance.bundle$broad.tolerance
      }

      output[output.row, 2:ncol(output)] <- c(peak.bundle$peak.preference,
                                              peak.bundle$peak.response,
                                              final.tol,
                                              str.resp.bundle$hd.strength,
                                              str.resp.bundle$hi.strength,
                                              str.resp.bundle$responsiveness,
                                              round(smoothing.parameter, 3))
    }

    if (summary.out != FALSE) {
      write.csv(output, summary.out, row.names = FALSE, na = "")
    }

    if (graph.out != FALSE) {
      dev.off()
    }

    if (points.out != FALSE){
      write.csv(predicted.points, points.out, row.names = FALSE)
    }
  }
}


# Secondary Function Definitions:

CheckValues <- function(input.data, k, blocklist, allowlist,
                        sp.assign, max.sp, min.sp,
                        peak.within, forgui) {
  # This function checks over several of the settings and helps the user
  # resolve conflicts.
  if (ncol(input.data) < 2) {
    stop("    Input data file must have at least 2 columns.
      Check that you imported it correctly.")
  }

  if (nrow(input.data) <= 2) {
    stop("Input data must have at least 2 rows")
  }

  if (length(unique(input.data[,1])) < 10) {
    if (!forgui) {
      print("    Your input data has fewer than ten rows.
      As a result, the degrees of freedom for the smoothing
      function will be automatically adjusted (see 'choose.k').
      In most cases, this is not a problem.
      Consider lowering 'min.sp'.")
    }
    k <- length(unique(input.data[,1]))
  }

  if (is.vector(blocklist) & is.vector(allowlist)) {
    stop("You may not use block.list and allow.list simultaneously")
  }

  if (max.sp < min.sp) {
    stop("max.sp cannot be smaller than min.sp")
  }

  if (length(peak.within) == 1) {
    if (peak.within <= 0 | peak.within > 1) {
      stop("peak.within must be greater than 0 and less than or equal to 1")
    }
  }
  return(k)
}


CalculateMaxY <- function(input.data) {
  # Automatically calculates maximum y-value for graphs if none is specified
  max.resp <- max(input.data[ , 2:ncol(input.data)])
  min.resp <- min(input.data[ , 2:ncol(input.data)])
  resp.range <- max.resp - min.resp
  max.y <- max.resp + 0.02 * resp.range
  return(max.y)
}


BlockOrAllowList <- function(input.data, blocklist, allowlist) {
  # Excludes individuals listed in the blocklist OR includes only those
  # individuals in the allowlist. The two are not to be used simultaneously.
  # When they are, allowlist takes precedence.
  all.names <- names(input.data)
  all.columns <- 1:ncol(input.data)
  used.list <- all.columns
  first.column <- input.data[1]

  if (is.vector(blocklist)) {
    used.list <- blocklist
  }

  if (is.vector(allowlist)) {
    used.list <- allowlist
  }

  if (is.character(used.list)) {
    used.columns <- vector(length = length(used.list))
    for (i in 1:length(used.list)) {
      used.columns[i] <- which(all.names == used.list[i])
    }
  }
  else {
    used.columns <- used.list
  }

  if (is.vector(blocklist)) {
    used.columns <- which(!all.columns %in% used.columns)
  }

  input.data <- input.data[used.columns]
  if (!all.names[1] %in% names(input.data)) {
    input.data <- cbind(first.column, input.data)
  }
  return(input.data)
}


GeneratePredictionLocations <- function(input.data, n.predictions) {
  # Creates the x-axis values for the predict.gam function either at the
  # exact same x-axis values as the input data, or at a set of evenly spaced
  # points, the number of which is dictated by the variable "n.predictions"
  if (n.predictions < 2) {
    pred.x.vals <- input.data[1]
  } else {
    pred.x.vals <- data.frame(stimulus = seq(min(input.data[,1]),
                                             max(input.data[,1]),
                                             length.out=n.predictions))
  }
  return(pred.x.vals)
}


CheckForFlat <- function(input.data, response.column) {
  # Checks to see if all data points have the same y-values.
  if (sd(input.data[, response.column]) == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


Peak <- function(input.stimuli, preference.function, peak.within, is.flat) {
  # Finds the peak of a preference function.
  max.stim <- max(input.stimuli)
  min.stim <- min(input.stimuli)
  stim.range <- max.stim - min.stim
  if (length(peak.within) == 1) {
    end.caps <- ((1 - peak.within) / 2) * stim.range
    inner.max <- max.stim - end.caps
    inner.min <- min.stim + end.caps
  } else {
    inner.min <- min(peak.within)
    inner.max <- max(peak.within)
  }

  predicting.stimuli <- data.frame(stimulus = seq(min.stim, max.stim,
                                   length.out = 201))
  inner.max.index <- min(which(abs(predicting.stimuli - inner.max) ==
                         min(abs(predicting.stimuli - inner.max))))
  inner.min.index <- max(which(abs(predicting.stimuli - inner.min) ==
                         min(abs(predicting.stimuli - inner.min))))
  predicted.response1 <- predict.gam(preference.function, predicting.stimuli,
                                     se.fit = TRUE)

  if (is.flat == FALSE) {
    peak.response <- max(
      predicted.response1$fit[inner.min.index: inner.max.index])
    peak.response.index <- min(which(predicted.response1$fit == peak.response))

    if (peak.response.index == inner.min.index |
        peak.response.index == inner.max.index) {
      peak.response <- max(predicted.response1$fit)
      peak.response.index <- min(
        which(predicted.response1$fit == peak.response))
    }

    if (peak.response.index != 1 &
        peak.response.index != length(predicted.response1$fit)) {

      pred.stim2 <- data.frame(
        stimulus = seq(predicting.stimuli$stimulus[peak.response.index - 1],
                       predicting.stimuli$stimulus[peak.response.index + 1],
                       length.out = 201))
      pred.resp2 <- predict.gam(preference.function, pred.stim2)

      peak.response <- max(pred.resp2)
      peak.response.index <- min(which(pred.resp2 == peak.response))

      peak.preference <- pred.stim2[peak.response.index, ]
    } else {
      peak.preference <- predicting.stimuli$stimulus[peak.response.index]
    }
  }
  else {
    peak.preference <- NA
    peak.response <- mean(predicted.response1$fit)
    peak.response.index <- NA
  }

  peak.bundle <- list(peak.preference = peak.preference,
                      peak.response = peak.response,
                      peak.response.index = peak.response.index,
                      predicting.stimuli = predicting.stimuli,
                      predicted.response = predicted.response1$fit,
                      predicted.se = predicted.response1$se.fit,
                      max.stim = max.stim,
                      min.stim = min.stim)

  return(peak.bundle)
}


Tolerance <- function(drop, peak.bundle, is.flat, preference.function,
                      tol.floor) {
  # Finds the tolerance (the width of the curve at a given height) for a
  # preference function.
  submerged <- FALSE
  if(is.flat == TRUE){
    broad.tol = peak.bundle$max.stim - peak.bundle$min.stim
    strict.tol = broad.tol
    tolerance.height = mean(peak.bundle$predicted.response) * (1 - drop)
    cross.points = c(peak.bundle$min.stim, peak.bundle$max.stim)
    strict.lo = peak.bundle$min.stim
    strict.hi = peak.bundle$max.stim
  } else if(tol.floor >= peak.bundle$peak.response){
      submerged <- TRUE
      broad.tol <- 0
      strict.tol <- 0
      tolerance.height <- tol.floor
      cross.points <- vector()
      strict.lo <- NA
      strict.hi <- NA
  } else {
    predicting.stimuli.df <- peak.bundle$predicting.stimuli
    pred.stim <- predicting.stimuli.df$stimulus
    peak.pref <- peak.bundle$peak.preference
    predicted.response <- peak.bundle$predicted.response
    peak.response <- peak.bundle$peak.response
    tol.drop.amount <- (peak.response - tol.floor) * drop
    tolerance.height <- peak.response - tol.drop.amount
    shifted.points <- predicted.response - tolerance.height
    sign.shpt <- sign(shifted.points)
    delta.sign <- abs(diff(sign.shpt))
    cross.pt.ix <- which(delta.sign > 0)
    cross.points <- vector("numeric", 0)
    for (i in cross.pt.ix) {
      pred.stim2 <- data.frame(stimulus = seq(pred.stim[i], pred.stim[i + 1],
                               length.out = 101))
      pred.resp2 <- predict.gam(preference.function, pred.stim2)
      shifted.pts2 <- pred.resp2 - tolerance.height
      sign.shpt2 <- sign(shifted.pts2)
      for (j in 1:length(sign.shpt2)) {
        if (sign.shpt2[j] == 0 &
            !(pred.stim2$stimulus[j] %in% cross.points)) {
          cross.points <- append(cross.points, pred.stim2$stimulus[j])
        }
        if (!is.na(sign.shpt2[j + 1]) &
            sign.shpt2[j] == (sign.shpt2[j + 1] * -1) &
            !(pred.stim2$stimulus[j] %in% cross.points)) {
          cross.points <- append(cross.points, pred.stim2$stimulus[j])
        }
      }
    }
    if (sign.shpt[1] > -1 &
        !(pred.stim[1] %in% cross.points)) {
      cross.points <- append(cross.points, pred.stim[1])
    }
    if (sign.shpt[length(sign.shpt)] > -1 &
        !(pred.stim[length(sign.shpt)] %in% cross.points)) {
      cross.points <- append(cross.points, pred.stim[length(sign.shpt)])
    }
    cross.points <- sort(cross.points)

  # Calculate tolerance values
    stim.diffs <- cross.points - peak.pref
    if (peak.pref == max(cross.points)) {
      strict.hi <- peak.pref
    } else {
      strict.hi <- min(cross.points[which(stim.diffs > 0)])
    }
    if (peak.pref == min(cross.points)) {
      strict.lo <- peak.pref
    } else {
      strict.lo <- max(cross.points[which(stim.diffs < 0)])
    }
    strict.tol <- strict.hi - strict.lo
    if (length(cross.points) == 2) {
      broad.tol <- strict.tol
    } else {
      tol.bits <- vector('numeric', 0)
      for (i in seq(1, (length(cross.points) - 1), by = 2)) {
        tol.bits <- append(tol.bits, (cross.points[i + 1] - cross.points[i]))
        broad.tol <- sum(tol.bits)
      }
    }
  }
  tolerance.bundle <- list(broad.tolerance = broad.tol,
                           strict.tolerance = strict.tol,
                           tolerance.height = tolerance.height,
                           cross.points = cross.points,
                           strict.points = c(strict.lo, strict.hi),
                           submerged = submerged)
  return(tolerance.bundle)
}


StrResp <- function(input.d, predicted.r, is.flat, allfromsplines) {
  # Calculate strength and responsiveness of a preference function.
  if (allfromsplines == FALSE) {
    input <- input.d
  }
  else {
    input <- predicted.r
  }
  if (is.flat == TRUE){
    hd.strength <- 0
    hi.strength <- 0
  } else {
    hd.strength <- round((sd(input)/mean(input)) ^ 2, 3) #formerly cv strength
    hi.strength <- round(sd(input)/diff(range(input)), 3) #formerly Gray's str
  }
  responsiveness <- round (mean (input), 3)
  str.resp.bundle <- list(hd.strength = hd.strength,
                          hi.strength = hi.strength,
                          responsiveness = responsiveness)
  return(str.resp.bundle)
}


GraphSpline <- function (input.data, peak.bundle, tolerance.bundle,
                         response.column.name, response.column, graph.points,
                         graph.peak, graph.tol, graph.sp, tol.mode, max.y,
                         ghost, ghost.bundle, is.flat, smoothing.parameter,
                         orig.col.num, forgui=FALSE, group=FALSE,
                         graph.se=FALSE, graph.spline=TRUE) {
  peak.response <- peak.bundle$peak.response
  predicting.stimuli <- peak.bundle$predicting.stimuli
  predicted.response <- peak.bundle$predicted.response
  predicted.se <- peak.bundle$predicted.se
  peak.preference <- 0
  if (is.flat == FALSE) {
    peak.preference <- peak.bundle$peak.preference
  }

  tolerance.height <- tolerance.bundle$tolerance.height
  all.points <- tolerance.bundle$cross.points #change variable name?

  if (length(ghost.bundle) > 1) {
    ghost.peak.s <- ghost.bundle$d.peak.bundle$peak.preference
    ghost.peak.r <- ghost.bundle$d.peak.bundle$peak.response
    ghost.pred.r <- ghost.bundle$d.peak.bundle$predicted.response

    ghost.drop.r <- ghost.bundle$d.tolerance.bundle$tolerance.height
    ghost.all.p <- ghost.bundle$d.tolerance.bundle$cross.points
    ghost.strict.p <- ghost.bundle$d.tolerance.bundle$strict.points
  }

  if (!forgui) {
    plot(predicting.stimuli[, 1], predicted.response,
         main = paste(orig.col.num, " - ", response.column.name), xlab = "",
         ylab = "", ylim = c(0, max.y), type = "l", lwd = 1)
  } else if(forgui & graph.points & !group) {
    points(input.data[, 1], input.data[, response.column])
  } else if(forgui & graph.points & group) {
    n_constituents <- length(levels(as.factor(input.data$names)))
    for (i in 1:n_constituents) {
      current.subset.name <- levels(as.factor(input.data$names))[i]
      current.subset.rows <- which(input.data$names == current.subset.name)
      constituent.x.vals <- input.data[current.subset.rows, 2]
      constituent.y.vals <- input.data[current.subset.rows, 3]
      lines(constituent.x.vals, constituent.y.vals, col = gray(0.7))
    }
  }

  if (length(ghost.bundle) > 1) {
    lines(predicting.stimuli[, 1], ghost.pred.r, lwd = 1, col = gray(.8))
  }

  if (graph.sp == TRUE) {
    mtext(paste("sp = ", smoothing.parameter), 3, cex = .4, line = 0)
  }

  if (graph.peak == TRUE & is.flat == FALSE) {
    if (length(ghost.bundle) > 1) {
      lines(c(ghost.peak.s, ghost.peak.s), c(0, ghost.peak.r), col = "pink")
    }
    lines(c(peak.preference, peak.preference),
          c(0, peak.response), col = "red")
  }

  if (graph.tol == TRUE & tolerance.bundle$submerged == FALSE) {
    if (tol.mode == "strict") {
      if (length(ghost.bundle) > 1) {
        lines(ghost.strict.p,
              c(ghost.drop.r, ghost.drop.r), col = "light blue")
      }
      lines(strict.points,
            c(tolerance.height, tolerance.height), col = "blue")
    } else {
      if (length(ghost.bundle) > 1) {
        for (t in seq(1, (length(ghost.all.p) - 1), by = 2)) {
          lines(c(ghost.all.p[t], ghost.all.p[t + 1]),
                c(ghost.drop.r, ghost.drop.r), col = "light blue")
        }
      }
      for (t in seq(1, (length(all.points) - 1), by = 2)) {
        lines(c(all.points[t], all.points[t + 1]),
              c(tolerance.height, tolerance.height), col = "blue")
      }
    }
  }
  if (graph.se) {
    upper.se <- predicted.response + predicted.se
    lower.se <- predicted.response - predicted.se
    lines(predicting.stimuli[, 1], upper.se, lwd = 1, lty = 2, col = '#666666')
    lines(predicting.stimuli[, 1], lower.se, lwd = 1, lty = 2, col = '#666666')
  }
  if (graph.spline) {
    lines(predicting.stimuli[, 1], predicted.response,
          main = response.column.name, xlab = "",
          ylab = "", ylim = c(0, max.y), type = "l", lwd = 1)
  }
}


SPBinding <- function(smoothing.parameter, max.sp, min.sp) {
  # Restricts the smoothing parameter between two values: max.sp and min.sp
  if (smoothing.parameter > max.sp) {
    smoothing.parameter <- max.sp
  }
  if (smoothing.parameter < min.sp) {
    smoothing.parameter <- min.sp
  }
  return(smoothing.parameter)
}


SPAssign <- function(smoothing.parameter, orig.col.num, sp.assign){
  # Specify a matrix of two columns with target's column number in the first
  # column, and the target's desired smoothing parameter in the second. To
  # construct such a matrix, follow this example:
  # new.sp <- matrix(c(2, 1, 3, .5, 4, 0), ncol=2, byrow=T)
  # In the above example, the individual in input.data column 2 gets a sp of
  # 1, the individual in input.data column 3 gets an sp of .5. The one in
  # column 4 gets an sp of 0.
  if (orig.col.num %in% sp.assign[, 1]) {
    sp.assign.row <- which(sp.assign[, 1] == orig.col.num)
    if (sp.assign[sp.assign.row, 2] >= 0) {
      smoothing.parameter <- sp.assign[sp.assign.row, 2]
    }
  }
  return(smoothing.parameter)
}


Ghost <- function(input.data, response.column, stimulus, k,
                  preference.function, smoothing.parameter, orig.col.num,
                  sp.assign, input.stimuli, peak.within, drop, is.flat,
                  tol.floor) {
  # Plots the spline with the default smoothing parameter alongside the spline
  # with the current smoothing parameter in lighter colors.
  # Useful for comparing effects of changes in the SP.
  ghost.bundle <- list(NULL)
  if (orig.col.num %in% sp.assign[, 1]) {
    preference.function <- gam(input.data[, response.column] ~
                               s(stimulus, k = k), data = input.data,
                               scale = -1, sp = smoothing.parameter)
    default.pf <- preference.function
    default.sp <- smoothing.parameter
    d.peak.bundle <- Peak(input.stimuli, preference.function, peak.within,
                          is.flat)
    d.tolerance.bundle <- Tolerance(drop, d.peak.bundle, is.flat,
                                    preference.function, tol.floor)
    ghost.bundle <- list(default.pf = default.pf,
                         default.sp = default.sp,
                         d.peak.bundle = d.peak.bundle,
                         d.tolerance.bundle = d.tolerance.bundle
                         )
  }

    return(ghost.bundle)
}


Diagnose <- function(input.data, diagnose.col, diagnose.sp, peak.within, drop,
                     tol.mode, max.y, pred.x.vals, allfromsplines, k,
                     sp.binding, min.sp, max.sp, assign.sp, graph.points,
                     graph.se, forgui, tol.floor, points.out) {
  # A very useful function that allows users to view individual splines
  # without outputting any files. It is called by passing a number as the
  # second positional argument in PFunc() (the first being the name of the data
  # file). This number is the column number of the individual they want to
  # examine.
  # For example, to view the individual in column 3, call PFunc(mydata, 3).
  # To view that same spline with a smoothing parameter of 0.1, call
  # PFunc(mydata, 3, 0.1)
  names(input.data)[1] <- "stimulus"
  input.stimuli <- input.data[1]

  preference.function <- gam(input.data[, diagnose.col] ~
                             s(stimulus, k = k), data = input.data, scale = -1,
                             sp=diagnose.sp)
  pred.y.vals <- predict.gam(preference.function, pred.x.vals, se.fit = TRUE)
  predicted.points <- cbind(pred.x.vals, pred.y.vals)

  if (points.out != FALSE) {
    names(predicted.points)[2] <- names(input.data)[diagnose.col]
    names(predicted.points)[3] <- "se"
    write.csv(predicted.points, points.out, row.names = FALSE)
  }

  smoothing.parameter <- preference.function$sp

  if (diagnose.sp > 0) {
    smoothing.parameter <- diagnose.sp
  } else if (sp.binding == TRUE) {
    smoothing.parameter <- SPBinding(smoothing.parameter, max.sp, min.sp)
    preference.function <- gam(input.data[, diagnose.col] ~
                               s(stimulus, k = k), data = input.data,
                               scale = -1, sp=smoothing.parameter)
  }

  is.flat <- CheckForFlat(input.data, diagnose.col)

  peak.bundle <- Peak(input.stimuli, preference.function, peak.within, is.flat)

  tolerance.bundle <- Tolerance(drop, peak.bundle, is.flat,
                                preference.function, tol.floor)

  str.resp.bundle <- StrResp(input.data[, diagnose.col],
                             predicted.points[, 2],
                             is.flat, allfromsplines)

  peak.preference <- peak.bundle$peak.preference
  peak.response <- peak.bundle$peak.response
  predicting.stimuli <- peak.bundle$predicting.stimuli
  predicted.response <- peak.bundle$predicted.response
  predicted.se <- peak.bundle$predicted.se

  tolerance.height <- tolerance.bundle$tolerance.height
  all.points <- tolerance.bundle$cross.points
  strict.points <- tolerance.bundle$strict.points

  hd.strength <- str.resp.bundle$hd.strength
  hi.strength <- str.resp.bundle$hi.strength
  responsiveness <- str.resp.bundle$responsiveness

  if (forgui == FALSE) {
    plot(predicting.stimuli[, 1], predicted.response,
         main = paste(diagnose.col, " - ", names(input.data)[diagnose.col]),
         xlab = "", ylab = "", ylim = c(0, max.y), type = "l", lwd = 1)
    if (graph.points == TRUE){
      points(input.data[, 1], input.data[, diagnose.col])
    }

    lines(c(peak.preference, peak.preference),
          c(0, peak.response), col = "red")

    if (tol.mode == "strict") {
      lines(strict.points,
            c(tolerance.height, tolerance.height), col = "blue")
      final.tol <- tolerance.bundle$strict.tolerance #tolhere
    } else {
      for (t in seq(1, (length(all.points) - 1), by = 2)) {
        lines(c(all.points[t], all.points[t + 1]),
              c(tolerance.height, tolerance.height), col = "blue")
      }
      final.tol <- tolerance.bundle$broad.tolerance
    }

    return(c(paste("Peak Preference: ", round(peak.preference, 2), sep = ""),
      paste("Peak Height: ", round(peak.response, 2), sep = ""),
      paste("Tolerance: ", round(final.tol, 2), sep = ""),
      paste("HD Strength: ", round(hd.strength, 3), sep = ""),
      paste("HI Strength: ", round(hi.strength, 3), sep = ""),
      paste("Responsiveness: ", round(responsiveness, 3), sep = ""),
      paste("Smoothing Parameter: ", round(smoothing.parameter, 5), sep = "")))

  } else if (forgui == TRUE) {
    # if (tol.mode == 'strict') {
    #   final.tol <- tolerance.bundle$strict.tolerance
    #   tol.points <- tolerance.bundle$strict.points
    # } else if (tol.mode == 'broad' | tol.mode == 'broad') {
    #   final.tol <- tolerance.bundle$broad.tolerance
    #   tol.points <- tolerance.bundle$cross.points
    # }
    gui.bundle <- list(data.x = input.data[, 1],
                       data.y = input.data[, diagnose.col],
                       gam.object = preference.function,
                       stimulus = predicting.stimuli[, 1],
                       response = predicted.response,
                       se = predicted.se,
                       peak.preference = peak.preference,
                       peak.response = round(peak.response, 3),
                       broad.tol = tolerance.bundle$broad.tolerance,
                       strict.tol = tolerance.bundle$strict.tolerance,
                       broad.tol.points = tolerance.bundle$cross.points,
                       strict.tol.points = tolerance.bundle$strict.points,
                       tol.height = tolerance.bundle$tolerance.height,
                       hd.strength = hd.strength,
                       hi.strength = hi.strength,
                       responsiveness = responsiveness,
                       smoothing.parameter = smoothing.parameter,
                       is.flat = is.flat)
    return(gui.bundle)
  }
}


InCheck <- function (term, domain){
  # The sole purpose of this function is to circumvent issues with python
  # interpreting the "%" character.
  return(term %in% domain)
}


# Load confirmation:
print("PFunc version 1.0.3 (2022-08-21) successfully loaded.")
