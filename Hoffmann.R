# R code for generating Hoffmann plots
# Code authors: Dr. Daniel T. Holmes, MD and Dr. Kevin A Buhr, PhD
# ----------------------------------------
# Usage
# ----------------------------------------
# Hoff.QQ.plot <- function(x, alpha = 0.05, from = NA, to = NA, xrange = c(NA,NA))
# ----------------------------------------
# Arguments
# ----------------------------------------
# x         a vector of observations drawn from the mixure for which
#           decomposition is desired.
# from      a vector of length 1 to define the lower limit of a linear
#           section of interest in the QQ-plot of x. Both 'from' and
#           'to' must be supplied values to obtain a reference
#           interval estimate.
# to        a vector of length 1 to define the upper limit of a linear
#           section of interest in the QQ-plot of x.Both 'from' and
#           'to' must be supplied values to obtain a reference
#           interval estimate.
# alpha     alpha value reference interval estimate. 
# xlim      vector of length = 2 defining a plot range.
# ----------------------------------------
# Value
# ----------------------------------------
# Returns a named vector with the following values:
# lower:    the estimate of the lower limit of normal
# upper:    the estimate of the upper limit of normal

Hoff.QQ.plot <- function(x, alpha=0.05, from=NA, to=NA, xlim=range(x)) {
  x <- sort(x)
  qq.data <- as.data.frame(qqnorm(x, datax = TRUE, plot.it = FALSE))
  plot(y ~ x, 
       data = qq.data,
       type = "l",
       xlab = "Result",
       ylab = "Quantiles of the Normal Distibution",
       xlim = xlim)
  if (!is.na(from) & !is.na(to)) {
    linear <- subset(qq.data, x >= from & x <= to)
    lin.mod <- lm(y ~ x, data = linear)
    abline(lin.mod)
    RI <- (c(qnorm(alpha/2),qnorm(1-alpha/2))
      - lin.mod$coefficients[1])/lin.mod$coefficients[2]
    result <- c(RI[1], RI[2])
    names(result) <- c("lower", "upper")
    abline(h = c(qnorm(alpha/2), qnorm(1 - alpha/2)), col = "blue")
    abline(v = RI, col = "blue", lty = 2)
    return(result)
  }
}

# Example
# ----------------------------------------
# Built-in R data set for waiting times between eruptions of the Old
# Faithful geyser in Yellowstone National Park.

data("faithful")
hist(faithful$waiting)
Hoff.QQ.plot(faithful$waiting)
# linear sections are seen from ~ 48 to 53 min and ~ 77 to 89
Hoff.QQ.plot(faithful$waiting, from = 48, to = 54)
Hoff.QQ.plot(faithful$waiting, from = 77, to = 89)
