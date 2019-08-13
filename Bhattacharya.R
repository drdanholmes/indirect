# R code for generating Bhattacharya plots
# Code author: Dr. Daniel T. Holmes, MD
# ----------------------------------------
# Usage
# ----------------------------------------
# bhatt.plot(x, breaks = "Sturges", from = NA, to = NA,
# alpha = 0.05, xlim = range(x)
# ----------------------------------------
# Arguments
# ----------------------------------------
# x         a vector of observations drawn from the mixure for which decomposition is desired.
# breaks    generally a vector of breakpoints to define bins for Bhattacharya plot. 
#           Enter ?hist for other options. Defaults to same algorithm as hist function.
#           Breaks must not be so numerous that empty bins are created. Breaks must be evenly spaced.
# from      optional vector to define the starting bin number(s) for the regression of the 
#           linear section(s) of interest in the Bhattacharya plot. Must be of equal length to
#           to vector 'to'.
# to        optional vector to define the ending bin number(s) for the regression of the 
#           linear section(s) of interest in the Bhattacharya plot. Must be of equal length to
#           vector 'from'.
# alpha     alpha value reference interval estimate. 
# xlim      vector of length = 2 defining a plot range.
# plot.num  logical; if TRUE bin numbers are plotted beside points. Facilitates selection of
#           variables from and to.
# ----------------------------------------
# Value
# ----------------------------------------
# Returns dataframe with the following columns:
# mu:       the means of the fitted distribution(s)
# sigma:    the standard deviation(s) of the fitted distributon(s)
# p:        the proportion(s) of the fitted distributions
# lower:    mu - 1.96*sigma
# upper:    mu + 1.96*sigma

# textxy taken from the calibrate package
textxy <- function (X, Y, labs, m = c(0, 0), cex = 0.5, offset = 0.8, ...) 
{
  posXposY <- ((X >= m[1]) & ((Y >= m[2])))
  posXnegY <- ((X >= m[1]) & ((Y < m[2])))
  negXposY <- ((X < m[1]) & ((Y >= m[2])))
  negXnegY <- ((X < m[1]) & ((Y < m[2])))
  if (sum(posXposY) > 0) 
    text(X[posXposY], Y[posXposY], labs[posXposY],
         adj = c(0.5 - offset, 0.5 - offset), cex = cex, ...)
  if (sum(posXnegY) > 0) 
    text(X[posXnegY], Y[posXnegY], labs[posXnegY],
         adj = c(0.5 - offset, 0.5 + offset), cex = cex, ...)
  if (sum(negXposY) > 0) 
    text(X[negXposY], Y[negXposY], labs[negXposY],
         adj = c(0.5 + offset, 0.5 - offset), cex = cex, ...)
  if (sum(negXnegY) > 0) 
    text(X[negXnegY], Y[negXnegY], labs[negXnegY],
         adj = c(0.5 + offset, 0.5 + offset), cex = cex, ...)
}

bhatt.plot <- function(x, breaks = "Sturges", from = NA, to = NA, alpha = 0.05, xlim = range(x), plot.num = TRUE){
  if (is.numeric(breaks) && length(breaks) > 1) {
    # if explicit breaks, limit x to range covered by breaks
    x <- x[x >= min(breaks) & x <= max(breaks)]
  }
  dhist <- hist(x,
                breaks = breaks,
                plot = FALSE)
  if (any(dhist$counts==0)) {
    warning("some bin counts were zero: setting them to 0.1")
    dhist$counts[dhist$counts==0] <- 0.1
  }
  ly <- log(dhist$counts)
  dly <- diff(ly)
  df <- data.frame(xm = dhist$mids[-length(dhist$mids)],
                   ly = dly,
                   counts = dhist$counts[-length(dhist$mids)])
  h <- diff(df$xm)[1]
  plot(ly ~ xm,
       data = df,
       xlab = "Bin Midpoint, xm",
       ylab = expression(paste(Delta, log(y))),
       xlim = xlim)
  abline(h = 0)
  
  #number all the points if requested
  if(plot.num == TRUE){
    textxy(df$xm, df$ly, row.names(df),
           row.names(df),
           cex = 0.8,
           offset = 1,
           col = "blue")
  }
  
  if(!all(is.na(from)) & !all(is.na(to)) & length(from) == length(to)){
    mu <- rep(NA,length(from))
    sigma <- rep(NA,length(from))
    N <- rep(NA,length(from))
    p <- rep(NA,length(from))
    
    for(i in 1:length(from)){
      linear <- subset(df[from[i]:to[i],])
      lin.mod <- lm(ly ~ xm, data = linear, weights = linear$counts)
      abline(lin.mod)
      lambda <- -coef(lin.mod)[1]/coef(lin.mod)[2]
      mu[i] <- lambda + h/2
      sigma[i] <- sqrt(-h/coef(lin.mod)[2] - h^2/12)
      P <- pnorm((df$xm[from[i]:to[i]] + h/2 - mu[i])/sigma[i]) - 
        pnorm((df$xm[from[i]:to[i]] - h/2 - mu[i])/sigma[i])
      N[i] <- sum(df$counts[from[i]:to[i]])/sum(P)
    }
    p <- N/sum(N)
    result <- data.frame(mu,sigma,p, mu + qnorm(alpha/2)*sigma, mu + qnorm(1 - alpha/2)*sigma)
    names(result) <- c("mu","sigma","p","lower", "upper")
    return(result)
  }
}

# Example 1
# ----------------------------------------
# Data is taken from Bhattachary's original paper (Bhattacharya CG. A simple method of 
# resolution of a distribution into Gaussian components. Biometrics. 1967;1:115-35.)
# which in turn came from an earlier paper looking at the fork length of Porgy fish caught 
# in the East China sea: (Tanaka, S. A method of analysing of polymodal frequency 
# distribution and its application to the length distribution of the Porgy, Taius 
# tumifrons (J. and S.). J. Fish. Res. Bd. Can. 1962;19:1143-59.)

mids <- seq(from = 9.5, to = 29.5, by = 1)
counts <- c(509,2240,2341,623,476,1230,1439,921,448,512,
            719,673,445,341,310,228,168,140,114,64,22)
ex1.df <- data.frame(mids, counts)
# Reconstruct Bhattacharya's data to emulate the structure of raw data
ex1.data <- unlist(mapply(rep, ex1.df$mids, ex1.df$counts))
# Plot Bhattachary plot with no linear sections defined
bhatt.plot(ex1.data, breaks = 9:30)
# Once the Bhattachary plot is displayed the 'from' and 'to' vectors represent the four 
# obviously linear sections: points 1-3, points 5-7, points 10-12 and points 18-20
from <- c(1,5,10,18)
to <- c(3,7,12,20)
bhatt.plot(ex1.data, breaks = 9:30, from = from, to = to)

# Compare w Bhattacharya's graphically derived results of:
# mu       sigma    p
# 11.03    0.81     0.4065
# 15.28    1.13     0.3067
# 19.86    1.60     0.2087
# 26.62    1.47     0.0361

# Example 2
# ----------------------------------------
# Data is taken from Bhattacharya's paper as per example 1.

mids <- seq(from = 8.5, to = 25.5, by = 1)
counts <- c(31,532,2198,2297,685,494,1188,1479,938,486,537,702,664,431,192,59,12,2)
ex2.df <- data.frame(mids, counts)
ex2.data <- unlist(mapply(rep, ex2.df$mids, ex2.df$counts))
bhatt.plot(ex2.data, breaks = 8:26, from = c(1,6,11), to = c(4,8,17))

# Compare w Bhattacharya's graphically derived results of:
# mu        sigma       p
# 11.04     0.81        0.4311
# 15.31     1.16        0.3444
# 19.85     1.60        0.2245

# Example 3
# ----------------------------------------
# Built-in R data set for waiting times between eruptions of the Old Faithful geyser in 
# Yellowstone National Park.

data("faithful")
hist(faithful$waiting)
bhatt.plot(faithful$waiting)
fit <- bhatt.plot(faithful$waiting, from = c(1,7), to = c(3,10))
fit
# visualize resulting fit
hist(faithful$waiting,
     freq = FALSE,
     main = "Waiting times for Old Faithful",
     xlab = "time(min)")
curve(fit$p[1]*dnorm(x, mean = fit$mu[1], sd = fit$sigma[1]),
      from  = 0,
      to = 100,
      lwd = 2,
      col = "red",
      add = TRUE)
curve(fit$p[2]*dnorm(x, mean = fit$mu[2], sd = fit$sigma[2]),
      from  = 0,
      to = 100,
      lwd = 2,
      col = "green",
      add = TRUE)


