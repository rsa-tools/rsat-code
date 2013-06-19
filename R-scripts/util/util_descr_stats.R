################################################################
##
## Descriptive statistics
##
## Jacques van Helden <jvanheld@ucmb.ulb.ac.be>
##
## Univariate descriptive statistics
##

source(paste(dir.util, 'palette.R', sep='/'))
## source('util_plots.R')

################################################################
##
## Moment calculation
##

moment <- function(x=rnorm(1000),k=1,c=0) {
  n <- length(x)
  m <- sum(as.numeric((x-c)^k)/n)
  return(m)
}

moment.central <- function(x=rnorm(1000),k=1) {
  m <- moment(x,1,0) # mean
  return(moment(x,k,m))
}

################################################################
##
## Summary of descriptive parameters
##
descr.stats <- function(vect=rnorm(1000),print=F) {
  n <- length(vect)
  m <- moment(vect,1,0)
  s2 <- moment(vect,2,m)
  s2.pop <- s2 * n/(n-1)
  s <- sqrt(s2)
  s.pop <- sqrt(s2.pop)
  g1 <- moment(vect,3,m)/sqrt(moment(vect,2)^3)
  b1 <- g1^2
  b2 <- moment(vect,4,m)/moment(vect,2,m)^2
  g2 <- b2 - 3
  
  description <- data.frame(
                            n = n,
                            
                            mean = m, # mean

                            median = median(vect), # median

                            s2 = s2, # sample variance
                            s2.pop = s2.pop, # population variance (inferred)

                            s = s, # sample standard deviation
                            s.pop = s.pop, # population standard deviation (inferred)

                            V = s/abs(m), # coefficient of variation
                            min = min(vect,na.rm=T), # minimum
                            max = max(vect,na.rm=T), # maximum
                            range = max(vect,na.rm=T)-min(vect,na.rm=T), # range

                            b1 = b1, # Pearson dissimetry coefficient
                            g1 = g1, # Fisher dissimetry coefficient

                            b2 = b2, # Pearson kurtosis coefficient
                            g2 = g2 # Fisher kurtosis coefficient
                            )
  if (print) {
    write.table(t(description),col.names=F,quote=F,sep="\t")
  } 
  return(description)
}

################################################################
##
## Class grouping.
##
## This function is obsolete, it is easier to use the standard R
## function hist()
class.grouping <- function(vect=rnorm(1000),
                           breaks=NA,
                           class.number=NA,
                           distrib.min=NA,
                           distrib.max=NA,
                           class.interval=NA,
                           print=FALSE,
                           return.type="descr", # supported: "descr","hist"
                           display.plot=FALSE,
                           col="#ffbb77"
                           ) {

  vect <- na.omit(vect)

  if (!is.na(breaks)) {
    distrib.min <- min(breaks)
    distrib.max <- max(breaks)
    class.interval <- breaks[2] - breaks[1]
    class.number <- length(breaks) - 1
  }
  
  if (is.na(distrib.min)) {
    if (is.na(class.interval) ||
        is.na(class.number) ||
        is.na(distrib.max)) {
      distrib.min <- min(vect,na.rm=T)
      ##  distrib.max <- max(vect,na.rm=T)
    } else {
      distrib.min <- distrib.max - class.interval*(class.number-1)
    }
  }
  if (is.na(class.number)) {
    if (is.na(class.interval) || 
        is.na(distrib.min) || 
        is.na(distrib.max)) {
      class.number <- 25
    } else {
      class.number <- as.integer((distrib.max - distrib.min)/class.interval) + 1
    }
  }

  if (is.na(class.interval)) {
    if (is.na(distrib.max)) {
      distrib.max <- max(vect,na.rm=T)
      ## this is to make sure we include the highest value in one interval
      class.interval <- (distrib.max - distrib.min)/(class.number -1) 
    } else {
      class.interval <- (distrib.max - distrib.min)/class.number
    }
  }
  ## class.breaks <- distrib.min + class.interval*0:class.number
  class.breaks <- pretty(c(distrib.min,distrib.max),n=class.number)

  my.hist <- hist(vect,breaks=class.breaks,col=col,plot=display.plot)

  c   <- length(my.hist$mids)
  class.index <- 1:c
  class.min <- my.hist$breaks[1:c]
  class.max <- my.hist$breaks[2:(c+1)]
  class.mid <- my.hist$mid
  class.intervals <- my.hist$breaks[2:(c+1)] - my.hist$breaks[1:c]
  class.occ <- my.hist$counts
  class.occ.cum <- cumsum(class.occ)
  class.freq <- my.hist$counts/sum(my.hist$counts)
  class.freq.cum <- cumsum(class.freq)
  class.intensity <- class.freq/class.intervals
  
  class.descr <- data.frame(
                            "min" = class.min,
                            mid = class.mid,
                            "max" = class.max,
                            occ = class.occ,
                            occ.cum = class.occ.cum,
                            freq = class.freq,
                            freq.cum = class.freq.cum,
                            intensity = class.intensity
                            )

  distrib.summary <- data.frame(
                                p = class.number,
                                interv = class.interval,
                                min = distrib.min,
                                max = distrib.max
                                )

  if (return.type == "descr") {
    to.return <- class.descr
  } else if (return.type == "hist") {
    to.return <- my.hist
  } else if (return.type == "summary") {
    to.return <- distrib.summary
  } else {
    to.return <- NULL
  }


  if (print) {
    print(class.descr, quote=FALSE,row.names=FALSE)
  } else {
    return(to.return)
  }
}

