## ##############################################################
## Calculate the negative binomial, with the recursive formula given
## by Dagnelie
## source('~/rsa-tools/R-scripts/negbin.R')

## Density function
dnegbin <- function (x, p=NA, k=NA, m=NA, v=NA, log=FALSE, ...) {
    return(negbin(x, p, k, m, v, log, cum=FALSE, ...))
}

## Probability function
pnegbin <- function (q, p=NA, k=NA, m=NA, v=NA, ...) {
    return(negbin(q, p, k, m, v, log=FALSE, cum=TRUE, ...))
}

## generic function called by the other ones
negbin <- function (x, p=NA, k=NA, m=NA, v=NA, log=FALSE, cum=FALSE, verbose=F) {
    x.max <- max(x)

    if ((is.na(p)) && (is.na(k)) && (!is.na(m)) && (!is.na(v))) {
	p <- v/m -1
	k <- m/p
    } else if ((!is.na(p)) && (!is.na(k)) && (is.na(m)) && (is.na(v))) {
	m <- k*p
	v <- k*p*(1+p)
    } else {
	stop('You sould provide either p and k, or m and v')
    }
    q <- 1+p

    if (verbose) {
	print (paste("p", p))
	print (paste("k", k))
	print (paste("m", m))
	print (paste("v", v))
    }

    p.log <- log(p)
    q.log <- log(q)
    k.log <- log(k)

    ## value for x=0
    negbin.log <- -k*q.log
    negbin.series.log <- negbin.log
	
    ## recurse
    for (i in 1:x.max) {
	negbin.log <- negbin.log + p.log + log(k+i-1) - log(i) - q.log
	negbin.series.log <-  append(negbin.series.log, negbin.log)
    }

    if (cum) {
	return(cumsum(exp(negbin.series.log))[x+1])
    } else if (log) {
	return (negbin.series.log[x+1])
    } else {
	return (exp(negbin.series.log[x+1]))
    }
}


