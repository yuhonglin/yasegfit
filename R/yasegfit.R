##' Segfit a series
##'
##' Segfit a series
##' @title Segfit a series
##' @param data The input series to be segfitted.
##' @param seg.alg Segmentation algorithm, possible values: c("topdown", "dp")
##' @param fit.alg Fitting algorithm for each segment,
##'                possible values: c("abline", "powerabc", "const")
##' @param smp Controlling the number of segments when seg.alg="fit"
##' @param nseg The exact number of segments when seg.alg="topdown"
##' @return A yasegfit object containing various kinds of information.
##'         For each segment, it contains the head index (hi), tail index (ei),
##'         parameter $a$ (a) parameter $b$ (b), parameter $c$ (c), fitting order
##'         (order), fitted values (fit), and fitting residuals (residual).
##'         It also has a attribute "data" for the original data; A attribute
##'         "seg.alg" for dp algorithm used; A attribute "fit.alg" for fitting
##'         algorithm used for each segment; A attribute "smp" for the smp used
##'         in "dp" algorithm; A attribute "nseg" for number of segments.
##' @export
yasegfit <- function(data, seg.alg="dp", fit.alg="abline", smp=2.3, nseg=3) {
    p = .Call("yasegfitc", s=as.numeric(data), seg.alg, fit.alg, smp=as.numeric(smp))

    ## compute fitted value and residuals
    ret = list()
    ratio = max(data)/100

    for (i in 1:ncol(p)) {

        fitted   = rep(NA, p[2,i] - p[1,i] + 1)
        residual = rep(NA, p[2,i] - p[1,i] + 1)

        for (j in 1:(p[2,i] - p[1,i] + 1)) {
            if (fit.alg=="powerabc")
                fitted[j] = (p[3,i] * j^p[4,i] + p[5,i])*ratio
            if (fit.alg=="abline")
                fitted[j] = (p[3,i] * j + p[4,i])*ratio
            if (fit.alg=="const")
                fitted[j] = p[3,i]*ratio
        }


        if (p[6,i] == 1) {
            fitted = rev(fitted)
        }

        for (j in 1:(p[2,i] - p[1,i] + 1)) {
            if (i==1) {
                residual[j]  = data[j] - fitted[j]
            } else {
                residual[j]  = data[p[2,i-1]+j] - fitted[j]
            }
        }

        if (fit.alg=="powerabc")
            b = p[4,i]
        if (fit.alg=="abline")
            b = p[4,i]*ratio
        if (fit.alg=="const")
            b = NaN
        
        ret[[i]] = list(
            hi = p[1,i],
            ti = p[2,i],
            a = p[3,i]*ratio,
            b = b,
            c = p[5,i]*ratio,
            order = p[6,i],
            fit = fitted,
            res = residual
        )
    }

    ## also append the original sequence, use attr to
    ## make "length(ret) == #segments"
    attr(ret, "data") <- data
    ## other attributes
    attr(ret, "smp")  <- smp
    attr(ret, "seg.alg")  <- seg.alg    
    attr(ret, "fit.alg")  <- fit.alg
    
    attr(ret, "nseg") <- length(ret)

    class(ret) <- "yasegfit"
    return (ret);
}

##' Plot a yasegfit object
##'
##' Plot a yasegfit object
##' @title Plot a segfit object
##' @param x If x is a segfit object, input y is ignored.
##' @param y A yasegfit object, If not provided, x must
##'        be a yasegfit object.
##' @param col.data color of data
##' @param col.seg  color of segments
##' @param legend.pos legend positions
##' @export 
plot.yasegfit <- function(x, y='', col.data='black', col.seg='red', legend.pos="topleft", ...) {
    if (y == '') {
        y = x;
        x = 1:length(attr(y,"data"))
    }

    plot.default(x, attr(y, "data"), type="l", col=col.data, ...)
    
    for (i in 1:length(y)) {
        lines(x[y[[i]]$hi:y[[i]]$ti], y[[i]]$fit, col=col.seg)
    }

    legend(x=legend.pos, legend=c("data", "segments"), lty=c(1,1), col=c('black', 'red'))
}

##' Summarise a yasegfit object
##'
##' Summarise a yasegfit object
##' @title Summarise a yasegfit object
##' @param sf The yasegfit object to summarise
##' @export 
summary.yasegfit <- function(sf) {
    cat(sprintf('Number of segments: %d\n', length(sf)))
    cat(sprintf('smp    : %f\n', attr(sf, "smp")))
    cat(sprintf('seg.alg: %s\n', attr(sf, "seg.alg")))
    cat(sprintf('fit.alg: %s\n', attr(sf, "fit.alg")))        

    cat(sprintf('\n% 8s% 8s% 8s% 8s% 8s% 8s% 8s% 8s\n',
                'Segment', 'hi', 'ti', 'a', 'b', 'c', 'mse', 'var.res'))

    for (i in 1:length(sf)) {
        cat(sprintf('% 8d% 8d% 8d%8.3f%8.3f%8.3f%8.3f%8.3f\n',
                      i, sf[[i]]$hi, sf[[i]]$ti, sf[[i]]$a, sf[[i]]$b,
                      sf[[i]]$c, mean(sf[[i]]$res^2), var(sf[[i]]$res)))
    }
}

##' Extract the residuals from a yasegfit object
##'
##' Extract the residuals from a yasegfit object
##' @title Extract the residuals from a yasegfit object
##' @param sf The yasegfit object
##' @param concat Whether the residuals of different segments will be
##'               concatenated.
##' @return If concat==TRUE, the residuals of different segments
##'               will be concatenated into one vector. Otherwise,
##'               the returned value will be a list of residuals
##'               of each segment.
##' @export
residuals.yasegfit <- function(sf, concat=TRUE) {
    if (concat==TRUE) {
        ret = rep(NA, length(attr(sf, "data")))
        stidx = 1
        for (i in 1:length(sf)) {
            ret[stidx:(stidx + length(sf[[i]]$res) - 1)] = sf[[i]]$res
            stidx = stidx + length(sf[[i]]$res)
        }
        return (ret)
    } else {
        return (lapply(sf, function(x) {x$res}))
    }
}

##' Extract the fitted values from a yasegfit object
##'
##' Extract the fitted values from a yasegfit object
##' @title Extract the fitted values from a yasegfit object
##' @param sf The yasegfit object
##' @param concat Whether the fitted values of different segments will be
##'               concatenated.
##' @return If concat==TRUE, the fitted values of different segments
##'               will be concatenated into one vector. Otherwise,
##'               the returned value will be a list of fitted values
##'               of each segment.
##' @export
fitted.yasegfit <- function(sf, concat=TRUE) {
    if (concat==TRUE) {
        ret = rep(NA, length(attr(sf, "data")))
        stidx = 1
        for (i in 1:length(sf)) {
            ret[stidx:(stidx + length(sf[[i]]$fit) - 1)] = sf[[i]]$fit
            stidx = stidx + length(sf[[i]]$fit)
        }
        return (ret)
    } else {
        return (lapply(sf, function(x) {x$fit}))
    }    
}
