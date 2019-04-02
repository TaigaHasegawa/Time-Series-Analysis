##################################################################
#tsdiag2 returns the p-values for the 
#Ljung-Box test statistic for a range of lags
#H=p+q+1,...p+q+gof.lag, where p and q are the ARMA model orders.
#Note that your time series must be longer than gof.lag+p+q+1.
##################################################################

tsdiag2<-function (object,gof.lag = 10, ...)
{
oldpar <- par(mfrow = c(3, 1))
on.exit(par(oldpar))
rs <- object$residuals

p <- as.list(object$call$order)[[2]][1]
q <- as.list(object$call$order)[[4]][1]

stdres <- rs/sqrt(object$sigma2)
plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
abline(h = 0)
acf(object$residuals, plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
nlag <- gof.lag+p+q
pval <- c()
s<-(1+p+q)
for (i in s:nlag) pval[i] <- Box.test2(rs, i, type = "Ljung-Box",fitdf=p+q)$p.value
plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")
}


############################################
#when there is a seasonal component,       #
#we use tsdiagseasonal2                    #
############################################
tsdiagseasonal2<-function (object,gof.lag = 10, ...)
{
oldpar <- par(mfrow = c(3, 1))
on.exit(par(oldpar))
rs <- object$residuals

p <- as.list(object$call$order)[[2]][1]
q <- as.list(object$call$order)[[4]][1]

P<-as.list(object$call$seasonal$order)[[2]][1]
Q<-as.list(object$call$seasonal$order)[[4]][1]

stdres <- rs/sqrt(object$sigma2)
plot(stdres, type = "h", main = "Standardized Residuals", ylab = "")
abline(h = 0)
acf(object$residuals, plot = TRUE, main = "ACF of Residuals", na.action = na.pass)
nlag <- gof.lag+p+q+P+Q
pval <- c()
s<-(1+p+q+P+Q)
for (i in s:nlag) pval[i] <- Box.test2(rs, i, type = "Ljung-Box",fitdf=p+q+P+Q)$p.value
plot(1:nlag, pval, xlab = "lag", ylab = "p value", ylim = c(0, 1), main = "p values for Ljung-Box statistic")
abline(h = 0.05, lty = 2, col = "blue")
}

Box.test2<-function (x, lag = 1, type = c("Box-Pierce", "Ljung-Box"), fitdf = 0) 
{
    if (NCOL(x) > 1) 
        stop("x is not a vector or univariate time series")
    DNAME <- deparse(substitute(x))
    type <- match.arg(type)
    cor <- acf(x, lag.max = lag, plot = FALSE, na.action = na.pass)
    n <- sum(!is.na(x))
    PARAMETER <- lag - fitdf
    obs <- cor$acf[2:(lag + 1)]
    if (type == "Box-Pierce") {
        METHOD <- "Box-Pierce test"
        STATISTIC <- n * sum(obs^2)
        PVAL <- 1 - pchisq(STATISTIC, lag - fitdf)
    }
    else {
        METHOD <- "Box-Ljung test"
        STATISTIC <- n * (n + 2) * sum(1/seq.int(n - 1, n - lag) * 
            obs^2)
        PVAL <- 1 - pchisq(STATISTIC, lag - fitdf) # fitdf should be p+q instead of the default 0 in the original code
    }
    names(STATISTIC) <- "X-squared"
    names(PARAMETER) <- "df"
    structure(list(statistic = STATISTIC, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME), 
        class = "htest")
}


#Example: fish data
rec = scan("recruit.dat")
#Exact Maximum likelihood
rec.ML<-arima(rec,order=c(2,0,0),method="ML")
rec.ML
tsdiag2(rec.ML,10)

