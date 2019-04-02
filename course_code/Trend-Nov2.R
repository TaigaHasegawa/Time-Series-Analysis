###########################
#generate trend + noise 
###########################
a0<-1
a1<--0.01
a2<-0.0002

set.seed(100)
n<-500
data<-a0+a1*(1:n)+a2*((1:n)^2)+rnorm(n,0,1)
plot(data,type="l")

#remove the trend by fitting polynomials
#assuming the order of polynomials are known

X1<-1:n
X2<-(1:n)^2
fit<-lm(data~1+X1+X2)
residual<-fit$resid
hata0<-as.numeric(fit$coef[1])
hata1<-as.numeric(fit$coef[2])
hata2<-as.numeric(fit$coef[3])
par(mfrow=c(2,2))
plot(data,type="l",main="original data")
plot(hata0+hata1*X1+hata2*X2,main="fitted trend")
plot(residual,type="l",main="residual series")
acf(residual,lag.max=50,main="ACF of residual series")

#########################################
#apply smoothing to estimate the trend
#########################################
#q<-10
q<-10
estitrend<-rep(0,n)
for (j in (q+1):(n-q))
{
estitrend[j]<-mean(data[(j-q):(j+q)])
}
for (j in 1:q)
{
estitrend[j]<-mean(c(data[1:(j+q)],rep(data[1],q+1-j)))
}
for (j in (n-q+1):n)
{
estitrend[j]<-mean(c(data[(j-q):n],rep(data[n],2*q+1-(n-j+q+1))))
}
par(mfrow=c(2,2))

plot(data,type="l",main="original data")

plot(estitrend,type="l",main="Estimated trend using smoothing")

############
#plot the residual and its acf
##################
plot(data-estitrend,type="l")
acf(data-estitrend)

#############################################
#apply differencing to eliminate the trend
#####################
diffdata<-data[2:n]-data[1:(n-1)]
diff2data<-diffdata[2:(n-1)]-diffdata[1:(n-2)]
par(mfrow=c(2,3))
plot(data,type="l",main="original data")
plot(diffdata,type="l",main="1st differenced data")
plot(diff2data,type="l",main="2nd differenced data")
acf(data,40)
acf(diffdata,40)
acf(diff2data,40)

#################################################
#use the data from the book as illustrations
#################################################
#The weekly cardiovascular mortality series (in Los Angeles County)
#the strong seasonal components in all of the series, corresponding to
#winter-summer varations and the downward trend in the cardiovascular mortality
#over the 10-year period. (see Example 2.2. in the book for more details)

#Example 2.10
pdf(file="mortMA.pdf")
mort=scan("cmort.dat")
 t = 1:length(mort)
 ma5 = filter(mort, sides=2, rep(1,5)/5)
 ma53 = filter(mort, sides=2, rep(1,53)/53)
 plot(t, mort, xlab="week", ylab="mortality")
 lines(ma5)
 lines(ma53)
dev.off()

Example 2.11 

pdf(file="mortPara.pdf")
 t = 1:length(mort)
 t2 = t^2
 t3 = t^3
 c = cos(2*pi*t/52)
 s = sin(2*pi*t/52)
 fit1 = lm(mort~t + t2 + t3)
 fit2 = lm(mort~t + t2 + t3 + c + s)
 plot(t, mort)
 lines(fit1$fit)
 lines(fit2$fit)
par(mfrow=c(2,1))
plot(fit2$residual)
acf(fit2$residual)

dev.off()


Example 2.12 
pdf(file="mortSmooth.pdf")
 plot(t, mort)
 lines(ksmooth(t, mort, "normal", bandwidth=10))
 lines(ksmooth(t, mort, "normal", bandwidth=104))
dev.off()

Example 2.14
 plot(t, mort)
 lines(smooth.spline(t, mort, spar=.0000001))
 lines(smooth.spline(t, mort, spar=.1))

############################################
#Paleoclimatic Glacial Varves (Example 2.5)
#Valves (sedimentary deposits) is used as proxies
#for paleoclimatic parameters, such as temperature, because
#warm years corresponds to more sand and silt , and thick varve.
#The data contains the thicknesses of the yearly varves collected from
#one location in Massachusetts for 634 years, beginning 11,834 years ago
#(Nonstationarity in variance, long memory feature)
##################################
pdf(file="varve.pdf")
varve<-scan("varve.dat")
t = 1:length(varve)
par(mfrow=c(2,1))
plot(t,varve,type="l")
plot(t,log(varve),type="l")
dev.off()
#normality improved
par(mfrow=c(2,2))
hist(varve)
hist(log(varve))
qqnorm(varve)
qqnorm(log(varve))

##
#long memory?
###############
acf(log(varve),lag.max=100)


#######################
#######################


