## armaspec.s

armaspec=function(T=5000,arcoeff,macoeff,vare){

## Input
## T = resolution in frequency (length of time series)
## arcoeff = vector of coefficients of the AR component of the time series
## macoeff = vector of coefficients of the MA component
## vare = variance of the white noise error

## Output 
## specout = spectrum of the ARMA process defined from (0, pi]

## Note: This is the parameterization of the ARMA(p,q) process 
## that I will adopt in this program 
## X(t) = arcoeff[1]X(t-1) + ... + arcoeff[p]X(t-p) + 
## W(t) + macoeff[1]W(t-1) + ... + macoeff[q]W(t-q)
## I will assume that the process is already causal and 
## invertible, i.e., the user has done these checks before
## attempting to derive the spectrum.

freq = 2*pi*seq(0, T-1)/T;

arcoeff = c(1, -arcoeff);
macoeff = c(1, macoeff);

theta=c(1:T); phi=c(1:T);

for (k in 1:(T)){

theta.temp=0;
for (m in 1: (length(macoeff)))
{
theta.temp = theta.temp + macoeff[m]*complex(mod=1, arg=-1*m*freq[k])}
theta[k] = (abs(theta.temp))^2

phi.temp=0;
for (m in 1: (length(arcoeff))){
phi.temp = phi.temp + arcoeff[m]*complex(mod=1, arg=-1*m*freq[k])}
phi[k] = (abs(phi.temp))^2

}

spec = vare*(theta)/(phi)
specout = spec[2:(T/2 + 1)]
return(specout)
}

######################################### 
#look at the spectrum for ARMA models 
#########################################
#AR(1), phi=0.9
T = 1024;
par(mfrow=c(2,2))
spec1 = armaspec(T, arcoeff=c(0.9), macoeff=c(0), vare=1);
plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="AR1=0.9")
plot(arima.sim(n=T,list(ar=c(0.9))),type="l",main="AR1=0.9")

#AR(1), phi=-0.9
spec11 = armaspec(T, arcoeff=c(-0.9), macoeff=c(0), vare=1);
plot((1:(T/2))/T,spec11, xlab="",ylab="spectrum", type="l",main="AR1=-0.9")
plot(arima.sim(n=T,list(ar=c(-0.9))),type="l",main="AR1=-0.9")

#MA(1), theta=0.9
T = 1024;
par(mfrow=c(2,2))
spec1 = armaspec(T, arcoeff=c(0), macoeff=c(0.9), vare=1);
plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="MA1=0.9")
plot(arima.sim(n=T,list(ma=c(0.9))),type="l",main="MA1=0.9")

#MA(1), theta=-0.9
T = 1024;
spec1 = armaspec(T, arcoeff=c(0), macoeff=c(-0.9), vare=1);
plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="MA1=-0.9")
plot(arima.sim(n=T,list(ma=c(-0.9))),type="l",main="MA1=-0.9")

#AR(2), phi_1=1, phi_2=-0.9
#
T = 1024;
par(mfrow=c(2,2))
spec1 = armaspec(T, arcoeff=c(1,-0.9), macoeff=c(0), vare=1);
plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="AR2(1,-0.9)")
#observe a peak at w=0.16, corresponding to 0.16 cycles per point/or 6/7 points per cycle
set.seed(100)
simudata<-arima.sim(n=T,list(ar=c(1,-0.9)))
plot(simudata,type="l",main="AR2(1,-0.9)")
#roughly 8 cycles for the first  50 points
plot(simudata[1:50],type="l",main="AR2(1,-0.9)")
#roughly 16 cycles for the first 100 points
plot(simudata[1:100],type="l",main="AR2(1,-0.9)")
###################
#ACF for AR(1), complex roots, 
#rho(h)=a|z_1|^{-h}\cos(h\theta+b)
#where z_1 is one of the roots,


#############################################################
#average: low pass filter
#differencing: high pass filter
#as applied to white noise series (for demonstration purpose)
#############################################################
#average
T = 1024;
freq<-seq(1, T/2)/T;
par(mfrow=c(2,2))
spec1 = armaspec(T, arcoeff=c(0), macoeff=c(0), vare=1);
spec2 = (1+2*cos(2*pi*freq))^2/9


plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="White noise")
plot((1:(T/2))/T,spec2, xlab="",ylab="spectrum",type="l",main="Averaged white noise")
simudata<-arima.sim(n=T,list(order=c(0,0,0)))
plot(simudata,type="l",main="White noise")

aversimudata<-filter(simudata,sides=2,rep(1,3)/3)
plot(aversimudata,type="l",main="Averaged White noise")




#differencing

T = 1024;
par(mfrow=c(2,2))
spec1 = armaspec(T, arcoeff=c(0), macoeff=c(0), vare=1);
spec2 = 2*(1-cos(2*pi*freq))


plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="White noise")
plot((1:(T/2))/T,spec2, xlab="",ylab="spectrum",type="l",main="Differenced white noise")
simudata<-arima.sim(n=T,list(order=c(0,0,0)))
plot(simudata,type="l",main="White noise")

diffsimudata<-diff(simudata)
plot(diffsimudata,type="l",main="Differenced White noise")





##########################################################################
#################################################################################
################### Filtering ###################################################
#################################################################################


#Background: Monthly values of Southern Oscillation Index (SOI), 1950-1987.
The SOI measures changes in air pressure, related to sea surface temperatures
in the central Pacific. The central Pacific Ocean warms every three to seven years due
to the El Nino effect, 
#Explain the effect of low and high pass filter
soi<-scan("soi.dat")
n<-length(soi)
 par(mfrow=c(2,2))
 plot.ts(soi) # the data
acf(soi)
 plot.ts(diff(soi),main="First Difference") # first difference
#get the local average (1/3,1/3,1/3)
aversoi<-(soi[1:(n-2)]+soi[2:(n-1)]+soi[3:n])/3
plot.ts(aversoi,main="local average")

 windows() # open new graphics device - use x11() in unix
 par(mfrow=c(3,1))
 w = seq(0,1/2, length=1000)         #--  power transfer function
 PT1  = abs(1-exp(1i*2*pi*w))^2
 plot(w, PT1, type="l",xlab="Freq",main="difference")

 PT2<- ((1+2*cos(2*pi*w))/3)^2
 plot(w, PT2, type="l",xlab="Freq",main="local average")

 PT3<- abs(1-exp(1i*8*pi*w))^2
 plot(w, PT3, type="l",xlab="Freq",main="seasonal differencing")

###################################
#effect of seasonal differencing
###################################
#generate a AR(1)_4 model with \Phi_1=0.9

T = 1024;
freq<-seq(1, T/2)/T;
par(mfrow=c(3,2))
spec1 = armaspec(T, arcoeff=c(0,0,0,0.9), macoeff=c(0), vare=1);
spec2 = spec1*abs(1-exp(1i*8*pi*freq))^2

plot((1:(T/2))/T,spec1, xlab="",ylab="spectrum",type="l",main="SAR(1), Phi_1=0.9")
plot((1:(T/2))/T,spec2, xlab="",ylab="spectrum",type="l",main="Seasonally differenced series")
simudata<-arima.sim(n=T,list(order=c(4,0,0),ar=c(0,0,0,0.9)) )



plot(simudata,type="l",main="SAR(1), Phi_1=0.9")

Sdsimudata<-diff(simudata,4)
plot(Sdsimudata,type="l",main="Seasonally differenced series")

plot(simudata[1:200],type="l",main="SAR(1), Phi_1=0.9")

plot(Sdsimudata[1:200],type="l",main="Seasonally differenced series")


#######################
#what does FFT do
#######################
x<-c(10,2,36,40,3)
n<-length(x)
dft<-fft(x)
dft

dft2<-rep(0,n)
for (j in 1:n)
{
dft2[j]<-sum(x*complex(mod=1, arg=-2*pi*(1:n)*(j-1)/n+2*pi*(j-1)/n))
}

dft2
dft-dft2
##########################
#fft(x) returns $\sum_{t=1}^{n}x_t e^{-2\pi i w_j t} e^{i2\pi w_j}$ for $j=0,1,...n-1
#for the purpose of calculating periodogram, 
$it is okay to use abs(fft(x))^2/n, where $n$ is the length of the series x
##########
###########################################################################
#Periodogram 
###############################################
 t = 1:100
 x1 = 2*cos(2*pi*t*6/100) + 3*sin(2*pi*t*6/100)
 x2 = 4*cos(2*pi*t*10/100) + 5*sin(2*pi*t*10/100)
 x3 = 6*cos(2*pi*t*40/100) + 7*sin(2*pi*t*40/100)
 x = x1 + x2 + x3
 par(mfrow=c(2,2))
 plot.ts(x1, ylim=c(-16,16), main="freq= 6/100, amp^2=13")
 plot.ts(x2, ylim=c(-16,16), main="freq= 10/100, amp^2=41")
 plot.ts(x3, ylim=c(-16,16), main="freq= 40/100, amp^2=85")
 plot.ts(x, ylim=c(-16,16), main="sum")

#look at the periodogram
#symmetric with respect to 1/2

 P = abs(fft(x))^2/100
 f = 0:99/100
 plot(f, P[1:100], type="o", xlab="frequency", ylab="periodogram")

##########
#spectral ANOVA
#####################
x = c(1,2,3,2,1)
plot(x)

abs(fft(x))^2/5
sum((x-mean(x))^2)

#a similar cos signal
y=2+cos(2*pi*1/5*(1:5-3))
plot(y,type="l",ylim=c(min(c(x,y)),max(c(x,y))))
lines(x,col=2)
points(x,pch=1)
points(y,pch=2)

abs(fft(y))^2/5
###########################



