#remove previous records
rm(list=ls())

library(xlsx) # load xlsx 
all<-read.xlsx("/Users/Michelle/Desktop/LTdata/phyto CB3.3.xlsx",1) # every samples at time series
#phyto<-read.xlsx("/Users/Michelle/Desktop/LTdata/phyto ET5.2.xlsx",1) # 12-monthly avg
head(all)

###First steps###
### autocorrelation function for phyto species ###
indkv<-is.na(all$Karlodinium.monthlyavg.)
kv<-all$Karlodinium.monthlyavg.[!indkv]
indpr<-is.na(all$Prorocentrum.monthlyavg.)
pr<-all$Prorocentrum.monthlyavg.[!indpr]
date<-all$SampleDate[!indkv]
indCsmall<-is.na(all$Crypto..10.Microns.monthlyavg.)
Csmall<-all$Crypto..10.Microns.monthlyavg.[!indCsmall] # cryptomonoas <10 Microns
indMsmall<-is.na(all$microphyto..10.Microns.monthlyavg.)
Msmall<-all$microphyto..10.Microns.monthlyavg.[!indMsmall]# micro-phytoflagellate <10 Microns

par(mfrow =c(2,2))
acf(kv, main = "Karlodinium veneficum"); acf(pr, main = "Prorocentrum minimum"); 
acf(Csmall, main = "Cryptomonas <10 Microns"); acf(Msmall, main = "Micro-phytoflagellate <10 Microns")

### acf for temp, salinity and nutrients ###
library(forecast)
temp<-na.interp(all$WTEMP)
sal<-na.interp(all$Salinity)
DIN<-na.interp(all$DIN)
NH4<-na.interp(all$NH4F)
NO3<-na.interp(all$NO3F)
NO2<-na.interp(all$NO2F)
sc<-na.interp(all$Secchi)

DIP<-na.interp(all$PO4F)
TN<-na.interp(all$TN)
TP<-na.interp(all$TP)
DON<-na.interp(all$DON)
DOP<-na.interp(all$DOP)
flow<-all$Flow.Susquehanna
par(mfrow =c(2,2))
acf(temp); acf(sal); acf(DIN); acf(DIP); 
acf(DON); acf(DOP); acf(TN); acf(TP);acf(flow)

dinp<-DIN/DIP
tnp<-TN/TP
dntp<-DIN/TP
donp<-DON/DOP

acf(dinp); acf(tnp); acf(dntp); acf(donp)

#### wavk.test ####
### small p-value means systematic change that could be a trend or a seasonal/cyclical component

library(funtimes)
kvtest<-wavk.test(kv);
prtest<-wavk.test(pr);
Csmalltest<-wavk.test(Csmall);
Msmalltest<-wavk.test(Msmall);

Temptest<-wavk.test(temp)
Saltest<-wavk.test(sal)

PO4Ftest<-wavk.test(DIP)
DINtest<-wavk.test(DIN)
DONtest<-wavk.test(DON)
DOPtest<-wavk.test(DOP)
TNtest<-wavk.test(TN)
TPtest<-wavk.test(TP)
flowtest<-wavk.test(flow)
Result<-matrix ( c(kvtest$statistic, kvtest$p.value,
                   prtest$statistic, prtest$p.value,
                   Csmalltest$statistic, Csmalltest$p.value,
                   Msmalltest$statistic, Msmalltest$p.value,
                   Temptest$statistic, Temptest$p.value, 
                   Saltest$statistic, Saltest$p.value, 
                   PO4Ftest$statistic, PO4Ftest$p.value, 
                   DINtest$statistic, DINtest$p.value, 
                   DONtest$statistic, DONtest$p.value,
                   DOPtest$statistic, DOPtest$p.value, 
                   TNtest$statistic, TNtest$p.value, 
                   TPtest$statistic, TPtest$p.value), 12,2, byrow = T, 
                 dimnames = list( c("Karlodinium","Prorocentrum","Cryptomonas","Microphytoflagellate","Temperature","salinity","PO4","DIN","DON",
                                    "DOP","TN", "TP"), c("wavk.statistic","p-value")))

Result 
wavk.test(dinp)
wavk.test(tnp)
wavk.test(dntp)
wavk.test(donp)

### parameters  at p-value < 0.05 using wavk.test needs to conduct de-seasonality:

#################(4)PO4 ######################
ydip<-as.vector(DIP)
tdip<-1:length(DIP)
rdip<-lm(ydip~tdip+tdip^2+cos(2*pi*tdip/12)+sin(2*pi*tdip/12))
wavk.test(rdip$residuals)
ddip<-DIP-predict(rdip)
wavk.test(ddip)


#################TP#####################
ytp<-as.vector(TP)
ttp<-1:length(TP)
rtp<-lm(ytp~ttp+cos(2*pi*ttp/12)+sin(2*pi*ttp/12))
wavk.test(rtp$residuals)
dtp<-TP-predict(rtp)
wavk.test(dtp)

yNH4<-as.vector(NH4)
tNH4<-1:length(NH4)
rNH4<-lm(yNH4~cos(2*pi*NH4/12)+sin(2*pi*NH4/12))
wavk.test(rNH4$residuals)
dNH4<-NH4-predict(rNH4)
wavk.test(dNH4)

##########temperature##########
yt<-as.vector(temp)

tscs<-ts(yt,frequency = 12)
na.new <- function(x) ts(na.exclude(x), frequency = 12)
c<-stl(yt, na.action = na.new, s.window = "periodic")
plot(c)
atemp<-tscs-c$time.series[,1]
wavk.test(atemp)
acf(atemp)
##############Step3: select initial variables##############
library(astsa)
lag2.plot(pr,kv,10) # dpr(t-1: r = 0.77)
lag2.plot(Msmall,kv,10) #dMsmall (t-1, r=0.14)
lag2.plot(Csmall,kv,10) #Csmall (t,r=0.18)


#Csmall and Msmall have negative correlation with KV with a negative 
#time lags, indicating 
flow<-flow[!indkv]
lag2.plot(flow,kv,10) # flow(t,r=0.45)
Temperature<-atemp[!indkv]
lag2.plot(Temperature,kv,10) # temp (t, r=0.19)
Salinity<-sal[!indkv]
lag2.plot(Salinity,kv,10) # Salinity (t-5, r= 0.29)

DissolveIN<-DIN[!indkv]
lag2.plot(DissolveIN,kv,10) #DIN(t-2, r=0.13)
ccf(DissolveIN,kv,10) 
DissolveIP<-ddip[!indkv]
lag2.plot(DissolveIP,kv,10) #DIP(t-7,r=0.28)
ccf(DissolveIP,kv,10) 

DissolveON<-DON[!indkv]
lag2.plot(DissolveON,kv,10) # DON(t-9, r=0.11)
DissolveOP<-DOP[!indkv]
lag2.plot(DissolveOP,kv,10) # DOP(t-7, r=0.1)

TotalN<-TN[!indkv]
lag2.plot(TotalN,kv,10) # TN(t-3, r=0.11)
TotalP<-dtp[!indkv]
lag2.plot(TotalP,kv,10) # TP(t-7, r=0.21)

DINP<-dinp[!indkv]
lag2.plot(DINP,kv,10) #t0 r=-0.18; t-2 r=0.28
TNP<-tnp[!indkv]
lag2.plot(TNP,kv,10) #t0 r=-0.18; t-3 r=0.24
DNTP<-dntp[!indkv]
lag2.plot(DNTP,kv,10) #t0 r=-0.2; t-2 r=0.15
DONP<-donp[!indkv]
lag2.plot(DONP,kv,10) # t-7 r= -0.15


#####

D = cbind(kv,prlag1 = lag(pr,-1), sallag2 = lag(Salinity,-2), DINlag7=lag(DissolveIN,-7), DissolveIN, 
          DIPlag7 = lag(DissolveIP, -7), TNlag7= lag(TotalN,-7), TPlag7 = lag(TotalP,-7))
cor(D)
data = na.omit(D)
dat = data.frame (data)

Dr = cbind (kv, prlag1 = lag(pr,-1), sallag2 = lag(Salinity,-2), 
            DINPlag3 = lag(DINP,-3), DINPlag2 = lag(DINP,-2), TNPlag7 = lag(TNP,-7), TNPlag3 = lag(TNP,-3))

cor(Dr)
datadr = na.omit(Dr)
datdr = data.frame(datadr)

library(GGally)
library(Hmisc)
ggpairs(D,axisLabels = "internal", lab.cex = 0.8)
rcorr(as.matrix(D),type="pearson")

model= lm (kv~ prlag1+sallag2 + DissolveIN+ DIPlag7 + TNlag7 + TPlag7, data=dat)
model1 = lm (kv~ prlag1+ sallag2+DINPlag2 + TNPlag7  + TNPlag3, data=datdr)

library(MASS)
step<-stepAIC(model)
step<-stepAIC(model1)
step$anova

######step7 diagnostic checks######
model.lm = lm(kv~sallag2 + DissolveIN , data = dat)
modelr.lm = lm(kv~DINPlag2 ,data = datdr)
summary(model.lm)
summary(modelr.lm)

model.res = resid(model.lm, data=dat)

par(mfrow =c(2,2))
plot(density(model.res)) 
qqnorm(resid(model.lm))
qqline(resid(model.lm))
acf(model.res); pacf(model.res)
# a portmanteau test returns large p-value suggest the residuals are white noise
Box.test(model.res,lag = 24, fitdf = 4, type = "Ljung")
shapiro.test(model.lm$residuals)

wavk.test(resid(model1.lm))

wavk.test(resid(model1a.lm))

#### incorporate ARMA structure 
x<-predict(model1.lm)
fit<-arima(x,order = c(3,0,2))
summary(fit)
Box.test(resid(fit),lag = 24, fitdf = 5 , type = "Ljung")
qqnorm(resid(fit))
qqline(resid(fit))
acf(resid(fit));pacf(resid(fit))
wavk.test(resid(fit))

