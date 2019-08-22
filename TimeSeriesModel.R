#remove previous records
rm(list=ls())

library(xlsx) # load xlsx library
all<-read.xlsx("/Users/Michelle/Desktop/LTdata/phyto CB3.3.xlsx",1) # time series data
#phyto<-read.xlsx("/Users/Michelle/Desktop/LTdata/phyto CB3.3.xlsx",1) # 12-monthly average
head(all)

###Step 1: apply autocorrelation function (ACF) for phytoplankton species ###

# Data cleaning for ACF analysis 
indkv<-is.na(all$Karlodinium.monthlyavg.) # index for nan data of phytoplankton predator 
kv<-all$Karlodinium.monthlyavg.[!indkv] # get the number for non-nan data 
indpr<-is.na(all$Prorocentrum.monthlyavg.)
pr<-all$Prorocentrum.monthlyavg.[!indpr]
date<-all$SampleDate[!indkv]
indCsmall<-is.na(all$Crypto..10.Microns.monthlyavg.)
Csmall<-all$Crypto..10.Microns.monthlyavg.[!indCsmall] # cryptomonoas <10 Microns
indMsmall<-is.na(all$microphyto..10.Microns.monthlyavg.)
Msmall<-all$microphyto..10.Microns.monthlyavg.[!indMsmall]# micro-phytoflagellate <10 Microns

# Displying the ACF result for phytoplankton species, including predator and prey
par(mfrow =c(2,2))
acf(kv, main = "Karlodinium veneficum"); acf(pr, main = "Prorocentrum minimum"); 
acf(Csmall, main = "Cryptomonas <10 Microns"); acf(Msmall, main = "Micro-phytoflagellate <10 Microns")

# For independent varialbles: temp, salinity and nutrients, use na.interp to interpolate missing values in time series (sinces there has more frequent data points and very few missing data). Those data has strong seasonality, so a robust STL decomposition is first compurted. Then, a linear interpolation is applied to the seasonally adjusted data, and the seasonal component is added back.


library(forecast)
# independent variables: physical conditions of water environment   
temp<-na.interp(all$WTEMP); sal<-na.interp(all$Salinity); sc<-na.interp(all$Secchi); 
flow<-all$Flow.Susquehanna

# independet varilables: different nutrients forms 
DIN<-na.interp(all$DIN) # dissolved inorganic nitroent
NH4<-na.interp(all$NH4F) # ammonium; NO3<-na.interp(all$NO3F) # nitrate; NO2<-na.interp(all$NO2F) # nitrite
DIP<-na.interp(all$PO4F) # dissolved inorganic phosphorus
TN<-na.interp(all$TN) # total nitrogen; TP<-na.interp(all$TP) # total phosphate
DON<-na.interp(all$DON) # Dissolved organic nitrogen; DOP<-na.interp(all$DOP) # Dissolved organic phosphate

# Displying the ACF result for independent variables 
par(mfrow =c(2,2))
acf(temp); acf(sal); acf(DIN); acf(DIP); 
acf(DON); acf(DOP); acf(TN); acf(TP);acf(flow)

dinp<-DIN/DIP # inorganic N/P ratio; tnp<-TN/TP # total N/P ratio; 
dntp<-DIN/TP # inorganic N/ total P ratio; donp<-DON/DOP # dissolved organic N/P ratio

acf(dinp); acf(tnp); acf(dntp); acf(donp)

#### Step 2: Conducting wavk.test. if you have small p-value, which means systematic change that could be a trend or a seasonal/cyclical component

library(funtimes)
kvtest<-wavk.test(kv); prtest<-wavk.test(pr); Csmalltest<-wavk.test(Csmall); Msmalltest<-wavk.test(Msmall);

Temptest<-wavk.test(temp); Saltest<-wavk.test(sal);flowtest<-wavk.test(flow)

PO4Ftest<-wavk.test(DIP); DINtest<-wavk.test(DIN); DONtest<-wavk.test(DON); DOPtest<-wavk.test(DOP)
TNtest<-wavk.test(TN); TPtest<-wavk.test(TP)

# Summarize the Result
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
                   TPtest$statistic, TPtest$p.value), 12,2, byrow = T, dimnames = list( c("Karlodinium","Prorocentrum","Cryptomonas","Microphytoflagellate","Temperature","salinity","PO4","DIN","DON","DOP","TN", "TP"), c("wavk.statistic","p-value")))

# Result 
# wavk.test(dinp); wavk.test(tnp); wavk.test(dntp); wavk.test(donp)

# Notes: parameters at p-value < 0.05 using wavk.test needs to conduct de-seasonality, which you have conducted it individually as below sections

# De-seasonality Sections using lm function and re-do the wavk.test 
### PO4 
ydip<-as.vector(DIP); tdip<-1:length(DIP); rdip<-lm(ydip~tdip+tdip^2+cos(2*pi*tdip/12)+sin(2*pi*tdip/12))
wavk.test(rdip$residuals); ddip<-DIP-predict(rdip); wavk.test(ddip)
### TP
ytp<-as.vector(TP); ttp<-1:length(TP); rtp<-lm(ytp~ttp+cos(2*pi*ttp/12)+sin(2*pi*ttp/12))
wavk.test(rtp$residuals); dtp<-TP-predict(rtp); wavk.test(dtp)
### NH4
yNH4<-as.vector(NH4); tNH4<-1:length(NH4); rNH4<-lm(yNH4~cos(2*pi*NH4/12)+sin(2*pi*NH4/12))
wavk.test(rNH4$residuals); dNH4<-NH4-predict(rNH4); wavk.test(dNH4)

### temperature
yt<-as.vector(temp); tscs<-ts(yt,frequency = 12)
na.new <- function(x) ts(na.exclude(x), frequency = 12)
c<-stl(yt, na.action = na.new, s.window = "periodic")
plot(c)
atemp<-tscs-c$time.series[,1]
wavk.test(atemp)
acf(atemp)



### Step 3: Select initial variables, use astsa to gain lag-lag relationships among phytoplatnkon predator (dependent varilable: kv) and independent variables ###
library(astsa)
lag2.plot(pr,kv,10) # dpr(t-1: r = 0.77)
lag2.plot(Msmall,kv,10) #dMsmall (t-1, r=0.14); lag2.plot(Csmall,kv,10) #Csmall (t,r=0.18)
# Prey Csmall and Msmall have a zero or 1 day time-lag with predator karlodinium

flow<-flow[!indkv]; lag2.plot(flow,kv,10) # flow(t,r=0.45)
Temperature<-atemp[!indkv]; lag2.plot(Temperature,kv,10) # temp (t, r=0.19)
Salinity<-sal[!indkv]; lag2.plot(Salinity,kv,10) # Salinity (t-5, r= 0.29)

DissolveIN<-DIN[!indkv]; lag2.plot(DissolveIN,kv,10) #DIN(t-2, r=0.13) # ccf(DissolveIN,kv,10) 
DissolveIP<-ddip[!indkv]; lag2.plot(DissolveIP,kv,10) #DIP(t-7,r=0.28) # ccf(DissolveIP,kv,10) 

DissolveON<-DON[!indkv]; lag2.plot(DissolveON,kv,10) # DON(t-9, r=0.11)
DissolveOP<-DOP[!indkv]; lag2.plot(DissolveOP,kv,10) # DOP(t-7, r=0.1)

TotalN<-TN[!indkv]; lag2.plot(TotalN,kv,10) # TN(t-3, r=0.11)
TotalP<-dtp[!indkv]; lag2.plot(TotalP,kv,10) # TP(t-7, r=0.21)

DINP<-dinp[!indkv]; lag2.plot(DINP,kv,10) #t0 r=-0.18; t-2 r=0.28
TNP<-tnp[!indkv]; lag2.plot(TNP,kv,10) #t0 r=-0.18; t-3 r=0.24
DNTP<-dntp[!indkv]; lag2.plot(DNTP,kv,10) #t0 r=-0.2; t-2 r=0.15
DONP<-donp[!indkv]; lag2.plot(DONP,kv,10) # t-7 r= -0.15


### Step 4: Combine the variables and check correlations again

D = cbind(kv,prlag1 = lag(pr,-1), sallag2 = lag(Salinity,-2), DINlag7=lag(DissolveIN,-7), DissolveIN, 
          DIPlag7 = lag(DissolveIP, -7), TNlag7= lag(TotalN,-7), TPlag7 = lag(TotalP,-7))
cor(D)
data = na.omit(D); dat = data.frame (data)

Dr = cbind (kv, prlag1 = lag(pr,-1), sallag2 = lag(Salinity,-2), 
            DINPlag3 = lag(DINP,-3), DINPlag2 = lag(DINP,-2), TNPlag7 = lag(TNP,-7), TNPlag3 = lag(TNP,-3))

cor(Dr)
datadr = na.omit(Dr); datdr = data.frame(datadr)

# check if there any multicollinearity, with |r|>0.8
library(GGally); library(Hmisc)
ggpairs(D,axisLabels = "internal", lab.cex = 0.8) 
rcorr(as.matrix(D),type="pearson")

### Step 5: Estimate a full model 

model= lm (kv~ prlag1+sallag2 + DissolveIN+ DIPlag7 + TNlag7 + TPlag7, data=dat)
model1 = lm (kv~ prlag1+ sallag2+DINPlag2 + TNPlag7  + TNPlag3, data=datdr)

### Step 6: Shrink model by selecing variables based on AIC

library(MASS)
step<-stepAIC(model)
step<-stepAIC(model1)
step$anova

### Step 7 : diagnostic checks
model.lm = lm(kv~sallag2 + DissolveIN , data = dat)
modelr.lm = lm(kv~DINPlag2 ,data = datdr)
summary(model.lm); summary(modelr.lm)

model.res = resid(model.lm, data=dat)

par(mfrow =c(2,2))
plot(density(model.res)); qqnorm(resid(model.lm)); qqline(resid(model.lm))
acf(model.res); pacf(model.res)

# a portmanteau test returns large p-value suggest the residuals are white noise
Box.test(model.res,lag = 24, fitdf = 4, type = "Ljung")
shapiro.test(model.lm$residuals)

wavk.test(resid(model1.lm)); wavk.test(resid(model1a.lm))

### Optinal: incorporate ARMA structure 
x<-predict(model1.lm)
fit<-arima(x,order = c(3,0,2))
summary(fit)
Box.test(resid(fit),lag = 24, fitdf = 5 , type = "Ljung")
qqnorm(resid(fit))
qqline(resid(fit))
acf(resid(fit));pacf(resid(fit))
wavk.test(resid(fit))

