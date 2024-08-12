# setting
getwd()
setwd("D:/Biological_Sciences_IC/R/Code")
# packages needed
library(ggplot2)
library(ohvbd)
library(dplyr)
library(minpack.lm)
library(cowplot)

# Applying to an insect population (Aedes aegypti)

## Growth

# if step 1 was ran before, read in data
df<-read.csv("../Output/Aa_at0_ogm.csv")
### 1. get data ##########
# pulled from VT dataset (https://www.vectorbyte.org/)
# reference doi: 10.1093/jmedent/27.5.892
# Mass(mg)-time(day), under a range of temperatures

# mass data
m <- search_vt_smart("DatasetID", "eq", "359") %>% 
  get_vt_byid() %>%
  extract_vt_data(
    cols = c(
      "OriginalTraitValue",
      "OriginalTraitUnit",
      "Interactor1Stage", 
      "Interactor1Temp"),
  )
# subset for only mass data
# exclude adult because their mass decreases
m<-m %>% filter(OriginalTraitUnit == "mg", Interactor1Stage != "adult")

# time data
t <- search_vt_smart("DatasetID", "eq", "355") %>% 
  get_vt_byid() %>%
  extract_vt_data(
    cols = c(
      "OriginalTraitValue",
      "Interactor1Stage", 
      "Interactor1Temp"),
  )
t<-t %>% filter(Interactor1Stage != "L1", 
                Interactor1Stage != "juvenile (not inc egg stage)")

# new dataframe
time<-t$OriginalTraitValue
# convert mass unit, mg>g
mass<-(m$OriginalTraitValue)/1000
Temp<-m$Interactor1Temp
stage<-t$Interactor1Stage
df<-data.frame(stage, time, mass, Temp)
# save the data
write.csv(df, "../Output/Aa_at0_ogm.csv")


### 2. calculate a(T) #########
# only 4 data points for each temperature, a bit insufficient?

  # subset data into different temperatures ####
t15<-subset(df, Temp==15)
t20<-subset(df, Temp==20)
t25<-subset(df, Temp==25)
t27<-subset(df, Temp==27)
t30<-subset(df, Temp==30)
t34<-subset(df, Temp==34)
  # get cumulative time columns ####
t15$Time<-cumsum(t15$time)
t20$Time<-cumsum(t20$time)
t25$Time<-cumsum(t25$time)
t27$Time<-cumsum(t27$time)
t30$Time<-cumsum(t30$time)
t34$Time<-cumsum(t34$time)
  # plot, check that model (growth curve) suits the shape ####
par(mfrow = c(2, 3))
plot(t15$Time, t15$mass)
plot(t20$Time, t20$mass)
plot(t25$Time, t25$mass) # flatting
plot(t27$Time, t27$mass) # flatting
plot(t30$Time, t30$mass)
plot(t34$Time, t34$mass)
dev.off()

  # model: growth equation (10) from Olivia (https://docs.google.com/document/d/1woAKI_FgjuBd3x44OBA4plzuH2DkRFsojzo4GsUKyo8/edit) ####
#Ontogentic Growth Model. 
ogm <- function(t, a, t0, m0 = 0.0004, mA = 0.0035){ #t: time in days from above function, Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT_0^2
  k <- 8.62*(10^(-5))
  m <- mA * ( 1- (1-((m0/mA)^0.25)) * (exp(-a*(t+t0)/(4*(mA^0.25))) )) ^4
  return(m) #mass reached after amount of time t
}
  # fitting ####
ogm_nls_15<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t15, start=list(a = 0.01, t0 = 0) )
ogm_nls_20<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t20, start=list(a = 0.01, t0 = 0) )
ogm_nls_25<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t25, start=list(a = 0.01, t0 = 0) )
ogm_nls_27<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t27, start=list(a = 0.01, t0 = 0) )
ogm_nls_30<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t30, start=list(a = 0.01, t0 = 0) )
ogm_nls_34<- nls(mass ~ ogm(t = Time, a, t0),
                 data=t34, start=list(a = 0.01, t0 = 0) )

  # cehck fitting with plot ####
p15<-ggplot(aes(Time, mass), data = t15)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_15)))
plot(p15)

p20<-ggplot(aes(Time, mass), data = t20)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_20)))
p25<-ggplot(aes(Time, mass), data = t25)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_25)))
p27<-ggplot(aes(Time, mass), data = t27)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_27)))
p30<-ggplot(aes(Time, mass), data = t30)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_30)))
p34<-ggplot(aes(Time, mass), data = t34)+
  geom_point()+
  geom_line(aes( y=predict(ogm_nls_34)))

plot_grid(p15,p20,p25,p27,p30,p34)

  # get a at each temperature ####
#Boltzmann's constant
k<- 8.62*(10^(-5)) 
#parameter a from fits to data
a<-c(coef(ogm_nls_15)[1], coef(ogm_nls_20)[1], coef(ogm_nls_25)[1], coef(ogm_nls_27)[1], coef(ogm_nls_30)[1], coef(ogm_nls_34)[1])
#Time elapsed
t0<-c(coef(ogm_nls_15)[2], coef(ogm_nls_20)[2], coef(ogm_nls_25)[2], coef(ogm_nls_27)[2], coef(ogm_nls_30)[2], coef(ogm_nls_34)[2])
#Temperature
temps<-unique(df$Temp)
ats<-as.data.frame(cbind(a, t0, temps))
# check by plot a~temp (should like TPC curve, not quite?)
plot(ats$temps, ats$a)


### 3. get at0 ##########
#Inverse temperature
ats$inversetemp<- (1/( k * (ats$temps +273))) - (1/( k * (273)))
# linear regratiom
lm(log(ats$a)~ats$inversetemp)
aiT_fit<-lm(log(ats$a)~ats$inversetemp)
# visualize
plot(ats$inversetemp, log(ats$a))
abline(lm(log(ats$a)~ats$inversetemp), col="blue")
# get at0
ait<-exp(coef(aiT_fit)[1])
ait

### 4. growth prediction #######
M_Simulator<-function()

#######