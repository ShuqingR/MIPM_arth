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
df_ogm<-read.csv("../Output/Aa_at0_ogm.csv")
### 1. get data ##########
# pulled from VT dataset (https://www.vectorbyte.org/)
# reference doi: 10.1093/jmedent/27.5.892
  # L1 mass predicyion from head width (../Code/1_) ####

  # Mass(mg)-time(day), under a range of temperatures ####

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

# read in L1 prediction data
L1_pre<-read.csv("../Output/L1_ogm_pre.csv")
# combine
df_ogm<-rbind(L1_pre, df)

# save the data
write.csv(df_ogm, "../Output/Aa_at0_ogm.csv")



### 2. calculate a(T) #########
# only 4 data points for each temperature, a bit insufficient?
# so prediction of L1 mass was made using linear regression of head_width~mass

  # set L1 time to 0
df_ogm$time[1:6]<-0
  # subset data into different temperatures ####
t15<-subset(df_ogm, Temp==15)
t20<-subset(df_ogm, Temp==20)
t25<-subset(df_ogm, Temp==25)
t27<-subset(df_ogm, Temp==27)
t30<-subset(df_ogm, Temp==30)
t34<-subset(df_ogm, Temp==34)
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
ogm <- function(t, a, t0, m0 = 0.0002, mA){ #t: time in days from above function, Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT_0^2
  k <- 8.62*(10^(-5))
  m <- mA * ( 1- (1-((m0/mA)^0.25)) * (exp(-a*(t+t0)/(4*(mA^0.25))) )) ^4
  return(m) #mass reached after amount of time t
}
  # fitting ####
# with estimated mA from plot
ogm_nls_15<- nls(mass ~ ogm(t = Time, a, t0, mA=0.0093),
                 data=t15, start=list(a = 0.01, t0 = 0) )
ogm_nls_20<- nls(mass ~ ogm(t = Time, a, t0, mA=0.0075),
                 data=t20, start=list(a = 0.01, t0 = 0) )
ogm_nls_25<- nls(mass ~ ogm(t = Time, a, t0, mA=0.005),
                 data=t25, start=list(a = 0.01, t0 = 0) )
ogm_nls_27<- nls(mass ~ ogm(t = Time, a, t0, mA=0.0045),
                 data=t27, start=list(a = 0.01, t0 = 0) )
ogm_nls_30<- nls(mass ~ ogm(t = Time, a, t0, mA=0.0043),
                 data=t30, start=list(a = 0.01, t0 = 0) )
ogm_nls_34<- nls(mass ~ ogm(t = Time, a, t0, mA=0.0041),
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

  # more check (growth curve prediction) ####
# model plotting was zig-zag, try with more Time points
test<-data.frame(Time=0:80)
test15<-predict(ogm_nls_15, newdata=test)
test20<-predict(ogm_nls_20, newdata=test)
test25<-predict(ogm_nls_25, newdata=test)
test27<-predict(ogm_nls_27, newdata=test)
test30<-predict(ogm_nls_30, newdata=test)
test34<-predict(ogm_nls_34, newdata=test)

# plot mass against Time
plot(test$Time,test15, type = "l", col="blue", lwd=2,
     xlab="Time (day)", ylab="mass (g)")
lines(test$Time,test20, col="cyan3", lwd=2)
lines(test$Time,test25, col="green3", lwd=2)
lines(test$Time,test27, col="pink2", lwd=2)
lines(test$Time,test30, col="orange2", lwd=2)
lines(test$Time,test34, col="red3", lwd=2)
leg<-legend("bottomright",
            legend=c("15","20","25","27","30","34"),
            col=c("blue","cyan3","green3","pink2","orange2","red3"),
            pch=15,
            ncol=2)
dev.off()

# plot log(mass) against Time
plot(test$Time[1:40],log(test15[1:40]), type = "l", col="blue", lwd=2,
     xlab="Time (day)", ylab="ln(mass)")
lines(test$Time[1:40],log(test20[1:40]), col="cyan3", lwd=2)
lines(test$Time[1:40],log(test25[1:40]), col="green3", lwd=2)
lines(test$Time[1:40],log(test27[1:40]), col="pink2", lwd=2)
lines(test$Time[1:40],log(test30[1:40]), col="orange2", lwd=2)
lines(test$Time[1:40],log(test34[1:40]), col="red3", lwd=2)
leg<-legend("bottomright",
            legend=c("15","20","25","27","30","34"),
            col=c("blue","cyan3","green3","pink2","orange2","red3"),
            pch=15,
            ncol=2)
dev.off()

  # get a at each temperature ####
#Boltzmann's constant
k<- 8.62*(10^(-5)) 
#parameter a from fits to data
a<-c(coef(ogm_nls_15)[1], coef(ogm_nls_20)[1], coef(ogm_nls_25)[1], coef(ogm_nls_27)[1], coef(ogm_nls_30)[1], coef(ogm_nls_34)[1])
#Time elapsed
t0<-c(coef(ogm_nls_15)[2], coef(ogm_nls_20)[2], coef(ogm_nls_25)[2], coef(ogm_nls_27)[2], coef(ogm_nls_30)[2], coef(ogm_nls_34)[2])
#Temperature
temps<-unique(df_ogm$Temp)
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
at0<-exp(coef(aiT_fit)[1])
at0


### 4. growth curve simulation ####

# constant parameter
k = 8.62*(10^(-5))
#a. Development time from mass NEW
development_time <- function(m, Tc, at0, E=0.45, m0 = 0.0002, M = 0.0093){ #m: mass(g), Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT20
  alpha = E/(k*(273^2))
  t <-  log( (1-( (m0/M)^(0.25) ))/ ( 1- ((m/M)^(0.25)) ) ) * (4 * (M^0.25) / ( at0 * exp(alpha* ( Tc/(1+(Tc/273))) )) ) 
  return(t)   #time in days to reach mass m
}

#b. Mass reached over developmental time NEW
development_mass <- function(t, Tc, at0, E = 0.45, m0 = 0.0002, M = 0.0093){ #t: time in days from above function, Tc: temp degrees celsius, t0: normalisation constant, a: slope, alpha=-E/kT_0^2
  k <- 8.62*(10^(-5))
  alpha = E/(k*(273^2)) 
  a_temp<- at0 * exp(alpha* (Tc/(1+(Tc/273))) )
  m <- M * ( 1- (1-((m0/M)^0.25)) * (exp(-a_temp*t/(4*(M^0.25))) )) ^4
  return(m) #mass reached after amount of time t
}

# simulatior function -> t-m, fixed temp.
G_sim<-function(Tc, tEnd = 80, mStart = 0.00022, at0 = 0.01317569){
  # empty lists storing simulation values
  timeVec<-c(0)
  massVec<-c(0.0002)
  # initial time
  t<-0
  
  # set simulation loop limit
  while(t<tEnd){
    # calculate new values
    dt<-development_time(mStart, Tc, at0)
    t<-dt+t
    m<-development_mass(t, Tc, at0)
    # replace old variable
    mStart<-m
    # record new values
    timeVec<-c(timeVec, t)
    massVec<-c(massVec, m)
    }
  result<-data.frame(timeVec, massVec)
  return(result)
  }

# simulate under a range of temp, plot
s10<-G_sim(10)
s15<-G_sim(15)
s20<-G_sim(20)
s25<-G_sim(25)
s30<-G_sim(30)
s35<-G_sim(35)
# check
plot(s10$timeVec, log(s10$massVec))
# plot
plot(s35$timeVec,log(s35$massVec), col="red3", lwd=2, type = "l",
     xlab="Time (day)", ylab="ln_mass (g)")
lines(s15$timeVec,log(s15$massVec), col="cyan3", lwd=2)
lines(s20$timeVec,log(s20$massVec), col="green3", lwd=2)
lines(s25$timeVec,log(s25$massVec), col="pink2", lwd=2)
lines(s30$timeVec,log(s30$massVec), col="orange2", lwd=2)
lines(s10$timeVec,log(s10$massVec), col="blue", lwd=2)
leg<-legend("bottomright",
            legend=c("10","15","20","25","30","35"),
            col=c("blue","cyan3","green3","pink2","orange2","red3"),
            pch=15,
            ncol=2)
dev.off()


#######