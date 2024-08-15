library(ggplot2)
library(ohvbd)
library(dplyr)

## A.aegypti mass data

# scaling
# get mass from head width?
# try to predict L1 mass from its head width, 
# add another data point for ogm model
# linear regression of L2-L4, mean head width~mean wet mass

size<-read.csv("../Output/Aa_DevSize.csv")
# load data ####
# size (width + mass_mg)
size <- search_vt_smart("DatasetID", "eq", "359") %>% 
  get_vt_byid() %>%
  extract_vt_data(
    cols = c(
      "OriginalTraitDef",
      "OriginalTraitValue",
      "Interactor1Stage", 
      "Interactor1Temp"),
  )
size<-size %>% filter(OriginalTraitDef != "mean body length", 
                Interactor1Stage != "adult",
                Interactor1Stage != "pupal")
write.csv(size, "../Output/Aa_DevSize.csv")

# time_day
t_L1 <- search_vt_smart("DatasetID", "eq", "355") %>% 
  get_vt_byid() %>%
  extract_vt_data(
    cols = c(
      "OriginalTraitValue",
      "Interactor1Stage", 
      "Interactor1Temp"),
  )
t_L1<-t_L1 %>% filter(Interactor1Stage == "L1", 
                Interactor1Stage != "juvenile (not inc egg stage)")


# ln & predict ####
# seperate into width & mass
lwidth<-size %>% subset(OriginalTraitDef == "mean head width") %>%
  subset(Interactor1Stage != "L1")
lmass<-subset(size, OriginalTraitDef == "mean wet mass")

wm<-data.frame(width=lwidth$OriginalTraitValue,
               mass=lmass$OriginalTraitValue,
               ln_mass=log(lmass$OriginalTraitValue))
plot(wm$width, wm$ln_mass)
# is this a linear relationship?

# try the linear model
wm_fit<-lm(ln_mass~width, data=wm)

#plot(wm_fit)
wm_p<-ggplot(aes(width, ln_mass), data = wm)+
  geom_point()+
  geom_line(aes( y=predict(wm_fit)))
plot(wm_p)
summary(wm_fit)

# predict L1 mass
L1<-size %>% subset(OriginalTraitDef == "mean head width") %>%
  subset(Interactor1Stage == "L1")
w_L1<-data.frame(width=L1$OriginalTraitValue)
ln_mass<-predict(wm_fit, newdata=w_L1)
wm_L1<-cbind(w_L1,ln_mass)
L1$mass<-exp(wm_L1$ln_mass)

# output L1 data frame ####
# convert mass unit, mg>g
L1_pre<-data.frame(stage=L1$Interactor1Stage,
                   time=t_L1$OriginalTraitValue,
                   mass=(L1$mass)/1000,
                   Temp=t_L1$Interactor1Temp)
write.csv(L1_pre, "../Output/L1_ogm_pre.csv")
