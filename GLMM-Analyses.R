# R version 4.2.3 (2023-03-15) -- "Shortstop Beagle"

# Required packages
library(tidyr) #replace NAs with zeros
library(dplyr) #group_by
library(ggplot2)
library(lme4) #glmm
library(lmtest) #lrtest #not used currently; better to use anova() for nested models; though results are identical in this case
library(glmmTMB) #just in case there's non-convergence w/ lme4
library(MuMIn) #pseudo r-squared
library(merTools) #predictInterval
library(boot) #inverse logit backtransform
library(emmeans) #tukey's test
library(effects)

################################################################*
# 1. JC seed experiment####
# Research questions:
# Is seed survival density or distance dependent?
# Does human impact affect seed survival?
# Data: seed fate in dispersal manipulation experiment. Observed after two months.

ebony1 <- read.csv(file.path("data", "ExtFile_Seed_fate_vs_density&distance","ExtFile_Seed_fate_vs_density&distance_Data.csv"),   header = T, na.strings=c("NA", "NULL", "", "."))
names(ebony1)<-tolower(names(ebony1))

# organize data
ebony1$locality<-as.factor(ebony1$locality)
ebony1$group<-as.factor(ebony1$group)
ebony1$density<-as.factor(ebony1$density)
ebony1$remained<- ifelse(ebony1$removed == 1, 0, 1) # the opposite of removed is remained 
# distance as numeric b/c there's three distances
# in this analysis #1, models with spatialgroup/group as random effect did not converge.  So continuing without it. 

# global model (w/ interaction)  
global<-glmer(remained~density+distance*locality+(1|group), family=binomial, data=ebony1)
summary(global) #everything w/ good support except density

# global minus density ->BEST FIT
lm1<-glmer(remained~distance*locality+(1|group), family=binomial, data=ebony1)
summary(lm1) #everything well supported
anova(global, lm1) #models are NOT different ...so density doesn't control any variation 

# extracting fixed effect parameter estimates from best-fit model for table
summary_table <- summary(lm1)
fixed_effects_table1 <- summary_table$coefficients

# extracting LRT results
lrt<-anova(global, lm1)
lrt.g1<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# above w/o locality or interaction 
lm2<-glmer(remained~distance+(1|group), family=binomial, data=ebony1)
summary(lm2) #everything well supported

anova(lm1, lm2)
AIC(lm1, lm2)

# extracting LRT results
lrt<-anova(lm1,lm2)
lrt.12<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# additive model
lm3<-glmer(remained~distance+locality+(1|group), family=binomial, data=ebony1)
summary(lm3) #

anova(lm3, lm1) # models are different - interaction improves things considerably vs. additive
AIC(lm3, lm1) 
lrt<-anova(lm1,lm3)
lrt.13<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# null model for comparison
null<-glmer(remained~(1|group), family=binomial, data=ebony1)
anova(lm1, null) #yes, best-fit is much better than the null
lrt<-anova(lm1,null)
lrt.1n<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# cbinding all the LRT results & AICs and exporting for supp table
selection.tab<-data.frame(rbind(lrt.g1, lrt.12, lrt.13, lrt.1n))
names(selection.tab)<-c("x1","Chisq", "x2", "df", "x3","P")
selection.tab<-subset(selection.tab, select = c("Chisq", "df", "P"))
write.csv(selection.tab, file.path("outputs", "selection.tab1.csv"))
aic<-AIC(global, lm1, lm2, lm3, null)
write.csv(aic, file.path("outputs", "aic1.csv"))

# checking univariate locality model
lm.locality<-glmer(remained~locality+(1|group), family=binomial, data=ebony1)
summary(lm.locality) # marginal additive effect of locality

# checking an additive model w/ distance as categorical
ebony1$distance.cat<-as.factor(ebony1$distance)
lm4<-glmer(remained~distance.cat+locality+(1|group), family=binomial, data=ebony1)
summary(lm4) 

anova(lm4, lm1) #fits much nicer as a categorical!  We might want to explore this?
AIC(lm4, lm1) #no, the categorical variable doesn't add anything vs. the best-fit interactive model


# predicted values

# finding a non-influential random effect (group) to base predictions on
ranef(lm1)
pred.mx <- expand.grid(locality = c("Bouamir","Kompia"),
                       distance = c(2, 3, 5, 7, 10, 15, 20, 30, 40, 50),  #2, 20, and 50 are measured values
                       group = "VD0324")# least influential random effect that won't skew predictions

preds <- predictInterval(lm1, newdata = pred.mx, level=0.95, n.sims = 999, which="fixed",include.resid.var=FALSE) #excludes random effect variance
preds<-cbind(pred.mx, preds)

# back transform, because  family=binomial specifies the logit transformation
preds$upr<-inv.logit(preds$upr)
preds$lwr<-inv.logit(preds$lwr)
preds$fit<-inv.logit(preds$fit)

write.csv(preds, file.path("outputs", "preds1.csv"))

# double-checking the param estimates ***NOT CURRENTLY USED
effect.factor1 <-emmeans::emmeans(model5,tukey~landscape,type="response")
effect.factor1
eff.factor1 <- data.frame(effects::effect("distance:locality", lm1.int,
                                          xlevels=list(distance=c(2,20,50))))

# Pseudo R-squared on best-fit 
r.squaredGLMM(lm1)

# on univariate models
r.squaredGLMM(lm2) 
r.squaredGLMM(lm.locality) 

################################################################*
# 2. JC seedling experiment####
# Research questions:
# Is seedling survival density or distance dependent?
# Does human impact affect seedling survival?
# Data: seedling fate in dispersal manipulation experiment. Observed after two years.

ebony2 <- read.csv(file.path("data", "JC_w_spatialgroup","ExtFile_Seedling_fate_vs_density&distance_Data.csv"),   header = T, na.strings=c("NA", "NULL", "", ".")) 
names(ebony2)<-tolower(names(ebony2))

# organize data
ebony2$locality<-as.factor(ebony2$locality)
ebony2$group<-as.factor(ebony2$group)
ebony2$spatialgroup<-as.factor(ebony2$spatialgroup)
ebony2$density<-as.factor(ebony2$density)
# distance as numeric b/c there's three distances

# global w/ interaction
global<-glmer(survived~density+distance*locality+(1|spatialgroup/group), family=binomial, data=ebony2)
summary(global) # pretty good support except for the interaction

# global model minus iteraction *BEST FIT based on marginal significance
global.add<-glmer(survived~density+distance+locality+(1|spatialgroup/group), family=binomial, data=ebony2)
summary(global.add) # everything very significant except for density (which is marginal)

anova(global.add, global)   # P = 0.848 
AIC(global.add, global)  # Interactive model is worse - no evidence for interaction

# extracting fixed effect parameter estimates from best-fit model for table
summary_table <- summary(global.add)
fixed_effects_table2 <- summary_table$coefficients


lrt<-anova(global, global.add)
lrt.g.ga<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# removing density 
lm1<-glmer(survived~distance+locality+(1|spatialgroup/group), family=binomial, data=ebony2)
summary(lm1) # everything with strong support


anova(global.add, lm1) 
AIC(global.add, lm1) # no difference. No evidence for density
lrt<-anova( global.add, lm1)
lrt.ga.1<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# ditching locality
lm2<-glmer(survived~distance+(1|spatialgroup/group), family=binomial, data=ebony2)
summary(lm2) 

anova(lm2, lm1)  
AIC(lm2, lm1)

lrt<-anova( lm1,lm2)
lrt.12<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# compare lm1 to null
null<-glmer(survived~1+(1|spatialgroup/group), family=binomial, data=ebony2)
anova(null, lm1)   #P = <.0001
AIC(null, lm1)  # lm1 is better than NULL

lrt<-anova( lm1,null)
lrt.1n<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# compare global.add to null
anova(null, global.add)   #P = <.0001
AIC(null, lm1)  # lm1 is much better than NULL

lrt<-anova( lm1,null)
lrt.1n<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

# cbinding all the LRT results & AICs and exporting for supp table 
selection.tab2<-data.frame(rbind(lrt.g.ga, lrt.ga.1, lrt.12, lrt.1n))
names(selection.tab2)<-c("x1","Chisq", "x2", "df", "x3","P")
selection.tab2<-subset(selection.tab2, select = c("Chisq", "df", "P"))
write.csv(selection.tab2, file.path("outputs", "selection.tab2.csv"))
aic<-AIC(global, global.add, lm1, lm2, null)
write.csv(aic, file.path("outputs", "aic2.csv"))

lm4<-glmer(survived~locality+(1|spatialgroup/group), family=binomial, data=ebony2)
summary(lm4)
anova(lm4, null)

# Pseudo R-squared 
r.squaredGLMM(lm1)

# Pseudo R-squared 
r.squaredGLMM(global.add)

anova(lm2, global.add)

# predicted values lm1
ranef(lm1) #viewing random effect (RE) params
pred.mx <- expand.grid(locality = c("Bouamir","Kompia"),
                       distance = c(2, 3, 5, 7, 10, 15, 20, 30, 40, 50),
                       group = "VD0402",
                       spatialgroup="SpatialGroup_seeds_07")# least influential REs that won't skew predictions

preds <- predictInterval(lm1, newdata = pred.mx, level=0.95, n.sims = 999, which="fixed",include.resid.var=FALSE) #excludes random effect variance
preds<-cbind(pred.mx, preds)

# back transform, because  family=binomial specifies the logit transformation
preds$upr<-inv.logit(preds$upr)
preds$lwr<-inv.logit(preds$lwr)
preds$fit<-inv.logit(preds$fit)
write.csv(preds, file.path("outputs", "preds2.csv"))

# predicted values, global.add
ranef(global.add) #viewing random effect (RE) params
pred.mx <- expand.grid(locality = c("Bouamir","Kompia"),
                       distance = c(2, 3, 5, 7, 10, 15, 20, 30, 40, 50),
                       density = c("2","25"),
                       group = "VD0402",
                       spatialgroup="SpatialGroup_seeds_07") # least influential REs that won't skew predictions

preds <- predictInterval(global.add, newdata = pred.mx, level=0.95, n.sims = 999, which="fixed",include.resid.var=FALSE) #excludes random effect variance
preds<-cbind(pred.mx, preds)

# back transform, because  family=binomial specifies the logit transformation
preds$upr<-inv.logit(preds$upr)
preds$lwr<-inv.logit(preds$lwr)
preds$fit<-inv.logit(preds$fit)
write.csv(preds, file.path("outputs", "preds2_globaladd.csv"))

################################################################*
# 3. Dung experiment (in cage)####
# Research questions:
# Does Environment (fruit, ground, dung, dung+gut passage) affect seed germination rate?
# Data: seeds in cage
# Model: germinated~treatment+(1|spatialgroup/group)

ebony3 <- read.csv(file.path("data", "ExtFile_Seed_fate_vs_dispersal","ExtFile_Seed_fate_vs_dispersal_Data.csv"),   header = T, na.strings=c("NA", "NULL", "", ".")) #calls file from subfolder(s) within WD

# organize data
names(ebony3)<-tolower(names(ebony3))
ebony3$year<-as.factor(ebony3$year)
ebony3$treatment<-as.factor(ebony3$treatment)
ebony3$group<-as.factor(ebony3$group) #tree
ebony3$spatialgroup<-as.factor(ebony3$spatialgroup) #watershed

# subsetting
ebony4<-subset(ebony3, cage==0) #for analysis 4 (out of cage only)
ebony3<-subset(ebony3, cage==1) #for this analysis (in cage only)


# global model ->BEST FIT MODEL
global<-glmer(germinated~treatment+(1|spatialgroup/group), family=binomial, data=ebony3)
summary(global) #

null<-glmer(germinated~1+(1|spatialgroup/group), family=binomial, data=ebony3)

# extracting fixed effect parameter estimates from best-fit model for table
summary_table <- summary(global)
fixed_effects_table3 <- summary_table$coefficients

anova(null, global)   #P = <.0001
AIC(null, global)  #global is better than NULL

lrt<-anova( null, global)
lrt.ng<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)
write.csv(lrt.ng, file.path("outputs", "selection.tab3.csv"))
aic3<-AIC(global, null)
write.csv(aic3, file.path("outputs", "aic3.csv"))

# Pseudo R-squared on best-fit
r.squaredGLMM(global)

# predicted values
ranef(global) # viewing random effect (RE) params
pred.mx <- expand.grid(treatment = c("Dung","Dung+Elephant", "Ground", "Fruit"),
                       group = "Group_55_2022", 
                       spatialgroup="SpatialGroup_05_2022" )# least influential REs that won't skew predictions

preds <- predictInterval(global, newdata = pred.mx, level=0.95, n.sims = 999, which="fixed",include.resid.var=FALSE) # excludes random effect variance
preds<-cbind(pred.mx, preds)

# back transform, because  family=binomial specifies the logit transformation
preds$upr<-inv.logit(preds$upr)
preds$lwr<-inv.logit(preds$lwr)
preds$fit<-inv.logit(preds$fit)
write.csv(preds, file.path("outputs", "preds3.csv"))

# Pairwise treatment group comparisons
effect.factor1 <-emmeans::emmeans(global,tukey~treatment,type="response")
effect.factor1
write.csv(effect.factor1$contrasts, file.path("outputs", "Tukeys3.csv"))

################################################################*
# 4. Dung experiment (out of cage)####
# Research questions:
# Environment (fruit, ground, dung) affect seed predation rate?
# Data: seeds out of cage (two years of data)
# Model: remained~treatment+year+(1|spatialgroup/group)

# the opposite of removed is remained 
ebony4$remained<- ifelse(ebony4$removed == 1, 0, 1)

# global model 
global<-glmer(remained~treatment+year+(1|spatialgroup/group), family=binomial, data=ebony4)
summary(global) #

# exclude year -> BEST FIT MODEL
lm1<-glmer(remained~treatment+(1|spatialgroup/group), family=binomial, data=ebony4)
summary(lm1) 

# extracting fixed effect parameter estimates from best-fit model for table
summary_table <- summary(lm1)
fixed_effects_table4 <- summary_table$coefficients

# binding all four of the fixed effect tables together and exporting
# Inserting the experiment name and number into the table
e1<-c("J-C seed experiment", "", "", "")
e2<-c("J-C seedling experiment", "", "", "")
e3<-c("Dung experiment (in cage)", "", "", "")
e4<-c("Dung experiment (out of cage)", "", "", "")
comment<-c("All best fit-models above include nested random effects of group within watershed, except the J-C seed experiment, which includes only the random effect of group.", "", "", "")
fixed.tab<-data.frame(rbind(e1, fixed_effects_table1, e2, fixed_effects_table2, e3, fixed_effects_table3,e4, fixed_effects_table4, comment))
write.csv(fixed.tab, file.path("outputs", "fixed.table.csv"))

anova(lm1, global)   # no effect of year
lrt<-anova(global, lm1)
lrt.g1<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)
AIC(lm1, global)

lrt.g1<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

null<-glmer(remained~1+(1|spatialgroup/group), family=binomial, data=ebony4)
summary(null) #

anova(null, lm1)   # P = <.0001
AIC(null, lm1)  # global is better than NULL
lrt<-anova(lm1,null)
lrt.1n<-c(lrt$Chisq, lrt$Df, lrt$`Pr(>Chisq)`)

selection.tab4<-data.frame(rbind(lrt.g1, lrt.1n))
names(selection.tab4)<-c("x1","Chisq", "x2", "df", "x3","P")
selection.tab4<-subset(selection.tab4, select = c("Chisq", "df", "P"))
write.csv(selection.tab4, file.path("outputs", "selection.tab4.csv"))

aic4<-AIC(global, lm1, null)
write.csv(aic4, file.path("outputs", "aic4.csv"))

# Pseudo R-squared on best-fit
r.squaredGLMM(lm1)

# predicted values
ranef(lm1) # viewing random effect (RE) params
pred.mx <- expand.grid(treatment = c("Dung", "Ground", "Fruit"),
                       group = "Group_07_2021", 
                       spatialgroup="SpatialGroup_01_2021" )# least influential REs that won't skew predictions

preds <- predictInterval(lm1, newdata = pred.mx, level=0.95, n.sims = 999, which="fixed",include.resid.var=FALSE) # excludes random effect variance
preds<-cbind(pred.mx, preds)

# back transform, because  family=binomial specifies the logit transformation
preds$upr<-inv.logit(preds$upr)
preds$lwr<-inv.logit(preds$lwr)
preds$fit<-inv.logit(preds$fit)

write.csv(preds, file.path("outputs", "preds4.csv"))

# Pairwise treatment group comparisons
effect.factor1 <-emmeans::emmeans(lm1,tukey~treatment,type="response")
effect.factor1
# also produces estimates for each factor
write.csv(effect.factor1$contrasts, file.path("outputs", "Tukeys4.csv"))