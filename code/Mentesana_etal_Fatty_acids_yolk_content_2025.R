################################################################################
# Author: 
# Lucia Mentesana (lucia.mentesana@fcien.edu.uy) wrote the code.

# Script created in May 2025

################################################################################
# Description of script and Instructions
################################################################################

# Manuscript titled "Variation in seed abundance predicts yolk fatty acid 
# composition in a wild population of birds"


# In this script, I:
# 1) Imported the database and prepared all variables
# 2) Prepared the final data base for analysis (a - I removed from the general
# data base 11 eggs that were analyzed separately; b - prepared covariates and
# fixed factors for statistical models).
# 2) Run linear-mixed models to test for the differences in fatty acid composition
# between years.

################################################################################
# Packages needed
################################################################################

# install.packages('stringr')
# install.packages('ggplot2')
# install.packages('data.table')
# install.packages('lme4')
# install.packages('arm')
# install.packages('blmeco')
# install.packages('reshape2')
# install.packages('BBmisc')
# install.packages('plyr')
# install.packages('rptR')
# install.packages('Rtools')
# install.packages("DHARMa")


require(knitr)
require(tidyr)
require(dplyr)
require(ggplot2)
require(broom)
require(gridExtra)
require(MCMCglmm)
require(arm)
require(lme4)
require(reshape2)
require(DHARMa)

################################################################################
# 1) Importing the data 
################################################################################
all.eggs<-read.csv("Mentesana_etal_Fatty_acid_yolk_content_2025.csv", header = TRUE)
str(all.eggs)

################################################################################
# Organizing variables
################################################################################
all.eggs$Nest <- as.factor(all.eggs$Nest)
all.eggs$Year <- as.factor(all.eggs$Year)
all.eggs$Clutch.number <- as.factor(all.eggs$Clutch.number)
all.eggs$Prop.Total.SFA <- as.numeric(all.eggs$Prop.Total.SFA)
all.eggs$Prop.Total.MUFA <- as.numeric(all.eggs$Prop.Total.MUFA)
all.eggs$Prop.Total.n.3.PUFA <- as.numeric(all.eggs$Prop.Total.n.3.PUFA)
all.eggs$Prop.Total.n.6.PUFA <- as.numeric(all.eggs$Prop.Total.n.6.PUFA)
all.eggs$Female <- as.factor(all.eggs$Female)

################################################################################
# 2a) Preparing final data 
################################################################################

# In 2017, we analyzed a subset of eggs collected in 2015 (N = 11 entire clutches) 
# because we wanted to study variation within clutch in yolk components 
# (i.e., hormones, # antioxidants and fatty acids). The results from these 
# analyses are already published: Mentesana et al. (2019): "Female variation in 
# allocation of steroid hormones, antioxidants and fatty acids: a multilevel 
# analysis in a wild passerine bird" (https://doi.org/10.1111/jav.01859)

# To address the question of the current manuscript, in 2019 we analyzed the 
# remaining eggs collected in 2015 and all eggs collected in 2016. 
# Because eggs were analyzed over different years, it is 
# very likley that differences between batches arise not because the egg content
# was different but because they were analyzed by different people. Therefore, 
# wetook a conservative approach and for the current analysis we only worked with
# eggs analyzed in 2019 (N = 112 eggs).

# To subset the database, I first need to create a unique variable that consideres
# Nest ID, year and clutch number. This is because the same nest box could have
# been used for more than one clutch within the same year, and between years.
# By creating this unique variable, I make sure that the final subset of the 
# data excludes eggs analyzed in 2017.

all.eggs$nest_year <- paste(all.eggs$Nest, all.eggs$Year, sep = "_")
all.eggs$nest_year_clutch <- paste(all.eggs$nest_year, 
                                   all.eggs$Clutch.number, sep = "_")

all.eggs$nest_year_clutch <- as.factor(all.eggs$nest_year_clutch)

all.eggs.2015_2016 <- droplevels(subset(all.eggs, all.eggs$nest_year_clutch!="95_2015_1" &
                                all.eggs$nest_year_clutch!="97_2015_1"&
                                all.eggs$nest_year_clutch!="100_2015_1"&
                                all.eggs$nest_year_clutch!="129_2015_1"&
                                all.eggs$nest_year_clutch!="136_2015_1"&
                                all.eggs$nest_year_clutch!="138_2015_1"&
                                all.eggs$nest_year_clutch!="163_2015_1"&
                                all.eggs$nest_year_clutch!="170_2015_1"&
                                all.eggs$nest_year_clutch!="177_2015_1"&
                                all.eggs$nest_year_clutch!="189_2015_1"&
                                all.eggs$nest_year_clutch!="198_2015_1"))

nrow(all.eggs.2015_2016)
# 112 

all.eggs.2015 <- droplevels(subset(all.eggs.2015_2016, all.eggs.2015_2016$Year == "2015"))
table(all.eggs.2015$nest_year)
table(all.eggs.2015$Clutch.number)

all.eggs.2016 <- droplevels(subset(all.eggs.2015_2016, all.eggs.2015_2016$Year == "2016"))
table(all.eggs.2016$Clutch.number)
################################################################################
# 2b) Organizing variables: covariates and fixed factors for statistical models
################################################################################

# Temperature before capture
str(all.eggs.2015_2016$Temp.3.days.before)
# continuous variable
all.eggs.2015_2016$Temp.3.days.before_scale <- scale(all.eggs.2015_2016$Temp.3.days.before)

# Date of collection
str(all.eggs.2015_2016$Date.collection.num)
# continuous variable
all.eggs.2015_2016$Date.collection.num_scale <- scale(all.eggs.2015_2016$Date.collection.num)

# Transform proportions into logit
all.eggs.2015_2016$SFA.prop <- (all.eggs.2015_2016$Prop.Total.SFA/100)
all.eggs.2015_2016$SFA.logit <- log ((all.eggs.2015_2016$SFA.prop)/(1 - all.eggs.2015_2016$SFA.prop))


all.eggs.2015_2016$MUFA.prop <- (all.eggs.2015_2016$Prop.Total.MUFA/100)
all.eggs.2015_2016$MUFA.logit <- log ((all.eggs.2015_2016$MUFA.prop)/(1 - all.eggs.2015_2016$MUFA.prop))


all.eggs.2015_2016$n6.PUFA.prop <- (all.eggs.2015_2016$Prop.Total.n.6.PUFA/100)
all.eggs.2015_2016$n6.PUFA.logit <- log((all.eggs.2015_2016$n6.PUFA.prop)/(1 - all.eggs.2015_2016$n6.PUFA.prop))

all.eggs.2015_2016$n3.PUFA.prop <- (all.eggs.2015_2016$Prop.Total.n.3.PUFA/100)
all.eggs.2015_2016$n3.PUFA.logit <- log((all.eggs.2015_2016$n3.PUFA.prop)/(1 - all.eggs.2015_2016$n3.PUFA.prop))


################################################################################
# 3) Statistical models 
# Aim: test for year differences in fatty acid composition between years
################################################################################

# SFA 
sfa.year<-lmer(SFA.logit ~ Year + Date.collection.num_scale + Temp.3.days.before_scale +
                          (1|Female) + (1|nest_year), data = all.eggs.2015_2016)
summary(sfa.year)


# 1 - Tukey-Anscombe Plot
par(mfrow = c(2, 2))
# fixed effects
scatter.smooth(fitted(sfa.year), resid(sfa.year)); 
abline(h = 0, lty = 2)

# 2 - QQ of residuals
qqnorm(resid(sfa.year))
qqline(resid(sfa.year))

# 3 - res.vari vs fitted
scatter.smooth(fitted(sfa.year), sqrt(abs(resid(sfa.year))))

# 4 random effects
qqnorm(ranef(sfa.year)$Female[, 1])
qqline(ranef(sfa.year)$Female[, 1])

qqnorm(ranef(sfa.year)$nest_year[, 1])
qqline(ranef(sfa.year)$nest_year[, 1])

# I calculate posterior probabilities:
nsim <- 10000
bsim.sfa.year <- arm::sim(sfa.year, n.sim = nsim)
apply(bsim.sfa.year@fixef, 2, mean) 
apply(bsim.sfa.year@fixef, 2, quantile, prob = c(0.025, 0.975))

# I want to see the size of the effects for each fixed factor or covariate:

mean(bsim.sfa.year@fixef[, 2] < 0)
# 99%

mean(bsim.sfa.year@fixef[, 3] < 0)
# 99%

mean(bsim.sfa.year@fixef[, 4] < 0)
# 33%


# Random factors:
# Median and 95% CrI
round(quantile (apply(bsim.sfa.year@ranef$Female[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

round(quantile (apply(bsim.sfa.year@ranef$nest_year[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

# Residual variance
# Median and 95% CrI
round(as.vector(quantile(bsim.sfa.year@sigma ^ 2,
                         c(0.025, 0.5, 0.975))), 3) 


# To create the plot
newdat.sfa.year <- expand.grid(Year = levels(all.eggs.2015_2016$Year),
                               Date.collection.num_scale = mean(all.eggs.2015_2016$Date.collection.num_scale),
                               Temp.3.days.before_scale = mean(all.eggs.2015_2016$Temp.3.days.before_scale))

Xmat.sfa.year <- model.matrix(~ Year + Date.collection.num_scale + 
                                Temp.3.days.before_scale, 
                                data = newdat.sfa.year)
head(Xmat.sfa.year)

fitmat.mod.sfa.year <- matrix(ncol = nsim,nrow = nrow(newdat.sfa.year))

for(i in 1:nsim) fitmat.mod.sfa.year[, i] <- 
  (1/(1 + exp(-Xmat.sfa.year%*%bsim.sfa.year@fixef[i, ])))

newdat.sfa.year$lower <- apply(fitmat.mod.sfa.year, 1, 
                               quantile, prob = 0.025)
newdat.sfa.year$upper <- apply(fitmat.mod.sfa.year, 1, 
                               quantile, prob = 0.975)
newdat.sfa.year$fit <- (1/(1 + exp(-Xmat.sfa.year%*%fixef(sfa.year))))


plot_sfa.year <- ggplot(newdat.sfa.year, aes(x = Year, y = fit)) + 
  geom_point(data = newdat.sfa.year, 
              alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.sfa.year, aes(ymax = upper, ymin = lower),
             alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_sfa.year_final <- plot_sfa.year +  
  geom_point(data = all.eggs.2015_2016, aes(x = Year, y = SFA.prop),
             color = "grey69", size = 3.5, alpha = 0.3, position = position_nudge(x = 0.1)) +  
  theme(panel.background = element_blank(),legend.key=element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        legend.title = element_text(colour = "black",size = 15),
        strip.text.y = element_text(colour = "black",size = 5, face = "bold",angle = 65),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.position = "right", 
        legend.direction = "vertical") +   
  labs(x = "Year", y = "SFA") 



################################################################################
# MUFA 
MUFA.year<-lmer(MUFA.logit ~ Year + Date.collection.num_scale + Temp.3.days.before_scale +
                 (1|nest_year), data = all.eggs.2015_2016)

# Female explains almost no variance and can therefore be excluded from final 
# model. 

summary(MUFA.year)

# to check differences between years:
ggplot(all.eggs.2015_2016, aes(x = Date.collection.num, y = MUFA.logit)) +
geom_point() +
facet_wrap(~Year, ncol = 1, scales = "fixed") +  # ensures identical axes
labs(x = "Date", y = "MUFA", title = "Comparison by Year") +
theme_minimal() 

# Checking residuals: 

# 1 - Tukey-Anscombe Plot
par(mfrow = c(2, 2))
# fixed effects
scatter.smooth(fitted(MUFA.year), resid(MUFA.year)); 
abline(h = 0, lty = 2)

# 2 - QQ of residuals
qqnorm(resid(MUFA.year))
qqline(resid(MUFA.year))

# 3 - res.vari vs fitted
scatter.smooth(fitted(MUFA.year), sqrt(abs(resid(MUFA.year))))

# 4 random effects
qqnorm(ranef(MUFA.year)$nest_year[, 1])
qqline(ranef(MUFA.year)$nest_year[, 1])


# I calculate posterior probabilities:
nsim <- 10000
bsim.MUFA.year <- arm::sim(MUFA.year, n.sim = nsim)
apply(bsim.MUFA.year@fixef, 2, mean) 
apply(bsim.MUFA.year@fixef, 2, quantile, prob = c(0.025, 0.975))

# I want to see the size of the effects for each fixed factor or covariate:

mean(bsim.MUFA.year@fixef[, 2] < 0)
# 1%

mean(bsim.MUFA.year@fixef[, 3] < 0)
# 0%

mean(bsim.MUFA.year@fixef[, 4] < 0)
# 82%


# Random factors:
# Median and 95% CrI
round(quantile (apply(bsim.MUFA.year@ranef$nest_year[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

# Residual variance
# Median and 95% CrI
round(as.vector(quantile(bsim.MUFA.year@sigma ^ 2,
                         c(0.025, 0.5, 0.975))), 3) 



# To create the plot
newdat.mufa.year <- expand.grid(Year = levels(all.eggs.2015_2016$Year),
                               Date.collection.num_scale = mean(all.eggs.2015_2016$Date.collection.num_scale),
                               Temp.3.days.before_scale = mean(all.eggs.2015_2016$Temp.3.days.before_scale))


Xmat.mufa.year <- model.matrix(~ Year + Date.collection.num_scale + 
                                Temp.3.days.before_scale, 
                              data = newdat.mufa.year)
head(Xmat.mufa.year)

fitmat.mod.mufa.year <- matrix(ncol = nsim,nrow = nrow(newdat.mufa.year))

for(i in 1:nsim) fitmat.mod.mufa.year[, i] <- 
  (1/(1 + exp(-Xmat.mufa.year%*%bsim.MUFA.year@fixef[i, ])))

newdat.mufa.year$lower <- apply(fitmat.mod.mufa.year, 1, 
                               quantile, prob = 0.025)
newdat.mufa.year$upper <- apply(fitmat.mod.mufa.year, 1, 
                               quantile, prob = 0.975)
newdat.mufa.year$fit <- (1/(1 + exp(-Xmat.mufa.year%*%fixef(MUFA.year))))


plot_mufa.year <- ggplot(newdat.mufa.year, aes(x = Year, y = fit)) + 
  geom_point(data = newdat.mufa.year, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.mufa.year, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_mufa.year_final <- plot_mufa.year +  
  geom_point(data = all.eggs.2015_2016, aes(x = Year, y = MUFA.prop),
             color = "grey69", size = 3.5, alpha = 0.3, position = position_nudge(x = 0.1)) +  
  theme(panel.background = element_blank(),legend.key=element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        legend.title = element_text(colour = "black",size = 15),
        strip.text.y = element_text(colour = "black",size = 5, face = "bold",angle = 65),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.position = "right", 
        legend.direction = "vertical") +   
  labs(x = "Year", y = "MUFA") 


##################################################################################
# n6 PUFA
n6.PUFA.year<-lmer(n6.PUFA.logit ~ Year + Date.collection.num_scale + Temp.3.days.before_scale +
                 (1|nest_year), data = all.eggs.2015_2016)

summary(n6.PUFA.year)

# 1 - Tukey-Anscombe Plot
par(mfrow = c(2, 2))
# fixed effects
scatter.smooth(fitted(n6.PUFA.year), resid(n6.PUFA.year)); 
abline(h = 0, lty = 2)

# 2 - QQ of residuals
qqnorm(resid(n6.PUFA.year))
qqline(resid(n6.PUFA.year))

# 3 - res.vari vs fitted
scatter.smooth(fitted(n6.PUFA.year), sqrt(abs(resid(n6.PUFA.year))))

# 4 random effects
qqnorm(ranef(n6.PUFA.year)$nest_year[, 1])
qqline(ranef(n6.PUFA.year)$nest_year[, 1])


# I calculate posterior probabilities:
nsim <- 10000
bsim.n6.PUFA.year <- arm::sim(n6.PUFA.year, n.sim = nsim)
apply(bsim.n6.PUFA.year@fixef, 2, mean) 
apply(bsim.n6.PUFA.year@fixef, 2, quantile, prob = c(0.025, 0.975))

# I want to see the size of the effects for each fixed factor or covariate:

mean(bsim.n6.PUFA.year@fixef[, 2] < 0)
# 1%

mean(bsim.n6.PUFA.year@fixef[, 3] < 0)
# 0%

mean(bsim.n6.PUFA.year@fixef[, 4] < 0)
# 97%


# Random factors:
# Median and 95% CrI
round(quantile (apply(bsim.n6.PUFA.year@ranef$nest_year[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

# Residual variance
# Median and 95% CrI
round(as.vector(quantile(bsim.n6.PUFA.year@sigma ^ 2,
                         c(0.025, 0.5, 0.975))), 3) 



# To create the plot
newdat.n6.PUFA.year <- expand.grid(Year = levels(all.eggs.2015_2016$Year),
                                Date.collection.num_scale = mean(all.eggs.2015_2016$Date.collection.num_scale),
                                Temp.3.days.before_scale = mean(all.eggs.2015_2016$Temp.3.days.before_scale))


Xmat.n6.PUFA.year <- model.matrix(~ Year + Date.collection.num_scale + 
                                 Temp.3.days.before_scale, 
                               data = newdat.n6.PUFA.year)
head(Xmat.n6.PUFA.year)

fitmat.mod.n6.PUFA.year <- matrix(ncol = nsim,nrow = nrow(newdat.n6.PUFA.year))

for(i in 1:nsim) fitmat.mod.n6.PUFA.year[, i] <- 
  (1/(1 + exp(-Xmat.n6.PUFA.year%*%bsim.n6.PUFA.year@fixef[i, ])))

newdat.n6.PUFA.year$lower <- apply(fitmat.mod.n6.PUFA.year, 1, 
                                quantile, prob = 0.025)
newdat.n6.PUFA.year$upper <- apply(fitmat.mod.n6.PUFA.year, 1, 
                                quantile, prob = 0.975)
newdat.n6.PUFA.year$fit <- (1/(1 + exp(-Xmat.n6.PUFA.year%*%fixef(n6.PUFA.year))))


plot_n6.PUFA.year <- ggplot(newdat.n6.PUFA.year, aes(x = Year, y = fit)) + 
  geom_point(data = newdat.n6.PUFA.year, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.n6.PUFA.year, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_n6.PUFA.year_final <- plot_n6.PUFA.year +  
  geom_point(data = all.eggs.2015_2016, aes(x = Year, y = n6.PUFA.prop),
             color = "grey69", size = 3.5, alpha = 0.3, position = position_nudge(x = 0.1)) +  
  theme(panel.background = element_blank(),legend.key=element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        legend.title = element_text(colour = "black",size = 15),
        strip.text.y = element_text(colour = "black",size = 5, face = "bold",angle = 65),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.position = "right", 
        legend.direction = "vertical") +   
  labs(x = "Year", y = "n6.PUFA") 


##################################################################################
# n3 PUFA

n3.PUFA.year<-lmer(n3.PUFA.logit ~ Year + Date.collection.num_scale + 
                     Temp.3.days.before_scale +
                     (1|nest_year), data = all.eggs.2015_2016)

# Female explains no variance. I removed it. 
summary(n3.PUFA.year)


# 1 - Tukey-Anscombe Plot
par(mfrow = c(2, 2))
# fixed effects
scatter.smooth(fitted(n3.PUFA.year), resid(n3.PUFA.year)); 
abline(h = 0, lty = 2)

# 2 - QQ of residuals
qqnorm(resid(n3.PUFA.year))
qqline(resid(n3.PUFA.year))

# 3 - res.vari vs fitted
scatter.smooth(fitted(n3.PUFA.year), sqrt(abs(resid(n3.PUFA.year))))

# 4 random effects
qqnorm(ranef(n3.PUFA.year)$nest_year[, 1])
qqline(ranef(n3.PUFA.year)$nest_year[, 1])


# I calculate posterior probabilities:
nsim <- 10000
bsim.n3.PUFA.year <- arm::sim(n3.PUFA.year, n.sim = nsim)
apply(bsim.n3.PUFA.year@fixef, 2, mean) 
apply(bsim.n3.PUFA.year@fixef, 2, quantile, prob = c(0.025, 0.975))

# I want to see the size of the effects for each fixed factor or covariate:

mean(bsim.n3.PUFA.year@fixef[, 2] < 0)
# 1%

mean(bsim.n3.PUFA.year@fixef[, 3] < 0)
# 0%

mean(bsim.n3.PUFA.year@fixef[, 4] < 0)
# 3%


# Random factors:
# Median and 95% CrI
round(quantile (apply(bsim.n3.PUFA.year@ranef$nest_year[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

# Residual variance
# Median and 95% CrI
round(as.vector(quantile(bsim.n3.PUFA.year@sigma ^ 2,
                         c(0.025, 0.5, 0.975))), 3) 



# To create the plot
newdat.n3.PUFA.year <- expand.grid(Year = levels(all.eggs.2015_2016$Year),
                                   Date.collection.num_scale = mean(all.eggs.2015_2016$Date.collection.num_scale),
                                   Temp.3.days.before_scale = mean(all.eggs.2015_2016$Temp.3.days.before_scale))


Xmat.n3.PUFA.year <- model.matrix(~ Year + Date.collection.num_scale + 
                                    Temp.3.days.before_scale, 
                                  data = newdat.n3.PUFA.year)
head(Xmat.n3.PUFA.year)

fitmat.mod.n3.PUFA.year <- matrix(ncol = nsim,nrow = nrow(newdat.n3.PUFA.year))

for(i in 1:nsim) fitmat.mod.n3.PUFA.year[, i] <- 
  (1/(1 + exp(-Xmat.n3.PUFA.year%*%bsim.n3.PUFA.year@fixef[i, ])))

newdat.n3.PUFA.year$lower <- apply(fitmat.mod.n3.PUFA.year, 1, 
                                   quantile, prob = 0.025)
newdat.n3.PUFA.year$upper <- apply(fitmat.mod.n3.PUFA.year, 1, 
                                   quantile, prob = 0.975)
newdat.n3.PUFA.year$fit <- (1/(1 + exp(-Xmat.n3.PUFA.year%*%fixef(n3.PUFA.year))))


plot_n3.PUFA.year <- ggplot(newdat.n3.PUFA.year, aes(x = Year, y = fit)) + 
  geom_point(data = newdat.n3.PUFA.year, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.n3.PUFA.year, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_n3.PUFA.year_final <- plot_n3.PUFA.year +  
  geom_point(data = all.eggs.2015_2016, aes(x = Year, y = n3.PUFA.prop),
             color = "grey69", size = 3.5, alpha = 0.3, position = position_nudge(x = 0.1)) +  
  theme(panel.background = element_blank(),legend.key=element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        legend.title = element_text(colour = "black",size = 15),
        strip.text.y = element_text(colour = "black",size = 5, face = "bold",angle = 65),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.position = "right", 
        legend.direction = "vertical") +   
  labs(x = "Year", y = "n3.PUFA") 



##################################################################################
# Ratio omega 6/omega 3

ratio.n6.n3.year<-lmer(Prop.Ratio.n6.n3 ~ Year + Date.collection.num_scale + 
                     Temp.3.days.before_scale + 
                     (1|nest_year), data = all.eggs.2015_2016)

# Female explains no variance. I removed it. 
summary(ratio.n6.n3.year)


# 1 - Tukey-Anscombe Plot
par(mfrow = c(2, 2))
# fixed effects
scatter.smooth(fitted(ratio.n6.n3.year), resid(ratio.n6.n3.year)); 
abline(h = 0, lty = 2)

# 2 - QQ of residuals
qqnorm(resid(ratio.n6.n3.year))
qqline(resid(ratio.n6.n3.year))

# 3 - res.vari vs fitted
scatter.smooth(fitted(ratio.n6.n3.year), sqrt(abs(resid(ratio.n6.n3.year))))

# 4 random effects
qqnorm(ranef(ratio.n6.n3.year)$nest_year[, 1])
qqline(ranef(ratio.n6.n3.year)$nest_year[, 1])


# I calculate posterior probabilities:
nsim <- 10000
bsim.ratio.n6.n3.year <- arm::sim(ratio.n6.n3.year, n.sim = nsim)
apply(bsim.ratio.n6.n3.year@fixef, 2, mean) 
apply(bsim.ratio.n6.n3.year@fixef, 2, quantile, prob = c(0.025, 0.975))

# I want to see the size of the effects for each fixed factor or covariate:

mean(bsim.ratio.n6.n3.year@fixef[, 2] < 0)
# 1%

mean(bsim.ratio.n6.n3.year@fixef[, 3] < 0)
# 0%

mean(bsim.ratio.n6.n3.year@fixef[, 4] < 0)
# 84%


# Random factors:
# Median and 95% CrI
round(quantile (apply(bsim.ratio.n6.n3.year@ranef$nest_year[ , , 1], 1, var),
                prob = c(0.025, 0.5, 0.975)), 3) 

# Residual variance
# Median and 95% CrI
round(as.vector(quantile(bsim.ratio.n6.n3.year@sigma ^ 2,
                         c(0.025, 0.5, 0.975))), 3) 



# To create the plot
newdat.ratio.n6.n3.year <- expand.grid(Year = levels(all.eggs.2015_2016$Year),
                                   Date.collection.num_scale = mean(all.eggs.2015_2016$Date.collection.num_scale),
                                   Temp.3.days.before_scale = mean(all.eggs.2015_2016$Temp.3.days.before_scale))


Xmat.ratio.n6.n3.year <- model.matrix(~ Year + Date.collection.num_scale + 
                                    Temp.3.days.before_scale, 
                                  data = newdat.ratio.n6.n3.year)
head(Xmat.ratio.n6.n3.year)

fitmat.mod.ratio.n6.n3.year <- matrix(ncol = nsim,nrow = nrow(newdat.ratio.n6.n3.year))

for(i in 1:nsim) fitmat.mod.ratio.n6.n3.year[, i] <- 
  Xmat.ratio.n6.n3.year%*%bsim.ratio.n6.n3.year@fixef[i, ]

newdat.ratio.n6.n3.year$lower <- apply(fitmat.mod.ratio.n6.n3.year, 1, 
                                   quantile, prob = 0.025)
newdat.ratio.n6.n3.year$upper <- apply(fitmat.mod.ratio.n6.n3.year, 1, 
                                   quantile, prob = 0.975)
newdat.ratio.n6.n3.year$fit <- Xmat.ratio.n6.n3.year%*%fixef(ratio.n6.n3.year)


plot_ratio.n6.n3.year <- ggplot(newdat.ratio.n6.n3.year, aes(x = Year, y = fit)) + 
  geom_point(data = newdat.ratio.n6.n3.year, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.ratio.n6.n3.year, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_ratio.n6.n3.year_final <- plot_ratio.n6.n3.year +  
  geom_point(data = all.eggs.2015_2016, aes(x = Year, y = Prop.Ratio.n6.n3),
             color = "grey69", size = 3.5, alpha = 0.3, position = position_nudge(x = 0.1)) +  
  theme(panel.background = element_blank(),legend.key=element_blank(),
        axis.line = element_line(colour = "black",size = 1),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(colour = "black", size = 20),
        axis.text.y = element_text(colour = "black",size = 10),
        axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black",size = 15),
        legend.title = element_text(colour = "black",size = 15),
        strip.text.y = element_text(colour = "black",size = 5, face = "bold",angle = 65),
        legend.text = element_text(colour = "black", size = 10, face = "bold"),
        legend.position = "right", 
        legend.direction = "vertical") +   
  labs(x = "Year", y = "Ratio n6/n3") 


t.test(all.eggs.2015$Fledgling.number, all.eggs.2016$Fledgling.number)
t.test(all.eggs.2015$Clutch.size, all.eggs.2016$Clutch.size)
t.test(all.eggs.2015$Hatching.number, all.eggs.2016$Hatching.number)
