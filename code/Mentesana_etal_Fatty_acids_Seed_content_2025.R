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
# 1) Imported the database and prepared the variable of interest
# 2) Run linear-mixed models to test for the differences in seed fatty acid 
# composition between years.

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
fa_seeds<-read.csv("Mentesana_etal_Fatty_acid_Seed_content_2025.csv", header = TRUE)
str(fa_seeds)

################################################################################
# Organizing variables
################################################################################
fa_seeds$Tree_spp <- as.factor(fa_seeds$Tree_spp)


# Transform proportions into logit
fa_seeds$SFA.prop <- (fa_seeds$Total_SFA/100)
fa_seeds$SFA.logit <- log ((fa_seeds$SFA.prop)/(1 - fa_seeds$SFA.prop))


fa_seeds$MUFA.prop <- (fa_seeds$Total_MUFA/100)
fa_seeds$MUFA.logit <- log ((fa_seeds$MUFA.prop)/(1 - fa_seeds$MUFA.prop))


fa_seeds$n6.PUFA.prop <- (fa_seeds$Total_n6_PUFA/100)
fa_seeds$n6.PUFA.logit <- log((fa_seeds$n6.PUFA.prop)/(1 - fa_seeds$n6.PUFA.prop))

fa_seeds$n3.PUFA.prop <- (fa_seeds$Total_n3_PUFA/100)
fa_seeds$n3.PUFA.logit <- log((fa_seeds$n3.PUFA.prop)/(1 - fa_seeds$n3.PUFA.prop))

################################################################################
# 2) Statistical models 
# Aim: test for year differences in seed fatty acid composition between years
################################################################################

# SFA 
sfa.year_seeds<-lm(SFA.logit ~ Tree_spp, data = fa_seeds)
summary(sfa.year_seeds)

# Checking residuals
par(mfrow = c(2, 2))
plot(sfa.year_seeds)

# I calculate posterior probabilities:
nsim <- 10000
bsim.sfa.year_seeds <- arm::sim(sfa.year_seeds, n.sim = nsim)
apply(bsim.sfa.year_seeds@coef, 2, mean) 
apply(bsim.sfa.year_seeds@coef, 2, quantile, prob = c(0.025, 0.975))


# To create the plot
newdat.sfa.year_seeds <- expand.grid(Tree_spp = levels(fa_seeds$Tree_spp))

Xmat.sfa.year_seeds <- model.matrix(~ Tree_spp, 
                              data = newdat.sfa.year_seeds)
head(Xmat.sfa.year_seeds)

fitmat.mod.sfa.year_seeds <- matrix(ncol = nsim,nrow = nrow(newdat.sfa.year_seeds))

for(i in 1:nsim) fitmat.mod.sfa.year_seeds[, i] <- 
  (1/(1 + exp(-Xmat.sfa.year_seeds%*%bsim.sfa.year_seeds@coef[i, ])))

newdat.sfa.year_seeds$lower <- apply(fitmat.mod.sfa.year_seeds, 1, 
                               quantile, prob = 0.025)
newdat.sfa.year_seeds$upper <- apply(fitmat.mod.sfa.year_seeds, 1, 
                               quantile, prob = 0.975)
newdat.sfa.year_seeds$fit <- (1/(1 + exp(-Xmat.sfa.year_seeds%*%coef(sfa.year_seeds))))


plot_sfa.year_seeds <- ggplot(newdat.sfa.year_seeds, aes(x = Tree_spp, y = fit)) + 
  geom_point(data = newdat.sfa.year_seeds, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.sfa.year_seeds, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_sfa.year_seeds_final <- plot_sfa.year_seeds +  
  geom_point(data = fa_seeds, aes(x = Tree_spp, y = SFA.prop),
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
  labs(x = "Tree_spp", y = "SFA") 



################################################################################
# MUFA

MUFA.year_seeds<-lm(MUFA.logit ~ Tree_spp, data = fa_seeds)
summary(MUFA.year_seeds)

# Checking residuals
par(mfrow = c(2, 2))
plot(MUFA.year_seeds)

# I calculate posterior probabilities:
nsim <- 10000
bsim.MUFA.year_seeds <- arm::sim(MUFA.year_seeds, n.sim = nsim)
apply(bsim.MUFA.year_seeds@coef, 2, mean) 
apply(bsim.MUFA.year_seeds@coef, 2, quantile, prob = c(0.025, 0.975))


# To create the plot
newdat.MUFA.year_seeds <- expand.grid(Tree_spp = levels(fa_seeds$Tree_spp))

Xmat.MUFA.year_seeds <- model.matrix(~ Tree_spp, 
                                    data = newdat.MUFA.year_seeds)
head(Xmat.MUFA.year_seeds)

fitmat.mod.MUFA.year_seeds <- matrix(ncol = nsim,nrow = nrow(newdat.MUFA.year_seeds))

for(i in 1:nsim) fitmat.mod.MUFA.year_seeds[, i] <- 
  (1/(1 + exp(-Xmat.MUFA.year_seeds%*%bsim.MUFA.year_seeds@coef[i, ])))

newdat.MUFA.year_seeds$lower <- apply(fitmat.mod.MUFA.year_seeds, 1, 
                                     quantile, prob = 0.025)
newdat.MUFA.year_seeds$upper <- apply(fitmat.mod.MUFA.year_seeds, 1, 
                                     quantile, prob = 0.975)
newdat.MUFA.year_seeds$fit <- (1/(1 + exp(-Xmat.MUFA.year_seeds%*%coef(MUFA.year_seeds))))


plot_MUFA.year_seeds <- ggplot(newdat.MUFA.year_seeds, aes(x = Tree_spp, y = fit)) + 
  geom_point(data = newdat.MUFA.year_seeds, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.MUFA.year_seeds, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_MUFA.year_seeds_final <- plot_MUFA.year_seeds +  
  geom_point(data = fa_seeds, aes(x = Tree_spp, y = MUFA.prop),
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
  labs(x = "Tree_spp", y = "MUFA") 


################################################################################
# n6_PUFA

n6.PUFA.year_seeds<-lm(n6.PUFA.logit ~ Tree_spp, data = fa_seeds)
summary(n6.PUFA.year_seeds)

# Checking residuals
par(mfrow = c(2, 2))
plot(n6.PUFA.year_seeds)

# I calculate posterior probabilities:
nsim <- 10000
bsim.n6.PUFA.year_seeds <- arm::sim(n6.PUFA.year_seeds, n.sim = nsim)
apply(bsim.n6.PUFA.year_seeds@coef, 2, mean) 
apply(bsim.n6.PUFA.year_seeds@coef, 2, quantile, prob = c(0.025, 0.975))


# To create the plot
newdat.n6.PUFA.year_seeds <- expand.grid(Tree_spp = levels(fa_seeds$Tree_spp))

Xmat.n6.PUFA.year_seeds <- model.matrix(~ Tree_spp, 
                                     data = newdat.n6.PUFA.year_seeds)
head(Xmat.n6.PUFA.year_seeds)

fitmat.mod.n6.PUFA.year_seeds <- matrix(ncol = nsim,nrow = nrow(newdat.n6.PUFA.year_seeds))

for(i in 1:nsim) fitmat.mod.n6.PUFA.year_seeds[, i] <- 
  (1/(1 + exp(-Xmat.n6.PUFA.year_seeds%*%bsim.n6.PUFA.year_seeds@coef[i, ])))

newdat.n6.PUFA.year_seeds$lower <- apply(fitmat.mod.n6.PUFA.year_seeds, 1, 
                                      quantile, prob = 0.025)
newdat.n6.PUFA.year_seeds$upper <- apply(fitmat.mod.n6.PUFA.year_seeds, 1, 
                                      quantile, prob = 0.975)
newdat.n6.PUFA.year_seeds$fit <- (1/(1 + exp(-Xmat.n6.PUFA.year_seeds%*%coef(n6.PUFA.year_seeds))))


plot_n6.PUFA.year_seeds <- ggplot(newdat.n6.PUFA.year_seeds, aes(x = Tree_spp, y = fit)) + 
  geom_point(data = newdat.n6.PUFA.year_seeds, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.n6.PUFA.year_seeds, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_n6.PUFA.year_seeds_final <- plot_n6.PUFA.year_seeds +  
  geom_point(data = fa_seeds, aes(x = Tree_spp, y = n6.PUFA.prop),
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
  labs(x = "Tree_spp", y = "n6.PUFA") 



################################################################################
# n3_PUFA

n3.PUFA.year_seeds<-lm(n3.PUFA.logit ~ Tree_spp, data = fa_seeds)
summary(n3.PUFA.year_seeds)

# Checking residuals
par(mfrow = c(2, 2))
plot(n3.PUFA.year_seeds)

# I calculate posterior probabilities:
nsim <- 10000
bsim.n3.PUFA.year_seeds <- arm::sim(n3.PUFA.year_seeds, n.sim = nsim)
apply(bsim.n3.PUFA.year_seeds@coef, 2, mean) 
apply(bsim.n3.PUFA.year_seeds@coef, 2, quantile, prob = c(0.025, 0.975))


# To create the plot
newdat.n3.PUFA.year_seeds <- expand.grid(Tree_spp = levels(fa_seeds$Tree_spp))

Xmat.n3.PUFA.year_seeds <- model.matrix(~ Tree_spp, 
                                        data = newdat.n3.PUFA.year_seeds)
head(Xmat.n3.PUFA.year_seeds)

fitmat.mod.n3.PUFA.year_seeds <- matrix(ncol = nsim,nrow = nrow(newdat.n3.PUFA.year_seeds))

for(i in 1:nsim) fitmat.mod.n3.PUFA.year_seeds[, i] <- 
  (1/(1 + exp(-Xmat.n3.PUFA.year_seeds%*%bsim.n3.PUFA.year_seeds@coef[i, ])))

newdat.n3.PUFA.year_seeds$lower <- apply(fitmat.mod.n3.PUFA.year_seeds, 1, 
                                         quantile, prob = 0.025)
newdat.n3.PUFA.year_seeds$upper <- apply(fitmat.mod.n3.PUFA.year_seeds, 1, 
                                         quantile, prob = 0.975)
newdat.n3.PUFA.year_seeds$fit <- (1/(1 + exp(-Xmat.n3.PUFA.year_seeds%*%coef(n3.PUFA.year_seeds))))


plot_n3.PUFA.year_seeds <- ggplot(newdat.n3.PUFA.year_seeds, aes(x = Tree_spp, y = fit)) + 
  geom_point(data = newdat.n3.PUFA.year_seeds, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.n3.PUFA.year_seeds, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_n3.PUFA.year_seeds_final <- plot_n3.PUFA.year_seeds +  
  geom_point(data = fa_seeds, aes(x = Tree_spp, y = n3.PUFA.prop),
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
  labs(x = "Tree_spp", y = "n3.PUFA") 



################################################################################
# Ratio omega 6 / omega 3
fa_seeds$Ratio_n6_n3 <- (fa_seeds$Total_n6_PUFA/fa_seeds$Total_n3_PUFA)

ratio.n6.n3.PUFA.year_seeds<-lm(Ratio_n6_n3 ~ Tree_spp, data = fa_seeds)
summary(ratio.n6.n3.PUFA.year_seeds)

# Checking residuals
par(mfrow = c(2, 2))
plot(ratio.n6.n3.PUFA.year_seeds)

# I calculate posterior probabilities:
nsim <- 10000
bsim.ratio.n6.n3.PUFA.year_seeds <- arm::sim(ratio.n6.n3.PUFA.year_seeds, n.sim = nsim)
apply(bsim.ratio.n6.n3.PUFA.year_seeds@coef, 2, mean) 
apply(bsim.ratio.n6.n3.PUFA.year_seeds@coef, 2, quantile, prob = c(0.025, 0.975))


# To create the plot
newdat.ratio.n6.n3.PUFA.year_seeds <- expand.grid(Tree_spp = levels(fa_seeds$Tree_spp))

Xmat.ratio.n6.n3.PUFA.year_seeds <- model.matrix(~ Tree_spp, 
                                        data = newdat.ratio.n6.n3.PUFA.year_seeds)
head(Xmat.ratio.n6.n3.PUFA.year_seeds)

fitmat.mod.ratio.n6.n3.PUFA.year_seeds <- matrix(ncol = nsim,nrow = nrow(newdat.ratio.n6.n3.PUFA.year_seeds))

for(i in 1:nsim) fitmat.mod.ratio.n6.n3.PUFA.year_seeds[, i] <- 
  Xmat.ratio.n6.n3.PUFA.year_seeds%*%bsim.ratio.n6.n3.PUFA.year_seeds@coef[i, ]

newdat.ratio.n6.n3.PUFA.year_seeds$lower <- apply(fitmat.mod.ratio.n6.n3.PUFA.year_seeds, 1, 
                                         quantile, prob = 0.025)
newdat.ratio.n6.n3.PUFA.year_seeds$upper <- apply(fitmat.mod.ratio.n6.n3.PUFA.year_seeds, 1, 
                                         quantile, prob = 0.975)
newdat.ratio.n6.n3.PUFA.year_seeds$fit <- Xmat.ratio.n6.n3.PUFA.year_seeds%*%coef(ratio.n6.n3.PUFA.year_seeds)


plot_ratio.n6.n3.PUFA.year_seeds <- ggplot(newdat.ratio.n6.n3.PUFA.year_seeds, aes(x = Tree_spp, y = fit)) + 
  geom_point(data = newdat.ratio.n6.n3.PUFA.year_seeds, 
             alpha = 0.90, color = "black", size = 5) +
  geom_errorbar(data = newdat.ratio.n6.n3.PUFA.year_seeds, aes(ymax = upper, ymin = lower),
                alpha = 0.90, color = "black", size = 1, width = 0.1)  

plot_ratio.n6.n3.PUFA.year_seeds_final <- plot_ratio.n6.n3.PUFA.year_seeds +  
  geom_point(data = fa_seeds, aes(x = Tree_spp, y = Ratio_n6_n3),
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
  labs(x = "Tree_spp", y = "ratio.n6.n3.PUFA") 
