###############################
###############################
## R code to accompany:      ##
## Schoemann (2015)          ##
## SEM with lavaan           ##
## Presented at SPSP, 2015   ##
###############################
###############################

library(lavaan)
library(psych)
library(semTools)
library(semPlot)

## CFA model
mod <- '
energetic =~ active + vigorous + wakeful + lively 
negative =~  jittery + nervous + scared
'
#Default is a marker variable method of scale setting
fit <- cfa(mod, data = msq)
summary(fit)
#Get path diagram
semPaths(fit, "est", nCharNodes = 0)

#sd.lv changes to a fixex factor method of scale setting
fit1 <- cfa(mod, data = msq, std.lv = TRUE)
summary(fit1)

semPaths(fit1, "est", nCharNodes = 0)

## Multiple group CFA

fitg <- cfa(mod, data = msq, group = "scale")
summary(fitg)

#Set loadings equal across groups
fitgW <- cfa(mod, data = msq, group = "scale",
            group.equal = "loadings")
#Compare models with anova
anova(fitg, fitgW)

#Test for invariance across groups
measurementInvariance(mod, data = msq, group = "scale")

## Cateogrical indicators
ability <- data.frame(ability)

modCat <- '
reason =~ reason.4 + reason.16 + reason.17 + reason.19
letter =~ letter.7 + letter.33 + letter.34 + letter.58
matrix =~ matrix.45 + matrix.46 + matrix.47 + matrix.55
rotate =~ rotate.3 + rotate.4 + rotate.6 + rotate.8
'

fitCat <- cfa(modCat, data = ability, std.lv = TRUE, 
              ordered = names(ability))
summary(fitCat)
