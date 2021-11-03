#Importing relevant libraries
library(car)
library(MASS)
library(robustbase)
library(MPV)

#Importing data
data.full = read.table(file.choose(),header=T)

#Dividing data into training and test sets
set.seed(0767950)
d.test <- sample(1:dim(data.full)[1], 200 )
data.test <- data.full[d.test, ]
data.training <- data.full[-d.test, ]

#1) EXPLORATORY DATA ANALYSIS
sum(is.na(data.training)) #0. No missing values.
summary(data.training)
attach(data.training)
boxplot(SWI, SWF, temperature, size, management, duration)
detach(data.training)
data.training.standardized <- data.frame(scale(data.training))
summary(data.training.standardized)
attach(data.training.standardized)
boxplot(SWI, SWF, temperature, size, management, duration)
detach(data.training.standardized)
#Standardized data is examined.
#No significant problem with outliers and skewness.

round(cor(data.training), 4)
#High positive correlations - SWI&SWF and temperature&duration

pairs(data.training) #temperature&duration pair shows some interesting properties.
plot(duration~temperature, data=data.training)
abline(v=25, col="blue")
#Irregularity on duration, duration variable is discarded for parametric analyses.

#2) FITTING A FIRST-ORDER LINEAR REGRESSION MODEL

#(a) Model fit
lm1 <- lm(SWI~SWF+temperature+size+management, data=data.training)
summary(lm1) #size isn't a significant variable.
par(mfrow = c(2,2))
plot(lm1) #There may be some problems.
par(mfrow = c(1,1))

#(b) Gauss-Markov conditions
gauss_markov <- function(linear_model) {
        #First Gauss-Markov condition. H0: E[residuals] = 0.
        gm1 <- t.test(linear_model$residuals, mu=0)
        #Second Gauss-Markov condition. H0: Homoscedasticity.
        gm2 <- ncvTest(linear_model)
        #Third Gauss-Markov condition. H0: Residuals are uncorrelated.
        gm3 <- durbinWatsonTest(linear_model)
        #All tests should have p-values above 0.05.
        results <- c(gm1$p.value > 0.05, gm2$p > 0.05, gm3$p > 0.05)
        if(sum(results) == 3) {
                return("All Gauss-Markov conditions are satisfied.")
                } else {
                print("Not all Gauss-Markov conditions are satisfied.")
                if(results[1] == F){
                        print("Residuals are not zero, on average.")
                        }
                if(results[2] == F){
                        print("Residuals are not homoscedastic.")
                        }
                if(results[3] == F){
                        print("Residuals are not uncorrelated.")
                        }
        }
}
gauss_markov(lm1) #Conditions aren't satisfied.

#(c) Multicollinearity
vifs <- vif(lm1)
sum(vifs > 10) #0
#No multicollinearity according to variation inflation factors.

rxx <- cor(data.training[,-c(1,6)])
eigen_rxx <- eigen(rxx)$values
condition <- sqrt(eigen_rxx[1]/eigen_rxx[2:4])
sum(condition > 30) #0
#No multicollinearity according to condition numbers.

#(d) Influential outliers
#Studentized residuals
lm1_studres <- studres(lm1)
plot(lm1_studres)
abline(h=c(-2.5,2.5))
#Only one influential point.

#DFFITS values
p=4; n=200
sum(abs(dffits(lm1)) > 2*sqrt(p/n))
#15 influential points.

#Cook's distances
sum(cooks.distance(lm1) > 1)
#No influential points.

#DFBETAS values
dfbetas <- (abs(dfbetas(lm1)) > 2/sqrt(n))[,2:5]
sum(rowSums(dfbetas) > 1)
#15 influential points.

#General influence table
influence.measures(lm1)
#3 influential points. (1,27,43)

#Robust diagnostics
lts1 <- ltsReg(SWI~SWF+temperature+size+management, data=data.training)
summary(lts1)
plot(lts1, which="rdiag")
#7 good leverage points and 1 vertical outlier.
#Vertical outlier is obs. 355.

#3) FINDING A BETTER MODEL
#Backward elimination
step_lm_backward <- stepAIC(lm1, direction = "back")
summary(step_lm_backward)
#Forward selection
step_lm_forward <- stepAIC(lm(SWI~1,data=data.training), direction = "forward", scope=formula(lm1))
summary(step_lm_forward)
#Stepwise regression
step_lm_stepwise <- stepAIC(lm(SWI~1,data=data.training), direction = "both", scope=formula(lm1))
summary(step_lm_stepwise)
#All 3 give the same model. (SWI ~ SWF + temperature + management)

#Updating the model
lm2 <- step_lm_stepwise
summary(lm2)
gauss_markov(lm2) #Conditions aren't satisfied.

#Creating new models using transformation
lm3 <- lm(sqrt(SWI)~SWF+temperature+management, data=data.training) #Square-root transformation
summary(lm3)
gauss_markov(lm3) #Conditions are satisfied.

lm4 <- lm(log(SWI)~SWF+temperature+management, data=data.training) #Logarithmic transformation
summary(lm4)
gauss_markov(lm4) #Conditions aren't satisfied.

lm_boxcox <- boxcox(lm1, lambda = seq(-2,2,0.001), plotit = TRUE) #Finding lambda for Box-Cox transformation
best_lambda <- lm_boxcox$x[which(lm_boxcox$y == max(lm_boxcox$y))]
best_lambda #0.645
range(lm_boxcox$x[lm_boxcox$y > max(lm_boxcox$y)-qchisq(0.95,1)/2]) #Zero is outside the interval. Lambda can be used.
data.training$boxcoxSWI <- (data.training$SWI^best_lambda-1)/best_lambda #Box-Cox transformation
lm5 <- lm(boxcoxSWI ~ SWF+temperature+management, data=data.training)
summary(lm5)
gauss_markov(lm5) #Conditions are satisfied.

lm6 <- lm(SWI ~ SWF + I(SWF^2)+temperature+management, data=data.training) #Square transformation
summary(lm6)
gauss_markov(lm6) #Conditions are satisfied.

lm7 <- lm(sqrt(SWI) ~ SWF + I(SWF^2)+temperature+management, data=data.training)
summary(lm7)
gauss_markov(lm7) #Conditions are satisfied.

lm8 <- lm(boxcoxSWI ~ SWF + I(SWF^2)+temperature+management, data=data.training)
summary(lm8)
gauss_markov(lm8) #Conditions are satisfied.

#Weighted regression models

stdev <- lm(abs(residuals(lm3))~SWF+temperature+management, data=data.training)
w <- 1/stdev$fitted.values^2 #Weights
lm3w <- lm(sqrt(SWI)~SWF+temperature+management, data=data.training, weights=w)
summary(lm3w)
gauss_markov(lm3w) #Conditions are satisfied.

stdev <- lm(abs(residuals(lm5))~SWF+temperature+management, data=data.training)
w <- 1/stdev$fitted.values^2 #Weights
lm5w <- lm(boxcoxSWI ~ SWF+temperature+management, data=data.training, weights=w)
summary(lm5w)
gauss_markov(lm5w) #Conditions are satisfied.

stdev <- lm(abs(residuals(lm6))~SWF+I(SWF^2)+temperature+management, data=data.training)
w <- 1/stdev$fitted.values^2 #Weights
lm6w <- lm(SWI ~ SWF+I(SWF^2)+temperature+management, data=data.training, weights=w)
summary(lm6w)
gauss_markov(lm6w) #Conditions are satisfied.

stdev <- lm(abs(residuals(lm7))~SWF+I(SWF^2)+temperature+management, data=data.training)
w <- 1/stdev$fitted.values^2 #Weights
lm7w <- lm(sqrt(SWI) ~ SWF+I(SWF^2)+temperature+management, data=data.training, weights=w)
summary(lm7w)
gauss_markov(lm7w) #Conditions are satisfied.

stdev <- lm(abs(residuals(lm8))~SWF+I(SWF^2)+temperature+management, data=data.training)
w <- 1/stdev$fitted.values^2 #Weights
lm8w <- lm(boxcoxSWI ~ SWF+I(SWF^2)+temperature+management, data=data.training, weights=w)
summary(lm8w)
gauss_markov(lm8w) #Conditions are satisfied.

#TEN MODELS TO COMPARE
AIC(lm3); AIC(lm5); AIC(lm6); AIC(lm7); AIC(lm8) #Unweighted AIC
AIC(lm3w); AIC(lm5w); AIC(lm6w); AIC(lm7w); AIC(lm8w) #Weighted AIC
PRESS(lm3);PRESS(lm5); PRESS(lm6); PRESS(lm7); PRESS(lm8) #Unweighted PRESS
PRESS(lm3w);PRESS(lm5w); PRESS(lm6w); PRESS(lm7w); PRESS(lm8w) #Weighted PRESS

#CONCLUSION: BEST MODELS ARE lm3 & lm7.
#lm3: sqrt(SWI)~I(SWF^2)+temperature+management
#lm7: sqrt(SWI)~SWF+temperature+management
summary(lm3)
summary(lm7)

#4) FITTING MODELS TO VALIDATION DATA
lm3_val <- lm(sqrt(SWI)~SWF+temperature+size+management, data=data.test)
summary(lm3_val)
gauss_markov(lm3_val)

lm7_val <- lm(sqrt(SWI)~SWF+I(SWF^2)+temperature+management, data=data.test)
summary(lm7_val)
gauss_markov(lm7_val)

#Performance measures
AIC(lm3_val); AIC(lm7_val)
PRESS(lm3_val); PRESS(lm7_val)

#CONCLUSION: lm7 IS (SLIGHTLY) BETTER.

#5) FITTING THE ULTIMATE MODEL ON THE FULL DATA
lm7_full <- lm(sqrt(SWI)~SWF+I(SWF^2)+temperature+management, data=data.full)
summary(lm7_full)
gauss_markov(lm7_full)
AIC(lm7_full); PRESS(lm7_full)

# Model combinations for partial contributions
lm7_w_swf <- lm(sqrt(SWI)~SWF, data=data.full); summary(lm7_w_swf)
lm7_w_swf2 <- lm(sqrt(SWI)~I(SWF^2), data=data.full); summary(lm7_w_swf2)
lm7_w_temp <- lm(sqrt(SWI)~temperature, data=data.full); summary(lm7_w_temp)
lm7_w_mgmt <- lm(sqrt(SWI)~management, data=data.full); summary(lm7_w_mgmt)
lm7_ww1 <- (lm(sqrt(SWI)~SWF+I(SWF^2), data=data.full)); summary(lm7_ww1)
lm7_ww2 <- (lm(sqrt(SWI)~SWF+temperature, data=data.full)); summary(lm7_ww2)
lm7_ww3 <- (lm(sqrt(SWI)~SWF+management, data=data.full)); summary(lm7_ww3)
lm7_ww4 <- (lm(sqrt(SWI)~I(SWF^2)+temperature, data=data.full)); summary(lm7_ww4)
lm7_ww5 <- (lm(sqrt(SWI)~I(SWF^2)+management, data=data.full)); summary(lm7_ww5)
lm7_ww6 <- (lm(sqrt(SWI)~temperature+management, data=data.full)); summary(lm7_ww6)
lm7_wo_swf <- lm(sqrt(SWI)~I(SWF^2)+temperature+management, data=data.full); summary(lm7_wo_swf)
lm7_wo_swf2 <- lm(sqrt(SWI)~SWF+temperature+management, data=data.full); summary(lm7_wo_swf2)
lm7_wo_temp <- lm(sqrt(SWI)~SWF+I(SWF^2)+management, data=data.full); summary(lm7_wo_temp)
lm7_wo_mgmt <- lm(sqrt(SWI)~SWF+I(SWF^2)+temperature, data=data.full); summary(lm7_wo_mgmt)

#CONCLUSION: Variable importance is SWF^2 > SWF > temperature > management

#6) NONPARAMETRIC AND QUADRATIC MODELS

#(a) Non-parametric model
attach(data.training)
nonp1 <- loess(duration~temperature, span=0.25, degree=1); nonp1
nonp2 <- loess(duration~temperature, span=0.5, degree=1); nonp2
nonp3 <- loess(duration~temperature, span=0.75, degree=1); nonp3
nonp4 <- loess(duration~temperature, span=0.25, degree=2); nonp4
nonp5 <- loess(duration~temperature, span=0.5, degree=2); nonp5
nonp6 <- loess(duration~temperature, span=0.75, degree=2); nonp6

#Comparing with the null model
anova(nonp1, update(nonp1,span=1))
anova(nonp2, update(nonp2,span=1))
anova(nonp3, update(nonp3,span=1))
anova(nonp4, update(nonp4,span=1))
anova(nonp5, update(nonp5,span=1))
anova(nonp6, update(nonp6,span=1))
#All models are better than the null model.

#MSE and MAPE
mean(nonp1$residuals^2); mean(abs(nonp1$residuals)/duration)
mean(nonp2$residuals^2); mean(abs(nonp2$residuals)/duration)
mean(nonp3$residuals^2); mean(abs(nonp3$residuals)/duration)
mean(nonp4$residuals^2); mean(abs(nonp4$residuals)/duration)
mean(nonp5$residuals^2); mean(abs(nonp5$residuals)/duration)
mean(nonp6$residuals^2); mean(abs(nonp6$residuals)/duration)

#CONCLUSION: BEST MODEL IS nonp4.

#(b) Quadratic linear model
lmq1 <- lm(duration~I(temperature^2), data=data.training); summary(lmq1); gauss_markov(lmq1)
lmq2 <- lm(duration~temperature+I(temperature^2), data=data.training); summary(lmq2); gauss_markov(lmq2)
lmq3 <- lm(duration~I(temperature^3), data=data.training); summary(lmq3); gauss_markov(lmq3)
lmq4 <- lm(duration~temperature+I(temperature^3), data=data.training); summary(lmq4); gauss_markov(lmq4)
lmq5 <- lm(duration~I(temperature^2)+I(temperature^3), data=data.training); summary(lmq5); gauss_markov(lmq5)
lmq6 <- lm(duration~temperature+I(temperature^2)+I(temperature^3), data=data.training); summary(lmq6); gauss_markov(lmq6)
#Performance criteria of the quadratic models
AIC(lmq1);AIC(lmq2);AIC(lmq3);AIC(lmq4);AIC(lmq5);AIC(lmq6)
PRESS(lmq1);PRESS(lmq2);PRESS(lmq3);PRESS(lmq4);PRESS(lmq5);PRESS(lmq6)

#CONCLUSION: THE BEST MODEL IS lmq6.

#(c) Plot and fits
plot(duration~temperature)
abline(v=25, lty=3)
lines(loess.smooth(temperature, duration, span=0.25, degree=2), col="blue") #Non-parametric fit
lines(temperature, 24.132803-1.495634*temperature+0.191738*temperature^2-0.004434*temperature^3, col="red") #Quadratic fit
legend(4,42, legend=c("Non-parametric", "Quadratic"), col=c("blue","red"), lty=1:2)

#Checking the validity of the non-parametric model
nonp_null <- update(nonp4, span=1)
scatter.smooth(temperature, residuals(nonp4), span=0.25, degree=2) #No problem. Around 0.
scatter.smooth(temperature, residuals(nonp4), span=1, degree=1) #No problem. Around 0.
qqnorm(residuals(nonp4)); qqline(residuals(nonp4)) #Seems normal.
shapiro.test(residuals(nonp4)) #p-value 0.8635
abline(h=0, lty=2)

#(d) Non-parametric vs. Quadratic comparison test
#Res.Std., MSE, MAD, MAPE
nonp4; summary(lmq6)
mean(nonp4$residuals^2); mean(lmq6$residuals^2)
mean(abs(nonp4$residuals)); mean(abs(lmq6$residuals))
mean(abs(nonp4$residuals)/duration); mean(abs(lmq6$residuals)/duration)
