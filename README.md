# Time-Series-Analysis

1- prepare all the time series data (add the P variable(s) and the time series for unidentified species)
2- check trend presence in each variable (e.g., use ACF plots [linear decay might be a signal of a trend] or WAVK test [small p-values tell about systematic change that could be a trend or a seasonal/cyclical component]). If you notice trend presence somewhere let me know (send the graphs and test results) and we can discuss what to do.

If no obvious trend found, proceed with the following steps.

3- to select initial variables and form of relationship for the regression, in replace to the current CCF plots, use 
library(astsa)
lag2.plot(Xvariable, Yvariable, Klag)

to explore the K lags of all your X-variables (change X-variable one by one). This function also reports the linear correlation coefficients in the top right corner that should be the same/similar as calculated by CCF. The difference with CCF is that i) you will see the actual form of relationship from the plots; ii) you will not see the future lags, which is okay, since you cannot know future values of the X-variables anyway, without predicting them separately and introducing even more uncertainty into the model for you Y-variable.
4- Combine the variables (Y, lagged X, and non-lagged X that you found relevant) in a matrix (e.g., call it "D") and check correlations again, check if there any multicollinearity, with |r|>0.8
library(GGally)
library(Hmisc)
ggpairs(D, axisLabels="internal", lab.cex=0.8)
rcorr(as.matrix(D), type="pearson")
5 - Use all the variables and estimate a 'full model':

model =lm(Y~all pre-selected variables)

6-Use 
library(MASS)
stepAIC(model)
to automatically shrink the model by selecting variables based on AIC.

7- Do diagnostics checks for the new model's residuals (ACF, PACF, Q-Q plots, etc.). If needed, you can incorporate ARMA structure at this stage (based on ACF and PACF), as well as a different distribution (based on the Q-Q and residuals vs. time plots). 