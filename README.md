## Time series models
* [General Info](#general-info)
* [Technologies](#technologies)
* [Steps](#steps)

## General Info

This project is use time series dataset to model and describe relationships among the abundance of targeted phytoplankton predator, different nutrient types (i.e., phosphate, nitrogen) and its prey in Chesapeake Bay. Here this project aims to advance knowledge of phytoplankton predator bloom ecology in three ways:
	
i) Use cross-correlation functions to model phytoplankton predator abundance (i.e., the dependent variables) as a potential functions of past lags of phytoplankton predator abundance, and current and past lags of nutrients and prey (i.e., the independent variables).
ii)Explore lagged relationships along with non-lagged relationships in a full transfer function model for identifying variables that might be useful predictor of phytoplankton bloom.
iii)Evaluate interactions among predator, nutrients and preys simultaneously.


The work is published at Harmful Algae Journal. 2018 Mar;73:110-118. DOI: 10.1016/j.hal.2018.02.002


## Technologies 

Project is created with:
* MySQL 
* R Studio

## Steps of statistical procedures

STEP 1- Prepare all the time series data (add the P variable(s) and the time series for targeted phytoplankton species). Here I use MySQL as a tool to preprocess the CB datasets (such as data cleaning and transformation; https://datahub.chesapeakebay.net/). 

STEP 2- Check trend presence in each variable (e.g., use ACF plots [linear decay might be a signal of a trend] or WAVK test [small p-values tell about systematic change that could be a trend or a seasonal/cyclical component]).
If no obvious trend found, proceed with the following steps.

STEP 3.1- Select initial variables and form of relationship for the regression, in replace to the current CCF plots, use 

library(astsa)
lag2.plot(Xvariable, Yvariable, Klag)

STEP 3.2- Explore the K lags of all your X-variables (change X-variable one by one). This function also reports the linear correlation coefficients in the top right corner that should be the same/similar as calculated by CCF. The difference with CCF is that i) you will see the actual form of relationship from the plots; ii) you will not see the future lags, which is okay, since you cannot know future values of the X-variables anyway, without predicting them separately and introducing even more uncertainty into the model for you Y-variable.

STEP4- Combine the variables (Y, lagged X, and non-lagged X that you found relevant) in a matrix (e.g., call it "D") and check correlations again, check if there any multicollinearity, with |r|>0.8, use

library(GGally)
library(Hmisc)
ggpairs(D, axisLabels="internal", lab.cex=0.8)
rcorr(as.matrix(D), type="pearson")

STEP 5 - Use all the variables and estimate a 'full model':

model =lm(Y~all pre-selected variables)

STEP 6-Use 

library(MASS)
stepAIC(model)
to automatically shrink the model by selecting variables based on AIC.

STEP 7- Do diagnostics checks for the new model's residuals (ACF, PACF, Q-Q plots, etc.). If needed, you can incorporate ARMA structure at this stage (based on ACF and PACF), as well as a different distribution (based on the Q-Q and residuals vs. time plots). 