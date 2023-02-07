# R function description
These R functions are utilised for evaluating clinical prediction performance measures for binary outcomes for complete data (or single imputed data using median or mean imputation) using a logistic regression model. Also computed (on request) are optimism adjusted performance measures obtained via penalised (regularised) logistic regression model, with their corresponding Bias Corrected accelarated (BCa) percentile confidence intervals following case resampling (blocked) bootstrap.


# Required R packages
The entire "cprlogit_boot_function.R" should be saved. For the code to run, the following R packages should be Installed (and loaded)Install (and load) the following
 
 * tidyverse
 * glmnet
 * InformationValue
 * stats
 * rsample
 * pROC
 * fmsb
 * coxed
 
 You may find it easier to mass install all the R packages using the following R code.

```{r eval = FALSE, echo = FALSE}

lib.load.list <- c("tidyverse",
                   "glmnet",
                   "InformationValue",
                   "stats",
                   "rsample",
                   "pROC",
                   "fmsb",
                   "coxed")


lapply(lib.load.list,
       character.only=TRUE,
       library)
```

 
# User input

Head to head comparison using pairwise deleted data may b performed by running "logit.penalised.comparison()" function. Below is a demonstration using data (not shared) where three clinical variables (oxygen saturation, age and sex) and an outcome (met.end.point.updated) are used. Lasso regression is fitted using the R package glmnet  with alpha=1 (ridge =0 and elastic net any value between 0 and 1). The $\lambda$ value is the regularisation parameter estimated using cross validation. Please visit [link] (https://glmnet.stanford.edu/articles/glmnet.html) for more information.

```{r eval = FALSE, echo = FALSE}

################################################################################
############ FIT REGULARISED AND LOGICTIC REGRESSION
################################################################################
# Lasso (L1) regression

Need.o2_x <- data.matrix(test_data.use%>%
                         dplyr::select(vsoxy, ageyr, sex, CRP))

Need.o2_y <- data.matrix(test_data.use%>%
                        dplyr::select(met.end.point.updated))


# Get minimum value of lambda

cvfit_Need.o2 <-  glmnet::cv.glmnet(x=Need.o2_x,
                                   y=Need.o2_y,
                                   type.measure ="class",
                                   alpha=1,
                                   nfolds=10)
cvfit_Need.o2$lambda.min

# Run lasso binary regression
lasso_Need.o2_mod <-  glmnet::glmnet(x=Need.o2_x,
                                    y=Need.o2_y,
                                    alpha = 1, 
                                    lambda=cvfit_Need.o2$lambda.min,
                                    family = "binomial")

# Run logistic regression
logistic_Need.o2_mod <- stats::glm(met.end.point.updated~ vsoxy + ageyr + sex + CRP,
                                   data=test_data.use,
                                   family=binomial(link="logit")) 


summary(logistic_Need.o2_mod)


# Head to head comparison between logistic and regularised regression model

logit.penalised.comparison(logit.model=logistic_Need.o2_mod,
                           reg.model=lasso_Need.o2_mod,
                           lambda.min=cvfit_Need.o2$lambda.min,
                           data=test_data.use,
                           outcome="met.end.point.updated",
                           round.off.digits=4,
                           prob.class=0.25)
```
Case resampling bootstrapping is performed using the "bootstrap.sample()" function. For blocked bootsrap, a grouping variable may be specified using "strata=group" argument in the function. The object from the function is utilised by the "glm.net.cc.fit()" function to generate a data frame of the performance measures for each of the bootstrap sample.  
 
The predictors (at least one for regularised regression) and outcome variable must be specified along with number bootstrap performed. More than one probability threshold for classification are permitted, with default of 0.5 used if none is provided. For the models, the following are available - 1 (regularised/ penalised logistic regression with alpha.param set to 0=ridge,  1=lasso and any value in-between being elastic net regression), 2 (logistic regression) and 3 (both regularised and logistic).

```{r eval = FALSE, echo = FALSE}

################################################################################
#################FIT REGULARISED AND LOGICTIC REGRESSION WITH BCa CIs
################################################################################
# Bootstrap data
block.boot.cc.samples <- bootstrap.samples(data=test_data.use,
                                           samples=100)

# Bootstrap CI
tests.reg <- glm.net.cc.fit(data.boot.complete=block.boot.cc.samples,
                            predictors=c("vsoxy", "ageyr", "sex", "CRP"),
                            outcome="met.end.point.updated",
                            data.set.no=100,
                            thres.prob.classifier=c(0.25,0.3),
                            model.options=2,
                            alpha.param=1,
                            seed.input=7879690)
                            
```

The different performance measures of overall model performance, calibration and discrimination are provided for the different models specified.

```{r eval = FALSE, echo = FALSE}

tests.reg.summary <- parameter_estimates.cc(results=tests.reg,
                                            alpha=0.05)

tests.reg.summary
 