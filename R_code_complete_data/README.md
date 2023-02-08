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
 
 You may find it easier to mass install  (and) load all the R packages using the pacman R package.
 
 ```{r eval = FALSE, echo = FALSE}
 
 ### Not run
 
 if(!(require("pacman"))){install.packages("pacman")}
 
 pacman::p_load(tidyverse,
                glmnet,
                stats,
                rsample,
                pROC,
                fmsb,
                coxed)


``` 
Note that  *InformationValue* R package is not supported in CRAN and should thus be installed from the archive here [archive](https://cran.r-project.org/src/contrib/Archive/InformationValue/). If the packages are already installed, you may mass load them without using pacman R package as follows

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

Head to head comparison using pairwise deleted data may be performed by running "logit.penalised.comparison()" function. Below is a demonstration using data (not shared) where three clinical variables (oxygen saturation, age and sex) and an outcome (met.end.point.updated) are used. Lasso regression is fitted using the R package glmnet  with alpha=1 (alternatively, set ridge = 0 or  any value between 0 and 1 for elastic net). The $\lambda$ value is the regularisation parameter estimated using cross validation. Please visit [link](https://glmnet.stanford.edu/articles/glmnet.html) for more information.

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
                                   
```


```{r eval = FALSE, echo = FALSE}
summary(logistic_Need.o2_mod)


# Head to head comparison between logistic and regularised regression model

logit.penalised.comparison(logit.model=logistic_Need.o2_mod,
                           reg.model=lasso_Need.o2_mod,
                           lambda.min=cvfit_Need.o2$lambda.min,
                           data=test_data.use,
                           outcome="met.end.point.updated",
                           round.off.digits=4,
                           prob.class=0.25)

$parameters
            logistic.model penalised.logistic.model shrinkage
(Intercept)        40.7473                  37.2067        NA
vsoxy              -0.4368                  -0.3992    0.9140
ageyr               0.0012                   0.0000    0.0000
sex                -0.0559                   0.0000    0.0000
CRP                 0.0041                   0.0035    0.8582

$overall.measures.of.performance
              logit.model penalised.logit.model
R2.Nagelkerke      0.1522                0.0986
Brier.score        0.2919                0.2895

$measures.of.discrimination
    Statistic logit.model penalised.logit.model
1         auc      0.7268                0.7279
2 concordance      0.7268                0.7279
3 sensitivity      0.6292                0.6180
4 specificity      0.7295                0.7356
5 brier.score      0.2919                0.2895
6         ppv      0.3862                0.3873
7         npv      0.8791                0.8768
8         LR+      2.3260                2.3369
9         LR-      0.5083                0.5194

$measures.of.calibration
               Statistic logit.model penalised.logit.model
1 n.predicted/n.outcomes      1.6292                1.5955
2          cox.intercept      0.0000                0.1382
3              cox.slope      1.0000                1.1188

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

$logistic.regression
$logistic.regression[[1]]
                       PT_0.25_Mean BCa.LL_2.5% BCa.UL_97.5%
Intercept                   41.0862     28.1488      59.1972
Beta_vsoxy                  -0.4408     -0.6266      -0.3121
Beta_ageyr                   0.0020     -0.0159       0.0157
Beta_sex                    -0.0969     -0.6453       0.3426
Beta_CRP                     0.0040      0.0011       0.0062
auc                          0.7331      0.6802       0.7878
mse                         -0.0446     -0.1180       0.0259
R2.Nagelkerke                0.1608      0.1054       0.2619
brier.score                  0.2810      0.2273       0.3325
concordance                  0.7331      0.6802       0.7878
n.predicted/n.outcomes       1.5280      1.2820       1.7357
cox.int                      0.0000      0.0000       0.0000
cox.slope                    1.0000      1.0000       1.0000
sensitivity                  0.5987      0.4609       0.7363
specificity                  0.7496      0.6800       0.8285
ppv                          0.3913      0.3300       0.4463
npv                          0.8756      0.8447       0.9101
LR+                          2.4281      1.8282       3.1908
LR-                          0.5342      0.3667       0.7053

$logistic.regression[[2]]
                       PT_0.3_Mean BCa.LL_2.5% BCa.UL_97.5%
Intercept                  41.0862     28.1488      59.1972
Beta_vsoxy                 -0.4408     -0.6266      -0.3121
Beta_ageyr                  0.0020     -0.0159       0.0157
Beta_sex                   -0.0969     -0.6453       0.3426
Beta_CRP                    0.0040      0.0011       0.0062
auc                         0.7331      0.6802       0.7878
mse                        -0.0446     -0.1180       0.0259
R2.Nagelkerke               0.1608      0.1054       0.2619
brier.score                 0.2502      0.1980       0.2895
concordance                 0.7331      0.6802       0.7878
n.predicted/n.outcomes      1.1122      0.8277       1.3395
cox.int                     0.0000      0.0000       0.0000
cox.slope                   1.0000      1.0000       1.0000
sensitivity                 0.4633      0.2778       0.6524
specificity                 0.8248      0.7711       0.8847
ppv                         0.4123      0.2822       0.4893
npv                         0.8529      0.8070       0.8932
LR+                         2.6749      1.6222       3.8026
LR-                         0.6493      0.4341       0.8602

```