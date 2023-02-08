# R function description
These R functions are utilised for evaluating clinical prediction performance measures for binary outcomes for Multiply Imputed data using a logistic regression model. Also computed (on request) are optimism adjusted performance measures obtained via penalised (regularised) logistic regression model, with their corresponding Bias Corrected accelarated (BCa) percentile confidence intervals.


# Required R packages
The entire "cprlogit_mi_no_boot_function.R" should be saved. For the code to run, the following R packages should be Installed (and loaded) 

* tidyverse
* glmnet
* InformationValue
* stats
* qpcR
* pROC
* fmsb
* coxed

You may find it easier to mass install  (and) load all the R packages using the pacman R package.
 
 ```{r eval = FALSE, echo = FALSE}
 
 ### Not run
 
 if(!(require("pacman"))){install.packages("pacman")}
 
 pacman::p_load(tidyverse,
                glmnet,
                InformationValue,
                stats,
                qpcR,
                pROC,
                fmsb,
                coxed)


``` 
Note that  *InformationValue* R package is not supported in CRAN and should thus be installed from the archive  [here](https://cran.r-project.org/src/contrib/Archive/InformationValue/).
If the packages are already installed, you may mass load them without using *pacman* R package as follows

```{r eval = FALSE, echo = FALSE}

lib.load.list <- c("tidyverse",
                   "glmnet",
                   "InformationValue",
                   "stats",
                   "qpcR",
                   "pROC",
                   "fmsb",
                   "coxed")


lapply(lib.load.list,
       character.only=TRUE,
       library)

                                       

```

# User input

Multiply imputed datasets from R package "mice" converted to long format are required for the evaluation of the clinical prediction model. The predictors (at least one for regularised regression) and outcome variable must be specified along with number of imputations. More than one probability threshold for classification are permitted, with default of 0.5 used if none is provided. For the models, the following are available - 1 (regularised/ penalised logistic regression with alpha.param set to 0=ridge,  1=lasso and any value in-between for elastic net regression), 2 (logistic regression) and 3 (both regularised and logistic). Please visit [link](https://glmnet.stanford.edu/articles/glmnet.html) for more information on the model types supported by *glmnet* R package.

```{r eval = FALSE, echo = FALSE}

test_results <- glm.net.no_boot.m.fit(data.boot.mi_complete=data, # multiply imputed dataframe 
                                      predictors=c("predictor1","predictor2",....),
                                      outcome="outcome_variable",
                                      data.set.no=5, # number of imputations
                                      thres.prob.classifier=c(0.05,0.075), 
                                      model.options=1,
                                      alpha.param=1,
                                      seed.input=8754654)
                                       

```

The different performance measures of overall model performance, calibration and discrimination are provided for the different models specified.

```{r eval = FALSE, echo = FALSE}
test_results_summary <- parameter_estimates.mi_no_boot(results=test_results,
                                                       alpha=0.05)
test_results_summary

```