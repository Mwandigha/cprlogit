################################################################################
################################################################################
### USEFUL OWN WRITTEN FUNCTIONS
################################################################################
################################################################################
### BRIER SCORE

##
#Journal of the American Medical Informatics Association, 27(4), 2020, 621–633
#doi: 10.1093/jamia/ocz228
#Advance Access Publication Date: 27 February 2020

# Equation 2

brier.score <- function(outcome,
                        pred,
                        threshold)
{
  pred.outcome <- ifelse(pred>=threshold,1,0)
  
  bs <- (sum((outcome-pred.outcome)^2))/(length(outcome))
  
  return(bs)
  
}
## COX INTERCEPT AND SLOPE
#  Equation 7
# NB: The values obtained with function match those from R package rms

calibration.int.slope <- function(outcome,
                                  pred)
{
  #####
  soln.model <- glm(formula=outcome~log(pred/(1-pred)),
                    family=binomial(link="logit")) 
  
  #####
  soln.model.df <- as.data.frame(soln.model$coefficients)
  
  int.model <- soln.model.df["(Intercept)",]
  slope.model <- soln.model.df["log(pred/(1 - pred))", ]
  #####
  output <- c(int.model,slope.model)
  
  return(output)
  
}

### N.OUTCOME vs. N.EXPECTED RATIO
outcome_expected.scaled.diff <- function(outcome,
                                         pred,
                                         threshold.prob)
{
  #####
  n.outcomes <- sum(outcome)
  #####
  expected <- ifelse((pred >= threshold.prob),
                     1,
                     0)
  
  n.expected <- sum(expected)
  #####
  outcome_expected.ratio <- (n.expected/n.outcomes)
  
  return(outcome_expected.ratio)
  
}

################################################################################
################################################################################
### R FUNCTION TO COMPUTE OPTIMISM ADJUSTED PEROMANCE MEASURES WITH MI DATA


################################################################################
################################################################################

regularisation.no_boot.mi.logit <- function(data.boot.mi_complete,
                                            predictors,
                                            outcome,
                                            data.set.no,
                                            thres.prob.classifier=NULL,
                                            model.options=NULL,
                                            alpha.param,
                                            seed.input)
{
  ################################################################################
  #  (PENALISED)  MODEL WITH LOGIT LINK
  ################################################################################
  
  ################################################################################
  
  set.seed(seed=seed.input)
  
  if(!is.null(thres.prob.classifier))
  {
    if(min(thres.prob.classifier) <0 | max(thres.prob.classifier) >1)
    { 
      stop("The probability threshold for classification MUST be between 0 and 1")
    }
    
  }
  
  ####
  if(is.null(model.options)){model.options <- 2}
  
  if(!(!is.null(model.options)&model.options%in%c(1,2,3)))
  {
    cat("Note: model options argument can assume only 1, 2 or 3","\n")
  }
  ###
  if(!(!is.null(model.options)&model.options%in%c(1,2,3))){cat("Note: model options argument defaulted to 2:(logistic model only)","\n")}
  
  if(!(!is.null(model.options)&model.options%in%c(1,2,3))){model.options <- 2}
  
  if(!is.null(model.options)&(model.options==1|model.options==3)&length(predictors)==1)
  {
    cat("Note: Regularisation not possible with one predictor. Model options defaulted to 2: (logistic model only)","\n")
  }
  
  if(!is.null(model.options)&(model.options==1|model.options==3)&length(predictors)==1)
  {
    model.options <- 2
  }
  
  
 
  
  ################################################################################
  ####### PRELIMINARY DATA PREPARATION
  ################################################################################
  dfs.mi_no_boot <- data.boot.mi_complete%>%
                      dplyr::group_split(.imp)
  
  
  # Dummy code ALL the variables if required
  
  pred.use <- predictors
  outcome.use <- outcome
  
  dfs.boot <- purrr::map(seq(data.set.no),
                         function(k)
                         {
                           na.omit(dfs.mi_no_boot[[k]]%>% 
                                     dplyr::select(all_of(pred.use),
                                                   all_of(outcome)))
                         })
  
  ################################################################################
  ####### Regularised regression only
  ################################################################################
  
  if(model.options==1)
  {
    pred.vars <- purrr::map(seq(data.set.no),
                            function(j)
                            {
                              data.matrix(dfs.boot[[j]] %>%
                                            dplyr::select(all_of(pred.use)))})
    
    outcome.var <- purrr::map(seq(data.set.no),
                              function(j)
                              {
                                as.matrix(dfs.boot[[j]] %>%
                                            dplyr::select(all_of(outcome.use)))})
    
    # For each data, find lambda that minimize the cross-validation prediction error rate. 
    
    
    cv.lasso.dev <- lapply(seq(data.set.no),
                           function(j)
                           {
                             glmnet::cv.glmnet(x=pred.vars[[j]],
                                               y=outcome.var[[j]], 
                                               alpha = alpha.param,
                                               nfolds=10,
                                               type.measure="deviance",
                                               family = "binomial")})
    
    best.lambda.dev <- lapply(seq(data.set.no),
                              function(k)
                              {
                                cv.lasso.dev[[k]]$lambda.min})
    
    
    # Fit the final models
    lasso.model <- lapply(seq(data.set.no),
                          function(s)
                          {
                            glmnet::glmnet(x=pred.vars[[s]],
                                           y=outcome.var[[s]], 
                                           alpha = alpha.param,
                                           family = "binomial",
                                           lambda = best.lambda.dev[[s]])})
    
    
    ##### COEFFICIENTS
    lasso.coefs <- lapply(seq(data.set.no),
                          function(s)
                          {
                            
                            predict(object = lasso.model[[s]],
                                    s = best.lambda.dev[[s]],
                                    type = "coefficients")
                          })
    
    
    
    
    lasso.coef.matrix <- lapply(seq(data.set.no),
                                function(s)
                                {
                                  as.matrix(lasso.coefs[[s]])
                                })
    
    matrix.coefs <- matrix(nrow=length(lasso.coef.matrix), 
                           ncol=length(lasso.coef.matrix[[1]]))
    for (j in seq(data.set.no))
    {
      
      matrix.coefs[j,] <- as.matrix(lasso.coefs[[j]])
      
    }
    
    
    lasso.coefficients <- data.frame(matrix.coefs)
    colnames(lasso.coefficients) <- c("Intercept",paste( "Beta",pred.use,sep="_"))
    lasso.coefficients$Imputation_no <- seq(data.set.no)
    
    lasso.coefficients.df <- lasso.coefficients[, c(ncol(lasso.coefficients),
                                                    seq(ncol(lasso.coefficients)-1))]
    
    
    ######### MEASURES OF PEFORMANCE
    gofstats.lasso <- lapply(seq(data.set.no),
                             function(s)
                             {
                               t(glmnet::assess.glmnet(object=lasso.model[[s]],
                                                       newx=pred.vars[[s]],
                                                       newy=outcome.var[[s]],
                                                       family="binomial",
                                                       s= best.lambda.dev[[s]]))
                             })
    
    
    matrix.gofs.lasso <- matrix(nrow=length(gofstats.lasso),
                                ncol=5)
    
    for (j in seq(data.set.no))
    {
      
      matrix.gofs.lasso[j,] <- as.matrix.data.frame(gofstats.lasso[[j]])
      
    }
    
    
    df.gof.lasso <- as.data.frame(matrix.gofs.lasso)
    colnames(df.gof.lasso) <- c("deviance","class","auc","mse","mae")
    
    # Nagalkerke's R2
    R2.Nagelkerke <- lapply(seq(data.set.no),
                            function(k)
                            {
                              lasso.model[[k]]$dev.ratio[which(lasso.model[[k]]$lambda==min(lasso.model[[k]]$lambda))]
                            })
    
    
    
    
    
    
    df.gof.lasso$R2.Nagelkerke <- as.vector(unlist(R2.Nagelkerke)) 
    df.gof.lasso$Imputation_no <- seq(data.set.no)
    
    
    
    #####
    
    df.gof.stats.lasso <- df.gof.lasso[, c(ncol(df.gof.lasso),
                                           seq(ncol(df.gof.lasso)-1))]
    
    
    df.lasso.results <- merge(x=lasso.coefficients.df,
                              y=df.gof.stats.lasso,
                              by="Imputation_no")
    
    df.lasso.results.output <- df.lasso.results %>%
                                   dplyr::select(-deviance,-class,-mae)
    
    
    # predictions
    preds.lasso <- lapply(seq(data.set.no),
                          function(k)
                          {
                            stats::predict(lasso.model[[k]], 
                                           newx = pred.vars[[k]], 
                                           type = 'response',
                                           s= best.lambda.dev[[k]])
                            
                            
                          })
    
    
    
    # optimal threshold probability for classification
    
    if(is.null(thres.prob.classifier))
    {
      optimal.lasso.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     InformationValue::optimalCutoff(dfs.boot[[k]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[k]])
                                     
                                   })
    }
    
    if(!is.null(thres.prob.classifier))
    {
      optimal.lasso.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     thres.prob.classifier
                                   })
    }
    
    
    # Brier score
    brier.score.lasso.lists <- purrr::map(seq(data.set.no),
                                          function(x)
                                          {
                                            brier.score(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                        pred=as.vector(preds.lasso[[x]]),
                                                        threshold=optimal.lasso.list[[x]])
                                            
                                          })
    
    
    
    # concordance  
    concordance.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::Concordance(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]])
                                       
                                     })
    
    # sensitivity
    sensitivity.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::sensitivity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]],
                                                                     optimal.lasso.list[[s]])
                                       
                                     })
    # outcome vs expected scaled difference
    outcome.vs.expected.ratio.lasso <- lapply(seq(data.set.no),
                                              function(x)
                                              {
                                                outcome_expected.scaled.diff(outcome=dfs.boot[[x]]%>%dplyr::select(all_of(outcome)),
                                                                             pred =preds.lasso[[x]],
                                                                             threshold.prob =optimal.lasso.list[[x]])
                                              })
    
    # cox intercept and slope
    cox.int.slope.lasso <- lapply(seq(data.set.no), 
                                  function(x)
                                  {
                                    
                                    calibration.int.slope(outcome=as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome))),
                                                          pred=preds.lasso[[x]])                                  
                                  })
    
    # specificity
    specificity.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::specificity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]],
                                                                     threshold=optimal.lasso.list[[s]])
                                       
                                     })
    
    
    
    # ppv 
    ppv.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::precision(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                           preds.lasso[[s]],
                                                           threshold=optimal.lasso.list[[s]])
                               
                             })
    
    # npv
    npv.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::npv(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                     preds.lasso[[s]],
                                                     threshold=optimal.lasso.list[[s]])
                               
                             })
    
    # plr
    
    plr.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               sensitivity.lasso.list[[s]]/(1-specificity.lasso.list[[s]])
                               
                             })
    
    # nlr
    
    nlr.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               (1-sensitivity.lasso.list[[s]])/specificity.lasso.list[[s]]
                               
                             })
    
    
    matrix.brier.score.lasso <- matrix(nrow=length(brier.score.lasso.lists))
    matrix.concordance.lasso <- matrix(nrow=length(concordance.lasso.list))
    matrix.n.O.vs.n.E.lasso <- matrix(nrow=length(outcome.vs.expected.ratio.lasso))
    matrix.cox.int.slope.lasso <- matrix(nrow=length(cox.int.slope.lasso),
                                         ncol=2)
    matrix.optimal.lasso <- matrix(nrow=length(optimal.lasso.list),
                                   ncol=length(optimal.lasso.list[[1]]))
    matrix.sensitivity.lasso <- matrix(nrow=length(sensitivity.lasso.list))
    matrix.specificity.lasso <- matrix(nrow=length(sensitivity.lasso.list))
    matrix.ppv.lasso <- matrix(nrow=length(ppv.lasso.list))
    matrix.npv.lasso <- matrix(nrow=length(npv.lasso.list))
    matrix.plr.lasso <- matrix(nrow=length(plr.lasso.list))
    matrix.nlr.lasso <- matrix(nrow=length(nlr.lasso.list))
    
    
    for (k in seq(data.set.no))
    {
      matrix.brier.score.lasso[k,] <- brier.score.lasso.lists[[k]]
      matrix.concordance.lasso[k,] <- concordance.lasso.list[[k]]$Concordance
      matrix.n.O.vs.n.E.lasso[k,] <- outcome.vs.expected.ratio.lasso[[k]]
      matrix.cox.int.slope.lasso[k,] <- cox.int.slope.lasso[[k]]
      matrix.optimal.lasso[k,] <- optimal.lasso.list[[k]]
      matrix.sensitivity.lasso[k,] <- sensitivity.lasso.list[[k]]
      matrix.specificity.lasso[k,] <- specificity.lasso.list[[k]]
      matrix.ppv.lasso[k,]<- ppv.lasso.list[[k]]
      matrix.npv.lasso[k,]<- npv.lasso.list[[k]]
      matrix.plr.lasso[k,]<- plr.lasso.list[[k]]
      matrix.nlr.lasso[k,]<- nlr.lasso.list[[k]]
    }
    
    perf.measures.lasso.df <- data.frame(seq(data.set.no),
                                         matrix.brier.score.lasso,
                                         matrix.concordance.lasso,
                                         matrix.n.O.vs.n.E.lasso,
                                         matrix.cox.int.slope.lasso,
                                         matrix.optimal.lasso,
                                         matrix.sensitivity.lasso,
                                         matrix.specificity.lasso,
                                         matrix.ppv.lasso,
                                         matrix.npv.lasso,
                                         matrix.plr.lasso,
                                         matrix.nlr.lasso)
    
    colnames(perf.measures.lasso.df) <-c("Imputation_no",
                                         "brier.score",
                                         "concordance",
                                         "n.predicted/n.outcomes",
                                         "cox.int",
                                         "cox.slope",
                                         "threshold.prob",
                                         "sensitivity",
                                         "specificity",
                                         "ppv",
                                         "npv",
                                         "LR+",
                                         "LR-")
    
    
    
    lasso.results.output <- merge(df.lasso.results.output,
                                  perf.measures.lasso.df,
                                  by="Imputation_no")
    
    
    lasso.results.output <- lasso.results.output%>%
                             dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                                dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
    
    return(list("regularised.regression"=dplyr::as_tibble(lasso.results.output),
                "model.option"=model.options))
    
    
  }
  
  ################################################################################
  ####### Logistic regression only
  ################################################################################
  
  if(model.options==2)
  {
    
    
    ################################################################################
    #  LOGISTIC MODEL WITH LOGIT LINK
    ################################################################################
    
    formula.model <- as.formula(paste(outcome.use, 
                                      paste(pred.use, 
                                            collapse=" + "),
                                      sep=" ~ "))
    
    logit.results.list <- lapply(seq(data.set.no),
                                 function(p)
                                 {
                                   glm(formula.model,
                                       family=binomial(link=logit),
                                       data=dfs.boot[[p]])
                                 })
    
    
    
    logit.results.matrix <- matrix(nrow=length(logit.results.list),
                                   ncol=length(logit.results.list[[1]]$coefficients))
    
    for (i in seq(data.set.no))
    {
      logit.results.matrix[i,] <- logit.results.list[[i]]$coefficients
    }
    
    logit.results.df <- as.data.frame(logit.results.matrix)
    colnames(logit.results.df) <- c("Intercept",paste( "Beta",pred.use,sep="_"))
    
    logit.results.df$Imputation_no <- seq(data.set.no)
    
    logit.results.df <- logit.results.df[,c(ncol(logit.results.df),
                                            seq(ncol(logit.results.df)-1))]
    
    
    # Performance measures
    
    if(!(require(pROC))){install.packages("pROC", dependencies=TRUE)}
    if(!(require(fmsb))){install.packages("fmsb", dependencies=TRUE)}
    # auc
    auc.logit <- lapply(seq(data.set.no),
                        function(s)
                        {
                          suppressMessages(pROC::auc(response=as.vector(as.matrix(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)))),
                                                     predictor=logit.results.list[[s]]$fitted.values))
                          #suppressMessages() prevents progress status from being printed
                        })
    
    
    
    
    # mse
    mse.logit <- lapply(seq(data.set.no),
                        function(s)
                        {
                          (sum(logit.results.list[[s]]$residuals)/
                             (nrow(dfs.boot[[s]])-length(pred.use)))
                        })
    
    # R2 Nagelkerke
    R2.Nagelkerke.logit.lists <- lapply(seq(data.set.no),
                                        function(s)
                                        {
                                          fmsb::NagelkerkeR2(logit.results.list[[s]])
                                        })
    
    
    # concordance  
    concordance.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::Concordance(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values)
                                       
                                     })
    
    # cox intercept and slope
    cox.int.slope.logit <- lapply(seq(data.set.no), 
                                  function(x)
                                  {
                                    
                                    calibration.int.slope(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                          pred=logit.results.list[[x]]$fitted.values)                                  
                                  })
    # threshold probability
    # optimal threshold probability for classification
    if(is.null(thres.prob.classifier))
    {
      optimal.logit.list <- lapply(seq(data.set.no),
                                   function(s)
                                   {
                                     InformationValue::optimalCutoff(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values)
                                   })
      
    }
    
    if(!is.null(thres.prob.classifier))
    {
      optimal.logit.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     thres.prob.classifier
                                   })
    }
    
    # Brier score
    
    brier.score.logit.lists <- purrr::map(seq(data.set.no),
                                          function(x)
                                          {
                                            brier.score(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                        pred=as.vector(logit.results.list[[x]]$fitted.values),
                                                        threshold=optimal.logit.list[[x]])
                                            
                                          })
    
    # sensitivity
    sensitivity.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::sensitivity(as.vector(as.matrix(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)))),
                                                                     as.vector(logit.results.list[[s]]$fitted.values),
                                                                     threshold=optimal.logit.list[[s]])
                                       
                                     })
    
    # specificity
    specificity.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::specificity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values,
                                                                     threshold=optimal.logit.list[[s]])
                                       
                                     })
    
    
    # outcome vs expected scaled difference
    outcome.vs.expected.ratio.logit <- lapply(seq(data.set.no),
                                              function(x)
                                              {
                                                outcome_expected.scaled.diff(outcome=dfs.boot[[x]]%>%dplyr::select(all_of(outcome)),
                                                                             pred =logit.results.list[[x]]$fitted.values,
                                                                             threshold.prob =optimal.logit.list[[x]])
                                              })
    
    
    # ppv 
    ppv.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::precision(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                           logit.results.list[[s]]$fitted.values,
                                                           threshold=optimal.logit.list[[s]])
                               
                             })
    
    # npv
    npv.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::npv(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                     logit.results.list[[s]]$fitted.values,
                                                     threshold=optimal.logit.list[[s]])
                               
                             })
    
    # plr
    
    plr.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               sensitivity.logit.list[[s]]/(1-specificity.logit.list[[s]])
                               
                             })
    
    # nlr
    
    nlr.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               (1-sensitivity.logit.list[[s]])/specificity.logit.list[[s]]
                               
                             })
    
    
    
    matrix.auc.logit <- matrix(nrow=length(auc.logit))
    matrix.mse.logit <- matrix(nrow=length(mse.logit))
    matrix.R2.Nagelkerke.logit <- matrix(nrow=length(R2.Nagelkerke.logit.lists))
    matrix.brier.score.logit <- matrix(nrow=length(brier.score.logit.lists))
    matrix.concordance.logit <- matrix(nrow=length(concordance.logit.list))
    matrix.n.O.vs.n.E.logit <- matrix(nrow=length(outcome.vs.expected.ratio.logit))
    matrix.cox.int.slope.logit <- matrix(nrow=length(cox.int.slope.logit),
                                         ncol=2)
    matrix.optimal.logit <- matrix(nrow=length(optimal.logit.list))
    matrix.sensitivity.logit <- matrix(nrow=length(sensitivity.logit.list))
    matrix.specificity.logit <- matrix(nrow=length(sensitivity.logit.list))
    matrix.ppv.logit <- matrix(nrow=length(ppv.logit.list))
    matrix.npv.logit <- matrix(nrow=length(npv.logit.list))
    matrix.plr.logit <- matrix(nrow=length(plr.logit.list))
    matrix.nlr.logit <- matrix(nrow=length(nlr.logit.list))
    
    
    for (k in seq(data.set.no))
    {
      matrix.auc.logit[k,] <- as.numeric(gsub("Area under the curve:",
                                              "", auc.logit[[k]]))
      matrix.mse.logit[k,] <- mse.logit[[k]]
      matrix.R2.Nagelkerke.logit[k,] <- R2.Nagelkerke.logit.lists[[k]]$R2
      matrix.brier.score.logit[k,] <- brier.score.logit.lists[[k]]
      matrix.concordance.logit[k,] <- concordance.logit.list[[k]]$Concordance
      matrix.n.O.vs.n.E.logit[k,] <- outcome.vs.expected.ratio.logit[[k]]
      matrix.cox.int.slope.logit[k,] <- cox.int.slope.logit[[k]]
      matrix.optimal.logit[k,] <- optimal.logit.list[[k]]
      matrix.sensitivity.logit[k,] <- sensitivity.logit.list[[k]]
      matrix.specificity.logit[k,] <- specificity.logit.list[[k]]
      matrix.ppv.logit[k,] <- ppv.logit.list[[k]]
      matrix.npv.logit[k,] <- npv.logit.list[[k]]
      matrix.plr.logit[k,]<- plr.logit.list[[k]]
      matrix.nlr.logit[k,]<- nlr.logit.list[[k]]
      
    }
    
    perf.measures.logit.df <- data.frame(seq(data.set.no),
                                         matrix.auc.logit,
                                         matrix.mse.logit,
                                         matrix.R2.Nagelkerke.logit,
                                         matrix.brier.score.logit,
                                         matrix.concordance.logit,
                                         matrix.n.O.vs.n.E.logit,
                                         matrix.cox.int.slope.logit,
                                         matrix.optimal.logit,
                                         matrix.sensitivity.logit,
                                         matrix.specificity.logit,
                                         matrix.ppv.logit,
                                         matrix.npv.logit,
                                         matrix.plr.logit,
                                         matrix.nlr.logit)
    
    colnames(perf.measures.logit.df) <-c("Imputation_no",
                                         "auc",
                                         "mse",
                                         "R2.Nagelkerke",
                                         "brier.score",
                                         "concordance",
                                         "n.predicted/n.outcomes",
                                         "cox.int",
                                         "cox.slope",
                                         "threshold.prob",
                                         "sensitivity",
                                         "specificity",
                                         "ppv",
                                         "npv",
                                         "LR+",
                                         "LR-")
    
    
    
    logit.results.output <- merge(logit.results.df,
                                  perf.measures.logit.df,
                                  by="Imputation_no")
    
    logit.results.output <- logit.results.output%>%
                             dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                               dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
    
    return(list("logistic.regression"=dplyr::as_tibble(logit.results.output),
                "model.option"=model.options))
    
  }
  
  ################################################################################
  ####### Both Regularised  and logistic regression
  ################################################################################
  
  if(model.options==3)
  {
    
    pred.vars <- purrr::map(seq(data.set.no),
                            function(j)
                            {
                              data.matrix(dfs.boot[[j]] %>%
                                            dplyr::select(all_of(pred.use)))})
    
    outcome.var <- purrr::map(seq(data.set.no),
                              function(j)
                              {
                                as.matrix(dfs.boot[[j]] %>%
                                            dplyr::select(all_of(outcome.use)))})
    
    # For each data, find lambda that minimize the cross-validation prediction error rate. 
    
    
    cv.lasso.dev <- lapply(seq(data.set.no),
                           function(j)
                           {
                             glmnet::cv.glmnet(x=pred.vars[[j]],
                                               y=outcome.var[[j]], 
                                               alpha = alpha.param,
                                               nfolds=10,
                                               type.measure="deviance",
                                               family = "binomial")})
    
    best.lambda.dev <- lapply(seq(data.set.no),
                              function(k)
                              {
                                cv.lasso.dev[[k]]$lambda.min})
    
    
    # Fit the final models
    lasso.model <- lapply(seq(data.set.no),
                          function(s)
                          {
                            glmnet::glmnet(x=pred.vars[[s]],
                                           y=outcome.var[[s]], 
                                           alpha = alpha.param,
                                           family = "binomial",
                                           lambda = best.lambda.dev[[s]])})
    
    
    ##### COEFFICIENTS
    lasso.coefs <- lapply(seq(data.set.no),
                          function(s)
                          {
                            
                            predict(object = lasso.model[[s]],
                                    s = best.lambda.dev[[s]],
                                    type = "coefficients")
                          })
    
    
    
    
    lasso.coef.matrix <- lapply(seq(data.set.no),
                                function(s)
                                {
                                  as.matrix(lasso.coefs[[s]])
                                })
    
    matrix.coefs <- matrix(nrow=length(lasso.coef.matrix), 
                           ncol=length(lasso.coef.matrix[[1]]))
    for (j in seq(data.set.no))
    {
      
      matrix.coefs[j,] <- as.matrix(lasso.coefs[[j]])
      
    }
    
    
    lasso.coefficients <- data.frame(matrix.coefs)
    colnames(lasso.coefficients) <- c("Intercept",paste( "Beta",pred.use,sep="_"))
    lasso.coefficients$Imputation_no <- seq(data.set.no)
    
    lasso.coefficients.df <- lasso.coefficients[, c(ncol(lasso.coefficients),
                                                    seq(ncol(lasso.coefficients)-1))]
    
    
    ######### MEASURES OF PEFORMANCE
    gofstats.lasso <- lapply(seq(data.set.no),
                             function(s)
                             {
                               t(glmnet::assess.glmnet(object=lasso.model[[s]],
                                                       newx=pred.vars[[s]],
                                                       newy=outcome.var[[s]],
                                                       family="binomial",
                                                       s= best.lambda.dev[[s]]))
                             })
    
    
    matrix.gofs.lasso <- matrix(nrow=length(gofstats.lasso),
                                ncol=5)
    
    for (j in seq(data.set.no))
    {
      
      matrix.gofs.lasso[j,] <- as.matrix.data.frame(gofstats.lasso[[j]])
      
    }
    
    
    df.gof.lasso <- as.data.frame(matrix.gofs.lasso)
    colnames(df.gof.lasso) <- c("deviance","class","auc","mse","mae")
    
    # Nagalkerke's R2
    R2.Nagelkerke <- lapply(seq(data.set.no),
                            function(k)
                            {
                              lasso.model[[k]]$dev.ratio[which(lasso.model[[k]]$lambda==min(lasso.model[[k]]$lambda))]
                            })
    
    
    
    
    
    
    df.gof.lasso$R2.Nagelkerke <- as.vector(unlist(R2.Nagelkerke)) 
    df.gof.lasso$Imputation_no <- seq(data.set.no)
    
    
    
    #####
    
    df.gof.stats.lasso <- df.gof.lasso[, c(ncol(df.gof.lasso),
                                           seq(ncol(df.gof.lasso)-1))]
    
    
    df.lasso.results <- merge(x=lasso.coefficients.df,
                              y=df.gof.stats.lasso,
                              by="Imputation_no")
    
    df.lasso.results.output <- df.lasso.results %>%
      dplyr::select(-deviance,-class,-mae)
    
    
    # predictions
    preds.lasso <- lapply(seq(data.set.no),
                          function(k)
                          {
                            stats::predict(lasso.model[[k]], 
                                           newx = pred.vars[[k]], 
                                           type = 'response',
                                           s= best.lambda.dev[[k]])
                            
                            
                          })
    
    
    
    # optimal threshold probability for classification
    
    if(is.null(thres.prob.classifier))
    {
      optimal.lasso.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     InformationValue::optimalCutoff(dfs.boot[[k]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[k]])
                                     
                                   })
    }
    
    if(!is.null(thres.prob.classifier))
    {
      optimal.lasso.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     thres.prob.classifier
                                   })
    }
    
    
    # Brier score
    brier.score.lasso.lists <- purrr::map(seq(data.set.no),
                                          function(x)
                                          {
                                            brier.score(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                        pred=as.vector(preds.lasso[[x]]),
                                                        threshold=optimal.lasso.list[[x]])
                                            
                                          })
    
    
    
    # concordance  
    concordance.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::Concordance(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]])
                                       
                                     })
    
    # sensitivity
    sensitivity.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::sensitivity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]],
                                                                     optimal.lasso.list[[s]])
                                       
                                     })
    # outcome vs expected scaled difference
    outcome.vs.expected.ratio.lasso <- lapply(seq(data.set.no),
                                              function(x)
                                              {
                                                outcome_expected.scaled.diff(outcome=dfs.boot[[x]]%>%dplyr::select(all_of(outcome)),
                                                                             pred =preds.lasso[[x]],
                                                                             threshold.prob =optimal.lasso.list[[x]])
                                              })
    
    # cox intercept and slope
    cox.int.slope.lasso <- lapply(seq(data.set.no), 
                                  function(x)
                                  {
                                    
                                    calibration.int.slope(outcome=as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome))),
                                                          pred=preds.lasso[[x]])                                  
                                  })
    
    # specificity
    specificity.lasso.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::specificity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     preds.lasso[[s]],
                                                                     threshold=optimal.lasso.list[[s]])
                                       
                                     })
    
    
    
    # ppv 
    ppv.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::precision(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                           preds.lasso[[s]],
                                                           threshold=optimal.lasso.list[[s]])
                               
                             })
    
    # npv
    npv.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::npv(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                     preds.lasso[[s]],
                                                     threshold=optimal.lasso.list[[s]])
                               
                             })
    
    # plr
    
    plr.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               sensitivity.lasso.list[[s]]/(1-specificity.lasso.list[[s]])
                               
                             })
    
    # nlr
    
    nlr.lasso.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               (1-sensitivity.lasso.list[[s]])/specificity.lasso.list[[s]]
                               
                             })
    
    
    matrix.brier.score.lasso <- matrix(nrow=length(brier.score.lasso.lists))
    matrix.concordance.lasso <- matrix(nrow=length(concordance.lasso.list))
    matrix.n.O.vs.n.E.lasso <- matrix(nrow=length(outcome.vs.expected.ratio.lasso))
    matrix.cox.int.slope.lasso <- matrix(nrow=length(cox.int.slope.lasso),
                                         ncol=2)
    matrix.optimal.lasso <- matrix(nrow=length(optimal.lasso.list),
                                   ncol=length(optimal.lasso.list[[1]]))
    matrix.sensitivity.lasso <- matrix(nrow=length(sensitivity.lasso.list))
    matrix.specificity.lasso <- matrix(nrow=length(sensitivity.lasso.list))
    matrix.ppv.lasso <- matrix(nrow=length(ppv.lasso.list))
    matrix.npv.lasso <- matrix(nrow=length(npv.lasso.list))
    matrix.plr.lasso <- matrix(nrow=length(plr.lasso.list))
    matrix.nlr.lasso <- matrix(nrow=length(nlr.lasso.list))
    
    
    for (k in seq(data.set.no))
    {
      matrix.brier.score.lasso[k,] <- brier.score.lasso.lists[[k]]
      matrix.concordance.lasso[k,] <- concordance.lasso.list[[k]]$Concordance
      matrix.n.O.vs.n.E.lasso[k,] <- outcome.vs.expected.ratio.lasso[[k]]
      matrix.cox.int.slope.lasso[k,] <- cox.int.slope.lasso[[k]]
      matrix.optimal.lasso[k,] <- optimal.lasso.list[[k]]
      matrix.sensitivity.lasso[k,] <- sensitivity.lasso.list[[k]]
      matrix.specificity.lasso[k,] <- specificity.lasso.list[[k]]
      matrix.ppv.lasso[k,]<- ppv.lasso.list[[k]]
      matrix.npv.lasso[k,]<- npv.lasso.list[[k]]
      matrix.plr.lasso[k,]<- plr.lasso.list[[k]]
      matrix.nlr.lasso[k,]<- nlr.lasso.list[[k]]
    }
    
    perf.measures.lasso.df <- data.frame(seq(data.set.no),
                                         matrix.brier.score.lasso,
                                         matrix.concordance.lasso,
                                         matrix.n.O.vs.n.E.lasso,
                                         matrix.cox.int.slope.lasso,
                                         matrix.optimal.lasso,
                                         matrix.sensitivity.lasso,
                                         matrix.specificity.lasso,
                                         matrix.ppv.lasso,
                                         matrix.npv.lasso,
                                         matrix.plr.lasso,
                                         matrix.nlr.lasso)
    
    colnames(perf.measures.lasso.df) <-c("Imputation_no",
                                         "brier.score",
                                         "concordance",
                                         "n.predicted/n.outcomes",
                                         "cox.int",
                                         "cox.slope",
                                         "threshold.prob",
                                         "sensitivity",
                                         "specificity",
                                         "ppv",
                                         "npv",
                                         "LR+",
                                         "LR-")
    
    
    
    lasso.results.output <- merge(df.lasso.results.output,
                                  perf.measures.lasso.df,
                                  by="Imputation_no")
    
    
    lasso.results.output <- lasso.results.output%>%
                              dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                                dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
    
    
    ################################################################################
    #  LOGISTIC MODEL WITH LOGIT LINK
    ################################################################################
    
    formula.model <- as.formula(paste(outcome.use, 
                                      paste(pred.use, 
                                            collapse=" + "),
                                      sep=" ~ "))
    
    logit.results.list <- lapply(seq(data.set.no),
                                 function(p)
                                 {
                                   glm(formula.model,
                                       family=binomial(link=logit),
                                       data=dfs.boot[[p]])
                                 })
    
    
    
    logit.results.matrix <- matrix(nrow=length(logit.results.list),
                                   ncol=length(logit.results.list[[1]]$coefficients))
    
    for (i in seq(data.set.no))
    {
      logit.results.matrix[i,] <- logit.results.list[[i]]$coefficients
    }
    
    logit.results.df <- as.data.frame(logit.results.matrix)
    colnames(logit.results.df) <- c("Intercept",paste( "Beta",pred.use,sep="_"))
    
    logit.results.df$Imputation_no <- seq(data.set.no)
    
    logit.results.df <- logit.results.df[,c(ncol(logit.results.df),
                                            seq(ncol(logit.results.df)-1))]
    
    
    # Performance measures
    
    if(!(require(pROC))){install.packages("pROC", dependencies=TRUE)}
    if(!(require(fmsb))){install.packages("fmsb", dependencies=TRUE)}
    # auc
    auc.logit <- lapply(seq(data.set.no),
                        function(s)
                        {
                          suppressMessages(pROC::auc(response=as.vector(as.matrix(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)))),
                                                     predictor=logit.results.list[[s]]$fitted.values))
                          #suppressMessages() prevents progress status from being printed
                        })
    
    
    
    
    # mse
    mse.logit <- lapply(seq(data.set.no),
                        function(s)
                        {
                          (sum(logit.results.list[[s]]$residuals)/
                             (nrow(dfs.boot[[s]])-length(pred.use)))
                        })
    
    # R2 Nagelkerke
    R2.Nagelkerke.logit.lists <- lapply(seq(data.set.no),
                                        function(s)
                                        {
                                          fmsb::NagelkerkeR2(logit.results.list[[s]])
                                        })
    
    
    # concordance  
    concordance.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::Concordance(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values)
                                       
                                     })
    
    # cox intercept and slope
    cox.int.slope.logit <- lapply(seq(data.set.no), 
                                  function(x)
                                  {
                                    
                                    calibration.int.slope(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                          pred=logit.results.list[[x]]$fitted.values)                                  
                                  })
    # threshold probability
    # optimal threshold probability for classification
    if(is.null(thres.prob.classifier))
    {
      optimal.logit.list <- lapply(seq(data.set.no),
                                   function(s)
                                   {
                                     InformationValue::optimalCutoff(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values)
                                   })
      
    }
    
    if(!is.null(thres.prob.classifier))
    {
      optimal.logit.list <- lapply(seq(data.set.no),
                                   function(k)
                                   {
                                     thres.prob.classifier
                                   })
    }
    
    # Brier score
    
    brier.score.logit.lists <- purrr::map(seq(data.set.no),
                                          function(x)
                                          {
                                            brier.score(outcome=as.vector(as.matrix(dfs.boot[[x]]%>%dplyr::select(all_of(outcome)))),
                                                        pred=as.vector(logit.results.list[[x]]$fitted.values),
                                                        threshold=optimal.logit.list[[x]])
                                            
                                          })
    
    # sensitivity
    sensitivity.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::sensitivity(as.vector(as.matrix(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)))),
                                                                     as.vector(logit.results.list[[s]]$fitted.values),
                                                                     threshold=optimal.logit.list[[s]])
                                       
                                     })
    
    # specificity
    specificity.logit.list <- lapply(seq(data.set.no),
                                     function(s)
                                     {
                                       InformationValue::specificity(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                                     logit.results.list[[s]]$fitted.values,
                                                                     threshold=optimal.logit.list[[s]])
                                       
                                     })
    
    
    # outcome vs expected scaled difference
    outcome.vs.expected.ratio.logit <- lapply(seq(data.set.no),
                                              function(x)
                                              {
                                                outcome_expected.scaled.diff(outcome=dfs.boot[[x]]%>%dplyr::select(all_of(outcome)),
                                                                             pred =logit.results.list[[x]]$fitted.values,
                                                                             threshold.prob =optimal.logit.list[[x]])
                                              })
    
    
    # ppv 
    ppv.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::precision(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                           logit.results.list[[s]]$fitted.values,
                                                           threshold=optimal.logit.list[[s]])
                               
                             })
    
    # npv
    npv.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               InformationValue::npv(dfs.boot[[s]]%>%dplyr::select(all_of(outcome)),
                                                     logit.results.list[[s]]$fitted.values,
                                                     threshold=optimal.logit.list[[s]])
                               
                             })
    
    # plr
    
    plr.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               sensitivity.logit.list[[s]]/(1-specificity.logit.list[[s]])
                               
                             })
    
    # nlr
    
    nlr.logit.list <- lapply(seq(data.set.no),
                             function(s)
                             {
                               (1-sensitivity.logit.list[[s]])/specificity.logit.list[[s]]
                               
                             })
    
    
    
    matrix.auc.logit <- matrix(nrow=length(auc.logit))
    matrix.mse.logit <- matrix(nrow=length(mse.logit))
    matrix.R2.Nagelkerke.logit <- matrix(nrow=length(R2.Nagelkerke.logit.lists))
    matrix.brier.score.logit <- matrix(nrow=length(brier.score.logit.lists))
    matrix.concordance.logit <- matrix(nrow=length(concordance.logit.list))
    matrix.n.O.vs.n.E.logit <- matrix(nrow=length(outcome.vs.expected.ratio.logit))
    matrix.cox.int.slope.logit <- matrix(nrow=length(cox.int.slope.logit),
                                         ncol=2)
    matrix.optimal.logit <- matrix(nrow=length(optimal.logit.list))
    matrix.sensitivity.logit <- matrix(nrow=length(sensitivity.logit.list))
    matrix.specificity.logit <- matrix(nrow=length(sensitivity.logit.list))
    matrix.ppv.logit <- matrix(nrow=length(ppv.logit.list))
    matrix.npv.logit <- matrix(nrow=length(npv.logit.list))
    matrix.plr.logit <- matrix(nrow=length(plr.logit.list))
    matrix.nlr.logit <- matrix(nrow=length(nlr.logit.list))
    
    
    for (k in seq(data.set.no))
    {
      matrix.auc.logit[k,] <- as.numeric(gsub("Area under the curve:",
                                              "", auc.logit[[k]]))
      matrix.mse.logit[k,] <- mse.logit[[k]]
      matrix.R2.Nagelkerke.logit[k,] <- R2.Nagelkerke.logit.lists[[k]]$R2
      matrix.brier.score.logit[k,] <- brier.score.logit.lists[[k]]
      matrix.concordance.logit[k,] <- concordance.logit.list[[k]]$Concordance
      matrix.n.O.vs.n.E.logit[k,] <- outcome.vs.expected.ratio.logit[[k]]
      matrix.cox.int.slope.logit[k,] <- cox.int.slope.logit[[k]]
      matrix.optimal.logit[k,] <- optimal.logit.list[[k]]
      matrix.sensitivity.logit[k,] <- sensitivity.logit.list[[k]]
      matrix.specificity.logit[k,] <- specificity.logit.list[[k]]
      matrix.ppv.logit[k,] <- ppv.logit.list[[k]]
      matrix.npv.logit[k,] <- npv.logit.list[[k]]
      matrix.plr.logit[k,]<- plr.logit.list[[k]]
      matrix.nlr.logit[k,]<- nlr.logit.list[[k]]
      
    }
    
    perf.measures.logit.df <- data.frame(seq(data.set.no),
                                         matrix.auc.logit,
                                         matrix.mse.logit,
                                         matrix.R2.Nagelkerke.logit,
                                         matrix.brier.score.logit,
                                         matrix.concordance.logit,
                                         matrix.n.O.vs.n.E.logit,
                                         matrix.cox.int.slope.logit,
                                         matrix.optimal.logit,
                                         matrix.sensitivity.logit,
                                         matrix.specificity.logit,
                                         matrix.ppv.logit,
                                         matrix.npv.logit,
                                         matrix.plr.logit,
                                         matrix.nlr.logit)
    
    colnames(perf.measures.logit.df) <-c("Imputation_no",
                                         "auc",
                                         "mse",
                                         "R2.Nagelkerke",
                                         "brier.score",
                                         "concordance",
                                         "n.predicted/n.outcomes",
                                         "cox.int",
                                         "cox.slope",
                                         "threshold.prob",
                                         "sensitivity",
                                         "specificity",
                                         "ppv",
                                         "npv",
                                         "LR+",
                                         "LR-")
    
    
    
    logit.results.output <- merge(logit.results.df,
                                  perf.measures.logit.df,
                                  by="Imputation_no")
    
    logit.results.output <- logit.results.output%>%
                               dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                                   dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
    
    return(list("regularised.regression"=dplyr::as_tibble(lasso.results.output),
                "logistic.regression"=dplyr::as_tibble(logit.results.output),
                "model.option"=model.options))
  }
}


## This allows multiple probability threshold to be specified by the user

glm.net.no_boot.mi.fit <- function(data.boot.mi_complete,
                                   predictors,
                                   outcome,
                                   data.set.no,
                                   thres.prob.classifier=NULL,
                                   model.options,
                                   alpha.param,
                                   seed.input)
{
  #####
  if (is.null(thres.prob.classifier))
  {
    output <- regularisation.no_boot.mi.logit(data.boot.mi_complete=data.boot.mi_complete,
                                              predictors=predictors,
                                              outcome=outcome,
                                              data.set.no=data.set.no,
                                              thres.prob.classifier=NULL,
                                              model.options=model.options,
                                              alpha.param=alpha.param,
                                              seed.input=seed.input) 
  }
  #####
  if (!is.null(thres.prob.classifier))
  {
    
    
    output <-  purrr::map(thres.prob.classifier, 
                          function(x)
                          {
                            regularisation.no_boot.mi.logit(data.boot.mi_complete=data.boot.mi_complete,
                                                            predictors=predictors,
                                                            outcome=outcome,
                                                            data.set.no=data.set.no,
                                                            thres.prob.classifier=x,
                                                            model.options=model.options,
                                                            alpha.param=alpha.param,
                                                            seed.input=seed.input) 
                          })
    
  }
  ######
  
  return(output)
}


################################################################################
################################################################################
#  FUNCTION FOR AGGREGATING SOLUTION FROM regularisation.fit.logit.link()
################################################################################
################################################################################

parameter_estimates.mi_no_boot <-  function(results,
                                            alpha=NULL)
{
  ##
  if(is.null(alpha)){alpha <- 0.05}
  
  if(alpha < 0 | alpha >1)
  { 
    stop("The type I error  MUST be between 0 and 1")
  }
  
  
  
  LL.null <- "BCa.LL_"
  UL.null <- "BCa.UL_"
  
  ##############################################################################
  ####### When probability threshold not specified
  ##############################################################################
  
  if (purrr::vec_depth(results)==3)
  {
    
    ############# model.option=1
    if(results$model.option==1)
    {
      ############# Regularised model
      reg.df <- results$regularised.regression %>% 
                      dplyr::select(-Imputation_no)
      ###
      results.df.regularised <- reg.df %>%
        dplyr::summarise_if(is.numeric,
                            mean,
                            na.rm = TRUE)
      ###
      results.df.regularised.bca <-  map(seq(ncol(reg.df)),
                                         function(x)
                                         {
                                           coxed::bca(as.matrix(na.omit(reg.df[,x])), 
                                                      conf.level =1-alpha)
                                         })
      
      ###
      soln.reg.df <- round(cbind(t(results.df.regularised),
                                 data.frame(matrix(unlist(results.df.regularised.bca), 
                                                   nrow=length(results.df.regularised.bca), 
                                                   byrow=TRUE),
                                            stringsAsFactors=FALSE)),4)
      
      colnames(soln.reg.df) <- c(noquote(paste0("Mean")),
                                 noquote(paste0(LL.null,
                                                (alpha/2)*100,
                                                "%",
                                                collapse="")),
                                 noquote(paste0(UL.null,
                                                (1-alpha/2)*100,
                                                "%",
                                                collapse="")))
      
      ####
      output <- list("regularised.regression"=soln.reg.df)
      
      return(output)
    }
    
    ############# model.option=2
    if(results$model.option==2)
    {
      ##################### Logistic model
      log.df <- results$logistic.regression %>% 
        dplyr::select(-Imputation_no)
      ###
      results.df.logistic <- log.df %>%
        dplyr::summarise_if(is.numeric,
                            mean,
                            na.rm = TRUE)
      ###
      results.df.logistic.bca <-  map(seq(ncol(log.df)),
                                      function(x)
                                      {
                                        coxed::bca(as.matrix(na.omit(log.df[,x])), 
                                                   conf.level =1-alpha)
                                      })
      
      ###
      soln.log.df <- round(cbind(t(results.df.logistic),
                                 data.frame(matrix(unlist(results.df.logistic.bca), 
                                                   nrow=length(results.df.logistic.bca), 
                                                   byrow=TRUE),
                                            stringsAsFactors=FALSE)),4)
      
      colnames(soln.log.df) <- c(noquote(paste0("Mean")),
                                 noquote(paste0(LL.null,
                                                (alpha/2)*100,
                                                "%",
                                                collapse="")),
                                 noquote(paste0(UL.null,
                                                (1-alpha/2)*100,
                                                "%",
                                                collapse="")))
      
      output <- list("logistic.regression"=soln.log.df)
      
      return(output)
    }  
    ############# model.option=3
    if(results$model.option==3)
    {
      ############# Regularised model
      reg.df <- results$regularised.regression %>% 
        dplyr::select(-Imputation_no)
      ###
      results.df.regularised <- reg.df %>%
        dplyr::summarise_if(is.numeric,
                            mean,
                            na.rm = TRUE)
      ###
      results.df.regularised.bca <-  map(seq(ncol(reg.df)),
                                         function(x)
                                         {
                                           coxed::bca(as.matrix(na.omit(reg.df[,x])), 
                                                      conf.level =1-alpha)
                                         })
      
      ###
      soln.reg.df <- round(cbind(t(results.df.regularised),
                                 data.frame(matrix(unlist(results.df.regularised.bca), 
                                                   nrow=length(results.df.regularised.bca), 
                                                   byrow=TRUE),
                                            stringsAsFactors=FALSE)),4)
      
      colnames(soln.reg.df) <- c(noquote(paste0("Mean")),
                                 noquote(paste0(LL.null,
                                                (alpha/2)*100,
                                                "%",
                                                collapse="")),
                                 noquote(paste0(UL.null,
                                                (1-alpha/2)*100,
                                                "%",
                                                collapse="")))
      
      ##################### Logistic model
      log.df <- results$logistic.regression %>% 
        dplyr::select(-Imputation_no)
      ###
      results.df.logistic <- log.df %>%
                             dplyr::summarise_if(is.numeric,
                                                mean,
                                                na.rm = TRUE)
      ###
      results.df.logistic.bca <-  map(seq(ncol(log.df)),
                                      function(x)
                                      {
                                        coxed::bca(as.matrix(na.omit(log.df[,x])), 
                                                   conf.level =1-alpha)
                                      })
      
      ###
      soln.log.df <- round(cbind(t(results.df.logistic),
                                 data.frame(matrix(unlist(results.df.logistic.bca), 
                                                   nrow=length(results.df.logistic.bca), 
                                                   byrow=TRUE),
                                            stringsAsFactors=FALSE)),4)
      
      colnames(soln.log.df) <- c(noquote(paste0("Mean")),
                                 noquote(paste0(LL.null,
                                                (alpha/2)*100,
                                                "%",
                                                collapse="")),
                                 noquote(paste0(UL.null,
                                                (1-alpha/2)*100,
                                                "%",
                                                collapse="")))
      
      
      output <- list("regularised.regression"=soln.reg.df,
                     "logistic.regression"=soln.log.df)
      
      return(output)
    }
    
    
  }
  
  ##############################################################################
  ####### When probability threshold(s) is(are) specified
  ##############################################################################
  
  
  if (purrr::vec_depth(results)==4)
  {
    #### R functions for Bca CI
    coxbca.l <- function(x) 
    {
      output <- coxed::bca(as.matrix(na.omit(x)),
                           conf.level =1-alpha)
      
      output.l <- output[1]
      
      return(output.l)
    }
    
    coxbca.u <- function(x) 
    {
      output <- coxed::bca(as.matrix(na.omit(x)),
                           conf.level =1-alpha)
      
      output.u <- output[2]
      
      return(output.u)
    }
    
    ######### model.option=1
    if(results[[1]]$model.option==1)
    {
      ##### Regularised model
      ###
      results.df.regularised <- t(purrr::map_df(.x=seq(length(results)),
                                                .f=function(x)
                                                {
                                                  results[[x]]$regularised.regression %>%
                                                    dplyr::select(-Imputation_no)%>%
                                                    dplyr::group_by(threshold.prob)%>%
                                                    dplyr::summarise_if(is.numeric,
                                                                        mean, 
                                                                        na.rm = TRUE)
                                                  
                                                }))
      
      ###
      reg.bca.l <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$regularised.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$regularised.regression)[-1],
                                                 coxbca.l)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      ###
      reg.bca.u <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$regularised.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$regularised.regression)[-1],
                                                 coxbca.u)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      
      ###
      # # replace NaN with NA_real_
      # 
      results.df.regularised <- as.data.frame(results.df.regularised)%>%
                                  dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                                   dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      #
      reg.bca.l <- as.data.frame(reg.bca.l)%>%
                     dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                         dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      #
      reg.bca.u <- as.data.frame(reg.bca.u)%>%
                       dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                           dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))

      ##
      
      soln.reg.df <- list()
      unique.prob.threshold <- results.df.regularised[1,]
      rownames.use <- list()
      rownames.vec <- list()
      soln.reg.df.output <- list()
      
      for (i in seq(length(results)))
      {
        
        library(qpcR)
        soln.reg.df[[i]] <- round(qpcR:::cbind.na(results.df.regularised[-1,i],
                                             reg.bca.l[!(row.names(reg.bca.l) %in% c("threshold.prob")),i],
                                             reg.bca.u[!(row.names(reg.bca.l) %in% c("threshold.prob")),i]),4)
        
        rownames.use[[i]] <- colnames(results[[i]]$regularised.regression)
        
        rownames.vec[[i]] <- rownames.use[[i]][!rownames.use[[i]]%in%(c("Imputation_no","threshold.prob"))]
        
        
        
        colnames(soln.reg.df[[i]]) <- c(noquote(paste0("PT_", 
                                                       unique.prob.threshold[i],
                                                       "_Mean",
                                                       collapse="")),
                                        noquote(paste0(LL.null,
                                                       (alpha/2)*100,
                                                       "%",
                                                       collapse="")),
                                        noquote(paste0(UL.null,
                                                       (1-alpha/2)*100,
                                                       "%",
                                                       collapse="")))
        
        soln.reg.df.output[[i]] <-  cbind.data.frame(rownames.vec[[i]],soln.reg.df[[i]])
        
        colnames(soln.reg.df.output[[i]]) <- c("Parameter",colnames(soln.reg.df[[i]]))
        
      }
      
      output <- list("regularised.regression"=soln.reg.df.output)
      
      return(output)
      
    }
    
    ######### model.option=2
    if(results[[1]]$model.option==2)
    {
      
      ##### Logistic model
      ###
      
      results.df.logistic <- t(purrr::map_df(.x=seq(length(results)),
                                             .f=function(x)
                                             {
                                               results[[x]]$logistic.regression %>%
                                                 dplyr::select(-Imputation_no)%>%
                                                 dplyr::group_by(threshold.prob)%>%
                                                 dplyr::summarise_if(is.numeric,
                                                                     mean, 
                                                                     na.rm = TRUE)
                                               
                                             }))
      
      ###
      log.bca.l <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$logistic.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$logistic.regression)[-1],
                                                 coxbca.l)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      ###
      log.bca.u <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$logistic.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$logistic.regression)[-1],
                                                 coxbca.u)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      
      ###
      # # replace NaN with NA_real_
      # 
      results.df.logistic <- as.data.frame(results.df.logistic)%>%
                                dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                                dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      #
      log.bca.l <- as.data.frame(log.bca.l)%>%
                      dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                        dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      #
      log.bca.u <- as.data.frame(log.bca.u)%>%
                       dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                           dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      ##
      
      
      soln.log.df <- list()
      unique.prob.threshold <- results.df.logistic[1,]
      rownames.use <- list()
      rownames.vec <- list()
      soln.log.df.output <- list()
      for (i in seq(length(results)))
      {
        
        library(qpcR)
       
        
        soln.log.df[[i]] <- round(qpcR:::cbind.na(results.df.logistic[-1,i],
                                                  log.bca.l[!(row.names(log.bca.l) %in% c("threshold.prob")),i],
                                                  log.bca.u[!(row.names(log.bca.l) %in% c("threshold.prob")),i]),4)
        
        
        rownames.use[[i]] <- colnames(results[[i]]$logistic.regression)
        
        rownames.vec[[i]] <- rownames.use[[i]][!rownames.use[[i]]%in%(c("Imputation_no","threshold.prob"))]
        
        
        
        colnames(soln.log.df[[i]]) <- c(noquote(paste0("PT_", 
                                                       unique.prob.threshold[i],
                                                       "_Mean",
                                                       collapse="")),
                                        noquote(paste0(LL.null,
                                                       (alpha/2)*100,
                                                       "%",
                                                       collapse="")),
                                        noquote(paste0(UL.null,
                                                       (1-alpha/2)*100,
                                                       "%",
                                                       collapse="")))
        
        soln.log.df.output[[i]] <-  cbind.data.frame(rownames.vec[[i]],soln.log.df[[i]])
        
        colnames(soln.log.df.output[[i]]) <- c("Parameter",colnames(soln.log.df[[i]]))
        
        
        
      }
      
      output <- list("logistic.regression"=soln.log.df.output)
      
      return(output)
      
      
    } 
    
    ######### model.option=3
    if(results[[1]]$model.option==3)
    {
      ##### Regularised model
      ###
      results.df.regularised <- t(purrr::map_df(.x=seq(length(results)),
                                                .f=function(x)
                                                {
                                                  results[[x]]$regularised.regression %>%
                                                    dplyr::select(-Imputation_no)%>%
                                                    dplyr::group_by(threshold.prob)%>%
                                                    dplyr::summarise_if(is.numeric,
                                                                        mean, 
                                                                        na.rm = TRUE)
                                                  
                                                }))
      
      ###
      reg.bca.l <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$regularised.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$regularised.regression)[-1],
                                                 coxbca.l)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      ###
      reg.bca.u <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$regularised.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$regularised.regression)[-1],
                                                 coxbca.u)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      
      ###
      # # replace NaN with NA_real_
      # 
      results.df.regularised <- as.data.frame(results.df.regularised)%>%
        dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
        dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      #
      reg.bca.l <- as.data.frame(reg.bca.l)%>%
        dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
        dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      #
      reg.bca.u <- as.data.frame(reg.bca.u)%>%
        dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
        dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      ##
      
      soln.reg.df <- list()
      unique.prob.threshold <- results.df.regularised[1,]
      rownames.use <- list()
      rownames.vec <- list()
      soln.reg.df.output <- list()
      
      for (i in seq(length(results)))
      {
        
        library(qpcR)
        soln.reg.df[[i]] <- round(qpcR:::cbind.na(results.df.regularised[-1,i],
                                                  reg.bca.l[!(row.names(reg.bca.l) %in% c("threshold.prob")),i],
                                                  reg.bca.u[!(row.names(reg.bca.l) %in% c("threshold.prob")),i]),4)
        
        rownames.use[[i]] <- colnames(results[[i]]$regularised.regression)
        
        rownames.vec[[i]] <- rownames.use[[i]][!rownames.use[[i]]%in%(c("Imputation_no","threshold.prob"))]
        
        
        
        colnames(soln.reg.df[[i]]) <- c(noquote(paste0("PT_", 
                                                       unique.prob.threshold[i],
                                                       "_Mean",
                                                       collapse="")),
                                        noquote(paste0(LL.null,
                                                       (alpha/2)*100,
                                                       "%",
                                                       collapse="")),
                                        noquote(paste0(UL.null,
                                                       (1-alpha/2)*100,
                                                       "%",
                                                       collapse="")))
        
        soln.reg.df.output[[i]] <-  cbind.data.frame(rownames.vec[[i]],soln.reg.df[[i]])
        
        colnames(soln.reg.df.output[[i]]) <- c("Parameter",colnames(soln.reg.df[[i]]))
        
      }
      
      ##### Logistic model
      ###
     
      results.df.logistic <- t(purrr::map_df(.x=seq(length(results)),
                                             .f=function(x)
                                             {
                                               results[[x]]$logistic.regression %>%
                                                 dplyr::select(-Imputation_no)%>%
                                                 dplyr::group_by(threshold.prob)%>%
                                                 dplyr::summarise_if(is.numeric,
                                                                     mean, 
                                                                     na.rm = TRUE)
                                               
                                             }))
      
      ###
      log.bca.l <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$logistic.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$logistic.regression)[-1],
                                                 coxbca.l)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      ###
      log.bca.u <- t(map_df(.x=seq(length(results)),
                            .f=function(x)
                            {
                              results[[x]]$logistic.regression %>%
                                dplyr::select(-Imputation_no)%>%
                                dplyr::group_by(threshold.prob)%>%
                                dplyr::mutate_at(colnames(results[[x]]$logistic.regression)[-1],
                                                 coxbca.u)%>%
                                dplyr::distinct(.keep_all=TRUE)
                            }))
      
      ###
      # # replace NaN with NA_real_
      # 
      results.df.logistic <- as.data.frame(results.df.logistic)%>%
                              dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                               dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      #
      log.bca.l <- as.data.frame(log.bca.l)%>%
                      dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                       dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      #
      log.bca.u <- as.data.frame(log.bca.u)%>%
                    dplyr::mutate_all(~ifelse(is.nan(.), NA, .))%>%
                       dplyr::mutate_if(is.numeric, list(~na_if(., Inf)))
      
      ##
      
      
      soln.log.df <- list()
      unique.prob.threshold <- results.df.logistic[1,]
      rownames.use <- list()
      rownames.vec <- list()
      soln.log.df.output <- list()
      for (i in seq(length(results)))
      {
        
        library(qpcR)
        
        
        soln.log.df[[i]] <- round(qpcR:::cbind.na(results.df.logistic[-1,i],
                                                  log.bca.l[!(row.names(log.bca.l) %in% c("threshold.prob")),i],
                                                  log.bca.u[!(row.names(log.bca.l) %in% c("threshold.prob")),i]),4)
        
        
        rownames.use[[i]] <- colnames(results[[i]]$logistic.regression)
        
        rownames.vec[[i]] <- rownames.use[[i]][!rownames.use[[i]]%in%(c("Imputation_no","threshold.prob"))]
        
        
        
        colnames(soln.log.df[[i]]) <- c(noquote(paste0("PT_", 
                                                       unique.prob.threshold[i],
                                                       "_Mean",
                                                       collapse="")),
                                        noquote(paste0(LL.null,
                                                       (alpha/2)*100,
                                                       "%",
                                                       collapse="")),
                                        noquote(paste0(UL.null,
                                                       (1-alpha/2)*100,
                                                       "%",
                                                       collapse="")))
        
        soln.log.df.output[[i]] <-  cbind.data.frame(rownames.vec[[i]],soln.log.df[[i]])
        
        colnames(soln.log.df.output[[i]]) <- c("Parameter",colnames(soln.log.df[[i]]))
        
      }
      
      
      ####
      output <- list("regularised.regression"=soln.reg.df.output,
                     "logistic.regression"=soln.log.df.output)
      
      return(output)
      
    }
  }
  
  
  
}