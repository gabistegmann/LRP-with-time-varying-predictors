lrp_tvp = function (method, nlme.model = NULL, randomFormula, fixedFormula = NULL, 
                         data, start, group, rPartFormula, weight = NULL, R = NULL, 
                         min.dev = NULL, control = rpart.control(), use_parallel = FALSE,
                         time_var_covs = NULL) 
{
  if (method == "lme") {
    lmeFormula <- fixedFormula
  }
  if (is.null(min.dev) == FALSE) {
    if (method == "lme") {
      mod <- lme(lmeFormula, data = data, random = randomFormula, 
                 correlation = R, na.action = na.omit)
    }
    else {
      mod <- nlme(model = nlme.model, fixed = fixedFormula, 
                  data = data, random = randomFormula, correlation = R, 
                  na.action = na.omit, start = start, groups = group)
    }
    control$cp <- 1 - (-2 * mod$logLik - min.dev)/(-2 * mod$logLik)
  }
  if (method == "lme") {
    groupingName = attr(terms(splitFormula(randomFormula, 
                                           "|")[[2]]), "term.labels")
    responseName = attr(terms(getResponseFormula(lmeFormula)), 
                        "term.labels")
    groupingFactor = data[, names(data) == groupingName]
    terms = attr(terms(lmeFormula), "term.labels")
    continuous = !is.factor(data[, names(data) == terms[1]])
    evaluation <- function(y, wt, parms) {
      model = lme(lmeFormula, data = parms[groupingFactor %in% 
                                             y, ], random = randomFormula, correlation = R, 
                  na.action = na.omit)
      if (continuous) {
        slope = model$coefficients$fixed[2]
      }
      else {
        levels = length(levels(data[, names(data) == 
                                      terms[1]]))
        y = model$coefficients$fixed[1:levels]
        x = 1:levels
        slope = lm(y ~ x)$coefficients[2]
      }
      list(label = slope, deviance = -2 * (model$logLik))
    }
  }
  else if (method == "nlme") {
    groupingName = attr(terms(splitFormula(group, "~")[[1]]), 
                        "term.labels")
    responseName = attr(terms(getResponseFormula(nlme.model)), 
                        "term.labels")
    groupingFactor = data[, names(data) == groupingName]
    evaluation <- function(y, wt, parms) {
      model = nlme(model = nlme.model, fixed = fixedFormula, 
                   data = parms[groupingFactor %in% y, ], random = randomFormula, 
                   correlation = R, na.action = na.omit, start = start, 
                   groups = group)
      slope = 1
      list(label = slope, deviance = -2 * (model$logLik))
    }
  }
  
  
  
  if(use_parallel == FALSE){
    
    split <- function(y, wt, x, parms, continuous) {
      print(paste("splitting:", length(unique(x)), "values"))
      dev = vector()
      xUnique = unique(x)
      if (method == "lme") {
        rootDev = lme(lmeFormula, data = parms[groupingFactor %in% 
                                                 y, ], random = randomFormula, correlation = R, 
                      na.action = na.omit)$logLik
      }
      else if (method == "nlme") {
        rootDev = nlme(model = nlme.model, fixed = fixedFormula, 
                       data = parms[groupingFactor %in% y, ], random = randomFormula, 
                       correlation = R, na.action = na.omit, start = start, 
                       groups = group)$logLik
      }
      if (continuous) {
        for (i in xUnique) {
          yLeft = y[x <= i]
          yRight = y[x > i]
          if (length(yLeft) < control$minbucket || length(yRight) < 
              control$minbucket) {
            dev = c(dev, 0)
          }
          else {
            if (method == "lme") {
              modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                             yLeft, ], random = randomFormula, correlation = R, 
                                  na.action = na.omit), silent = TRUE)
              modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                              yRight, ], random = randomFormula, correlation = R, 
                                   na.action = na.omit), silent = TRUE)
            }
            else if (method == "nlme") {
              modelLeft = try(nlme(model = nlme.model, 
                                   fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                        yLeft, ], random = randomFormula, correlation = R, 
                                   na.action = na.omit, start = start, groups = group), 
                              silent = TRUE)
              modelRight = try(nlme(model = nlme.model, 
                                    fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                         yRight, ], random = randomFormula, correlation = R, 
                                    na.action = na.omit, start = start, groups = group), 
                               silent = TRUE)
            }
            if (any(class(modelLeft) == "lme") && any(class(modelRight) == 
                                                      "lme")) {
              dev = c(dev, modelLeft$logLik + modelRight$logLik)
            }
            else {
              dev = c(dev, 0)
            }
          }
        }
        good = rep(0, length(x))
        for (i in 1:length(xUnique)) {
          good[x == xUnique[i]] = dev[i]
        }
        good = good[1:(length(good) - 1)]
        list(goodness = good + abs(rootDev) * (good != 0) * 
               2, direction = rep(-1, length(good)))
      }
      else {
        order = rep(0, length(xUnique))
        response = parms[, names(parms) == responseName]
        for (i in 1:length(xUnique)) {
          order[i] = mean(response[x == xUnique[i]], na.rm = TRUE)
        }
        dir = sort(order, index.return = TRUE)$ix
        for (i in 1:(length(dir) - 1)) {
          yLeft = y[x %in% dir[1:i]]
          yRight = y[x %in% dir[(i + 1):length(dir)]]
          if (length(yLeft) < control$minbucket || length(yRight) < 
              control$minbucket) {
            dev = c(dev, 0)
          }
          else {
            if (method == "lme") {
              modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                             yLeft, ], random = randomFormula, correlation = R, 
                                  na.action = na.omit), silent = TRUE)
              modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                              yRight, ], random = randomFormula, correlation = R, 
                                   na.action = na.omit), silent = TRUE)
            }
            else if (method == "nlme") {
              modelLeft = try(nlme(model = nlme.model, 
                                   fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                        yLeft, ], random = randomFormula, correlation = R, 
                                   na.action = na.omit, start = start, groups = group), 
                              silent = TRUE)
              modelRight = try(nlme(model = nlme.model, 
                                    fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                         yRight, ], random = randomFormula, correlation = R, 
                                    na.action = na.omit, start = start, groups = group), 
                               silent = TRUE)
            }
            if ((any(class(modelLeft) == "lme") | any(class(modelLeft) == 
                                                      "nlme")) && (any(class(modelRight) == "lme") | 
                                                                   any(class(modelRight) == "nlme"))) {
              dev = c(dev, modelLeft$logLik + modelRight$logLik)
            }
            else {
              dev = c(dev, 0)
            }
          }
        }
        list(goodness = dev + abs(rootDev) * (dev != 0) * 
               2, direction = dir)
      }
    }
    
  }else if(use_parallel == TRUE){
    
    ######IF use_parallel == TRUE:
    
    library(parallel)
    # Calculate the number of cores
    no_cores <- detectCores() - 1
    # Initiate cluster
    CL <- makeCluster(no_cores)
    
    
    split <- function(y, wt, x, parms, continuous) {
      
      print(paste("splitting:", length(unique(x)), "values"))
      dev = vector()
      xUnique = unique(x)
      if (method == "lme") {
        rootDev = lme(lmeFormula, data = parms[groupingFactor %in% 
                                                 y, ], random = randomFormula, correlation = R, 
                      na.action = na.omit)$logLik
      }
      else if (method == "nlme") {
        rootDev = nlme(model = nlme.model, fixed = fixedFormula, 
                       data = parms[groupingFactor %in% y, ], random = randomFormula, 
                       correlation = R, na.action = na.omit, start = start, 
                       groups = group)$logLik
      }
      if (continuous) {
        
        #for (i in xUnique) {
        clusterExport(cl=CL, list("y", "x", "parms", "nlme", "lme", "lrp_tvp"), 
                      envir = environment())
        
        dev = parSapply(CL, xUnique, 
                        
                        function(i){
                          
                          
                          yLeft = y[x <= i]
                          yRight = y[x > i]
                          if (length(yLeft) < control$minbucket || length(yRight) < 
                              control$minbucket) {
                            dev = c(dev, 0)
                          }
                          else {
                            if (method == "lme") {
                              modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                                             yLeft, ], random = randomFormula, correlation = R, 
                                                  na.action = na.omit), silent = TRUE)
                              modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                                              yRight, ], random = randomFormula, correlation = R, 
                                                   na.action = na.omit), silent = TRUE)
                            }
                            else if (method == "nlme") {
                              modelLeft = try(nlme(model = nlme.model, 
                                                   fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                                        yLeft, ], random = randomFormula, correlation = R, 
                                                   na.action = na.omit, start = start, groups = group), 
                                              silent = TRUE)
                              modelRight = try(nlme(model = nlme.model, 
                                                    fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                                         yRight, ], random = randomFormula, correlation = R, 
                                                    na.action = na.omit, start = start, groups = group), 
                                               silent = TRUE)
                            }
                            if (any(class(modelLeft) == "lme") && any(class(modelRight) == 
                                                                      "lme")) {
                              dev = c(dev, modelLeft$logLik + modelRight$logLik)
                            }
                            else {
                              0
                            }
                          }
                        }
                        
                        
        )
        
        good = rep(0, length(x))
        for (i in 1:length(xUnique)) {
          good[x == xUnique[i]] = dev[i]
        }
        good = good[1:(length(good) - 1)]
        list(goodness = good + abs(rootDev) * (good != 0) * 
               2, direction = rep(-1, length(good)))
      }
      else {
        order = rep(0, length(xUnique))
        response = parms[, names(parms) == responseName]
        for (i in 1:length(xUnique)) {
          order[i] = mean(response[x == xUnique[i]], na.rm = TRUE)
        }
        dir = sort(order, index.return = TRUE)$ix
        # for (i in 1:(length(dir) - 1)) {
        
        clusterExport(cl=CL, list("y", "x", "parms", "nlme", "lme", "lrp_tvp"), 
                      envir = environment())
        
        dev = parSapply(CL, 1:(length(dir) - 1), 
                        
                        function(i){
                          
                          yLeft = y[x %in% dir[1:i]]
                          yRight = y[x %in% dir[(i + 1):length(dir)]]
                          if (length(yLeft) < control$minbucket || length(yRight) < 
                              control$minbucket) {
                            dev = c(dev, 0)
                          }
                          else {
                            if (method == "lme") {
                              modelLeft = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                                             yLeft, ], random = randomFormula, correlation = R, 
                                                  na.action = na.omit), silent = TRUE)
                              modelRight = try(lme(lmeFormula, data = parms[groupingFactor %in% 
                                                                              yRight, ], random = randomFormula, correlation = R, 
                                                   na.action = na.omit), silent = TRUE)
                            }
                            else if (method == "nlme") {
                              modelLeft = try(nlme(model = nlme.model, 
                                                   fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                                        yLeft, ], random = randomFormula, correlation = R, 
                                                   na.action = na.omit, start = start, groups = group), 
                                              silent = TRUE)
                              modelRight = try(nlme(model = nlme.model, 
                                                    fixed = fixedFormula, data = parms[groupingFactor %in% 
                                                                                         yRight, ], random = randomFormula, correlation = R, 
                                                    na.action = na.omit, start = start, groups = group), 
                                               silent = TRUE)
                            }
                            if ((any(class(modelLeft) == "lme") | any(class(modelLeft) == 
                                                                      "nlme")) && (any(class(modelRight) == "lme") | 
                                                                                   any(class(modelRight) == "nlme"))) {
                              dev = c(dev, modelLeft$logLik + modelRight$logLik)
                            }
                            else {
                              0
                            }
                          }
                        }
        )
        list(goodness = dev + abs(rootDev) * (dev != 0) * 
               2, direction = dir)
      }
    }
    
    
  }
  
  
  
  initialize <- function(y, offset, parms = 0, wt) {
    list(y = y, parms = parms, numresp = 1, numy = 1, summary = function(yval, 
                                                                         dev, wt, ylevel, digits) {
      paste("deviance (-2logLik)", format(signif(dev), 
                                          3), "slope", signif(yval, 2))
    }, text = function(yval, dev, wt, ylevel, digits, n, 
                       use.n) {
      if (!use.n) {
        paste("m:", format(signif(yval, 1)))
      } else {
        paste("n:", n)
      }
    })
  }
  model <- list()
  model.rpart = rpart(paste(groupingName, c(rPartFormula)), 
                      method = list(eval = evaluation, split = split, init = initialize), 
                      control = control, data = data, parms = data)
  
  ### If Parallel, close clusters.
  if(use_parallel == TRUE){
    stopCluster(CL)
  }
  
  model$rpart_out <- model.rpart
  if (method == "lme") {
    model$lmeModel = lme(lmeFormula, data = data, random = randomFormula, 
                         correlation = R, na.action = na.omit)
    model$fixedFormula = lmeFormula
    model$lmeFormula = lmeFormula
  }
  else if (method == "nlme") {
    model$nlmeModel <- nlme(model = nlme.model, fixed = fixedFormula, 
                            data = data, random = randomFormula, correlation = R, 
                            na.action = na.omit, start = start, groups = group)
    model$fixedFormula <- fixedFormula
    model$nlme.model <- nlme.model
  }
  frame <- model.rpart$frame
  node2 = row.names(frame[frame[, "var"] == "<leaf>", ])
  model$node2 = node2
  model$leaf_node <- model.rpart$where
  summary = fixed_effects = var.corr = resid.var = list()
  
  data$resids = matrix(NA,nrow(data),1)
  data$fixed.fitted = matrix(NA,nrow(data),1)
  data$complete.fitted = matrix(NA,nrow(data),1)
  data$error = matrix(NA,nrow(data),1)
  
  
  for (j in 1:length(table(model.rpart$where))) {
    id <- names(table(model.rpart$where))[j] == model.rpart$where
    if (method == "lme") {
      model.out = lme(lmeFormula, data = data[id, ], random = randomFormula, 
                      correlation = R, na.action = na.omit)
      summary[[node2[j]]] <- summary(model.out)
      fixed_effects[[node2[j]]] <- fixed.effects(model.out)
      var.corr[[node2[j]]] <- model.out$modelStruct[[1]]
      resid.var[[node2[j]]] <- model.out$sigma
    }
    else if (method == "nlme") {
      model.out <- nlme(model = nlme.model, fixed = fixedFormula, 
                        data = data[id, ], random = randomFormula, correlation = R, 
                        na.action = na.omit, start = start, groups = group)
      summary[[node2[j]]] <- summary(model.out)
      fixed_effects[[node2[j]]] <- fixed.effects(model.out)
      var.corr[[node2[j]]] <- model.out$modelStruct[[1]]
      resid.var[[node2[j]]] <- model.out$sigma
      
      data[id,]$fixed.fitted = model.out$fitted[,1]
      data[id,]$complete.fitted = model.out$fitted[,2]
      data[id,]$resids = model.out$residuals[,2]
      data[id,]$error = model.out$residuals[,1]
      
    }
  }
  
  
  model$summary <- summary
  model$fixed_effects <- fixed_effects
  model$var.corr <- var.corr
  model$resid.var <- resid.var
  model$rpart_out <- model.rpart
  model$randomFormula = randomFormula
  model$R = R
  if (method == "nlme") {
    model$group = group
    model$start = start
  }
  model$method = method
#  model$data = data
  model$groupingName = groupingName
  model$rPartFormula = rPartFormula
  model$call <- match.call()
  class(model) <- "lrp_tvp"
  
  # Time-varying predictors predicting residuals
  
  tvp.rpart = rpart(paste("resids", c(time_var_covs)), data = data)
  
  #predicted values
  data$pred.resids.tvp = predict(tvp.rpart, data = data)
  data$final.tiv.tvp.pred = data$pred.resids.tvp + data$fixed.fitted

  model$data = data
  
  model$tvp.rpart = tvp.rpart
  
  model
}
