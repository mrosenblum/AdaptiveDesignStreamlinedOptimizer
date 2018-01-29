### Description: #################################################################
# This code is added by Tianchen Qian, in order to estimate R^2_W given a trial data set.
# Definition of R^2_W is given below:
#    R^2_W = [ Var{E(Y|W,A=1) | A=1} + Var{E(Y|W,A=0) | A=0} ] / [ Var{Y | A=1} + Var{Y | A=0} ]
# see page 13 of http://biostats.bepress.com/jhubiostat/paper285/.
#
# Relationship between R^2_W and relative efficency is also given in the paper as a main theorem.


### Major update log: ############################################################
#     2017.08.17: implemented main function.
#     2017.08.28: added wrapper to output relative efficiency


### Note: ########################################################################
#  1. Treatment indicator must be binary (taking values 0 and 1).
#  2. In the calculation, we use cross-validated predicted outcome value. Details are below:
#     To get the predicted value of mean of Y setting A=0 for a given participant's baseline vector W_i,
#     use \hat{E}(Y|A=0,W_i) where \hat{E} is fit using the leave-one-out data set omitting participant i.
#     This is a form of cross-validation to avoid overly optimistic evaluation of prediction accuracy.
#     (W: baseline variable. A: treatment indicator. Y: outcome)



### estimate.relative.efficiency #################################################
# This function estimates the relative efficiency between an efficient adjusted estimator
# and the unadjusted estimator for a given dataset.
# Input:
#     dataset: the data frame of interest.
#     varname.treatment: variable name of the treatment indicator.
#     varname.baseline.variable: variable name of the baseline variable (could be vector)
#     varname.outcome: variable name of the outcome variable.
#     outcome.type: one of "continuous", "binary"
# Output:
#     The estimated relative efficency, which is a scalar >=1.

estimate.relative.efficiency <- function(
    dataset,
    varname.treatment,
    varname.baseline.variable,
    varname.outcome,
    outcome.type = c("continuous", "binary")
) {
    rsq <- estimate.rsquared(dataset, varname.treatment, varname.baseline.variable, varname.outcome, outcome.type)
    return(1 / (1 - rsq))
}

### estimate.rsquared ############################################################
# This function estimates R^2_W for a given dataset.
# Input:
#     dataset: the data frame of interest.
#     varname.treatment: variable name of the treatment indicator.
#     varname.baseline.variable: variable name of the baseline variable (could be vector)
#     varname.outcome: variable name of the outcome variable.
#     outcome.type: one of "continuous", "binary"
# Output:
#     The estimated rsquared, which is a scalar.

estimate.rsquared <- function(
    dataset,
    varname.treatment,
    varname.baseline.variable,
    varname.outcome,
    outcome.type = c("continuous", "binary")
) {
    
    match.arg(outcome.type)
    
    if (outcome.type == "continuous") {
        outcome.family <- "gaussian"
    } else if (outcome.type == "binary") {
        outcome.family <- "binomial"
    }
    
    # check that treatment indicator only takes value 0 and 1
    stopifnot(sort(unique(dataset[, varname.treatment])) == c(0, 1))
    
    reg.formula <- as.formula(paste(varname.outcome, "~", paste(varname.baseline.variable, collapse = "+")))
    
    data.a1 <- dataset[dataset[, varname.treatment] == 1, ]
    data.a0 <- dataset[dataset[, varname.treatment] == 0, ]
    
    data.a1$predicted.outcome <- NA
    data.a0$predicted.outcome <- NA
    
    # predict outcome using leave-one-out regression
    
    for (i.patient in 1:nrow(data.a1)) {
        fit <- glm(reg.formula, data = data.a1[-i.patient, ], family = outcome.family)
        data.a1$predicted.outcome[i.patient] <- predict(fit, newdata = data.a1[i.patient, ], type = "response")
    }
    
    for (i.patient in 1:nrow(data.a0)) {
        fit <- glm(reg.formula, data = data.a0[-i.patient, ], family = outcome.family)
        data.a0$predicted.outcome[i.patient] <- predict(fit, newdata = data.a0[i.patient, ], type = "response")
    }
    
    # compute R^2_W
    
    rsquared <- ( var(data.a1$predicted.outcome) + var(data.a0$predicted.outcome) ) /
        ( var(data.a1[, varname.outcome]) + var(data.a0[, varname.outcome]) )
    
    return(rsquared)
}

# examples

if (0) {
    
    set.seed(123)
    
    n <- 500
    W1 <- rnorm(n)
    W2 <- rnorm(n)
    A <- rbinom(n, 1, 0.5)
    error <- rnorm(n)
    
    # continuous outcome Y
    Y <- W1 + W2 + A + error
    data.example <- data.frame(W1 = W1, W2 = W2, A = A, Y = Y)
    
    estimate.rsquared(dataset = data.example,
                      varname.treatment = "A",
                      varname.baseline.variable = c("W1", "W2"),
                      varname.outcome = "Y",
                      outcome.type = "continuous") # true R^2_W = 2/3
    
    estimate.relative.efficiency(dataset = data.example,
                                 varname.treatment = "A",
                                 varname.baseline.variable = c("W1", "W2"),
                                 varname.outcome = "Y",
                                 outcome.type = "continuous")
    
    estimate.rsquared(dataset = data.example,
                      varname.treatment = "A",
                      varname.baseline.variable = "W1",
                      varname.outcome = "Y",
                      outcome.type = "continuous") # true R^2_W = 1/3
    
    estimate.relative.efficiency(dataset = data.example,
                                 varname.treatment = "A",
                                 varname.baseline.variable = "W1",
                                 varname.outcome = "Y",
                                 outcome.type = "continuous")
    
    # binary outcome Y
    expit <- function(x) {exp(x) / (1 + exp(x))}
    Y <- rbinom(n, 1, expit(W1 + W2 + A))
    data.example <- data.frame(W1 = W1, W2 = W2, A = A, Y = Y)
    
    estimate.rsquared(dataset = data.example,
                      varname.treatment = "A",
                      varname.baseline.variable = c("W1", "W2"),
                      varname.outcome = "Y",
                      outcome.type = "binary")
    
    estimate.relative.efficiency(dataset = data.example,
                                 varname.treatment = "A",
                                 varname.baseline.variable = c("W1", "W2"),
                                 varname.outcome = "Y",
                                 outcome.type = "binary")
}