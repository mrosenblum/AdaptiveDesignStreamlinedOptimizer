#######
# This code is added by Tianchen Qian, in order to calculate
# bias, variance, MSE, CI coverage, CI inflation factor of a given design,
# using the simulated test statistics.
#######

#######
# Major update log:
#     2017.07.18: added function prototype, implemented compute.performance().
#     2017.07.19: implemented compute.performance.under.single.scenario()
#                 to compute bias, variance, mse, ci.coverage.
#     2017.08.01: added ci.inflation.factor output in compute.performance.under.single.scenario().
#     2017.08.01: added ci.level argument in compute.performance() to allow user-specified confidence level.
#     2017.08.01: added n.treatment.excluding.control and n.population arguments in
#                 compute.performance() to allow for different number of treatments/populations.
#     2017.08.11: add argument restrict.enrollment in compute.performance() to allow
#                 only enrolling subpopulation 1 in the first stage. 
#                 ASSUMING ONLY TWO STAGES (HARD CODED FOR NOW)
#     2017.08.17: add argument convert.to.hazard.ratio.scale and ni.margin in compute.performance()
#                 to convert the estimator of log hazard ratio (minus log non-inferiority-margin) 
#                 to hazard ratio scale, when dealing with survival outcome. In this case,
#                 bias, variance, and mse will be computed on hazard ratio scale.
#                 (ci.coverage and ci.inflation.factor are still on log hazard ratio scale.)
#     2017.09.11: now also handles one stage trials
#######

#######
# Note:
#  1. I use the word "scenario" to denote a arm*population combination.
#     So we are considering L * J scenarios: L arms * J populations.
#  2. The ordering in the test statistics is: "arm1.pop1", ..., "arm1.popJ", ..., "armL.pop1", ..., "armL.popJ".
#     This is the same ordering as in 2populations2arms.R.
#  3. When estimating survival outcome, hazard ratio is understood as (hazard.rate.treatment/hazard.rate.control).
#######


# This function computes performance measures of simulated trials under a single scenario
# (for a single treatment arm and a single population).
# Input:
#     test.statistics: a matrix of test statistics with dimension: n.sim * n.stages.
#     cov.var.matrix: the (true) covariance matrix of the test statistics.
#     stage.decision: a vector of length n.sim, indicating at which stage the trial stops for
#         this single scenario.
#     information.level: a vector of information levels at all K stages. That is, the inverse 
#         of variance of the estimators at all K stages.
#     true.estimand.value: the true estimand value (the basis of bias calculation)
#     ci.level: the confidence level of confidence interval (used in computing ci.coverage and ci.inflation.factor).
#         Default ci.level = 0.95.
#     restrict.enrollment: if TRUE, only enroll subpopulation 1 in the first stage. In this case
#         subpopulation 2 may not be enrolled even at the second stage (in which case stage.decision = 0).
#     convert.to.hazard.ratio.scale: if TRUE (only in the case of survival outcomes), convert estimators and
#         estimand from log hazard ratio (minus log ni.margin) to hazard ratio scale.
#         Bias, variance, and mse will be computed on hazard ratio scale.
#         (ci.coverage and ci.inflation.factor are still on log hazard ratio scale.)
#     ni.margin: the non-inferiority margin. Used to convert to hazard ratio scale for survival outcome with
#         non-inferiority trials.
# Output: A list consisting of four elements:
#     bias:        a number.
#     variance:    a number.
#     mse:         a number.
#     ci.coverage: a number.
#     ci.inflation.factor: a number, the inflation factor for ci to have correct coverage.

compute.performance.under.single.scenario <- function(
    test.statistics,
    cov.var.matrix,
    stage.decision,
    information.level,
    true.estimand.value,
    ci.level = 0.95,
    restrict.enrollment,
    convert.to.hazard.ratio.scale = FALSE,
    ni.margin = NULL
) {
    
    # for single stage trials, convert test.statistics from vector to matrix
    single.stage.trial <- is.vector(test.statistics)
    if (single.stage.trial) {
        test.statistics <- matrix(test.statistics, ncol = 1)
        cov.var.matrix <- matrix(cov.var.matrix, ncol = 1)
    }
    
    # checks on dimensions
    stopifnot(nrow(test.statistics) == length(stage.decision))
    stopifnot(ncol(test.statistics) == nrow(cov.var.matrix))
    stopifnot(ncol(test.statistics) == length(information.level))
    stopifnot(length(true.estimand.value) == 1)
    
    # check on valid argument values
    stopifnot((ci.level > 0 & ci.level <= 1))
    
    # check that stage.decision can be 0 only if restrict.enrollment = TRUE
    stopifnot((!is.null(restrict.enrollment) && restrict.enrollment) | all(stage.decision > 0))
    enrolled.sim.index <- which(stage.decision > 0)
    test.statistics <- test.statistics[enrolled.sim.index, ]
    stage.decision <- stage.decision[enrolled.sim.index]
    
    # number of MC evaluations (excluding simulations where the corresponding subpopulation is never enrolled)
    if (single.stage.trial) {
        n.sim <- length(test.statistics)
    } else {
        n.sim <- nrow(test.statistics)
    }
    
    # diagonal matrix of estimator variances
    if (single.stage.trial) {
        D <- matrix(information.level^(-1), ncol = 1, nrow = 1)
    } else {
        D <- diag(information.level^(-1))
    }
    
    # obtain estimator at the end of each simulated trial
    if (single.stage.trial) {
        estimators <- matrix(test.statistics, ncol = 1)
    } else {
        estimators <- test.statistics
    }
    for (istage in 1:nrow(D)) {
        # using for loop to avoid the impact of NA in stage 1 on stage 2, when restrict.enrollment = TRUE
        estimators[, istage] <- estimators[, istage] * D[istage, istage]^(1/2)
    }
    estimators.at.end.of.trial <- rep(NA, n.sim)
    for (i.sim in 1:n.sim) {
        estimators.at.end.of.trial[i.sim] <- estimators[i.sim, stage.decision[i.sim]]
    }
    
    # extract known variance of the estimator at the end of each trial
    variance.at.end.of.trial <- information.level[stage.decision]^(-1)
    
    # compute ci.coverage
    normal.quantile <- qnorm(1 - (1 - ci.level) / 2)
    ci.left     <- estimators.at.end.of.trial - normal.quantile * sqrt(variance.at.end.of.trial)
    ci.right    <- estimators.at.end.of.trial + normal.quantile * sqrt(variance.at.end.of.trial)
    ci.coverage <- mean((true.estimand.value >= ci.left) & (true.estimand.value <= ci.right))

    # compute ci.inflation.factor
    ci.inflation.factor <- tryCatch(
        {
            binary.search(function(inflation.factor) {
                ci.left.inflated <- estimators.at.end.of.trial - inflation.factor * normal.quantile * sqrt(variance.at.end.of.trial)
                ci.right.inflated <- estimators.at.end.of.trial + inflation.factor * normal.quantile * sqrt(variance.at.end.of.trial)
                ci.coverage.inflated <- mean((true.estimand.value >= ci.left.inflated) & (true.estimand.value <= ci.right.inflated))
                return(ci.coverage.inflated - ci.level)
            }, interval = c(0, 20))
        },
        error = function(cond) {
            message("\nCatched error in binary.search() when computing ci.inflation.factor:")
            message(cond)
            message("\n")
            return("Cannot find ci.inflation.factor using binary.search()")
        })
    
    # compute bias, variance, mse
    
    if (convert.to.hazard.ratio.scale) {
        if (is.null(ni.margin)) {
            ni.margin <- 1
        }
        # convert estimator to hazard ratio scale (hazard.rate.treatment/hazard.rate.control)
        estimators.at.end.of.trial <- exp(log(ni.margin) - estimators.at.end.of.trial)
        true.estimand.value <- exp(log(ni.margin) - true.estimand.value)
    }
    
    bias     <- mean(estimators.at.end.of.trial - true.estimand.value)
    variance <- var(estimators.at.end.of.trial)
    mse      <- mean((estimators.at.end.of.trial - true.estimand.value)^2)
    
    return(list(bias = bias,
                variance = variance,
                mse = mse,
                ci.coverage = ci.coverage,
                ci.inflation.factor = ci.inflation.factor))
}

# An example of simulating 10 trials, each with 2 stages.
# This example is to showcase that the function runs.
# Not a reasonable example in a real trial setting.
if (0) {
    n.sim <- 100
    test.statistics <- cbind(rnorm(n.sim), rnorm(n.sim))
    cov.var.matrix <- diag(c(1, 1))
    stage.decision <- sample(c(1,2), size = n.sim, replace = T)
    information.level <- c(10, 20)
    true.estimand.value <- 0
    compute.performance.under.single.scenario(
        test.statistics, cov.var.matrix, stage.decision, information.level, true.estimand.value)
    compute.performance.under.single.scenario(
        test.statistics, cov.var.matrix, stage.decision, information.level, true.estimand.value, ci.level = 0.80)
    compute.performance.under.single.scenario(
        test.statistics, cov.var.matrix, stage.decision, information.level, true.estimand.value, ci.level = -1)
    
}






# This function computes performance measures of simulated trials under all scenarios
# (for all treatment comparisons and all populations).
# Input:
#     test.statistics: a matrix of test statistics with dimension: n.sim * (n.stages * n.arms * n.pops).
#         Order of the test statistics within a stage is given by
#         ("arm1.pop1", ..., "arm1.popJ", ..., "armL.pop1", ..., "armL.popJ").
#     cov.var.matrix: the (true) covariance matrix of the test statistics.
#     stage.decision: a matrix with dimension n.sim * (n.arms * n.pops), indicating at which stage the decision
#         is made, in the order of ("arm1.pop1", ..., "arm1.popJ", ..., "armL.pop1", ..., "armL.popJ").
#     information.level: a matrix with dimension: n.stages * (n.arms * n.pops). The k-th row is the
#         information levels (inverse of estimator's variance) for stage k, with the order
#         ("arm1.pop1", ..., "arm1.popJ", ..., "armL.pop1", ..., "armL.popJ").
#         This is the output "information.vector" from construct.test.statistics.joint.distribution()
#         in the file "Backend2Populations2Arms.R" on master branch.
#     true.estimand.value: a vector of the true estimand value (the basis of bias calculation),
#         in the order of ("arm1.pop1", ..., "arm1.popJ", ..., "armL.pop1", ..., "armL.popJ").
#     ci.level: the confidence level of confidence interval (used in computing ci.coverage and ci.inflation.factor).
#         Default = 0.95.
#     n.treatment.excluding.control: number of comparator treatments (excluding control). Default = 2.
#     n.population: number of subpopulations. Default = 2.
#     restrict.enrollment: if TRUE, only enroll subpopulation 1 in the first stage. In this case
#         subpopulation 2 may not be enrolled even at the second stage (in which case stage.decision = 0).
#     convert.to.hazard.ratio.scale: if TRUE (only in the case of survival outcomes), convert estimators and
#         estimand from log hazard ratio (minus log ni.margin) to hazard ratio scale.
#         Bias, variance, and mse will be computed on hazard ratio scale.
#         (ci.coverage and ci.inflation.factor are still on log hazard ratio scale.)
#     ni.margin: the non-inferiority margin. Used to convert to hazard ratio scale for survival outcome with
#         non-inferiority trials.
# Output: A list of lists (one list for each scenario), each list consisting of five elements:
#     bias, variance, mse, ci.coverage, ci.inflation.factor.

compute.performance <- function(
    test.statistics,
    cov.var.matrix,
    stage.decision,
    information.level,
    true.estimand.value,
    ci.level = 0.95,
    n.treatment.excluding.control = 2,
    n.population = 2,
    restrict.enrollment = FALSE,
    convert.to.hazard.ratio.scale = FALSE,
    ni.margin = NULL
) {
    # Hard coded modification when restrict.enrollment = TRUE
    # browser()
    if (!is.null(restrict.enrollment) && restrict.enrollment) {
        if (n.population != 2 | ncol(test.statistics) != 3) {
            stop("restrict.enrollment = TRUE only allows for 2 subpopulation 2 stage trials.")
        }
        test.statistics <- cbind(test.statistics[, 1], rep(NA, nrow(test.statistics)), test.statistics[, 2:3])
        cov.var.matrix <- cbind(cov.var.matrix[, 1], rep(NA, nrow(cov.var.matrix)), cov.var.matrix[, 2:3])
        cov.var.matrix <- rbind(cov.var.matrix[1, ], rep(NA, ncol(cov.var.matrix)), cov.var.matrix[2:3, ])
    }
    
    # messages dealing with survival outcome
    #if (convert.to.hazard.ratio.scale) {
    #    if (is.null(ni.margin)) {
    #        message("Superiority trial for survival outcome.\n",
    #                "    Performance metrics reported are on hazard ratio scale.\n",
    #                "    Estimators are internally converted from log(hazard ratio) to hazard ratio.\n",
    #                "    ci.coverage and ci.inflation.factor are still computed on log(hazard ratio) scale.")
    #    } else {
    #        message("Non-inferiority trial for survival outcome.\n",
    #                "    Performance metrics reported are on hazard ratio scale.\n",
    #                "    Estimators are internally converted from log(hazard ratio) - log(ni.margin) to hazard ratio.\n",
    #                "    ci.coverage and ci.inflation.factor are still computed on log(hazard ratio) scale.")
    #    }
    #}
    
    # number of scenarios (arm*population combinations)
    n.arm.pop <- n.treatment.excluding.control * n.population
    
    # number of stages
    n.stages <- ncol(test.statistics) / n.arm.pop
    
    # scenarios <- c("arm1.pop1", "arm1.pop2", "arm2.pop1", "arm2.pop2")
    scenarios <- as.vector(matrix(outer(
        paste0("arm", 1:n.treatment.excluding.control), paste0("pop", 1:n.population), paste, sep = "."
    ), nrow = 1, byrow = TRUE))
    
    # a list to store the results under all scenarios
    perf.metric.list <- list()
    for (i.scn in 1:length(scenarios)) {
        # an index to extract test statistics under a single scenario
        index <- seq(from = i.scn, by = n.arm.pop, length = n.stages)
        
        perf.metric <- compute.performance.under.single.scenario(
            test.statistics[, index],
            cov.var.matrix[index, index],
            stage.decision[, i.scn],
            information.level[, i.scn],
            true.estimand.value[i.scn],
            ci.level,
            restrict.enrollment = restrict.enrollment,
            convert.to.hazard.ratio.scale = convert.to.hazard.ratio.scale,
            ni.margin = ni.margin
        )
        
        # concatenate the result under single scenario to the overall list
        perf.metric.list <- c(perf.metric.list, list(perf.metric))
    }
    
    names(perf.metric.list) <- scenarios
    
    return(perf.metric.list)
}


# An example of simulating 10 trials, each with 2 stages, 2 comparator amrs, 2 populations
# This example is to showcase that the function runs.
# Not a reasonable example in a real trial setting.
if (0) {
    n.sim <- 100
    test.statistics <- cbind(rnorm(n.sim), rnorm(n.sim),
                             rnorm(n.sim), rnorm(n.sim),
                             rnorm(n.sim), rnorm(n.sim),
                             rnorm(n.sim), rnorm(n.sim))
    cov.var.matrix <- diag(rep(1, 8))
    stage.decision <- cbind(sample(c(1,2), size = n.sim, replace = T),
                            sample(c(1,2), size = n.sim, replace = T),
                            sample(c(1,2), size = n.sim, replace = T),
                            sample(c(1,2), size = n.sim, replace = T))
    information.level <- cbind(c(10, 20),
                               c(10, 20),
                               c(10, 20),
                               c(10, 20))
    true.estimand.value <- c(0, 1, 0, 1)
    compute.performance(
        test.statistics, cov.var.matrix, stage.decision, information.level, true.estimand.value)
    compute.performance(
        test.statistics, cov.var.matrix, stage.decision, information.level, true.estimand.value, ci.level = 0.80)
}
