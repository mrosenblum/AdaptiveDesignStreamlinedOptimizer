### Streamlined Design Optimization Pipeline ##########
# Authors: Josh Betz (jbetz@jhu.edu) and Michael Rosenblum
library(knitr)
library(mvtnorm)
library(plyr)
### Fixed Parameters ###########################################################
min.n.per.arm <- 25       # For Continuous/Binary Outcomes
min.enrollment.period <- 0.5    # For Survival Outcomes
simulated.annealing.parameter.function.scale <- 1
simulated.annealing.parameter.n.scale <- 100
simulated.annealing.parameter.period.scale <- 2
simulated.annealing.parameter.max.iterations <- 1000
# simulated.annealing.parameter.max.iterations <- 5e4 # Use for production
simulated.annealing.parameter.n.simulations <- 1e4
simulated.annealing.parameter.means.temperature <- 100
simulated.annealing.parameter.survival.temperature <- 10
simulated.annealing.parameter.evals.per.temp <- 10
simulated.annealing.parameter.report.iteration <- 1
simulated.annealing.parameter.power.penalty <- 100000
simulated.annealing.parameter.boundary.to.enroll <- 1

#for local testing only:
#code.dir <- "."
#optimizer.file <- "design_optimizer.R"
#performance.file <- "ComputePerformanceMetrics.R"
#binsearch.file <- "Utility_BinarySearch.R"
#OneTreatmentArm.file <- "DesignClass2Subpopulations1TreatmentVsControl.R"
#TwoTreatmentArms.file <- "DesignClass2Subpopulations2TreatmentVsControl.R"

# Read in bash arguments from command line - extract parameters
bash.args <- commandArgs(trailingOnly=TRUE)
if(length(bash.args)>0) {
  for(i in 1:length(bash.args)) eval(parse(text=bash.args[i]))
}

# Read in optimizer code
source(file.path(code.dir, optimizer.file))
source(file.path(code.dir, binsearch.file))
source(file.path(code.dir, performance.file))
source(file.path(code.dir, OneTreatmentArm.file))
source(file.path(code.dir, TwoTreatmentArms.file))
source(file.path(code.dir, GroupSequentialOneTreatmentArm.file))

# Load parameters from user interface
load(file.path(data.dir, "parameters", "ui.parameters.rda"))


# Load parameters from user interface
#load(file.path(data.dir, "parameters", "ui.parameters.rda"))

# TWO ARM EXAMPLES
# Continuous Example
#load("michael_09252017_1335.rda"); ui.n.arms <- 2;
# Binary Example
#load("michael_09252017_1508.rda"); ui.n.arms <- 2;
# Survival Example
# load("michael_09252017_1551.rda")

# THREE ARM EXAMPLES
# Continuous Example
#load("michael_09252017_2052.rda")
# Binary Example
#load("michael_09252017_2100.rda")
# Survival Example
#load("michael_09252017_2112.rda")

# Get start time
isa.start.time <- proc.time()

### NOTE: restricted to two subpopulations ###
n.subpopulations <- 2 
n.arms <- ui.n.arms
if(n.arms==3){simulated.annealing.parameter.max.iterations <- 10} # Since 3 arm designs require substantially more computation time
ui.subpopulation.sizes <- c(ui.subpopulation.1.size, 1-ui.subpopulation.1.size)
# If random seed is supplied, specify seeds. Otherwise pseudorandom seeds
# are chosen based on the initial RNG state.
if(!exists("initial.seed")){
  initial.seed <- sample(x=1:1e8, size=1)
}

# Set random seed
set.seed(initial.seed)

# Source design.evaluation code corresponding to number of arms in trial
if(n.arms==2){
  # Computes distribution of test statistics in a given scenario,
  # using canonical joint distribution
  construct.joint.distribution.of.test.statistics <- 
    function(...){
      construct.joint.distribution.of.test.statistics.OneTreatmentArm(...)
    }
  # Computes efficacy stopping boundaries
  generate.efficacy.boundaries <-
    function(...){
      get.eff.bound.OneTreatmentArm(...)
    }
  # Evaluates performance of simulated trials 
  design.evaluate <- 
    function(...){
      design.evaluate.OneTreatmentArm(...)
    }
} else if(n.arms==3){
  # Computes distribution of test statistics in a given scenario,
  # using canonical joint distribution
  construct.joint.distribution.of.test.statistics <- 
    function(...){
      construct.joint.distribution.of.test.statistics.TwoTreatmentArms(...)
    }
  # Computes efficacy stopping boundaries
  generate.efficacy.boundaries <-
    function(...){
      get.eff.bound.TwoTreatmentArms(...)
    }
  # Evaluates performance of simulated trials 
  design.evaluate <- 
    function(...){
      design.evaluate.TwoTreatmentArms(...)
    }
}
# Set functions for computing design features and design evaluation
##
## Format User Inputs from Graphical User Interface
##

if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  if(n.arms==2){
    if(ui.type.of.outcome.data=="binary") {
      ui.outcome.mean <- subset(ui.population.parameters,select=c(2,4,1,3))
      ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
    } else{
      ui.outcome.mean <- cbind(array(0,c(nrow(ui.population.parameters),2)),subset(ui.population.parameters,select=c(1,2)))
      ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(4,6,3,5)))
    }
  } else if(n.arms==3){
    if(ui.type.of.outcome.data=="binary") {
      ui.outcome.mean <- ui.population.parameters
      ui.outcome.sd <- sqrt(ui.outcome.mean*(1-ui.outcome.mean))
    } else{
      ui.outcome.mean <- subset(ui.population.parameters,select=c(1,3,5,7,9,11))
      ui.outcome.sd <- sqrt(subset(ui.population.parameters,select=c(2,4,6,8,10,12)))
      }
  }
  arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
  colnames(ui.outcome.sd) <- colnames(ui.outcome.mean) <-
    paste0(rep(arm.names, each=n.subpopulations),
           rep(1:n.subpopulations, n.arms))
} else { # Survival Cases
    ui.hazard.rate <- ui.population.parameters
    if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}
    arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
    colnames(ui.hazard.rate) <- 
      paste0(rep(arm.names, each=n.subpopulations),
             rep(1:n.subpopulations, n.arms))
}

## Run optimizations, starting with 1 stage
n.stages <- 1 # Single Stage
number.of.alpha.allocation.components <- n.stages*n.subpopulations
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
 # if(ui.optimization.target=="ESS") {
  #   Switch Objective Function and Parameters
  # }
  max.enrollment.period <- (ui.max.duration-ui.followup.length)
  max.possible.accrual <- ui.accrual.yearly.rate*max.enrollment.period
  #Binary search to minimize sample size over feasible designs
  feasible.n.per.arm <- osea.result <- NULL
  n.per.arm.upper.bound <- min(ui.max.size, max.possible.accrual)/n.arms
  n.per.arm.lower.bound <- min.n.per.arm
  #Increase upper bound if necessary to meet power requirements
  #repeat{
  #  osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type=ui.type.of.outcome.data,
  #                                                      n.per.arm=floor(n.per.arm.upper.bound),
  #                                                      n.arms=n.arms,
  #                                                      accrual.rate=ui.accrual.yearly.rate,
  #                                                      delay=ui.followup.length,
  #                                                      subpopulation.sizes=ui.subpopulation.sizes,
  #                                                      interim.info.times=NULL,
  #                                                      outcome.mean=ui.outcome.mean,
  #                                                      outcome.sd=ui.outcome.sd,
  #                                                      mcid=ui.mcid,
  #                                                      futility.boundaries=NULL,
  #                                                      relative.efficiency=ui.relative.efficiency,        #                                                      n.simulations=simulated.annealing.parameter.n.simulations,
  #                                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
  #                                                                           number.of.alpha.allocation.components),
  #                                                      total.alpha=ui.total.alpha)
  #  discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
  #  feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
  #                                  discrepancy.between.desired.power.empirical.power<=0)
  #  if(feasibility.indicator){
  #    break
  #  }
  #  n.per.arm.upper.bound <- 2*n.per.arm.upper.bound
  #}
  while(n.per.arm.upper.bound-n.per.arm.lower.bound>0.1){
    candidate.n.per.arm <- mean(c(n.per.arm.lower.bound,n.per.arm.upper.bound))
    osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type=ui.type.of.outcome.data,
                                      n.per.arm=floor(candidate.n.per.arm),
                                      n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      delay=ui.followup.length,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      interim.info.times=NULL,
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency,
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                           number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate)
    discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
    feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
                                        discrepancy.between.desired.power.empirical.power<=0)
    if(feasibility.indicator==TRUE){#Current sample size is feasible; store results and explore smaller sample sizes
      #osea.result <- osea.design.performance.evaluation;
      feasible.n.per.arm <- candidate.n.per.arm;
      n.per.arm.upper.bound <- candidate.n.per.arm} else{ #Current sample size is infeasible; explore larger sample sizes
      n.per.arm.lower.bound <- candidate.n.per.arm  
    }
  }
  if(is.null(feasible.n.per.arm)){feasible.n.per.arm <- min(ui.max.size, max.possible.accrual)/n.arms} #If none feasible, consider max allowed sample size
  #Placeholder to force output into format expected by .Rnw file for report building
  osea.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=floor(feasible.n.per.arm)),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x){floor(feasible.n.per.arm)}),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      delay=ui.followup.length,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type=ui.type.of.outcome.data,
                                      interim.info.times=NULL,
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency,    
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                           number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=simulated.annealing.parameter.n.scale,
                max.iterations=2,
                temperature=simulated.annealing.parameter.means.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
  
  } else {#Survival Outcome 
  feasible.enrollment.period <- osea.design.performance.evaluation <- NULL
  enrollment.period.upper.bound <- min(ui.max.duration,
                                        ui.max.size/ui.accrual.yearly.rate)
  enrollment.period.lower.bound <- min.enrollment.period
  feasible.max.duration <- ui.max.duration
  #repeat{
  #  osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type='survival',
  #                                                        enrollment.period=enrollment.period.upper.bound,
  #                                                        n.arms=n.arms,
  #                                                        accrual.rate=ui.accrual.yearly.rate,
  #                                                        subpopulation.sizes=ui.subpopulation.sizes,
  #                                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
  #                                                        hazard.rate=ui.hazard.rate,
  #                                                        time=max(ui.max.duration,enrollment.period.upper.bound),
  #                                                        max.follow=Inf,
  #                                                        censoring.rate=ui.time.to.event.censoring.rate,
  #                                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
  #                                                        restrict.enrollment=FALSE,
  #                                                        mcid=ui.mcid,
  #                                                        futility.boundaries=NULL,
  #                                                        relative.efficiency=ui.relative.efficiency, 
  #                                                         n.simulations=simulated.annealing.parameter.n.simulations,
  #                                                        alpha.allocation=
  #                                                          rep(1/number.of.alpha.allocation.components,
  #                                                              number.of.alpha.allocation.components),
  #                                                        total.alpha=ui.total.alpha)
  #  discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
  #  feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
  #                                  discrepancy.between.desired.power.empirical.power<=0)
  #  if(feasibility.indicator){
  #    break
  #  }
  #  enrollment.period.upper.bound <- 2*enrollment.period.upper.bound
  #}
  while(enrollment.period.upper.bound-enrollment.period.lower.bound>0.01){
    candidate.enrollment.period <- mean(c(enrollment.period.lower.bound,enrollment.period.upper.bound))
    osea.design.performance.evaluation <- triage.based.on.outcome.type(outcome.type='survival',
                                      enrollment.period=candidate.enrollment.period,
                                      n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=max(ui.max.duration,candidate.enrollment.period),
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=
                                        rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate)
    discrepancy.between.desired.power.empirical.power <- max(ui.desired.power-cbind(osea.design.performance.evaluation$empirical.power,osea.design.performance.evaluation$conj.power),na.rm=TRUE)
    feasibility.indicator <- ifelse(is.na(discrepancy.between.desired.power.empirical.power),TRUE,
                                    discrepancy.between.desired.power.empirical.power<=0)
    if(feasibility.indicator){#Current sample size is feasible; store current results and explore smaller sample sizes
      #osea.result <- osea.design.performance.evaluation;
      feasible.enrollment.period <- candidate.enrollment.period;
      feasible.max.duration <- max(ui.max.duration,feasible.enrollment.period)
      enrollment.period.upper.bound <- candidate.enrollment.period} else{
        #Current sample size is infeasible; explore larger sample sizes
      enrollment.period.lower.bound <- candidate.enrollment.period  
    }
  }
  if(is.null(feasible.enrollment.period)){feasible.enrollment.period <- min(ui.max.duration,
                                                                            ui.max.size/ui.accrual.yearly.rate)} #if no feasible solution, use maximum duration
  #Placeholder to force output into format expected by .Rnw file for report building
  osea.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=feasible.enrollment.period),
                search.transforms=
                  list(enrollment.period=function(x){feasible.enrollment.period}),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=max(feasible.enrollment.period,ui.max.duration),
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=
                                        rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=simulated.annealing.parameter.period.scale,
                max.iterations=2,
                temperature=simulated.annealing.parameter.survival.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
}

## 1SOA 1 stage optimized alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  osoa.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                       alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components)
                       ),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x) 
                    ceiling(
                      squash(x, 
                             min.n.per.arm,
                             min(ui.max.size, max.possible.accrual)/n.arms)),
                    alpha.allocation=reals.to.probability
                                      ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      delay=ui.followup.length,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type=ui.type.of.outcome.data,
                                      interim.info.times=NULL,
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,number.of.alpha.allocation.components)),
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.means.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,                 
                optimization.target=ui.optimization.target)
  
} else { # Survival Cases
  osoa.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                         min(c(feasible.max.duration,
                              ui.max.size/ui.accrual.yearly.rate))),
                       alpha.allocation=
                         rep(1/number.of.alpha.allocation.components,
                             number.of.alpha.allocation.components)
                       ),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,feasible.enrollment.period),
                    alpha.allocation=reals.to.probability
                    ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=feasible.max.duration,
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=c(simulated.annealing.parameter.period.scale,rep(1,number.of.alpha.allocation.components)),
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.survival.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
}

## Two stage design
n.stages <- 2 # Two Stage
number.of.alpha.allocation.components <- n.stages*n.subpopulations

## 2SEA 2 stage equal alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  # if(ui.optimization.target=="ESS") {
  #   Switch Objective Function and Parameters
  # }
  tsea.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual)),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x) 
                    ceiling(
                      squash(x, 
                             min.n.per.arm,
                             min(ui.max.size, max.possible.accrual)/n.arms))
                  ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      delay=ui.followup.length,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      interim.info.times=c(1/2,1),
                                      outcome.type=ui.type.of.outcome.data,
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      mcid=ui.mcid,
                                      futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                           number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=simulated.annealing.parameter.n.scale,
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.means.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
} else { # Survival Cases
  if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}

  
    tsea.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                min(feasible.max.duration,
                                                     ui.max.size/ui.accrual.yearly.rate))),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           min(ui.max.duration,
                               ui.max.size/ui.accrual.yearly.rate))),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=c(feasible.max.duration/2,feasible.max.duration),
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      alpha.allocation=
                                        rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=simulated.annealing.parameter.period.scale,
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.survival.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
}

## 2SOA 2 stage optimized alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
   tsoa.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                       interim.info.times=c(1/2,1),
                       futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                       alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components)
                       ),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x) 
                    ceiling(
                      squash(x, 
                             min.n.per.arm,
                             min(ui.max.size, max.possible.accrual)/n.arms)),
                      interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                      alpha.allocation=reals.to.probability
                  ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      delay=ui.followup.length,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type=ui.type.of.outcome.data,
                                      outcome.mean=ui.outcome.mean,
                                      outcome.sd=ui.outcome.sd,
                                      mcid=ui.mcid,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate
                                      ),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.means.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
  
} else { # Survival Cases
  tsoa.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                min(feasible.max.duration,
                                                     ui.max.size/ui.accrual.yearly.rate)),
                       time=c(feasible.max.duration/2,feasible.max.duration),
                       futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                       alpha.allocation=
                         rep(1/number.of.alpha.allocation.components,
                             number.of.alpha.allocation.components)
                  ),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           min(ui.max.duration,
                               ui.max.size/ui.accrual.yearly.rate)),
                    time=function(t){t1 <- squash(t[1],min.enrollment.period,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))}, 
                    alpha.allocation=reals.to.probability
                    ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      relative.efficiency=ui.relative.efficiency, 
                                      n.simulations=simulated.annealing.parameter.n.simulations,
                                      total.alpha=ui.total.alpha,
                                      construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                      generate.efficacy.boundaries=generate.efficacy.boundaries,
                                      design.evaluate=design.evaluate),
                create.object=triage.based.on.outcome.type,
                evaluate.object=power.penalized.weighted,
                function.scale=simulated.annealing.parameter.function.scale,
                parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                max.iterations=simulated.annealing.parameter.max.iterations,
                temperature=simulated.annealing.parameter.survival.temperature,
                evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                report.iteration=simulated.annealing.parameter.report.iteration,
                scenario.weights=ui.scenario.weights,
                power.penalty=simulated.annealing.parameter.power.penalty,
                power.constraints=ui.desired.power,
                optimization.target=ui.optimization.target)
}

####
# Optimize 2 Stage, Group Sequential Design, for 2 arm trials
####
if(n.arms==2){
  # Computes distribution of test statistics in a given scenario,
  # using canonical joint distribution
  construct.joint.distribution.of.test.statistics <-
    function(...){
      construct.joint.distribution.of.test.statistics.GroupSequential.OneTreatmentArm(...)
    }
  # Computes efficacy stopping boundaries
  generate.efficacy.boundaries <-
    function(...){
      get.eff.bound.GroupSequential.OneTreatmentArm(...)
    }
  # Evaluates performance of simulated trials
  design.evaluate <-
    function(...){
      design.evaluate.GroupSequential.OneTreatmentArm(...)
    }
  ## Optimize:
  if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
    group.sequential.tsoa.result <-
      sa.optimize(search.parameters=
                    list(n.per.arm=ifelse(!is.null(feasible.n.per.arm),feasible.n.per.arm,max.possible.accrual),
                         interim.info.times=c(1/2,1),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                              number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    # Cap sample size at minimum of the maximum specified size
                    # and the accrual rate x maximum allowable duration
                    list(n.per.arm=function(x)
                      ceiling(
                        squash(x,
                               min.n.per.arm,
                               min(ui.max.size, max.possible.accrual)/n.arms)),
                      interim.info.times=function(x){c(squash(x[1],0.1,0.9),1)},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        delay=ui.followup.length,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type=ui.type.of.outcome.data,
                                        outcome.mean=ui.outcome.mean,
                                        outcome.sd=ui.outcome.sd,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate
                  ),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.means.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
    
  } else { # Survival Cases
    group.sequential.tsoa.result <-
      sa.optimize(search.parameters=
                    list(enrollment.period=ifelse(!is.null(feasible.enrollment.period),feasible.enrollment.period,
                                                  min(feasible.max.duration,
                                                      ui.max.size/ui.accrual.yearly.rate)),
                         time=c(feasible.max.duration/2,feasible.max.duration),
                         futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                         alpha.allocation=
                           rep(1/number.of.alpha.allocation.components,
                               number.of.alpha.allocation.components)
                    ),
                  search.transforms=
                    list(enrollment.period=function(x)
                      squash(x, min.enrollment.period,
                             min(ui.max.duration,
                                 ui.max.size/ui.accrual.yearly.rate)),
                      time=function(t){t1 <- squash(t[1],0.01,feasible.max.duration-0.02); t2<-squash(t[2],t1+0.01,feasible.max.duration); return(c(t1,t2))},
                      alpha.allocation=reals.to.probability
                    ),
                  fixed.parameters=list(n.arms=n.arms,
                                        accrual.rate=ui.accrual.yearly.rate,
                                        subpopulation.sizes=ui.subpopulation.sizes,
                                        outcome.type='survival',
                                        non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                        hazard.rate=ui.hazard.rate,
                                        max.follow=Inf,
                                        censoring.rate=ui.time.to.event.censoring.rate,
                                        ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                        restrict.enrollment=FALSE,
                                        mcid=ui.mcid,
                                        relative.efficiency=ui.relative.efficiency,
                                        n.simulations=simulated.annealing.parameter.n.simulations,
                                        total.alpha=ui.total.alpha,
                                        construct.joint.distribution.of.test.statistics=construct.joint.distribution.of.test.statistics,
                                        generate.efficacy.boundaries=generate.efficacy.boundaries,
                                        design.evaluate=design.evaluate),
                  create.object=triage.based.on.outcome.type,
                  evaluate.object=power.penalized.weighted,
                  function.scale=simulated.annealing.parameter.function.scale,
                  parameter.scale=c(simulated.annealing.parameter.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                  max.iterations=simulated.annealing.parameter.max.iterations,
                  temperature=simulated.annealing.parameter.survival.temperature,
                  evals.per.temp=simulated.annealing.parameter.evals.per.temp,
                  report.iteration=simulated.annealing.parameter.report.iteration,
                  scenario.weights=ui.scenario.weights,
                  power.penalty=simulated.annealing.parameter.power.penalty,
                  power.constraints=ui.desired.power,
                  optimization.target=ui.optimization.target)
  }
}

setwd(file.path(data.dir,"results"))
if(n.arms==2){
save(osea.result,osoa.result,tsea.result,tsoa.result,group.sequential.tsoa.result,file="optimizer_output.rda")} else {
  save(osea.result,osoa.result,tsea.result,tsoa.result,file="optimizer_output.rda")
}
knit(file.path(code.dir,report.generator.file),output="trial_design_performance_report.tex")
