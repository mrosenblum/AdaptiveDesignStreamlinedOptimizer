### Streamlined Design Optimization Pipeline ##########
# Authors: Josh Betz (jbetz@jhu.edu) and Michael Rosenblum
# 
# 

### Fixed Parameters ###########################################################
# Default Parameter Values - ## CONFIRM ## That they are reasonable
library(knitr)
default.pct.of.max <- 0.8 # Start at 80% of max sample size/enrollment period

min.n.per.arm <- 25       # For Continuous/Binary Outcomes
min.enrollment.period <- 0.5    # For Survival Outcomes

default.function.scale <- 1
default.n.scale <- 100
default.period.scale <- 2
default.max.iterations <- 2 # Use for testing
# default.max.iterations <- 5e4 # Use for production
default.n.simulations <- 1e4
default.means.temperature <- 100
default.survival.temperature <- 10
default.evals.per.temp <- 10
default.report.iteration <- 1
default.power.penalty <- 1000
default.boundary.to.enroll <- 1

#code.dir <- "."
#optimizer.file <- "design_optimizer.R"
#performance.file <- "ComputePerformanceMetrics.R"
#binary.search.file <- "Utility_BinarySearch.R"
#backend.1tvc.file <- "Backend2Populations1Arm.R"
#backend.2tvc.file <- "Backend2Populations2Arms.R"

# Read in bash arguments from command line - extract parameters
bash.args <- commandArgs(trailingOnly=TRUE)
if(length(bash.args)>0) {
  for(i in 1:length(bash.args)) eval(parse(text=bash.args[i]))
}

# Read in optimizer code
source(file.path(code.dir, optimizer.file))
source(file.path(code.dir, binsearch.file))
source(file.path(code.dir, performance.file))

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
ui.subpopulation.sizes <- c(ui.subpopulation.1.size, 1-ui.subpopulation.1.size)
# If random seed is supplied, specify seeds. Otherwise pseudorandom seeds
# are chosen based on the initial RNG state.
if(!exists("initial.seed")){
  initial.seed <- sample(x=1:1e8, size=1)
}

# Set random seed
set.seed(initial.seed)

if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  n.stages <- 1 # Single Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations
  if(n.arms==2){
    source(file.path(code.dir, backend.1tvc.file))
    if(ui.type.of.outcome.data=="binary") {
      ui.outcome.mean <- subset(ui.population.parameters,select=c(2,4,1,3))
      ui.outcome.sd <- ui.outcome.mean*(1-ui.outcome.mean)
    } else{
      ui.outcome.mean <- cbind(array(0,c(nrow(ui.population.parameters),2)),subset(ui.population.parameters,select=c(1,2)))
      ui.outcome.sd <- subset(ui.population.parameters,select=c(4,6,3,5))
    }
  } else if(n.arms==3){
    source(file.path(code.dir, backend.2tvc.file))
    if(ui.type.of.outcome.data=="binary") {
      ui.outcome.mean <- ui.population.parameters
      ui.outcome.sd <- ui.outcome.mean*(1-ui.outcome.mean)
    } else{
      ui.outcome.mean <- subset(ui.population.parameters,select=c(1,3,5,7,9,11))
      ui.outcome.sd <- subset(ui.population.parameters,select=c(2,4,6,8,10,12))
    }
    
  }

  # Computes distribution of test statistics in a given scenario,
  # using canonical joint distribution
  get.z.distribution <- 
    function(...){
      construct.test.statistics.joint.distribution(...)
    }
  # Computes efficacy stopping boundaries
  get.efficacy.dunnett <-
    function(...){
      get.eff.bound(...)
    }
  # Evaluates performance of simulated trials 
  evaluate.design.dunnett <- 
    function(...){
      design.evaluate(...)
    }
  

  arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
  colnames(ui.outcome.sd) <- colnames(ui.outcome.mean) <-
    paste0(rep(arm.names, each=n.subpopulations),
           rep(1:n.subpopulations, n.arms))
  
  # if(ui.optimization.target=="ESS") {
  #   Switch Objective Function and Parameters
  # }
  
  max.enrollment.period <- (ui.max.duration-ui.followup.length)
  max.possible.accrual <- ui.accrual.yearly.rate*max.enrollment.period
  
  
  osea.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=default.pct.of.max*pmin(ui.max.size, max.possible.accrual)),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x) 
                    ceiling(
                      squash(x, 
                             min.n.per.arm,
                             pmin(ui.max.size, max.possible.accrual)))
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
                                      n.simulations=default.n.simulations,
                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                           number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=default.n.scale,
                max.iterations=default.max.iterations,
                temperature=default.means.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
  
  
} else { # Survival Cases
  n.stages <- 1 # Single Stage
  ui.hazard.rate <- ui.population.parameters
  number.of.alpha.allocation.components <- n.stages*n.subpopulations
  if(ui.include.designs.start.subpop.1){
    number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}
    
  if(n.arms==2){
    source(file.path(code.dir, backend.1tvc.file))
  } else if(n.arms==3){
    source(file.path(code.dir, backend.2tvc.file))
  }
  

  # Computes distribution of test statistics in a given scenario,
  # using canonical joint distribution
  get.z.distribution <- 
    function(...){
      construct.test.statistics.joint.distribution(...)
    }
  # Computes efficacy stopping boundaries
  get.efficacy.dunnett <-
    function(...){
      get.eff.bound(...)
    }
  # Evaluates performance of simulated trials 
  evaluate.design.dunnett <- 
    function(...){
      design.evaluate(...)
    }
  
  arm.names <- c(LETTERS[3], LETTERS[1:n.arms][-3])[1:n.arms]
  colnames(ui.hazard.rate) <- 
    paste0(rep(arm.names, each=n.subpopulations),
           rep(1:n.subpopulations, n.arms))

  osea.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=default.pct.of.max*
                         pmin(ui.max.duration,
                              ui.max.size/ui.accrual.yearly.rate)),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           pmin(ui.max.duration,
                                ui.max.size/ui.accrual.yearly.rate))),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=ui.max.duration,
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      n.simulations=default.n.simulations,
                                      alpha.allocation=
                                        rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=default.period.scale,
                max.iterations=default.max.iterations,
                temperature=default.survival.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
}

## 1SOA 1 stage optimized alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  osoa.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=default.pct.of.max*
                         pmin(ui.max.size, max.possible.accrual),
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
                             pmin(ui.max.size, max.possible.accrual))),
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
                                      n.simulations=default.n.simulations,
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=c(default.n.scale,rep(1,number.of.alpha.allocation.components)),
                max.iterations=default.max.iterations,
                temperature=default.means.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
  
  
} else { # Survival Cases
  osoa.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=default.pct.of.max*
                         pmin(ui.max.duration,
                              ui.max.size/ui.accrual.yearly.rate),
                       alpha.allocation=
                         rep(1/number.of.alpha.allocation.components,
                             number.of.alpha.allocation.components)
                       ),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           pmin(ui.max.duration,
                                ui.max.size/ui.accrual.yearly.rate)),
                    alpha.allocation=reals.to.probability
                    ),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=ui.max.duration,
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=NULL,
                                      n.simulations=default.n.simulations,
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=c(default.period.scale,rep(1,number.of.alpha.allocation.components)),
                max.iterations=default.max.iterations,
                temperature=default.survival.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
}

## 2SEA 2 stage equal alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
  n.stages <- 2 # Two Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations

  # if(ui.optimization.target=="ESS") {
  #   Switch Objective Function and Parameters
  # }
  
  tsea.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=default.pct.of.max*
                         pmin(ui.max.size, max.possible.accrual)),
                search.transforms=
                  # Cap sample size at minimum of the maximum specified size
                  # and the accrual rate x maximum allowable duration 
                  list(n.per.arm=function(x) 
                    ceiling(
                      squash(x, 
                             min.n.per.arm,
                             pmin(ui.max.size, max.possible.accrual)))
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
                                      n.simulations=default.n.simulations,
                                      alpha.allocation=rep(1/number.of.alpha.allocation.components,
                                                           number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=default.n.scale,
                max.iterations=default.max.iterations,
                temperature=default.means.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
  
  
} else { # Survival Cases
  n.stages <- 2 # Single Stage
  number.of.alpha.allocation.components <- n.stages*n.subpopulations
  if(ui.include.designs.start.subpop.1){
      number.of.alpha.allocation.components <- number.of.alpha.allocation.components - (n.subpopulations-1)}

  
    tsea.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=default.pct.of.max*
                         pmin(ui.max.duration,
                              ui.max.size/ui.accrual.yearly.rate)),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           pmin(ui.max.duration,
                                ui.max.size/ui.accrual.yearly.rate))),
                fixed.parameters=list(n.arms=n.arms,
                                      accrual.rate=ui.accrual.yearly.rate,
                                      subpopulation.sizes=ui.subpopulation.sizes,
                                      outcome.type='survival',
                                      non.inferiority=ifelse(ui.time.to.event.trial.type=="non-inferiority",TRUE,FALSE),
                                      hazard.rate=ui.hazard.rate,
                                      time=c(ui.max.duration/2,ui.max.duration),
                                      max.follow=Inf,
                                      censoring.rate=ui.time.to.event.censoring.rate,
                                      ni.margin=ui.time.to.event.non.inferiority.trial.margin,
                                      restrict.enrollment=FALSE,
                                      mcid=ui.mcid,
                                      futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                                      n.simulations=default.n.simulations,
                                      alpha.allocation=
                                        rep(1/number.of.alpha.allocation.components,
                                            number.of.alpha.allocation.components),
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=default.period.scale,
                max.iterations=default.max.iterations,
                temperature=default.survival.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
}

## 2SOA 2 stage optimized alpha
if(ui.type.of.outcome.data!="time-to-event"){ # Continuous and Binary Cases
   tsoa.result <-
    sa.optimize(search.parameters=
                  list(n.per.arm=default.pct.of.max*
                         pmin(ui.max.size, max.possible.accrual),
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
                             pmin(ui.max.size, max.possible.accrual))),
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
                                      n.simulations=default.n.simulations,
                                      total.alpha=ui.total.alpha
                                      ),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=c(default.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                max.iterations=default.max.iterations,
                temperature=default.means.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
  
  
} else { # Survival Cases
  tsoa.result <-
    sa.optimize(search.parameters=
                  list(enrollment.period=default.pct.of.max*
                         pmin(ui.max.duration,
                              ui.max.size/ui.accrual.yearly.rate),
                       time=c(ui.max.duration/2,ui.max.duration),
                       futility.boundaries=rep(-3,(n.arms-1)*n.subpopulations),
                       alpha.allocation=
                         rep(1/number.of.alpha.allocation.components,
                             number.of.alpha.allocation.components)
                  ),
                search.transforms=
                  list(enrollment.period=function(x)
                    squash(x, min.enrollment.period,
                           pmin(ui.max.duration,
                                ui.max.size/ui.accrual.yearly.rate)),
                    time=function(t){t1 <- squash(t[1],0.01,ui.max.duration-0.01); t2<-squash(t[2],t1+0.01,ui.max.duration); return(c(t1,t2))}, 
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
                                      n.simulations=default.n.simulations,
                                      total.alpha=ui.total.alpha),
                create.object=dunnett.wrapper,
                evaluate.object=power.penalized.weighted.ess,
                function.scale=default.function.scale,
                parameter.scale=c(default.n.scale,rep(1,n.stages+(n.arms-1)*n.subpopulations+number.of.alpha.allocation.components)),
                max.iterations=default.max.iterations,
                temperature=default.survival.temperature,
                evals.per.temp=default.evals.per.temp,
                report.iteration=default.report.iteration,
                power.penalty=default.power.penalty,
                power.constraints=ui.desired.power)
}

#save(osea.result,osoa.result,tsea.result,tsoa.result,file=file.path(data.dir,"results","optimizer_output.rda"))
knit(file.path(code.dir, report.generator.file),output=file.path(data.dir,"results","trial_design_performance_report.tex"))
