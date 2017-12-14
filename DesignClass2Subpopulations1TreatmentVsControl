# Uses library(mvtnorm)

# Comparing a single treatment to control in two disjoint sub-populations

# This R file creates the necessary backend files for the optimizer to call
# There are three key functions 
# 1) construct.joint.distribution.of.test.statistics.OneTreatmentArm creates mean and covariance 
# matrices associated with the vector of statistics
# 2) get.eff.bound calculates the efficacy boundaries for the design
# 3) design.evaluate for different vectors of test statistics calculates
# which hypothesis are rejected and at which stage.

# Throughout the sequence of test statistics, is given by blocks of stages and within a block 
# of a stage k the vector of test statistics is given by 
# (Z_{1,k}, Z_{2,k}), where the first subscript
# indicates treatment and the second stage.



# This function calculates the covariate matrix for a binary and a continuous outcome
# It assumes one treatments and 
# a control, two sub-populations, and arbitrary number of stages
# prop.samp.vec.pop.1 
# Inputs: var.vec.pop.1: variance vector for population 1
#  var.vec.pop.2: variance vector for population 2
#  prop.samp.vec.pop.1: the proportion of total number of subjects in sub-population one
#  which are enrolled at each stage. E.g. c(0.5, 0.5) means 50% of the obs are enrolled 
#  at stage one and 50% at stage 2.
#  prop.samp.vec.pop.1: the proportion of total number of subjects in sub-population one
#  which are enrolled at each stage. E.g. c(0.5, 0.5) means 50 \% of the obs are enrolled 
#  at stage one and 50% at stage 2.

# Output: The covariance matrix associated with the vector of test statistics

cov.mat.cont.bin.OneTreatmentArm = function(var.vec.pop.1, var.vec.pop.2,
                            prop.samp.vec.pop.1, prop.samp.vec.pop.2){
  
  # K is the number of stages
  K = length(prop.samp.vec.pop.1)
  
  # Proportion of sample size enrolled up to and until stage k in each population
  cumsum.prop.samp.pop.1 = cumsum(prop.samp.vec.pop.1)
  cumsum.prop.samp.pop.2 = cumsum(prop.samp.vec.pop.2)
  cov.mat = matrix(0, nrow = 2 * K, ncol = 2 * K)
  
  # Filling in the covariance matrix
  
  # Filling in by blocks
  for(i in 1:K){
    for(j in 1:K){
      min.max.term.pop.1 =
        sqrt(min(cumsum.prop.samp.pop.1[i], cumsum.prop.samp.pop.1[j])/
               max(cumsum.prop.samp.pop.1[i], cumsum.prop.samp.pop.1[j]))
      
      min.max.term.pop.2 =
        sqrt(min(cumsum.prop.samp.pop.2[i], cumsum.prop.samp.pop.2[j])/
               max(cumsum.prop.samp.pop.2[i], cumsum.prop.samp.pop.2[j]))
      # (Z_{1,1}, Z_{2,1}, Z_{1,2}, Z_{2,2})
      cov.mat[(i-1)*2 + 1, (j-1)*2 + 1] = min.max.term.pop.1
      cov.mat[(i-1)*2 + 2, (j-1)*2 + 2] = min.max.term.pop.2
      #      cov.mat[(i-1)*2 + 3, (j-1)*2 + 3] = min.max.term.pop.1
      #      cov.mat[(i-1)*2 + 4, (j-1)*2 + 4] = min.max.term.pop.2
    }
  }
  return(cov.mat)
}


# This function calculates the covariance matrix associatec with 
# a survival outcome. 
# Input d.l.j, l = 0,1, and j = 1,2. Here, d.l.j is a vector of length
# K where element k is the expected number of deaths at or before analysis
# k in subpopulation j and treatment l.

# Output covariance matrix associate with vector of test statistics

cov.mat.surv.OneTreatmentArm = function(d.0.1, d.1.1, d.0.2, d.1.2){
  
  K = length(d.0.1)
  cov.mat = matrix(0, nrow = 2 * K, ncol = 2 * K)
  
  # Filling in by blocks
  for(i in 1:K){
    for(j in 1:K){
      min.max.term.pop.1.treatment.1 = sqrt((d.0.1[min(i,j)] + d.1.1[min(i,j)])/(d.0.1[max(i,j)] + d.1.1[max(i,j)]))
      min.max.term.pop.2.treatment.1 = sqrt((d.0.2[min(i,j)] + d.1.2[min(i,j)])/(d.0.2[max(i,j)] + d.1.2[max(i,j)]))
      
      # Both treatment 1 and subpopulation 1 different stages
      cov.mat[(i-1) * 2 + 1, (j-1) * 2 + 1] = min.max.term.pop.1.treatment.1
      
      # Both treatment 1 and subpopulation 2 different stages
      cov.mat[(i-1) * 2 + 2, (j-1) * 2 + 2] = min.max.term.pop.2.treatment.1
    }
  }
  # Throughout the sequence of test statistics is given by blocks of stages 
  # and within a block  of a stage k the vector of test statistics is given by 
  # (Z_{1,1,k}, Z_{1,2,k}, Z_{2,1,k}, Z_{2,2,k}), where the first subscript
  # indicates treatment, the second sub-populaiton and the third stage.
  
  return(cov.mat)
}



# This function calculates the covariance matrix associatec with 
# a survival outcome when just sub-population one is enrolled at beginning of trial. 
# Inputs: d.i.j a vector corresponding to the of expected number of deaths in treatment i 
#         and subpopulation j at each stage
# Output: Covarince matrix

cov.mat.surv.restrict.yes.OneTreatmentArm = function(d.0.1, d.1.1, d.0.2, d.1.2){
  
  # Creating covariance matrix (note we delete a row and a column later)
  K = length(d.0.1)
  cov.mat = matrix(0, nrow = 2 * K , ncol = 2 * K)
  
  # Filling in by blocks for subpopulation 1
  for(i in 1:K){
    for(j in 1:K){
      min.max.term.pop.1.treatment.1 = sqrt((d.0.1[min(i,j)] + d.1.1[min(i,j)])/(d.0.1[max(i,j)] + d.1.1[max(i,j)]))
      # Both treatment 1 and subpopulation 1 different stages
      cov.mat[(i-1) * 2 + 1, (j-1) * 2 + 1] = min.max.term.pop.1.treatment.1
    }
  }
  
  # Remove second row and second column as Z_{2,1} not available
  cov.mat = cov.mat[, -2]
  cov.mat = cov.mat[-2,]
  
  # Filling in by blocks for subpopulation 1
  for(i in 1:(K-1)){
    for(j in 1:(K-1)){
      min.max.term.pop.2.treatment.1 = sqrt((d.0.2[min(i,j)] + d.1.2[min(i,j)])/(d.0.2[max(i,j)] + d.1.2[max(i,j)]))
      # Both treatment 1 and subpopulation 2 different stages
      cov.mat[((i-1) * 2 + 3), ((j-1) * 2 + 3)] = min.max.term.pop.2.treatment.1
    }
  }
  
  return(cov.mat)
}


# This function creates the covariance matrix and mean vector 
# associated with the test statistic
# Inputs: #   analytic.n.per.stage - [K x J(L+1)] matrix: patients with primary outcome
#             at each interim analysis. For a survival outcome this should be total patients enrolled which 
#.            equals analytic.n.per.stage if dealy is set to zero
#             stage 1: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#             stage 2: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#                ...
#             stage k: T0S1 T0S2 ... T0SJ ... TLS1 TLS2 ... TLSJ
#         outcome.type: type of outcome, one of continuous, binary or survival
#         mean.sub.pop.1: the assumed means assocated with each treatment in sub-population 1
#         the mean vector is input in the order (control, treatment) with
#         mean.sub.pop.2: the assumed means assocated with each treatment in sub-population 2
#         the mean vector is input in the order (control, treatment) with
#         var.vec.pop.1: the variance vector associated with each treatment in sub-population one
#         var.vec.pop.2: the variance vector associated with each treatment in sub-population two
#         prop.pop.1: The proportion of subjects in population one. Assumed known.
#         NOTE: For a survival outcome the code is setup such that if the hazard rate in the treatment 
#         group is smaller than in the control group that leads to large value of the test statistics
#         and rejection of the null-hypothesis.
#         NI: A logical variable if a non-inferiority test is done
#         ni.margin: the non-inferiority margin. Only for a non-inferiority test.
#         NOTE: Non-inferiority only works for a survival outcome and the non-inferiority margin
#         is \lambda_{Treatment}/\lambda_{control} and is therefore greater than 1.
#         max.follow: For survival outcome, how long each participant is followed up for
#         enrollment.period: For survival outcome, the maximum time participants are enrolled
#         hazard.rate.pop.1: For survival outcome, hazard rate for subpopulation 1 in 
#         the order (control, treatment)
#         hazard.rate.pop.2: For survival outcome, hazard rate for subpopulation 2 in 
#         the order (control, treatment)
#         time: for survival outcomes is the timing of all analysis, only needs to be specified for
#         a survival outcome
#         restrict.enrollment: for a survival outcome, if enrollment at the start data of the trial 
#         is restricted only to sub-population 1
#         censoring.rate: For a survival outcome only. It is the proportion of participants that are not
#         administratively censored which drop out of the study. For example, if 100 events are expected
#         without any dropout then setting censoring rate to 0.5 means that 100*0.5 events are expected
#         relative.efficiency: ratio of the asymptotic variance of unadjusted estimator to 
#         asymptotic variance of adjusted estimator. Should be greater than one.



# Output: A list with three elements:
#         cov.mat.used: Covariance matrix associated with vector of test statistic.
#         non.centrality.parameter.vec = The mean vector associated with each test statistic.
#         information.vector =  for a given stage k elements [(1+ (k-1) * 2):(2 + (k-1) * 2)] are 
#         (var(\beta_{1,k}, var(\beta_{2,k})) where the
#         first subscript indicates sub-populaiton and the second stage. 
#         beta is the estimator of the treatment effect

# Connection to Betz et al. paper: "Comparison of Adaptive Randomized Trial Designs for Time-to-Event Outcomes that Expand Versus Restrict Enrollment Criteria, to Test Non-Inferiority Design Optimization"
#  The 3 design classes defined in Section 5.1 of that paper are:
#  1. D_{ONE-STAGE} [standard, single stage designs]
#  2. D_{ADAPTIVE,START-BOTH} [adaptive, starts by enrolling both subpops in stage 1, can restrict enrollment at interim analyses]
#  3. D_{ADAPTIVE,START-SUBPOP.1} [adaptive, starts by enrolling only subpop. 1; can initiate enrollment of subpop 2 based on decision at first interim analysis. At interim analyses after the first, the decisions are to continue enrollment or restrict.] 
#  To get the design class (1), use a single stage and set restrict.enrollment = FALSE.
#  To get the design class (2), set restrict.enrollment = FALSE.
#  To get the design class (3), set restrict.enrollment = TRUE.

construct.joint.distribution.of.test.statistics.OneTreatmentArm <- function(analytic.n.per.stage,
                                                         mean.sub.pop.1=NULL,
                                                         mean.sub.pop.2=NULL,
                                                         var.vec.pop.1=NULL,
                                                         var.vec.pop.2=NULL,
                                                         outcome.type,
                                                         prop.pop.1,
                                                         NI = FALSE,
                                                         ni.margin = NULL,
                                                         max.follow = NULL,
                                                         enrollment.period = NULL,
                                                         hazard.rate.pop.1 = NULL,
                                                         hazard.rate.pop.2 = NULL,
                                                         time = NULL,
                                                         restrict.enrollment = FALSE,
                                                         censoring.rate = NULL,
                                                         relative.efficiency = NULL){
  
  # Number of stages
  K <- nrow(analytic.n.per.stage)
  
  # Calculating the total number of subjects at each analysis for both sub-populations
  # Note: We assume that an equal number is enrolled to treatment and control.
  n.pop.1 = analytic.n.per.stage[, 1]
  n.pop.2 = analytic.n.per.stage[, 2]
  
  # Calculating the proportion of observation sampled at each stage for the
  # Two treatments
  prop.samp.vec.pop.1 = diff(c(0, n.pop.1))/n.pop.1[K]
  prop.samp.vec.pop.2 = diff(c(0, n.pop.2))/n.pop.2[K]
  
  
  # Creating storage space for mean vector
  mean.vec = rep(NA, 2 * K)
  
  # Do the calculations seperately depending on the type of outcome
  if(outcome.type == "continuous"){
    for(i in 1:K){
      mean.vec[((i-1)*2+1):((i-1)*2+2)] = 
        c(sqrt(n.pop.1[i])*(mean.sub.pop.1[2]-mean.sub.pop.1[1])/
            sqrt(var.vec.pop.1[2]+var.vec.pop.1[1]), 
          sqrt(n.pop.2[i])*(mean.sub.pop.2[2] - mean.sub.pop.2[1])/
            sqrt(var.vec.pop.2[2]+var.vec.pop.2[1]))    
    }
    
    cov.mat.used = cov.mat.cont.bin.OneTreatmentArm(var.vec.pop.1,
                                    var.vec.pop.2,
                                    prop.samp.vec.pop.1,
                                    prop.samp.vec.pop.2)
  }
  
  if(outcome.type == "binary"){
    var.vec.pop.1 = mean.sub.pop.1*(1 - mean.sub.pop.1)
    var.vec.pop.2 = mean.sub.pop.2*(1 - mean.sub.pop.2)
    
    for(i in 1:K){
      mean.vec[((i-1)*2+1):((i-1)*2+2)] = 
        c(sqrt(n.pop.1[i])*(mean.sub.pop.1[2]-mean.sub.pop.1[1])/
            sqrt(var.vec.pop.1[2]+var.vec.pop.1[1]), 
          sqrt(n.pop.2[i])*(mean.sub.pop.2[2] - mean.sub.pop.2[1])/
            sqrt(var.vec.pop.2[2]+var.vec.pop.2[1]))    
    }
    
    cov.mat.used = cov.mat.cont.bin.OneTreatmentArm(var.vec.pop.1,
                                    var.vec.pop.2,
                                    prop.samp.vec.pop.1,
                                    prop.samp.vec.pop.2)
  }
  
  
  # Create information vector on the extimator scale
  # for a given stage k elements [(1+ (k-1) * 4):(4 + (k-1) * 4)] are 
  # (var(\beta_{1,1,k}, var(\beta_{1,2,k}), var(\beta_{2,1,k}), var(beta_{2,2,k})
  #, where the first subscript indicates treatment, the second sub-populaiton and the third stage.
  
  # For a continuous outcome 
  if(outcome.type == "continuous"){
    # Initialize the vector
    information.vector.inv = rep(NA, 2 * K)
    
    for(i in 1:K){
      information.vector.inv[((i-1)*2+1):((i-1)*2+2)] = 
        c(1/n.pop.1[i] * (var.vec.pop.1[2]+var.vec.pop.1[1]), 
          1/n.pop.2[i] * (var.vec.pop.2[2]+var.vec.pop.2[1]))
    }
  } # End if outcome.type is continuous
  
  # For a binary outcome 
  if(outcome.type == "binary"){
    # Initialize the vector
    information.vector.inv = rep(NA, 2 * K)
    
    var.vec.pop.1 = mean.sub.pop.1*(1 - mean.sub.pop.1)
    var.vec.pop.2 = mean.sub.pop.2*(1 - mean.sub.pop.2)
    
    for(i in 1:K){
      information.vector.inv[((i-1)*2 + 1):((i-1)*2 + 2)] = 
        c(1/n.pop.1[i] * (var.vec.pop.1[2]+var.vec.pop.1[1]), 
          1/n.pop.2[i] * (var.vec.pop.2[2]+var.vec.pop.2[1]))
    }
  } # End if outcome.type is binary
  
  if(outcome.type != "survival"){
    # Make information vector fit Tianchens code
    information.vector.inv.matrix = matrix(information.vector.inv, nrow = K, byrow = TRUE)
  }
  
  
  
  if(outcome.type == "survival"){
    
    mean.sub.pop.1 = 1/hazard.rate.pop.1
    mean.sub.pop.2 = 1/hazard.rate.pop.2
    
    
    # Do modification for NI test
    if(NI == FALSE){
      # The ratios of log-rank tests
      theta = -c(log(hazard.rate.pop.1[2]/hazard.rate.pop.1[1]), log(hazard.rate.pop.2[2]/hazard.rate.pop.2[1]))
    }
    
    if(NI == TRUE){
      # The ratios of log-rank tests
      theta = -c(log(hazard.rate.pop.1[2]/hazard.rate.pop.1[1]), log(hazard.rate.pop.2[2]/hazard.rate.pop.2[1])) + log(ni.margin)
    }
    
    
    if(restrict.enrollment == FALSE){
      # non-centraility parameter in same order as test statistics mentioned above
      mean.vec = rep(NA, 2 * K)
      
      d.0.1 = rep(0,K) # number of deaths in control group pop 1
      d.1.1 = rep(0,K) # number of deaths in treatment group pop 1
      d.0.2 = rep(0,K) # number of deaths in control group pop 2
      d.1.2 = rep(0,K) # number of deaths in treatment group pop 2
      
      for(i in 1:K){
        
        # Calculating the expected number of deaths for each treatment + sub-population combination 
        # at interim analys i
        
        # We cycle through the 6 different cases
        if(enrollment.period >= time[i] & time[i] >= max.follow){
          d.0.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[1] * max.follow))/(max.follow * hazard.rate.pop.1[1]))
          d.1.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[2] * max.follow))/(max.follow * hazard.rate.pop.1[2]))
          d.0.2[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.2[1] * max.follow))/(max.follow * hazard.rate.pop.2[1]))
          d.1.2[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.2[2] * max.follow))/(max.follow * hazard.rate.pop.2[2]))
        }
        
        if(enrollment.period >= max.follow & max.follow >= time[i]){
          d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
          d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time[i]))/(time[i] * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time[i]))/(time[i] * hazard.rate.pop.2[2])
        }
        
        if(time[i] >= enrollment.period & enrollment.period >= max.follow){
          k = time[i] - enrollment.period
          d.0.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * k) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - k) * hazard.rate.pop.1[1]))
          d.1.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * k) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - k) * hazard.rate.pop.1[2]))
          d.0.2[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[1] * k) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - k) * hazard.rate.pop.2[1]))
          d.1.2[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[2] * k) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - k) * hazard.rate.pop.2[2]))
        }
        
        if(max.follow >= enrollment.period & enrollment.period >= time[i]){
          d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
          d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time[i]))/(time[i] * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time[i]))/(time[i] * hazard.rate.pop.2[2])
        }
        
        if(max.follow >= time[i] & time[i] >= enrollment.period){
          d.0.1[i] = 1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * time[i]))/(enrollment.period * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * time[i]))/(enrollment.period * hazard.rate.pop.1[2])
          d.0.2[i] = 1 - (exp(-hazard.rate.pop.2[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[1] * time[i]))/(enrollment.period * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (exp(-hazard.rate.pop.2[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[2] * time[i]))/(enrollment.period * hazard.rate.pop.2[2])
        }
        
        if(time[i] >= max.follow & max.follow >= enrollment.period){
          d.0.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[1]))
          d.1.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[2]))
          d.0.2[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.2[1]))
          d.1.2[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.2[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.2[2]))
        }
        
        d.0.1[i] = d.0.1[i] * n.pop.1[i] * (1 - censoring.rate)
        d.1.1[i] = d.1.1[i] * n.pop.1[i] * (1 - censoring.rate)
        d.0.2[i] = d.0.2[i] * n.pop.2[i] * (1 - censoring.rate)
        d.1.2[i] = d.1.2[i] * n.pop.2[i] * (1 - censoring.rate)
        
        # Calculating the information and covariance matrix
        mean.vec[((i-1) * 2 + 1):((i-1) * 2 + 2)] = theta * sqrt(c((d.0.1[i] + d.1.1[i])/4, (d.0.2[i] + d.1.2[i])/4))
      }
      
      cov.mat.used = cov.mat.surv.OneTreatmentArm(d.0.1, d.1.1, d.0.2, d.1.2)
      
      
      # Initialize the vector
      information.vector.inv = rep(NA, 2 * K)
      
      for(i in 1:K){
        information.vector.inv[((i-1)*2+1):((i-1)*2+2)] = 
          c(4/((d.0.1[i] + d.1.1[i])), 
            4/((d.0.2[i] + d.1.2[i])))
      }
      # Make information vector fit Tianchens code
      information.vector.inv.matrix = matrix(information.vector.inv, nrow = K, byrow = TRUE)
      
    } # End if restrict enrollment 
    
    # If enrollment at beginning of trial is restricted to sub-population 1
    if(restrict.enrollment == TRUE){
      
      # information vector in same order as the vector test statistics mentioned above.
      mean.vec = rep(NA, 2 * (K-1) + 1)
      
      # Initialize the vector
      information.vector.inv = rep(NA, 2 * K - 1)
      
      # Calculate the expected number of deaths in each treatment arm and subpopulations at each stage
      d.0.1 = rep(0,K)
      d.1.1 = rep(0,K)
      d.0.2 = rep(0,K-1)
      d.1.2 = rep(0,K-1)
      
      # Cycle through sub-population 1
      for(i in 1:K){
        
        # Calculating the expected number of deaths for each treatment in subpopulation 1 
        # at interim analys i
        
        # Note first we calcualte the probability for an individual observation and then at the end
        # multiply by the sample size to get the expectation
        
        # We cycle through the 6 different cases
        if(enrollment.period >= time[i] & time[i] >= max.follow){
          d.0.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[1] * max.follow))/(max.follow * hazard.rate.pop.1[1]))
          d.1.1[i] = (time[i] - max.follow)/time[i] * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + max.follow/time[i] * (1 - (1-exp(-hazard.rate.pop.1[2] * max.follow))/(max.follow * hazard.rate.pop.1[2]))
        }
        
        if(enrollment.period >= max.follow & max.follow >= time[i]){
          d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
        }
        
        if(time[i] >= enrollment.period & enrollment.period >= max.follow){
          k = time[i] - enrollment.period
          d.0.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * k) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - k) * hazard.rate.pop.1[1]))
          d.1.1[i] = (enrollment.period + k -max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (max.follow - k)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * k) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - k) * hazard.rate.pop.1[2]))
        }
        
        if(max.follow >= enrollment.period & enrollment.period >= time[i]){
          d.0.1[i] = 1 - (1-exp(-hazard.rate.pop.1[1] * time[i]))/(time[i] * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (1-exp(-hazard.rate.pop.1[2] * time[i]))/(time[i] * hazard.rate.pop.1[2])
        }
        
        if(max.follow >= time[i] & time[i] >= enrollment.period){
          d.0.1[i] = 1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * time[i]))/(enrollment.period * hazard.rate.pop.1[1])
          d.1.1[i] = 1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * time[i]))/(enrollment.period * hazard.rate.pop.1[2])
        }
        
        if(time[i] >= max.follow & max.follow >= enrollment.period){
          d.0.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[1] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[1] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[1] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[1]))
          d.1.1[i] = (time[i] - max.follow)/enrollment.period * (1 - exp(-hazard.rate.pop.1[2] * max.follow)) + (enrollment.period - time[i] + max.follow)/enrollment.period * (1 - (exp(-hazard.rate.pop.1[2] * (time[i]-enrollment.period)) - exp(-hazard.rate.pop.1[2] * max.follow))/((max.follow - time[i] + enrollment.period) * hazard.rate.pop.1[2]))
        }
        
        # Multiplying by sample size to get expecations
        d.0.1[i] = d.0.1[i] * n.pop.1[i] * (1 - censoring.rate)
        d.1.1[i] = d.1.1[i] * n.pop.1[i] * (1 - censoring.rate)
        
        # Calculating the information vector for sub-population1
        if(i == 1){
          mean.vec[1] = theta[1] * sqrt((d.0.1[i] + d.1.1[i])/4)
          information.vector.inv[1] = 4/((d.0.1[1] + d.1.1[1]))
        }
        if(i != 1){
          mean.vec[(i-2) *2 + 2] = theta[1] * sqrt((d.0.1[i] + d.1.1[i])/4)
          information.vector.inv[(i-2) *2 + 2] =  4/((d.0.1[i] + d.1.1[i]))
        }
      }
      
      
      
      # Cycle through sub-population 2
      for(i in 1:(K-1)){
        
        # Calculating the expected number of deaths for each treatment in subpopulation 1 
        # at interim analys i
        
        # Note first we calcualte the probability for an individula observation and then at the end
        # multiply by the sample size to get the expectation
        
        # As enrollemnt starts at time[1] lenght of study for subpopulation 1 is
        enrollment.period.2 = enrollment.period - time[1]
        # and time of interim analysis since start of enrollment is
        time.2 = (time - time[1])[-1]
        
        
        # We cycle through the 6 different cases described in the pdf file
        if(enrollment.period.2 >= time.2[i] & time.2[i] >= max.follow){
          d.0.2[i] = (time.2[i] - max.follow)/time.2[i] * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + max.follow/time.2[i] * (1 - (1-exp(-hazard.rate.pop.2[1] * max.follow))/(max.follow * hazard.rate.pop.2[1]))
          d.1.2[i] = (time.2[i] - max.follow)/time.2[i] * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + max.follow/time.2[i] * (1 - (1-exp(-hazard.rate.pop.2[2] * max.follow))/(max.follow * hazard.rate.pop.2[2]))
        }
        
        if(enrollment.period.2 >= max.follow & max.follow >= time.2[i]){
          d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time.2[i]))/(time.2[i] * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time.2[i]))/(time.2[i] * hazard.rate.pop.2[2])
        }
        
        if(time.2[i] >= enrollment.period.2 & enrollment.period.2 >= max.follow){
          k = time.2[i] - enrollment.period.2
          d.0.2[i] = (enrollment.period.2 + k -max.follow)/enrollment.period.2 * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (max.follow - k)/enrollment.period.2 * (1 - (exp(-hazard.rate.pop.2[1] * k) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - k) * hazard.rate.pop.2[1]))
          d.1.2[i] = (enrollment.period.2 + k -max.follow)/enrollment.period.2 * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (max.follow - k)/enrollment.period.2 * (1 - (exp(-hazard.rate.pop.2[2] * k) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - k) * hazard.rate.pop.2[2]))
        }
        
        if(max.follow >= enrollment.period.2 & enrollment.period.2 >= time.2[i]){
          d.0.2[i] = 1 - (1-exp(-hazard.rate.pop.2[1] * time.2[i]))/(time.2[i] * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (1-exp(-hazard.rate.pop.2[2] * time.2[i]))/(time.2[i] * hazard.rate.pop.2[2])
        }
        
        if(max.follow >= time.2[i] & time.2[i] >= enrollment.period.2){
          d.0.2[i] = 1 - (exp(-hazard.rate.pop.2[1] * (time.2[i]-enrollment.period.2)) - exp(-hazard.rate.pop.2[1] * time.2[i]))/(enrollment.period.2 * hazard.rate.pop.2[1])
          d.1.2[i] = 1 - (exp(-hazard.rate.pop.2[2] * (time.2[i]-enrollment.period.2)) - exp(-hazard.rate.pop.2[2] * time.2[i]))/(enrollment.period.2 * hazard.rate.pop.2[2])
        }
        
        if(time.2[i] >= max.follow & max.follow >= enrollment.period.2){
          d.0.2[i] = (time.2[i] - max.follow)/enrollment.period.2 * (1 - exp(-hazard.rate.pop.2[1] * max.follow)) + (enrollment.period.2 - time.2[i] + max.follow)/enrollment.period.2 * (1 - (exp(-hazard.rate.pop.2[1] * (time.2[i]-enrollment.period.2)) - exp(-hazard.rate.pop.2[1] * max.follow))/((max.follow - time.2[i] + enrollment.period.2) * hazard.rate.pop.2[1]))
          d.1.2[i] = (time.2[i] - max.follow)/enrollment.period.2 * (1 - exp(-hazard.rate.pop.2[2] * max.follow)) + (enrollment.period.2 - time.2[i] + max.follow)/enrollment.period.2 * (1 - (exp(-hazard.rate.pop.2[2] * (time.2[i]-enrollment.period.2)) - exp(-hazard.rate.pop.2[2] * max.follow))/((max.follow - time.2[i] + enrollment.period.2) * hazard.rate.pop.2[2]))
        }
        
        # Multiplying by sample size to get expecations
        # Need to adjust sample size to account for late enrollment
        d.0.2[i] = d.0.2[i] * n.pop.2[i+1] * (1 - censoring.rate)
        d.1.2[i] = d.1.2[i] * n.pop.2[i+1] * (1 - censoring.rate)
        
        mean.vec[((i-1) * 2 + 3)] = theta[2] * sqrt((d.0.2[i] + d.1.2[i])/4)
        
        # Calculating the information vector for sub-population 2
        information.vector.inv[((i-1) * 2 + 3)] = 4/((d.0.2[i] + d.1.2[i]))
      } 
      # Make information vector fit Tianchens code
      information.vector.inv.matrix = matrix(c(information.vector.inv[1], NA, information.vector.inv[2:length(information.vector.inv)]), nrow = K, byrow = TRUE)
      
      
      cov.mat.used = cov.mat.surv.restrict.yes.OneTreatmentArm(d.0.1, d.1.1, d.0.2, d.1.2)
      
    } # End restrict enrollment
  } # End if outcome = survival
  
  if(!is.null(relative.efficiency)){
    mean.vec = mean.vec * sqrt(relative.efficiency)
    information.vector.inv.matrix = information.vector.inv.matrix/relative.efficiency
  }
  
  # Need to calculate the inverse to get the information vector
  return(list(cov.mat.used=cov.mat.used,
              non.centrality.parameter.vec = mean.vec,
              information.vector = 1/information.vector.inv.matrix))
}



# This function calculates the efficacy boundaries
# Input: alpha.alloc = Vector of alpha allocations:
#        The alpha allocation vector is of length 2 * number of stages
#        the first two elements are the alpha allocations at stage one to each subpopulation
#        the next two elements are the alpha allocations at the next stage and so on
#        cov.mat.used: covariance matrix under the scenario
#        err.tol = how precise is the binary search

# Output:eff.bound: The vector z_{j,k} with blocks corresponding to stages and (z_{1,k}, z_{2,k}) within stage
#        eff.bound.alpha: K blocks where each block is (\tilde z_{1,k}, \tilde z_{2,k}) where block
#        k corresponds to alpha reallocated if both treatments in other sub-pop are rejected at stage k


get.eff.bound.OneTreatmentArm = function(alpha.allocation, cov.mat.used, err.tol = 10^-4, restrict.enrollment = FALSE){
  
  if(restrict.enrollment == FALSE){ 
    # Number of stages
    K = length(alpha.allocation)/2
    
    # Getting index corresponding to which sup-population is being used in each 
    # alpha allocation
    index.sub.pop = rep(c(1, 2), K)
    
    # eff.bound is the vector of efficacy boundaries with the elements corresponding to the
    # same stage and subpopulation combinations as in the alpha.allocation vector
    eff.bound = rep(NA, 2 * K)
    
    # Cumulative alpha allocation with subpopulation 1 and 2
    cum.alpha.1 = cumsum(alpha.allocation[which(index.sub.pop == 1)])
    cum.alpha.2 = cumsum(alpha.allocation[which(index.sub.pop == 2)])
    
    # Calculating the first elements of the efficacy boundary u_{j,1} corresponding to each 
    # subpopulation
    eff.bound[1] = qnorm(1-alpha.allocation[1])
    eff.bound[2] = qnorm(1-alpha.allocation[2])
    
    # A function that calculates the cumulative type one error corresponding 
    # to a sub-population j
    # effecacy boundaries is the efficacy boundary
    # cov.mat.used is the covariance matrix
    # j is the subpopulation
    sign.lev = function(eff.bound.used, cov.mat.used, sub.pop.numb){
      
      numb.stages = length(eff.bound.used)
      
      # Index which test-statistic belongs to population and treatment, respectivly
      index.sub.pop.eff = rep(c(1, 2), numb.stages)    
      
      cov.mat.eff.bound =
        cov.mat.used[1:(2* numb.stages), 1:(2* numb.stages)][which(index.sub.pop.eff == sub.pop.numb), which(index.sub.pop.eff == sub.pop.numb)]
      # Calculating the overall type one error under the global null
      type.1.err = 1- pmvnorm(mean = rep(0, numb.stages), sigma= cov.mat.eff.bound,lower = rep(-Inf, numb.stages), upper= eff.bound.used, algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))[1]
      
      return(type.1.err)
    }
    
    
    # Start calculating the efficacy boundaries associated with population 1
    # Cycling through the stages after and calculating the efficacy boundary z_{1,k} for k = 2, \ldots, K
    if(K > 1){
      for(i in 2:K){
        
        # Start by doing binary search for z_{1, k}
        # upper and lower values of interval
        upper = 20
        lower = -20
        length.int = upper - lower
        
        # Initial guess
        upper.bound.term = mean(c(upper, lower))
        eff.bound.j = eff.bound[index.sub.pop == 1]
        
        while(length.int > err.tol){
          eff.bound.j[i] = upper.bound.term
          alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 1)
          
          if(alpha.used < cum.alpha.1[i]){
            upper = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, lower))
          }
          if(alpha.used >= cum.alpha.1[i]){
            lower = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, upper))
          }
          length.int = upper - lower
        }
        
        # "Rounding" up to preserve type 1 error
        upper.bound.term = upper.bound.term + length.int
        eff.bound[(i-1)*2 + 1] = upper.bound.term
      }
      
      # Calculate the efficacy boundaries associated with population 2
      # Cycling through the stages and calculating the efficacy boundary 
      # z_{2,k} for k = 1, \ldots, K
      for(i in 2:K){
        # upper and lower values of interval
        upper = 20
        lower = -20
        length.int = upper - lower
        
        # Initial guess
        upper.bound.term = mean(c(upper, lower))
        eff.bound.j = eff.bound[index.sub.pop == 2]
        
        while(length.int > err.tol){
          eff.bound.j[i] = upper.bound.term
          alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 2)
          
          if(alpha.used < cum.alpha.2[i]){
            upper = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, lower))
          }
          if(alpha.used >= cum.alpha.2[i]){
            lower = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, upper))
          }
          length.int = upper - lower
        }
        
        # "Rounding" up to preserve type 1 error
        upper.bound.term = upper.bound.term + length.int
        eff.bound[(i-1)*2 + 2] = upper.bound.term
        
      } # End if K >1 statement
    }
    
    # efficacy boundaries associated with alphar reallocation
    # eff.bound.alpha[2*(k-1) +1] is \tilde z_{1,K} if both H_0 in subpopulation two 
    # are rejected at stage k
    # eff.bound.alpha[2*k] is \tilde u_{2,K} if both H_0 in subpopulation one 
    # are rejected at stage k
    
    eff.bound.alpha = rep(NA, 2 * K)
    
    for(k in 1:K){
      # Start with sub-population 1
      # Now we calculate the efficacy boundaries for pop 1 if null hypothesis corresponding to
      # pop 2 is rejected at stage k
      
      # Start a binary search for \tilde u_{1,K}
      # upper and lower values of interval
      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound[index.sub.pop == 1]
      # The cumulative alpha level that the last stage is allowed to test at (note not \alpha_{1,K})
      # \sum_{j=1}^K \alpha_{1,j} + \sum_{j=k}^K \alpha_{2,j}
      alpha.allowed = cum.alpha.1[K] + (cum.alpha.2[K] - c(0,cum.alpha.2)[k])
      
      while(length.int > err.tol){
        eff.bound.j[K] = upper.bound.term
        alpha.used = sign.lev(eff.bound.j, cov.mat.used, 1)
        
        if(alpha.used < alpha.allowed){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= alpha.allowed){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound.alpha[(k-1) * 2 + 1] = upper.bound.term
    }
    
    for(k in 1:K){
      # Now sub-population 2
      # Calculate the efficacy boundaries for pop 2 if both null hypothesis corresponding to
      # pop 1 are rejected at stage k
      
      # Start a binary search for \tilde u
      # upper and lower values of interval
      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound[index.sub.pop == 2]
      # The alpha level that the last stage is allowed to test at
      # \alpha_{1,K} + \sum_{j=k}^K \alpha_{2,j}
      alpha.allowed = cum.alpha.2[K] + (cum.alpha.1[K] - c(0,cum.alpha.1)[k])
      
      while(length.int > err.tol){
        eff.bound.j[K] = upper.bound.term
        alpha.used = sign.lev(eff.bound.j, cov.mat.used, 2)
        
        if(alpha.used < alpha.allowed){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= alpha.allowed){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound.alpha[k*2] = upper.bound.term
    } # end for k loop
  }
  
  if(restrict.enrollment == TRUE){
    
    # Number of stages
    K = (length(alpha.allocation) + 1)/2
    
    # Getting index corresponding to which sup-population is being used in each 
    # alpha allocation
    index.sub.pop = c(1,rep(c(1, 2), K-1))
    
    # eff.bound is the vector of efficacy boundaries with the elements corresponding to the
    # same stage and subpopulation combinations as in the alpha.allocation vector
    eff.bound = rep(NA, 2 * (K-1) + 1)
    
    # Cumulative alpha allocation with subpopulation 1 and 2
    cum.alpha.1 = cumsum(alpha.allocation[which(index.sub.pop == 1)])
    cum.alpha.2 = cumsum(alpha.allocation[which(index.sub.pop == 2)])
    
    # Calculating the first elements of the efficacy boundary u_{j,1} corresponding to each 
    # subpopulation
    eff.bound[1] = qnorm(1-alpha.allocation[1])
    eff.bound[3] = qnorm(1-alpha.allocation[3])
    
    # A function that calculates the cumulative type one error corresponding 
    # to a sub-population j
    # effecacy boundaries is the efficacy boundary
    # cov.mat.used is the covariance matrix
    # sub.pop.numb is the subpopulation
    sign.lev = function(eff.bound.used, cov.mat.used, sub.pop.numb){
      
      numb.stages = length(eff.bound.used)
      
      # Index which test-statistic belongs to population and treatment, respectivly
      index.sub.pop.eff = c(1, rep(c(1, 2), numb.stages -1))    
      
      cov.mat.eff.bound =
        cov.mat.used[1:(2* (numb.stages -1) +1), 1:(2* (numb.stages -1) +1)][which(index.sub.pop.eff == sub.pop.numb), which(index.sub.pop.eff == sub.pop.numb)]
      
      numb.stages = sum(index.sub.pop.eff == sub.pop.numb)
      # Calculating the overall type one error under the global null
      type.1.err = 1- pmvnorm(mean = rep(0, numb.stages), sigma= cov.mat.eff.bound,lower = rep(-Inf, numb.stages), upper= eff.bound.used, algorithm=GenzBretz(abseps = 0.000000001,maxpts=100000))[1]
      
      return(type.1.err)
    }
    
    
    # Start calculating the efficacy boundaries associated with population 1
    # Cycling through the stages after and calculating the efficacy boundary z_{1,k} for k = 2, \ldots, K
    if(K > 1){
      for(i in 2:K){
        
        # Start by doing binary search for z_{1, k}
        # upper and lower values of interval
        upper = 20
        lower = -20
        length.int = upper - lower
        
        # Initial guess
        upper.bound.term = mean(c(upper, lower))
        eff.bound.j = eff.bound[index.sub.pop == 1]
        
        while(length.int > err.tol){
          eff.bound.j[i] = upper.bound.term
          alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 1)
          
          if(alpha.used < cum.alpha.1[i]){
            upper = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, lower))
          }
          if(alpha.used >= cum.alpha.1[i]){
            lower = upper.bound.term
            upper.bound.term = mean(c(upper.bound.term, upper))
          }
          length.int = upper - lower
        }
        
        # "Rounding" up to preserve type 1 error
        upper.bound.term = upper.bound.term + length.int
        eff.bound[(i-2)*2 + 2] = upper.bound.term
      }
      
      # Calculate the efficacy boundaries associated with population 2
      # Cycling through the stages and calculating the efficacy boundary 
      # z_{2,k} for k = 1, \ldots, K
      if(K > 2){
        for(i in 2:(K-1)){
          # upper and lower values of interval
          upper = 20
          lower = -20
          length.int = upper - lower
          
          # Initial guess
          upper.bound.term = mean(c(upper, lower))
          eff.bound.j = eff.bound[index.sub.pop == 2]
          
          while(length.int > err.tol){
            eff.bound.j[i] = upper.bound.term
            alpha.used = sign.lev(eff.bound.j[1:i], cov.mat.used, 2)
            
            if(alpha.used < cum.alpha.2[i]){
              upper = upper.bound.term
              upper.bound.term = mean(c(upper.bound.term, lower))
            }
            if(alpha.used >= cum.alpha.2[i]){
              lower = upper.bound.term
              upper.bound.term = mean(c(upper.bound.term, upper))
            }
            length.int = upper - lower
          }
          
          # "Rounding" up to preserve type 1 error
          upper.bound.term = upper.bound.term + length.int
          eff.bound[(i-2)*2 + 5] = upper.bound.term
        } # End K > 2 statement
      } # End if K >1 statement
    }
    
    # efficacy boundaries associated with alphar reallocation
    # eff.bound.alpha[2*(k-1) +1] is \tilde z_{1,K} if both H_0 in subpopulation two 
    # are rejected 
    # eff.bound.alpha[2*k] is \tilde u_{2,K} if both H_0 in subpopulation one 
    # are rejected 
    
    eff.bound.alpha = rep(NA, 2)
    
    # Start with sub-population 1
    # Now we calculate the efficacy boundaries for pop 1 if poplation 1 is not 
    # strted or the null hypothesis corresponding to
    # pop 2 is rejected
    
    # Start a binary search for \tilde u_{1,K}
    # upper and lower values of interval
    upper = 20
    lower = -20
    length.int = upper - lower
    
    # Initial guess
    upper.bound.term = mean(c(upper, lower))
    eff.bound.j = eff.bound[index.sub.pop == 1]
    # The cumulative alpha level that the last stage is allowed to test at (note not \alpha_{1,K})
    # \sum_{j=1}^K \alpha_{1,j} + \sum_{j=k}^K \alpha_{2,j}
    alpha.allowed = cum.alpha.1[K] + cum.alpha.2[K-1]
    
    while(length.int > err.tol){
      eff.bound.j[K] = upper.bound.term
      alpha.used = sign.lev(eff.bound.j, cov.mat.used, 1)
      
      if(alpha.used < alpha.allowed){
        upper = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, lower))
      }
      if(alpha.used >= alpha.allowed){
        lower = upper.bound.term
        upper.bound.term = mean(c(upper.bound.term, upper))
      }
      length.int = upper - lower
    }
    
    # "Rounding" up to preserve type 1 error
    upper.bound.term = upper.bound.term + length.int
    eff.bound.alpha[1] = upper.bound.term
    
    # Now sub-population 2
    # Calculate the efficacy boundaries for pop 2 if both null hypothesis corresponding to
    # pop 1 are rejected
    
    
    if(K >2){
      # Start a binary search for \tilde u
      # upper and lower values of interval
      upper = 20
      lower = -20
      length.int = upper - lower
      
      # Initial guess
      upper.bound.term = mean(c(upper, lower))
      eff.bound.j = eff.bound[index.sub.pop == 2]
      # The alpha level that the last stage is allowed to test at
      # \alpha_{1,K} + \sum_{j=k}^K \alpha_{2,j}
      alpha.allowed = cum.alpha.2[K-1] + cum.alpha.1[K]
      
      while(length.int > err.tol){
        eff.bound.j[K-1] = upper.bound.term
        alpha.used = sign.lev(eff.bound.j, cov.mat.used, 2)
        
        if(alpha.used < alpha.allowed){
          upper = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, lower))
        }
        if(alpha.used >= alpha.allowed){
          lower = upper.bound.term
          upper.bound.term = mean(c(upper.bound.term, upper))
        }
        length.int = upper - lower
      }
      
      # "Rounding" up to preserve type 1 error
      upper.bound.term = upper.bound.term + length.int
      eff.bound.alpha[2] = upper.bound.term
    }
    # If only two stages
    if(K == 2){
      alpha.allowed = cum.alpha.2[K-1] + cum.alpha.1[K]
      eff.bound.alpha[2] = qnorm(1-alpha.allowed)
    }
  }  
  
  return(list(eff.bound = eff.bound, eff.bound.alpha = eff.bound.alpha))
}



# This function evaluates the performance of a given design
# Inputs: test.statistics: A matrix of test statistics where each row is a vector of test statistics
#         efficacy.boundary: all the different effacacy boundaries outputted from get.eff.bound
#         futility.boundary: A vector of futility boundaries with blocks corresponding to stages
#         and within a block the null hypothesis are ordered as in cov.mat.bin
#         alpha.allocation: The alpha allocation to each stage
#         restrict.enrollment: An indicator if enrollment at beginning of trial is restricted to 
#         sub-population 1
#         boundary.to.enroll: The boundary that the test statistic for sup-population 1
#         needs to cross at stage one to start enrolling from sub-population 2 at stage 1

# Output: A list consisting of three elements. The first one is a matrix of dim
#         n.sim times 2 where each column corresponds to if H_0 is rejected where the 
#         hypothesis are in the same order as in the covariance matrix. One means rejected 
#         and zero means not rejected.
#         The second element is of the same nature as the first element where
#         each column indicates at what stage the decision to reject or not reject
#         the corresponding hypothesis is made.
#         The third vector is the list of efficacy boundaries as outputted by get.eff.bound

# Connection to Betz et al. paper: "Comparison of Adaptive Randomized Trial Designs for Time-to-Event Outcomes that Expand Versus Restrict Enrollment Criteria, to Test Non-Inferiority Design Optimization"
#  The 3 design classes defined in Section 5.1 of that paper are:
#  1. D_{ONE-STAGE} [standard, single stage designs]
#  2. D_{ADAPTIVE,START-BOTH} [adaptive, starts by enrolling both subpops in stage 1, can restrict enrollment at interim analyses]
#  3. D_{ADAPTIVE,START-SUBPOP.1} [adaptive, starts by enrolling only subpop. 1; can initiate enrollment of subpop 2 based on decision at first interim analysis. At interim analyses after the first, the decisions are to continue enrollment or restrict.] 
#  To get the design class (1), use a single stage and set restrict.enrollment = FALSE.
#  To get the design class (2), set restrict.enrollment = FALSE.
#  To get the design class (3), set restrict.enrollment = TRUE.

design.evaluate.OneTreatmentArm <- function(test.statistics,
                            efficacy.boundary,
                            futility.boundary,
                            alpha.allocation,
                            restrict.enrollment = FALSE,
                            boundary.to.enroll = NULL){
  
  if(restrict.enrollment == FALSE){
    # Number of stages
    K = length(alpha.allocation)/2
    
    # Number of MC evaluations
    n.sim <- nrow(test.statistics)
    
    # z_{j,k} boundaries
    est.eff.bound= efficacy.boundary$eff.bound
    
    # The alpha re-allocated z_{j,k}
    est.eff.bound.alpha = efficacy.boundary$eff.bound.alpha
    
    # Going through each hypothesis being tested finding out if rejected or not
    # and when rejected/stopped
    
    reject.hyp = matrix(0, nrow = n.sim, ncol = 2)
    stage.decision = matrix(NA, nrow = n.sim, ncol = 2)
    
    for(i in 1:n.sim){
      # Start by looking at subpopulation 1
      index.sub.pop = rep(c(1, 2), K)
      index.stage = rep(1:K, each = 2)
      
      still.enrolled.1 = TRUE
      still.enrolled.2 = TRUE
      
      # Cycling through stages 
      for(k in 1:K){
        # Do the testing first in subpopulation 1
        if(still.enrolled.1 & test.statistics[i, 2* (k-1) + 1] > est.eff.bound[2* (k-1) + 1]){
          reject.hyp[i, 1] = 1
          stage.decision[i, 1] = k
          still.enrolled.1 = FALSE
        }
        
        # Doing futility stopping for sub-pop 1
        if(still.enrolled.1 & test.statistics[i, 2* (k-1) + 1] < futility.boundary[2* (k-1) + 1]){
          reject.hyp[i, 1] = 0
          stage.decision[i, 1] = k
          if(k != K){
            still.enrolled.1 = FALSE       
          }
        }
        
        # Do the testing in subpopulation 2
        if(still.enrolled.2 & test.statistics[i, 2* (k-1) + 2] > est.eff.bound[2* (k-1) + 2]){
          reject.hyp[i, 2] = 1
          stage.decision[i, 2] = k
          still.enrolled.2 = FALSE
        }
        
        # Doing futility stopping for sub-pop 2
        if(still.enrolled.2 & test.statistics[i, 2* (k-1) + 2] < futility.boundary[2* (k-1) + 2]){
          reject.hyp[i, 2] = 0
          if(k != K){
            stage.decision[i, 2] = k
            still.enrolled.2 = FALSE       
          }
        }
      }
      
      # All H_0 not already stopped for efficacy or futility for futility at stage K
      stage.decision[i, which(is.na(stage.decision[i, ]))] = K  
      
      
      # Alpha re-allocation
      # Cycle through stages
      for(k in 1:K){
        
        # Start with subpopulation 1
        # if subpop 2 is rejected at k
        if(stage.decision[i, 2] == k & reject.hyp[i, 2] == 1){
          
          # If still enrolled at stage K
          if(stage.decision[i, 1] == K & test.statistics[i, 2* (K-1) + 1] > est.eff.bound.alpha[1]){
            reject.hyp[i, 1] = 1
            stage.decision[i, 1] = K
          }
        } # end if loop if  subpop 2 is rejected
        
        # Now with subpopulation 2
        # if subpop 1 is rejected at stage k
        if(stage.decision[i, 1] == k & reject.hyp[i, 1] == 1){
          
          if(stage.decision[i, 2] == K & test.statistics[i, 2* (K-1) + 2] > est.eff.bound.alpha[2]){
            reject.hyp[i, 2] = 1
            stage.decision[i, 2] = K
          }
        } # end if loop if  subpop 2 is rejected
      } # end for k loop
      
    } # end for i loop
  }
  
  if(restrict.enrollment == TRUE){
    
    # Number of stages
    K = (length(alpha.allocation) + 1)/2
    
    # Number of MC evaluations
    n.sim <- nrow(test.statistics)
    
    # z_{j,k} boundaries
    est.eff.bound= efficacy.boundary$eff.bound
    
    # The alpha re-allocated z_{j,k}
    est.eff.bound.alpha = efficacy.boundary$eff.bound.alpha
    
    # Going through each hypothesis being tested finding out if rejected or not
    # and when rejected/stopped
    
    reject.hyp = matrix(0, nrow = n.sim, ncol = 2)
    stage.decision = matrix(NA, nrow = n.sim, ncol = 2)
    
    for(i in 1:n.sim){
      # Start by looking at subpopulation 1
      index.sub.pop = c(1, rep(c(1, 2), K-1))
      index.stage = c(1, rep(2:K, each = 2))
      
      still.enrolled.1 = TRUE
      still.enrolled.2 = TRUE
      
      # Cycling through stages 
      for(k in 1:K){
        # Do the testing first in subpopulation 1
        if(still.enrolled.1 & test.statistics[i, which(index.sub.pop ==1)][k] > est.eff.bound[which(index.sub.pop ==1)][k]){
          reject.hyp[i, 1] = 1
          stage.decision[i, 1] = k
          still.enrolled.1 = FALSE
        }
        
        # Doing futility stopping for sub-pop 1
        if(still.enrolled.1 & test.statistics[i, which(index.sub.pop ==1)][k] < futility.boundary[which(index.sub.pop ==1)][k]){
          reject.hyp[i, 1] = 0
          stage.decision[i, 1] = k
          if(k != K){
            still.enrolled.1 = FALSE       
          }
        }
        
        # See if enrollment does not start for population 2
        if(test.statistics[i, 1] <= boundary.to.enroll){
          still.enrolled.2 == FALSE
          stage.decision[i, 2] = 1
          reject.hyp[i, 2] = 0
        }
        
        # If enrollment starts for sub-population 2 at stage 1
        if(test.statistics[i, 1] > boundary.to.enroll){
          if(k >1){
            if(still.enrolled.2 & test.statistics[i, which(index.sub.pop ==2)][k-1] > est.eff.bound[which(index.sub.pop ==2)][k-1]){
              reject.hyp[i, 2] = 1
              stage.decision[i, 2] = k
              still.enrolled.2 = FALSE
            }
            
            # Doing futility stopping for sub-pop 2
            if(still.enrolled.2 & test.statistics[i, which(index.sub.pop ==2)][k-1] < futility.boundary[which(index.sub.pop ==2)][k-1]){
              reject.hyp[i, 2] = 0
              stage.decision[i, 2] = k
              if(k != K){
                still.enrolled.2 = FALSE       
              }
            }
          }  # End if k > 1 statment
        } # End if test.statistics[i, 1] > boundary.to.enroll statment
      }
      
      # All H_0 not already stopped for efficacy or futility for futility at stage K
      stage.decision[i, which(is.na(stage.decision[i, ]))] = K  
    } # End i loop
    
    # Alpha re-allocation
    # Cycle through stages
    for(i in 1:n.sim){
      # Start with sub-population 1
      if((reject.hyp[i, 2] == 1 & reject.hyp[i, 1] == 0)){
        # If still enrolled at stage K
        if(stage.decision[i, 1] == K & test.statistics[i, 2* (K-1)] > est.eff.bound.alpha[1]){
          reject.hyp[i, 1] = 1
        }
      } # end if loop if  subpop 2 is rejected
      
      # Now sub-population 2
      if(reject.hyp[i, 1] == 1 & reject.hyp[i, 2] == 0 & stage.decision[i, 2] != 0){
        # If still enrolled at stage K
        if(stage.decision[i, 2] == K & test.statistics[i, 2* (K-1) + 1] > est.eff.bound.alpha[2]){
          reject.hyp[i, 2] = 1
          stage.decision[i, 2] = K
        }
      } 
    } # end i loop 
  } # End if restrict.enrollment is TRUE
  
  colnames(reject.hyp) = c("A1", "A2")
  colnames(stage.decision) = c("A1", "A2")

  return(list(rejection.matrix = reject.hyp,
              stage.decision = stage.decision,
              eff.boundaries = efficacy.boundary))
}
