
#racing algorithm setup
library(irace)
target_runner <- function(experiment, scenario){
  library(ppso)
  N = 18
  a = 0
  b = 500
  u = seq(from = a, to = b, length.out = N)
  v0 = 60
  v1 = 294
  v2 = 25
  v02 = 60
  v12 = 340
  v22 = 107.14
  v03 = 49.62
  v13 = 290.51
  v23 = 150
  v33 = 45.51
  
  p1 = 2
  p2 = 3
  p3 = p2
  p4 = 4
  
  F1i = list()
  F2i = list()
  F3i = list()
  F4i = list()
  
  for(j in 1:N){
    d = u[j]
    f1 = matrix(c(1, d), nrow = 1)
    f2 = matrix(c(1, d / (v2 + d) - v1 * d / (v2 + d)^2), nrow = 1)
    f3 = matrix(c(1, d / (v22 + d) - v12 * d / (v22 +d)^2), nrow = 1)
    g = exp((v23 - d) / v33)
    f4 = matrix(c(1, 
                  1 / (1 + g), -v13 * g / v33 / ((1 + g)^2),
                  v13 * g * (v23 - d) / (v33^2) / ((1 + g)^2)))
    F1i[[j]] = t(f1) %*% f1
    F2i[[j]] = t(f2) %*% f2
    F3i[[j]] = t(f3) %*% f3
    F4i[[j]] = t(f4) %*% f4
  }
  
  
  library(metaheuristicOpt)
  algos = c('PSO',
            'ALO',
            'GWO',
            'DA',
            'FFA',
            'GA',
            'GOA',
            'HS',
            'MFO',
            'SCA',
            'WOA',
            'CLONALG',
            'DE',
            'SFL',
            # 'CSO',
            'ABC',
            # 'KH',
            'CS',
            'BA',
            'GBS',
            'BHO')
  rangeV = rbind(rep(0, N),
                 rep(1, N))
  
  x = c(.5, rep(0, N-2), .5)
  cvx_func = function(x){
    A1 = matrix(0, p1, p1)
    for(j in 1:N){
      A1 = A1 +F1i[[j]] * x[j]
    }
    obj = -det(A1)
    penalty = ifelse(sum(x) > 1, 
                     1000000000000000000000,
                     0)
    return(obj + penalty)
  }
  
  start = proc.time()
  configuration = experiment$configuration
  
  a = optim_pso(objective_function = cvx_func,
                number_of_parameters = N,
                # initial_estimates = data.frame((rangeV[1,])),
                parameter_bounds = t(rangeV),
                max_number_of_iterations = 10000,
                max_number_function_calls = 10000,
                w = configuration[['w']],
                C1 = configuration[['Ca']],
                C2 = configuration[['Cb']]
  )
  time1 = proc.time() - start
  return(list(cost = a$value,
              time = time1[1]))
}
# time1 = proc.time() - start
# unregister_dopar <- function() {
#   env <- foreach:::.foreachGlobals
#   rm(list=ls(name=env), pos=env)
# }
# unregister_dopar()
# return(list(#cost = a$value,
#   cost = results@best_cost,
#   # cost = a$optimumValue,
#   time = time1[1]))
library(tictoc)
tic()
scenario <- list(targetRunner = target_runner,
                 instances = NULL,
                 maxExperiments = 1000,
                 # parallel = 20,
                 # mpi = 1,
                 # Do not create a logFile
                 logFile = "amgen_pso_race.Rdata"
)
parameters_table <- '
 w "" r (1, 20000)
 Ca "" r (1, 20000)
 Cb "" r (1, 20000)
 '

## We use the irace function readParameters to read this table:
parameters <- readParameters(text = parameters_table)
setwd('C:/Users/Admin/Desktop/Github/tuning')
checkIraceScenario(scenario, parameters = parameters)
ss = irace(scenario = scenario, parameters = parameters)
toc()