#Optimize 
library(tictoc)


#Metaheuristicopt
# library(metaheuristicOpt)
# algos = c('PSO',
#           'ALO',
#           'GWO',
#           'DA',
#           'FFA',
#           'GA',
#           'GOA',
#           'HS',
#           'MFO',
#           'SCA',
#           'WOA',
#           'CLONALG',
#           'DE',
#           'SFL',
#           # 'CSO',
#           'ABC',
#           # 'KH',
#           'CS',
#           'BA',
#           'GBS',
#           'BHO')
# algos = c('PSO')
# rangeV = rbind(rep(0, N),
#                rep(1, N))
# ans_a =   metaOpt(FUN = cvx_func,
#                   optimType = "MIN",
#                   algorithm = "PSO",
#                   numVar    = N,
#                   rangeVar = rangeV,
#                   control = list(maxIter = 1000))

tic()
indo =   metaOpt(FUN = cvx_func_1,
                 optimType = "MIN",
                 algorithm = "PSO",
                 numVar    = N,
                 rangeVar = rangeV,
                 control = list(maxIter = 38))
toc()

#EmiR
# library(EmiR)
# params = rbind(rep(0, N),
#                rep(1, N))
# params = parameters(params)
# conf_algo <- config_algo(algorithm_id = 'BAT', 
#                          iterations = 10000, 
#                          population_size = 100)
# ans_b = minimize(algorithm_id = list("BAT"), 
#                  obj_func = cvx_func, 
#                  parameters = params,
#                  config = conf_algo,
#                  save_pop_history = T)

# library(ppso)
# a = optim_pso(objective_function = cvx_func,
#           number_of_parameters = N,
#           initial_estimates = data.frame((rangeV[1,])),
#           parameter_bounds = t(rangeV),
#           max_number_of_iterations = 1000,
#           max_number_function_calls = 1000
#           )

