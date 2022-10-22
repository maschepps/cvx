# Journal of Computational and Graphical Statistics
# Using CVX to Construct Optimal Designs for
# Biomedical Studies with Multiple Objectives
# Weng Kee Wong & Julie Zhou
# Department of Biostatistics, University of California, Los Angeles, CA;
# Department of Mathematics and Statistics, University of Victoria, Victoria, BC,Canada

# Table 5 
# D-Optimal answer
# Right-most column
# Can alter for both other columns.

N = 201
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
  f2 = matrix(c(1, d / (v2 + d), -v1 * d / (v2 + d)^2), nrow = 1)
  f3 = matrix(c(1, d / (v22 + d), -v12 * d / (v22 +d)^2), nrow = 1)
  g = exp((v23 - d) / v33)
  f4 = matrix(c(1, 
                1 / (1 + g), -v13 * g / v33 / ((1 + g)^2),
                v13 * g * (v23 - d) / (v33^2) / ((1 + g)^2)))
  F1i[[j]] = t(f1) %*% f1
  F2i[[j]] = t(f2) %*% f2
  F3i[[j]] = t(f3) %*% f3
  F4i[[j]] = t(f4) %*% f4
}

# Create objective functions
cvx_func_1 = function(x){
  x = round(x, digits = 3)
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
# Optimal answer for cvx_func_1 
# x = c(.5, rep(0, N-2), .5)


cvx_func_2 = function(x){
  A2 = matrix(0, p2, p2)
  for(j in 1:N){
    A2 = A2 +F2i[[j]] * x[j]
  }
  obj = -det(A2)
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}
# So called Optimal Answer for cvx_func_2
# x = c(1/3, rep(0, 8), 1/3, rep(0, 201 - 10), 1/3)

cvx_func_3 = function(x){
  A3 = matrix(0, p3, p3)
  for(j in 1:N){
    A3 = A3 +F3i[[j]] * x[j]
  }
  obj = -det(A3)
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}

cvx_func_4 = function(x){
  A4 = matrix(0, p4, p4)
  for(j in 1:N){
    A4 = A4 +F4i[[j]] * x[j]
  }
  obj = -det(A4)
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}



#multiobjective optimization.

################
##################
###################
##################
#End of problem 1.




