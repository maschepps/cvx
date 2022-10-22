N = 201
a = -6.91
b = 6.91
u = seq(from = a, to = b, length.out = N)
# Constraints
# Will need to change to complete rest of table 3.
e2 = 0.73
e3 = 0.76

#v0 = 60
v1 = 1.563
v2 = 0.825
v3 = 0.653
v4 = 0.137
p1 = 4

F1i = list()
for(j in 1:N){
  d = u[j]
  g = exp((v2 * d) + v3)
  f1 = matrix(c(1/(1+g), -v1*d*g/((1+g)^2), -v1*g/((1+g)^2), 1))
  F1i[[j]] = t(f1) %*% f1
  
}

cvx_func_1 = function(x){
  A1 = matrix(0, p1, p1)
  for(j in 1:N){
    A1 = A1 +F1i[[j]] * x[j]
  }
  obj = -det(A1)
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}
# Optimal solution
# x = 1/4 * u[c(1, 72, 108, 201)]

c1 = matrix(c(0, v3 / (v2 ^ 2), -1 / v2, 0), ncol = 1)
cvx_func_2 = function(x){
  A1 = matrix(0, p1, p1)
  for(j in 1:N){
    A1 = A1 +F1i[[j]] * x[j]
  }
  obj = t(c1) %*% A1 %*% c1
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}
# Optimal solution
# c(.2568, .1164, .3836, .2432) * u[1, 88, 89, 201]


c2 = matrix(c(-1/v1-1/v2, (v3+log(v1-1))/v2^2, -1 / v2, 0), ncol = 1)
cvx_func_3 = function(x){
  A1 = matrix(0, p1, p1)
  for(j in 1:N){
    A1 = A1 +F1i[[j]] * x[j]
  }
  obj = t(c2) %*% A1 %*% c2
  x = round(x, digits = 3)
  penalty = ifelse(sum(x) > 1, 
                   1000000000000000000000,
                   0)
  return(obj + penalty)
}
# Optimal solution
# c(.5, .0752, .4006, .0242) * u[1, 96, 97, 201]

#Multiobjective

