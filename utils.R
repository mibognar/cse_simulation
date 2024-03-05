library(MASS)
# Function to create a 3x3 covariance matrix
covmatrix = function(sd_a,sd_b,sd_c,cor_AB,cor_AC,cor_BC){
  covarAB = cor_AB*sd_a*sd_b
  covarAC = cor_AC*sd_a*sd_b
  covarBC = cor_BC*sd_b*sd_c
  
  cov_matrix = matrix(c(sd_a^2, covarAB, covarAC,
                        covarAB, sd_b^2, covarBC,
                        covarAC, covarBC,sd_c^2),
                      3,3)
  
  return(cov_matrix)
}
# function to use 3x3 covariance matrix with the MASS::mvnorm distribution
generate_condition = function(mean_a,mean_b,mean_c, cov_matrix){
  means = c(mean_a, mean_b, mean_c)
  data = MASS::mvrnorm(n = 1,  
                       mu = means,  
                       Sigma = cov_matrix
  )
  return(data)
}

# combine covariance matrix and mvnorm 
generate_subject_condition = function(mean_a,mean_b,mean_c,sd_a,sd_b,sd_c,cor_AB,cor_AC,cor_BC){
  covariance_matrix = covmatrix(sd_a,sd_b,sd_c,cor_AB,cor_AC,cor_BC)
  value_list = generate_condition(mean_a,mean_b,mean_c,covariance_matrix)
  return(value_list)
}