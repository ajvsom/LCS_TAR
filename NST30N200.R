# Timepoint and Sample Size Conditions: T = 30; N = 200 ------------------------------------------------

# Libraries
library(plm)
library (lavaan)
library(lcsm)

# Data Generating Process  ------------------------------------------------
# y_{it} = b_i + a*y{it-1} + e_it 
makedata = function(N, T, a = .75){
  count = 0
  y = c()
  Pmat <- matrix(nrow=(N*T), ncol=4)
  for(i in 1:N) {
    int = rnorm(1, 0, 1)
    for(j in 1:T) {
      count = count + 1
      if (j == 1){
        y[j] = int
      }
      else if (j > 1){
        y[j] = int + a*y[j-1] + rnorm(1, 0, .1)
      }
      Pmat[count,1] = i
      Pmat[count,2] = j
      Pmat[count,3] = y[j]
      Pmat[count,4] = int
    }
  }
  panel_df = data.frame(Pmat)
  names(panel_df)<-c("Unit", "Time", "Y", "Int")
  return(panel_df)
}


# Test Dataset ------------------------------------------------------------
dfn200u0t5_a = makedata(N = 200, T = 30, a = .5)
dfn200u0t5_a = pdata.frame(dfn200u0t5_a,index=c("Unit","Time")) 
dfn200u0t5_a$Y_lagged = lag(dfn200u0t5_a$Y, lag = 1)   
fixed <- plm(Y ~ Y_lagged, data=dfn200u0t5_a, 
             index=c("Unit", "Time"), model="within")
summary(fixed)

# Negative Nonstationary Trajectories ------------------------------------------------
# AR(1) = -1.5 ------------------------------------------------------------
# FE Model
K = 500
fe_n15 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n15 = makedata(N = 200, T = 30, a = -1.5)
  df_fe_n15 = pdata.frame(df_fe_n15,index=c("Unit","Time")) 
  df_fe_n15$Y_lagged = lag(df_fe_n15$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n15, 
               index=c("Unit", "Time"), model="within")
  fe_n15[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n15[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n15[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n15 = as.data.frame(fe_n15)
names(fe_n15) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n15 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n15 = makedata(N = 200, T = 30, a = -1.5)
  df_lcs_n15 = pdata.frame(df_lcs_n15,index=c("Unit","Time")) 
  df_lcs_n15_wide <- reshape(df_lcs_n15, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n15_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n15_lcsm_fit <- sem(df_lcs_n15_lcsm, data= df_lcs_n15_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n15_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n15[i,1] = 
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n15[i,2] =
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n15[i,3] =
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n15[i,4] =
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n15[i,5] =
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n15[i,6] =
      round(parameterEstimates(df_lcs_n15_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n15[i,7] =
      round(fitmeasures(df_lcs_n15_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n15[i,8] =
      round(fitmeasures(df_lcs_n15_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n15[i,9] =
      round(fitmeasures(df_lcs_n15_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n15 = as.data.frame(lcsm_n15)
names(lcsm_n15) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')


# AR(1) = -1.4 ------------------------------------------------------------
#FE Model
K = 500
fe_n14 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n14 = makedata(N = 200, T = 30, a = -1.4)
  df_fe_n14 = pdata.frame(df_fe_n14,index=c("Unit","Time")) 
  df_fe_n14$Y_lagged = lag(df_fe_n14$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n14, 
               index=c("Unit", "Time"), model="within")
  fe_n14[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n14[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n14[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n14 = as.data.frame(fe_n14)
names(fe_n14) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n14 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n14 = makedata(N = 200, T = 30, a = -1.4)
  df_lcs_n14 = pdata.frame(df_lcs_n14,index=c("Unit","Time")) 
  df_lcs_n14_wide <- reshape(df_lcs_n14, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n14_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n14_lcsm_fit <- sem(df_lcs_n14_lcsm, data= df_lcs_n14_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n14_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n14[i,1] = 
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n14[i,2] =
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n14[i,3] =
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n14[i,4] =
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n14[i,5] =
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n14[i,6] =
      round(parameterEstimates(df_lcs_n14_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n14[i,7] =
      round(fitmeasures(df_lcs_n14_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n14[i,8] =
      round(fitmeasures(df_lcs_n14_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n14[i,9] =
      round(fitmeasures(df_lcs_n14_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n14 = as.data.frame(lcsm_n14)
names(lcsm_n14) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -1.3 ------------------------------------------------------------
# FE Model 
K = 500
fe_n13 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n13 = makedata(N = 200, T = 30, a = -1.3)
  df_fe_n13 = pdata.frame(df_fe_n13,index=c("Unit","Time")) 
  df_fe_n13$Y_lagged = lag(df_fe_n13$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n13, 
               index=c("Unit", "Time"), model="within")
  fe_n13[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n13[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n13[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n13 = as.data.frame(fe_n13)
names(fe_n13) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n13 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n13 = makedata(N = 200, T = 30, a = -1.3)
  df_lcs_n13 = pdata.frame(df_lcs_n13,index=c("Unit","Time")) 
  df_lcs_n13_wide <- reshape(df_lcs_n13, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n13_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n13_lcsm_fit <- sem(df_lcs_n13_lcsm, data= df_lcs_n13_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n13_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n13[i,1] = 
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n13[i,2] =
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n13[i,3] =
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n13[i,4] =
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n13[i,5] =
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n13[i,6] =
      round(parameterEstimates(df_lcs_n13_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n13[i,7] =
      round(fitmeasures(df_lcs_n13_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n13[i,8] =
      round(fitmeasures(df_lcs_n13_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n13[i,9] =
      round(fitmeasures(df_lcs_n13_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n13 = as.data.frame(lcsm_n13)
names(lcsm_n13) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -1.2 ------------------------------------------------------------
# FE Model 
K = 500
fe_n12 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n12 = makedata(N = 200, T = 30, a = -1.2)
  df_fe_n12 = pdata.frame(df_fe_n12,index=c("Unit","Time")) 
  df_fe_n12$Y_lagged = lag(df_fe_n12$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n12, 
               index=c("Unit", "Time"), model="within")
  fe_n12[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n12[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n12[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n12 = as.data.frame(fe_n12)
names(fe_n12) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n12 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n12 = makedata(N = 200, T = 30, a = -1.2)
  df_lcs_n12 = pdata.frame(df_lcs_n12,index=c("Unit","Time")) 
  df_lcs_n12_wide <- reshape(df_lcs_n12, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n12_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n12_lcsm_fit <- sem(df_lcs_n12_lcsm, data= df_lcs_n12_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n12_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n12[i,1] = 
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n12[i,2] =
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n12[i,3] =
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n12[i,4] =
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n12[i,5] =
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n12[i,6] =
      round(parameterEstimates(df_lcs_n12_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n12[i,7] =
      round(fitmeasures(df_lcs_n12_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n12[i,8] =
      round(fitmeasures(df_lcs_n12_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n12[i,9] =
      round(fitmeasures(df_lcs_n12_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n12 = as.data.frame(lcsm_n12)
names(lcsm_n12) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')


# AR(1) = -1.1 ------------------------------------------------------------
# FE Model 
K = 500
fe_n11 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n11 = makedata(N = 200, T = 30, a = -1.1)
  df_fe_n11 = pdata.frame(df_fe_n11,index=c("Unit","Time")) 
  df_fe_n11$Y_lagged = lag(df_fe_n11$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n11, 
               index=c("Unit", "Time"), model="within")
  fe_n11[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n11[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n11[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n11 = as.data.frame(fe_n11)
names(fe_n11) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n11 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n11 = makedata(N = 200, T = 30, a = -1.1)
  df_lcs_n11 = pdata.frame(df_lcs_n11,index=c("Unit","Time")) 
  df_lcs_n11_wide <- reshape(df_lcs_n11, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n11_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n11_lcsm_fit <- sem(df_lcs_n11_lcsm, data= df_lcs_n11_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n11_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n11[i,1] = 
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n11[i,2] =
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n11[i,3] =
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n11[i,4] =
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n11[i,5] =
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n11[i,6] =
      round(parameterEstimates(df_lcs_n11_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n11[i,7] =
      round(fitmeasures(df_lcs_n11_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n11[i,8] =
      round(fitmeasures(df_lcs_n11_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n11[i,9] =
      round(fitmeasures(df_lcs_n11_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n11 = as.data.frame(lcsm_n11)
names(lcsm_n11) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -1.0 ------------------------------------------------------------
# FE Model 
K = 500
fe_n10 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n10 = makedata(N = 200, T = 30, a = -1.0)
  df_fe_n10 = pdata.frame(df_fe_n10,index=c("Unit","Time")) 
  df_fe_n10$Y_lagged = lag(df_fe_n10$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n10, 
               index=c("Unit", "Time"), model="within")
  fe_n10[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n10[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n10[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n10 = as.data.frame(fe_n10)
names(fe_n10) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n10 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n10 = makedata(N = 200, T = 30, a = -1.0)
  df_lcs_n10 = pdata.frame(df_lcs_n10,index=c("Unit","Time")) 
  df_lcs_n10_wide <- reshape(df_lcs_n10, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n10_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n10_lcsm_fit <- sem(df_lcs_n10_lcsm, data= df_lcs_n10_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n10_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n10[i,1] = 
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n10[i,2] =
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n10[i,3] =
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n10[i,4] =
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n10[i,5] =
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n10[i,6] =
      round(parameterEstimates(df_lcs_n10_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n10[i,7] =
      round(fitmeasures(df_lcs_n10_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n10[i,8] =
      round(fitmeasures(df_lcs_n10_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n10[i,9] =
      round(fitmeasures(df_lcs_n10_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n10 = as.data.frame(lcsm_n10)
names(lcsm_n10) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# Positive Nonstationary Processes ------------------------------------------------------------
# AR(1) = 1.0 ------------------------------------------------------------
# FE Model 
K = 500
fe_10 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_10 = makedata(N = 200, T = 30, a = 1.0)
  df_fe_10 = pdata.frame(df_fe_10,index=c("Unit","Time")) 
  df_fe_10$Y_lagged = lag(df_fe_10$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_10, 
               index=c("Unit", "Time"), model="within")
  fe_10[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_10[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_10[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_10 = as.data.frame(fe_10)
names(fe_10) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_10 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_10 = makedata(N = 200, T = 30, a = 1.0)
  df_lcs_10 = pdata.frame(df_lcs_10,index=c("Unit","Time")) 
  df_lcs_10_wide <- reshape(df_lcs_10, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_10_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_10_lcsm_fit <- sem(df_lcs_10_lcsm, data= df_lcs_10_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_10_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_10[i,1] = 
      round(parameterEstimates(df_lcs_10_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_10[i,2] =
      round(parameterEstimates(df_lcs_10_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_10[i,3] =
      round(parameterEstimates(df_lcs_10_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_10[i,4] =
      round(parameterEstimates(df_lcs_10_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_10[i,5] =
      round(parameterEstimates(df_lcs_10_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_10[i,6] =
      round(parameterEstimates(df_lcs_10_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_10[i,7] =
      round(fitmeasures(df_lcs_10_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_10[i,8] =
      round(fitmeasures(df_lcs_10_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_10[i,9] =
      round(fitmeasures(df_lcs_10_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_10 = as.data.frame(lcsm_10)
names(lcsm_10) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')



# AR(1) = 1.1 ------------------------------------------------------------
# FE Model 
K = 500
fe_11 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_11 = makedata(N = 200, T = 30, a = 1.1)
  df_fe_11 = pdata.frame(df_fe_11,index=c("Unit","Time")) 
  df_fe_11$Y_lagged = lag(df_fe_11$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_11, 
               index=c("Unit", "Time"), model="within")
  fe_11[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_11[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_11[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_11 = as.data.frame(fe_11)
names(fe_11) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_11 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_11 = makedata(N = 200, T = 30, a = 1.1)
  df_lcs_11 = pdata.frame(df_lcs_11,index=c("Unit","Time")) 
  df_lcs_11_wide <- reshape(df_lcs_11, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_11_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_11_lcsm_fit <- sem(df_lcs_11_lcsm, data= df_lcs_11_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_11_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_11[i,1] = 
      round(parameterEstimates(df_lcs_11_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_11[i,2] =
      round(parameterEstimates(df_lcs_11_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_11[i,3] =
      round(parameterEstimates(df_lcs_11_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_11[i,4] =
      round(parameterEstimates(df_lcs_11_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_11[i,5] =
      round(parameterEstimates(df_lcs_11_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_11[i,6] =
      round(parameterEstimates(df_lcs_11_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_11[i,7] =
      round(fitmeasures(df_lcs_11_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_11[i,8] =
      round(fitmeasures(df_lcs_11_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_11[i,9] =
      round(fitmeasures(df_lcs_11_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_11 = as.data.frame(lcsm_11)
names(lcsm_11) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 1.2 ------------------------------------------------------------
# FE Model 
K = 500
fe_12 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_12 = makedata(N = 200, T = 30, a = 1.2)
  df_fe_12 = pdata.frame(df_fe_12,index=c("Unit","Time")) 
  df_fe_12$Y_lagged = lag(df_fe_12$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_12, 
               index=c("Unit", "Time"), model="within")
  fe_12[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_12[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_12[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_12 = as.data.frame(fe_12)
names(fe_12) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_12 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_12 = makedata(N = 200, T = 30, a = 1.2)
  df_lcs_12 = pdata.frame(df_lcs_12,index=c("Unit","Time")) 
  df_lcs_12_wide <- reshape(df_lcs_12, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_12_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_12_lcsm_fit <- sem(df_lcs_12_lcsm, data= df_lcs_12_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_12_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_12[i,1] = 
      round(parameterEstimates(df_lcs_12_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_12[i,2] =
      round(parameterEstimates(df_lcs_12_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_12[i,3] =
      round(parameterEstimates(df_lcs_12_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_12[i,4] =
      round(parameterEstimates(df_lcs_12_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_12[i,5] =
      round(parameterEstimates(df_lcs_12_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_12[i,6] =
      round(parameterEstimates(df_lcs_12_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_12[i,7] =
      round(fitmeasures(df_lcs_12_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_12[i,8] =
      round(fitmeasures(df_lcs_12_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_12[i,9] =
      round(fitmeasures(df_lcs_12_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_12 = as.data.frame(lcsm_12)
names(lcsm_12) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 1.3 ------------------------------------------------------------
# FE Model 
K = 500
fe_13 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_13 = makedata(N = 200, T = 30, a = 1.3)
  df_fe_13 = pdata.frame(df_fe_13,index=c("Unit","Time")) 
  df_fe_13$Y_lagged = lag(df_fe_13$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_13, 
               index=c("Unit", "Time"), model="within")
  fe_13[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_13[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_13[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_13 = as.data.frame(fe_13)
names(fe_13) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_13 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_13 = makedata(N = 200, T = 30, a = 1.3)
  df_lcs_13 = pdata.frame(df_lcs_13,index=c("Unit","Time")) 
  df_lcs_13_wide <- reshape(df_lcs_13, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_13_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_13_lcsm_fit <- sem(df_lcs_13_lcsm, data= df_lcs_13_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_13_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_13[i,1] = 
      round(parameterEstimates(df_lcs_13_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_13[i,2] =
      round(parameterEstimates(df_lcs_13_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_13[i,3] =
      round(parameterEstimates(df_lcs_13_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_13[i,4] =
      round(parameterEstimates(df_lcs_13_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_13[i,5] =
      round(parameterEstimates(df_lcs_13_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_13[i,6] =
      round(parameterEstimates(df_lcs_13_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_13[i,7] =
      round(fitmeasures(df_lcs_13_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_13[i,8] =
      round(fitmeasures(df_lcs_13_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_13[i,9] =
      round(fitmeasures(df_lcs_13_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_13 = as.data.frame(lcsm_13)
names(lcsm_13) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 1.4 ------------------------------------------------------------
#FE Model
K = 500
fe_14 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_14 = makedata(N = 200, T = 30, a = 1.4)
  df_fe_14 = pdata.frame(df_fe_14,index=c("Unit","Time")) 
  df_fe_14$Y_lagged = lag(df_fe_14$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_14, 
               index=c("Unit", "Time"), model="within")
  fe_14[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_14[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_14[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_14 = as.data.frame(fe_14)
names(fe_14) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_14 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_14 = makedata(N = 200, T = 30, a = 1.4)
  df_lcs_14 = pdata.frame(df_lcs_14,index=c("Unit","Time")) 
  df_lcs_14_wide <- reshape(df_lcs_14, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_14_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_14_lcsm_fit <- sem(df_lcs_14_lcsm, data= df_lcs_14_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_14_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_14[i,1] = 
      round(parameterEstimates(df_lcs_14_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_14[i,2] =
      round(parameterEstimates(df_lcs_14_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_14[i,3] =
      round(parameterEstimates(df_lcs_14_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_14[i,4] =
      round(parameterEstimates(df_lcs_14_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_14[i,5] =
      round(parameterEstimates(df_lcs_14_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_14[i,6] =
      round(parameterEstimates(df_lcs_14_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_14[i,7] =
      round(fitmeasures(df_lcs_14_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_14[i,8] =
      round(fitmeasures(df_lcs_14_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_14[i,9] =
      round(fitmeasures(df_lcs_14_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_14 = as.data.frame(lcsm_14)
names(lcsm_14) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 1.5 ------------------------------------------------------------
# FE Model
K = 500
fe_15 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_15 = makedata(N = 200, T = 30, a = 1.5)
  df_fe_15 = pdata.frame(df_fe_15,index=c("Unit","Time")) 
  df_fe_15$Y_lagged = lag(df_fe_15$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_15, 
               index=c("Unit", "Time"), model="within")
  fe_15[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_15[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_15[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_15 = as.data.frame(fe_15)
names(fe_15) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_15 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_15 = makedata(N = 200, T = 30, a = 1.5)
  df_lcs_15 = pdata.frame(df_lcs_15,index=c("Unit","Time")) 
  df_lcs_15_wide <- reshape(df_lcs_15, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_15_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_15_lcsm_fit <- sem(df_lcs_15_lcsm, data= df_lcs_15_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_15_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_15[i,1] = 
      round(parameterEstimates(df_lcs_15_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_15[i,2] =
      round(parameterEstimates(df_lcs_15_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_15[i,3] =
      round(parameterEstimates(df_lcs_15_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_15[i,4] =
      round(parameterEstimates(df_lcs_15_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_15[i,5] =
      round(parameterEstimates(df_lcs_15_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_15[i,6] =
      round(parameterEstimates(df_lcs_15_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_15[i,7] =
      round(fitmeasures(df_lcs_15_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_15[i,8] =
      round(fitmeasures(df_lcs_15_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_15[i,9] =
      round(fitmeasures(df_lcs_15_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_15 = as.data.frame(lcsm_15)
names(lcsm_15) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')


# Calculating Bias for FE -------------------------------------------------
#Combine all FE datasets
all_fe_data = rbind(fe_n15, fe_n14, fe_n13, 
                    fe_n12, fe_n11, fe_n10,
                    fe_15, fe_14, fe_13, 
                    fe_12, fe_11, fe_10)
cond = c(rep(-1.5, 500), rep(-1.4, 500), rep(-1.3, 500),
         rep(-1.2, 500), rep(-1.1, 500), rep(-1.0, 500),
         rep(1.5, 500), rep(1.4, 500), rep(1.3, 500),
         rep(1.2, 500), rep(1.1, 500), rep(1.0, 500))
all_fe_outputs = cbind(cond, all_fe_data)

#Compute average FE Estimate
Ave_fe_est = rbind(mean(subset(all_fe_outputs, cond == -1.5)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -1.4)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -1.3)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -1.2)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -1.1)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -1.0)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.0)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.1)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.2)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.3)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.4)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 1.5)$`Est(b)`))

#Compute bias in average FE Estimate
Bias_fe_est = rbind(mean(abs(subset(all_fe_outputs, cond == -1.5)$`Est(b)` - (-1.5))),
                    mean(abs(subset(all_fe_outputs, cond == -1.4)$`Est(b)`  - (-1.4))),
                    mean(abs(subset(all_fe_outputs, cond == -1.3)$`Est(b)`  - (-1.3))),
                    mean(abs(subset(all_fe_outputs, cond == -1.2)$`Est(b)`  - (-1.2))),
                    mean(abs(subset(all_fe_outputs, cond == -1.1)$`Est(b)`  - (-1.1))),
                    mean(abs(subset(all_fe_outputs, cond == -1.0)$`Est(b)`  - (-1.0))),
                    mean(abs(subset(all_fe_outputs, cond == 1.0)$`Est(b)`  - (1.0))),
                    mean(abs(subset(all_fe_outputs, cond == 1.1)$`Est(b)`  - (1.1))),
                    mean(abs(subset(all_fe_outputs, cond == 1.2)$`Est(b)`  - (1.2))),
                    mean(abs(subset(all_fe_outputs, cond == 1.3)$`Est(b)`  - (1.3))),
                    mean(abs(subset(all_fe_outputs, cond == 1.4)$`Est(b)`  - (1.4))),
                    mean(abs(subset(all_fe_outputs, cond == 1.5)$`Est(b)`  - (1.5))))

# RMSE_fe_est = rbind(mean((subset(all_fe_outputs, cond == -1.5)$`Est(b)`-(-1.5))^2),
#                     mean((subset(all_fe_outputs, cond == -1.0)$`Est(b)`-(-1))^2),
#                     mean((subset(all_fe_outputs, cond == -0.5)$`Est(b)`-(-0.5))^2),
#                     mean((subset(all_fe_outputs, cond == 0)$`Est(b)`-0)^2),
#                     mean((subset(all_fe_outputs, cond == 0.5)$`Est(b)`-0.5)^2),
#                     mean((subset(all_fe_outputs, cond == 1.0)$`Est(b)`-1.0)^2),
#                     mean((subset(all_fe_outputs, cond == 1.5)$`Est(b)`-1.5))^2)


#Compute average FE Standard Errors
Ave_fe_se = rbind(mean(subset(all_fe_outputs, cond == -1.5)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -1.4)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -1.3)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -1.2)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -1.1)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -1.0)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.0)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.1)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.2)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.3)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.4)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 1.5)$`SE(b)`))



#Calculate Type 2 errors - what percentage were found not significant
for (i in 1:length(all_fe_outputs$cond)){
  if (all_fe_outputs$`pval(b)`[i]>.05){
    all_fe_outputs$Type_2[i] = 1
  }
  else {
    all_fe_outputs$Type_2[i] = 0
  }
}

T2_fe_percent = rbind(mean(subset(all_fe_outputs, cond == -1.5)$Type_2),
                      mean(subset(all_fe_outputs, cond == -1.4)$Type_2),
                      mean(subset(all_fe_outputs, cond == -1.3)$Type_2),
                      mean(subset(all_fe_outputs, cond == -1.2)$Type_2),
                      mean(subset(all_fe_outputs, cond == -1.1)$Type_2),
                      mean(subset(all_fe_outputs, cond == -1.0)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.0)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.1)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.2)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.3)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.4)$Type_2),
                      mean(subset(all_fe_outputs, cond == 1.5)$Type_2))

#Calculate the percentage of estimates that fall within the 95% confidence interval
#Create Lower CI
all_fe_outputs$CIL = all_fe_outputs$`Est(b)`-1.96*all_fe_outputs$`SE(b)`
#Create Upper CI
all_fe_outputs$CIU = all_fe_outputs$`Est(b)`+1.96*all_fe_outputs$`SE(b)`
#If statement 
all_fe_outputs$CI = ifelse(all_fe_outputs$cond < all_fe_outputs$CIU & 
                             all_fe_outputs$cond > all_fe_outputs$CIL, 1, 0)

CI_fe_percent = rbind(mean(subset(all_fe_outputs, cond == -1.5)$CI),
                      mean(subset(all_fe_outputs, cond == -1.4)$CI),
                      mean(subset(all_fe_outputs, cond == -1.3)$CI),
                      mean(subset(all_fe_outputs, cond == -1.2)$CI),
                      mean(subset(all_fe_outputs, cond == -1.1)$CI),
                      mean(subset(all_fe_outputs, cond == -1.0)$CI),
                      mean(subset(all_fe_outputs, cond == 1.0)$CI),
                      mean(subset(all_fe_outputs, cond == 1.1)$CI),
                      mean(subset(all_fe_outputs, cond == 1.2)$CI),
                      mean(subset(all_fe_outputs, cond == 1.3)$CI),
                      mean(subset(all_fe_outputs, cond == 1.4)$CI),
                      mean(subset(all_fe_outputs, cond == 1.5)$CI))




a_coef = c(seq(-1.5, -1.0, .1), seq(1.0, 1.5, .1))
out_fe = as.data.frame(cbind(a_coef, Ave_fe_est, Bias_fe_est, Ave_fe_se, T2_fe_percent, CI_fe_percent))
names(out_fe) = c('a', 'Ave_fe_est', 'Bias_fe_est', 'Ave_fe_se', 'error_rate', 'CI_error_rate') 

# Calculating Bias for LCSM -------------------------------------------------
#Combine LCSM datasets
all_lcsm_data = rbind(lcsm_n15, lcsm_n14, lcsm_n13, 
                      lcsm_n12, lcsm_n11, lcsm_n10,
                      lcsm_15, lcsm_14, lcsm_13, 
                      lcsm_12, lcsm_11, lcsm_10)
lcsm_cond = c(rep(-2.5, 500), rep(-2.4, 500), rep(-2.3, 500),
              rep(-2.2, 500), rep(-2.1, 500), rep(-2.0, 500),
              rep(.5, 500), rep(.4, 500), rep(.3, 500),
              rep(.2, 500), rep(.1, 500), rep(0, 500))
all_lcsm_outputs = cbind(lcsm_cond, all_lcsm_data)

#Compute average LCSM Estimate
Ave_lcsm_est = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == 0)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == .1)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == .2)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == .3)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == .4)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == .5)$`Est(b)`, na.rm = T))
#Compute Bias in LCSM Estimate
Bias_lcsm_est = rbind(mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.5)$`Est(b)`- (-2.5)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.4)$`Est(b)` - (-2.4)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.3)$`Est(b)` - (-2.3)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.2)$`Est(b)` - (-2.2)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.1)$`Est(b)` - (-2.1)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == -2.0)$`Est(b)` - (-2.0)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == 0)$`Est(b)` - (0)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == .1)$`Est(b)` - (.1)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == .2)$`Est(b)` - (.2)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == .3)$`Est(b)` - (.3)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == .4)$`Est(b)` - (.4)), na.rm = T),
                      mean(abs(subset(all_lcsm_outputs, lcsm_cond == .5)$`Est(b)` - (.5)), na.rm = T))

#Compute average LCSM Standard error
Ave_lcsm_se = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == 0)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == .1)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == .2)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == .3)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == .4)$`SE(b)`, na.rm = T),
                    mean(subset(all_lcsm_outputs, lcsm_cond == .5)$`SE(b)`, na.rm = T))


# RMSE_lcsm_est = rbind(mean((subset(all_lcsm_outputs, lcsm_cond == -2.5)$`Est(b)`-(-2.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -2.0)$`Est(b)`-(-2.0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -1.5)$`Est(b)`-(-1.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -1.0)$`Est(b)`-(-1.0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -0.5)$`Est(b)`-(-0.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == 0)$`Est(b)`-(0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == 0.5)$`Est(b)`-(.5))^2, na.rm = T))

#Compute average fit index
Ave_CFI = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == 0)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == .1)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == .2)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == .3)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == .4)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == .5)$cfi, na.rm = T))

Ave_RMSEA = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == 0)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == .1)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == .2)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == .3)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == .4)$rmsea, na.rm = T),
                  mean(subset(all_lcsm_outputs, lcsm_cond == .5)$rmsea, na.rm = T))

Ave_SRMR = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == 0)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == .1)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == .2)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == .3)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == .4)$srmr, na.rm = T),
                 mean(subset(all_lcsm_outputs, lcsm_cond == .5)$srmr, na.rm = T))

#Compute percent significant
all_lcsm_outputs[is.na(all_lcsm_outputs)] = -99

for (i in 1:length(all_lcsm_outputs$lcsm_cond)){
  if (all_lcsm_outputs$`pval(b)`[i]>.05){
    all_lcsm_outputs$Type_2[i] = 1
  }
  else if (all_lcsm_outputs$`pval(b)`[i]==-99){
    all_lcsm_outputs$Type_2[i] = NA
  }
  else {
    all_lcsm_outputs$Type_2[i] = 0
  }
}

all_lcsm_outputs[all_lcsm_outputs == -99] = NA

T2_lcsm_percent = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == 0)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .1)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .2)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .3)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .4)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .5)$Type_2, na.rm = T))

#Calculate the percentage of estimates that fall within the 95% confidence interval
#Create Lower CI
all_lcsm_outputs$CIL = all_lcsm_outputs$`Est(b)`-1.96*all_lcsm_outputs$`SE(b)`
#Create Upper CI
all_lcsm_outputs$CIU = all_lcsm_outputs$`Est(b)`+1.96*all_lcsm_outputs$`SE(b)`
#If statement 
all_lcsm_outputs$CI = ifelse(is.na(all_lcsm_outputs$`Est(b)`), 0, 
                             ifelse(all_lcsm_outputs$lcsm_cond < all_lcsm_outputs$CIU & 
                                      all_lcsm_outputs$lcsm_cond > all_lcsm_outputs$CIL, 1, 0))

CI_lcsm_percent = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -2.5)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.4)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.3)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.2)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.1)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -2.0)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == 0)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .1)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .2)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .3)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .4)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == .5)$CI, na.rm = T))


p_coef = c(seq(-2.5, -2.0, .1), seq(0, 0.5, .1))
out_lcsm = as.data.frame(cbind(p_coef, Ave_lcsm_est, Bias_lcsm_est, Ave_lcsm_se, T2_lcsm_percent, CI_lcsm_percent,
                               Ave_CFI, Ave_RMSEA, Ave_SRMR))
names(out_lcsm) = c('p', 'Ave_est', 'Bias_est', 'Ave_se', 'error_rate', 'CI_error_rate', 'Ave_CFI', 'Ave_RMSEA', 'Ave_SRMR')                

#final table

bias_table_out = cbind(out_fe, out_lcsm)

