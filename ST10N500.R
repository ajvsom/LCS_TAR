# Timepoint and Sample Size Conditions: T = 10; N = 500 ------------------------------------------------

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
dfn200u0t5_a = makedata(N = 500, T = 10, a = .5)
dfn200u0t5_a = pdata.frame(dfn200u0t5_a,index=c("Unit","Time")) 
dfn200u0t5_a$Y_lagged = lag(dfn200u0t5_a$Y, lag = 1)   
fixed <- plm(Y ~ Y_lagged, data=dfn200u0t5_a, 
             index=c("Unit", "Time"), model="within")
summary(fixed)

# Negative Stationary Trajectories ------------------------------------------------
# AR(1) = -0.9 ------------------------------------------------------------
# FE Model
K = 500
fe_n09 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n09 = makedata(N = 500, T = 10, a = -0.9)
  df_fe_n09 = pdata.frame(df_fe_n09,index=c("Unit","Time")) 
  df_fe_n09$Y_lagged = lag(df_fe_n09$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n09, 
               index=c("Unit", "Time"), model="within")
  fe_n09[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n09[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n09[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n09 = as.data.frame(fe_n09)
names(fe_n09) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n09 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n09 = makedata(N = 500, T = 10, a = -0.9)
  df_lcs_n09 = pdata.frame(df_lcs_n09,index=c("Unit","Time")) 
  df_lcs_n09_wide <- reshape(df_lcs_n09, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n09_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n09_lcsm_fit <- sem(df_lcs_n09_lcsm, data= df_lcs_n09_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n09_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n09[i,1] = 
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n09[i,2] =
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n09[i,3] =
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n09[i,4] =
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n09[i,5] =
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n09[i,6] =
      round(parameterEstimates(df_lcs_n09_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n09[i,7] =
      round(fitmeasures(df_lcs_n09_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n09[i,8] =
      round(fitmeasures(df_lcs_n09_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n09[i,9] =
      round(fitmeasures(df_lcs_n09_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n09 = as.data.frame(lcsm_n09)
names(lcsm_n09) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.8 ------------------------------------------------------------
# FE Model
K = 500
fe_n08 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n08 = makedata(N = 500, T = 10, a = -0.8)
  df_fe_n08 = pdata.frame(df_fe_n08,index=c("Unit","Time")) 
  df_fe_n08$Y_lagged = lag(df_fe_n08$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n08, 
               index=c("Unit", "Time"), model="within")
  fe_n08[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n08[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n08[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n08 = as.data.frame(fe_n08)
names(fe_n08) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n08 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n08 = makedata(N = 500, T = 10, a = -0.8)
  df_lcs_n08 = pdata.frame(df_lcs_n08,index=c("Unit","Time")) 
  df_lcs_n08_wide <- reshape(df_lcs_n08, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n08_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n08_lcsm_fit <- sem(df_lcs_n08_lcsm, data= df_lcs_n08_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n08_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n08[i,1] = 
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n08[i,2] =
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n08[i,3] =
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n08[i,4] =
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n08[i,5] =
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n08[i,6] =
      round(parameterEstimates(df_lcs_n08_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n08[i,7] =
      round(fitmeasures(df_lcs_n08_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n08[i,8] =
      round(fitmeasures(df_lcs_n08_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n08[i,9] =
      round(fitmeasures(df_lcs_n08_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n08 = as.data.frame(lcsm_n08)
names(lcsm_n08) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.7 ------------------------------------------------------------
# FE Model
K = 500
fe_n07 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n07 = makedata(N = 500, T = 10, a = -0.7)
  df_fe_n07 = pdata.frame(df_fe_n07,index=c("Unit","Time")) 
  df_fe_n07$Y_lagged = lag(df_fe_n07$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n07, 
               index=c("Unit", "Time"), model="within")
  fe_n07[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n07[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n07[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n07 = as.data.frame(fe_n07)
names(fe_n07) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n07 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n07 = makedata(N = 500, T = 10, a = -0.7)
  df_lcs_n07 = pdata.frame(df_lcs_n07,index=c("Unit","Time")) 
  df_lcs_n07_wide <- reshape(df_lcs_n07, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n07_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n07_lcsm_fit <- sem(df_lcs_n07_lcsm, data= df_lcs_n07_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n07_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n07[i,1] = 
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n07[i,2] =
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n07[i,3] =
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n07[i,4] =
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n07[i,5] =
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n07[i,6] =
      round(parameterEstimates(df_lcs_n07_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n07[i,7] =
      round(fitmeasures(df_lcs_n07_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n07[i,8] =
      round(fitmeasures(df_lcs_n07_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n07[i,9] =
      round(fitmeasures(df_lcs_n07_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n07 = as.data.frame(lcsm_n07)
names(lcsm_n07) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.6 ------------------------------------------------------------
# FE Model
K = 500
fe_n06 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n06 = makedata(N = 500, T = 10, a = -0.6)
  df_fe_n06 = pdata.frame(df_fe_n06,index=c("Unit","Time")) 
  df_fe_n06$Y_lagged = lag(df_fe_n06$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n06, 
               index=c("Unit", "Time"), model="within")
  fe_n06[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n06[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n06[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n06 = as.data.frame(fe_n06)
names(fe_n06) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n06 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n06 = makedata(N = 500, T = 10, a = -0.6)
  df_lcs_n06 = pdata.frame(df_lcs_n06,index=c("Unit","Time")) 
  df_lcs_n06_wide <- reshape(df_lcs_n06, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n06_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n06_lcsm_fit <- sem(df_lcs_n06_lcsm, data= df_lcs_n06_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n06_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n06[i,1] = 
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n06[i,2] =
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n06[i,3] =
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n06[i,4] =
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n06[i,5] =
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n06[i,6] =
      round(parameterEstimates(df_lcs_n06_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n06[i,7] =
      round(fitmeasures(df_lcs_n06_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n06[i,8] =
      round(fitmeasures(df_lcs_n06_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n06[i,9] =
      round(fitmeasures(df_lcs_n06_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n06 = as.data.frame(lcsm_n06)
names(lcsm_n06) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.5 ------------------------------------------------------------
# FE Model
K = 500
fe_n05 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n05 = makedata(N = 500, T = 10, a = -0.5)
  df_fe_n05 = pdata.frame(df_fe_n05,index=c("Unit","Time")) 
  df_fe_n05$Y_lagged = lag(df_fe_n05$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n05, 
               index=c("Unit", "Time"), model="within")
  fe_n05[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n05[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n05[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n05 = as.data.frame(fe_n05)
names(fe_n05) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n05 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n05 = makedata(N = 500, T = 10, a = -0.5)
  df_lcs_n05 = pdata.frame(df_lcs_n05,index=c("Unit","Time")) 
  df_lcs_n05_wide <- reshape(df_lcs_n05, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n05_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n05_lcsm_fit <- sem(df_lcs_n05_lcsm, data= df_lcs_n05_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n05_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n05[i,1] = 
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n05[i,2] =
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n05[i,3] =
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n05[i,4] =
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n05[i,5] =
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n05[i,6] =
      round(parameterEstimates(df_lcs_n05_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n05[i,7] =
      round(fitmeasures(df_lcs_n05_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n05[i,8] =
      round(fitmeasures(df_lcs_n05_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n05[i,9] =
      round(fitmeasures(df_lcs_n05_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n05 = as.data.frame(lcsm_n05)
names(lcsm_n05) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.4 ------------------------------------------------------------
# FE Model
K = 500
fe_n04 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n04 = makedata(N = 500, T = 10, a = -0.4)
  df_fe_n04 = pdata.frame(df_fe_n04,index=c("Unit","Time")) 
  df_fe_n04$Y_lagged = lag(df_fe_n04$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n04, 
               index=c("Unit", "Time"), model="within")
  fe_n04[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n04[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n04[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n04 = as.data.frame(fe_n04)
names(fe_n04) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n04 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n04 = makedata(N = 500, T = 10, a = -0.4)
  df_lcs_n04 = pdata.frame(df_lcs_n04,index=c("Unit","Time")) 
  df_lcs_n04_wide <- reshape(df_lcs_n04, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n04_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n04_lcsm_fit <- sem(df_lcs_n04_lcsm, data= df_lcs_n04_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n04_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n04[i,1] = 
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n04[i,2] =
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n04[i,3] =
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n04[i,4] =
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n04[i,5] =
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n04[i,6] =
      round(parameterEstimates(df_lcs_n04_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n04[i,7] =
      round(fitmeasures(df_lcs_n04_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n04[i,8] =
      round(fitmeasures(df_lcs_n04_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n04[i,9] =
      round(fitmeasures(df_lcs_n04_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n04 = as.data.frame(lcsm_n04)
names(lcsm_n04) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.3 ------------------------------------------------------------
# FE Model
K = 500
fe_n03 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n03 = makedata(N = 500, T = 10, a = -0.3)
  df_fe_n03 = pdata.frame(df_fe_n03,index=c("Unit","Time")) 
  df_fe_n03$Y_lagged = lag(df_fe_n03$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n03, 
               index=c("Unit", "Time"), model="within")
  fe_n03[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n03[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n03[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n03 = as.data.frame(fe_n03)
names(fe_n03) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n03 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n03 = makedata(N = 500, T = 10, a = -0.3)
  df_lcs_n03 = pdata.frame(df_lcs_n03,index=c("Unit","Time")) 
  df_lcs_n03_wide <- reshape(df_lcs_n03, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n03_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n03_lcsm_fit <- sem(df_lcs_n03_lcsm, data= df_lcs_n03_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n03_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n03[i,1] = 
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n03[i,2] =
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n03[i,3] =
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n03[i,4] =
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n03[i,5] =
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n03[i,6] =
      round(parameterEstimates(df_lcs_n03_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n03[i,7] =
      round(fitmeasures(df_lcs_n03_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n03[i,8] =
      round(fitmeasures(df_lcs_n03_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n03[i,9] =
      round(fitmeasures(df_lcs_n03_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n03 = as.data.frame(lcsm_n03)
names(lcsm_n03) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.2 ------------------------------------------------------------
# FE Model
K = 500
fe_n02 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n02 = makedata(N = 500, T = 10, a = -0.2)
  df_fe_n02 = pdata.frame(df_fe_n02,index=c("Unit","Time")) 
  df_fe_n02$Y_lagged = lag(df_fe_n02$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n02, 
               index=c("Unit", "Time"), model="within")
  fe_n02[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n02[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n02[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n02 = as.data.frame(fe_n02)
names(fe_n02) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n02 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n02 = makedata(N = 500, T = 10, a = -0.2)
  df_lcs_n02 = pdata.frame(df_lcs_n02,index=c("Unit","Time")) 
  df_lcs_n02_wide <- reshape(df_lcs_n02, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n02_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n02_lcsm_fit <- sem(df_lcs_n02_lcsm, data= df_lcs_n02_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n02_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n02[i,1] = 
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n02[i,2] =
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n02[i,3] =
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n02[i,4] =
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n02[i,5] =
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n02[i,6] =
      round(parameterEstimates(df_lcs_n02_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n02[i,7] =
      round(fitmeasures(df_lcs_n02_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n02[i,8] =
      round(fitmeasures(df_lcs_n02_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n02[i,9] =
      round(fitmeasures(df_lcs_n02_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n02 = as.data.frame(lcsm_n02)
names(lcsm_n02) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = -0.1 ------------------------------------------------------------
# FE Model
K = 500
fe_n01 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_n01 = makedata(N = 500, T = 10, a = -0.1)
  df_fe_n01 = pdata.frame(df_fe_n01,index=c("Unit","Time")) 
  df_fe_n01$Y_lagged = lag(df_fe_n01$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_n01, 
               index=c("Unit", "Time"), model="within")
  fe_n01[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_n01[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_n01[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_n01 = as.data.frame(fe_n01)
names(fe_n01) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_n01 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_n01 = makedata(N = 500, T = 10, a = -0.1)
  df_lcs_n01 = pdata.frame(df_lcs_n01,index=c("Unit","Time")) 
  df_lcs_n01_wide <- reshape(df_lcs_n01, timevar = "Time", 
                             idvar = "Unit", direction = "wide")
  df_lcs_n01_lcsm <- specify_uni_lcsm(timepoints = 5,
                                      var = "Y.",  
                                      change_letter = "g",
                                      model = list(alpha_constant = TRUE, 
                                                   beta = TRUE, 
                                                   phi = TRUE))
  df_lcs_n01_lcsm_fit <- sem(df_lcs_n01_lcsm, data= df_lcs_n01_wide, 
                             missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_n01_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_n01[i,1] = 
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_n01[i,2] =
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_n01[i,3] =
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_n01[i,4] =
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_n01[i,5] =
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_n01[i,6] =
      round(parameterEstimates(df_lcs_n01_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_n01[i,7] =
      round(fitmeasures(df_lcs_n01_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_n01[i,8] =
      round(fitmeasures(df_lcs_n01_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_n01[i,9] =
      round(fitmeasures(df_lcs_n01_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_n01 = as.data.frame(lcsm_n01)
names(lcsm_n01) = c('Est(b)', 'SE(b)', 'pval(b)', 
                    'Est(u)', 'SE(u)', 'pval(u)', 
                    'rmsea', 'cfi', 'srmr')

# AR(1) = 0 ------------------------------------------------------------
# FE Model
K = 500
fe_0 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_0 = makedata(N = 500, T = 10, a = 0)
  df_fe_0 = pdata.frame(df_fe_0,index=c("Unit","Time")) 
  df_fe_0$Y_lagged = lag(df_fe_0$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_0, 
               index=c("Unit", "Time"), model="within")
  fe_0[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_0[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_0[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_0 = as.data.frame(fe_0)
names(fe_0) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_0 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_0 = makedata(N = 500, T = 10, a = 0)
  df_lcs_0 = pdata.frame(df_lcs_0,index=c("Unit","Time")) 
  df_lcs_0_wide <- reshape(df_lcs_0, timevar = "Time", 
                           idvar = "Unit", direction = "wide")
  df_lcs_0_lcsm <- specify_uni_lcsm(timepoints = 5,
                                    var = "Y.",  
                                    change_letter = "g",
                                    model = list(alpha_constant = TRUE, 
                                                 beta = TRUE, 
                                                 phi = TRUE))
  df_lcs_0_lcsm_fit <- sem(df_lcs_0_lcsm, data= df_lcs_0_wide, 
                           missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_0_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_0[i,1] = 
      round(parameterEstimates(df_lcs_0_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_0[i,2] =
      round(parameterEstimates(df_lcs_0_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_0[i,3] =
      round(parameterEstimates(df_lcs_0_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_0[i,4] =
      round(parameterEstimates(df_lcs_0_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_0[i,5] =
      round(parameterEstimates(df_lcs_0_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_0[i,6] =
      round(parameterEstimates(df_lcs_0_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_0[i,7] =
      round(fitmeasures(df_lcs_0_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_0[i,8] =
      round(fitmeasures(df_lcs_0_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_0[i,9] =
      round(fitmeasures(df_lcs_0_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_0 = as.data.frame(lcsm_0)
names(lcsm_0) = c('Est(b)', 'SE(b)', 'pval(b)', 
                  'Est(u)', 'SE(u)', 'pval(u)', 
                  'rmsea', 'cfi', 'srmr')

# AR(1) = 0.1 ------------------------------------------------------------
# FE Model
K = 500
fe_01 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_01 = makedata(N = 500, T = 10, a = 0.1)
  df_fe_01 = pdata.frame(df_fe_01,index=c("Unit","Time")) 
  df_fe_01$Y_lagged = lag(df_fe_01$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_01, 
               index=c("Unit", "Time"), model="within")
  fe_01[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_01[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_01[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_01 = as.data.frame(fe_01)
names(fe_01) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_01 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_01 = makedata(N = 500, T = 10, a = 0.1)
  df_lcs_01 = pdata.frame(df_lcs_01,index=c("Unit","Time")) 
  df_lcs_01_wide <- reshape(df_lcs_01, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_01_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_01_lcsm_fit <- sem(df_lcs_01_lcsm, data= df_lcs_01_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_01_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_01[i,1] = 
      round(parameterEstimates(df_lcs_01_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_01[i,2] =
      round(parameterEstimates(df_lcs_01_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_01[i,3] =
      round(parameterEstimates(df_lcs_01_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_01[i,4] =
      round(parameterEstimates(df_lcs_01_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_01[i,5] =
      round(parameterEstimates(df_lcs_01_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_01[i,6] =
      round(parameterEstimates(df_lcs_01_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_01[i,7] =
      round(fitmeasures(df_lcs_01_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_01[i,8] =
      round(fitmeasures(df_lcs_01_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_01[i,9] =
      round(fitmeasures(df_lcs_01_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_01 = as.data.frame(lcsm_01)
names(lcsm_01) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.2 ------------------------------------------------------------
# FE Model
K = 500
fe_02 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_02 = makedata(N = 500, T = 10, a = 0.2)
  df_fe_02 = pdata.frame(df_fe_02,index=c("Unit","Time")) 
  df_fe_02$Y_lagged = lag(df_fe_02$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_02, 
               index=c("Unit", "Time"), model="within")
  fe_02[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_02[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_02[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_02 = as.data.frame(fe_02)
names(fe_02) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_02 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_02 = makedata(N = 500, T = 10, a = 0.2)
  df_lcs_02 = pdata.frame(df_lcs_02,index=c("Unit","Time")) 
  df_lcs_02_wide <- reshape(df_lcs_02, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_02_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_02_lcsm_fit <- sem(df_lcs_02_lcsm, data= df_lcs_02_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_02_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_02[i,1] = 
      round(parameterEstimates(df_lcs_02_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_02[i,2] =
      round(parameterEstimates(df_lcs_02_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_02[i,3] =
      round(parameterEstimates(df_lcs_02_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_02[i,4] =
      round(parameterEstimates(df_lcs_02_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_02[i,5] =
      round(parameterEstimates(df_lcs_02_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_02[i,6] =
      round(parameterEstimates(df_lcs_02_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_02[i,7] =
      round(fitmeasures(df_lcs_02_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_02[i,8] =
      round(fitmeasures(df_lcs_02_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_02[i,9] =
      round(fitmeasures(df_lcs_02_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_02 = as.data.frame(lcsm_02)
names(lcsm_02) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.3 ------------------------------------------------------------
# FE Model
K = 500
fe_03 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_03 = makedata(N = 500, T = 10, a = 0.3)
  df_fe_03 = pdata.frame(df_fe_03,index=c("Unit","Time")) 
  df_fe_03$Y_lagged = lag(df_fe_03$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_03, 
               index=c("Unit", "Time"), model="within")
  fe_03[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_03[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_03[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_03 = as.data.frame(fe_03)
names(fe_03) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_03 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_03 = makedata(N = 500, T = 10, a = 0.3)
  df_lcs_03 = pdata.frame(df_lcs_03,index=c("Unit","Time")) 
  df_lcs_03_wide <- reshape(df_lcs_03, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_03_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_03_lcsm_fit <- sem(df_lcs_03_lcsm, data= df_lcs_03_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_03_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_03[i,1] = 
      round(parameterEstimates(df_lcs_03_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_03[i,2] =
      round(parameterEstimates(df_lcs_03_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_03[i,3] =
      round(parameterEstimates(df_lcs_03_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_03[i,4] =
      round(parameterEstimates(df_lcs_03_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_03[i,5] =
      round(parameterEstimates(df_lcs_03_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_03[i,6] =
      round(parameterEstimates(df_lcs_03_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_03[i,7] =
      round(fitmeasures(df_lcs_03_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_03[i,8] =
      round(fitmeasures(df_lcs_03_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_03[i,9] =
      round(fitmeasures(df_lcs_03_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_03 = as.data.frame(lcsm_03)
names(lcsm_03) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.4 ------------------------------------------------------------
# FE Model
K = 500
fe_04 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_04 = makedata(N = 500, T = 10, a = 0.4)
  df_fe_04 = pdata.frame(df_fe_04,index=c("Unit","Time")) 
  df_fe_04$Y_lagged = lag(df_fe_04$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_04, 
               index=c("Unit", "Time"), model="within")
  fe_04[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_04[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_04[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_04 = as.data.frame(fe_04)
names(fe_04) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_04 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_04 = makedata(N = 500, T = 10, a = 0.4)
  df_lcs_04 = pdata.frame(df_lcs_04,index=c("Unit","Time")) 
  df_lcs_04_wide <- reshape(df_lcs_04, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_04_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_04_lcsm_fit <- sem(df_lcs_04_lcsm, data= df_lcs_04_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_04_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_04[i,1] = 
      round(parameterEstimates(df_lcs_04_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_04[i,2] =
      round(parameterEstimates(df_lcs_04_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_04[i,3] =
      round(parameterEstimates(df_lcs_04_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_04[i,4] =
      round(parameterEstimates(df_lcs_04_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_04[i,5] =
      round(parameterEstimates(df_lcs_04_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_04[i,6] =
      round(parameterEstimates(df_lcs_04_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_04[i,7] =
      round(fitmeasures(df_lcs_04_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_04[i,8] =
      round(fitmeasures(df_lcs_04_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_04[i,9] =
      round(fitmeasures(df_lcs_04_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_04 = as.data.frame(lcsm_04)
names(lcsm_04) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.5 ------------------------------------------------------------
# FE Model
K = 500
fe_05 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_05 = makedata(N = 500, T = 10, a = 0.5)
  df_fe_05 = pdata.frame(df_fe_05,index=c("Unit","Time")) 
  df_fe_05$Y_lagged = lag(df_fe_05$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_05, 
               index=c("Unit", "Time"), model="within")
  fe_05[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_05[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_05[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_05 = as.data.frame(fe_05)
names(fe_05) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_05 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_05 = makedata(N = 500, T = 10, a = 0.5)
  df_lcs_05 = pdata.frame(df_lcs_05,index=c("Unit","Time")) 
  df_lcs_05_wide <- reshape(df_lcs_05, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_05_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_05_lcsm_fit <- sem(df_lcs_05_lcsm, data= df_lcs_05_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_05_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_05[i,1] = 
      round(parameterEstimates(df_lcs_05_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_05[i,2] =
      round(parameterEstimates(df_lcs_05_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_05[i,3] =
      round(parameterEstimates(df_lcs_05_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_05[i,4] =
      round(parameterEstimates(df_lcs_05_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_05[i,5] =
      round(parameterEstimates(df_lcs_05_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_05[i,6] =
      round(parameterEstimates(df_lcs_05_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_05[i,7] =
      round(fitmeasures(df_lcs_05_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_05[i,8] =
      round(fitmeasures(df_lcs_05_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_05[i,9] =
      round(fitmeasures(df_lcs_05_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_05 = as.data.frame(lcsm_05)
names(lcsm_05) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.6 ------------------------------------------------------------
# FE Model
K = 500
fe_06 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_06 = makedata(N = 500, T = 10, a = 0.6)
  df_fe_06 = pdata.frame(df_fe_06,index=c("Unit","Time")) 
  df_fe_06$Y_lagged = lag(df_fe_06$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_06, 
               index=c("Unit", "Time"), model="within")
  fe_06[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_06[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_06[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_06 = as.data.frame(fe_06)
names(fe_06) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_06 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_06 = makedata(N = 500, T = 10, a = 0.6)
  df_lcs_06 = pdata.frame(df_lcs_06,index=c("Unit","Time")) 
  df_lcs_06_wide <- reshape(df_lcs_06, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_06_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_06_lcsm_fit <- sem(df_lcs_06_lcsm, data= df_lcs_06_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_06_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_06[i,1] = 
      round(parameterEstimates(df_lcs_06_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_06[i,2] =
      round(parameterEstimates(df_lcs_06_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_06[i,3] =
      round(parameterEstimates(df_lcs_06_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_06[i,4] =
      round(parameterEstimates(df_lcs_06_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_06[i,5] =
      round(parameterEstimates(df_lcs_06_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_06[i,6] =
      round(parameterEstimates(df_lcs_06_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_06[i,7] =
      round(fitmeasures(df_lcs_06_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_06[i,8] =
      round(fitmeasures(df_lcs_06_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_06[i,9] =
      round(fitmeasures(df_lcs_06_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_06 = as.data.frame(lcsm_06)
names(lcsm_06) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.7 ------------------------------------------------------------
# FE Model
K = 500
fe_07 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_07 = makedata(N = 500, T = 10, a = 0.7)
  df_fe_07 = pdata.frame(df_fe_07,index=c("Unit","Time")) 
  df_fe_07$Y_lagged = lag(df_fe_07$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_07, 
               index=c("Unit", "Time"), model="within")
  fe_07[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_07[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_07[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_07 = as.data.frame(fe_07)
names(fe_07) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_07 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_07 = makedata(N = 500, T = 10, a = 0.7)
  df_lcs_07 = pdata.frame(df_lcs_07,index=c("Unit","Time")) 
  df_lcs_07_wide <- reshape(df_lcs_07, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_07_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_07_lcsm_fit <- sem(df_lcs_07_lcsm, data= df_lcs_07_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_07_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_07[i,1] = 
      round(parameterEstimates(df_lcs_07_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_07[i,2] =
      round(parameterEstimates(df_lcs_07_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_07[i,3] =
      round(parameterEstimates(df_lcs_07_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_07[i,4] =
      round(parameterEstimates(df_lcs_07_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_07[i,5] =
      round(parameterEstimates(df_lcs_07_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_07[i,6] =
      round(parameterEstimates(df_lcs_07_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_07[i,7] =
      round(fitmeasures(df_lcs_07_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_07[i,8] =
      round(fitmeasures(df_lcs_07_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_07[i,9] =
      round(fitmeasures(df_lcs_07_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_07 = as.data.frame(lcsm_07)
names(lcsm_07) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.8 ------------------------------------------------------------
# FE Model
K = 500
fe_08 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_08 = makedata(N = 500, T = 10, a = 0.8)
  df_fe_08 = pdata.frame(df_fe_08,index=c("Unit","Time")) 
  df_fe_08$Y_lagged = lag(df_fe_08$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_08, 
               index=c("Unit", "Time"), model="within")
  fe_08[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_08[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_08[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_08 = as.data.frame(fe_08)
names(fe_08) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_08 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_08 = makedata(N = 500, T = 10, a = 0.8)
  df_lcs_08 = pdata.frame(df_lcs_08,index=c("Unit","Time")) 
  df_lcs_08_wide <- reshape(df_lcs_08, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_08_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_08_lcsm_fit <- sem(df_lcs_08_lcsm, data= df_lcs_08_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_08_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_08[i,1] = 
      round(parameterEstimates(df_lcs_08_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_08[i,2] =
      round(parameterEstimates(df_lcs_08_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_08[i,3] =
      round(parameterEstimates(df_lcs_08_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_08[i,4] =
      round(parameterEstimates(df_lcs_08_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_08[i,5] =
      round(parameterEstimates(df_lcs_08_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_08[i,6] =
      round(parameterEstimates(df_lcs_08_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_08[i,7] =
      round(fitmeasures(df_lcs_08_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_08[i,8] =
      round(fitmeasures(df_lcs_08_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_08[i,9] =
      round(fitmeasures(df_lcs_08_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_08 = as.data.frame(lcsm_08)
names(lcsm_08) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# AR(1) = 0.9 ------------------------------------------------------------
# FE Model
K = 500
fe_09 = matrix(nrow = K, ncol = 3)
set.seed(123)
for (i in 1:K){
  df_fe_09 = makedata(N = 500, T = 10, a = 0.9)
  df_fe_09 = pdata.frame(df_fe_09,index=c("Unit","Time")) 
  df_fe_09$Y_lagged = lag(df_fe_09$Y, lag = 1)   
  
  fixed <- plm(Y ~ Y_lagged, data=df_fe_09, 
               index=c("Unit", "Time"), model="within")
  fe_09[i,1] = round(coef(summary(fixed))[1], 2) #est (b)
  fe_09[i,2] = round(coef(summary(fixed))[2], 3) #se (b)
  fe_09[i,3] = round(coef(summary(fixed))[4], 4) #p (b)
  print(i)
}
fe_09 = as.data.frame(fe_09)
names(fe_09) = c('Est(b)', 'SE(b)', 'pval(b)')

#LCS
K = 500
lcsm_09 = matrix(nrow = K, ncol = 9)
set.seed(123)
for (i in 1:K){
  df_lcs_09 = makedata(N = 500, T = 10, a = 0.9)
  df_lcs_09 = pdata.frame(df_lcs_09,index=c("Unit","Time")) 
  df_lcs_09_wide <- reshape(df_lcs_09, timevar = "Time", 
                            idvar = "Unit", direction = "wide")
  df_lcs_09_lcsm <- specify_uni_lcsm(timepoints = 5,
                                     var = "Y.",  
                                     change_letter = "g",
                                     model = list(alpha_constant = TRUE, 
                                                  beta = TRUE, 
                                                  phi = TRUE))
  df_lcs_09_lcsm_fit <- sem(df_lcs_09_lcsm, data= df_lcs_09_wide, 
                            missing = "fiml", estimator ="ml")
  if (lavInspect(df_lcs_09_lcsm_fit, 'converged') == 'TRUE'){
    lcsm_09[i,1] = 
      round(parameterEstimates(df_lcs_09_lcsm_fit)[49,5], 2) #estimate pc (b)
    lcsm_09[i,2] =
      round(parameterEstimates(df_lcs_09_lcsm_fit)[49,6], 3) #se pc (b)
    lcsm_09[i,3] =
      round(parameterEstimates(df_lcs_09_lcsm_fit)[49,8], 4) #pval pc (b)
    lcsm_09[i,4] =
      round(parameterEstimates(df_lcs_09_lcsm_fit)[6,5], 2) #estimate constant change (u)
    lcsm_09[i,5] =
      round(parameterEstimates(df_lcs_09_lcsm_fit)[6,6], 3) #se constant change (u)
    lcsm_09[i,6] =
      round(parameterEstimates(df_lcs_09_lcsm_fit)[6,8], 4) #pval constant change (u)
    lcsm_09[i,7] =
      round(fitmeasures(df_lcs_09_lcsm_fit, c('rmsea')), 3) #rmsea
    lcsm_09[i,8] =
      round(fitmeasures(df_lcs_09_lcsm_fit, c('cfi')), 3) #cfi
    lcsm_09[i,9] =
      round(fitmeasures(df_lcs_09_lcsm_fit, c('srmr')), 3) #srmr
  }
  else{
    next
  }
  print(i)
}
lcsm_09 = as.data.frame(lcsm_09)
names(lcsm_09) = c('Est(b)', 'SE(b)', 'pval(b)', 
                   'Est(u)', 'SE(u)', 'pval(u)', 
                   'rmsea', 'cfi', 'srmr')

# Calculating Bias for FE -------------------------------------------------
#Combine all FE datasets
all_fe_data = rbind(fe_n09, fe_n08, fe_n07, 
                    fe_n06, fe_n05, fe_n04,
                    fe_n03, fe_n02, fe_n01, 
                    fe_0, fe_01, fe_02,
                    fe_03, fe_04, fe_05, 
                    fe_06, fe_07, fe_08,
                    fe_09)
cond = c(rep(-0.9, 500), rep(-0.8, 500), rep(-0.7, 500),
         rep(-0.6, 500), rep(-0.5, 500), rep(-0.4, 500),
         rep(-0.3, 500), rep(-0.2, 500), rep(-0.1, 500),
         rep(0, 500),    rep(0.1, 500), rep(0.2, 500),
         rep(0.3, 500),  rep(0.4, 500), rep(0.5, 500),
         rep(0.6, 500),  rep(0.7, 500), rep(0.8, 500),
         rep(0.9, 500))
all_fe_outputs = cbind(cond, all_fe_data)

#Compute average FE Estimate
Ave_fe_est = rbind(mean(subset(all_fe_outputs, cond == -0.9)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.8)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.7)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.6)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.5)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.4)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.3)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.2)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == -0.1)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.0)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.1)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.2)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.3)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.4)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.5)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.6)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.7)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.8)$`Est(b)`),
                   mean(subset(all_fe_outputs, cond == 0.9)$`Est(b)`))

#Compute bias in average FE Estimate
Bias_fe_est = rbind(mean(abs(subset(all_fe_outputs, cond == -0.9)$`Est(b)`  - (-0.9))),
                    mean(abs(subset(all_fe_outputs, cond == -0.8)$`Est(b)`  - (-0.8))),
                    mean(abs(subset(all_fe_outputs, cond == -0.7)$`Est(b)`  - (-0.7))),
                    mean(abs(subset(all_fe_outputs, cond == -0.6)$`Est(b)`  - (-0.6))),
                    mean(abs(subset(all_fe_outputs, cond == -0.5)$`Est(b)`  - (-0.5))),
                    mean(abs(subset(all_fe_outputs, cond == -0.4)$`Est(b)`  - (-0.4))),
                    mean(abs(subset(all_fe_outputs, cond == -0.3)$`Est(b)`  - (-0.3))),
                    mean(abs(subset(all_fe_outputs, cond == -0.2)$`Est(b)`  - (-0.2))),
                    mean(abs(subset(all_fe_outputs, cond == -0.1)$`Est(b)`  - (-0.1))),
                    mean(abs(subset(all_fe_outputs, cond == 0)$`Est(b)`     - (0))),
                    mean(abs(subset(all_fe_outputs, cond == 0.1)$`Est(b)`  - (0.1))),
                    mean(abs(subset(all_fe_outputs, cond == 0.2)$`Est(b)`  - (0.2))),
                    mean(abs(subset(all_fe_outputs, cond == 0.3)$`Est(b)`  - (0.3))),
                    mean(abs(subset(all_fe_outputs, cond == 0.4)$`Est(b)`  - (0.4))),
                    mean(abs(subset(all_fe_outputs, cond == 0.5)$`Est(b)`  - (0.5))),
                    mean(abs(subset(all_fe_outputs, cond == 0.6)$`Est(b)`  - (0.6))),
                    mean(abs(subset(all_fe_outputs, cond == 0.7)$`Est(b)`  - (0.7))),
                    mean(abs(subset(all_fe_outputs, cond == 0.8)$`Est(b)`  - (0.8))),
                    mean(abs(subset(all_fe_outputs, cond == 0.9)$`Est(b)`  - (0.9))))


# RMSE_fe_est = rbind(mean((subset(all_fe_outputs, cond == -1.5)$`Est(b)`-(-1.5))^2),
#                     mean((subset(all_fe_outputs, cond == -1.0)$`Est(b)`-(-1))^2),
#                     mean((subset(all_fe_outputs, cond == -0.5)$`Est(b)`-(-0.5))^2),
#                     mean((subset(all_fe_outputs, cond == 0)$`Est(b)`-0)^2),
#                     mean((subset(all_fe_outputs, cond == 0.5)$`Est(b)`-0.5)^2),
#                     mean((subset(all_fe_outputs, cond == 1.0)$`Est(b)`-1.0)^2),
#                     mean((subset(all_fe_outputs, cond == 1.5)$`Est(b)`-1.5))^2)


#Compute average FE Standard Errors
Ave_fe_se = rbind(mean(subset(all_fe_outputs, cond == -0.9)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.8)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.7)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.6)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.5)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.4)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.3)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.2)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == -0.1)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.1)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.2)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.3)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.4)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.5)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.6)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.7)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.8)$`SE(b)`),
                  mean(subset(all_fe_outputs, cond == 0.9)$`SE(b)`))



#Calculate Type 2 errors - what percentage were found not significant
for (i in 1:length(all_fe_outputs$cond)){
  if (all_fe_outputs$`pval(b)`[i]>.05){
    all_fe_outputs$Type_2[i] = 1
  }
  else {
    all_fe_outputs$Type_2[i] = 0
  }
}

T2_fe_percent = rbind(mean(subset(all_fe_outputs, cond == -0.9)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.8)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.7)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.6)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.5)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.4)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.3)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.2)$Type_2),
                      mean(subset(all_fe_outputs, cond == -0.1)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.1)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.2)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.3)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.4)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.5)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.6)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.7)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.8)$Type_2),
                      mean(subset(all_fe_outputs, cond == 0.9)$Type_2))

#Calculate the percentage of estimates that fall within the 95% confidence interval
#Create Lower CI
all_fe_outputs$CIL = all_fe_outputs$`Est(b)`-1.96*all_fe_outputs$`SE(b)`
#Create Upper CI
all_fe_outputs$CIU = all_fe_outputs$`Est(b)`+1.96*all_fe_outputs$`SE(b)`
#If statement 
all_fe_outputs$CI = ifelse(all_fe_outputs$cond < all_fe_outputs$CIU & 
                             all_fe_outputs$cond > all_fe_outputs$CIL, 1, 0)

CI_fe_percent = rbind(mean(subset(all_fe_outputs, cond == -0.9)$CI),
                      mean(subset(all_fe_outputs, cond == -0.8)$CI),
                      mean(subset(all_fe_outputs, cond == -0.7)$CI),
                      mean(subset(all_fe_outputs, cond == -0.6)$CI),
                      mean(subset(all_fe_outputs, cond == -0.5)$CI),
                      mean(subset(all_fe_outputs, cond == -0.4)$CI),
                      mean(subset(all_fe_outputs, cond == -0.3)$CI),
                      mean(subset(all_fe_outputs, cond == -0.2)$CI),
                      mean(subset(all_fe_outputs, cond == -0.1)$CI),
                      mean(subset(all_fe_outputs, cond == 0)$CI),
                      mean(subset(all_fe_outputs, cond == 0.1)$CI),
                      mean(subset(all_fe_outputs, cond == 0.2)$CI),
                      mean(subset(all_fe_outputs, cond == 0.3)$CI),
                      mean(subset(all_fe_outputs, cond == 0.4)$CI),
                      mean(subset(all_fe_outputs, cond == 0.5)$CI),
                      mean(subset(all_fe_outputs, cond == 0.6)$CI),
                      mean(subset(all_fe_outputs, cond == 0.7)$CI),
                      mean(subset(all_fe_outputs, cond == 0.8)$CI),
                      mean(subset(all_fe_outputs, cond == 0.9)$CI))

a_coef = seq(-0.9, 0.9, .1)
out_fe = as.data.frame(cbind(a_coef, Ave_fe_est, Bias_fe_est, Ave_fe_se, T2_fe_percent, CI_fe_percent))
names(out_fe) = c('a', 'Ave_fe_est', 'Bias_fe_est', 'Ave_fe_se', 'error_rate', 'CI_fe_est') 

# Calculating Bias for LCSM -------------------------------------------------
#Combine LCSM datasets
all_lcsm_data = rbind(lcsm_n09, lcsm_n08, lcsm_n07, 
                      lcsm_n06, lcsm_n05, lcsm_n04,
                      lcsm_n03, lcsm_n02, lcsm_n01, 
                      lcsm_0,   lcsm_01,  lcsm_02,
                      lcsm_03,  lcsm_04,  lcsm_05, 
                      lcsm_06,  lcsm_07,  lcsm_08,
                      lcsm_09)
lcsm_cond = c(rep(-1.9, 500), rep(-1.8, 500), rep(-1.7, 500),
              rep(-1.6, 500), rep(-1.5, 500), rep(-1.4, 500),
              rep(-1.3, 500), rep(-1.2, 500), rep(-1.1, 500),
              rep(-1.0, 500), rep(-0.9, 500), rep(-0.8, 500),
              rep(-0.7, 500), rep(-0.6, 500), rep(-0.5, 500),
              rep(-0.4, 500), rep(-0.3, 500), rep(-0.2, 500),
              rep(-0.1, 500))

all_lcsm_outputs = cbind(lcsm_cond, all_lcsm_data)

#Compute average LCSM Estimate
Ave_lcsm_est = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$`Est(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$`Est(b)`, na.rm = T))
#Compute Bias in LCSM Estimate
Bias_lcsm_est =rbind(mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.9)$`Est(b)` - (-1.9)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.8)$`Est(b)` - (-1.8)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.7)$`Est(b)` - (-1.7)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.6)$`Est(b)` - (-1.6)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.5)$`Est(b)` - (-1.5)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.4)$`Est(b)` - (-1.4)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.3)$`Est(b)` - (-1.3)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.2)$`Est(b)` - (-1.2)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.1)$`Est(b)` - (-1.1)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -1.0)$`Est(b)` - (-1.0)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.9)$`Est(b)` - (-0.9)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.8)$`Est(b)` - (-0.8)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.7)$`Est(b)` - (-0.7)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.6)$`Est(b)` - (-0.6)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.5)$`Est(b)` - (-0.5)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.4)$`Est(b)` - (-0.4)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.3)$`Est(b)` - (-0.3)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.2)$`Est(b)` - (-0.2)), na.rm = T),
                     mean(abs(subset(all_lcsm_outputs, lcsm_cond == -0.1)$`Est(b)` - (-0.1)), na.rm = T))

#Compute average standard error
Ave_lcsm_se =  rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$`SE(b)`, na.rm = T),
                     mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$`SE(b)`, na.rm = T))


# RMSE_lcsm_est = rbind(mean((subset(all_lcsm_outputs, lcsm_cond == -2.5)$`Est(b)`-(-2.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -2.0)$`Est(b)`-(-2.0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -1.5)$`Est(b)`-(-1.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -1.0)$`Est(b)`-(-1.0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == -0.5)$`Est(b)`-(-0.5))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == 0)$`Est(b)`-(0))^2, na.rm = T),
#                       mean((subset(all_lcsm_outputs, lcsm_cond == 0.5)$`Est(b)`-(.5))^2, na.rm = T))

#Compute average fit index estimate
Ave_CFI = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$cfi, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$cfi, na.rm = T))

Ave_RMSEA=rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$rmsea, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$rmsea, na.rm = T))

Ave_SRMR =rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$srmr, na.rm = T),
                mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$srmr, na.rm = T))

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

T2_lcsm_percent = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$Type_2, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$Type_2, na.rm = T))


#Calculate the percentage of estimates that fall within the 95% confidence interval
#Create Lower CI
all_lcsm_outputs$CIL = all_lcsm_outputs$`Est(b)`-1.96*all_lcsm_outputs$`SE(b)`
#Create Upper CI
all_lcsm_outputs$CIU = all_lcsm_outputs$`Est(b)`+1.96*all_lcsm_outputs$`SE(b)`
#If statement 
all_lcsm_outputs$CI = ifelse(is.na(all_lcsm_outputs$CIU), 0, 
                             ifelse(all_lcsm_outputs$lcsm_cond < all_lcsm_outputs$CIU & 
                                      all_lcsm_outputs$lcsm_cond > all_lcsm_outputs$CIL, 1, 0))

CI_lcsm_percent = rbind(mean(subset(all_lcsm_outputs, lcsm_cond == -1.9)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.8)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.7)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.6)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.5)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.4)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.3)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.2)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.1)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -1.0)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.9)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.8)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.7)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.6)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.5)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.4)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.3)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.2)$CI, na.rm = T),
                        mean(subset(all_lcsm_outputs, lcsm_cond == -0.1)$CI, na.rm = T))



p_coef = seq(-1.9, -0.1, .1)
out_lcsm = as.data.frame(cbind(p_coef, Ave_lcsm_est, Bias_lcsm_est, Ave_lcsm_se, T2_lcsm_percent,
                               CI_lcsm_percent, Ave_CFI, Ave_RMSEA, Ave_SRMR))
names(out_lcsm) = c('p', 'Ave_est', 'Bias_est', 'Ave_se', 'error_rate', 'CI_lcsm_est', 'Ave_CFI', 'Ave_RMSEA', 'Ave_SRMR')                

#final table

bias_table_out = cbind(out_fe, out_lcsm)

