# ANOCOVA Corn Yield Mapping 2019
# Created By A.J. Brown
# Updated 06 April 2020


library(ggplot2)
library(gridExtra)
library(ggpubr)
library(hydroGOF)
library(psych)
library(minpack.lm)
library(caret)
library(nlstools)
library(robustbase)

obs_df = read.csv(file.choose()) #use "1819_ObsOnly.csv"
corn_yield1 = read.csv(file.choose()) #use "1819_allFields.csv"
corn_yield1$ln1.0m = log(corn_yield1$CV.1.0m)
corn_yield1$ln0.5m = log(corn_yield1$CV.0.5m)
corn_yield1$lnkg = log(corn_yield1$kg.ha)
tapply(corn_yield1$kg.ha, corn_yield1$Field, describe)

#head(corn_yield1)
typeof(corn_yield1$kg.ha) #should be 'double'
corn_yield1$Field = as.factor(corn_yield1$Field)
obs_df$Field = as.factor(obs_df$Field)
levels(corn_yield1$Field)
corn_yield1$Field = factor(corn_yield1$Field,levels(corn_yield1$Field)[c(5,1,2,3,4,6,7,8,9,10)])
levels(corn_yield1$Field)

#Removing fields b/c salt was not dominant yield impact factor
#Rule: remove if max ECe <1.7 for corn, or pearsons correlation not significant (p<0.05)
tapply(corn_yield1$ECe_Deep, corn_yield1$Field, describe)
tapply(corn_yield1$ECe_Gyp_Call, corn_yield1$Field, describe)
df = subset(corn_yield1, Field=="G2")
cor.test(df$ECe_Gyp_Call, df$Yr_., short=FALSE) #remove G2 b/c pearsons not significant
corn_yield = corn_yield1
corn_yield = droplevels(subset(corn_yield1, Field!='M23' & Field!='M19' & Field!='G8' & Field!='G2'))
obs_df = droplevels(subset(obs_df, Field!='M23' & Field!='M19' & Field!='G8' & Field!='G2'))
levels(obs_df$Field)
levels(corn_yield$Field)

#pre-generation of columns needed for proper calculations later
estimates = function(){
  corn_yield$relative_yield <<- c(1:length(corn_yield$ID))
  corn_yield$ece_relative_yield <<- c(1:length(corn_yield$ID))
  corn_yield$eca_relative_yield <<- c(1:length(corn_yield$ID))
  corn_yield$anocova_eca_relative_yield <<- c(1:length(corn_yield$ID))
  corn_yield$maas_relative_yield <<- c(1:length(corn_yield$ID))
  corn_yield$maas_yield <<- c(1:length(corn_yield$ID))
  obs_df$relative_yield <<- c(1:length(obs_df$ID))
  obs_df$ece_relative_yield <<- c(1:length(obs_df$ID))
  obs_df$eca_relative_yield <<- c(1:length(obs_df$ID))
  obs_df$anocova_eca_relative_yield <<- c(1:length(obs_df$ID))
  obs_df$maas_relative_yield <<- c(1:length(obs_df$ID))
  obs_df$maas_yield <<- c(1:length(obs_df$ID))
}
estimates()

#Summary statistics for 1m and 0.5m ECa directly to yield models via ANOCOVA model
Summary_stats = function() {
  #Correlation of ECa to Yield
  cor1.0m = cor(corn_yield$kg.ha, corn_yield$CV.1.0m)
  cor0.5m = cor(corn_yield$kg.ha, corn_yield$CV.0.5m)

  #Summary of Statistics - to get into bu/ac, change bu_kg to equal 0.0159
  Summary = data.frame(
    Model = c('1.5m', '0.75m'),
    Corr = c(cor1.0m, cor0.5m)
  )
  Summary
}
Summary_stats()

##################Relative Yield, Yr with ECe, and ECa with linear and nonlinear
corn_yield$ECe = corn_yield$ECe_Deep#/0.7819-.4013
#corn_yield$ECeg = ifelse(corn_yield$ECe<3.08, corn_yield$ECe, (corn_yield$ECe*0.9131)-0.7256)
corn_yield$ECeg = corn_yield$ECe_Gyp_Call #Run for Callaghan (2016) ECeg
obs_df$ECe = obs_df$ECe_Deep#/0.7819-.4013
#obs_df$ECeg = ifelse(obs_df$ECe<3.08, obs_df$ECe, (obs_df$ECe*0.9131)-0.7256)
obs_df$ECeg = obs_df$ECe_Gyp_Call
#corn_yield$ECeg = corn_yield$ECe_Gyp_Bern #Run for Selective Dilution ECeg
range(corn_yield$ECe) #should be: 1.215487 - 10.198949
range(corn_yield$ECeg) #should be: 1.293638 - 7.739676
range(obs_df$ECe) #should be: 1.266667 - 11.206667
range(obs_df$ECeg) #should be: 1.266667 - 9.507207
######################################separation of testing and calibration data
sample_size = 20 #either 38, 20, 12, or 6
corn_yield_cal = subset(corn_yield, Field.ID<sample_size) #use to change sample size
corn_yield_test = subset(corn_yield, Field.ID>=sample_size) #use to change sample size
################################################################################
#calculation of observed relative yield df for model comparisons
max_yield_df = function(df, avgtop) {
  #Ideas on theoretical max yield:
  # -why average the top 3 yields for a "max" in each field? Maybe We should average all of the yields
  #  below the salinity thresholds (1.7, and then try 3.7 for gypsum shift).  This would average out
  #  the other potential factors reducing/bolstering yields. PROS: This would allow our 'scatter' to
  #  shift outward, showing the potential gypsum effects. CONS: We will have values above 100% relative
  #  yield (but this could be okay in a field setting such as ours, which isn't a controlled
  #  environment).

  #calculation of average top 3 yields in each field for relative yield calcs
  field_list = c()
  top_yld = c() #Top averaged yield in each field
  top_ece = c() #Top averaged ECe level in each field
  top_eca = c() #Top averaged ECa level in each field
  top_eceg = c() #Top averaged ECeg level in each field
  half_yld = c() #observed yield closest to top_yld/2
  half_ece = c() #ECe level at half yield location
  half_eca = c() #ECa level at half yield location
  half_eceg = c() #ECeg level at half yield location
  avg_top = avgtop #number of locations that we'd like to average to obtain "top"
  #note: I had to put this^^^ to 1 in order to work with 6 sample size sets

  for(i in levels(df$Field)){
    df2 = subset(df, Field == i)
    top_yld1 = mean(tail(sort(df2$kg.ha), avg_top))
    top_ece1 = mean(tail(sort(df2$ECe), avg_top))
    top_eca1 = mean(tail(sort(df2$CV.1.0m), avg_top))
    top_eceg1 = mean(tail(sort(df2$ECeg), avg_top))
    half_yld1 = df2$kg.ha[which.min(abs(df2$kg.ha - top_yld1/2))]
    half_ece1 = subset(df2, kg.ha == half_yld1)$ECe
    half_eca1 = subset(df2, kg.ha == half_yld1)$CV.1.0m
    half_eceg1 = subset(df2, kg.ha == half_yld1)$ECeg
    top_yld = append(top_yld, top_yld1)
    top_ece = append(top_ece, top_ece1)
    top_eca = append(top_eca, top_eca1)
    top_eceg = append(top_eceg, top_eceg1)
    half_yld = append(half_yld, half_yld1)
    half_ece = append(half_ece, half_ece1)
    half_eca = append(half_eca, half_eca1)
    half_eceg = append(half_eceg, half_eceg1)
    field_list = append(field_list, i)
  }
  max_yields <<- data.frame(Field = field_list,
                            max_yield = top_yld,
                            max_ECe = top_ece,
                            max_ECeg = top_eceg,
                            max_ECa = top_eca,
                            max_lnECa = log(top_eca),
                            half_yield = half_yld,
                            half_ECe = half_ece,
                            half_ECeg = half_eceg,
                            half_ECa = half_eca,
                            half_lnECa = log(half_eca))
  print(max_yields)
  return(max_yields)
}
max_yields_all = max_yield_df(corn_yield, avgtop = 1)
max_yields_obs = max_yield_df(obs_df, avgtop = 1)
max_yields_cal = max_yield_df(corn_yield_cal, avgtop = 1) #use avgtop=3 for 12+ sample design
#max_yields_cal = max_yield_df(corn_yield_cal, avgtop = 1) #for 6- sample design

#creating RY% columns in full and calibration df's
actual_rel_yield = function(df, graph, max) {
  #calculation of actual relative yield % for dataframe
  for(i in 1:length(df$ID)) {
    x = df$Field
    y = df$kg.ha
    z = df$CV.1.0m
    df$relative_yield[i] <- (y[i]/max$max_yield[max$Field == x[i]])*100
    df$relative_ECa[i] <- (z[i]/max$max_ECa[max$Field == x[i]])*100
  }
  df$ln_relative_yield <- log(df$relative_yield)
  if (graph == T){
    a=qplot(CV.1.0m, relative_yield, data=df, col=Field)+stat_smooth(method='lm', se=F)
    b=qplot(ln1.0m, log(relative_yield), data=df, col=Field)+stat_smooth(method='lm', se=F)
    c=qplot(ECe, relative_yield, data=df, col=Field)+stat_smooth(method='lm', se=F)
    d=qplot(ECeg, relative_yield, data=df, col=Field)+stat_smooth(method='lm', se=F)
    grid.arrange(a,b,c,d, nrow=2, ncol=2)
  }
  return(df)
}
corn_yield = actual_rel_yield(df = corn_yield, graph = F, max=max_yields_all)
obs_df = actual_rel_yield(df = obs_df, graph = F, max=max_yields_obs)
corn_yield_cal = actual_rel_yield(df = corn_yield_cal, graph = F, max=max_yields_cal)
corn_yield_test = actual_rel_yield(df = corn_yield_test, graph = F, max=max_yields_cal)

#test for correlation between ECe and Yr
by(corn_yield[,c("ECe_Deep","relative_yield")], corn_yield$Field, corr.test)

#predictions of RY and Y from linear and ANOCOVA models with ECe and ECa
#max variable is the max_yields df from above
linear_fit = function(df, test_df, full_data, max) {
  #calculation of ECe relative yield (with linear regression ECe fitting)
  ece_model <<- lm(relative_yield ~ ECe, data = df)
    print(summary(ece_model))
    print("Confidence Intervals:#############################")
    print(confint(ece_model))
    print("Confidence Intervals:#############################")
    df$ece_relative_yield <- predict(ece_model, data = df)
    test_df$ece_relative_yield_test <- predict(ece_model, newdata = test_df)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$ece_relative_yield
      df$ece_yield_pred_from_rel[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$ece_relative_yield_test
      test_df$ece_yield_pred_from_rel_test[i] <- ((y2[i]/100)*max$max_yield[max$Field == x2[i]])
    }

  eceg_model <<- lm(relative_yield ~ ECeg, data = df)
    print(summary(eceg_model))
    print("Confidence Intervals:#############################")
    print(confint(eceg_model))
    print("Confidence Intervals:#############################")
    df$eceg_relative_yield <- predict(eceg_model, data = df)
    test_df$eceg_relative_yield_test <- predict(eceg_model, newdata = test_df)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$eceg_relative_yield
      df$eceg_yield_pred_from_rel[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$eceg_relative_yield_test
      test_df$eceg_yield_pred_from_rel_test[i] <- ((y2[i]/100)*max$max_yield[max$Field == x2[i]])
    }

  #calculation of ECa relative yield (with linear regression ECa fitting)
  eca_model <<- lm(relative_yield ~ CV.1.0m, data = df)
    print(summary(eca_model))
    print("Confidence Intervals:#############################")
    print(confint(eca_model))
    print("Confidence Intervals:#############################")
    df$eca_relative_yield <- predict(eca_model, data = df)
    test_df$eca_relative_yield_test <- predict(eca_model, newdata = test_df)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$eca_relative_yield
      df$eca_yield_pred_from_rel[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$eca_relative_yield_test
      test_df$eca_yield_pred_from_rel_test[i] <- ((y2[i]/100)*max$max_yield[max$Field == x2[i]])
    }

  #calculation of ECa relative yield (with ANOCOVA regression ECa fitting)
  anocova_eca_model <<- lm(relative_yield ~ CV.1.0m+Field, data = df)
    print(summary(anocova_eca_model))
    #print(confint(anocova_eca_model))
    df$anocova_eca_relative_yield <- predict(anocova_eca_model, data = df)
    test_df$anocova_eca_relative_yield_test <- predict(anocova_eca_model, newdata = test_df)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$anocova_eca_relative_yield
      df$anocova_eca_yield_pred_from_rel[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$anocova_eca_relative_yield_test
      test_df$anocova_eca_yield_pred_from_rel_test[i] <- ((y2[i]/100)*max$max_yield[max$Field == x2[i]])
    }

  #calculation of relative yield thresholds (with Maas and Hoffman, 1977 method)
  salt_threshold = 1.7
  salt_slope = 12
  df$maas_relative_yield <- ifelse(df$ECe >= salt_threshold, 100 - salt_slope*(df$ECe - salt_threshold), 100)
  test_df$maas_relative_yield_test <- ifelse(test_df$ECe >= salt_threshold, 100 - salt_slope*(test_df$ECe - salt_threshold), 100)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$maas_relative_yield
      df$maas_yield[i] <- (y[i]/100*max_yields$max_yield[max_yields$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$maas_relative_yield_test
      test_df$maas_yield_test[i] <- (y2[i]/100*max_yields$max_yield[max_yields$Field == x2[i]])
    }

  #Applying ECeg to Maas Hoffman eq.
  df$eceg_maas_relative_yield <- ifelse(df$ECeg >= salt_threshold, 100 - salt_slope*(df$ECeg - salt_threshold), 100)
  test_df$eceg_maas_relative_yield_test <- ifelse(test_df$ECeg >= salt_threshold, 100 - salt_slope*(test_df$ECeg - salt_threshold), 100)
    for(i in 1:length(df$ID)) {
      x = df$Field
      y = df$eceg_maas_relative_yield
      df$eceg_maas_yield[i] <- (y[i]/100*max_yields$max_yield[max_yields$Field == x[i]])
    }
    for(i in 1:length(test_df$ID)) {
      x2 = test_df$Field
      y2 = test_df$eceg_maas_relative_yield_test
      test_df$eceg_maas_yield_test[i] <- (y2[i]/100*max_yields$max_yield[max_yields$Field == x2[i]])
    }

  #RMSE
  RMSE_ECa_model = rmse(df$eca_yield_pred_from_rel, df$kg.ha)
  RMSE_ECe_model = rmse(df$ece_yield_pred_from_rel, df$kg.ha)
  RMSE_ECeg_model = rmse(df$eceg_yield_pred_from_rel, df$kg.ha)
  RMSE_Maas_ECe_model = rmse(df$maas_yield, df$kg.ha)
  RMSE_Maas_ECeg_model = rmse(df$eceg_maas_yield, df$kg.ha)
  RMSE_ECa_anocova = rmse(df$anocova_eca_yield_pred_from_rel, df$kg.ha)

  RMSE_Rel_ECa_model = rmse(df$relative_yield, df$eca_relative_yield)
  RMSE_Rel_ECe_model = rmse(df$relative_yield, df$ece_relative_yield)
  RMSE_Rel_ECeg_model = rmse(df$relative_yield, df$eceg_relative_yield)
  RMSE_Rel_Maas_ECe_model = rmse(df$relative_yield, df$maas_relative_yield)
  RMSE_Rel_Maas_ECeg_model = rmse(df$relative_yield, df$eceg_maas_relative_yield)
  RMSE_Rel_ECa_anocova = rmse(df$anocova_eca_relative_yield, df$relative_yield)

  #NRMSE
  NRMSE_ECa_model = nrmse(df$eca_yield_pred_from_rel, df$kg.ha)
  NRMSE_ECe_model = nrmse(df$ece_yield_pred_from_rel, df$kg.ha)
  NRMSE_ECeg_model = nrmse(df$eceg_yield_pred_from_rel, df$kg.ha)
  NRMSE_Maas_ECe_model = nrmse(df$maas_yield, df$kg.ha)
  NRMSE_Maas_ECeg_model = nrmse(df$eceg_maas_yield, df$kg.ha)
  NRMSE_ECa_anocova = nrmse(df$anocova_eca_yield_pred_from_rel, df$kg.ha)

  NRMSE_Rel_ECa_model = nrmse(df$relative_yield, df$eca_relative_yield)
  NRMSE_Rel_ECe_model = nrmse(df$relative_yield, df$ece_relative_yield)
  NRMSE_Rel_ECeg_model = nrmse(df$relative_yield, df$eceg_relative_yield)
  NRMSE_Rel_Maas_ECe_model = nrmse(df$relative_yield, df$maas_relative_yield)
  NRMSE_Rel_Maas_ECeg_model = nrmse(df$relative_yield, df$eceg_maas_relative_yield)
  NRMSE_Rel_ECa_anocova = nrmse(df$anocova_eca_relative_yield, df$relative_yield)

  #Correlation Statistic
  Cor_ECa_model = cor(df$kg.ha, df$eca_yield_pred_from_rel)
  Cor_ECe_model = cor(df$kg.ha, df$ece_yield_pred_from_rel)
  Cor_ECeg_model = cor(df$kg.ha, df$eceg_yield_pred_from_rel)
  Cor_Maas_ECe_model = cor(df$kg.ha, df$maas_yield)
  Cor_Maas_ECeg_model = cor(df$kg.ha, df$eceg_maas_yield)
  Cor_ECa_anocova = cor(df$kg.ha, df$anocova_eca_yield_pred_from_rel)

  Cor_Rel_ECa_model = cor(df$relative_yield, df$eca_relative_yield)
  Cor_Rel_ECe_model = cor(df$relative_yield, df$ece_relative_yield)
  Cor_Rel_ECeg_model = cor(df$relative_yield, df$eceg_relative_yield)
  Cor_Rel_Maas_ECe_model = cor(df$relative_yield, df$maas_relative_yield)
  Cor_Rel_Maas_ECeg_model = cor(df$relative_yield, df$eceg_maas_relative_yield)
  Cor_Rel_ECa_anocova = cor(df$relative_yield,df$anocova_eca_relative_yield)

  #index of agreement (IOE) https://www.rforge.net/doc/packages/hydroGOF/d.html
  IOE_ECa_model = d(df$eca_yield_pred_from_rel, df$kg.ha)
  IOE_ECe_model = d(df$ece_yield_pred_from_rel, df$kg.ha)
  IOE_ECeg_model = d(df$eceg_yield_pred_from_rel, df$kg.ha)
  IOE_Maas_ECe_model = d(df$maas_yield, df$kg.ha)
  IOE_Maas_ECeg_model = d(df$eceg_maas_yield, df$kg.ha)
  IOE_ECa_anocova = d(df$anocova_eca_yield_pred_from_rel, df$kg.ha)

  IOE_Rel_ECa_model = d(df$relative_yield, df$eca_relative_yield)
  IOE_Rel_ECe_model = d(df$relative_yield, df$ece_relative_yield)
  IOE_Rel_ECeg_model = d(df$relative_yield, df$eceg_relative_yield)
  IOE_Rel_Maas_ECe_model = d(df$relative_yield, df$maas_relative_yield)
  IOE_Rel_Maas_ECeg_model = d(df$relative_yield, df$eceg_maas_relative_yield)
  IOE_Rel_ECa_anocova = d(df$relative_yield, df$anocova_eca_relative_yield)

  #Jack-Knifed Mean Squared Prediction Error - mspe
  if (full_data == T) {
    #Linear ECa models
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$CV.1.0m
      y = test$relative_yield
      mdl = lm(y ~ x)
      ydf = data.frame(x=lofo$CV.1.0m, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df = data.frame(Left_out_Field = field_list,
                            ECa_Linear_RMSPE = rmspe_list)
    #Linear ECe models
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$ECe
      y = test$relative_yield
      mdl = lm(y ~ x)
      ydf = data.frame(x=lofo$ECe, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df1 = data.frame(Left_out_Field = field_list,
                            ECe_Linear_RMSPE = rmspe_list)

    #Linear ECeg models
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$ECeg
      y = test$relative_yield
      mdl = lm(y ~ x)
      ydf = data.frame(x=lofo$ECeg, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df2 = data.frame(Left_out_Field = field_list,
                            ECeg_Linear_RMSPE = rmspe_list)

      #mergine RMSPE results
      #master_df = merge(final_df, final_df1, final_df2, by='Left_out_Field')
      master_df <- Reduce(function(x, y) merge(x, y, by='Left_out_Field', all=TRUE), list(final_df, final_df1, final_df2))
      mergedData <- merge(master_df, max,
                          by.x=c('Left_out_Field'), by.y=c('Field'))
      mergedData$Lin_ECa_RMSPE_Y = (mergedData$ECa_Linear_RMSPE/100)*mergedData$max_yield
      mergedData$Lin_ECe_RMSPE_Y = (mergedData$ECe_Linear_RMSPE/100)*mergedData$max_yield
      mergedData$Lin_ECeg_RMSPE_Y = (mergedData$ECeg_Linear_RMSPE/100)*mergedData$max_yield
      keeps <- c('Left_out_Field',
                 'ECa_Linear_RMSPE',
                 'ECe_Linear_RMSPE',
                 'ECeg_Linear_RMSPE',
                 'Lin_ECa_RMSPE_Y',
                 'Lin_ECe_RMSPE_Y',
                 'Lin_ECeg_RMSPE_Y')
      print(mergedData[keeps])
    }
  else {
    #RMSE
    RMSE_ECa_model_test = rmse(test_df$eca_yield_pred_from_rel_test, test_df$kg.ha)
    RMSE_ECe_model_test = rmse(test_df$ece_yield_pred_from_rel_test, test_df$kg.ha)
    RMSE_ECeg_model_test = rmse(test_df$eceg_yield_pred_from_rel_test, test_df$kg.ha)
    RMSE_Maas_ECe_model_test = rmse(test_df$maas_yield_test, test_df$kg.ha)
    RMSE_Maas_ECeg_model_test = rmse(test_df$eceg_maas_yield_test, test_df$kg.ha)
    RMSE_ECa_anocova_test = rmse(test_df$anocova_eca_yield_pred_from_rel_test, test_df$kg.ha)

    RMSE_Rel_ECa_model_test = rmse(test_df$relative_yield, test_df$eca_relative_yield_test)
    RMSE_Rel_ECe_model_test = rmse(test_df$relative_yield, test_df$ece_relative_yield_test)
    RMSE_Rel_ECeg_model_test = rmse(test_df$relative_yield, test_df$eceg_relative_yield_test)
    RMSE_Rel_Maas_ECe_model_test = rmse(test_df$relative_yield, test_df$maas_relative_yield_test)
    RMSE_Rel_Maas_ECeg_model_test = rmse(test_df$relative_yield, test_df$eceg_maas_relative_yield_test)
    RMSE_Rel_ECa_anocova_test = rmse(test_df$relative_yield, test_df$anocova_eca_relative_yield_test)

    #NRMSE
    NRMSE_ECa_model_test = nrmse(test_df$eca_yield_pred_from_rel_test, test_df$kg.ha)
    NRMSE_ECe_model_test = nrmse(test_df$ece_yield_pred_from_rel_test, test_df$kg.ha)
    NRMSE_ECeg_model_test = nrmse(test_df$eceg_yield_pred_from_rel_test, test_df$kg.ha)
    NRMSE_Maas_ECe_model_test = nrmse(test_df$maas_yield_test, test_df$kg.ha)
    NRMSE_Maas_ECeg_model_test = nrmse(test_df$eceg_maas_yield_test, test_df$kg.ha)
    NRMSE_ECa_anocova_test = nrmse(test_df$anocova_eca_yield_pred_from_rel_test, test_df$kg.ha)

    NRMSE_Rel_ECa_model_test = nrmse(test_df$relative_yield, test_df$eca_relative_yield_test)
    NRMSE_Rel_ECe_model_test = nrmse(test_df$relative_yield, test_df$ece_relative_yield_test)
    NRMSE_Rel_ECeg_model_test = nrmse(test_df$relative_yield, test_df$eceg_relative_yield_test)
    NRMSE_Rel_Maas_ECe_model_test = nrmse(test_df$relative_yield, test_df$maas_relative_yield_test)
    NRMSE_Rel_Maas_ECeg_model_test = nrmse(test_df$relative_yield, test_df$eceg_maas_relative_yield_test)
    NRMSE_Rel_ECa_anocova_test = nrmse(test_df$relative_yield, test_df$anocova_eca_relative_yield_test)

    #Pearson's Correlation Statistic
    Cor_ECa_model_test = cor(test_df$kg.ha, test_df$eca_yield_pred_from_rel_test)
    Cor_ECe_model_test = cor(test_df$kg.ha, test_df$ece_yield_pred_from_rel_test)
    Cor_ECeg_model_test = cor(test_df$kg.ha, test_df$eceg_yield_pred_from_rel_test)
    Cor_Maas_ECe_model_test = cor(test_df$kg.ha, test_df$maas_yield_test)
    Cor_Maas_ECeg_model_test = cor(test_df$kg.ha, test_df$eceg_maas_yield_test)
    Cor_ECa_anocova_test = cor(test_df$kg.ha, test_df$anocova_eca_yield_pred_from_rel_test)

    Cor_Rel_ECa_model_test = cor(test_df$relative_yield, test_df$eca_relative_yield_test)
    Cor_Rel_ECe_model_test = cor(test_df$relative_yield, test_df$ece_relative_yield_test)
    Cor_Rel_ECeg_model_test = cor(test_df$relative_yield, test_df$eceg_relative_yield_test)
    Cor_Rel_Maas_ECe_model_test = cor(test_df$relative_yield, test_df$maas_relative_yield_test)
    Cor_Rel_Maas_ECeg_model_test = cor(test_df$relative_yield, test_df$eceg_maas_relative_yield_test)
    Cor_Rel_ECa_anocova_test = cor(test_df$relative_yield,test_df$anocova_eca_relative_yield_test)

    #index of agreement (IOE) https://www.rforge.net/doc/packages/hydroGOF/d.html
    IOE_ECa_model_test = d(test_df$eca_yield_pred_from_rel_test, test_df$kg.ha)
    IOE_ECe_model_test = d(test_df$ece_yield_pred_from_rel_test, test_df$kg.ha)
    IOE_ECeg_model_test = d(test_df$eceg_yield_pred_from_rel_test, test_df$kg.ha)
    IOE_Maas_ECe_model_test = d(test_df$maas_yield_test, test_df$kg.ha)
    IOE_Maas_ECeg_model_test = d(test_df$eceg_maas_yield_test, test_df$kg.ha)
    IOE_ECa_anocova_test = d(test_df$anocova_eca_yield_pred_from_rel_test, test_df$kg.ha)

    IOE_Rel_ECa_model_test = d(test_df$relative_yield, test_df$eca_relative_yield_test)
    IOE_Rel_ECe_model_test = d(test_df$relative_yield, test_df$ece_relative_yield_test)
    IOE_Rel_ECeg_model_test = d(test_df$relative_yield, test_df$eceg_relative_yield_test)
    IOE_Rel_Maas_ECe_model_test = d(test_df$relative_yield, test_df$maas_relative_yield_test)
    IOE_Rel_Maas_ECeg_model_test = d(test_df$relative_yield, test_df$eceg_maas_relative_yield_test)
    IOE_Rel_ECa_anocova_test = d(test_df$relative_yield, test_df$anocova_eca_relative_yield_test)
    #Summary of Statistics - to get into bu/ac, change bu_kg to equal 0.0159
    print("test dataset statistics:###########################################")
    bu_kg = 1
    Summary_df = data.frame(
      Model = c('ECa_Linear_model', 'ECe_Linear_model', 'ECeg_Linear_model','ECa_ANOCOVA_model', 'Maas_ECe_model', 'Maas_ECeg_model'),
      NRMSE = c(NRMSE_ECa_model_test, NRMSE_ECe_model_test, NRMSE_ECeg_model_test, NRMSE_ECa_anocova_test, NRMSE_Maas_ECe_model_test, NRMSE_Maas_ECeg_model_test),
      RMSE = c(RMSE_ECa_model_test*bu_kg, RMSE_ECe_model_test*bu_kg, RMSE_ECeg_model_test*bu_kg, RMSE_ECa_anocova_test*bu_kg, RMSE_Maas_ECe_model_test*bu_kg, RMSE_Maas_ECeg_model_test*bu_kg),
      Corr = c(Cor_ECa_model_test, Cor_ECe_model_test, Cor_ECeg_model_test, Cor_ECa_anocova_test, Cor_Maas_ECe_model_test, Cor_Maas_ECeg_model_test),
      IOE = c(IOE_ECa_model_test, IOE_ECe_model_test, IOE_ECeg_model_test, IOE_ECa_anocova_test, IOE_Maas_ECe_model_test, IOE_Maas_ECeg_model_test)
    )
    print(Summary_df)
    Summary_df1 = data.frame(
      Model = c('ECa_Linear_Relative', 'ECe_Linear_Relative', 'ECeg_Linear_Relative', 'ECa_ANOCOVA_Relative', 'Maas_Relative', 'ECeg_Maas_Relative'),
      NRMSE = c(NRMSE_Rel_ECa_model_test, NRMSE_Rel_ECe_model_test, NRMSE_Rel_ECeg_model_test, NRMSE_Rel_ECa_anocova_test, NRMSE_Rel_Maas_ECe_model_test, NRMSE_Rel_Maas_ECeg_model_test),
      RMSE = c(RMSE_Rel_ECa_model_test, RMSE_Rel_ECe_model_test, RMSE_Rel_ECeg_model_test, RMSE_Rel_ECa_anocova_test, RMSE_Rel_Maas_ECe_model_test, RMSE_Rel_Maas_ECeg_model_test),
      Corr = c(Cor_Rel_ECa_model_test, Cor_Rel_ECe_model_test, Cor_Rel_ECeg_model_test, Cor_Rel_ECa_anocova_test, Cor_Rel_Maas_ECe_model_test, Cor_Rel_Maas_ECeg_model_test),
      IOE = c(IOE_Rel_ECa_model_test, IOE_Rel_ECe_model_test, IOE_Rel_ECeg_model_test, IOE_Rel_ECa_anocova_test,  IOE_Rel_Maas_ECe_model_test, IOE_Rel_Maas_ECeg_model_test)
    )
    print(Summary_df1)
  }
  #Summary of Statistics - to get into bu/ac, change bu_kg to equal 0.0159
  print("full dataset statistics:###########################################")
  bu_kg = 1
  Summary = data.frame(
    Model = c('ECa_Linear_model', 'ECe_Linear_model', 'ECeg_Linear_model','ECa_ANOCOVA_model', 'Maas_ECe_model', 'Maas_ECeg_model'),
    NRMSE = c(NRMSE_ECa_model, NRMSE_ECe_model, NRMSE_ECeg_model, NRMSE_ECa_anocova, NRMSE_Maas_ECe_model, NRMSE_Maas_ECeg_model),
    RMSE = c(RMSE_ECa_model*bu_kg, RMSE_ECe_model*bu_kg, RMSE_ECeg_model*bu_kg, RMSE_ECa_anocova*bu_kg, RMSE_Maas_ECe_model*bu_kg, RMSE_Maas_ECeg_model*bu_kg),
    Corr = c(Cor_ECa_model, Cor_ECe_model, Cor_ECeg_model, Cor_ECa_anocova, Cor_Maas_ECe_model, Cor_Maas_ECeg_model),
    IOE = c(IOE_ECa_model, IOE_ECe_model, IOE_ECeg_model, IOE_ECa_anocova, IOE_Maas_ECe_model, IOE_Maas_ECeg_model)
  )
  print(Summary)
  Summary1 = data.frame(
    Model = c('ECa_Linear_Relative', 'ECe_Linear_Relative', 'ECeg_Linear_Relative', 'ECa_ANOCOVA_Relative', 'Maas_Relative', 'ECeg_Maas_Relative'),
    NRMSE = c(NRMSE_Rel_ECa_model, NRMSE_Rel_ECe_model, NRMSE_Rel_ECeg_model, NRMSE_Rel_ECa_anocova, NRMSE_Rel_Maas_ECe_model, NRMSE_Rel_Maas_ECeg_model),
    RMSE = c(RMSE_Rel_ECa_model, RMSE_Rel_ECe_model, RMSE_Rel_ECeg_model, RMSE_Rel_ECa_anocova, RMSE_Rel_Maas_ECe_model, RMSE_Rel_Maas_ECeg_model),
    Corr = c(Cor_Rel_ECa_model, Cor_Rel_ECe_model, Cor_Rel_ECeg_model, Cor_Rel_ECa_anocova, Cor_Rel_Maas_ECe_model, Cor_Rel_Maas_ECeg_model),
    IOE = c(IOE_Rel_ECa_model, IOE_Rel_ECe_model, IOE_Rel_ECeg_model, IOE_Rel_ECa_anocova,  IOE_Rel_Maas_ECe_model, IOE_Rel_Maas_ECeg_model)
  )
  print(Summary1)
  return(df)
}
max_yields_all = max_yield_df(corn_yield, avgtop = 1) # running this again here helps for some reason; else it gives error
corn_yield=linear_fit(df = corn_yield, test_df = corn_yield, max = max_yields_all, full_data = T)
obs_df=linear_fit(df = obs_df, test_df = obs_df, max = max_yields_obs, full_data = T)
corn_yield_cal=linear_fit(df = corn_yield_cal, test_df = corn_yield_test, max = max_yields_cal,  full_data = F)
#Get RMSPE from limited subset to compare:
corn_yield_cal=linear_fit(df = corn_yield_cal, test_df = corn_yield_test, max = max_yields_cal,  full_data = T)

#Calculation of Yield and Relative Yield% Sigmoidal Curve
#max variable is the max_yields df from above
sigmoid_fit = function(df, test_df, full_data, max) {
  #sigmoidal 4pl model generation
  #Equation form: Yr = d+(a-d)/(1+(C/c)^b)
  # a = upper Y-intercept
  # b = curve steepness?
  # c = mid slope?
  # d = lower Y-intercept

  #sigmoid_start_eca = list(a=100, b=4, c=125, d=25) #works with 38 and 20 sample size
  sigmoid_start_eca = list(a=91, b=4.6, c=116, d=25) #works with 12 sample size
  eca_sig_model <<- nlsLM(relative_yield~I(d+(a-d)/(1+(CV.1.0m/c)^b)), data=df, start=sigmoid_start_eca) #ECa, Yr
    df$eca_sig_Yr <- predict(eca_sig_model, data = df)
    test_df$eca_sig_Yr_test <- predict(eca_sig_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$eca_sig_Yr
       df$eca_sig_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$eca_sig_Yr_test
       test_df$eca_sig_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(eca_sig_model))
  #print(confint(eca_sig_model)) #use for nlrob> , method = "Wald"))

  sigmoid_start_ece = list(a=85, b=4, c=5.19, d=9.5)
  ece_sig_model <<- nlsLM(relative_yield~I(d+(a-d)/(1+(ECe/c)^b)), data=df, start=sigmoid_start_ece) #ECe, Yr
    df$ece_sig_Yr <- predict(ece_sig_model, data = df)
    test_df$ece_sig_Yr_test <- predict(ece_sig_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$ece_sig_Yr
       df$ece_sig_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$ece_sig_Yr_test
       test_df$ece_sig_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(ece_sig_model))
  #print(confint(ece_sig_model)) #use for nlrob> , method = "Wald"))

  sigmoid_start_eceg = list(a=80, b=4, c=4, d=12)
  #sigmoid_start_eceg = list(a=90, b=4, c=4) #test with custom "d"
  eceg_sig_model <<- nlsLM(relative_yield~I(d+(a-d)/(1+(ECeg/c)^b)), data=df, start=sigmoid_start_eceg) #ECe, Yr
  #eceg_sig_model <<- nlsLM(relative_yield~I(12+(a-12)/(1+(ECeg/c)^b)), data=df, start=sigmoid_start_eceg) #ECe, Yr
    df$eceg_sig_Yr <- predict(eceg_sig_model, data = df)
    test_df$eceg_sig_Yr_test <- predict(eceg_sig_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$eceg_sig_Yr
       df$eceg_sig_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$eceg_sig_Yr_test
       test_df$eceg_sig_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(eceg_sig_model))
  #print(confint(eceg_sig_model)) #use for nlrob> , method = "Wald"))

  #van genucthen relative yield models - ECe and ECa
  #Equation form: Yr = nlsLM(Yr ~ Ym/(1+(C/C_50)^exp(s*C_50))
  genuchten_start_eca = list(s=0.00536691, p=mean(max$half_ECa))
  #half_ECa = mean(max$half_ECa)
  #half_ECa = median(max$half_ECa)
  eca_genuchten_model <<- nlsLM(relative_yield~I(100/(1+(CV.1.0m/p)^exp(s*p))), data=df, start = genuchten_start_eca)
  #eca_genuchten_model <<- nlrob(relative_yield~I(100/(1+(CV.1.0m/p)^exp(s*p))), data=df, start = genuchten_start_eca, method = 'mtl')
    df$eca_gen_Yr <- predict(eca_genuchten_model, data = df)
    test_df$eca_gen_Yr_test <- predict(eca_genuchten_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$eca_gen_Yr
       df$eca_gen_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$eca_gen_Yr_test
       test_df$eca_gen_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(eca_genuchten_model))
  print(confint(eca_genuchten_model))
  half_ECa <<- coef(eca_genuchten_model)['p']
  s_ECa <<- coef(eca_genuchten_model)['s']

  genuchten_start_ece = list(s=0.01, p=mean(max$half_ECe))
  #half_ECe = mean(max$half_ECe)
  #half_ECe = median(max$half_ECe)
  #half_ECe = 5.867 #maas and hoffman derived c50 value. worked best on all except 6 sample design.
  ece_genuchten_model <<- nlsLM(relative_yield~I(100/(1+(ECe/p)^exp(s*p))), data=df, start = genuchten_start_ece)
  #ece_genuchten_model <<- nlrob(relative_yield~I(100/(1+(ECe/p)^exp(s*p))), data=df, start = genuchten_start_ece, method = 'mtl')
    df$ece_gen_Yr <- predict(ece_genuchten_model, data = df)
    test_df$ece_gen_Yr_test <- predict(ece_genuchten_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$ece_gen_Yr
       df$ece_gen_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$ece_gen_Yr_test
       test_df$ece_gen_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(ece_genuchten_model))
  print(confint(ece_genuchten_model))
  half_ECe <<- coef(ece_genuchten_model)['p']
  s_ECe <<- coef(ece_genuchten_model)['s']

  #genuchten_start_eceg = list(s=0.01, p=4.5)
  genuchten_start_eceg = list(s=0.01, p=mean(max$half_ECe))
  #half_ECe = mean(max$half_ECe)
  #half_ECeg = median(max$half_ECeg)
  #half_ECe = 5.867 #maas and hoffman derived c50 value. worked best on all except 6 sample design.
  eceg_genuchten_model <<- nlsLM(relative_yield~I(100/(1+(ECeg/p)^exp(s*p))), data=df, start = genuchten_start_eceg)
  # eceg_genuchten_model <<- nlrob(relative_yield~I(100/(1+(ECeg/p)^exp(s*p))), data=df,
  #                                method = 'M')
    df$eceg_gen_Yr <- predict(eceg_genuchten_model, data = df)
    test_df$eceg_gen_Yr_test <- predict(eceg_genuchten_model, newdata = test_df)
     for(i in 1:length(df$ID)) {
       x = df$Field
       y = df$eceg_gen_Yr
       df$eceg_gen_Y[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
     for(i in 1:length(test_df$ID)) {
       x = test_df$Field
       y = test_df$eceg_gen_Yr_test
       test_df$eceg_gen_Y_test[i] <- ((y[i]/100)*max$max_yield[max$Field == x[i]])
     }
  print(summary(eceg_genuchten_model))
  print(confint(eceg_genuchten_model)) #use for nlrob> , method = "Wald"))
  half_ECeg <<- coef(eceg_genuchten_model)['p']
  s_ECeg <<- coef(eceg_genuchten_model)['s']

  #Goodness of fit statistics
  mse_sig_Y_ECa = mse(df$eca_sig_Y, df$kg.ha)  #Y from ECa - 4pl
  mse_sig_Y_ECe = mse(df$ece_sig_Y, df$kg.ha)  #Y from ECe - 4pl
  mse_sig_Y_ECeg = mse(df$eceg_sig_Y, df$kg.ha)  #Y from ECeg - 4pl
  mse_sig_Yr_ECa = mse(df$eca_sig_Yr, df$relative_yield) #Yr from ECa - 4pl
  mse_sig_Yr_ECe = mse(df$ece_sig_Yr, df$relative_yield) #Yr from ECe - 4pl
  mse_sig_Yr_ECeg = mse(df$eceg_sig_Yr, df$relative_yield) #Yr from ECeg - 4pl

  mse5 = mse(df$ece_gen_Y, df$kg.ha)  #Y from ECe - genuchten
  mse6 = mse(df$eca_gen_Y, df$kg.ha)  #Y from ECa - genuchten
  mse_g = mse(df$eceg_gen_Y, df$kg.ha)  #Y from ECeg - genuchten
  mse7 = mse(df$ece_gen_Yr, df$relative_yield) #Yr from ECe - genuchten
  mse8 = mse(df$eca_gen_Yr, df$relative_yield) #Yr from ECa - genuchten
  mse_gr = mse(df$eceg_gen_Yr, df$relative_yield) #Yr from ECeg - genuchten

  rmse_sig_Y_ECa = rmse(df$eca_sig_Y, df$kg.ha)  #Y from ECa - 4pl
  rmse_sig_Y_ECe = rmse(df$ece_sig_Y, df$kg.ha)  #Y from ECe - 4pl
  rmse_sig_Y_ECeg = rmse(df$eceg_sig_Y, df$kg.ha)  #Y from ECeg - 4pl
  rmse_sig_Yr_ECa = rmse(df$eca_sig_Yr, df$relative_yield) #Yr from ECa - 4pl
  rmse_sig_Yr_ECe = rmse(df$ece_sig_Yr, df$relative_yield) #Yr from ECe - 4pl
  rmse_sig_Yr_ECeg = rmse(df$eceg_sig_Yr, df$relative_yield) #Yr from ECeg - 4pl
  rmse5 = rmse(df$ece_gen_Y, df$kg.ha)
  rmse6 = rmse(df$eca_gen_Y, df$kg.ha)
  rmse_g = rmse(df$eceg_gen_Y, df$kg.ha)  #Y from ECeg - genuchten
  rmse7 = rmse(df$ece_gen_Yr, df$relative_yield)
  rmse8 = rmse(df$eca_gen_Yr, df$relative_yield)
  rmse_gr = rmse(df$eceg_gen_Yr, df$relative_yield) #Yr from ECeg - genuchten

  nrmse_sig_Y_ECa = nrmse(df$eca_sig_Y, df$kg.ha)  #Y from ECa - 4pl
  nrmse_sig_Y_ECe = nrmse(df$ece_sig_Y, df$kg.ha)  #Y from ECe - 4pl
  nrmse_sig_Y_ECeg = nrmse(df$eceg_sig_Y, df$kg.ha)  #Y from ECeg - 4pl
  nrmse_sig_Yr_ECa = nrmse(df$eca_sig_Yr, df$relative_yield) #Yr from ECa - 4pl
  nrmse_sig_Yr_ECe = nrmse(df$ece_sig_Yr, df$relative_yield) #Yr from ECe - 4pl
  nrmse_sig_Yr_ECeg = nrmse(df$eceg_sig_Yr, df$relative_yield) #Yr from ECeg - 4pl
  nrmse5 = nrmse(df$ece_gen_Y, df$kg.ha)
  nrmse6 = nrmse(df$eca_gen_Y, df$kg.ha)
  nrmse_g = nrmse(df$eceg_gen_Y, df$kg.ha)  #Y from ECeg - genuchten
  nrmse7 = nrmse(df$ece_gen_Yr, df$relative_yield)
  nrmse8 = nrmse(df$eca_gen_Yr, df$relative_yield)
  nrmse_gr = nrmse(df$eceg_gen_Yr, df$relative_yield) #Yr from ECeg - genuchten

  cor_sig_Y_ECa = cor(df$eca_sig_Y, df$kg.ha)  #Y from ECa - 4pl
  cor_sig_Y_ECe = cor(df$ece_sig_Y, df$kg.ha)  #Y from ECe - 4pl
  cor_sig_Y_ECeg = cor(df$eceg_sig_Y, df$kg.ha)  #Y from ECeg - 4pl
  cor_sig_Yr_ECa = cor(df$eca_sig_Yr, df$relative_yield) #Yr from ECa - 4pl
  cor_sig_Yr_ECe = cor(df$ece_sig_Yr, df$relative_yield) #Yr from ECe - 4pl
  cor_sig_Yr_ECeg = cor(df$eceg_sig_Yr, df$relative_yield) #Yr from ECeg - 4pl
  cor5 = cor(df$ece_gen_Y, df$kg.ha)
  cor6 = cor(df$eca_gen_Y, df$kg.ha)
  cor_g = cor(df$eceg_gen_Y, df$kg.ha)  #Y from ECeg - genuchten
  cor7 = cor(df$ece_gen_Yr, df$relative_yield)
  cor8 = cor(df$eca_gen_Yr, df$relative_yield)
  cor_gr = cor(df$eceg_gen_Yr, df$relative_yield) #Yr from ECeg - genuchten

  d_sig_Y_ECa = d(df$eca_sig_Y, df$kg.ha)  #Y from ECa - 4pl
  d_sig_Y_ECe = d(df$ece_sig_Y, df$kg.ha)  #Y from ECe - 4pl
  d_sig_Y_ECeg = d(df$eceg_sig_Y, df$kg.ha)  #Y from ECeg - 4pl
  d_sig_Yr_ECa = d(df$eca_sig_Yr, df$relative_yield) #Yr from ECa - 4pl
  d_sig_Yr_ECe = d(df$ece_sig_Yr, df$relative_yield) #Yr from ECe - 4pl
  d_sig_Yr_ECeg = d(df$eceg_sig_Yr, df$relative_yield) #Yr from ECeg - 4pl
  d5 = d(df$ece_gen_Y, df$kg.ha)
  d6 = d(df$eca_gen_Y, df$kg.ha)
  d_g = d(df$eceg_gen_Y, df$kg.ha)  #Y from ECeg - genuchten
  d7 = d(df$ece_gen_Yr, df$relative_yield)
  d8 = d(df$eca_gen_Yr, df$relative_yield)
  d_gr = d(df$eceg_gen_Yr, df$relative_yield) #Yr from ECeg - genuchten

  if (full_data == T) {
    #for Yr sigmoidal 4pl ECa
     rmspe_list = c()
     field_list = c()
     for(i in unique(df$Field)) {
       test = subset(df, Field != i)
       lofo = subset(df, Field == i)
       x = test$CV.1.0m
       y = test$relative_yield
       mdl = nlsLM(y~d+(a-d)/(1+(x/c)^b), start=sigmoid_start_eca, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
       ydf = data.frame(x=lofo$CV.1.0m, y=lofo$relative_yield)
       y_pred = predict(mdl, newdata = ydf)
       rmspe = sqrt(mean((ydf$y - y_pred)^2))
       rmspe_list = append(rmspe_list, rmspe)
       field_list = append(field_list, i)
     }
       final_sig_df = data.frame(Left_out_Field = field_list,
                             ECa_Sig_RMSPE = rmspe_list)

     #for Yr sigmoidal 4pl ECe
     rmspe_list = c()
     field_list = c()
     for(i in unique(df$Field)){
       test = subset(df, Field != i)
       lofo = subset(df, Field == i)
       x = test$ECe
       y = test$relative_yield
       mdl = nlsLM(y~d+(a-d)/(1+(x/c)^b), start=sigmoid_start_ece, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
       ydf = data.frame(x=lofo$ECe, y=lofo$relative_yield)
       y_pred = predict(mdl, newdata = ydf)
       rmspe = sqrt(mean((ydf$y - y_pred)^2))
       rmspe_list = append(rmspe_list, rmspe)
       field_list = append(field_list, i)
       }
       final_sig_df1 = data.frame(Left_out_Field = field_list,
                             ECe_Sig_RMSPE = rmspe_list)

     #for Yr sigmoidal 4pl ECeg
     rmspe_list = c()
     field_list = c()
     for(i in unique(df$Field)){
       test = subset(df, Field != i)
       lofo = subset(df, Field == i)
       x = test$ECeg
       y = test$relative_yield
       mdl = nlsLM(y~d+(a-d)/(1+(x/c)^b), start=sigmoid_start_eceg, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
       #mdl = nlsLM(y~12+(a-12)/(1+(x/c)^b), start=sigmoid_start_eceg, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
       ydf = data.frame(x=lofo$ECeg, y=lofo$relative_yield)
       y_pred = predict(mdl, newdata = ydf)
       rmspe = sqrt(mean((ydf$y - y_pred)^2))
       rmspe_list = append(rmspe_list, rmspe)
       field_list = append(field_list, i)
       }
       final_sig_df2 = data.frame(Left_out_Field = field_list,
                             ECeg_Sig_RMSPE = rmspe_list)

    #for Yr van genuchten ECa
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$CV.1.0m
      y = test$relative_yield
      mdl = nlsLM(y~100/(1+(x/p)^exp(s*p)), start = genuchten_start_eca, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
      ydf = data.frame(x=lofo$CV.1.0m, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df = data.frame(Left_out_Field = field_list,
                            ECa_Gen_RMSPE = rmspe_list)
    #print(final_df2)
    #for Yr van genuchten ECe
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$ECe
      y = test$relative_yield
      mdl = nlsLM(y~100/(1+(x/p)^exp(s*p)), start = genuchten_start_ece, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
      ydf = data.frame(x=lofo$ECe, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df1 = data.frame(Left_out_Field = field_list,
                            ECe_Gen_RMSPE = rmspe_list)

    #for Yr van genuchten ECeg
    rmspe_list = c()
    field_list = c()
    for(i in unique(df$Field)){
      test = subset(df, Field != i)
      lofo = subset(df, Field == i)
      x = test$ECeg
      y = test$relative_yield
      mdl = nlsLM(y~100/(1+(x/p)^exp(s*p)), start = genuchten_start_eceg, control = nls.lm.control(maxiter = 1024, maxfev = 1024))
      ydf = data.frame(x=lofo$ECeg, y=lofo$relative_yield)
      y_pred = predict(mdl, newdata = ydf)
      rmspe = sqrt(mean((ydf$y - y_pred)^2))
      rmspe_list = append(rmspe_list, rmspe)
      field_list = append(field_list, i)
      }
      final_df2 = data.frame(Left_out_Field = field_list,
                            ECeg_Gen_RMSPE = rmspe_list)
    #print(final_df3)

    master_df <- Reduce(function(x, y) merge(x, y, by='Left_out_Field', all=TRUE), list(final_df, final_df1, final_df2, final_sig_df, final_sig_df1, final_sig_df2))
    mergedData <- merge(master_df, max,
                        by.x=c('Left_out_Field'), by.y=c('Field'))
    mergedData$Gen_ECa_RMSPE_Y = (mergedData$ECa_Gen_RMSPE/100)*mergedData$max_yield
    mergedData$Gen_ECe_RMSPE_Y = (mergedData$ECe_Gen_RMSPE/100)*mergedData$max_yield
    mergedData$Gen_ECeg_RMSPE_Y = (mergedData$ECeg_Gen_RMSPE/100)*mergedData$max_yield
    mergedData$Sig_ECa_RMSPE_Y = (mergedData$ECa_Sig_RMSPE/100)*mergedData$max_yield
    mergedData$Sig_ECe_RMSPE_Y = (mergedData$ECe_Sig_RMSPE/100)*mergedData$max_yield
    mergedData$Sig_ECeg_RMSPE_Y = (mergedData$ECeg_Sig_RMSPE/100)*mergedData$max_yield
    keeps <- c('Left_out_Field',
               'ECa_Gen_RMSPE',
               'ECe_Gen_RMSPE',
               'ECeg_Gen_RMSPE',
               'ECa_Sig_RMSPE',
               'ECe_Sig_RMSPE',
               'ECeg_Sig_RMSPE',
               'Gen_ECa_RMSPE_Y',
               'Gen_ECe_RMSPE_Y',
               'Gen_ECeg_RMSPE_Y',
               'Sig_ECa_RMSPE_Y',
               'Sig_ECe_RMSPE_Y',
               'Sig_ECeg_RMSPE_Y')
    print(mergedData[keeps])
  }
  else {
    #Goodness of fit statistics
    mse_test_sig_Y_ECa = mse(test_df$eca_sig_Y, test_df$kg.ha)  #Y from ECa - 4pl
    mse_test_sig_Y_ECe = mse(test_df$ece_sig_Y, test_df$kg.ha)  #Y from ECe - 4pl
    mse_test_sig_Y_ECeg = mse(test_df$eceg_sig_Y, test_df$kg.ha)  #Y from ECeg - 4pl
    mse_test_sig_Yr_ECa = mse(test_df$eca_sig_Yr, test_df$relative_yield) #Yr from ECa - 4pl
    mse_test_sig_Yr_ECe = mse(test_df$ece_sig_Yr, test_df$relative_yield) #Yr from ECe - 4pl
    mse_test_sig_Yr_ECeg = mse(test_df$eceg_sig_Yr, test_df$relative_yield) #Yr from ECeg - 4pl
    mse5_test = mse(test_df$ece_gen_Y_test, test_df$kg.ha)  #Y from ECe - genuchten
    mse6_test = mse(test_df$eca_gen_Y_test, test_df$kg.ha)  #Y from ECa - genuchten
    mse_g_test = mse(test_df$eceg_gen_Y_test, test_df$kg.ha)  #Y from ECeg - genuchten
    mse7_test = mse(as.numeric(test_df$ece_gen_Yr_test), test_df$relative_yield) #Yr from ECe - genuchten
    mse8_test = mse(as.numeric(test_df$eca_gen_Yr_test), test_df$relative_yield) #Yr from ECa - genuchten
    mse_gr_test = mse(as.numeric(test_df$eceg_gen_Yr_test), test_df$relative_yield) #Yr from ECeg - genuchten

    rmse_test_sig_Y_ECa = rmse(test_df$eca_sig_Y, test_df$kg.ha)  #Y from ECa - 4pl
    rmse_test_sig_Y_ECe = rmse(test_df$ece_sig_Y, test_df$kg.ha)  #Y from ECe - 4pl
    rmse_test_sig_Y_ECeg = rmse(test_df$eceg_sig_Y, test_df$kg.ha)  #Y from ECeg - 4pl
    rmse_test_sig_Yr_ECa = rmse(test_df$eca_sig_Yr, test_df$relative_yield) #Yr from ECa - 4pl
    rmse_test_sig_Yr_ECe = rmse(test_df$ece_sig_Yr, test_df$relative_yield) #Yr from ECe - 4pl
    rmse_test_sig_Yr_ECeg = rmse(test_df$eceg_sig_Yr, test_df$relative_yield) #Yr from ECeg - 4pl
    rmse5_test = rmse(test_df$ece_gen_Y_test, test_df$kg.ha)
    rmse6_test = rmse(test_df$eca_gen_Y_test, test_df$kg.ha)
    rmse_g_test = rmse(test_df$eceg_gen_Y_test, test_df$kg.ha)  #Y from ECeg - genuchten
    rmse7_test = rmse(as.numeric(test_df$ece_gen_Yr_test), test_df$relative_yield)
    rmse8_test = rmse(as.numeric(test_df$eca_gen_Yr_test), test_df$relative_yield)
    rmse_gr_test = rmse(as.numeric(test_df$eceg_gen_Yr_test), test_df$relative_yield) #Yr from ECeg - genuchten

    nrmse_test_sig_Y_ECa = nrmse(test_df$eca_sig_Y, test_df$kg.ha)  #Y from ECa - 4pl
    nrmse_test_sig_Y_ECe = nrmse(test_df$ece_sig_Y, test_df$kg.ha)  #Y from ECe - 4pl
    nrmse_test_sig_Y_ECeg = nrmse(test_df$eceg_sig_Y, test_df$kg.ha)  #Y from ECeg - 4pl
    nrmse_test_sig_Yr_ECa = nrmse(test_df$eca_sig_Yr, test_df$relative_yield) #Yr from ECa - 4pl
    nrmse_test_sig_Yr_ECe = nrmse(test_df$ece_sig_Yr, test_df$relative_yield) #Yr from ECe - 4pl
    nrmse_test_sig_Yr_ECeg = nrmse(test_df$eceg_sig_Yr, test_df$relative_yield) #Yr from ECeg - 4pl
    nrmse5_test = nrmse(test_df$ece_gen_Y_test, test_df$kg.ha)
    nrmse6_test = nrmse(test_df$eca_gen_Y_test, test_df$kg.ha)
    nrmse_g_test = nrmse(test_df$eceg_gen_Y_test, test_df$kg.ha)  #Y from ECeg - genuchten
    nrmse7_test = nrmse(as.numeric(test_df$ece_gen_Yr_test), test_df$relative_yield)
    nrmse8_test = nrmse(as.numeric(test_df$eca_gen_Yr_test), test_df$relative_yield)
    nrmse_gr_test = nrmse(as.numeric(test_df$eceg_gen_Yr_test), test_df$relative_yield) #Yr from ECeg - genuchten

    cor_test_sig_Y_ECa = cor(test_df$eca_sig_Y, test_df$kg.ha)  #Y from ECa - 4pl
    cor_test_sig_Y_ECe = cor(test_df$ece_sig_Y, test_df$kg.ha)  #Y from ECe - 4pl
    cor_test_sig_Y_ECeg = cor(test_df$eceg_sig_Y, test_df$kg.ha)  #Y from ECeg - 4pl
    cor_test_sig_Yr_ECa = cor(test_df$eca_sig_Yr, test_df$relative_yield) #Yr from ECa - 4pl
    cor_test_sig_Yr_ECe = cor(test_df$ece_sig_Yr, test_df$relative_yield) #Yr from ECe - 4pl
    cor_test_sig_Yr_ECeg = cor(test_df$eceg_sig_Yr, test_df$relative_yield) #Yr from ECeg - 4pl
    cor5_test = cor(test_df$ece_gen_Y_test, test_df$kg.ha)
    cor6_test = cor(test_df$eca_gen_Y_test, test_df$kg.ha)
    cor_g_test = cor(test_df$eceg_gen_Y_test, test_df$kg.ha)  #Y from ECeg - genuchten
    cor7_test = cor(as.numeric(test_df$ece_gen_Yr_test), test_df$relative_yield)
    cor8_test = cor(as.numeric(test_df$eca_gen_Yr_test), test_df$relative_yield)
    cor_gr_test = cor(as.numeric(test_df$eceg_gen_Yr_test), test_df$relative_yield) #Yr from ECeg - genuchten

    d_test_sig_Y_ECa = d(test_df$eca_sig_Y, test_df$kg.ha)  #Y from ECa - 4pl
    d_test_sig_Y_ECe = d(test_df$ece_sig_Y, test_df$kg.ha)  #Y from ECe - 4pl
    d_test_sig_Y_ECeg = d(test_df$eceg_sig_Y, test_df$kg.ha)  #Y from ECeg - 4pl
    d_test_sig_Yr_ECa = d(test_df$eca_sig_Yr, test_df$relative_yield) #Yr from ECa - 4pl
    d_test_sig_Yr_ECe = d(test_df$ece_sig_Yr, test_df$relative_yield) #Yr from ECe - 4pl
    d_test_sig_Yr_ECeg = d(test_df$eceg_sig_Yr, test_df$relative_yield) #Yr from ECeg - 4pl
    d5_test = d(test_df$ece_gen_Y_test, test_df$kg.ha)
    d6_test = d(test_df$eca_gen_Y_test, test_df$kg.ha)
    d_g_test = d(test_df$eceg_gen_Y_test, test_df$kg.ha)  #Y from ECeg - genuchten
    d7_test = d(as.numeric(test_df$ece_gen_Yr_test), test_df$relative_yield)
    d8_test = d(as.numeric(test_df$eca_gen_Yr_test), test_df$relative_yield)
    d_gr_test = d(as.numeric(test_df$eceg_gen_Yr_test), test_df$relative_yield) #Yr from ECeg - genuchten

    print("test dataset statistics:###########################################")
    bu_kg = 1 #a value of 1 here means units are kg/ha
    #summary of absolute yield predictions
    Summary = data.frame(
      Model = c('ECa_Genuchten_Yield', 'ECe_Genuchten_Yield', 'ECeg_Genuchten_Yield', 'ECa_Sigmoid_Yield', 'ECe_Sigmoid_Yield', 'ECeg_Sigmoid_Yield'),
      NRMSE = c(nrmse6_test, nrmse5_test, nrmse_g_test, nrmse_test_sig_Y_ECa, nrmse_test_sig_Y_ECe, nrmse_test_sig_Y_ECeg),
      MSE = c(mse6_test, mse5_test, mse_g_test, mse_test_sig_Y_ECa, mse_test_sig_Y_ECe, mse_test_sig_Y_ECeg),
      RMSE = c(rmse6_test*bu_kg, rmse5_test*bu_kg, rmse_g_test*bu_kg, rmse_test_sig_Y_ECa*bu_kg, rmse_test_sig_Y_ECe*bu_kg, rmse_test_sig_Y_ECeg*bu_kg),
      Corr = c(cor6_test, cor5_test, cor_g_test, cor_test_sig_Y_ECa, cor_test_sig_Y_ECe, cor_test_sig_Y_ECeg),
      IOA = c(d6_test, d5_test, d_g_test, d_test_sig_Y_ECa, d_test_sig_Y_ECe, d_test_sig_Y_ECeg)
    )
    print(Summary)
    #summary of relative yield predictions
    Summary1 = data.frame(
      Model = c('ECa_Genuchten_Relative', 'ECe_Genuchten_Relative', 'ECeg_Genuchten_Relative', 'ECa_Sigmoid_Relative', 'ECe_Sigmoid_Relative', 'ECeg_Sigmoid_Relative'),
      NRMSE = c(nrmse8_test, nrmse7_test, nrmse_gr_test, nrmse_test_sig_Yr_ECa, nrmse_test_sig_Yr_ECe, nrmse_test_sig_Yr_ECeg),
      MSE = c(mse8_test, mse7_test, mse_gr_test, mse_test_sig_Yr_ECa, mse_test_sig_Yr_ECe, mse_test_sig_Yr_ECeg),
      RMSE = c(rmse8_test, rmse7_test, rmse_gr_test*bu_kg, rmse_test_sig_Yr_ECa*bu_kg, rmse_test_sig_Yr_ECe*bu_kg, rmse_test_sig_Yr_ECeg*bu_kg),
      Corr = c(cor8_test, cor7_test, cor_gr_test, cor_test_sig_Yr_ECa, cor_test_sig_Yr_ECe, cor_test_sig_Yr_ECeg),
      IOA = c(d8_test, d7_test, d_gr_test, d_test_sig_Yr_ECa, d_test_sig_Yr_ECe, d_test_sig_Yr_ECeg)
    )
    print(Summary1)
  }
  print("full dataset statistics:###########################################")
  bu_kg = 1 #a value of 1 here means units are kg/ha
  #summary of absolute yield predictions
  Summary = data.frame(
    Model = c('ECa_Genuchten_Yield', 'ECe_Genuchten_Yield', 'ECeg_Genuchten_Yield', 'ECa_Sigmoid_Yield', 'ECe_Sigmoid_Yield', 'ECeg_Sigmoid_Yield'),
    NRMSE = c(nrmse6, nrmse5, nrmse_g, nrmse_sig_Y_ECa, nrmse_sig_Y_ECe, nrmse_sig_Y_ECeg),
    MSE = c(mse6, mse5, mse_g, mse_sig_Y_ECa, mse_sig_Y_ECe, mse_sig_Y_ECeg),
    RMSE = c(rmse6*bu_kg, rmse5*bu_kg, rmse_g*bu_kg, rmse_sig_Y_ECa*bu_kg, rmse_sig_Y_ECe*bu_kg, rmse_sig_Y_ECeg*bu_kg),
    Corr = c(cor6, cor5, cor_g, cor_sig_Y_ECa, cor_sig_Y_ECe, cor_sig_Y_ECeg),
    IOA = c(d6, d5, d_g, d_sig_Y_ECa, d_sig_Y_ECe, d_sig_Y_ECeg)
  )
  print(Summary)

  #summary of relative yield predictions
  Summary1 = data.frame(
    Model = c('ECa_Genuchten_Relative', 'ECe_Genuchten_Relative', 'ECeg_Genuchten_Relative', 'ECa_Sigmoid_Relative', 'ECe_Sigmoid_Relative', 'ECeg_Sigmoid_Relative'),
    NRMSE = c(nrmse8, nrmse7, nrmse_gr, nrmse_sig_Yr_ECa, nrmse_sig_Yr_ECe, nrmse_sig_Yr_ECeg),
    MSE = c(mse8, mse7, mse_gr, mse_sig_Yr_ECa, mse_sig_Yr_ECe, mse_sig_Yr_ECeg),
    RMSE = c(rmse8, rmse7, rmse_gr, rmse_sig_Yr_ECa*bu_kg, rmse_sig_Yr_ECe*bu_kg, rmse_sig_Yr_ECeg*bu_kg),
    Corr = c(cor8, cor7, cor_gr, cor_sig_Yr_ECa, cor_sig_Yr_ECe, cor_sig_Yr_ECeg),
    IOA = c(d8, d7, d_gr, d_sig_Yr_ECa, d_sig_Yr_ECe, d_sig_Yr_ECeg)
  )
  print(Summary1)
  return(df)
}
#Fit models using all predicted EC data
corn_yield=sigmoid_fit(df = corn_yield, test_df = corn_yield, max = max_yields_all, full_data = T)
#Fit model using observed EC data only
obs_df=sigmoid_fit(df = obs_df, test_df = obs_df, max = max_yields_all, full_data = T)
#Use left out data to get test GOF stats
corn_yield_cal=sigmoid_fit(df = corn_yield_cal, test_df = corn_yield_test, max = max_yields_cal,  full_data = F)
#Get RMSPE from limited subset to compare:
corn_yield_cal=sigmoid_fit(df = corn_yield_cal, test_df = corn_yield_test, max = max_yields_cal,  full_data = T)

kfold_gen = function(df, form, start, folds, reps) {
  modelInfo <- list(label = "Modified Discount Function",
                   library = "minpack.lm",
                   type = "Regression",
                   parameters = data.frame(parameter = "parameter",
                                           class = "character",
                                           label = "parameter"),
                   grid = function(x, y, len = NULL, search = "grid")
                     data.frame(parameter = "none"),
                   loop = NULL,
                   fit = function(x, y, data=df, ...) {
                     genuchten_start = start
                     eca_genuchten_model <- nlsLM(
                       form,
                       data=data,
                       start = genuchten_start)
                   },
                   predict = function(modelFit, newdata, submodels = NULL) {
                     predict(modelFit, newdata)
                   },
                   prob = NULL,
                   predictors = function(x, ...) {
                     unique(as.vector(variable.names(x)))
                   },
                   levels = NULL,
                   sort = function(x) x)
  set.seed(307)
  mboost_resamp <- train(x=df,
                        y=df$relative_yield,
                        method = modelInfo,
                        trControl = trainControl(method = "repeatedcv",
                                                 repeats = reps,
                                                 number = folds,
                                                 allowParallel = T)
                        )

  print(mboost_resamp$results)
  print(mboost_resamp$finalModel)
  # a = histogram(mboost_resamp, metric = 'RMSE')
  # b = densityplot(mboost_resamp, metric = 'RMSE')
  # print(a, position = c(0, 0, 0.5, 1), more = TRUE)
  # print(b, position = c(0.5, 0, 1, 1))
  #densityplot(mboost_resamp, pch = "|")
}
kfold_gen(df=corn_yield,
          form = formula(eca_genuchten_model),
          start = list(s=s_ECa, p=half_ECa),
          reps = 50,
          folds = 10)
kfold_gen(df=corn_yield,
          form = formula(ece_genuchten_model),
          start = list(s=s_ECe, p=half_ECe),
          reps = 50,
          folds = 10)
kfold_gen(df=corn_yield,
          form = formula(eceg_genuchten_model),
          start = list(s=s_ECeg, p=half_ECeg),
          reps = 50,
          folds = 10)
kfold_gen(df=corn_yield,
          form = formula(eca_sig_model),
          start = list(a=85, b=4, c=70, d=25),
          reps = 50,
          folds = 10)
kfold_gen(df=corn_yield,
          form = formula(ece_sig_model),
          start = list(a=100, b=4, c=4, d=25),
          reps = 50,
          folds = 10)
kfold_gen(df=corn_yield,
          form = formula(eceg_sig_model),
          start = list(a=90, b=4, c=4, d=25),
          reps = 50,
          folds = 10)

#Graphing Results
summary_graph = function(df, max) {
  df$Year = as.factor(df$Year)
  half_ECa <- coef(eca_genuchten_model)['p']
  half_ECe <- coef(ece_genuchten_model)['p']
  half_ECeg <- coef(eceg_genuchten_model)['p']
  s_ECa <- coef(eca_genuchten_model)['s']
  s_ECe <- coef(ece_genuchten_model)['s']
  s_ECeg <- coef(eceg_genuchten_model)['s']
  a_ECa <- coef(eca_sig_model)['a']
  a_ECe <- coef(ece_sig_model)['a']
  a_ECeg <- coef(eceg_sig_model)['a']
  b_ECa <- coef(eca_sig_model)['b']
  b_ECe <- coef(ece_sig_model)['b']
  b_ECeg <- coef(eceg_sig_model)['b']
  c_ECa <- coef(eca_sig_model)['c']
  c_ECe <- coef(ece_sig_model)['c']
  c_ECeg <- coef(eceg_sig_model)['c']
  d_ECa <- coef(eca_sig_model)['d']
  d_ECe <- coef(ece_sig_model)['d']
  d_ECeg <- coef(eceg_sig_model)['d']
  #d_ECeg <- 12 #test to fit sigmoid better
  y = df$kg.ha
  y1 = df$relative_yield
  x = df$CV.1.0m
  x1 = df$ECe
  xg = df$ECeg
  z = df$Field
  start = list(s=0.01)
  #plots comparing Relative Yield %
  se_ece = 1.177697#exp(summary(deep_model)$coefficients[17]) #for slope std. error, should be 1.177697
  se_yld = sd(df$relative_yield)/sqrt(length(df$relative_yield))
  formatter100 = function(x){ # <- function to scale ECa to dS/m
    x/100
  }
    ECa = ggplot(df, aes(x, y1))+#, col=Field)) +
      #geom_point(col=df$Year) +
      geom_point(col='gray') +
      geom_point(data = obs_df, aes(CV.1.0m, relative_yield), shape=17, size=2) +
      stat_function(fun = function(.x) 100/(1+(.x/half_ECa)^exp(s_ECa*half_ECa)), color='black') +
      stat_function(fun = function(.x) d_ECa+(a_ECa-d_ECa)/(1+(.x/c_ECa)^b_ECa), color='red') +
      #geom_smooth(method='lm', se=F, fullrange=T, color='red') +
      xlab(expression(EC[a]~dS~m^{-1})) +
      scale_x_continuous(labels = formatter100, limits = c(0,300)) +
      ylab(expression(Y[r]~'%')) +
      ylim(0,105) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), text = element_text(size=25))


    ECe = ggplot(df, aes(x1, y1))+#, col=Field)) +
      #geom_point(aes(col=Year)) +
      geom_point(col='gray') +
      geom_point(data = obs_df, aes(ECe, relative_yield, color = "black"), shape=17, size=2) +
      stat_function(fun = function(x1) 100/(1+(x1/half_ECe)^exp(s_ECe*half_ECe)), aes(linetype = 'Modified Discount Function')) +
      stat_function(fun = function(x1) d_ECe+(a_ECe-d_ECe)/(1+(x1/c_ECe)^b_ECe), aes(linetype = 'Custom Sigmoidal Model', color = 'red')) +
      #geom_smooth(method='lm', se=F, fullrange=T, aes(color = 'Linear Regression')) +
      stat_function(fun = function(x) ifelse(x>=1.7, 100 - 12*(x - 1.7), 100), aes(linetype = 'Threshold-Slope Function'), size = 1) +
      scale_linetype_manual(name='',
                            breaks=c('Modified Discount Function', 'Custom Sigmoidal Model', 'Threshold-Slope Function'),
                            values = c('solid', 'solid', 'dashed')) +
      scale_color_identity(guide = "legend",
                           name='',
                           labels=c('Observed Points')) +
      xlab(expression(EC[e]~dS~m^{-1})) +
      xlim(0,15) +
      ylab(expression(Y[r]~'%')) +
      ylim(0,105) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position=c('bottom'),
            text=element_text(size=20), legend.key.size = unit(3,"line"))


    ECeg = ggplot(df, aes(ECeg, relative_yield))+#, col=Field)) +
      #geom_point(col=df$Year) +
      geom_point(col='gray') +
      geom_point(data = obs_df, aes(ECeg, relative_yield), shape=17, size=2) +
      stat_function(fun = function(xg) 100/(1+(xg/half_ECeg)^exp(s_ECeg*half_ECeg)), color='black') +
      stat_function(fun = function(xg) d_ECeg+(a_ECeg-d_ECeg)/(1+(xg/c_ECeg)^b_ECeg), color='red') +
      #geom_smooth(method='nls', method.args=list(formula = y~100/I((1+(x/p)^exp(s*p))),
      #            start=list(s=s_ECeg, p=half_ECeg)), color='black') +
      #geom_smooth(method='lm', se=F, fullrange=T, color = 'red') +
      stat_function(fun = function(x) ifelse(x>=1.7, 100 - 12*(x - 1.7), 100), size = 1, linetype='dashed') +
      xlab(expression(EC[eg]~dS~m^{-1})) +
      xlim(0,15) +
      ylab(expression(Y[r]~'%')) +
      ylim(0,105) + #disregard 24 point omission warning due to this line of code
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5), legend.position=c('bottom'),
            text=element_text(size=25))

    ggarrange(ECa,ECe,ECeg, nrow = 1, ncol = 3, common.legend = TRUE, legend="bottom")
}
#run the sigmoid fxn again to define variables for graph
corn_yield=sigmoid_fit(df = corn_yield, test_df = corn_yield, max = max_yields_all, full_data = T)
summary_graph(df= corn_yield, max=max_yields_all)
#rerun the sigmoid fxn again to define variables for obs graph
obs_df=sigmoid_fit(df = obs_df, test_df = obs_df, max = max_yields_all, full_data = T)
summary_graph(df= obs_df, max=max_yields_obs)

summary_graph(df=corn_yield_cal, max=max_yields_cal)

#Residual plots
plot(nlsResiduals(eca_genuchten_model))
plot(nlsResiduals(ece_genuchten_model))
plot(nlsResiduals(eceg_genuchten_model))

#Variance by EC category
variance_fxn = function(df, obs) {
  salt_df1 = subset(df, ECe <= 2)
  salt_df2 = subset(df, ECe <= 4)
  salt_df3 = subset(df, ECe <= 6)
  salt_df4 = subset(df, ECe <= 8)
  salt_df5 = subset(df, ECe > 8)
  var1 = var(salt_df1$relative_yield)
  var2 = var(salt_df2$relative_yield)
  var3 = var(salt_df3$relative_yield)
  var4 = var(salt_df4$relative_yield)
  var5 = var(salt_df5$relative_yield)
  salt_obs1 = subset(obs, ECe_Deep <= 2)
  salt_obs2 = subset(obs, ECe_Deep <= 4)
  salt_obs3 = subset(obs, ECe_Deep <= 6)
  salt_obs4 = subset(obs, ECe_Deep <= 8)
  salt_obs5 = subset(obs, ECe_Deep > 8)
  var11 = var(salt_obs1$relative_yield)
  var22 = var(salt_obs2$relative_yield)
  var33 = var(salt_obs3$relative_yield)
  var44 = var(salt_obs4$relative_yield)
  var55 = var(salt_obs5$relative_yield)
  Summary = data.frame(
    Salinity = c('0-2 dS/m','2-4 dS/m','4-6 dS/m','6-8 dS/m','>8 dS/m'),
    Yr_Variance = c(var1, var2, var3, var4, var5),
    Obs_Variance = c(var11, var22, var33, var44, var55)
  )
  print(Summary)
}
variance_fxn(corn_yield, obs_df)
corr.test(corn_yield$CV.1.0m, corn_yield$relative_yield)
corr.test(corn_yield$ECe, corn_yield$relative_yield)
corr.test(corn_yield$ECeg, corn_yield$relative_yield)

prediction_csv = function(df, name) {
  ################ WRITE OUT SIGMOID RESULTS TO .CSV
  id = df$ID
  dfx2 = data.frame(id,
                    'X'=df$X,
                    'Y'=df$Y,
                    'Y_obs'=df$kg.ha,
                    'Yr_obs'=df$relative_yield,
                    'ECa_obs'=df$CV.1.0m,
                    'ECe_pred'=df$ECe,
                    'Field_ID'=df$Field,
                    'Y_ANOCOVA_ECa'=df$anocova_eca_yield_pred_from_rel,
                    'Yr_ANOCOVA_ECa'=df$anocova_eca_relative_yield,
                    #"Y_sig_ECa"=df$eca_sig_Y,
                    #"Yr_sig_ECa"=df$eca_sig_Yr,
                    #"Y_sig_ECe"=df$ece_sig_Y,
                    #"Yr_sig_ECe"=df$ece_sig_Yr,
                    "Y_Gen_ECe"=df$ece_gen_Y,
                    "Y_Gen_ECa"=df$eca_gen_Y,
                    "Yr_Gen_ECe"=df$ece_gen_Yr,
                    "Yr_Gen_ECa"=df$eca_gen_Yr,
                    'Y_lin_ECa'=df$eca_yield_pred_from_rel,
                    'Y_lin_ECe'=df$ece_yield_pred_from_rel,
                    'Yr_lin_ECa'=df$eca_relative_yield,
                    'Yr_lin_ECe'=df$ece_relative_yield,
                    'Y_Maas'=df$maas_yield,
                    'Yr_Maas'=df$maas_relative_yield
                    )
  #dfx2
  getwd()
  setwd() # set custom save directory here
  write.csv(dfx2, file = name)
}
prediction_csv(corn_yield, name='38_ECa_ECe_results.csv')
prediction_csv(corn_yield_cal, name='20_ECa_ECe_results.csv')
prediction_csv(corn_yield_cal, name='12_ECa_ECe_results.csv')
prediction_csv(corn_yield_cal, name='6_ECa_ECe_results.csv')
prediction_csv(corn_yield, name='2018_ECa_ECe_results.csv')

#Uses existing model to predict yield on master ANOCOVA dataset
master_df = read.csv(file.choose()) #use "Master ANOCOVA ECe Model Predictions_ver2020" or "extrap_template.csv"
prediction_csv2 = function(df, name, max) {
  dfx2 = data.frame('id'=df$Site_ID,
                    'X'=df$X,
                    'Y'=df$Y,
                    ##################################
                    #used only for extrap_template.csv
                    #'Y_obs'=corn_yield$kg.ha,
                    #'Yr_obs'=corn_yield$relative_yield,
                    ##################################
                    'ECa_obs'=df$CV.1.0m,
                    'ECe_pred'=df$ECe,
                    'Field_ID'=df$Field,
                    ##'Y_ANOCOVA_ECa'=df$anocova_eca_yield_pred_from_rel,
                    #'Yr_ANOCOVA_ECa'=predict(anocova_eca_model, newdata=df),
                    ##"Y_sig_ECa"=df$eca_sig_Y,
                    #"Yr_sig_ECa"=predict(eca_sig_model, newdata=df),
                    ##"Y_sig_ECe"=df$ece_sig_Y,
                    #"Yr_sig_ECe"=predict(ece_sig_model, newdata=df),
                    ##"Y_Gen_ECe"=df$ece_gen_Y,
                    ##"Y_Gen_ECa"=df$eca_gen_Y,
                    "Yr_Gen_ECa"=predict(eca_genuchten_model, newdata=df),
                    "Yr_Gen_ECe"=predict(ece_genuchten_model, newdata=df),
                    ##'Y_lin_ECa'=df$eca_yield_pred_from_rel,
                    ##'Y_lin_ECe'=df$ece_yield_pred_from_rel,
                    'Yr_lin_ECa'=predict(eca_model, newdata=df),
                    'Yr_lin_ECe'=predict(ece_model, newdata=df),
                    ##'Y_Maas_ECe'=df$maas_yield,
                    'Yr_Maas_ECe'= ifelse(df$ECe >= 1.7, 100 - 12*(df$ECe - 1.7), 100)
  )
  #dfx2
  dfx2$Y_Gen_ECa <- c(0)
  dfx2$Y_Gen_ECe <- c(0)
  dfx2$Y_lin_ECa <- c(0)
  dfx2$Y_lin_ECe <- c(0)
  for(i in 1:length(dfx2$Field)) {
    x2 = dfx2$Field
    y2 = dfx2$Yr_Gen_ECa
    y3 = dfx2$Yr_Gen_ECe
    y4 = dfx2$Yr_lin_ECa
    y5 = dfx2$Yr_lin_ECe
    dfx2$Y_Gen_ECa[i] <- ((y2[i]/100)*max$max_yield[max$Field == x2[i]])
    dfx2$Y_Gen_ECe[i] <- ((y3[i]/100)*max$max_yield[max$Field == x2[i]])
    dfx2$Y_lin_ECa[i] <- ((y4[i]/100)*max$max_yield[max$Field == x2[i]])
    dfx2$Y_lin_ECe[i] <- ((y5[i]/100)*max$max_yield[max$Field == x2[i]])
  }
  getwd()
  setwd() # set custom save directory here
  write.csv(dfx2, file = name)
}
prediction_csv2(df=master_df, name='38_Master_ECa_ECe_results.csv')
prediction_csv2(df=master_df, name='20_Master_ECa_ECe_results.csv')
prediction_csv2(df=master_df, name='12_Master_ECa_ECe_results.csv')
prediction_csv2(df=master_df, name='6_Master_ECa_ECe_results.csv')

prediction_csv2(df=master_df, name='6_38extrap_ECa_ECe_results.csv', max = max_yields_cal)
prediction_csv2(df=master_df, name='12_38extrap_ECa_ECe_results.csv', max = max_yields_cal)
prediction_csv2(df=master_df, name='20_38extrap_ECa_ECe_results.csv', max = max_yields_cal)
########################################################################Graphing
  #Plot of ln(ECa) v. ln(kg/ha of dry biomass)
  plot(corn_yield$kg.ha~ corn_yield$CV.1.0m)
  plot(corn_yield$lnkg~ corn_yield$ln1.0m)
  g2 = subset(corn_yield1, Field=='G2' )
  mdl = lm(kg.ha~ECe_Deep, data = g2)
  summary(mdl)
  mdl2 = lm(kg.ha~CV.1.0m, data = g2)
  summary(mdl2)

  #Plot of all fields with linear regression lines in same plot
  trend_plot = function(){
    a = ggplot(corn_yield, aes(CV.1.0m, kg.ha, color=factor(Field))) +
        geom_point() +
        stat_smooth(method = "lm", se = F) +
        xlab('1.0m ECa, dS/m') +
        ylab('Dry Biomass, kg/ha') +
        ggtitle("ECa 1.0m v. Corn Yield for 2019") +
        theme(plot.title = element_text(hjust = 0.5))
    b = ggplot(corn_yield, aes(CV.0.5m, kg.ha, color=factor(Field))) +
        geom_point() +
        stat_smooth(method = "lm", se = F) +
        xlab('0.5m ECa, dS/m') +
        ylab('Dry Biomass, kg/ha') +
        ggtitle("ECa 0.5m v. Corn Yield for 2019") +
        theme(plot.title = element_text(hjust = 0.5))
    c = ggplot(corn_yield, aes(ln1.0m, lnkg, color=factor(Field))) +
        geom_point() +
        stat_smooth(method = "lm", se = F) +
        xlab('ln(1.0m ECa, dS/m)') +
        ylab('ln(Dry Biomass, kg/ha)') +
        ggtitle("ln(ECa 1.0m) v. ln(Corn Yield) for 2018") +
        theme(plot.title = element_text(hjust = 0.5))
    d = ggplot(corn_yield, aes(ln0.5m, lnkg, color=factor(Field))) +
        geom_point() +
        stat_smooth(method = "lm", se = F) +
        xlab('ln(0.5m ECa, dS/m)') +
        ylab('ln(Dry Biomass, kg/ha)') +
        ggtitle("ln(ECa 0.5m) v. ln(Corn Yield) for 2018") +
        theme(plot.title = element_text(hjust = 0.5))
   grid.arrange(a, b, nrow=2)
  }
  trend_plot()

  #plots of each field separately
  facets = function(df){
    a = qplot(CV.1.0m, kg.ha, data=df, geom='point',xlab="ECa (1.5m)",
             ylab="kg/ha",facets=Field~.) + stat_smooth(method = "lm", se = FALSE)
    a2 = qplot(ECe_Deep, kg.ha, data=df, geom='point',xlab="ECe (1.2m)",
              ylab="kg/ha",facets=Field~.) + stat_smooth(method = "lm", se = FALSE)
    b = qplot(CV.0.5m, kg.ha, data=df, geom='point',xlab="ECa (0.75m)",
              ylab="kg/ha",facets=Field~.) + stat_smooth(method = "lm", se = FALSE)
    c = qplot(ln1.0m, lnkg, data=df, geom='point',xlab="ln(ECa; 1.5m)",
        ylab="ln(kg/ha)",facets=Field~.) + stat_smooth(method = "lm", se = FALSE)

    d = qplot(ln0.5m, lnkg, data=df, geom='point',xlab="ln(ECa; 0.75m)",
        ylab="ln(kg/ha)",facets=Field~.) + stat_smooth(method = "lm", se = FALSE)
    grid.arrange(a, a2, nrow=1, ncol=2)
    }
  facets(df=corn_yield1)

  #ggplot of ECe and ECa v. Relative Yield % and Yield
  linear_plot = function() {
    a = ggplot(corn_yield, aes(ECe, relative_yield)) +
      geom_point() +
      stat_smooth(method = "lm", se = F) +
      xlab('ECe, dS/m') +
      ylab('Relative Yield %') +
      ggtitle("ECe v. Relative Corn Yield % for 2018") +
      theme(plot.title = element_text(hjust = 0.5))
    b = ggplot(corn_yield, aes(ln1.0m, ln_relative_yield)) +
      geom_point() +
      stat_smooth(method = "lm", se = F) +
      xlab('ln(ECa, mS/m)') +
      ylab('ln(Relative Yield %') +
      ggtitle("ln(ECa) v. ln(Relative Corn Yield %) for 2018") +
      theme(plot.title = element_text(hjust = 0.5))
    c = ggplot(corn_yield, aes(ECe, kg.ha, color=factor(Field))) +
      geom_point() +
      stat_smooth(method = "lm", se = F) +
      xlab('ECe, dS/m') +
      ylab('Yield, kg/ha') +
      ggtitle("ECe v. Corn Yield for 2018") +
      theme(plot.title = element_text(hjust = 0.5))
    d = ggplot(corn_yield, aes(ln1.0m, lnkg, color=factor(Field))) +
      geom_point() +
      stat_smooth(method = "lm", se = F) +
      xlab('ln(ECa, mS/m)') +
      ylab('ln(Yield, kg/ha)') +
      ggtitle("ln(ECa) v. ln(Corn Yield) for 2018") +
      theme(plot.title = element_text(hjust = 0.5))
    grid.arrange(a,b,c,d, nrow = 2, ncol = 2)
  }
  linear_plot()

  #1:1 Plots of Yield
  one_to_one = function(){
    #Linear ECa Model
    a = qplot(kg.ha, eca_yield_pred_from_rel, data = corn_yield, geom = 'point', color = Field) +
          geom_abline(slope = 1) +
          coord_fixed(xlim = c(4000, 21000), ylim = c(4000, 21000)) +
          ggtitle("1:1 Plot of Observed and \n ECa Linear Model Predicted Yield") +
          xlab("Observed Corn Yield, kg/ha") +
          ylab("Predicted Corn Yield, kg/ha") +
          theme(plot.title = element_text(hjust = 0.5))

    #Linear ECe model
    b = qplot(kg.ha, ece_yield_pred_from_rel, data = corn_yield, geom = 'point', color = Field) +
          geom_abline(slope = 1) +
          coord_fixed(xlim = c(4000, 21000), ylim = c(4000, 21000)) +
          ggtitle("1:1 Plot of Observed and \n ECe Linear Model Predicted Yield") +
          xlab("Observed Corn Yield, kg/ha") +
          ylab("Predicted Corn Yield, kg/ha") +
          theme(plot.title = element_text(hjust = 0.5))

    #ANOCOVA ECa Model
    c = qplot(kg.ha, anocova_eca_yield_pred_from_rel, data = corn_yield, geom = 'point', color = Field) +
      geom_abline(slope = 1) +
      coord_fixed(xlim = c(4000, 21000), ylim = c(4000, 21000)) +
      ggtitle("1:1 Plot of Observed and \n ECa ANOCOVA Model Predicted Yield") +
      xlab("Observed Corn Yield, kg/ha") +
      ylab("Predicted Corn Yield, kg/ha") +
      theme(plot.title = element_text(hjust = 0.5))

    #Linear Maas Model
    d = qplot(kg.ha, maas_yield, data = corn_yield, geom = 'point', color = Field) +
          geom_abline(slope = 1) +
          coord_fixed(xlim = c(4000, 21000), ylim = c(4000, 21000)) +
          ggtitle("1:1 Plot of Observed and (Maas and Hoffman, 1977) \n Model Predicted Yield") +
          xlab("Observed Corn Yield, kg/ha") +
          ylab("Predicted Corn Yield, kg/ha") +
          theme(plot.title = element_text(hjust = 0.5))

    grid.arrange(a, b, c, d, nrow = 2, ncol = 2)
  }
  one_to_one()
