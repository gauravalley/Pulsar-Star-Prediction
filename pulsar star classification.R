#CLEARING THE ENVIRONMENT
rm(list = ls())

#IMPORTING THE DATA
setwd('C:\\Users\\Shekhar Lamba\\Documents\\Datasets')
star <- read.csv('pulsar_stars.csv')
head(star)
str(star)

#DATA PREPARATION
table(star$target_class)
star$target_class <- ifelse(star$target_class == 1, 'P', 'NP')
star$target_class <- factor(star$target_class)
summary(star)

#EDA
my_stats <- function(x) {
  n <- length(x)
  m <- mean(x)
  s <- sd(x)
  mini <- min(x)
  p1 <- quantile(x, .01)
  p3 <- quantile(x, .03)
  p5 <- quantile(x, .05)
  p10 <- quantile(x, .1)
  q1 <- quantile(x, .25)
  q2 <- quantile(x, .5)
  q3 <- quantile(x, .75)
  p90 <- quantile(x, .90)
  p95 <- quantile(x, .95)
  p97 <- quantile(x, .97)
  p99 <- quantile(x, .99)
  maxi <- max(x)
  lc <- m - (3 * s)
  uc <- m + (3 * s)
  return(c(No_of_values = n, Mean = m, Std_dev = s, LC = lc, Min = mini, P1 = p1, P3 = p3, P5 = p5, P10 = p10, 
           Q1 = q1, Q2 = q2, Q3 = q3, P90 = p90, P95 = p95, P97 = p97, P99 = p99, Max = maxi, UC = uc))
}
num_var <- sapply(star, is.numeric)
diag_stats <- t(data.frame(apply(star[num_var], 2, my_stats)))
View(diag_stats)

#GRAPHICAL ANALYSIS
library(ggplot2)
library(tidyverse)
star %>%
  gather(key = 'variable', value = 'value', c(1:4)) %>%
  ggplot(aes(x = reorder(variable, value, FUN = median), y = value, fill = variable)) +
  geom_boxplot(show.legend = F) +
  labs(x = element_blank(), y = element_blank(), title = 'Boxplot for Variables') +
  theme_bw() + coord_flip()
star %>%
  gather(key = 'variable', value = 'value', c(5:8)) %>%
  ggplot(aes(x = reorder(variable, value, FUN = median), y = value, fill = variable)) +
  geom_boxplot(show.legend = F) +
  labs(x = element_blank(), y = element_blank(), title = 'Boxplot for Variables') +
  theme_bw() + coord_flip()

#OUTLIER TREATMENT
star$Mean.of.the.integrated.profile[star$Mean.of.the.integrated.profile < 34.121162] <- 34.121162
star$Standard.deviation.of.the.integrated.profile[star$Standard.deviation.of.the.integrated.profile < 26.019963] <- 26.019963
star$Mean.of.the.integrated.profile[star$Mean.of.the.integrated.profile > 188.038774] <- 188.038774
star$Standard.deviation.of.the.integrated.profile[star$Standard.deviation.of.the.integrated.profile > 67.079100] <- 67.079100
star$Excess.kurtosis.of.the.integrated.profile[star$Excess.kurtosis.of.the.integrated.profile > 3.669976] <- 3.669976
star$Skewness.of.the.integrated.profile[star$Skewness.of.the.integrated.profile > 20.274019] <- 20.274019
star$Mean.of.the.DM.SNR.curve[star$Mean.of.the.DM.SNR.curve > 101.033091] <- 101.033091
star$Standard.deviation.of.the.DM.SNR.curve[star$Standard.deviation.of.the.DM.SNR.curve > 84.738232] <- 84.738232
star$Excess.kurtosis.of.the.DM.SNR.curve[star$Excess.kurtosis.of.the.DM.SNR.curve > 21.821832] <- 21.821832
star$Skewness.of.the.DM.SNR.curve[star$Skewness.of.the.DM.SNR.curve > 424.401327] <- 424.401327

#MORE GRAPHICAL INSIGHTS
library(GGally)
library(corrplot)
ggplot(star, aes(x = target_class, fill = target_class)) + geom_bar() + labs(title = 'Density Comparison of Non-Pulsar Star and Pulsar star')
corrplot(corr = cor(star[num_var], use = 'pairwise.complete.obs'), method = 'number',number.cex = 0.7, tl.cex = 0.55)
ggpairs(star, legend = 1, columns = 1:8, mapping = aes(color = target_class),
        lower = list(continuous = wrap('smooth', alpha = .5, size = .5))) + theme(legend.position = 'bottom')

#DIVIDING THE DATA SET INTO TRAIN AND TEST
set.seed(1234)
train_ind <- sample(1:nrow(star), size = floor(.8 * nrow(star)))
training <- star[train_ind, ]
testing <- star[-train_ind, ]
prop.table(table(training$target_class)) 
prop.table(table(testing$target_class))

#PREPARING DATA FOR H2O
library(h2o)
h2o.init(nthreads = -1, max_mem_size = '2G')
h2o.removeAll()
train <- as.h2o(training)
test <- as.h2o(testing)
training_data <- h2o.assign(train, 'train.hex')
testing_data <- h2o.assign(test, 'test.hex')

#BAGGING___________________________________________________________________________________________________

#INTIAL RANDOM FOREST MODEL
rf_model_init <- h2o.randomForest(training_frame = training_data, validation_frame = testing_data,
                                  x = 1:8, y = 9, model_id = 'rf_model_init', seed = 1234)
rf_model_init@model$validation_metrics
h2o.varimp_plot(rf_model_init)

#INITIAL RF PREDICTION AND ACCURACY
library(caret)
library(ROSE)
rf_predict_init <- h2o.predict(object = rf_model_init, newdata = testing_data) %>% as.data.frame() %>% pull(predict)
confusionMatrix(rf_predict_init, testing$target_class, positive = 'P') #Acc = 98.1%, Rec = 88.18%, Pre = 90.94%
roc.curve(testing$target_class, rf_predict_init, plotit = T) #AUC = 0.936

#MODIFIED RANDOM FOREST MODEL
rf_model_mod <- h2o.randomForest(training_frame = training_data, validation_frame = testing_data,
                                 x = 1:8, y = 9, model_id = 'rf_model_mod', seed = 1234,
                                 ntrees = 600, mtries = 4, stopping_rounds = 30, stopping_metric = 'AUC',
                                 stopping_tolerance = 0.000001, balance_classes = T, nfolds = 10, 
                                 score_each_iteration = T)
rf_model_mod@model$cross_validation_metrics
h2o.varimp_plot(rf_model_mod)

#FINAL RF PREDICTION AND ACCURACY
rf_predict_final <- h2o.predict(object = rf_model_mod, newdata = testing_data) %>% as.data.frame() %>% pull(predict)
confusionMatrix(rf_predict_final, testing$target_class, positive = 'P') #Acc = 97.91%, Rec = 89.39%, Pre = 88.06%
roc.curve(testing$target_class, rf_predict_final, plotit = T) #AUC = 0.941

#BOOSTING_________________________________________________________________________________________________

#INITIAL GRADIENT BOOST MODEL
gbm_model_init <- h2o.gbm(training_frame = training_data, validation_frame = testing_data,
                          x = 1:8, y = 9, model_id = 'gbm_model_init', seed = 1234)
gbm_model_init@model$validation_metrics
h2o.varimp_plot(gbm_model_init)

#INITIAL GBM PRECICTION AND ACCURACY
gbm_predict_init <- h2o.predict(object = gbm_model_init, newdata = testing_data) %>% as.data.frame() %>% pull(predict)
confusionMatrix(gbm_predict_init, testing$target_class, positive = 'P') #Acc = 98.07%, Rec = 88.18%, Pre = 90.65%
roc.curve(testing$target_class, gbm_predict_init, plotit = T) #AUC = 0.936

#MODIFIED GRADIENT BOOST MODEL
gbm_model_mod <- h2o.gbm(training_frame = training_data, validation_frame = testing_data,
                         x = 1:8, y = 9, model_id = 'gbm_model_mod', seed = 1234, ntrees = 1000,
                         max_depth = 2, learn_rate = 0.1, sample_rate = 0.8, col_sample_rate = 0.8,
                         stopping_metric = 'AUC', stopping_rounds = 6, stopping_tolerance = 1e-5,
                         balance_classes = T, score_tree_interval = 10, nfolds = 10)
gbm_model_mod@model$cross_validation_metrics
h2o.varimp_plot(gbm_model_mod)

#FINAL RF PREDICTION AND ACCURACY
gbm_predict_final <- h2o.predict(object = gbm_model_mod, newdata = testing_data) %>% as.data.frame() %>% pull(predict)
confusionMatrix(gbm_predict_final, testing$target_class, positive = 'P') #Acc = 98.16%, Rec = 89.70%, Pre = 90.24%
roc.curve(testing$target_class, gbm_predict_final, plotit = T) #AUC = 0.944
