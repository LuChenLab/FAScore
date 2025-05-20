## Model training

library(caret)
library(pROC)
library(randomForest)
library(ggplot2)
library(caret)

# load data
df_unlabeled <- read.table("Unlabeled_dataset.txt", header = T)
df_positive <- read.table("Positive_dataset.txt", header = T)


# train
set.seed(202312)
u_fold <- createMultiFolds(y=df_unlabeled$Label, k = 10, times = 1)
p_fold <- createMultiFolds(y=df_positive$Label, k = 10, times = 1)


res <- list()
rf_roc <- list()
auc_value <- list()
auPR <- list()
F1score <- list()
rankInx <- list()

for(i in 1:10){
  
  train_u <- df_unlabeled[ u_fold[[i]],] 
  test_u <- df_unlabeled[ -u_fold[[i]],] 
  
  train_p <- df_positive[ p_fold[[i]],] 
  test_p <- df_positive[-p_fold[[i]],] 
  
  df_test <- rbind(test_u, test_p)
  
  set.seed(202313)
  index <- replicate(100, sample(seq_len(nrow(train_u)), replace = F, size = nrow(train_p))) 
  
  lapply(1:100, function(x){
    print(paste0("x = ", x))
    
    train_n <- train_u[index[,x],]
    train_n$Label <- 0
    train_p$Label <- 1
    
    ### Initial training
    train_np <- rbind(train_p, train_n )
    train_np <- train_np[,-1]
    train_np$Label <- as.factor( train_np$Label)
    
    set.seed(202313)
    Classifier <- randomForest(Label ~ ., data = train_np, mtry =7, ntree = 1000, importance =T)
    
    
    ### Initial predicting
    pred_rf <- predict(Classifier, newdata = train_u[,-c(1,21)], type = "prob")
    pred_rf <- as.data.frame(pred_rf)
    
    tr_1 = 0.9
    train_p <- rbind(train_p, train_u[ pred_rf$`1` >= tr_1, ])
    train_p$Label <- 1
    
    tr_2 <- quantile(pred_rf$`1`, probs = seq(0,1,0.1))[3][[1]]
    train_n <- train_u[ pred_rf$`1` <= tr_2, ]
    train_n$Label <- 0
    
    idx = c( which(pred_rf$`1` >= tr_1), which(pred_rf$`1` <= tr_2))
    
    train_np <- rbind(train_p, train_n )
    train_np <- train_np[,-1]
    train_np$Label <- as.factor( train_np$Label)
    
    set.seed(202313)
    Classifier <- randomForest(Label ~ ., data = train_np, mtry =7, ntree = 1000, importance =T)
    
    
    ### again training
    train_q <- train_u[-idx, ]
    n = nrow(train_q)
    while ( n > 0 ){
      print(paste0("n = ", n))
      
      pred_rf <- predict(Classifier, newdata = train_q[,-c(1,21)], type = "prob")
      pred_rf <- as.data.frame(pred_rf)
      
      tr_1 = 0.9
      train_p <- rbind(train_p, train_q[ pred_rf$`1` >= tr_1, ])
      train_p$Label <- 1
      
      tr_2 <- quantile(pred_rf$`1`, probs = seq(0,1,0.1))[3][[1]]
      train_n <- rbind(train_n, train_q[ pred_rf$`1` <= tr_2, ])
      train_n$Label <- 0
      
      idx = c( which(pred_rf$`1` >= tr_1), which(pred_rf$`1` <= tr_2))
      
      train_np <- rbind(train_p, train_n )
      train_np$Label <- as.factor( train_np$Label)
      
      set.seed(202313)
      Classifier <- randomForest(Label ~ ., data = train_np[,-1], mtry =7, ntree = 1000, importance =T)
      
      train_q <- train_q[-idx, ]
      n = nrow(train_q)
      
    }
    return( Classifier )
  }) -> res[[i]]
  
  
  ## performance evaluation
  ### roc and auc
  lapply(1:100, function(y) {
    temp <- predict(res[[i]][[y]], newdata = df_test[, -c(1, 21)], type = "prob") %>% as.data.frame()
    temp[,2]
  }) %>% do.call(cbind, .)  -> tmp0
  
  obs_p_rf <- data.frame(pred = rowMeans(tmp0),
                         obs = df_test$Label)
  
  rf_roc[[i]] <- roc(obs ~ pred, obs_p_rf, levels = c("0", "1" ))
  auc_value[[i]] <- as.numeric(auc(rf_roc[[i]]))
  
  
  ### auPR
  fg <- obs_p_rf[obs_p_rf$obs == 1, "pred"]
  bg <- obs_p_rf[obs_p_rf$obs == 0, "pred"]
  pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  
  auPR[[i]] <- pr$auc.integral
  
  
  ### rank index
  obs_p_rf <- obs_p_rf[order(obs_p_rf$pred, decreasing = T), ]
  obs_p_rf$Rank <- 1:nrow(obs_p_rf)
  SP <- subset(obs_p_rf, obs == 1)
  
  rankInx[[i]] <- (1/nrow(SP)) * (sum(SP$Rank)/nrow(obs_p_rf))
  
  
  ### F1 score
  lapply(1:100, function(y) {
    
    temp <- predict(res[[i]][[y]], newdata = test_p[,-c(1,21)])
    temp <- as.character(temp) %>% as.numeric()
    
  }) %>% do.call(cbind, .) %>% rowMeans(.) -> tmp0
  
  
  tmp <- confusionMatrix((ifelse(tmp0 >= 0.5, "1", "0") %>% as.factor()), 
                         factor(rep("1",length(tmp0)), levels = c(0,1)),
                         mode = "everything",
                         positive = "1")
  
  F1score[[i]] <- tmp$byClass["F1"]
}

which.max(unlist(auc_value))
which.max(unlist(auPR))
which.max(unlist(F1score))
which.min(unlist(rankInx))

finalModel <- res[[which.max(unlist(auc_value))]]

save(auc_value, auPR, F1score, rankInx, file = "Model_evaluation.RData")
saveRDS(finalModel, file = "finalModel.Rds")




