## Model training

library(e1071)
library(caret)
library(pROC)
library(randomForest)
library(ggplot2)
library(PRROC)
library(dplyr)
library(parallel)

# load data
df_unlabeled <- read.table("Unlabeled_dataset.txt", header = T)
df_positive <- read.table("Positive_dataset.txt", header = T)


set.seed(202312)
p_fold <- createMultiFolds(y=df_positive$Label, k = 10, times = 1)


rf_roc <- list()
auc_value <- list()
pr <- list()
auPR <- list()
F1score <- list()
FNR <- list()
rankInx <- list()
finalPU <- list()
res_val <- vector("list", 10)
res_model <- list()

for(i in 1:10){
  
  train_u <- df_unlabeled
  
  train_p <- df_positive[ p_fold[[i]],] 
  valid_p <- df_positive[-p_fold[[i]],] 
  
  set.seed(i)
  index <- replicate(100, sample(seq_len(nrow(train_u)), replace = F, size = nrow(train_p))) 
  
  mclapply(1:100, function(x){
    print(paste0("x = ", x))
    
    train_n <- train_u[index[,x],]
    train_n$Label <- 0
    train_p$Label <- 1
    
    ### Initial training
    train_np <- rbind(train_p, train_n )
    train_np <- train_np[,-1]
    train_np$Label <- as.factor( train_np$Label)
    
    set.seed(202313)
    Classifier <- svm(Label~.,data = train_np, kernel="radial", probability = T)
    
    
    ### Initial predicting
    pred_svm <- predict(Classifier, newdata = train_u[,-c(1,21)], probability = T)
    pred_svm <- attr(pred_svm, "probabilities") %>% as.data.frame()
    
    tr_1 = 0.8
    train_p <- rbind(train_p, train_u[ pred_svm$`1` >= tr_1, ])
    train_p$Label <- 1
    
    tr_2 <- 0.2
    
    train_n <- train_u[pred_svm$`1` <= tr_2, ]
    if (nrow(train_n) > 0) {
      train_n$Label <- 0
    } else {
      tr_2 <- 0.3
      train_n <- train_u[pred_svm$`1` <= tr_2, ]
      train_n$Label <- 0
    }
    
    idx = c( which(pred_svm$`1` >= tr_1), which(pred_svm$`1` <= tr_2))
    
    train_np <- rbind(train_p, train_n )
    train_np <- train_np[,-1]
    train_np$Label <- as.factor( train_np$Label)
    
    set.seed(202313)
    Classifier <- svm(Label~.,data = train_np, kernel="radial", probability = T)
    
    
    ### again training
    train_q <- train_u[-idx, ]
    n = nrow(train_q)
    while ( length(idx) > 0 & n > 0 ){
      if (nrow(train_q) == 0) break
      print(paste0("n = ", n))
      
      pred_svm <- predict(Classifier, newdata = train_q[,-c(1,21)], probability = T)
      pred_svm <- attr(pred_svm, "probabilities") %>% as.data.frame()
      
      tr_1 <- 0.9
      train_p <- rbind(train_p, train_q[ pred_svm$`1` >= tr_1, ])
      train_p$Label <- 1
      
      tr_2 <- 0.1
      train_n <- rbind(train_n, train_q[ pred_svm$`1` <= tr_2, ])
      train_n$Label <- 0
      
      idx = c( which(pred_svm$`1` >= tr_1), which(pred_svm$`1` <= tr_2))
      
      train_np <- rbind(train_p, train_n )
      train_np$Label <- as.factor( train_np$Label)
      
      set.seed(202313)
      Classifier <- svm(Label~.,data = train_np[,-1], kernel="radial", probability = T)
      
      train_q <- train_q[-idx, ]
      n = nrow(train_q)
      
    }
    
    temp <- predict(Classifier, newdata = valid_p[, -c(1, 21)])
    TPR <- sum(temp == 1)/nrow(valid_p)
    resval <- (TPR >= 5/6) 
    
    res <- train_np[,c(1,21)]
    colnames(res)[2] <- x
    
    return( list(val = resval, PU = res) )
  }, mc.cores = 6) -> res_tmp
  
  res_PU <- list()
  for(x in 1:100) {
    res_val[[i]][[x]] <- res_tmp[[x]]$val
    res_PU[[x]] <- res_tmp[[x]]$PU
  }
  
  ## filter replicates by TPR and calculate final PN datasets
  finalPU[[i]] <- Reduce(function(x,y){merge(x,y, all = T)}, res_PU[unlist(res_val[[i]])])
  finalPU[[i]] <- finalPU[[i]][ rowSums(is.na(finalPU[[i]][,-1]))/ncol(finalPU[[i]][,-1]) <= 0.2, ]
  finalPU[[i]] <- finalPU[[i]] %>% mutate(across(-1, as.character)) %>% mutate(across(-1, as.numeric))
  finalPU[[i]]$Mean <- rowMeans(finalPU[[i]][,-1], na.rm = T)
  finalPU[[i]]$Label <- ifelse(finalPU[[i]]$Mean >= 0.5, 1, 0)
  
  
  ## training RF
  finalPU[[i]] <- rbind(train_p,train_u)[,-21] %>% merge(finalPU[[i]][,c("ASID", "Label")], ., by = "ASID", all.x = T)
  finalPU[[i]]$Label <-  as.factor( finalPU[[i]]$Label)
  finalPU[[i]] <- finalPU[[i]][,c(1, 3:21, 2)]
  
  set.seed(123)
  train_index <- createDataPartition(
    y = finalPU[[i]]$Label,     
    p = 0.8,        
    list = FALSE,      
    times = 1          
  )
  
  train_data <- finalPU[[i]][train_index, ]  
  test_data  <- finalPU[[i]][-train_index, ] 
  
  res_model[[i]] <- randomForest(Label ~ ., data = train_data[,-1], mtry =7, ntree = 1000, importance =T)
  resPred <- predict(res_model[[i]], newdata = test_data[,-c(1,21)], type = "prob") %>% as.data.frame()
  
  print("performance evaluation")
  
  ## performance evaluation
  ### roc and auc
  obs_p_rf <- data.frame(pred = resPred[,2], obs = test_data$Label)
  rf_roc[[i]] <- roc(obs ~ pred, obs_p_rf, levels = c("0", "1" ))
  auc_value[[i]] <- as.numeric(auc(rf_roc[[i]]))
  
  ### auPR
  fg <- obs_p_rf[obs_p_rf$obs == 1, "pred"]
  bg <- obs_p_rf[obs_p_rf$obs == 0, "pred"]
  pr[[i]] <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
  auPR[[i]] <- pr[[i]]$auc.integral
  
  ### rank index
  obs_p_rf <- obs_p_rf[order(obs_p_rf$pred, decreasing = T), ]
  obs_p_rf$Rank <- 1:nrow(obs_p_rf)
  SP <- subset(obs_p_rf, obs == 1)
  rankInx[[i]] <- (1/nrow(SP)) * (sum(SP$Rank)/nrow(obs_p_rf))
  
  ### F1 score and false negative ratio
  tmp <- confusionMatrix((ifelse(obs_p_rf$pred >= 0.5, "1", "0") %>% as.factor()),
                         obs_p_rf$obs,
                         mode = "everything",
                         positive = "1")
  F1score[[i]] <- tmp$byClass["F1"]
  recall <- tmp$byClass["Recall"]
  FNR[[i]] <- 1 - as.numeric(recall)
  
}

save(res_val, finalPU, rf_roc, auc_value, pr, auPR, F1score, FNR, rankInx, file = "Model_evaluation.RData")
saveRDS(res_model, file = "finalModel.Rds")


