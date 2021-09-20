# cmccmpredfunc prediction function
#
#
library(dplyr)
require(caret)
library(MASS)

#gender 性别.女0.男1..
#ecog ECOG评分.0.小于等于2分.1.大于等于3分.
#pathological_types 病理分型2.0.非腺癌.1.腺癌.
#ganyu 本虚3.肝郁脾虚.0否.1是.
#pishen 本虚1.脾肾亏虚.0否.1是.
#line 门诊初诊时西医治疗线数.0.小于等于一线.1.二线及以上.
#left 左右半分类.左0.右1.
#transfer 转移部位1.0.肝.腹腔.卵巢.脑.腹膜.1.以上2.3合并.

model_file_path = '/home/opencpu/test2/model/cmccm_2group.rda'
csv_data_path = 'cmccm_2group.csv'
buildmodel = function(csv_data_path){

  data_model = read.csv(csv_data_path, fileEncoding = 'GBK',stringsAsFactors = F)
  data_model <- data_model[,c(2,3,4,5,6,7,8,9,10)] # 选出有重要影响的协变量(按需调整)
  for (i in c(1,2,3,4,5,6,7,8,9)) {
    data_model[,i] <- as.factor(data_model[,i])} # 将分类变量转换成因子型
  # 修改列名
  names(data_model)<-c("group","gender","ecog","pathological_types","ganyu","pishen","line","left","transfer")
  # 人群分类
  data_model[,"group"] = ifelse(data_model[,"group"]==1,0,1)
  data_model[,"gender"] = ifelse(data_model[,"gender"]==1,0,1)
  data_model[,"ecog"] = ifelse(data_model[,"ecog"]==1,0,1)
  data_model[,"pathological_types"] = ifelse(data_model[,"pathological_types"]==1,0,1)
  data_model[,"ganyu"] = ifelse(data_model[,"ganyu"]==1,0,1)
  data_model[,"pishen"] = ifelse(data_model[,"pishen"]==1,0,1)
  data_model[,"line"] = ifelse(data_model[,"line"]==1,0,1)
  data_model[,"left"] = ifelse(data_model[,"left"]==1,0,1)
  data_model[,"transfer"] = ifelse(data_model[,"transfer"]==1,0,1)

  #建立模型
  set.seed(7)
  require(caret)
  folds <- createFolds(y=data_model[,1],k=5) # 五折交叉验证

  new_data1 <- as.data.frame(lapply(data_model, as.numeric))

  #构建for循环，得5次交叉验证的测试集精确度、训练集精确度
  max=0
  num=0
  for(i in 1:5){

    fold_test <- new_data1[folds[[i]],]   # 取folds[[i]]作为测试集
    fold_train <- new_data1[-folds[[i]],]   # 剩下的数据作为训练集

    print("***组号***")

    fold_preq <- qda(group~.,data=fold_train) # 训练集训练模型
    fold_predictq <- predict(fold_preq,newdata=fold_test[-1]) # 在测试集预测
    ## 注意这里fold_predictq结果中1相当于group中的0；fold_predictq结果中2相当于group中的1
    newGroupq1 <- fold_predictq$class # 测试集预测结果
    tabq1 <- table(fold_test[,1],newGroupq1) # 测试集真实值和预测值的混淆矩阵
    tabq1
    fold_accuracyq = sum(diag(prop.table(tabq1))) # 测试集精确度
    print(i)
    print("***测试集精确度***")
    print(fold_accuracyq)
    print(tabq1)
    print("***训练集精确度***")
    fold_predictq2 <- predict(fold_preq,newdata=fold_train[-1]) # 在训练集预测
    newGroupq2 <- fold_predictq2$class # 训练集预测结果
    tabq2 <- table(fold_train[,1],newGroupq2) # 训练集真实值和预测值的混淆矩阵
    tabq2
    fold_accuracyq2 = sum(diag(prop.table(tabq2))) # 训练集精确度
    print(fold_accuracyq2)
    print(tabq2)


    if(fold_accuracyq>max)
    {
      max=fold_accuracyq
      num=i
      best_model = fold_preq
    }

  }

  save(best_model, file=model_file_path)


  ### 最终结果汇总
  fold_test <- data_model[folds[[num]],]   # 最终结果的测试集
  fold_train <- data_model[-folds[[num]],]  # 最终结果训练集
  fold_predictq <- predict(best_model,newdata=fold_test[-1])
  fold_predictq2 <- predict(best_model,newdata=fold_train[-1])
  fold_test['pre_group'] = fold_predictq
  test_posterior = fold_predictq$posterior # 测试集的后验概率
  fold_train['pre_group'] = fold_predictq2
  train_posterior = fold_predictq2$posterior # 训练集的后验概率
}

predictfunc = function(input){
  #print('cmccmpredfunc')
  library(MASS)
  #"gender","ecog","pathological_types","ganyu","pishen","line","left","transfer"
  New_patient = as.data.frame(lapply(input,as.numeric))#as.data.frame(input)
  New_patient[,"gender"] = New_patient[,"f1"]
  New_patient[,"ecog"] = New_patient[,"f2"]
  New_patient[,"pathological_types"] = New_patient[,"f3"]
  New_patient[,"ganyu"] = New_patient[,"f4"]
  New_patient[,"pishen"] = New_patient[,"f5"]
  New_patient[,"line"] = sqrt(New_patient[,"f6"])
  New_patient[,"left"] = New_patient[,"f7"]
  New_patient[,"transfer"] = New_patient[,"f8"]

  load(model_file_path)
  predmodel.train.new_patient = predict(best_model,New_patient)
  #predmodel.train.new_patient
  #print(predmodel.train.new_patient)
  #decide = ifelse(predmodel.train.new_patient$class=='优势人群','属于','不属于')
  #result = paste0('此患者属于优势人群概率为：',round(predmodel.train.new_patient$posterior[,1],2),'; ',decide,'优势人群。')
  #print(result)
  #as.character(result)
  predmodel.train.new_patient
}

testpredict = function(){
  #1    0                  1     0      1    0    1        0
  New_patient = data.frame('f1'=1,'f2'=0,'f3'=1,'f4'=0,"f5"=1,"f6"=0,"f7"=1,"f8"=0)
  print(predictfunc(New_patient))
  #0    1                  1     1      0    0    1        0
  New_patient = data.frame('f1'=0,'f2'=1,'f3'=1,'f4'=1,"f5"=0,"f6"=0,"f7"=1,"f8"=0)
  print(predictfunc(New_patient))
}
