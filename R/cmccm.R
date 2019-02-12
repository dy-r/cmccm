# cmccmpredfunc prediction function
#
#
library(MASS)

#
model_file_path = '/home/opencpu/test2/model/cmccm.rda'
csv_data_path = '西苑数据-判别分析原始变量.csv'
cmccmbuildmodel = function(csv_data_path){

  data_model = read.csv(csv_data_path, fileEncoding = 'GBK',stringsAsFactors = F)

  # 人群分类
  data_model[,"event"] = ifelse(data_model[,"死亡日期"]=='',0,1)
  data_model[,"group"] <- ifelse(data_model[,"event"]==1&data_model[,"IV期生存时间"] <= 365,'劣势人群',ifelse((data_model[,"基因汇总"]=='野生型'&data_model[,"IV期生存时间"] > 30*30)|
                                                                                           (data_model[,"基因汇总"]=='BRAF突变型'&data_model[,"IV期生存时间"] > 18*30)|
                                                                                           (data_model[,"基因汇总"]=='RAS突变型'&data_model[,"IV期生存时间"] > 24*30)|
                                                                                           (data_model[,"基因汇总"]=='未查'&data_model[,"IV期生存时间"] > 18*30),'优势人群','中间人群'))

  data_model[,"肿瘤原发部位_左半肠"] = ifelse(data_model[,"肿瘤原发部位汇总"]=="左半肠",1,0)
  data_model[,"BRAF突变型"] = ifelse(data_model[,"基因汇总"]=='BRAF突变型',1,0)
  data_model[,"RAS突变型"] = ifelse(data_model[,"基因汇总"]=='RAS突变型',1,0)


  #变量处理1-连续变量开根号处理-KPS评分、中医临床症状评分总分值、近3个月主症评分总分、症状妨碍生活评分总分
  data_model[,"KPS评分"]=sqrt(data_model[,"KPS评分"])
  data_model[,"中医临床症状评分总分值"]=sqrt(data_model[,"中医临床症状评分总分值"])
  data_model[,"近3个月主症评分总分"]=sqrt(data_model[,"近3个月主症评分总分"])
  data_model[,"症状妨碍生活评分总分"]=sqrt(data_model[,"症状妨碍生活评分总分"])

  #变量处理2-分类变量重编码为0-1变量
  data_model[,"肝转移"]=ifelse(data_model[,"肝转移"]=='是',1,0)
  data_model[,"肺转移"]=ifelse(data_model[,"肺转移"]=='是',1,0)

  data_model[,"肺"]=ifelse(data_model[,"肺"]=='有',1,0)
  data_model[,"肾"]=ifelse(data_model[,"肾"]=='有',1,0)
  data_model[,"肝"]=ifelse(data_model[,"肝"]=='有',1,0)
  data_model[,"脾"]=ifelse(data_model[,"脾"]=='有',1,0)

  #建立模型
  data_lda = data_model[,c('group','KPS评分','肝转移','肺转移','肺',"肾","肝","脾",
                           "中医临床症状评分总分值","近3个月主症评分总分","症状妨碍生活评分总分","肿瘤原发部位_左半肠",
                           "BRAF突变型","RAS突变型")]

  qda.model <- qda( group~ .,data=na.omit(data_lda))

  save(qda.model, file=model_file_path)
}

cmccmpredfunc = function(input){
  #print('cmccmpredfunc')
  library(MASS)
  New_patient = as.data.frame(lapply(input,as.numeric))#as.data.frame(input)
  New_patient[,"肿瘤原发部位_左半肠"] = New_patient[,"f1"]
  New_patient[,"RAS突变型"] = New_patient[,"f2"]
  New_patient[,"BRAF突变型"] = New_patient[,"f3"]
  New_patient[,"肝转移"] = New_patient[,"f4"]
  New_patient[,"肺转移"] = New_patient[,"f5"]
  New_patient[,"KPS评分"] = sqrt(New_patient[,"f6"])
  New_patient[,"肝"] = New_patient[,"f7"]
  New_patient[,"肺"] = New_patient[,"f8"]
  New_patient[,"肾"] = New_patient[,"f9"]
  New_patient[,"脾"] = New_patient[,"f10"]
  New_patient[,"中医临床症状评分总分值"] = sqrt(New_patient[,"f11"])
  New_patient[,"近3个月主症评分总分"] = sqrt(New_patient[,"f12"])
  New_patient[,"症状妨碍生活评分总分"] = sqrt(New_patient[,"f13"])

  load(model_file_path)
  predmodel.train.new_patient = predict(qda.model,New_patient)
  #predmodel.train.new_patient
  #print(predmodel.train.new_patient)
  decide = ifelse(predmodel.train.new_patient$class=='优势人群','属于','不属于')
  result = paste0('此患者属于优势人群概率为：',round(predmodel.train.new_patient$posterior[,2],2),'; ',decide,'优势人群。')
  #print(result)
  as.character(result)
}

testpredict = function(){
  New_patient = data.frame('f4'=70,'f5'=1,'f6'=0,'f7'=1,"f8"=0,"f9"=1,"f10"=0,
                           "f11"=24,"f12"=7,"f13"=11,
                           "f1"=0,"f2"=0,"f3"=1)
  print(cmccmpredfunc(New_patient))
  New_patient = data.frame('f4'=sqrt(70),'f5'=1,'f6'=0,'f7'=1,"f8"=0,"f9"=1,"f10"=0,
                           "f11"=sqrt(24),"f12"=sqrt(7),"f13"=sqrt(11),
                           "f1"=0,"f2"=0,"f3"=1)
  print(cmccmpredfunc(New_patient))
}
