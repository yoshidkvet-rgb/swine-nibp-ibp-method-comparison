source(file.path("analysis","00_setup.R"))
metrics<-function(truth,pred){TP<-sum(truth&pred);FN<-sum(truth&!pred);FP<-sum(!truth&pred);TN<-sum(!truth&!pred);data.frame(TP=TP,FN=FN,FP=FP,TN=TN,sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),PPV=TP/(TP+FP),NPV=TN/(TN+FN))}
out<-rbind(cbind(target="IBP<70",cutoff="NIBP<70",metrics(primary$ibp_same_map<70,primary$nibp_map<70)),cbind(target="IBP>=120",cutoff="NIBP>=120",metrics(primary$ibp_same_map>=120,primary$nibp_map>=120)),cbind(target="IBP>=120",cutoff="NIBP>=87 exploratory",metrics(primary$ibp_same_map>=120,primary$nibp_map>=87)))
readr::write_csv(out,file.path("output","threshold_2x2.csv"))
readr::write_csv(data.frame(target=c("IBP<70","IBP>=120"),AUC=c(as.numeric(pROC::auc(pROC::roc(primary$ibp_same_map<70,-primary$nibp_map,quiet=TRUE))),as.numeric(pROC::auc(pROC::roc(primary$ibp_same_map>=120,primary$nibp_map,quiet=TRUE))))),file.path("output","threshold_AUC.csv"))
