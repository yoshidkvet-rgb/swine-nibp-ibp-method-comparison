source(file.path("analysis","00_setup.R"))
# Deming: preserve first-occurrence pig order so fixed seed matches the verified run.
pig_order<-unique(primary$canonical_pig_id)
pig_stats<-primary|>dplyr::group_by(canonical_pig_id)|>dplyr::summarise(n=dplyr::n(),sx=sum(ibp_same_map),sy=sum(nibp_map),sx2=sum(ibp_same_map^2),sy2=sum(nibp_map^2),sxy=sum(ibp_same_map*nibp_map),.groups="drop")
pig_stats<-pig_stats[match(pig_order,pig_stats$canonical_pig_id),]
deming_w<-function(w,lambda=1){N<-sum(w*pig_stats$n);sx<-sum(w*pig_stats$sx);sy<-sum(w*pig_stats$sy);sx2<-sum(w*pig_stats$sx2);sy2<-sum(w*pig_stats$sy2);sxy<-sum(w*pig_stats$sxy);mx<-sx/N;my<-sy/N;Sxx<-sx2/N-mx^2;Syy<-sy2/N-my^2;Sxy<-sxy/N-mx*my;slope<-(Syy-lambda*Sxx+sqrt((Syy-lambda*Sxx)^2+4*lambda*Sxy^2))/(2*Sxy);c(intercept=my-slope*mx,slope=slope)}
set.seed(20261104);B<-3000;npig<-nrow(pig_stats);dem<-matrix(NA_real_,B,2,dimnames=list(NULL,c("intercept","slope")))
for(b in seq_len(B)){sampled<-sample.int(npig,npig,replace=TRUE);w<-tabulate(sampled,nbins=npig);dem[b,]<-deming_w(w)}
readr::write_csv(as.data.frame(dem),file.path("output","Deming_bootstrap_3000.csv"))
readr::write_csv(data.frame(estimate=deming_w(rep(1,npig))["slope"],CI_low=quantile(dem[,"slope"],.025),bootstrap_median=quantile(dem[,"slope"],.5),CI_high=quantile(dem[,"slope"],.975)),file.path("output","Deming_lambda1_final_summary.csv"))
# Threshold fast pig-cluster bootstrap.
pig_thr<-primary|>dplyr::group_by(canonical_pig_id)|>dplyr::summarise(low_TP=sum(ibp_same_map<70&nibp_map<70),low_FN=sum(ibp_same_map<70&nibp_map>=70),low_FP=sum(ibp_same_map>=70&nibp_map<70),low_TN=sum(ibp_same_map>=70&nibp_map>=70),h120_TP=sum(ibp_same_map>=120&nibp_map>=120),h120_FN=sum(ibp_same_map>=120&nibp_map<120),h120_FP=sum(ibp_same_map<120&nibp_map>=120),h120_TN=sum(ibp_same_map<120&nibp_map<120),h87_TP=sum(ibp_same_map>=120&nibp_map>=87),h87_FN=sum(ibp_same_map>=120&nibp_map<87),h87_FP=sum(ibp_same_map<120&nibp_map>=87),h87_TN=sum(ibp_same_map<120&nibp_map<87),.groups="drop")
calc4<-function(TP,FN,FP,TN)c(sensitivity=TP/(TP+FN),specificity=TN/(TN+FP),PPV=TP/(TP+FP),NPV=TN/(TN+FN))
set.seed(20261104);npig<-nrow(pig_thr);mat<-matrix(NA_real_,3000,12)
for(b in 1:3000){sampled<-sample.int(npig,npig,replace=TRUE);w<-tabulate(sampled,nbins=npig);s<-colSums(pig_thr[,-1]*w);mat[b,]<-c(calc4(s["low_TP"],s["low_FN"],s["low_FP"],s["low_TN"]),calc4(s["h120_TP"],s["h120_FN"],s["h120_FP"],s["h120_TN"]),calc4(s["h87_TP"],s["h87_FN"],s["h87_FP"],s["h87_TN"]))}
colnames(mat)<-c("low_sens","low_spec","low_PPV","low_NPV","high120_sens","high120_spec","high120_PPV","high120_NPV","high87_sens","high87_spec","high87_PPV","high87_NPV")
readr::write_csv(as.data.frame(mat),file.path("output","threshold_bootstrap_3000.csv"));ci<-apply(mat,2,quantile,c(.025,.5,.975),na.rm=TRUE);write.csv(ci,file.path("output","threshold_bootstrap_CI.csv"))
