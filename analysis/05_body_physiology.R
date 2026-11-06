source(file.path("analysis","00_setup.R"))
bw_data <- primary |> dplyr::mutate(diff=nibp_map-ibp_same_map,meanp=(nibp_map+ibp_same_map)/2,mean_c10=(meanp-mean(meanp))/10,bw10=body_weight_kg_registry/10,year_f=factor(year),pig_id=as.integer(factor(canonical_pig_id))) |> dplyr::arrange(pig_id,session_path,source_row_index_0based)
bw<-geepack::geeglm(diff~mean_c10*bw10+year_f,id=pig_id,data=bw_data,family=gaussian,corstr="exchangeable")
readr::write_csv(data.frame(term=names(coef(bw)),estimate=coef(bw),se=sqrt(diag(vcov(bw))),p=summary(bw)$coefficients[,4]),file.path("output","body_weight_model.csv"))
if(file.exists(wp7_path)){
 w<-readr::read_csv(wp7_path,show_col_types=FALSE)|>dplyr::mutate(diff=nibp_map-ibp_same_map,meanp=(nibp_map+ibp_same_map)/2,session_id=as.integer(factor(session_path)))
 plausible<-w|>dplyr::filter(HR>0,HR<=250,EtCO2>0,EtCO2<=80,expAG>=0.3,expAG<=4,temperature1>=30,temperature1<=42)|>dplyr::mutate(HR_c10=(HR-median(HR))/10,EtCO2_c5=(EtCO2-median(EtCO2))/5,expAG_c=expAG-median(expAG),temp_c=temperature1-median(temperature1),time_c_h=time_since_recording_start_h-median(time_since_recording_start_h))|>dplyr::arrange(session_id,source_row_index_0based)
 m1<-geepack::geeglm(diff~I(meanp/10)+HR_c10+EtCO2_c5+expAG_c+temp_c+time_c_h,id=session_id,data=plausible,family=gaussian,corstr="exchangeable")
 m2<-geepack::geeglm(diff~I(meanp/10)+HR_c10+EtCO2_c5+expAG_c+temp_c+time_c_h+factor(source_program_group),id=session_id,data=plausible,family=gaussian,corstr="exchangeable")
 readr::write_csv(data.frame(model=c("physiology","physiology_plus_source_program"),pressure_slope=c(coef(m1)[2],coef(m2)[2])),file.path("output","physiology_protocol_adjusted.csv"))
}
