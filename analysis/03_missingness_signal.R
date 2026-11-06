source(file.path("analysis","00_setup.R"))
reason <- function(s,m,d){if(is.na(s)&is.na(m)&is.na(d)) return("All IBP missing"); if(any(is.na(c(s,m,d)))) return("Partial IBP missing"); if(any(c(s,m,d)<=0)) return("Non-positive IBP"); if(!(s>m & m>d)) return("Ordering invalid"); "Valid"}
reasons<-mapply(reason,d$ibp_same_sys,d$ibp_same_map,d$ibp_same_dia)
readr::write_csv(as.data.frame(table(reasons)),file.path("output","same_row_IBP_availability_reasons.csv"))
p <- primary |> dplyr::mutate(diff=nibp_map-ibp_same_map,meanp=(nibp_map+ibp_same_map)/2,pp=ibp_same_sys-ibp_same_dia,
 pig_id=as.integer(factor(canonical_pig_id)),external_fail=(ibp_same_sys>200&pp<50)|(ibp_same_sys>150&ibp_same_sys<=200&pp<30)|(ibp_same_sys>=100&ibp_same_sys<=150&pp<15)|(ibp_same_sys<100&pp<10)|ibp_same_sys<=10|ibp_same_dia<=10|(ibp_same_sys==ibp_same_dia&ibp_same_dia==ibp_same_map))
fit_slope<-function(dat){dat<-dat|>dplyr::arrange(pig_id,session_path,source_row_index_0based);m<-geepack::geeglm(diff~I(meanp/10),id=pig_id,data=dat,family=gaussian,corstr="exchangeable");b<-coef(m)[2];se<-sqrt(diag(vcov(m)))[2];c(n=nrow(dat),pigs=dplyr::n_distinct(dat$pig_id),slope=b,lo=b-1.96*se,hi=b+1.96*se)}
sets<-list(primary=p,signal20=dplyr::filter(p,ibp_signal_quality_flag_20==1),signal10=dplyr::filter(p,ibp_signal_quality_flag_10==1),external=dplyr::filter(p,!external_fail),exclude_absdiff_gt120=dplyr::filter(p,abs(diff)<=120))
out<-do.call(rbind,lapply(names(sets),function(nm)c(analysis=nm,fit_slope(sets[[nm]]))))
readr::write_csv(as.data.frame(out),file.path("output","signal_extreme_sensitivity.csv"))
