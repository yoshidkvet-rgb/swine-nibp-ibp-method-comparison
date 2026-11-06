source(file.path("analysis","00_setup.R"))
p <- primary |>
  dplyr::mutate(diff=nibp_map-ibp_same_map, meanp=(nibp_map+ibp_same_map)/2,
                pig_id=as.integer(factor(canonical_pig_id))) |>
  dplyr::arrange(pig_id, session_path, source_row_index_0based)
m <- geepack::geeglm(diff~I(meanp/10), id=pig_id, data=p, family=gaussian, corstr="exchangeable")
b<-coef(m)[2]; se<-sqrt(diag(vcov(m)))[2]
out<-data.frame(slope_per_10mmHg=b,ci_low=b-1.96*se,ci_high=b+1.96*se,p_value=summary(m)$coefficients[2,4])
readr::write_csv(out,file.path("output","primary_pressure_slope.csv"))
