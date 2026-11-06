source(file.path("analysis","00_setup.R"))
base <- primary |>
 dplyr::mutate(diff=nibp_map-ibp_same_map, meanp=(nibp_map+ibp_same_map)/2,
               pig_id=as.integer(factor(canonical_pig_id)), session_id=as.integer(factor(session_path)))
p_pig <- base |> dplyr::arrange(pig_id,session_id,source_row_index_0based)
m_pig <- geepack::geeglm(diff~I(meanp/10),id=pig_id,data=p_pig,family=gaussian,corstr="exchangeable")
p_session <- base |> dplyr::arrange(session_id,event_number_in_session)
m_session <- geepack::geeglm(diff~I(meanp/10),id=session_id,data=p_session,family=gaussian,corstr="exchangeable")

# Fast Gaussian AR(1) GEE matching statsmodels Autoregressive(grid=TRUE) used in the frozen WP3 sensitivity.
p_ar <- p_session |> dplyr::mutate(x=meanp/10)
groups <- split(p_ar,p_ar$session_id)
X <- cbind(1,p_ar$x); y<-p_ar$diff
beta<-as.numeric(qr.solve(X,y)); rho<-0
for(iter in 1:100){
 lag0<-0;lag1<-0
 for(g in groups){r<-g$diff-(beta[1]+beta[2]*g$x);n<-length(r);if(n>1){lag1<-lag1+sum(r[-n]*r[-1])/(n-1);lag0<-lag0+sum(r^2)/n}}
 rho_new<-lag1/lag0; den<-1-rho_new^2; A<-matrix(0,2,2);bb<-numeric(2)
 for(g in groups){x<-g$x;yy<-g$diff;n<-length(x);if(n==1){qd<-1;qo<-numeric(0)}else{qd<-rep((1+rho_new^2)/den,n);qd[c(1,n)]<-1/den;qo<-rep(-rho_new/den,n-1)}
   A[1,1]<-A[1,1]+sum(qd)+2*sum(qo)
   if(n==1){A[1,2]<-A[1,2]+sum(qd*x);A[2,2]<-A[2,2]+sum(qd*x^2);bb[1]<-bb[1]+sum(qd*yy);bb[2]<-bb[2]+sum(qd*x*yy)} else {
     A[1,2]<-A[1,2]+sum(qd*x)+sum(qo*(x[-n]+x[-1])); A[2,2]<-A[2,2]+sum(qd*x^2)+2*sum(qo*x[-n]*x[-1]);
     bb[1]<-bb[1]+sum(qd*yy)+sum(qo*(yy[-n]+yy[-1])); bb[2]<-bb[2]+sum(qd*x*yy)+sum(qo*(x[-n]*yy[-1]+x[-1]*yy[-n]))}
 }
 A[2,1]<-A[1,2]; beta_new<-as.numeric(solve(A,bb)); if(max(abs(beta_new-beta))<1e-10 && abs(rho_new-rho)<1e-10){beta<-beta_new;rho<-rho_new;break};beta<-beta_new;rho<-rho_new
}
# within/between session decomposition; session-clustered
wb <- base |> dplyr::group_by(session_id) |>
 dplyr::mutate(session_mean_pressure=mean(meanp),within_pressure10=(meanp-session_mean_pressure)/10) |>
 dplyr::ungroup() |> dplyr::arrange(session_id,source_row_index_0based)
m_wb<-geepack::geeglm(diff~within_pressure10+I(session_mean_pressure/10),id=session_id,data=wb,family=gaussian,corstr="exchangeable")
out<-data.frame(analysis=c("animal_exchangeable","session_exchangeable","session_AR1","within_session","between_session"),
 value=c(coef(m_pig)[2],coef(m_session)[2],beta[2],coef(m_wb)[2],coef(m_wb)[3]),
 auxiliary=c(NA,NA,rho,NA,NA))
readr::write_csv(out,file.path("output","dependence_models.csv"))
