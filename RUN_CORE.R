scripts <- c("01_primary_pressure.R","02_dependence_models.R","03_missingness_signal.R","04_thresholds.R","05_body_physiology.R")
for(s in scripts){message("Running ",s);source(file.path("analysis",s),local=new.env(parent=globalenv()))}
writeLines(capture.output(sessionInfo()),file.path("output","sessionInfo.txt"))
