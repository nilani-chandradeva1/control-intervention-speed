#matamal modelling

#use Okell paper to triangulate estimates of prevalence from microscopy, PCR and RDT

input_micro <- 0.1
log_odds_micr_all_age <- qlogis(input_micro)
log_odds_RDT_prev = 0.108 + 0.907*log_odds_micr_all_age
RDT_prev <- plogis(log_odds_RDT_prev)


log_odds_PCR_prev = 0.108 + 0.907 * ((log_odds_RDT_prev-0.954)/0.868)
plogis(log_odds_PCR_prev)

#then

df_odin <- readRDS("C:/Users/nc1115/Documents/github/matamal-cluster-modelling/2.stat-analysis-model-params/output/df_odin_inputs.rds")


log_odds_PCR_prev <- qlogis(0.5) # get on log odds scale
log_odds_micr_all_age <- (log_odds_PCR_prev + 0.777) / 1.239
log_odds_micr_all_age
prev_micr <- plogis(log_odds_micr_all_age)
