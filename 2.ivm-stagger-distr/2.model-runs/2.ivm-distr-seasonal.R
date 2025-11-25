#model staggered and overnight distributions in seasonal setting
#Fatick, Senegal: highly seasonal setting
require(tidyverse)
devtools::load_all()

#we now pass this into the model

time_period <- 365*15 #long run to get to eqm
mda_int <- 30 #MDA starts every 30 days (for 3 months)


runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, rtol = 1e-8, atol = 1e-6)
  op<- mod$transform_variables(modx)
  return(op)
}

#distr plan
start <- (365*5)+200 #time for starting distr, help get to reach eqm, and to distr at right time when seasonality introduced

#parameters
Q0_in <- 0.92 #Human Blood Index of An. gambiae-like vector
init_EIR_in <- c(2, 100) #initial EIRs (low and high endemicity)
target_cov <- c(0.5, 0.7, 0.9) #target coverage of ivermectin

#create parameter set
df_var_all <- expand.grid(Q0 = Q0_in, ivm_cov = target_cov, init_EIR = init_EIR_in)

my_list_all <- list()
for (i in seq_len(nrow(df_var_all))){
  my_list_all[[i]] <- as.numeric(df_var_all[i,])
}


#we model the overnight distributions in both the original model (odin_model_endectocide) and model with staggered distributions (odin_model_endectocide_staggered_HS) to ensure that they are consistent
#for manuscript, we then simulate overnight distributions in the staggered distribution model
#this helps to ensure that any differences between the overnight and staggered distributions are due to the distribution alone

#read in HR curves
df_all <- readRDS("2.ivm-stagger-distr/output/HR_overnight.rds")
df_extended_10d_new_HR <- readRDS("2.ivm-stagger-distr/output/HR_10d_stagger.rds")
df_extended_20d_new_HR <- readRDS("2.ivm-stagger-distr/output/HR_20d_stagger.rds")

#first, set up parameters for modelling overnight distributions in the original model (odin_model_endectocide)
#all-in-one
ivm_parms_all <- ivRmectin::ivm_fun(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = c(start, start + 30, start + 60),
  #IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_all$HR_use_above1[1:23], #pass in HR 1 to 23
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

ivm_parms_all$IVRM_start[start:2055]

mod_all <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all$ttt,
    eff_len = ivm_parms_all$eff_len,
    haz = ivm_parms_all$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all$ivm_min_age,
    ivm_max_age = ivm_parms_all$ivm_max_age,
    IVRM_start = ivm_parms_all$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all <- function(){
  mod_out_list <- lapply(my_list_all, mod_all)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov,ivm_cov_par, slide_prev0to80, Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all <- my_sim_mod_all()


#checking against the all-in-one in the staggered model
ivm_parms_all_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  time_period = time_period,
  hazard_profile = df_all$HR_use_above1 , #pass in HR for
  prop_human_HR_threshold = df_all$prop_pop_cov, #still passing into ivm compartment but experience HR of 1
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

plot(df_all$Day, df_all$prop_pop_cov) #steps in the coverage. Goes to 0 when HR = 1

ivm_parms_all_stag$IVRM_start[2025]
ivm_parms_all$IVRM_start[2024]
ivm_parms_all_stag$prop_human_HR_threshold[2025]



mod_all_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all_stag$ttt,
    eff_len = ivm_parms_all_stag$eff_len,
    prop_human_HR_threshold = ivm_parms_all_stag$prop_human_HR_threshold,
    haz = ivm_parms_all_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_all_stag$ivm_min_age,
    ivm_max_age = ivm_parms_all_stag$ivm_max_age,
    IVRM_start = ivm_parms_all_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_all_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_all_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov,ivm_cov_par ,slide_prev0to80,Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "all-in-one-stag-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_all_stag <- my_sim_mod_all_stag()

checking <- df_mod_all_stag %>%
  filter(ref == 2)


ivm_parms_all_stag$IVRM_start[2025:2085]
checking$ivm_cov[2024:2055]
ivm_parms_all_stag$haz[2025:2055]
ivm_parms_all_stag$haz
checking$mv[2025:2048]

#compare treat all stag and normal code (i.e. not actually any staggering.)
df_mod_compare <- rbind(df_mod_all, df_mod_all_stag)

df_mod_compare <- df_mod_compare %>%
  mutate(init_EIR = case_when(ref == 1 ~ 2,
                              TRUE ~ 100)) %>%
  filter(init_EIR == 100)

check_mv <- df_mod_compare %>%
  filter(between(t, 2025, 2050)) %>%
  group_by(model_type) %>%
  summarise(t = t,
            mv = mv)%>%
  pivot_wider(names_from = model_type, values_from = mv)

df_mod_compare %>%
  #filter(between(t, 2025, 2050)) %>%
  group_by(model_type) %>%
  summarise(
    tot_mv = sum(mv))

#expect mv to match, but they do not: small discrepancy when ivermectin is turned off (from day 24 after start)
#I think driven by how the solver handles jumps - prop_human_HR is suddenly changing to 0 again?

unique(df_mod_compare$ref)

df_mod_compare %>%
  filter(ivm_cov_par == 0.7 & ref == 5) %>%
  ggplot(aes(x = t, y = mv, col = as.factor(model_type)))+
  geom_line()+
  xlim(2000, 2500)+
  theme_minimal()

df_mod_compare %>%
  group_by(model_type) %>%
  summarise(tot_mosq= sum(mv)) #marginal difference in total number of mosquitoes.

#then for the 10d and 20d distributions
ivm_parms_10d_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  prop_human_HR_threshold = df_extended_10d_new_HR$prop_pop_cov,
  time_period = time_period,
  hazard_profile = df_extended_10d_new_HR$HR_use_above1,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

range(ivm_parms_10d_stag$prop_human_HR_threshold)
plot(ivm_parms_10d_stag$prop_human_HR_threshold)


mod_10d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_10d_stag$ttt,
    eff_len = ivm_parms_10d_stag$eff_len,
    haz = ivm_parms_10d_stag$haz,
    prop_human_HR_threshold = ivm_parms_10d_stag$prop_human_HR_threshold,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_10d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_10d_stag$ivm_max_age,
    IVRM_start = ivm_parms_10d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_10d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_10d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, ivm_cov_par,slide_prev0to80,Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "10d-stagger-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_10d_stag <- my_sim_mod_10d_stag()


#then 20d stagger
ivm_parms_20d_stag <- ivRmectin::ivm_fun_stag_cov(#IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = start,
  prop_human_HR_threshold = df_extended_20d_new_HR$prop_pop_cov,
  time_period = time_period,
  hazard_profile = df_extended_20d_new_HR$HR_use_above1,
  #hazard_profile = hazzy,
  ivm_coverage=0.8, #gets updated in function
  ivm_min_age=5,
  ivm_max_age = 90)

mod_20d_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_20d_stag$ttt,
    prop_human_HR_threshold = ivm_parms_20d_stag$prop_human_HR_threshold,
    eff_len = ivm_parms_20d_stag$eff_len,
    haz = ivm_parms_20d_stag$haz,
    ivm_cov_par = ivm_cov_in,
    ivm_min_age = ivm_parms_20d_stag$ivm_min_age,
    ivm_max_age = ivm_parms_20d_stag$ivm_max_age,
    IVRM_start = ivm_parms_20d_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_20d_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_20d_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, ivm_cov_par,slide_prev0to80,Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "20d-stagger-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_20d_stag <- my_sim_mod_20d_stag()

#baseline
mod_baseline_stag <-  function(data_in){
  Q0_in <- data_in[1]
  ivm_cov_in <- data_in[2]
  init_EIR_in <- data_in[3]
  output <- ivRmectin:::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide_staggered_HS.R", package = "ivRmectin"),
    num_int = 1,
    #num_int = 2,
    #ITN_IRS_on = 100,
    #itn_cov = 0.75,
    #het_brackets = 5,
    #age = init_age,
    init_EIR = init_EIR_in,
    country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    admin2 = "Fatick",
    ttt = ivm_parms_all_stag$ttt,
    eff_len = ivm_parms_all_stag$eff_len,
    prop_human_HR_threshold = ivm_parms_all_stag$prop_human_HR_threshold,
    haz = ivm_parms_all_stag$haz,
    ivm_cov_par = 0, #set coverage to 0 - no impact
    ivm_min_age = ivm_parms_all_stag$ivm_min_age,
    ivm_max_age = ivm_parms_all_stag$ivm_max_age,
    IVRM_start = ivm_parms_all_stag$IVRM_start,
    Q0 = Q0_in
  )
  return(output)
}

my_sim_mod_baseline_stag <- function(){
  mod_out_list <- lapply(my_list_all, mod_baseline_stag)
  res_mod_out <- lapply(mod_out_list, runfun)
  mod_df <- do.call(rbind, sapply(1:(nrow(df_var_all)), function(x){
    df <- as.data.frame(res_mod_out[[x]])
    df2 <-  as.data.frame(dplyr::select(.data = df, t, mu, mv, mvx_dead,Q0, ivm_cov, ivm_cov_par,slide_prev0to80,Svtot, Evtot, Ivtot,EIR_tot, EIRout, clin_inc0to5))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, model_type = "baseline-Sen"))}, simplify = F))
  return(mod_df)
} #adding mvtot_1 and 2 and 3 so can rbind onto the rest

df_mod_baseline_stag <- my_sim_mod_baseline_stag()

df_mod_distr <- do.call("rbind", list(df_mod_10d_stag, df_mod_20d_stag, df_mod_all_stag, df_mod_all, df_mod_baseline_stag))

write_rds(df_mod_distr, file = "2.ivm-stagger-distr/output/df_distr_HR_3m_seasonal.rds")

