#getting model in right parameter space

devtools::load_all()
require(tidyverse)

init_EIR_vec <- c(12,14,16)
Q0_vec <- c(0.20, 0.82, 0.93)
bites_Bed_vec <- c(0.84, 0.95, 1)

acc_refusal_rate <- c(1-0.05, 1-0.21)
ivm_cov_vec <- c(0.609, 0.682, 0.747)

#getting lowest possible coverage and highest possible coverage of ivm
ivm_cov_vec_ref_adj <- c(ivm_cov_vec[1]*acc_refusal_rate[2], ivm_cov_vec[2], ivm_cov_vec[3])

years <- seq(2012, 2022, by = 1)
distr_years <- c(years[3], years[6], years[9])
distr_campaign <- (6*30)/365 #distributions are in June of the distr years

distr_campaign_days <- distr_campaign*365

#how many days into times is the first itn campaign?
time_period <- 365*13

net_seq <- seq(distr_campaign_days+(365*3), time_period, by = 3*365)
itn_on <- net_seq[1]

#the prevalence survey measuring ~14% qPCR all age prevalence and ~7.8% all age microscopy is
nov <- 30*11
y_2018 <- 365*7 # 8 years into simulation
prev_survey_date <- nov+y_2018

#put nets in the model: 55% phenotypic resistance acc Moss et al 2023
pyr_nets_df <- read.csv("C:/Users/nc1115/OneDrive - Imperial College London/PhD/PhD_malaria/data/ellie_net_efficacy/pyrethroid_only_nets.csv")
pyr_nets_res <- pyr_nets_df %>%
  filter(resistance == 0.55)

pyr_param_df <- data.frame(pyr_nets_res$dn0_med, pyr_nets_res$rn0_med, pyr_nets_res$gamman_med)
names(pyr_param_df) <- c("dn0", "rn0", "gamman")

mod_space_df <- expand.grid(init_EIR_vec, Q0_vec, bites_Bed_vec, ivm_cov_vec_ref_adj)
names(mod_space_df) <- c("init_EIR","Q0", "bites_Bed", "ivm_cov")

mod_space_df$dn0 <- pyr_param_df$dn0
mod_space_df$rn0 <- pyr_param_df$rn0
mod_space_df$gamman <- pyr_param_df$gamman
mod_space_df$itn_cov <- 0.75



mod_space_df <- mod_space_df %>%
  #filter(bites_Bed == bites_Bed_vec[2] & Q0 == Q0_vec[2]) %>%
  filter(bites_Bed == bites_Bed_vec[2])

mod_space_list <- list()

for(i in seq_len(nrow(mod_space_df))){
  mod_space_list[[i]] <- as.numeric(mod_space_df[i,])
}


#getting the ivermectin information
ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE)
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")
# Sourcing the extra functions required to generate the endectocide-specific parameters
source("R/mda_ivm_functions.R")



runfun <- function(mod_name){
  mod <- mod_name$generator$new(user= mod_name$state, use_dde = TRUE)
  modx <- mod$run(t = 1:time_period, tcrit = net_seq)
  op<- mod$transform_variables(modx)
  return(op)
}

#endectocide set up####
mda_int <- 30
eff_len <- 23
july <- 30*7
y_2021 <- 365*10
y_2022 <- 365*11
endec_y1_start <- july + y_2021 #start time of first MDA
endec_y2_start <- july+ y_2022

endec_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int,
                 endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)


#set IVM params for Hannah's mosq model with hazards#
ivm_parms1 <- ivRmectin::ivm_fun(
  #IVM_start_times = c(3120, 3150, 3180), #distribution every 3 months
  IVM_start_times = endec_start,
  time_period = time_period,
  hazard_profile = ivm_haz$IVM_300_3_HS[1:23],
  #hazard_profile = hazzy,
  ivm_coverage= 0.7, #gets updated later on!!
  ivm_min_age=5,
  ivm_max_age = 90)


eff_len <- 23


#########
#ento_space_mod <- function(itn_type_ivm_param){
#  init_EIR_in <- itn_type_ivm_param[1]
#  Q0_in  <- itn_type_ivm_param[2]
#  bites_Bed_in <- itn_type_ivm_param[3]
#  d_ITN0_in <- itn_type_ivm_param[4]
#  r_ITN0_in <- itn_type_ivm_param[5]
#  itn_half_life_in <- itn_type_ivm_param[6]
#  itn_cov_in <- itn_type_ivm_param[7]
#  output <- ivRmectin::create_r_model(
#    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
#    num_int = 2,
#    #num_int = 2, # number of vector control (IRS and ITN) population groups
#    #het_brackets = 5, # number of heterogeneous biting categories
#    #age = init_age, # the different age classes to be ran within the model
#    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
#    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
#    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
#    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
#    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
#    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
#    ivm_cov_par = 0, # proportion of popuulation receiving the endectocide
#    ivm_min_age = 5, # youngest age group receiving endectocide
#    ivm_max_age = 90, # oldest age group receiving endectocide
#    IVRM_start = ivm_parms1$IVRM_start,
#    Q0 = Q0_in,
#    bites_Bed = bites_Bed_in,
#    d_ITN0 = d_ITN0_in,
#    r_ITN0 = r_ITN0_in,
#    itn_half_life = itn_half_life_in*365,
#    itn_cov = itn_cov_in,
#    ITN_IRS_on = itn_on
#  )
#  return(output)
#}
#
#my_sim_ento_space_mod <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  out_list_antag <- lapply(mod_space_list, ento_space_mod)
#  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
#  out_df_antag <- do.call(rbind, sapply(1:(nrow(mod_space_df)), function(x){
#    df <- as.data.frame(res_out_antag[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
#                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
#  return(out_df_antag)
#
#}
#
#my_sim_ento_space_out <- my_sim_ento_space_mod()

#want to hit prev 7.8%

#ggplot(my_sim_ento_space_out, aes(x = t, y = slide_prev0to80))+
#  geom_line()+
#  facet_wrap(vars(ref))+
#  geom_vline(xintercept = prev_survey_date)
#
#my_sim_ento_space_out %>%
#  group_by(ref, Q0,bites_Bed) %>%
#  filter(t == prev_survey_date) %>%
#  filter(bites_Bed == bites_Bed_vec[2] & Q0 == Q0_vec[2]) %>%
#  summarise(prev = slide_prev0to80)
#
#target_prev <- 0.078
#
#get_prev <- function(model){
#  model %>%
#    filter(t == prev_survey_date) %>%
#    group_by(ref) %>%
#    summarise(prev = slide_prev0to80,
#              diff = (slide_prev0to80 - target_prev)^2)
#}
#
#get_prev(my_sim_ento_space_out)
#
#target_prev <- 0.078
#
## Function to compute squared difference for a given Q0, bites_Bed, and ref
#objective_fn <- function(ref, data = my_sim_ento_space_out, target_prev = 0.078) {
#  data %>%
#    dplyr::filter(ref == !!ref, t == prev_survey_date) %>%
#    dplyr::summarise(diff = (slide_prev0to80 - target_prev)^2, .groups = "drop") %>%
#    dplyr::pull(diff)
#}
#
## Find the best init_EIR (ref) within a given Q0–bites_Bed combo
#optimise_init_EIR <- function(Q0_in, bites_Bed_in, data = my_sim_ento_space_out) {
#  df_sub <- data %>%
#    dplyr::filter(Q0 == Q0_in, bites_Bed == bites_Bed_in, t == prev_survey_date) %>%
#    dplyr::mutate(diff = (slide_prev0to80 - target_prev)^2)
#
#  best <- df_sub %>%
#    dplyr::arrange(diff) %>%
#    dplyr::slice(1)
#
#  tibble(
#    Q0 = Q0_in,
#    bites_Bed = bites_Bed_in,
#    ref = best$ref,
#    best_prev = best$slide_prev0to80,
#    min_diff = best$diff
#  )
#}
#
## Run calibration grid
#param_grid <- expand.grid(Q0 = Q0_vec, bites_Bed = bites_Bed_vec)
#
#calibration_results <- param_grid %>%
#  purrr::pmap_dfr(~ optimise_init_EIR(..1, ..2, data = my_sim_ento_space_out))
#
#calibration_results
#
#best_ref <- calibration_results$ref
#
#init_EIR_vec[best_ref]
#
##so ref 14 and init_EIR_vec[2] give prev of 0.078 at the prev_survey_date

#ENDECTOCIDE MDA######

ivm_params_df <- mod_space_df %>%
  filter(init_EIR == 14) #may need to update this to init_EIR_vec[best_ref]]




######################
#ivm_cov_par <- ivm_parms1$ivm_cov_par*0.95 #if 5% refusal rate. or by ~0.8 if 20% refusal rate

#ivm_cov = ivm_cov_par*(exp(-ivm_parms1$ivm_min_age/21) - exp(-ivm_parms1$ivm_max_age/21))

#######################
#so 70% amongst eligibles goes to 54% in total population, or 51% w 5% refusal (for the median)


ivm_parms_list <- list()

for(i in seq_len(nrow(ivm_params_df))){
  ivm_parms_list[[i]] <- as.numeric(ivm_params_df[i,])
}




#then set up the endectocide model (model A) with the MDAs in 2021 and 2022


ivm_mda_mod <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in  <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  ivm_cov_in <- itn_type_ivm_param[4]
  d_ITN0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  itn_cov_in <- itn_type_ivm_param[8]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_endectocide.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = ivm_cov_in, # proportion of population receiving the endectocide
    ivm_min_age = 5, # youngest age group receiving endectocide
    ivm_max_age = 90, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in*365,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on
  )
  return(output)
}

my_sim_ivm_mda_mod <- function(){
  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
  out_list_antag <- lapply(ivm_parms_list, ivm_mda_mod)
  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_df)), function(x){
    df <- as.data.frame(res_out_antag[[x]])
    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, Ivtot))
    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "explore", int = "ITN-only"))}, simplify = F))
  return(out_df_antag)

}

my_sim_ivm_mda_out <- my_sim_ivm_mda_mod()
saveRDS(my_sim_ivm_mda_out, file = "analysis/matamal-runs/matamal-Q0-ivm-cov-combos.rds")

ggplot(my_sim_ivm_mda_out, aes(x = t/365, y = slide_prev0to80))+
  geom_line()
ggplot(my_sim_ivm_mda_out, aes(x = t/365, y = Ivtot))+
  geom_line()

##what value of endec_mu and wane give same total number of infectious vectors over the ivm killing period?

#need to run model C

wane_vec <- seq(0, 0.1, 0.01)
endec_mu_vec <- seq(0.05, 0.15, 0.01)

malsim_odin <- expand.grid(wane_vec, endec_mu_vec)
names(malsim_odin) <- c("wane", "endec_mu")


ivm_params_exp_decay_df <- dplyr::cross_join(ivm_params_df, malsim_odin)

ivm_params_exp_decay_list <- list()

#just one row
#ivm_params_exp_decay_df <- ivm_params_exp_decay_df[1:2,]

for(i in seq_len(nrow(ivm_params_exp_decay_df))){
  ivm_params_exp_decay_list[[i]] <- as.numeric(ivm_params_exp_decay_df[i,])
}


ivm_mda_mod_malsim <- function(itn_type_ivm_param){
  init_EIR_in <- itn_type_ivm_param[1]
  Q0_in  <- itn_type_ivm_param[2]
  bites_Bed_in <- itn_type_ivm_param[3]
  #ivm_cov_in <- itn_type_ivm_params[4]
  d_ITN0_in <- itn_type_ivm_param[5]
  r_ITN0_in <- itn_type_ivm_param[6]
  itn_half_life_in <- itn_type_ivm_param[7]
  itn_cov_in <- itn_type_ivm_param[8]
  wane_in <- itn_type_ivm_param[9]
  endec_mu_in <- itn_type_ivm_param[10]
  output <- ivRmectin::create_r_model(
    odin_model_path = system.file("extdata/odin_model_malariasim_exp_decay.R", package = "ivRmectin"),
    num_int = 2,
    #num_int = 2, # number of vector control (IRS and ITN) population groups
    #het_brackets = 5, # number of heterogeneous biting categories
    #age = init_age, # the different age classes to be ran within the model
    init_EIR = init_EIR_in, # the Entomological Innoculation Rate
    #country = "Senegal", # Country setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    #admin2 = "Fatick", # Admin 2 setting to be run - see admin_units_seasonal.rds in inst/extdata for more info
    ttt = ivm_parms1$ttt, # model specific parameter to control timing of endectocide delivery
    eff_len = ivm_parms1$eff_len, # number of days after receiving endectocide that HR is higher
    haz = ivm_parms1$haz, # hazard ratio for each off the eff_len number of days
    ivm_cov_par = 0.687, # redundant param in this model
    ivm_min_age = 5, # youngest age group receiving endectocide
    ivm_max_age = 90, # oldest age group receiving endectocide
    IVRM_start = ivm_parms1$IVRM_start,
    Q0 = Q0_in,
    bites_Bed = bites_Bed_in,
    d_ITN0 = d_ITN0_in,
    r_ITN0 = r_ITN0_in,
    itn_half_life = itn_half_life_in*365,
    itn_cov = itn_cov_in,
    ITN_IRS_on = itn_on,
    endec_mu = endec_mu_in,
    wane = wane_in
  )
  return(output)
}

#my_sim_ivm_mda_mod_malsim <- function(){
#  #pyr_out_list_antag_ITN <- purrr::map2(y, x, antag_ITN_cov_loop) #loop through all parameter values
#  out_list_antag <- lapply(ivm_params_exp_decay_list, ivm_mda_mod_malsim)
#  res_out_antag <- lapply(out_list_antag, runfun) #put these values into the model
#  out_df_antag <- do.call(rbind, sapply(1:(nrow(ivm_params_exp_decay_df)), function(x){
#    df <- as.data.frame(res_out_antag[[x]])
#    df2 <- as.data.frame(dplyr::select(.data = df,t, mu, mv, avhc, itn_cov, EIR_tot, slide_prev0to5, slide_prev0to80,
#                                       Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, Ivtot, wane, endec_mu))
#    df3 <- as.data.frame(dplyr::mutate(.data = df2, ref = x, scenario = "odin-malsim", int = "ITN-IVM"))}, simplify = F))
#  return(out_df_antag)
#
#}

my_sim_ivm_mda_mod_malsim <- function() {

  results_list <- vector("list", length(ivm_params_exp_decay_list))

  for (i in seq_along(ivm_params_exp_decay_list)) {
    params <- ivm_params_exp_decay_list[[i]]

    # Create and run model for this parameter set
    model_obj <- ivm_mda_mod_malsim(params)
    res <- runfun(model_obj)   # your simulation call

    # Keep only relevant variables (convert to small data frame)
    res_df <- as.data.frame(res) %>%
      dplyr::select(
        t, mu, mv, avhc, itn_cov, slide_prev0to5, slide_prev0to80,
        Q0, IVRM_sr, EIRout, clin_inc0to5, bites_Bed, Ivtot, wane, endec_mu
      ) %>%
      dplyr::mutate(ref = i, scenario = "odin-malsim", int = "ITN-IVM")

    results_list[[i]] <- res_df

    # Cleanup
    rm(model_obj, res, res_df)
    gc(verbose = FALSE)

    # Optional: write out every 50 runs to disk
    if (i %% 50 == 0) {
      saveRDS(results_list[1:i],
              file = file.path("analysis/matamal-runs/results_ivm_exp_decay", paste0("ivm_partial_", i, ".rds")))
    }
  }

  # Combine all results at the end
  out_df_antag <- dplyr::bind_rows(results_list)
  return(out_df_antag)
}

#need to do the last from 1050 I think because of the modulus thing.
#and bind that one


my_sim_ivm_mda_out_malsim <- my_sim_ivm_mda_mod_malsim()

saveRDS(my_sim_ivm_mda_out_malsim, file = "analysis/matamal-runs/matamal-Q0-ivm-cov-combos-malariasim-odin.rds")

length(unique(my_sim_ivm_mda_out_malsim$ref))

#beepr::beep(sound = 1)


odin_model <- my_sim_ivm_mda_out
exp_decay_model <- my_sim_ivm_mda_out_malsim

endec_t_start <- endec_start[1]
endec_t_end <- endec_start[6]+23

ivm_params_df
#ref 1-3 are low ivm cov, 4-6 are medium cov and 7-9 are high cov
#1,4,7 are low Q0, 2,5,8 are medium Q0 and 3,6,9 are high Q0

# Prepare odin model times split by scenario
# Prepare odin model times split by scenario
times_odin_model <- odin_model %>%
  mutate(ivm_cov = case_when(
    ref %in% c(1,2,3) ~ "low_cov",
    ref %in% c(4,5,6) ~ "medium_cov",
    ref %in% c(7,8,9) ~ "high_cov"
  )) %>%
  filter(between(t, endec_t_start, endec_t_end)) %>%
  select(t, Ivtot, ivm_cov, Q0)

# Get unique combinations of ivm_cov and Q0
unique_scenarios <- unique(times_odin_model[, c("ivm_cov", "Q0")])

# Split exp decay model by Q0
time_exp_decay_model <- exp_decay_model %>%
  filter(between(t, endec_t_start, endec_t_end)) %>%
  select(t, Ivtot, Q0, ref)

##Low Q0 and low coverage#####
Q0_vals <- unique(times_odin_model$Q0)
ivm_cov_vals <- unique(times_odin_model$ivm_cov)


exp_decay_low_Q0 <- time_exp_decay_model %>%
  filter(Q0 == Q0_vals[1])

odin_low_cov_low_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[1] & ivm_cov == ivm_cov_vals[1])

list_exp_decay_model_low_Q0 <- split(exp_decay_low_Q0, f = exp_decay_low_Q0$ref)

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_low_cov_low_Q0$Ivtot - list_exp_decay_model_low_Q0[[i]]$Ivtot)^2))
}
range(error)
error_low_cov_low_Q0 <- min(error)
index <- which.min(error)
best_fit_low_cov_low_Q0 <- malsim_odin[index,]

#low Q0 and medium coverage
odin_medium_cov_low_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[1] & ivm_cov == ivm_cov_vals[2])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_medium_cov_low_Q0$Ivtot - list_exp_decay_model_low_Q0[[i]]$Ivtot)^2))
}
range(error)
error_medium_cov_low_Q0 <- min(error)
index <- which.min(error)
best_fit_medium_cov_low_Q0 <- malsim_odin[index,]

#low Q0 and high coverage
odin_high_cov_low_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[1] & ivm_cov == ivm_cov_vals[3])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_high_cov_low_Q0$Ivtot - list_exp_decay_model_low_Q0[[i]]$Ivtot)^2))
}
range(error)
error_high_cov_low_Q0 <- min(error)
index <- which.min(error)
best_fit_high_cov_low_Q0 <- malsim_odin[index,]

#medium Q0 and low coverage #####

exp_decay_medium_Q0 <- time_exp_decay_model %>%
  filter(Q0 == Q0_vals[2])

list_exp_decay_model_medium_Q0 <- split(exp_decay_medium_Q0, f = exp_decay_medium_Q0$ref)

#copy
odin_low_cov_medium_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[2] & ivm_cov == ivm_cov_vals[1])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_low_cov_medium_Q0$Ivtot - list_exp_decay_model_medium_Q0[[i]]$Ivtot)^2))
}
range(error)
error_low_cov_medium_Q0 <- min(error)
index <- which.min(error)
best_fit_low_cov_medium_Q0 <- malsim_odin[index,]

#medium Q0 and medium cov

odin_medium_cov_medium_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[2] & ivm_cov == ivm_cov_vals[2])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_medium_cov_medium_Q0$Ivtot - list_exp_decay_model_medium_Q0[[i]]$Ivtot)^2))
}
range(error)
error_medium_cov_medium_Q0 <- min(error)
index <- which.min(error)
best_fit_medium_cov_medium_Q0 <- malsim_odin[index,]

#medium Q0 and high cov
odin_high_cov_medium_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[2] & ivm_cov == ivm_cov_vals[3])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_high_cov_medium_Q0$Ivtot - list_exp_decay_model_medium_Q0[[i]]$Ivtot)^2))
}
range(error)
error_high_cov_medium_Q0 <- min(error)
index <- which.min(error)
best_fit_high_cov_medium_Q0 <- malsim_odin[index,]

#high Q0 and low cov####
exp_decay_high_Q0 <- time_exp_decay_model %>%
  filter(Q0 == Q0_vals[3])

list_exp_decay_model_high_Q0 <- split(exp_decay_high_Q0, f = exp_decay_high_Q0$ref)

#copy
odin_low_cov_high_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[3] & ivm_cov == ivm_cov_vals[1])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_low_cov_high_Q0$Ivtot - list_exp_decay_model_high_Q0[[i]]$Ivtot)^2))
}
range(error)
error_low_cov_high_Q0 <- min(error)
index <- which.min(error)
best_fit_low_cov_high_Q0 <- malsim_odin[index,]

#medium cov and high Q0
odin_medium_cov_high_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[3] & ivm_cov == ivm_cov_vals[2])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_medium_cov_high_Q0$Ivtot - list_exp_decay_model_high_Q0[[i]]$Ivtot)^2))
}
range(error)
error_medium_cov_high_Q0 <- min(error)
index <- which.min(error)
best_fit_medium_cov_high_Q0 <- malsim_odin[index,]

#high cov and high Q0
odin_high_cov_high_Q0 <- times_odin_model %>%
  filter(Q0 == Q0_vals[3] & ivm_cov == ivm_cov_vals[3])

error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((odin_high_cov_high_Q0$Ivtot - list_exp_decay_model_high_Q0[[i]]$Ivtot)^2))
}
range(error)
error_high_cov_high_Q0 <- min(error)
index <- which.min(error)
best_fit_high_cov_high_Q0 <- malsim_odin[index,]

##then bind all together
scenario <- c("low_Q0_low_cov", "low_Q0_medium_cov", "low_Q0_high_cov",
              "medium_Q0_low_cov", "medium_Q0_medium_cov", "medium_Q0_high_cov",
              "high_Q0_low_cov", "high_Q0_medium_cov", "high_Q0_high_cov")

errors <- c(error_low_cov_low_Q0, error_low_cov_medium_Q0, error_low_cov_high_Q0,
            error_medium_cov_low_Q0, error_medium_cov_medium_Q0, error_medium_cov_high_Q0,
            error_high_cov_low_Q0, error_high_cov_medium_Q0, error_high_cov_high_Q0
            )

best_wane <- c(best_fit_low_cov_low_Q0$wane, best_fit_low_cov_medium_Q0$wane, best_fit_low_cov_high_Q0$wane,
               best_fit_medium_cov_low_Q0$wane, best_fit_medium_cov_medium_Q0$wane, best_fit_medium_cov_high_Q0$wane,
               best_fit_high_cov_low_Q0$wane, best_fit_high_cov_medium_Q0$wane, best_fit_high_cov_high_Q0$wane)

best_endec_mu <- c(best_fit_low_cov_low_Q0$endec_mu, best_fit_low_cov_medium_Q0$endec_mu, best_fit_low_cov_high_Q0$endec_mu,
               best_fit_medium_cov_low_Q0$endec_mu, best_fit_medium_cov_medium_Q0$endec_mu, best_fit_medium_cov_high_Q0$endec_mu,
               best_fit_high_cov_low_Q0$endec_mu, best_fit_high_cov_medium_Q0$endec_mu, best_fit_high_cov_high_Q0$endec_mu)

malariasim_endec_killing_inputs <- data.frame(scenario, errors, best_wane, best_endec_mu)

saveRDS(malariasim_endec_killing_inputs, file = "analysis/matamal-runs/malariasim_endec_killing_inputs.rds")
#need to generate plots: pipe through the exp_decay model and illustrate similar fit
#and impact of incidence, prevalence etc





###############################################################

times_odin_model <- odin_model %>% ##assign an id by group_by Q0 and ivm_cov
  filter(between(t, endec_t_start, endec_t_end)) %>%
  select(t, Ivtot)

time_exp_decay_model <- exp_decay_model %>%
  filter(between(t, endec_t_start, endec_t_end)) %>%
  select(t, Ivtot, ref)

list_exp_decay_model <- split(time_exp_decay_model, f = time_exp_decay_model$ref)
error <- numeric()
for(i in 1:nrow(malsim_odin)){
  #error <- c(error, sum(times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2)
  error <- c(error, sum((times_odin_model$Ivtot - list_exp_decay_model[[i]]$Ivtot)^2))
}
range(error)
index <- which.min(error)
best_fit <- malsim_odin[index,]

exp_decay_model %>%
  filter(wane == 0.02 & endec_mu == 0.05) %>%
  ggplot()+
  aes(x = t, y = Ivtot)+
  geom_line()+
  geom_line(data = odin_model, aes(x = t, y = Ivtot), col = "red")


