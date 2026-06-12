#' Create mass drug administration (MDA) intervention settings
#'
#' Constructs a list of parameters defining a mass drug administration
#' (MDA) intervention for use in the transmission model. User-specified
#' age limits are mapped to the nearest predefined age categories, and
#' up to nine treatment rounds are supported.
#'
#' @param ncc Integer indicating the number of consecutive treatment
#'   cycles or campaign components. Default is `2`.
#' @param MDA_times Numeric vector giving the times at which MDA rounds
#'   are administered. A maximum of nine rounds is supported.
#' @param MDA_grp_prop Proportion of the population belonging to the
#'   target treatment group. Default is `0.9`.
#' @param MDA_cov Proportion of the eligible target group receiving
#'   treatment (coverage). Default is `0.8`.
#' @param MDA_st_age Minimum eligible age (in years) for treatment.
#'   This value is mapped to the nearest predefined age category.
#'   Default is `0.5`.
#' @param MDA_en_age Maximum eligible age (in years) for treatment.
#'   This value is mapped to the nearest predefined age category.
#'   Default is `80`.
#' @param MDA_succ Probability that treatment is successful among
#'   treated individuals. Default is `0.95`.
#' @param MDA_drug_choice Integer identifying the drug regimen to use.
#'   The interpretation of this code depends on the model implementation.
#'
#' @return
#' A named list containing:
#' \describe{
#'   \item{MDA_t1--MDA_t9}{Times of up to nine MDA rounds. Unused rounds
#'   are filled with a large placeholder value (`1e5`).}
#'   \item{MDA_st_cat}{Index of the lower eligible age category.}
#'   \item{MDA_en_cat}{Index of the upper eligible age category.}
#'   \item{ncc}{Number of treatment cycles/components.}
#'   \item{MDA_grp_prop}{Target group proportion.}
#'   \item{MDA_cov}{Treatment coverage.}
#'   \item{MDA_succ}{Treatment success probability.}
#'   \item{mda_drug_choice}{Drug regimen identifier.}
#' }
#'
#' @export
mda_fun <- function(ncc = 2,
                    MDA_times = NULL,
                    MDA_grp_prop = 0.9,
                    MDA_cov = 0.8,
                    MDA_st_age = 0.5,
                    MDA_en_age = 80,
                    MDA_succ = 0.95,
                    MDA_drug_choice = 2
){

  # set up MDA times
  L1 = as.list(MDA_times)
  if(length(MDA_times) < 9) L1 = c(L1, as.list(rep(1e5, 9 - length(MDA_times))))
  names(L1) = paste0("MDA_t",1:9)
  if(length(MDA_times) >9)print("too many MDA times")

  # translate the user defined ages in to age categories
  MDA_st_cat = which.min(abs(init_age - MDA_st_age))
  MDA_en_cat = which.min(abs(init_age - MDA_en_age))
  if(!MDA_st_age %in% init_age) print(paste("Not a pre-defined age break point.",init_age[MDA_st_cat],"used instead"))
  if(!MDA_en_age %in% init_age) print(paste("Not a pre-defined age break point.",init_age[MDA_en_cat],"used instead"))

  op = c(L1,list(
    MDA_st_cat = MDA_st_cat,
    MDA_en_cat = MDA_en_cat,
    ncc = ncc,
    MDA_grp_prop = MDA_grp_prop,
    MDA_cov = MDA_cov,
    MDA_succ = MDA_succ,
    mda_drug_choice = MDA_drug_choice))

  return(op)

}

#' Ivermectin intervention function
#'
#' @param IVM_start_times Numeric vector of ivermectin start times.
#' @param time_period Numeric time period.
#' @param hazard_profile Hazard profile used in the simulation.
#' @param ivm_coverage Ivermectin coverage.
#' @param ivm_min_age Minimum eligible age.
#' @param ivm_max_age Maximum eligible age.
#' @param bites_per_3_days Number of bites per three days.
#'
#' @return Description of return value.
#' @export
ivm_fun <- function(IVM_start_times,
                    time_period,
                    hazard_profile,
                    ivm_coverage = 0.8,
                    ivm_min_age = 5,
                    ivm_max_age = 200,
                    bites_per_3_days = 1
){

  # function to make the ivermectin time profile
  make_IVRM = function(IVRM_st, eff_len,ttt){
    IVRM =  rep(max(ttt), length(ttt))
    for(i in 1:length(IVRM_st)){
      IVRM[ttt >= IVRM_st[i]-1 & ttt < IVRM_st[i] + eff_len -1] = IVRM_st[i]
    }
    return(IVRM)
  }

  ttt = 0:time_period
  eff_len = length(hazard_profile)
  IVRM_start = make_IVRM(IVM_start_times, eff_len, ttt)

  op = list(ttt = ttt,
            eff_len = eff_len,
            haz = hazard_profile,
            ivm_cov_par = ivm_coverage,
            ivm_min_age = ivm_min_age,
            ivm_max_age = ivm_max_age,
            B2 = bites_per_3_days,
            IVRM_start=IVRM_start
  )
  return(op)

}

#' Create ivermectin intervention settings with staggered high-risk coverage
#'
#' Constructs a time profile for an ivermectin intervention, including the
#' intervention start times, hazard profile, age eligibility limits, biting
#' rate parameter, and a staggered time-varying profile for the proportion
#' of humans above the high-risk threshold.
#'
#' @param IVM_start_times Numeric vector giving the start times of ivermectin
#'   intervention rounds.
#' @param time_period Numeric value giving the total simulation time period.
#'   The returned time vector runs from `0` to `time_period`.
#' @param hazard_profile Numeric vector giving the ivermectin hazard profile
#'   over the duration of the intervention effect.
#' @param prop_human_HR_threshold Numeric vector giving the time-varying
#'   proportion of humans above the high-risk threshold during each
#'   intervention period. Its length should be at least as long as the
#'   intervention effect duration.
#' @param ivm_coverage Numeric value giving ivermectin treatment coverage
#'   among the eligible population. Default is `0.8`.
#' @param ivm_min_age Numeric value giving the minimum eligible age for
#'   ivermectin treatment. Default is `5`.
#' @param ivm_max_age Numeric value giving the maximum eligible age for
#'   ivermectin treatment. Default is `200`.
#' @param bites_per_3_days Numeric value giving the number of bites per
#'   three-day period. Default is `1`.
#'
#' @return
#' A named list containing:
#' \describe{
#'   \item{ttt}{Integer time vector from `0` to `time_period`.}
#'   \item{eff_len}{Length of the ivermectin effect period, calculated as
#'   `length(hazard_profile)`.}
#'   \item{haz}{The supplied ivermectin hazard profile.}
#'   \item{ivm_cov_par}{Ivermectin treatment coverage parameter.}
#'   \item{ivm_min_age}{Minimum eligible treatment age.}
#'   \item{ivm_max_age}{Maximum eligible treatment age.}
#'   \item{B2}{Number of bites per three-day period.}
#'   \item{IVRM_start}{Time profile indicating the active ivermectin round
#'   start time for each time point.}
#'   \item{prop_human_HR_threshold}{Time profile for the proportion of humans
#'   above the high-risk threshold.}
#' }
#'
#' @export
ivm_fun_stag_cov <- function(IVM_start_times,
                             time_period,
                             hazard_profile,
                             prop_human_HR_threshold,
                             ivm_coverage = 0.8,
                             ivm_min_age = 5,
                             ivm_max_age = 200,
                             bites_per_3_days = 1
){

  # function to make the ivermectin time profile
  make_IVRM = function(IVRM_st, eff_len,ttt){
    IVRM =  rep(max(ttt), length(ttt))
    for(i in 1:length(IVRM_st)){
      IVRM[ttt >= IVRM_st[i]-1 & ttt < IVRM_st[i] + eff_len -1] = IVRM_st[i]
    }
    return(IVRM)
  }

  #make function for the prop_human_HR_threshold time profile
  make_stag_cov = function(IVRM_st, eff_len,ttt, prop_human_HR_threshold){
    prop <- rep(0, length(ttt))  # Start with 0s everywhere

    for (i in seq_along(IVRM_st)) {
      start_day <- IVRM_st[i]
      #idx <- which(ttt >= start_day-1 & ttt < start_day + eff_len-1) #days when intervention is on
      idx <- which(ttt >= start_day & ttt < start_day + eff_len-1) #days when intervention is on
      # Only assign if there's enough length in prop_human_HR_threshold
      if (length(idx) > 0 && length(prop_human_HR_threshold) >= length(idx)) {
        prop[idx] <- prop_human_HR_threshold[1:length(idx)]
      } else {
        warning("eff_len exceeds length of prop_human_HR_threshold or ttt range is too short")
      }
    }
    return(prop)
  }

  ttt = 0:time_period
  eff_len = length(hazard_profile)
  IVRM_start = make_IVRM(IVM_start_times, eff_len, ttt)
  #ivm_cov = ivm_coverage*(exp(-ivm_min_age/21) - exp(-ivm_max_age/21))
  prop_human_HR_threshold = make_stag_cov(IVM_start_times,eff_len, ttt, prop_human_HR_threshold)

  op = list(ttt = ttt,
            eff_len = eff_len,
            haz = hazard_profile,
            ivm_cov_par = ivm_coverage,
            ivm_min_age = ivm_min_age,
            ivm_max_age = ivm_max_age,
            B2 = bites_per_3_days,
            IVRM_start=IVRM_start,
            prop_human_HR_threshold = prop_human_HR_threshold

  )
  return(op)

}

#' Replace parameter values in a named list
#'
#' Updates a named list of parameters by replacing existing elements with
#' values supplied in another named list. Parameters in `replace_params`
#' overwrite elements of the same name in `statename`, while all other
#' elements are retained.
#'
#' @param replace_params A named list of parameter values to replace or add.
#'   Elements with names matching those in `statename` will overwrite the
#'   existing values.
#' @param statename A named list of existing parameter values.
#'
#' @return
#' A named list containing the parameters in `statename` with any matching
#' entries replaced by those in `replace_params`.
#'
#' @examples
#' defaults <- list(a = 1, b = 2, c = 3)
#' updates <- list(b = 10, d = 4)
#' new_params(updates, defaults)
#' # Returns: list(a = 1, c = 3, b = 10, d = 4)
#'
#' @export
new_params <- function(replace_params,
                       statename
){
  statename = statename[-which(names(statename) %in% names(replace_params))]
  statename_new = c(statename,replace_params)
  return(statename_new)
}

#' Calculate the effective vector mortality rate
#'
#' Computes the effective mortality rate across vector compartments by
#' combining baseline mortality and intervention-specific mortality,
#' weighted by the corresponding compartment populations and normalized
#' by the total vector population.
#'
#' @param op A list containing the model state variables and parameters.
#'   The list is expected to include mortality rates (`mu`, `mu_vi`),
#'   vector compartment populations (`Sv`, `Sv_F1`, `Sv_F2`, `Ev_F1`,
#'   `Ev_F2`, `Iv_F1`, `Iv_F2`, `Sx_F1`, `Sx_F2`, `Ex_F1`, `Ex_F2`,
#'   `Ix_F1`, `Ix_F2`), and the total vector population (`mv`).
#'
#' @return
#' A numeric vector giving the effective mortality rate at each time point
#' or spatial unit, calculated as the total mortality across all vector
#' compartments divided by the total vector population (`op$mv`).
#'
#' @export
effective_mortality <- function(op){
  x1 = op$mu*(op$Sv + op$Sv_F1 + op$Sv_F2 + rowSums(op$Ev_F1) + rowSums(op$Ev_F2) + op$Iv_F1 + op$Iv_F2)

  mmat = matrix(0, nrow = length(op$Sv), ncol = ncol(op$mu_vi))

  for(j in 1:ncol(op$mu_vi)){
    mmat[,j] = op$mu_vi[1,j] * (op$Sx_F1[,j] + op$Sx_F2[,j] + op$Ix_F1[,j] + op$Ix_F2[,j] + rowSums(op$Ex_F1[,,j]) + rowSums(op$Ex_F2[,,j]))
  }
  x2 = rowSums(mmat)
  return((x1 + x2)/op$mv)
}
