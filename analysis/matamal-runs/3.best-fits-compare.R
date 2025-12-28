malariasim_endec_killing_inputs <- readRDS("analysis/matamal-runs/malariasim_endec_killing_inputs.rds")

malariasim_endec_killing_inputs <- malariasim_endec_killing_inputs %>%
  mutate(Q0 = case_when(grepl("low_Q0", scenario) ~ "low",
                        grepl("medium_Q0", scenario) ~ "medium",
                        grepl("high_Q0", scenario) ~ "high"),
         ivm_cov = case_when(grepl("low_cov", scenario) ~ "low",
                             grepl("medium_cov", scenario) ~ "medium",
                             grepl("high_cov", scenario) ~ "high"))

ggplot(malariasim_endec_killing_inputs, aes(x = as.factor(Q0), y = as.factor(ivm_cov),
                                            fill = best_endec_mu))+
  geom_tile()+
  geom_text(
    aes(label = paste0("wane_endec = ", round(best_wane, 3))),
    size = 5, col = "black"
  )+
  scale_fill_distiller(palette = "RdPu", name = "mu_endec")

#odin model
my_sim_ivm_mda_out <- readRDS("analysis/matamal-runs/matamal-Q0-ivm-cov-combos.rds")
my_sim_ivm_mda_out_malsim <- readRDS("analysis/matamal-runs/matamal-Q0-ivm-cov-combos-malariasim-odin.rds")

odin_model <- my_sim_ivm_mda_out
exp_decay_model <- my_sim_ivm_mda_out_malsim

odin_model <- odin_model %>%
  mutate(ivm_cov = case_when(
    ref %in% c(1,2,3) ~ "low_cov",
    ref %in% c(4,5,6) ~ "medium_cov",
    ref %in% c(7,8,9) ~ "high_cov"
  ))

#low Q0
odin_low_Q0_low_cov <- odin_model %>%
  filter(Q0 == 0.20 & ivm_cov == "low_cov") %>%
  mutate(Q0_cov = "low_low")

odin_low_Q0_medium_cov <- odin_model %>%
  filter(Q0 == 0.20 & ivm_cov == "medium_cov") %>%
  mutate(Q0_cov = "low_medium")

odin_low_Q0_high_cov <- odin_model %>%
  filter(Q0 == 0.20 & ivm_cov == "high_cov") %>%
  mutate(Q0_cov = "low_high")

#medium Q0
odin_medium_Q0_low_cov <- odin_model %>%
  filter(Q0 == 0.82 & ivm_cov == "low_cov") %>%
  mutate(Q0_cov = "medium_low")

odin_medium_Q0_medium_cov <- odin_model %>%
  filter(Q0 == 0.82 & ivm_cov == "medium_cov") %>%
  mutate(Q0_cov = "medium_medium")

odin_medium_Q0_high_cov <- odin_model %>%
  filter(Q0 == 0.82 & ivm_cov == "high_cov") %>%
  mutate(Q0_cov = "medium_high")

#high Q0
odin_high_Q0_low_cov <- odin_model %>%
  filter(Q0 == 0.93 & ivm_cov == "low_cov") %>%
  mutate(Q0_cov = "high_low")

odin_high_Q0_medium_cov <- odin_model %>%
  filter(Q0 == 0.93 & ivm_cov == "medium_cov") %>%
  mutate(Q0_cov = "high_medium")

odin_high_Q0_high_cov <- odin_model %>%
  filter(Q0 == 0.93 & ivm_cov == "high_cov") %>%
  mutate(Q0_cov = "high_high")

odin_dfs <- do.call("rbind", list(odin_low_Q0_low_cov, odin_low_Q0_medium_cov, odin_low_Q0_high_cov,
                                  odin_medium_Q0_low_cov, odin_medium_Q0_medium_cov, odin_medium_Q0_high_cov,
                                  odin_high_Q0_low_cov,odin_high_Q0_medium_cov,odin_high_Q0_high_cov))

#exp decay model

#low Q0
exp_decay_low_Q0_low_cov <- exp_decay_model %>%
  filter(Q0 == 0.20 & endec_mu == 0.05 & wane == 0.06) %>%
  mutate(Q0_cov = "low_low")

exp_decay_low_Q0_medium_cov <- exp_decay_model %>%
  filter(Q0 == 0.20  & endec_mu == 0.08 & wane == 0.01) %>%
  mutate(Q0_cov = "low_medium")

exp_decay_low_Q0_high_cov <- exp_decay_model %>%
  filter(Q0 == 0.20  & endec_mu == 0.09 & wane == 0.00) %>%
  mutate(Q0_cov = "low_high")

#medium Q0
exp_decay_medium_Q0_low_cov <- exp_decay_model %>%
  filter(Q0 == 0.82  & endec_mu == 0.05 & wane == 0.02) %>%
  mutate(Q0_cov = "medium_low")

exp_decay_medium_Q0_medium_cov <- exp_decay_model %>%
  filter(Q0 == 0.82  & endec_mu == 0.10 & wane == 0.00) %>%
  mutate(Q0_cov = "medium_medium")

endec_mus <- unique(exp_decay_model$endec_mu)

exp_decay_medium_Q0_high_cov <- exp_decay_model %>%
  filter(Q0 == 0.82 & wane == 0.00 & endec_mu == endec_mus[8]) %>%
  mutate(Q0_cov = "medium_high")

#high Q0
exp_decay_high_Q0_low_cov <- exp_decay_model %>%
  filter(Q0 == 0.93  & endec_mu == 0.15 & wane == 0.00) %>%
  mutate(Q0_cov = "high_low")

exp_decay_high_Q0_medium_cov <- exp_decay_model %>%
  filter(Q0 == 0.93  & endec_mu == 0.11 & wane == 0.00) %>%
  mutate(Q0_cov = "high_medium")

exp_decay_high_Q0_high_cov <- exp_decay_model %>%
  filter(Q0 == 0.93  & endec_mu == 0.13 & wane == 0.00) %>%
  mutate(Q0_cov = "high_high")

exp_decay_df <- do.call("rbind", list(exp_decay_low_Q0_low_cov, exp_decay_low_Q0_medium_cov, exp_decay_low_Q0_high_cov,
                                      exp_decay_medium_Q0_low_cov, exp_decay_medium_Q0_medium_cov, exp_decay_medium_Q0_high_cov,
                                      exp_decay_high_Q0_low_cov, exp_decay_high_Q0_medium_cov, exp_decay_high_Q0_high_cov))


odin_dfs <- odin_dfs %>%
  select(t, EIRout, slide_prev0to80, clin_inc0to5, mv, Ivtot, Q0_cov, scenario)

exp_decay_df <- exp_decay_df %>%
  select(t, EIRout, slide_prev0to80, clin_inc0to5, mv, Ivtot, Q0_cov, scenario)

combo_df <- rbind(odin_dfs, exp_decay_df)

combo_df$Q0_cov <- factor(combo_df$Q0_cov, levels = c("high_high", "high_medium", "high_low",
                                                      "medium_high", "medium_medium", "medium_low",
                                                      "low_high", "low_medium", "low_low"))


combo_df <- combo_df %>%
  mutate(parameters = case_when(Q0_cov == "high_high" ~ "High Q0; High ivermectin coverage",
                                Q0_cov == "high_low" ~ "High Q0; Low ivermectin coverage",
                                Q0_cov == "high_medium" ~ "High Q0; Medium ivermectin coverage",
                                Q0_cov == "low_high" ~ "Low Q0; High ivermectin coverage",
                                Q0_cov == "low_low" ~ "Low Q0; Low ivermectin coverage",
                                Q0_cov == "low_medium" ~ "Low Q0; Medium ivermectin coverage",
                                Q0_cov == "medium_high" ~ "Medium Q0; High ivermectin coverage",
                                Q0_cov == "medium_low" ~ "Medium Q0; Low ivermectin coverage",
                                Q0_cov == "medium_medium" ~ "Medium Q0; Medium ivermectin coverage"))

pals_models <- c('#1b9e77', 'darkviolet')

itn_on <- 1275

mda_int <- 30
eff_len <- 23
july <- 30*7
y_2021 <- 365*10
y_2022 <- 365*11
endec_y1_start <- july + y_2021 #start time of first MDA
endec_y2_start <- july+ y_2022

endec_start <- c(endec_y1_start, endec_y1_start+mda_int, endec_y1_start+mda_int+mda_int,
                 endec_y2_start, endec_y2_start+mda_int, endec_y2_start+mda_int+mda_int)


iv_plot <- ggplot(combo_df, aes(x = ((t-itn_on)/365), y = Ivtot))+
  geom_line(size = 0.9, aes(col = as.factor(scenario)))+
  facet_wrap(vars(parameters), scales = "free")+
  theme_bw()+
  scale_colour_manual(values = pals_models, labels = c("Model A", "Model B3"),
                     name = "Scenario")+
  ylab("Infectious vectors")+
  theme(legend.position = c(0.9, 0.9))+
  xlab("Time (years) since first ITN campaign")+
  coord_cartesian(xlim = c(-0.1, 9))+
  geom_vline(xintercept  = (endec_start[1]-itn_on)/365, col = "pink", lty = "dashed")+
  geom_vline(xintercept  = (endec_start[6]+23-itn_on)/365, col = "pink", lty = "dashed")+



eir_plot <- ggplot(combo_df, aes(x = ((t-itn_on)/365), y = EIRout*1000, col = as.factor(scenario)))+
  geom_line(size = 0.9)+
  facet_wrap(vars(parameters), scales = "free")+
  theme_bw()+
  scale_colour_manual(values = pals_models, labels = c("Model A", "Model B3"),
                      name = "Scenario")+
  ylab("Daily EIR (per 1000 persons)")+
  theme(legend.position = c(0.9, 0.9))+
  xlab("Time (years) since first ITN campaign")+
  coord_cartesian(xlim = c(-0.1, 9))

prev_plot <- ggplot(combo_df, aes(x = ((t-itn_on)/365), y = slide_prev0to80*100, col = as.factor(scenario)))+
  geom_line(size = 0.9)+
  facet_wrap(vars(parameters), scales = "free")+
  theme_bw()+
  scale_colour_manual(values = pals_models, labels = c("Model A", "Model B3"),
                      name = "Scenario")+
  ylab("All-age slide prevalence (%)")+
  theme(legend.position = c(0.9, 0.92))+
  xlab("Time (years) since first ITN campaign")+
  coord_cartesian(xlim = c(-0.1, 9))

inc_plot <- ggplot(combo_df, aes(x = ((t-itn_on)/365), y = clin_inc0to5*1000, col = as.factor(scenario)))+
  geom_line(size = 0.9)+
  facet_wrap(vars(parameters), scales = "free")+
  theme_bw()+
  scale_colour_manual(values = pals_models, labels = c("Model A", "Model B3"),
                      name = "Scenario")+
  ylab("Clinical incidence in under 5-year-olds \n (per 1000 persons)")+
  theme(legend.position = c(0.8, 0.55))+
  xlab("Time (years) since first ITN campaign")+
  coord_cartesian(xlim = c(-0.1, 9))

ggsave(iv_plot, file = "analysis/matamal-runs/plots/iv_plot.pdf")
ggsave(eir_plot, file = "analysis/matamal-runs/plots/eir_plot.pdf")
ggsave(prev_plot, file = "analysis/matamal-runs/plots/prev_plot.pdf")
ggsave(inc_plot, file = "analysis/matamal-runs/plots/inc_plot.pdf")

