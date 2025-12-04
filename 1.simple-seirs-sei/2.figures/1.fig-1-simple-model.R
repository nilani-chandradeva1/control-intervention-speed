#figures for extreme density dependence and Logistic growth.

require(tidyverse)

df_1_emerge_T <- readRDS("1.simple-seirs-sei/output/model_results_base_df_constant_emergence_TRUE.rds")  #baseline constant emergence T
df_2_emerge_T <- readRDS("1.simple-seirs-sei/output/model_results_df_constant_emergence_TRUE.rds") #int constant emergence T

df_3_emerge_F <- readRDS("1.simple-seirs-sei/output/model_results_base_df_carrying_capacity.rds") #baseline Logistic growth
df_4_emerge_F <- readRDS("1.simple-seirs-sei/output/model_results_df_carrying_capacity.rds") #int Logistic growth

df_1_emerge_T <- df_1_emerge_T %>%
  select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D) %>%
  mutate(C0_daily = c(C0[1], diff(C0))) %>%
  ungroup()

df_2_emerge_T <- df_2_emerge_T %>%
  select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed) %>%
  mutate(C_daily = c(C[1], diff(C))) %>%
  ungroup()

df_3_emerge_F <- df_3_emerge_F %>%
  #select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D) %>%
  mutate(C0_daily = c(C0[1], diff(C0))) %>%
  ungroup()

df_4_emerge_F <- df_4_emerge_F %>%
  #select(-label) %>%
  mutate(prop_killed = delta_D/M0) %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed) %>%
  mutate(C_daily = c(C[1], diff(C))) %>%
  ungroup()

int_scenarios <- rbind(df_2_emerge_T, df_4_emerge_F)

delta_D_vec <- unique(int_scenarios$delta_D)
M0_vec <- unique(int_scenarios$M0)

int_scenarios_main <- int_scenarios %>%
  filter(delta_D == delta_D_vec[2] & M0 == M0_vec[2])



base_scenarios <- rbind(df_1_emerge_T, df_3_emerge_F) %>%
  rename(C = C0,
         C_daily = C0_daily)

base_scenarios_main <- base_scenarios %>%
  filter(M0 == M0_vec[2])

df_all_main <- rbind(int_scenarios_main, base_scenarios_main)
unique(df_all_main$delta_t)

start_int <- 100

scenario_pals <- c('#1b9e77','#d95f02','#7570b3') #colour for each product
scenario_pals2 <- c('#e7298a', scenario_pals) #baseline colour

mosq_killed_plot <- ggplot(df_all_main, aes(x = t-start_int, y = D, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(
    text = element_text(size = 14),
    legend.position = c(0.8, 0.5))+
  #guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Cumulative number of mosquitoes \n killed by intervention")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                       values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2,
                      labels = c("No intervention", "10-day killing period", "30-day killing period", "90-day killing period"),
                      name = "Scenario")  +
  coord_cartesian(xlim = c(-10, 250), ylim = c(0,2000))

daily_inc_plot <- ggplot(df_all_main, aes(x = t-start_int, y = C_daily, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(legend.position = c(0.8, 0.5),
        text = element_text(size = 14))+
  guides(col = "none", linetype = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Daily incidence")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 250))

prevalence_plot <- ggplot(df_all_main, aes(x = t-start_int, y = (I_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(
    text = element_text(size = 14))+
  ylab("Prevalence(%) in humans")+
  ylim(0, 24)+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  theme(legend.position = c(0.6, 0.5))+
  scale_colour_manual(values = scenario_pals2, labels = c("Baseline (no intervention)", "10 days", "30 days", "90 days"),
                      name = "Time taken to kill target mosquitoes")+
  coord_cartesian(xlim = c(-10, 250))+
  guides(col = "none", lty = "none")


df_all_main_baseline <- df_all_main %>%
  filter(delta_t == 0) %>%
  rename(baseline_delta_t = delta_t) %>%
  #select(I_h, N, baseline_delta_t) %>%
  mutate(scenario = "baseline",
         prev_baseline = I_h/N) %>%
  select(prev_baseline, baseline_delta_t,t, constant_emergence)

df_all_main_ints <- df_all_main %>%
  filter(delta_t != 0) %>%
  #select(I_h, N, delta_t) %>%
  mutate(scenario = "int",
         prev_int = I_h/N) %>%
  select(prev_int, delta_t,t, constant_emergence) %>%
  pivot_wider(names_from = delta_t, values_from = prev_int)

df_all_main_compare <- left_join(df_all_main_ints, df_all_main_baseline) %>%
  mutate(eff_prev_10d = ((prev_baseline - `10`)/prev_baseline)*100,
         eff_prev_30d = ((prev_baseline - `30`)/prev_baseline)*100,
         eff_prev_90d = ((prev_baseline - `90`)/prev_baseline)*100) %>%
  select(t, constant_emergence,eff_prev_10d,eff_prev_30d,eff_prev_90d) %>%
  pivot_longer(cols = c(eff_prev_10d,eff_prev_30d,eff_prev_90d),names_to = "delta_t", values_to = "eff_prev")

prev_eff_plot <- ggplot(df_all_main_compare, aes(x = t-start_int, y = eff_prev, col = as.factor(delta_t), lty = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(
    text = element_text(size = 14))+
  ylab("Percentage reduction (%) in \n prevalence")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  theme(legend.position = "none")+
  scale_colour_manual(values = scenario_pals2[2:4], labels = c("10 days", "30 days", "90 days"),
                      name = "Time taken to kill target mosquitoes")+
  coord_cartesian(xlim = c(-10, 250))

df_all_main_compare %>%
  filter(t == 110+28) %>%
  group_by(delta_t, constant_emergence) %>%
  summarise(eff_prev = eff_prev)


df_all_main_compare %>%
  group_by(delta_t, constant_emergence) %>%
  filter(eff_prev == max(eff_prev)) %>%
  select(delta_t, constant_emergence, t, eff_prev)

df_all_main_compare %>%
  group_by(delta_t, constant_emergence) %>%
  filter(eff_prev < 0) %>%
  slice_min(t, with_ties = FALSE) %>%
  select(delta_t, constant_emergence, t, eff_prev)



asprev_eff_plot2 <- ggplot(df_all_main_compare, aes(x = t-start_int, y = eff_prev, col = as.factor(delta_t), lty = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(
    text = element_text(size = 14))+
  ylab("Percentage reduction (%) in \n prevalence")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  theme(legend.position = c(0.6, 0.7))+
  scale_colour_manual(values = scenario_pals2[2:4], labels = c("Basline","10 days", "30 days", "90 days"),
                      name = "Time taken to kill target mosquitoes")+
  coord_cartesian(xlim = c(-10, 250))+
  guides(colour = "none", lty = "none")


mosq_pop_plot <- ggplot(df_all_main, aes(x = t-start_int, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(
    text = element_text(size = 14))+
  ylab("Mosquito population size")+
  guides(col = "none", linetype = "none")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2) +
  coord_cartesian(xlim = c(-10, 250), ylim = c(0,2000))


df_all_main %>%
  group_by(delta_t, constant_emergence) %>%
  summarise(min_M = min(M)) #lowest pop size reached if assume constant emergence - Logistic growth stops pop from dropping too low

Re_t_plot <- ggplot(df_all_main, aes(x = t-start_int, y = Re_t, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(legend.position = c(0.7, 0.7))+
  ylab(expression("Effective reproduction number " ~ R[e]))+
  #guides(col = "none", linetype = "none")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlim(-10,300)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2, labels = c("Baseline (no intervention)", "10 days", "30 days", "90 days"),
                      name = "Time taken to kill target mosquitoes")+
 coord_cartesian(ylim = c(0, 2.5))


#TC: move Re plot to the SM
#ggsave(Re_t_plot, file = "1.simple-seirs-sei/plots/Re_t_plot_SM.pdf")



figure_dynamics <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot,
                                      daily_inc_plot, prevalence_plot,
                                      prev_eff_plot,
                                      align = "v",
                                      labels = c("A", "B", "C", "D"))

#then epi impact plot for:
#different m0 (vector to host ratio - defines endemicity)
#different proportion of total vectors killed
#each faceted by the time of the measurement

#epi impact

#measure impact at different time periods

eqm_point <- 350

#assign time periods
time_ranges <- tibble(
  time_period = c("10d", "30d", "90d", "250d"), #number of days since start
  start = c(100, 100, 100, 100), #start times for measurements
  end = c(110, 130, 190, eqm_point) #end times for measurements
)


model_base_long <- base_scenarios_main %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

model_base_epi <- model_base_long %>%
  group_by(m0, M0, time_period, constant_emergence) %>% #for each m0, M0 and time period, sum total cases
  summarise(tot_cases_baseline = sum(C_daily), .groups = "drop")

model_int_long <- int_scenarios_main %>%
  select(t,delta_t, m0, M0, delta_D, prop_killed, C_daily, constant_emergence) %>%
  filter(prop_killed == 0.95) %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

model_int_summary <- model_int_long %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed, time_period, constant_emergence) %>%
  summarise(tot_cases_int = sum(C_daily), .groups = "drop")

#compare difference in incidence
summary_impact <- left_join(model_int_summary, model_base_epi) %>% #at each time period, for each scenario, what is rel diff in prev
  mutate(abs_diff_cases = tot_cases_baseline-tot_cases_int,
         rel_diff_cases = ((tot_cases_baseline-tot_cases_int)/tot_cases_baseline)*100,
         constant_emergence = case_when(constant_emergence == TRUE ~ "Constant emergence",
                                        constant_emergence == FALSE ~ "Logistic growth"))


impact_plot <- ggplot(summary_impact, aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")), y = log(rel_diff_cases+1), fill = as.factor(delta_t)))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(vars(constant_emergence))+
  theme_bw(base_size = 14)+
  #ylim(0,4)+
  scale_y_continuous(limits = c(0,4), labels = c(0, 10, 20, 30, 40))+
  ylab("Percentage (%) cases averted \n due to intervention")+
  xlab("Time period over which incidence measured since intervention start")+
  labs(fill = "Time to complete MDA (days)")+
  theme(legend.position = c(0.8, 0.9))+
  scale_fill_manual(values = scenario_pals)+
  guides(fill = "none")

##mosquito pop sizes at key times
model_base_long2 <- base_scenarios_main %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

d10_since_start <- int_scenarios_main %>%
  filter(t == 110) %>%
  filter(delta_t == 10) %>%
  group_by(delta_t, constant_emergence) %>%
  select(M) %>%
  mutate(M0 = 2000,
         prop_M_killed = ((M0-M)/M0)*100,
         t = "10d since start")

d30_since_start <- int_scenarios_main %>%
  filter(t == 140) %>%
  filter(delta_t == 30) %>%
  group_by(delta_t, constant_emergence) %>%
  select(M) %>%
  mutate(M0 = 2000,
         prop_M_killed = ((M0-M)/M0)*100,
         t = "30d since start")

d90_since_start <- int_scenarios_main %>%
  filter(t == 200) %>%
  filter(delta_t == 90) %>%
  group_by(delta_t, constant_emergence) %>%
  select(M) %>%
  mutate(M0 = 2000,
         prop_M_killed = ((M0-M)/M0)*100,
         t = "90d since start")

mosq_pop_tab <- do.call("rbind", list(d10_since_start, d30_since_start, d90_since_start))


model_int_long2 <- int_scenarios_main %>%
  select(t,delta_t, m0, M0, M, delta_D, prop_killed, constant_emergence) %>%
  filter(prop_killed == 0.95) %>%
  crossing(time_ranges) %>%
  filter(t >= start & t <= end)

model_int_summary <- model_int_long %>%
  group_by(delta_t, m0, M0, delta_D, prop_killed, time_period, constant_emergence) %>%
  summarise(tot_cases_int = sum(C_daily), .groups = "drop")

#compare difference in incidence
summary_impact <- left_join(model_int_summary, model_base_epi) %>% #at each time period, for each scenario, what is rel diff in prev
  mutate(abs_diff_cases = tot_cases_baseline-tot_cases_int,
         rel_diff_cases = ((tot_cases_baseline-tot_cases_int)/tot_cases_baseline)*100,
         constant_emergence = case_when(constant_emergence == TRUE ~ "Constant emergence",
                                        constant_emergence == FALSE ~ "Logistic growth"))

library(ggpattern)
impact_plot_main <- ggplot(summary_impact,
       aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")),
           #y = rel_diff_cases,
           y = sqrt(rel_diff_cases),
           fill = as.factor(delta_t),
           pattern = as.factor(constant_emergence))) +
  geom_bar_pattern(
    stat = "identity",
    position = position_dodge(),
    colour = "black",                # Border of bars
    pattern_colour = "black",        # Pattern line color
    pattern_fill = NA,               # Transparent so bar fill shows
    pattern_density = 0.4,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.5
  ) +
  theme_minimal() +
  scale_y_continuous(
  #  breaks = log(c(0, 10, 20, 30, 40) + 1),
    labels = c(0, 4, 16, 36)
  )+
  ylab("Percentage (%) cases averted \n due to intervention")+
  xlab("Time period (days) over \n which incidence measured since intervention started") +
  labs(fill = "Time to complete MDA (days)",
       pattern = "Adult emergence") +
  theme_bw(base_size = 14)+
  theme(legend.position = c(0.7, 0.8),
        text = element_text(size = 14)) +
  scale_fill_manual(values = scenario_pals) +
  scale_pattern_manual(values = c("Logistic growth" = "none",
                                  "Constant emergence" = "stripe"),
                       breaks = c("Logistic growth", "Constant emergence"))+
  guides(pattern_spacing = 0.5,
         pattern = guide_legend(
           override.aes = list(fill = "white"), # Force fill color in legend
         ),
         #fill = guide_legend(override.aes = list(pattern = "none")),
         fill = "none")






#impact_plot_main <- ggplot(summary_impact,
#                           aes(x = factor(time_period, levels = c("10d", "30d", "90d", "250d")),
#                               y = log(rel_diff_cases + 1),
#                               fill = as.factor(delta_t),
#                               pattern = as.factor(constant_emergence))) +
#  geom_bar_pattern(
#    stat = "identity",
#    position = position_dodge(),
#    colour = "black",                # Border of bars
#    pattern_colour = "black",        # Pattern line color
#    pattern_fill = NA,               # Transparent so bar fill shows
#    pattern_density = 0.4,
#    pattern_spacing = 0.05,
#    pattern_key_scale_factor = 0.5
#  ) +
#  theme_minimal() +
#  scale_y_continuous(limits = c(0, 4), labels = c(0, 10, 20, 30, 40)) +
#  ylab("Percentage (%) cases averted \n due to intervention")+
#  xlab("Time period (days) over \n which incidence measured since intervention started") +
#  labs(fill = "Time to complete MDA (days)",
#       pattern = "Adult emergence") +
#  theme_bw(base_size = 14)+
#  theme(legend.position = c(0.7, 0.8),
#        text = element_text(size = 14)) +
#  scale_fill_manual(values = scenario_pals) +
#  scale_pattern_manual(values = c("Logistic growth" = "none",
#                                  "Constant emergence" = "stripe"),
#                       breaks = c("Logistic growth", "Constant emergence"))+
#  guides(pattern_spacing = 0.5,
#         pattern = guide_legend(
#    override.aes = list(fill = "white"), # Force fill color in legend
#  ),
#  #fill = guide_legend(override.aes = list(pattern = "none")),
#  fill = "none")

summary_impact %>%
  filter(time_period == "250d")
0.7-0.310 #delta_t 10 Logistic growth - delta_t 90 Logistic growth
0.391-0.290 #delta_t 10 Logistic growth - delta_t 90 constant emergence




#figure_dynamic_impact <- cowplot::plot_grid(figure_dynamics, impact_plot_main,
#                                            labels = c("", "E"))
#
figure_dynamic_impact <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot,
                                           daily_inc_plot, prevalence_plot,
                                           prev_eff_plot, impact_plot_main,
                                           labels = c("A", "B", "C", "D", "E", "F"),
                                           nrow = 2, ncol = 3, align = "v")



ggsave(figure_dynamic_impact, file = "1.simple-seirs-sei/plots/fig_1_simple_model_plot.pdf")



figure_dynamic_impact_pres <- cowplot::plot_grid(mosq_killed_plot, mosq_pop_plot,
                                            prev_eff_plot2,
                                            impact_plot_main,
                                            labels = c("A", "B", "C", "D"),
                                            nrow = 2, ncol = 2, align = "v")
ggsave(figure_dynamic_impact_pres, file = "1.simple-seirs-sei/plots/fig_1_simple_model_plot_pres.pdf")
ggsave(figure_dynamic_impact_pres, file = "1.simple-seirs-sei/plots/fig_1_simple_model_plot_pres.png")
