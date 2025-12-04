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

ggplot(df_all_main, aes(x = t-start_int, y = M, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  #guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("mosquito density")+
  #scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
  #                     values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)#+
  #coord_cartesian(xlim = c(-10, 250), ylim = c(0,25))

inf_mosq <- ggplot(df_all_main, aes(x = t-start_int, y = (I_v/M)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Infectious vectors (%)")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 250), ylim = c(0,25))

susceptible_people <- ggplot(df_all_main, aes(x = t-start_int, y = (S_h/N)*100, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw()+
  theme(
    text = element_text(size = 14),
    legend.position = c(0.9, 0.6))+
  guides(col = "none", lty = "none")+
  labs(col = "Duration of killing (days)")+
  ylab("Susceptible humans (%)")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2)+
  coord_cartesian(xlim = c(-10, 250), ylim = c(0,10))

Re_t_plot <- ggplot(df_all_main, aes(x = t-start_int, y = Re_t, col = as.factor(delta_t), linetype = as.factor(constant_emergence)))+
  geom_line(linewidth = 0.9)+
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 1.1)+
  theme_bw(base_size = 14)+
  theme(legend.position = "right",
        legend.direction = "vertical")+
  ylab(expression("Effective reproduction number " ~ R[e]))+
  #guides(col = "none", linetype = "none")+
  scale_linetype_manual(name = "Adult emergence", labels = c("Logistic growth", "Constant emergence"),
                        values = c("solid", "dotdash"))+
  xlim(-10,250)+
  xlab("Time since intervention started (days)")+
  scale_colour_manual(values = scenario_pals2, labels = c("No intervention", "10-day killing period", "30-day killing period", "90-day killing period"),
                      name = "Scenario")+
  coord_cartesian(ylim = c(0, 2.5), xlim = c(-10, 250))

legend <- cowplot::get_legend(Re_t_plot + theme(legend.position = "right"))

Re_t_plot_no_legend <- Re_t_plot + theme(legend.position = "none")

sens_dynamics_epi <- cowplot::plot_grid(inf_mosq, susceptible_people, Re_t_plot_no_legend,
                                        legend,
                                        labels = c("A", "B", "C", ""),
                                        nrow = 2,
                                        align = "v")

ggsave(sens_dynamics_epi, file = "1.simple-seirs-sei/plots/sm_fig_1_simple_model_plot_epi.pdf")

