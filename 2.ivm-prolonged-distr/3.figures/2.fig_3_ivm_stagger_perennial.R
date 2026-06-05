#plots for impact of prolonged vs immediate ivermectin distributions in perennial settings
#high endemicity

require(tidyverse)
start <- (365*5)+200
df_distr_all <- readRDS("2.ivm-prolonged-distr/output/df_distr_HR_3m_perennial.rds")

model_types <-  unique(df_distr_all$model_type)

df_distr_all %>%
  filter(t == 1) %>% #EIR at t = 1 0.00407 or 0.203
  group_by(ref) %>%
  summarise(eir = EIRout) %>%
  distinct()  #so ref 1-3 is low EIR and 4-6 is high

#match to the initial EIR passed into the model
df_distr_all <- df_distr_all %>%
  mutate(init_EIR = case_when(ref <= 3 ~ 2,
                              TRUE ~ 100)) %>%
  filter(model_type != "all-in-one")

unique(df_distr_all$init_EIR)
unique(df_distr_all$model_type)

#distr_pals <- c('#66c2a5','#fc8d62','#8da0cb',)
distr_pals <- c('#a6d854', '#8da0cb', '#66c2a5', '#fc8d62')
covs <- unique(df_distr_all$ivm_cov_par)

df_distr <- df_distr_all %>%
  filter(init_EIR == 100) #filter for high endemicity

####
#measure efficacy at different time points, using MATAMAL/BOHEMIA trial protocols as examples

#MATAMA trial: prevalence survey takes place 4 weeks after last MDA
#BOHEMIA: incidence measured from first MDA, for 6 months
matamal_survey <- (30*4) + start
plot_matamal <- matamal_survey-start #diff between start and survey
2145-start
bohemia_inc_period <- start+(6*30)
plot_bohemia <- bohemia_inc_period- start

second_mda <- 30+start-start
third_mda <- 60+start-start

df_distr <- df_distr %>%
  mutate(ivm_cov_par = case_when(model_type == "baseline" ~ 0.7, ##just for plotting purposes
                                 TRUE ~ ivm_cov_par)) %>%
  filter(ivm_cov_par == 0.7)

#dotted arrow shows point of matamal prevalence survey
#bound box is bohemia incidence period
#full shaded area is area that I measure
mv_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = mv, col = as.factor(model_type))) +
  geom_line(size = 1.1) +
  theme_bw(base_size = 14) +
  ylim(0, 50) +
  ylab("Mosquito density")+
  theme(legend.position = "none") +

  scale_color_manual(values = distr_pals,
                     breaks = c("baseline", "all-in-one-stag", "10d-stagger", "20d-stagger"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  coord_cartesian(xlim = c(-0.25, 1))+
  xlab("Years since intervention started") +
  geom_segment(x = 0, y = 48, xend = 0, yend = 42, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = 48, xend = second_mda/365, yend = 42, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = 48, xend = third_mda/365, yend = 42, arrow = arrow(),
               col = "black", size = 1.1) #+
  #geom_segment(x = plot_matamal/365, y = 48, xend = plot_matamal/365, yend = 42,
  #             arrow = arrow(),
  #             col = "blue", size = 1.1) +
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 20, ymax = 42,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = 42,alpha = 0.1)

eir_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = EIRout, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  coord_cartesian(ylim = c(0, 0.3), xlim = c(-0.25, 1))+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline", "all-in-one-stag", "10d-stagger", "20d-stagger"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  guides(col = "none")+
  xlim(-0.25, 1) +
  ylab("Average number of infectious bites \n per person per day (daily EIR)")+
  xlab("Years since intervention started") +
  geom_segment(x = 0, y = 0.28, xend = 0, yend =0.23, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = 0.28, xend = second_mda/365, yend =0.23, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = 0.28, xend = third_mda/365, yend =0.23, arrow = arrow(),
               col = "black", size = 1.1)#+
  #geom_segment(x = plot_matamal/365, y = 0.28, xend = plot_matamal/365, yend = 0.23,
  #             arrow = arrow(),
  #             col = "blue", size = 1.1)+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = 0.23,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = 0.23,alpha = 0.1)

prev_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = slide_prev0to80*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  coord_cartesian(ylim = c(0, 75), , xlim = c(-0.25, 1))+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline", "all-in-one-stag", "10d-stagger", "20d-stagger"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  #guides(col = "none")+
  theme(legend.position = c(0.4, 0.2))+
  xlab("Years since intervention started")+
  ylab("All-age slide prevalence (%)") +
  geom_segment(x = 0, y = 57, xend = 0, yend =50, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = 57, xend = second_mda/365, yend =50, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = 57, xend = third_mda/365, yend =50, arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = 57, xend = plot_matamal/365, yend = 50,
               arrow = arrow(),
               col = "blue", size = 1.1)#+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = 64,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = 64,alpha = 0.1)


inc_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline", "all-in-one-stag", "10d-stagger", "20d-stagger"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  guides(col = "none")+
  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  coord_cartesian(xlim = c(-0.25, 1), ylim = c(0, 7.5))+
  xlab("Years since intervention started")+
  geom_segment(aes(x = 0, y = 7, xend = 0, yend =6), arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(aes(x = second_mda/365, y = 7, xend = second_mda/365, yend =6), arrow = arrow(),
               col = "black", size = 1.1)+
  geom_segment(aes(x = third_mda/365, y = 7, xend = third_mda/365, yend =6), arrow = arrow(),
               col = "black", size = 1.1)+
  #geom_segment(aes(x = plot_matamal/365, y = 7, xend = plot_matamal/365, yend = 6),
  #             arrow = arrow(),
  #             col = "blue", size = 1.1) +
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = 6,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = 6,alpha = 0.1)

dynamics <- cowplot::plot_grid(mv_plot, eir_plot, prev_plot, inc_plot, align = "v",
                               labels = c("A", "B", "C", "D"))

#then look at epi impact, at different times to understand the influence of the timing of the measurement

#incidence measurement
df_distr_wide_setting_ivm <- df_distr_all %>%
  select(t,ref, init_EIR, model_type, ivm_cov_par, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  }) %>%
  select(-baseline)

df_distr_wide_setting_baseline <- df_distr_all %>%
  filter(model_type == "baseline") %>%
  select(t,  ref,init_EIR, model_type, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  })

df_distr_wide_setting <- left_join(df_distr_wide_setting_ivm,df_distr_wide_setting_baseline, by = c("ref",
                                                                                                    "t",
                                                                                                    "init_EIR"))
#prevalence measurement
df_distr_wide_setting_prev_ivm<- df_distr_all %>%
  select(t, ref, init_EIR, model_type, ivm_cov_par, slide_prev0to80) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to80)
  }) %>%
  select(-baseline)

df_distr_wide_setting_prev_baseline <- df_distr_all %>%
  filter(model_type == "baseline") %>%
  select(t,ref,init_EIR, model_type, slide_prev0to80) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to80)
  })

df_distr_wide_setting_prev <- left_join(df_distr_wide_setting_prev_ivm,
                                        df_distr_wide_setting_prev_baseline,
                                        by = c("ref","t", "init_EIR"))

#overall impact measurement: cases averted in u5 from start to 1y later
impact_overall <- df_distr_wide_setting %>%
  group_by(ref, init_EIR, ivm_cov_par) %>%
  #group_by(ref) %>%
  filter(between(t, start, start + 365)) %>% #sum cases one year from start
  summarise(tot_cases_baseline = sum(baseline)*1000,
            tot_cases_10d = sum(`10d-stagger`)*1000,
            tot_cases_20d = sum(`20d-stagger`)*1000,
            tot_cases_all_in_stag = sum(`all-in-one-stag`)*1000) %>%
  mutate(
    impact_all_in_stag = ((tot_cases_baseline - tot_cases_all_in_stag)/tot_cases_baseline)*100,
    impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
    impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
    scenario = "Once year since start") #measuring overall impact


#overall impact acc. BOHEMIA protocol
impact_bohemia <- df_distr_wide_setting %>%
  group_by(ref, init_EIR, ivm_cov_par) %>%
  filter(between(t, start, bohemia_inc_period)) %>%
  summarise(tot_cases_baseline = sum(baseline)*1000,
            tot_cases_10d = sum(`10d-stagger`)*1000,
            tot_cases_20d = sum(`20d-stagger`)*1000,
            tot_cases_all_in_stag = sum(`all-in-one-stag`)*1000) %>%
  mutate(
    impact_all_in_stag = ((tot_cases_baseline - tot_cases_all_in_stag)/tot_cases_baseline)*100,
    impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
    impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
    scenario = "bohemia") #measuring impact across 6 months

#overall impact for MATAMAL protocol
impact_matamal <- df_distr_wide_setting_prev %>%
  group_by(ref, ivm_cov_par, init_EIR) %>%
  filter(t == matamal_survey) %>%
  summarise(prev_baseline = baseline,
            prev_10d =`10d-stagger`,
            prev_20d = `20d-stagger`,

            prev_all_in_stag = `all-in-one-stag`) %>%
  mutate(
    impact_all_in_stag = ((prev_baseline - prev_all_in_stag)/prev_baseline)*100,
    impact_20d = ((prev_baseline - prev_20d)/prev_baseline)*100,
    impact_10d = ((prev_baseline - prev_10d)/prev_baseline)*100,
    scenario = "matamal") #measuring point prevalence at 4 weeks after last MDA


impact_measurements <- do.call("rbind", list(impact_overall, impact_bohemia, impact_matamal)) %>%
  select(init_EIR, ivm_cov_par,impact_all_in_stag, impact_20d, impact_10d, scenario)

impact_measurements_long <- impact_measurements %>%
  pivot_longer(cols = c( impact_all_in_stag, impact_20d, impact_10d), names_to = "intervention", values_to = "impact")

levels_x <- unique(impact_measurements_long$scenario)

#distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
distr_pals <- c('#a6d854', '#8da0cb', '#66c2a5', '#fc8d62')
distr_pals2 <- distr_pals[2:4]


int_order <- c("impact_all_in_stag", "impact_10d", "impact_20d")
impact_measurements_long$intervention <- factor(impact_measurements_long$intervention,
                                                levels = int_order)

cov_error <- impact_measurements_long %>%
  ungroup() %>%
  select(-ref) %>%
  filter(init_EIR == 100 & ivm_cov_par %in% c(0.5, 0.9)) %>%
  pivot_wider(names_from = ivm_cov_par, values_from = impact) %>%
  rename(cov_low_0.5 = `0.5`,
         cov_high_0.9 = `0.9`)

impact_main_plot <- ggplot() +
  # Bars
  geom_bar(
    data = impact_measurements_long %>% filter(init_EIR == 100 & ivm_cov_par == 0.7),
    aes(x = factor(scenario), y = impact, fill = as.factor(intervention)),
    stat = "identity",
    position = position_dodge(width = 0.9)
  ) +
  # Error bars
  geom_errorbar(
    data = cov_error, size = 1.1,
    aes(
      x = factor(scenario),
      ymin = cov_low_0.5,
      ymax = cov_high_0.9,
      fill = as.factor(intervention)   # match fill to align dodging
    ),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.7, 0.8)) +
  guides(fill = "none")+
  scale_fill_manual(
    values = distr_pals2, name = "Scenario",
    breaks = c("impact_all_in_stag", "impact_10d", "impact_20d"),
    labels = c("Synchronised MDA", "10-day MDA", "20-day MDA")
    )+
  xlab("Time of measurement") +
  ylab("Efficacy (%)") +
  scale_x_discrete(labels = c(
    "bohemia" = "Incidence U5s \n (start to 6m later)",
    "matamal" = "All-age prevalence \n (1m after last MDA)",
    "Once year since start" = "Incidence U5s \n (start to 1y later)"
  ))

impact_measurements_long %>%
  filter(init_EIR == 100 & ivm_cov_par == 0.7 & scenario == "Once year since start")%>%
  group_by(scenario, intervention)

# all-in 15.7
#20d 19.1
cov_error %>%
  filter(init_EIR == 100 & scenario == "Once year since start" )  # all-in 13.9-16.8
#20d 16.8 to 20.6


dynamics_perennial <- cowplot::plot_grid(mv_plot, eir_plot, prev_plot, inc_plot,
                                         labels = c("A", "B", "C", "D"),
                                         align = "v")

plot_perennial <- cowplot::plot_grid(dynamics_perennial, impact_main_plot,
                                     labels = c("", "E"))

ggsave(plot_perennial, file = "2.ivm-prolonged-distr/plots/fig_3_plot_perennial.pdf", dpi = 100)
ggsave(plot_perennial, file = "2.ivm-prolonged-distr/plots/fig_3_plot_perennial.png")
