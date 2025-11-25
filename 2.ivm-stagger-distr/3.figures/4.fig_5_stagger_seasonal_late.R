#viz and epi impact of staggered vs all in one distribution
require(tidyverse)

start <- (365*5)+200
start <- start + 60 #update start time: 60 days later

df_distr_all <- readRDS("2.ivm-stagger-distr/output/df_distr_HR_3m_seasonal_late.rds")

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
  filter(model_type != "all-in-one-Sen")


distr_pals <- c('#a6d854', '#8da0cb', '#66c2a5', '#fc8d62')
covs <- unique(df_distr_all$ivm_cov_par)

df_distr <- df_distr_all %>%
  filter(init_EIR == 100) #filter for high endemicity

#Take measurements using different epi endpoints (prevalence or incidence u5) and at different time points
#use examples from recent cRCTs

#MATAMAL: 4 weeks after last MDA (prevalence)
#BOHEMIA: incidence from first MDA, for 6 months
matamal_survey <- (30*4) + start
plot_matamal <- matamal_survey-start #diff between start and survey
bohemia_inc_period <- start+(6*30)
plot_bohemia <- bohemia_inc_period- start
second_mda <- 30+start-start
third_mda <- 60+start-start

df_distr <- df_distr %>%
  mutate(ivm_cov_par = case_when(model_type == "baseline-Sen" ~ 0.7, ##just for plotting purposes
                                 TRUE ~ ivm_cov_par)) %>%
  filter(ivm_cov_par == 0.7)

#dotted arrow shows point of matamal prevalence survey
#bound box is bohemia incidence period
#full shaded area is area that I measure
max_mv <- max(df_distr$mv, na.rm = TRUE)

saveRDS(df_distr, file = "2.ivm-stagger-distr/output/df_distr_season_late.rds")


mv_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = mv, col = as.factor(model_type))) +
  geom_line(size = 1.1) +
  theme_bw(base_size = 14) +
  #ylim(0, 50) +
  ylab("Mosquito density")+
  #theme(legend.position = c(0.7, 0.3)) +
  guides(col = "none")+
  #coord_cartesian(ylim = c(0,200))+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 200))+
  xlab("Years since intervention started") +
  geom_segment(x = 0, y = max_mv+15, xend = 0, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #first MDA
  geom_segment(x = second_mda/365, y = max_mv+15, xend = second_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #second MDA
  geom_segment(x = third_mda/365, y = max_mv+15, xend = third_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #third MDA
  geom_segment(x = plot_matamal/365, y = max_mv+15, xend = plot_matamal/365, yend = max_mv,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_mv,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_mv,
           alpha = 0.05)

mv_plot_pres <- ggplot(df_distr, aes(x = (t - start)/365, y = mv, col = as.factor(model_type))) +
  geom_line(size = 1.1) +
  theme_bw(base_size = 18) +
  #ylim(0, 50) +
  ylab("Mosquito density")+
  #theme(legend.position = c(0.7, 0.3)) +
  guides(col = "none")+
  #coord_cartesian(ylim = c(0,200))+
  scale_color_manual(name = "Scenario", values = distr_pals, labels = c("10 days to complete MDA",
                                                                        "20 days to complete MDA",
                                                                        "1 day to complete MDA",
                                                                        "Baseline (no intervention)")) +
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 200))+
  xlab("Years since intervention started") +
  geom_segment(x = 0, y = max_mv+15, xend = 0, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #first MDA
  geom_segment(x = second_mda/365, y = max_mv+15, xend = second_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #second MDA
  geom_segment(x = third_mda/365, y = max_mv+15, xend = third_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #third MDA
  geom_segment(x = plot_matamal/365, y = max_mv+15, xend = plot_matamal/365, yend = max_mv,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_mv,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_mv,
           alpha = 0.05)




max_eir <- max(df_distr$EIRout, na.rm = TRUE)

eir_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = EIRout, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  guides(col = "none")+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 1))+
  ylab("Average number of infectious bites \n per person per day (daily EIR)")+
  xlab("Years since intervention started") +
  geom_segment(x = 0, y = max_eir+0.05, xend = 0, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_eir+0.05, xend = second_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_eir+0.05, xend = third_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_eir+0.05, xend = plot_matamal/365, yend = max_eir,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_eir,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_eir,
           alpha = 0.05)

max_prev <- max(df_distr$slide_prev0to80*100, na.rm = TRUE)

prev_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = slide_prev0to80*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  #ylim(0, 75)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  theme(legend.position = c(0.4, 0.2))+
  coord_cartesian(xlim = c(-0.25,1))+
  xlab("Years since intervention started")+
  ylab("All-age slide prevalence (%)") +
  geom_segment(x = 0, y = max_prev+5, xend = 0, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_prev+5, xend = second_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_prev+5, xend = third_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_prev+5, xend = plot_matamal/365, yend = max_prev,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_prev,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_prev,alpha = 0.05)

#distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#a6d854')

prev_plot_pres <- ggplot(df_distr, aes(x = (t - start)/365, y = slide_prev0to80*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 18)+
  #ylim(0, 75)+
  scale_color_manual(name = "Scenario", values = c("baseline-Sen" = distr_pals[4],
                                                   "all-in-one-stag-Sen" = distr_pals[3],
                                                   "10d-stagger-Sen" = distr_pals[1],
                                                   "20d-stagger-Sen" = distr_pals[2]),
                     labels = c("baseline-Sen" = "Baseline", "all-in-one-stag-Sen" = "1 day to complete MDA",
                                "10d-stagger-Sen" = "10 days to complete MDA", "20d-stagger-Sen" = "20 days to complete MDA"),
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"))+
  theme(legend.position = c(0.4, 0.2))+
  coord_cartesian(xlim = c(-0.25,1))+
  xlab("Years since intervention started")+
  ylab("Slide prevalence (%) in children \n under 5-years-old") +
  geom_segment(x = 0, y = max_prev+5, xend = 0, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_prev+5, xend = second_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_prev+5, xend = third_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_prev+5, xend = plot_matamal/365, yend = max_prev,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_prev,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_prev,alpha = 0.05)


max_inc <- df_distr %>%
  filter(model_type == "baseline-Sen") %>%
  summarise(max_inc = max(clin_inc0to5)*1000)
max_inc <- max_inc$max_inc

inc_plot <- ggplot(df_distr, aes(x = (t - start)/365, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  guides(col = "none")+

  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  coord_cartesian(xlim = c(-0.25,1))+
  xlab("Years since intervention started")+
  geom_segment(x = 0, y = max_inc+1, xend = 0, yend =max_inc, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_inc+1, xend = second_mda/365,
               yend =max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_inc+1, xend = third_mda/365, yend =max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_inc+1, xend = plot_matamal/365, yend = max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_inc,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_inc,alpha = 0.05)

dynamics <- cowplot::plot_grid(mv_plot, eir_plot, inc_plot, prev_plot, align = "v",
                               labels = c("A", "B", "C", "D"))


df_distr_wide_setting_ivm <- df_distr_all %>%
  select(t,ref, init_EIR, model_type, ivm_cov_par, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  }) %>%
  select(-`baseline-Sen`)

df_distr_wide_setting_baseline <- df_distr_all %>%
  filter(model_type == "baseline-Sen") %>%
  select(t,  ref,init_EIR, model_type, clin_inc0to5) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = clin_inc0to5)
  })


df_distr_wide_setting <- left_join(df_distr_wide_setting_ivm,df_distr_wide_setting_baseline,
                                   by = c("ref","t","init_EIR"))

df_distr_wide_setting_prev_ivm<- df_distr_all %>%
  select(t, ref, init_EIR, model_type, ivm_cov_par, slide_prev0to80) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to80)
  }) %>%
  select(-`baseline-Sen`)

df_distr_wide_setting_prev_baseline <- df_distr_all %>%
  filter(model_type == "baseline-Sen") %>%
  select(t,ref,init_EIR, model_type, slide_prev0to80) %>%
  group_by(ref) %>%
  group_modify(~ {
    df_wide <- .x %>%
      pivot_wider(names_from = model_type, values_from = slide_prev0to80)
  })

df_distr_wide_setting_prev <- left_join(df_distr_wide_setting_prev_ivm,df_distr_wide_setting_prev_baseline,
                                        by = c("ref","t","init_EIR"))

impact_overall <- df_distr_wide_setting %>%
  group_by(ref, init_EIR, ivm_cov_par) %>%
  filter(between(t, start, start + 365)) %>% #sum cases one year from start
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_20d = sum(`20d-stagger-Sen`)*1000,
            tot_cases_all_in_stag = sum(`all-in-one-stag-Sen`)*1000) %>%
  mutate(
    impact_all_in_stag = ((tot_cases_baseline - tot_cases_all_in_stag)/tot_cases_baseline)*100,
    impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
    impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
    scenario = "Once year since start") #measuring overall impact



impact_bohemia <- df_distr_wide_setting %>%
  group_by(ref, init_EIR, ivm_cov_par) %>%
  filter(between(t, start, bohemia_inc_period)) %>%
  summarise(tot_cases_baseline = sum(`baseline-Sen`)*1000,
            tot_cases_10d = sum(`10d-stagger-Sen`)*1000,
            tot_cases_20d = sum(`20d-stagger-Sen`)*1000,
            tot_cases_all_in_stag = sum(`all-in-one-stag-Sen`)*1000) %>%
  mutate(
    impact_all_in_stag = ((tot_cases_baseline - tot_cases_all_in_stag)/tot_cases_baseline)*100,
    impact_20d = ((tot_cases_baseline - tot_cases_20d)/tot_cases_baseline)*100,
    impact_10d = ((tot_cases_baseline - tot_cases_10d)/tot_cases_baseline)*100,
    scenario = "bohemia") #measuring impact across 6 months

impact_matamal <- df_distr_wide_setting_prev %>%
  group_by(ref, ivm_cov_par, init_EIR) %>%
  filter(t == matamal_survey) %>%
  summarise(prev_baseline = `baseline-Sen`,
            prev_10d =`10d-stagger-Sen`,
            prev_20d = `20d-stagger-Sen`,
            prev_all_in_stag = `all-in-one-stag-Sen`) %>%
  mutate(
    impact_all_in_stag = ((prev_baseline - prev_all_in_stag)/prev_baseline)*100,
    impact_20d = ((prev_baseline - prev_20d)/prev_baseline)*100,
    impact_10d = ((prev_baseline - prev_10d)/prev_baseline)*100,
    scenario = "matamal") #measuring point prevalence at 4 weeks after last MDA


impact_measurements <- do.call("rbind", list(impact_overall, impact_bohemia, impact_matamal)) %>%
  select(init_EIR, ivm_cov_par,impact_all_in_stag, impact_20d, impact_10d, scenario)

impact_measurements_long <- impact_measurements %>%
  pivot_longer(cols = c(impact_all_in_stag, impact_20d, impact_10d), names_to = "intervention", values_to = "impact")

levels_x <- unique(impact_measurements_long$scenario)

#distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
#distr_pals2 <- distr_pals[1:4]
impact_measurements_long$intervention <- factor(impact_measurements_long$intervention, levels = c("impact_all_in_stag",
                                                                                                  "impact_10d",
                                                                                                  "impact_20d"))
cov_error <- impact_measurements_long %>%
  ungroup() %>%
  select(-ref) %>%
  filter(init_EIR == 100 & ivm_cov_par %in% c(0.5, 0.9)) %>%
  pivot_wider(names_from = ivm_cov_par, values_from = impact) %>%
  rename(cov_low_0.5 = `0.5`,
         cov_high_0.9 = `0.9`)

distr_pals2 <- distr_pals[2:4]

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
    data = cov_error,
    aes(
      x = factor(scenario),
      ymin = cov_low_0.5,
      ymax = cov_high_0.9,
      fill = as.factor(intervention)   # match fill to align dodging
    ),
    width = 0.2,
    position = position_dodge(width = 0.9),
    size = 1.1
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = c(0.7, 0.8)) +
  scale_fill_manual(
    values = distr_pals2, name = "Scenario",
    labels = c(
      "10 days to complete monthly MDA",
      "20 days to complete monthly MDA",
      "1 day to complete monthly MDA (original model)",
      "1 day to complete monthly MDA (staggered model)"
    )
  ) +
  xlab("Time of measurement") +
  ylab("Efficacy (%)") +
  scale_x_discrete(labels = c(
    "bohemia" = "Incidence U5s (start to 6m later)",
    "matamal" = "All-age prevalence (1m after last MDA)",
    "Once year since start" = "Incidence U5s (start to 1y later)"
  ))+
  guides(fill = "none")


impact_measurements_long2 <- impact_measurements_long
impact_measurements_long2$intervention <- factor(impact_measurements_long2$intervention,
                                                 levels = c("impact_all_in_stag",
                                                            "impact_10d",
                                                            "impact_20d"))
cov_error2 <- cov_error
cov_error2$intervention <- factor(cov_error2$intervention,
                                  levels = c("impact_all_in_stag",
                                             "impact_10d",
                                             "impact_20d"))
distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

impact_main_plot_pres <- ggplot() +
  # Bars
  geom_bar(
    data = impact_measurements_long2 %>% filter(init_EIR == 100 & ivm_cov_par == 0.7),
    aes(x = factor(scenario), y = impact, fill = as.factor(intervention)),
    stat = "identity",
    position = position_dodge(width = 0.9)
  ) +
  # Error bars
  geom_errorbar(
    data = cov_error2,
    aes(
      x = factor(scenario),
      ymin = cov_low_0.5,
      ymax = cov_high_0.9,
      fill = as.factor(intervention)   # match fill to align dodging
    ),
    width = 0.2,
    position = position_dodge(width = 0.9),
    size = 1.1
  ) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.7, 0.8)) +
  scale_fill_manual(
    values = c("impact_all_in_stag" = distr_pals[3],
               "impact_10d" = distr_pals[1], "impact_20d" = distr_pals[2]),
    labels = c("impact_all_in_stag" = "1 day to complete MDA",
               "impact_10d" = "10 days to complete MDA",
               "impact_20d" = "20 days to complete MDA"),
    breaks = c("impact_all_in_stag", "impact_10d", "impact_20d"),
    name = "Scenario") +
  xlab("Time of measurement") +
  ylab("Efficacy (%)") +
  scale_x_discrete(labels = c(
    "bohemia" = "Incidence U5s \n (start to 6m later)",
    "matamal" = "Prevalence U5s \n (1m after last MDA)",
    "Once year since start" = "Incidence U5s \n (start to 1y later)"
  ))+
  guides(fill = "none")

impact_measurements_long %>%
  filter(init_EIR == 100 & ivm_cov_par == 0.7) %>%
  group_by(scenario) %>%
  summarise(mean_impact = mean(impact, na.rm = TRUE))
#for 1y since start, mean impact is 26.2%

dynamics_seasonal <- cowplot::plot_grid(mv_plot, eir_plot, prev_plot, inc_plot,
                                        labels = c("A", "B", "C", "D"),
                                        align = "v")

#impact_perennial <- cowplot::plot_grid(impact_cov_plot, impact_Q0_plot, labels = c("E", "F"),
#                                      nrow = 2, align = "v")




plot_seasonal <- cowplot::plot_grid(dynamics_seasonal, impact_main_plot,
                                    labels = c("", "E"))



ggsave(plot_seasonal, file = "2.ivm-stagger-distr/plots/fig_5_plot_seasonal_late.pdf")
ggsave(plot_seasonal, file = "2.ivm-stagger-distr/plots/fig_5_plot_seasonal_late.png")




dynamics_seasonal_pres <- cowplot::plot_grid(mv_plot_pres, prev_plot_pres,
                                             labels = c("A", "B"), align = "h")
plot_seasonal_pres <- cowplot::plot_grid(dynamics_seasonal_pres, impact_main_plot_pres,labels = c("", "C"))
ggsave(plot_seasonal_pres, file = "../glasgow-visit/plot_seasonal_pres_late.png")




#reporting stats
impact_measurements_long %>%
  filter(init_EIR == 100 & ivm_cov_par == 0.7 & scenario == "Once year since start")%>%
  group_by(scenario, intervention)


cov_error %>%
  filter(init_EIR == 100 & scenario == "Once year since start" )

