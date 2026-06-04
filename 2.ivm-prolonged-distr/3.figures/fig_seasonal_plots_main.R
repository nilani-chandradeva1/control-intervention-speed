require(tidyverse)
df_distr_seasonal <- readRDS("2.ivm-stagger-distr/output/df_distr_season_timed.rds")
df_distr_seasonal_late <- readRDS("2.ivm-stagger-distr/output/df_distr_season_late.rds")

distr_pals <- c('#a6d854', '#8da0cb', '#66c2a5', '#fc8d62')

start <- (365*5)+200
matamal_survey <- (30*4) + start
plot_matamal <- matamal_survey-start #diff between start and survey
bohemia_inc_period <- start+(6*30)
plot_bohemia <- bohemia_inc_period- start
second_mda <- 30+start-start
third_mda <- 60+start-start

max_mv <- max(df_distr_seasonal$mv, na.rm = TRUE)
mv_plot_on_time <- ggplot(df_distr_seasonal, aes(x = (t - start)/365, y = mv, col = as.factor(model_type))) +
  geom_line(size = 1.1) +
  theme_bw(base_size = 14) +
  #ylim(0, 50) +
  ylab("Mosquito density")+
  #theme(legend.position = c(0.7, 0.3)) +
  guides(col = "none")+
  #coord_cartesian(ylim = c(0,200))+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 200))+
  xlab("Years since start of on-time MDA") +
  geom_segment(x = 0, y = max_mv+25, xend = 0, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #first MDA
  geom_segment(x = second_mda/365, y = max_mv+25, xend = second_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+ #second MDA
  geom_segment(x = third_mda/365, y = max_mv+25, xend = third_mda/365, yend = max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1) #+ #third MDA
  #geom_segment(x = plot_matamal/365, y = max_mv+25, xend = plot_matamal/365, yend = max_mv,
  #             arrow = arrow(length = unit(0.3, "cm")),
  #             col = "blue", size = 1.1)+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_mv,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_mv,
  #         alpha = 0.05)

max_prev <- max(df_distr_seasonal$slide_prev0to80*100, na.rm = TRUE)

prev_plot_on_time <- ggplot(df_distr_seasonal, aes(x = (t - start)/365, y = slide_prev0to80*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  #ylim(0, 75)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  theme(legend.position = c(0.45, 0.2))+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(-15, 70))+
  xlab("Years since start of on-time MDA") +
  ylab("All-age slide prevalence (%)") +
  geom_segment(x = 0, y = max_prev+10, xend = 0, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_prev+10, xend = second_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_prev+10, xend = third_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_prev+10, xend = plot_matamal/365, yend = max_prev,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "blue", size = 1.1) #+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_prev,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_prev,alpha = 0.05)


#late plots
start <- (365*5)+200
start_late <- start + 60 #update start time: 60 days later
matamal_survey <- (30*4) + start_late
plot_matamal <- matamal_survey-start #diff between start and survey
bohemia_inc_period <- start_late+(6*30)
plot_bohemia <- bohemia_inc_period- start_late
second_mda <- 30+start_late-start
third_mda <- 60+start_late-start

max_mv <- max(df_distr_seasonal_late$mv, na.rm = TRUE)

mv_plot_late <- ggplot(df_distr_seasonal_late, aes(x = (t - start)/365, y = mv, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  guides(col = "none")+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 200))+
  ylab("Mosquito density")+
  xlab("Years since start of on-time MDA") +
  geom_segment(x = 60/365, y = max_mv+25, xend = 60/365, yend =max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_mv+25, xend = second_mda/365, yend =max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_mv+25, xend = third_mda/365, yend =max_mv, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)#+
#geom_segment(x = plot_matamal/365, y = max_eir+0.05, xend = plot_matamal/365, yend = max_eir,
#             arrow = arrow(length = unit(0.3, "cm")),
#             col = "blue", size = 1.1)+
#annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_eir,
#         fill = "white", alpha = 0.1, col = "black")+
#annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_eir,
#         alpha = 0.05)

max_prev <- max(df_distr_seasonal_late$slide_prev0to80*100, na.rm = TRUE)

prev_plot_late <- ggplot(df_distr_seasonal_late, aes(x = (t - start)/365, y = slide_prev0to80*100, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("No intervention", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Scenario") +
  guides(col = "none")+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(-15, 70))+
  ylab("All-age prevalence (%)")+
  xlab("Years since start of on-time MDA") +
  geom_segment(x = 60/365, y = max_prev+10, xend = 60/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_prev+10, xend = second_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_prev+10, xend = third_mda/365, yend =max_prev, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = plot_matamal/365, y = max_prev+10, xend = plot_matamal/365, yend = max_prev,
             arrow = arrow(length = unit(0.3, "cm")),
             col = "blue", size = 1.1)#+
#annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_eir,
#         fill = "white", alpha = 0.1, col = "black")+
#annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_eir,
#         alpha = 0.05)

dynamics <- cowplot::plot_grid(mv_plot_on_time, prev_plot_on_time,
                      mv_plot_late, prev_plot_late,
                      labels = c("A", "B", "C", "D"))


#impact
impact_on_time <- readRDS("2.ivm-stagger-distr/output/impact_seasonal_on_time.rds")
covs_error_on_time <- readRDS("2.ivm-stagger-distr/output/impact_seasonal_covs_error.rds") %>%
  mutate(seasonality = "seasonal-on-time"  )


#impact-late
impact_late <- readRDS("2.ivm-stagger-distr/output/impact_seasonal_late.rds")
cov_error_late <- readRDS("2.ivm-stagger-distr/output/impact_seasonal_late_covs_error.rds")%>%
  mutate(seasonality = "seasonal-on-late")


impact_df <- rbind(impact_on_time, impact_late)
covs_error_df <- rbind(covs_error_on_time, cov_error_late) %>%
  mutate(scenario = case_when(scenario == "Once year since start" ~ "One_year_since_start",
                              TRUE ~ scenario))

covs_error_df$intervention <- factor(covs_error_df$intervention, levels = c("impact_all_in_stag",
                                                                            "impact_10d",
                                                                            "impact_20d"))

covs_error_df$seasonality <- factor(covs_error_df$seasonality, levels = c("seasonal-on-time","seasonal-on-late"))

impact_df2 <- impact_df %>%
  filter(init_EIR == 100 & ivm_cov_par == 0.7) %>%
  mutate(scenario = case_when(scenario == "Once year since start" ~ "One_year_since_start",
                              TRUE ~ scenario))

impact_df2$seasonality <- factor(impact_df2$seasonality, levels = c("seasonal-on-time","seasonal-on-late"))
impact_df2$intervention <- factor(impact_df2$intervention, levels = c("impact_all_in_stag", "impact_10d", "impact_20d"))


require(ggpattern)


distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

#saveRDS(impact_df2, file = "2.ivm-stagger-distr/output/TB_impact_df2.rds")
#saveRDS(covs_error_df, file = "2.ivm-stagger-distr/output/TB_covs_error_df.rds")

impact_plot <- ggplot()  +
  geom_bar_pattern(data = impact_df2,
                   aes(
                     x = factor(scenario),
                     y = impact,
                     fill = as.factor(intervention),
                     pattern = factor(seasonality,
                                      levels = c("seasonal-on-time","seasonal-on-late"))
                   ),
                   stat = "identity",
                   position = position_dodge(preserve = "single", 0.9),
                   #position = position_dodge(),
                   color = "black",
                   pattern_colour = "black",
                   pattern_fill = NA,
                   pattern_density = 0.4,
                   pattern_spacing = 0.04,
                   pattern_key_scale_factor = 0.5
  ) +
  geom_errorbar(
    data = covs_error_df,
    aes(
      x = factor(scenario),
      #col = as.factor(intervention),
      ymin = cov_low_0.5,
      ymax = cov_high_0.9,
      #group = interaction(seasonality, intervention)  # match bars
      group = interaction(seasonality, intervention),
    ),
    width = 0.2,
    col = "black",
    #position = position_dodge(),
    position = position_dodge(preserve = "single", 0.9),
    size = 1.1
  ) +
  scale_fill_manual(
    values = c("impact_all_in_stag" = distr_pals[3],
               "impact_10d" = distr_pals[1],
               "impact_20d" = distr_pals[2]),
    breaks = c("impact_all_in_stag", "impact_10d", "impact_20d"),
    labels = c(
      "impact_all_in_stag"="Synchronised MDA",
      "impact_10d"="10 days to complete MDA",
      "impact_20d"="20 days to complete MDA"
    ),
    name = "Scenario"
  ) +
 #scale_colour_manual(
 #  values = c("impact_all_in_stag" = distr_pals[3],
 #             "impact_10d" = distr_pals[1],
 #             "impact_20d" = distr_pals[2]),
 #  breaks = c("impact_all_in_stag", "impact_10d", "impact_20d"),
 #  labels = c(
 #    "impact_all_in_stag"="1 day to complete MDA",
 #    "impact_10d"="10 days to complete MDA",
 #    "impact_20d"="20 days to complete MDA"
 #  ),
 #  name = "Scenario"
 #)+
  scale_pattern_manual(
    values = c("seasonal-on-time"="none",
               "seasonal-on-late"="stripe"),
    breaks = c("seasonal-on-time","seasonal-on-late"),
    labels = c("On-time", "Late (by 2 months)"),
    name = "Timing of MDA in seasonal setting"
  ) +
  guides(
    fill = "none",
    pattern = guide_legend(
      override.aes = list(
        fill = "white",
        pattern_fill = "white",
        pattern_colour = "black"
      )
    )
  )+
  theme_bw(base_size = 14)+
  theme(legend.position = c(c(0.5, 0.9)))+
  ylab("Efficacy (%)") +
  xlab("Time of measurement") +
  scale_x_discrete(labels = c(
    "bohemia" = "Incidence U5s \n (start to 6m later)",
    "matamal" = "All-age prevalence \n (1m after last MDA)",
    "One_year_since_start" = "Incidence U5s \n (start to 1y later)"
  ))

seasonal_figure <- cowplot::plot_grid(dynamics, impact_plot, labels = c("", "E"))





ggsave(seasonal_figure, file = "2.ivm-stagger-distr/plots/seasonality_main_figure.png")
#ggsave(seasonal_figure, file = "2.ivm-stagger-distr/plots/seasonality_main_figure.pdf")
