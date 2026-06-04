require(tidyverse)

df_distr_seasonal <- readRDS("2.ivm-stagger-distr/output/df_distr_season_timed.rds")
df_distr_seasonal_late <- readRDS("2.ivm-stagger-distr/output/df_distr_season_late.rds")

max_eir <- max(df_distr_seasonal$EIRout, na.rm = TRUE)

distr_pals <- c('#a6d854', '#8da0cb', '#66c2a5', '#fc8d62')

start <- (365*5)+200
matamal_survey <- (30*4) + start
plot_matamal <- matamal_survey-start #diff between start and survey
bohemia_inc_period <- start+(6*30)
plot_bohemia <- bohemia_inc_period- start
second_mda <- 30+start-start
third_mda <- 60+start-start

eir_plot_s_on_time <- ggplot(df_distr_seasonal, aes(x = (t - start)/365, y = EIRout, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  #guides(col = "none")+
  theme(legend.position = c(0.6, 0.5))+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 1))+
  ylab("Average number of infectious bites \n per person per day (daily EIR)")+
  xlab("Years since start of well-timed MDA implementation") +
  geom_segment(x = 0, y = max_eir+0.15, xend = 0, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_eir+0.15, xend = second_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_eir+0.15, xend = third_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)#+
  #geom_segment(x = plot_matamal/365, y = max_eir+0.15, xend = plot_matamal/365, yend = max_eir,
  #             arrow = arrow(length = unit(0.3, "cm")),
  #             col = "blue", size = 1.1)+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_eir,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_eir,
  #         alpha = 0.05)

max_inc <- df_distr_seasonal %>%
  filter(model_type == "baseline-Sen") %>%
  summarise(max_inc = max(clin_inc0to5)*1000)
max_inc <- max_inc$max_inc

inc_plot_s_on_time <- ggplot(df_distr_seasonal, aes(x = (t - start)/365, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  guides(col = "none")+

  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  coord_cartesian(xlim = c(-0.25,1))+
  xlab("Years since start of well-timed MDA implementation") +
  geom_segment(x = 0, y = max_inc+5, xend = 0, yend =max_inc, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_inc+5, xend = second_mda/365,
               yend =max_inc, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_inc+5, xend = third_mda/365, yend =max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  #geom_segment(x = plot_matamal/365, y = max_inc+5, xend = plot_matamal/365, yend = max_inc,
  #             arrow = arrow(length = unit(0.3, "cm")),
  #             col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_inc,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_inc,alpha = 0.05)



#late plots
start <- (365*5)+200
start_late <- start + 60 #update start time: 60 days later
matamal_survey <- (30*4) + start
#plot_matamal <- matamal_survey-start #diff between start and survey
bohemia_inc_period <- start_late+(6*30)
plot_bohemia <- bohemia_inc_period- start_late
second_mda <- 30+start_late-start
third_mda <- 60+start_late-start

max_eir <- max(df_distr_seasonal_late$EIRout, na.rm = TRUE)

eir_plot_late <- ggplot(df_distr_seasonal_late, aes(x = (t - start)/365, y = EIRout, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  guides(col = "none")+
  coord_cartesian(xlim = c(-0.25,1), ylim = c(0, 1))+
  ylab("Average number of infectious bites \n per person per day (daily EIR)")+
  xlab("Years since start of well-timed MDA implementation") +
  geom_segment(x = 60/365, y = max_eir+0.15, xend = 60/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_eir+0.15, xend = second_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_eir+0.15, xend = third_mda/365, yend =max_eir, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)#+
  #geom_segment(x = plot_matamal/365, y = max_eir+0.05, xend = plot_matamal/365, yend = max_eir,
  #             arrow = arrow(length = unit(0.3, "cm")),
  #             col = "blue", size = 1.1)+
  #annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_eir,
  #         fill = "white", alpha = 0.1, col = "black")+
  #annotate("rect", xmin = 0, xmax = 1, fill = "blue", ymin = 0, ymax = max_eir,
  #         alpha = 0.05)

max_inc <- df_distr_seasonal_late %>%
  filter(model_type == "baseline-Sen") %>%
  summarise(max_inc = max(clin_inc0to5)*1000)
max_inc <- max_inc$max_inc

inc_plot_late <- ggplot(df_distr_seasonal_late, aes(x = (t - start)/365, y = clin_inc0to5*1000, col = as.factor(model_type)))+
  geom_line(size = 1.1)+
  theme_bw(base_size = 14)+
  scale_color_manual(values = distr_pals,
                     breaks = c("baseline-Sen", "all-in-one-stag-Sen", "10d-stagger-Sen", "20d-stagger-Sen"),
                     labels = c("Baseline", "Synchronised MDA", "10-day MDA", "20-day MDA"),
                     name = "Time to complete MDA") +
  guides(col = "none")+

  ylab("Clinical incidence in children \n under 5-years-old, per 1000 persons")+
  coord_cartesian(xlim = c(-0.25,1))+
  xlab("Years since start of well-timed MDA implementation") +
  geom_segment(x = 60/365, y = max_inc+5, xend = 60/365, yend =max_inc, arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = second_mda/365, y = max_inc+5, xend = second_mda/365,
               yend =max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  geom_segment(x = third_mda/365, y = max_inc+5, xend = third_mda/365, yend =max_inc,
               arrow = arrow(length = unit(0.3, "cm")),
               col = "black", size = 1.1)+
  #geom_segment(x = plot_matamal/365, y = max_inc+1, xend = plot_matamal/365, yend = max_inc,
   #            arrow = arrow(length = unit(0.3, "cm")),
    #           col = "blue", size = 1.1)+
  annotate("rect", xmin = 0, xmax = plot_bohemia/365, ymin = 0, ymax = max_inc,
           fill = "white", alpha = 0.1, col = "black")+
  annotate("rect", xmin = 60/365, xmax = 1, fill = "blue", ymin = 0, ymax = max_inc,alpha = 0.05)

SM_seasonal_mda_plots <- cowplot::plot_grid(eir_plot_s_on_time, inc_plot_s_on_time,
                   eir_plot_late, inc_plot_late,
                   labels = c("A", "B", "C", "D"))






ggsave("2.ivm-stagger-distr/plots/SM_seasonal_mda_plots.png")
