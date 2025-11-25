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
  xlab("Years since start of well-timed MDA implementation") +
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
  xlab("Years since start of well-timed MDA implementation") +
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
  ylab("Average number of infectious bites \n per person per day (daily EIR)")+
  xlab("Years since start of well-timed MDA implementation") +
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
  xlab("Years since start of well-timed MDA implementation") +
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
