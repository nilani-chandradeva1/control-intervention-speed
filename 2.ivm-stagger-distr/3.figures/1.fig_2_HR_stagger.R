#plot of hazard ratio, proportion of treated groups that are lethal and time-varying proportion of population with letal dose of ivm


require(tidyverse)

facet_labels <- c("10d" = "10-day staggered MDA",
                  "20d" = "20-day staggered MDA",
                  "all_in_one_stag" = "Overnight MDA")

stag_HR <- readRDS("2.ivm-stagger-distr/output/HR_staggered.rds")
stag_cov <- readRDS("2.ivm-stagger-distr/output/prop_lethal_ivm.rds")

HR_plot <- ggplot(stag_HR, aes(x = Day, y = HR, group = group))+
  geom_line(aes(col = as.factor(group), lty = as.factor(group_lab)),size = 1.1,
            inherit.aes = TRUE, alpha = 0.6)+
  geom_point(size = 1.5, alpha = 0.8, aes(col = as.factor(group)))+
  facet_wrap(vars(stagger), labeller = as_labeller(facet_labels))+
  theme_bw(base_size = 14)+
  geom_hline(aes(yintercept = 1), lty = "dashed")+
  ylim(1, 16)+
  scale_colour_manual(values = c(
    "group1" = "#1b9e77",
    "group2" = "#d95f02",
    "group3" = "#7570b3",
    "HR_use_above1" = "black"), labels= c("Group 1", "Group 2", "Group 3", "Group average"),
    name = "Hazard ratio group"
  )+
  labs(linetype = "Average or group-level", y = "Hazard ratio")+
  #ggtitle("Staggered distribution")+
  #xlim(1, 106)+
  theme(legend.position = c(0.25, 0.7),
        legend.direction = "horizontal")+
  xlab("Hazard ratio")+
  xlim(-10, 107)+
  ylim(0, 24)


time_toxic_plot <- ggplot(stag_cov, aes(x = Day, y = prop_pop_cov))+
  geom_line(size = 1.1)+
  facet_wrap(vars(stagger), labeller = as_labeller(facet_labels))+
  ylab("Proportion of covered group \n with HR > 1")+
  theme_bw(base_size = 14)+
  #theme(text = element_text(size = 14))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.33, 0.66, 1))+
  xlim(-10, 107)+
  xlab("")

#cowplot::plot_grid(HR_plot, time_toxic_plot, nrow = 2)

df_time_covs <- df_distr %>%
  #select(t, ivm_cov, model_type) %>%
  filter(model_type %in% c("10d-stagger", "20d-stagger", "all-in-one-stag") & ivm_cov_par == 0.7) %>%
  select(t, model_type, ivm_cov_par, ivm_cov)

df_time_covs_error <- df_distr %>%
  select(t, model_type, ivm_cov_par, ivm_cov) %>%
  filter(model_type %in% c("10d-stagger", "20d-stagger", "all-in-one-stag") & ivm_cov_par %in% c(0.5, 0.9))%>%
  pivot_wider(names_from = ivm_cov_par, values_from=ivm_cov) %>%
  rename(lower_cov = `0.5`,
         upper_cov = `0.9`)

df_time_covs_plot <- left_join(df_time_covs, df_time_covs_error)
range(df_time_covs_error$lower_cov)

df_time_covs_plot

facet_labels2 <- c(
  "10d-stagger" = "10-day staggered MDA",
  "20d-stagger" = "20-day staggered MDA",
  "all-in-one-stag" = "Overnight MDA"
)

time_var_cov_plot <- ggplot(df_time_covs_plot, aes(x = t-start, y = ivm_cov*100))+
  geom_line(size = 1.1)+
  geom_ribbon(aes(ymin = lower_cov*100, ymax = upper_cov*100), fill = "lightblue", alpha = 0.4)+
  facet_wrap(vars(model_type), labeller = as_labeller(facet_labels2))+
  xlim(-10, 107)+
  xlab("Time (days) since intervention started")+
  theme_bw(base_size = 14)+
  ylab("Proportion of population with \n lethal dose of ivermectin (%)")

toxic_ivm_plot <- cowplot::plot_grid(HR_plot, time_toxic_plot, time_var_cov_plot,
                                     nrow = 3, align = "v",
                                     labels = c("A", "B", "C"))

ggsave(toxic_ivm_plot, file = "2.ivm-stagger-distr/plots/fig_2_toxic_ivm_plots.pdf")
