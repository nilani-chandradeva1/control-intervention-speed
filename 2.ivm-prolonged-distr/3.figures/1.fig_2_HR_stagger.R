#plot of hazard ratio, proportion of treated groups that are lethal and time-varying proportion of population with letal dose of ivm
require(tidyverse)

facet_labels <- c("10d" = "10-day MDA implementation",
                  "20d" = "20-day MDA implementation",
                  "all_in_one_stag" = "Synchronised MDA implementation")

stag_HR <- readRDS("2.ivm-stagger-distr/output/HR_staggered.rds")
stag_cov <- readRDS("2.ivm-stagger-distr/output/prop_lethal_ivm.rds")

stag_HR$stagger <- factor(stag_HR$stagger, levels = c("all_in_one_stag", "10d", "20d"))
stag_cov$stagger <- factor(stag_cov$stagger, levels = c("all_in_one_stag", "10d", "20d"))

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
  labs(linetype = "Average or group-level", y = "Hazard ratio (HR)")+
  #ggtitle("Staggered distribution")+
  #xlim(1, 106)+
  theme(legend.position = "top",
        legend.direction = "horizontal")+
  xlab("Hazard ratio")+
  xlim(-10, 107)+
  ylim(0, 10)


time_toxic_plot <- ggplot(stag_cov, aes(x = Day, y = prop_pop_cov))+
  geom_line(size = 1.1)+
  facet_wrap(vars(stagger), labeller = as_labeller(facet_labels))+
  ylab("Proportion of covered \n group with HR > 1")+
  theme_bw(base_size = 14)+
  #theme(text = element_text(size = 14))+
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.33, 0.66, 1))+
  xlim(-10, 107)+
  xlab("")



#read in one of the model runs to get correct coverage
df_distr_all <- readRDS("2.ivm-stagger-distr/output/df_distr_HR_3m_perennial.rds")
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
covs <- unique(df_distr_all$ivm_cov_par)

df_distr <- df_distr_all %>%
  filter(init_EIR == 100)




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
  "10d-stagger" = "10-day MDA implementation",
  "20d-stagger" = "20-day MDA implementation",
  "all-in-one-stag" = "Synchronised MDA implementation"
)
start <- (365*5)+200

df_time_covs_plot$model_type <- factor(df_time_covs_plot$model_type, levels = c("all-in-one-stag", "10d-stagger", "20d-stagger"))

time_var_cov_plot <- ggplot(df_time_covs_plot, aes(x = t-start, y = ivm_cov*100))+
  geom_line(size = 1.1)+
  geom_ribbon(aes(ymin = lower_cov*100, ymax = upper_cov*100), fill = "lightblue", alpha = 0.4)+
  facet_wrap(vars(model_type), labeller = as_labeller(facet_labels2))+
  xlim(-10, 107)+
  xlab("Time (days) since intervention started")+
  theme_bw(base_size = 14)+
  ylab("Population (%) with \n lethal dose of ivermectin")

toxic_ivm_plot <- cowplot::plot_grid(HR_plot, time_toxic_plot, time_var_cov_plot,
                                     nrow = 3, align = "v",
                                     labels = c("A", "B", "C"))


ggsave(toxic_ivm_plot, file = "2.ivm-stagger-distr/plots/SM_fig_toxic_ivm_plots.pdf")
