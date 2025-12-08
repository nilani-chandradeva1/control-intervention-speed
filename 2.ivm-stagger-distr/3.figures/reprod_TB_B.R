library(ggplot2)
library(ggpattern)

getwd()
distr_pals <- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')

impact_df2 <- readRDS("TB_impact_df2.rds")
covs_error_df <- readRDS("TB_covs_error_df.rds")

covs_error_df$seasonality = impact_df2$seasonality

ggplot()  +
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
      col = as.factor(intervention),
      ymin = cov_low_0.5,
      ymax = cov_high_0.9,
      #group = interaction(seasonality, intervention)  # match bars
      group = interaction(seasonality, intervention),
    ),
    width = 0.2,
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
      "impact_all_in_stag"="1 day to complete MDA",
      "impact_10d"="10 days to complete MDA",
      "impact_20d"="20 days to complete MDA"
    ),
    name = "Scenario"
  ) +
  scale_colour_manual(
    values = c("impact_all_in_stag" = distr_pals[3],
               "impact_10d" = distr_pals[1],
               "impact_20d" = distr_pals[2]),
    breaks = c("impact_all_in_stag", "impact_10d", "impact_20d"),
    labels = c(
      "impact_all_in_stag"="1 day to complete MDA",
      "impact_10d"="10 days to complete MDA",
      "impact_20d"="20 days to complete MDA"
    ),
    name = "Scenario"
  )+
  scale_pattern_manual(
    values = c("seasonal-on-time"="none",
               "seasonal-on-late"="stripe"),
    breaks = c("seasonal-on-time","seasonal-on-late"),
    labels = c("On time", "Late (by 2 months)"),
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
  #theme(legend.position = c(c(0.5, 0.9)))+
  ylab("Efficacy (%)") +
  xlab("Time of measurement") +
  scale_x_discrete(labels = c(
    "bohemia" = "Incidence U5s \n (start to 6m later)",
    "matamal" = "All-age prevalence \n (1m after last MDA)",
    "One_year_since_start" = "Incidence U5s \n (start to 1y later)"
  ))
