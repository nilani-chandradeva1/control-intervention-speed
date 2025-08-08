#generate HR curve for staggered distributions of ivermectin
#ivermectin MDA takes place once a month, for 3 months

#simulate 3 equally interspaced distributions per month, staggered across 10 days or 20 days

#modify the hazard ratio curve
require(tidyverse)

#treating pop once a month for three months
#each monthly MDA either completed in 11 days or 21 days.

ivm_haz <- read.table("IVM_derivation/ivermectin_hazards.txt", header=TRUE) #Smit Hazard Ratios
colnames(ivm_haz) = c("Day", "IVM_400_1_HS", "IVM_300_3_HS")

stag_10 <- 11
times_10_d <- seq(1, stag_10, length.out = 3)

ivm_haz <- ivm_haz %>%
  select(-IVM_400_1_HS) %>%
  filter(between(Day, 1, 23)) %>%
  mutate(IVM_300_3_HS = round(IVM_300_3_HS, 2))

ivm_haz_df <- as.data.frame(ivm_haz)

shift_days_10_MDA_3 <- times_10_d[3] - times_10_d[1]

max_day <- max(ivm_haz_df$Day)
new_days <- (max_day + 1):90 # whole thing will be completed within 90 days (3 months)

extra_rows_10d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_10d <- rbind(ivm_haz_df, extra_rows_10d)

group1_HR <- c(df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7))

df_extended_10d <- df_extended_10d %>%
  mutate(group1 = group1_HR)

#create new shifted column - group 3 - and fill with HR 1 first, no impact if not started yet

df_extended_10d$group3 <- 1 # HR1 - no impact

df_extended_10d$group3[(shift_days_10_MDA_3 + 1):nrow(df_extended_10d)] <- df_extended_10d$group1[1:(nrow(df_extended_10d) - shift_days_10_MDA_3)]

shift_days_10_MDA_2 <- times_10_d[2] - times_10_d[1]

df_extended_10d$group2 <- 1
df_extended_10d$group2[(shift_days_10_MDA_2 + 1):nrow(df_extended_10d)] <- df_extended_10d$group1[1:(nrow(df_extended_10d) - shift_days_10_MDA_2)]

df_extended_10d <- df_extended_10d %>%
  select(Day, IVM_300_3_HS, group1, group2, group3) #these are HRs for each day (repeated for each month)

#take mean HR across daily HR for groups, only if the HR is greater than 1 (i.e. only take mean if toxic)

df_extended_10d_new_HR <- df_extended_10d %>%
  rowwise() %>%
  mutate(HR_use_above1 = mean(c_across(group1:group3)[c_across(group1:group3) > 1], na.rm = TRUE), #rowwise mean, if greater than 1
         prop_pop_cov = sum(c_across(group1:group3) > 1)/3, #the proportion that are toxic (if HR is greater than 1, use in calc)
         stagger = "10d") %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  ungroup()

#for the 20d strategy

stag_20 <- 21
times_20_d <- seq(1, stag_20, length.out = 3)

shift_days_20_MDA_3 <- times_20_d[3] - times_20_d[1]

max_day <- max(ivm_haz_df$Day)
new_days <- (max_day + 1):123

extra_rows_20d <- data.frame(Day = new_days, IVM_300_3_HS = 1)

df_extended_20d <- rbind(ivm_haz_df, extra_rows_20d)

to_add <- nrow(df_extended_20d) - 90

group1_HR <- c(df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               df_extended_10d$IVM_300_3_HS[1:23], rep(1, 7),
               rep(1, to_add))

length(group1_HR)

df_extended_20d <- df_extended_20d %>%
  mutate(group1 = group1_HR)

#create new shifted column - group 3 - and fill with HR 1 first, no impact

df_extended_20d$group3 <- 1 #HR 1 - no impact

df_extended_20d$group3[(shift_days_20_MDA_3 + 1):nrow(df_extended_20d)] <- df_extended_20d$group1[1:(nrow(df_extended_20d) - shift_days_20_MDA_3)]

shift_days_20_MDA_2 <- times_20_d[2] - times_20_d[1]

df_extended_20d$group2 <- 1
df_extended_20d$group2[(shift_days_20_MDA_2 + 1):nrow(df_extended_20d)] <- df_extended_20d$group1[1:(nrow(df_extended_20d) - shift_days_20_MDA_2)]

df_extended_20d <- df_extended_20d %>%
  select(Day, IVM_300_3_HS, group1, group2, group3)

df_extended_20d <- df_extended_20d[1:113,]

df_extended_20d_new_HR <- df_extended_20d %>%
  rowwise() %>%
  mutate(HR_use_above1 = mean(c_across(group1:group3)[c_across(group1:group3) > 1], na.rm = TRUE), #rowwise mean, if greater than 1
         prop_pop_cov = sum(c_across(group1:group3) > 1)/3,
         stagger = "20d") %>%
  filter(if_all(everything(), ~ !is.na(.))) %>%
  ungroup()


#HR_all <- rbind(df_extended_10d_new_HR, df_extended_20d_new_HR)


#for the all-in-one distribution
#as a sanity check, check that implementing the treatment of target coverage works similarly in staggered model and original model
df_all <-df_extended_10d_new_HR %>%
  select(Day, group1) %>%
  rename(HR_use_above1 = group1) %>% #giving this name so easy to bind with others
  mutate(HR_use_above1 = case_when(HR_use_above1 == 1 ~ round(1),
                                   TRUE ~ HR_use_above1)) %>%
  mutate(stagger = "all_in_one_stag") %>%
  rowwise() %>%
  mutate(
    prop_pop_cov = sum(c_across(HR_use_above1) > 1)/1)

plot(df_all$Day[1:23], df_all$HR_use_above1[1:23]) #these will be inputs for all in one modelling

plot(df_all$Day, df_all$HR_use_above1) #these will be inputs for the staggered implementation of the all in one modelling

#confirming that all HRs are the same when repeated
df_all$HR_use_above1[1:23] ==  df_all$HR_use_above1[31:53]
df_all$HR_use_above1[1:23] ==  df_all$HR_use_above1[61:83]
df_all$HR_use_above1[31:53] == df_all$HR_use_above1[61:83]


df_all$prop_pop_cov[1:23] ==  df_all$prop_pop_cov[31:53]
df_all$prop_pop_cov[1:23] ==  df_all$prop_pop_cov[61:83]
df_all$prop_pop_cov[31:53] == df_all$prop_pop_cov[61:83]


df_all_HR <- df_all %>%
  mutate(
    group2 = NA,
    group3 = NA,
    group1 = HR_use_above1)

df_all_HR$IVM_300_3_HS <- df_extended_10d_new_HR$IVM_300_3_HS

df_all_HR <- df_all_HR %>%
  select(Day, IVM_300_3_HS, group1, group2, group3, HR_use_above1, prop_pop_cov,stagger)


HR_all <- do.call("rbind", list(df_extended_10d_new_HR, df_extended_20d_new_HR, df_all_HR))


head(HR_all)

HR_all_long <- HR_all %>%
  pivot_longer(cols = c(group1, group2, group3, HR_use_above1), names_to = "group",
               values_to = "HR") %>%
  mutate(group_lab = case_when(group == "HR_use_above1" ~ "average",
                               TRUE ~ "group-level"),
         group_combo = paste0(group_lab))

HR_plot <- ggplot(HR_all_long, aes(x = Day, y = HR, group = group))+
  geom_line(aes(col = as.factor(group), lty = as.factor(group_lab)),size = 1.1,
            inherit.aes = TRUE)+
  geom_point(size = 2,aes(col = as.factor(group)), alpha = 0.5)+
  facet_wrap(vars(stagger), labeller = label_both)+
  theme_bw()+
  geom_hline(aes(yintercept = 1), lty = "dashed")+
  ylim(1, 10)+
  scale_colour_manual(values = c(
    "group1" = "#1b9e77",
    "group2" = "#d95f02",
    "group3" = "#7570b3",
    "HR_use_above1" = "black"), labels= c("Group 1", "Group 2", "Group 3", "Group average"),
    name = "Hazard ratio group"
  )+
  labs(linetype = "Average or group-level", y = "Hazard ratio")+
  ggtitle("HR and distribution strategy")


time_cov_plot <- ggplot(HR_all, aes(x = Day, y = prop_pop_cov))+
  geom_line()+
  facet_wrap(vars(stagger), labeller = label_both)+
  ylab("Proportion of covered group with HR > 1")+
  theme_bw()+
  ylim(0,1)

write_rds(HR_all_long, file = "analysis/target-profiles-distrib-strat/chapter/output/test_HR_staggered_long.rds")
write_rds(HR_all, file = "analysis/target-profiles-distrib-strat/chapter/output/test_HR_staggered.rds")

cov_in <- 0.7 #if 70% of the pop are treated and everyone is eligible

prop_pop_toxic <- ggplot(HR_all, aes(x = Day, y = prop_pop_cov*cov_in))+
  geom_line()+
  facet_wrap(vars(stagger), labeller = label_both)+
  ylab("Proportion of population \n with lethal dose")+
  theme_bw()+
  ylim(0,1)

#check coverages and HR line up as expected
cowplot::plot_grid(HR_plot,time_cov_plot,prop_pop_toxic,
                   nrow = 3, align = "v")
