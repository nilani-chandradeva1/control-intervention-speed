require(tidyverse)

#Used to generate estimates of entomological efficacy in Table 1

model_out_constant <-readRDS("1.simple-seirs-sei/output/model_results_df_constant_emergence_TRUE.rds")
model_out_constant_base <- readRDS("1.simple-seirs-sei/output/model_results_base_df_constant_emergence_TRUE.rds")

start_base <- 0

time_periods <- unique(model_out_constant$delta_t)

unique(model_out_constant$M0)

base_10d <- model_out_constant_base %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & M0 == 2000) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[1]) %>% #for each time period, summarise the number of mosquitoes that have ever-lived, naturally died and were killed by the intervention
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_30d <- model_out_constant_base %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & M0 == 2000) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[2]) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_90d <- model_out_constant_base %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & M0 == 2000) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[3]) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_250d <- model_out_constant_base %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & M0 == 2000) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + 250) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

#model_out_constant %>%
#  mutate(prop_killed_all = delta_D/M0) %>%
#  filter(m0 == 2 & prop_killed_all == 0.95 & M0 == 2000) %>%
#  ggplot(aes(x = (t-start_int), y = M, col = as.factor(delta_t)))+
#  geom_line()+
#  coord_cartesian(xlim = c(0,10), ylim = c(0,2000))


2000 + ((2000*0.1)*10)

start_int <- 100

#measurement at 10 days after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & prop_killed_all == 0.95 & M0 == 2000) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_int + time_periods[1]) %>%
  summarise(#ever_lived = ever_lived,
            #nat_die = nat_die,
            D = round(D))%>%
  mutate(base_ever_lived = base_10d$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 30 days after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + time_periods[2]) %>%
  summarise(#ever_lived = ever_lived,
    #nat_die = nat_die,
    D = round(D))%>%
  mutate(base_ever_lived = base_30d$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 90 days after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + time_periods[3]) %>%
  summarise(
            D = round(D))%>%
  mutate(base_ever_lived = base_90d$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 250 days after start

model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 2 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + 250) %>%
  summarise(
            D = round(D)) %>%
  mutate(base_ever_lived = base_250d$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)


#at higher V-H ratio
#measurement at 10 days after start

base_10d_high_vh <- model_out_constant_base %>%
  #mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & M0 == 2e5) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[1]) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_30d_high_vh <- model_out_constant_base %>%
  #mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & M0 == 2e5) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[2]) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_90d_high_vh <- model_out_constant_base %>%
  #mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & M0 == 2e5) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + time_periods[3]) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)

base_250d_high_vh <- model_out_constant_base %>%
  #mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & M0 == 2e5) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_base + 250) %>%
  summarise(ever_lived = ever_lived,
            nat_die = nat_die,
            D = round(D),
            prop_killed = ((D/(ever_lived))*100),
            M = M,
            M0 = M0)
#add on stuff here

#measure at 10d after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & prop_killed_all == 0.95 & M0 == 2e5) %>%
  group_by(delta_t, m0, constant_emergence) %>%
  filter(t == start_int + time_periods[1]) %>%
  summarise(#ever_lived = ever_lived,
    #nat_die = nat_die,
    D = round(D))%>%
  mutate(base_ever_lived = base_10d_high_vh$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 30 days after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + time_periods[2]) %>%
  summarise(#ever_lived = ever_lived,
    #nat_die = nat_die,
    D = round(D))%>%
  mutate(base_ever_lived = base_30d_high_vh$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 90 days after start
model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + time_periods[3]) %>%
  summarise(#ever_lived = ever_lived,
    #nat_die = nat_die,
    D = round(D))%>%
  mutate(base_ever_lived = base_90d_high_vh$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

#measurement at 250 days after start

model_out_constant %>%
  mutate(prop_killed_all = delta_D/M0) %>%
  filter(m0 == 200 & prop_killed_all == 0.95) %>%
  group_by(delta_t, m0) %>%
  filter(t == start_int + 250) %>%
  summarise(#ever_lived = ever_lived,
    #nat_die = nat_die,
    D = round(D))%>%
  mutate(base_ever_lived = base_250d_high_vh$ever_lived,
         prop_killed_int = (D/base_ever_lived)*100)

