library(dplyr)
library(tidyr)
library(lmerTest)
library(lsmeans)
library(ggplot2)
library(knitr)
source("read_data.R")


phage %>%
  filter(ph=="6.5") -> phage_temperature

phage %>%
  filter(temperature=="6") -> phage_ph

#######################Temperature model#############################

lmer(log_count ~ temperature*Phage + day + Strain + Phage*Strain +
       milk_age +
       log_milk_apc +
       (1|milk_batch/cheese_make/cheese_plate) +
       (1|rep),
      data=phage_temperature) -> m_temperature

#Raw data used in model above
phage_temperature %>%
  select(Phage,temperature,Strain,day,experiment,log_count) %>%
  spread(key = Phage,value=log_count) %>%
  mutate(phage_effect=N-Y) %>%
  filter(!is.na(phage_effect)) %>%
  group_by(Strain,temperature,day) %>%
  summarize(mean_phage_effect=mean(phage_effect),
            se_phage_effect=sd(phage_effect)/sqrt(n())) %>%
  ggplot(aes(x=day,y=mean_phage_effect)) +
  geom_hline(aes(yintercept=0),color="red",alpha=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  labs(y="Phage Effect (Log)", x="Day") +
  geom_errorbar(aes(ymin=mean_phage_effect-2*se_phage_effect,
                    ymax=mean_phage_effect+2*se_phage_effect),
                width=.1, position=position_dodge(width=0.5)) +
  ylim(c(-3,5.5)) +
  facet_grid(Strain~temperature)
  
# Full plot for temperature model
lsm_temperature <- summary(lsmeans(m_temperature,~Phage+temperature+day))
lsm_temperature %>%
  data.frame() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=lsmean,color=Phage)) +
  geom_point(position=position_dodge(1)) +
  labs(y="Least Square Means", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_errorbar(width=1,
                aes(ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(~temperature)

#Raw data of plot above
phage_temperature %>%
  filter(!is.na(log_count)) %>%
  group_by(day,temperature,Phage) %>%
  summarize(mean_log_count=mean(log_count),
            se_log_count=sd(log_count)/sqrt(n())) %>%
  ungroup() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=mean_log_count,color=Phage)) +
  labs(y="Average Log Count", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=mean_log_count-2*se_log_count,
                    ymax=mean_log_count+2*se_log_count),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(~temperature)

# Full temperature model (same as above) separated by strain
tiff(filename = "Phage_Temp_Strain_Model.tiff", units = "in", width = 7, height = 5, res = 300)
lsm_temperature_w_strain <- summary(lsmeans(m_temperature,~Phage+temperature+day+Strain))
lsm_temperature_w_strain %>%
  data.frame() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=lsmean,color=Phage)) +
  labs(y="Least Square Means", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(Strain~temperature)

# Raw data of above plot
phage_temperature %>%
  filter(!is.na(log_count)) %>%
  group_by(day,temperature,Phage,Strain) %>%
  summarize(mean_log_count=mean(log_count),
            se_log_count=sd(log_count)/sqrt(n())) %>%
  ungroup() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=mean_log_count,color=Phage)) +
  labs(y="Average Log Count", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=mean_log_count-2*se_log_count,
                    ymax=mean_log_count+2*se_log_count),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(Strain~temperature)


################################pH model#####################################

lmer(log_count ~ ph*Phage + day + Strain + Phage*Strain +
       milk_age +
       log_milk_apc +
       (1|milk_batch/cheese_make/cheese_plate) +
       (1|rep),
     data=phage_ph) -> m_ph

phage_ph$ph <- factor(phage_ph$ph, levels = c(6.5, 6, 5.5))

#Raw data used in model above
phage_ph %>%
  select(Phage,ph,Strain,day,experiment,log_count) %>%
  spread(key = Phage,value=log_count) %>%
  mutate(phage_effect=N-Y) %>%
  filter(!is.na(phage_effect)) %>%
  group_by(Strain,ph,day) %>%
  summarize(mean_phage_effect=mean(phage_effect),
            se_phage_effect=sd(phage_effect)/sqrt(n())) %>%
  ggplot(aes(x=day,y=mean_phage_effect)) +
  geom_hline(aes(yintercept=0),color="red",alpha=0.5) +
  geom_point(position=position_dodge(width=0.5)) +
  labs(y="Phage Effect (Log)", x="Day") +
  geom_errorbar(aes(ymin=mean_phage_effect-2*se_phage_effect,
                    ymax=mean_phage_effect+2*se_phage_effect),
                width=.1, position=position_dodge(width=0.5)) +
  ylim(c(-3,5.5)) +
  facet_grid(Strain~ph)

# Full plot for pH model
lsm_ph <- summary(lsmeans(m_ph,~Phage+ph+day))
lsm_ph %>%
  data.frame() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=lsmean,color=Phage)) +
  labs(y="Least Square Means", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(~ph)

#Raw data from plot above
phage_ph %>%
  filter(!is.na(log_count)) %>%
  group_by(day,ph,Phage) %>%
  summarize(mean_log_count=mean(log_count),
            se_log_count=sd(log_count)/sqrt(n())) %>%
  ungroup() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=mean_log_count,color=Phage)) +
  labs(y="Average Log Counts", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=mean_log_count-2*se_log_count,
                    ymax=mean_log_count+2*se_log_count),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(~ph)

# Full pH model (same as above) separated by strain
lsm_ph_w_strain <- summary(lsmeans(m_ph,~Phage+ph+day+Strain))
lsm_ph_w_strain %>%
  data.frame() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=lsmean,color=Phage)) +
  labs(y="Least Square Means", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=lower.CL,ymax=upper.CL),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(Strain~ph)

# Raw data of above plot

phage_ph %>%
  filter(!is.na(log_count)) %>%
  group_by(day,ph,Phage,Strain) %>%
  summarize(mean_log_count=mean(log_count),
            se_log_count=sd(log_count)/sqrt(n())) %>%
  ungroup() %>%
  mutate(day=as.numeric(as.character(day))) %>%
  ggplot(aes(x=day,y=mean_log_count,color=Phage)) +
  labs(y="Average Log Count", x="Day", color = "Phage Added") +
  scale_color_manual(labels = c("No", "Yes"), values = c("red", "blue")) +
  geom_point(position=position_dodge(1)) +
  geom_errorbar(width=1,
                aes(ymin=mean_log_count-2*se_log_count,
                    ymax=mean_log_count+2*se_log_count),
                position=position_dodge(1)) +
  geom_line(position=position_dodge(1),aes(group=Phage)) +
  scale_x_continuous(breaks=c(1,7,14)) +
  ylim(3,11) +
  facet_grid(Strain~ph)