#Run this after the data script
#14/12/2022

# Load libraries ----------------------------------------------------------
library(tidyverse) # data wrangling tools
library(stringr) #working with strings
library(zoo) # roll functions
library(here) #directory management
library(cowplot) #composite plots for ggplot
library(openxlsx) # read in excel files xlsx
library(yardstick) #  MPE function
library(rstatix) # statistical analyses
library(caret) #classification analyses
library(viridis)
library(mgcv)
library(emmeans)


# Read in data files ------------------------------------------------------
#Read in Sr isotope data
strontium_readin <- read.xlsx(here('data','sutter_bypass_masterfile_juvenile_20230323.xlsx'), sheet="Sr8786")

#Read in Sr isoscape
water_isoscape_readin <- read.csv(here::here("data", "water_sr_isoscape_forBC.csv"))

#Growth data
growth_cagewildcomb_summary <- read.csv(here::here("output", "growth_cagewildcomb.csv"))


# Sr data: Clean and merge data ----------------------------------------------------
strontium_join <- strontium_readin%>%
  group_by(Sample_ID)%>%
  mutate(spot_start=lag(Distance_um),
         spot_stop= Distance_um)%>%
  ungroup()

data_merged <-oto_growth_lifetime %>%
  left_join (strontium_join, by="Sample_ID")%>%
  mutate(Sr87Sr86_assigned=case_when (Inc_distance >spot_start & 
                                        Inc_distance <= spot_stop ~Sr8786_norm)) %>%
  drop_na(Sr87Sr86_assigned)%>%
  filter(Sr87Sr86_assigned>0)

#Create the dataframe for further analysis and plotting
data_clean <- data_merged %>%
  dplyr::select (Sample_ID,Fish_ID, backdate,Tag=Tag.x,Inc_no, Inc_distance, Location, Type, Region,
                 FL_mm, Wt, growth_bimf, growth_bimf_lifetime,Spot_no, 
                 Distance_um, Sr87Sr86_assigned,SE2, Inc_no_max, Inc_distance_max)%>%
  mutate(Sr_habitat=case_when(
                              Sr87Sr86_assigned>= 0.70420 & Sr87Sr86_assigned< 0.70561~"Butte Creek watershed",
                              Sr87Sr86_assigned>= 0.70561 & Sr87Sr86_assigned< 0.70620~"Sacramento River",
                              TRUE ~"Other"))%>%
  mutate(Sr_habitat=factor(Sr_habitat, levels=c("Sacramento River","Butte Creek watershed", "Other")))%>%
  group_by(Fish_ID)%>%
  mutate (movement_class=case_when(max(Sr87Sr86_assigned)>=0.70561~"Recent Migrant",
          TRUE~"Butte Creek Resident"))%>%
  mutate (transition_incs=case_when(Sr_habitat=="Sacramento River" ~ Inc_no))%>%
  mutate (transition=max(transition_incs, na.rm=T))%>%
  ungroup()

#Clean isoscape
sr_isoscape <- water_isoscape_readin %>%
  dplyr::select(River, Sr8786,SD, MEAN,MIN,MAX)


# Strontium plots ---------------------------------------------------------

#Plot 87Sr86Sr spot data from the chemistry file
p_sr87sr86_spots<- ggplot(data=strontium_join%>% filter (!is.na(Sample_ID))) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.705610, ymax=0.70620, 
                      alpha = 0.4, fill = 'red')+ #Sac River
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.70420, ymax =0.705610, 
           alpha = 0.4, fill = 'orange')+ #Butte Creek
  annotate('text', x = 100,  y =0.7057100,alpha = 0.4, label = 'Sac River')+ #Sac River
  annotate('text', x = 100, y= 0.7047566,alpha = 0.4, label = 'Butte Creek watershed')+ #Butte Creek
  geom_line(aes(x=Distance_um,y=Sr8786_norm, group=Sample_ID))+
  geom_pointrange(aes(x=Distance_um,y=Sr8786_norm, group=Sample_ID, 
                      ymin=Sr8786_norm-SE2, ymax=Sr8786_norm+SE2))+
  scale_y_continuous(name=expression(""^"87"*"Sr/"^"86"*"Sr"))+
  scale_x_continuous(name=expression("Distance ["*mu*"m]"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(.~Sample_ID, ncol=3, scales="free_x")
p_sr87sr86_spots
ggsave(plot=p_sr87sr86_spots, "output/sr87sr86_spots.png", width = 12, height = 18)


# Plot 87Sr86Sr by increment no
p_sr87sr86<- ggplot(data=data_clean) +
  geom_line(aes(x= Inc_no,y=Sr87Sr86_assigned, group=Sample_ID))+
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.705610, ymax=0.70620, 
           alpha = 0.4, fill = 'red')+ #Sac River
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.70420, ymax =0.705610, 
           alpha = 0.4, fill = 'orange')+ #Butte Creek
  annotate('text', x = 100,  y =0.7057100,alpha = 0.4, label = 'Sac River')+ #Sac River
  annotate('text', x = 100, y= 0.7047566,alpha = 0.4, label = 'Butte Creek watershed')+ #Butte Creek
  # geom_smooth(aes(x= Inc_no,y=Sr87Sr86_assigned, group=Sample_ID))+
  scale_y_continuous(name=expression(""^"87"*"Sr/"^"86"*"Sr"))+
  scale_x_continuous(name=expression("Age"), expand=c(0,0))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(.~Sample_ID, ncol=5, scales="free_x")
p_sr87sr86
ggsave(plot=p_sr87sr86, "output/sr87sr86.png", width = 18, height = 20)

# Plot 87Sr86Sr by increment no and summary
p_sr87sr86_summary<- ggplot(data=data_clean) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.705610, ymax=0.70620, 
           alpha = 0.4, fill = 'grey50')+ #Sac River
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.70420, ymax =0.705610, 
           alpha = 0.4, fill = 'grey85')+ #Butte Creek
  annotate('text', x = 100,  y =0.7057100,alpha = 0.4, label = 'Sac River')+ #Sac River
  annotate('text', x = 100, y= 0.7047566,alpha = 0.4, label = 'Butte Creek watershed')+ #Butte Creek
  # geom_smooth(aes(x= Inc_no,y=Sr87Sr86_assigned, group=Sample_ID))+
  geom_smooth(aes(x= Inc_no,y=Sr87Sr86_assigned, group=Sample_ID, color=Type), se=F)+
  scale_y_continuous(name=expression(""^"87"*"Sr/"^"86"*"Sr"))+
  scale_x_continuous(name=expression("Age"), expand=c(0,0))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_sr87sr86_summary
ggsave(plot=p_sr87sr86_summary, "output/sr87sr86_summary.png", width = 5, height = 5)


#Lifetime growth reconstructions and place reconstructions for wild fish 
p_growth_place <-ggplot (data=data_clean%>%filter(Type=="Wild"))+
  geom_tile(aes(x=backdate, y=Fish_ID, fill=growth_bimf, group=Fish_ID))+
  scale_x_date(name="Date")+
  scale_y_discrete(name="Fish ID")+
  scale_fill_viridis(name="Growth rate (mm/d)")+
  theme_bw()+
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6))+
  facet_wrap(~Sr_habitat)
p_growth_place
ggsave(plot=p_growth_place, "output/growth_place.png", width = 8, height = 5)


p_growth_box <-ggplot (data=data_clean%>%filter(Type=="Wild"))+
  geom_boxplot(aes(x=movement_class, y=growth_bimf_lifetime, fill=movement_class), outlier.shape=NA, width=0.5)+
  geom_point(aes(x=movement_class, y=growth_bimf_lifetime),
             fill="grey", alpha=0.5, color="black", shape=21, size=2)+
  scale_x_discrete(name="Movement history")+
  scale_y_continuous(name="Growth rate (mm/d)")+
  scale_fill_manual(values=c("steelblue","darkgreen"))+
  theme_bw()+
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6))
p_growth_box
ggsave(plot=p_growth_box, "output/growth_place_box.png", width = 5, height = 5)

p_growth_box_gray <-ggplot (data=data_clean%>%filter(Type=="Wild"))+
  geom_boxplot(aes(x=movement_class, y=growth_bimf_lifetime, fill=movement_class), outlier.shape=NA, width=0.5)+
  geom_point(aes(x=movement_class, y=growth_bimf_lifetime),
             fill="grey", alpha=0.5, color="black", shape=21, size=2)+
  scale_x_discrete(name="Movement history")+
  scale_y_continuous(name="Growth rate (mm/d)")+
  scale_fill_manual(values=c("#D3D3D3","#D3D3D3"))+
  theme_bw()+
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8))
p_growth_box_gray
ggsave(plot=p_growth_box_gray, "output/growth_place_box_gray.png", width = 4, height = 5)

#FL by age
growth_cagewildcomb_summary_mov <- growth_cagewildcomb_summary %>%
  left_join(dplyr::select(data_clean,Fish_ID, movement_class), by="Fish_ID")%>%
  filter(Type=="Wild")%>%
  mutate(movement_class=case_when(is.na(movement_class)~"Unknown",
                                  TRUE~movement_class))%>%
  distinct(Fish_ID, Inc_no_max, Inc_distance_max, 
           FL_mm, Growth_cat, Growth_rate_mean, capture_location, movement_class)


p_FL_age_org <-ggplot ()+
  geom_smooth(data=growth_cagewildcomb_summary_mov,
              aes(x=Inc_no_max, y=FL_mm, color=capture_location), 
              method="lm", se=F, linetype="dashed")+
  geom_point(data=growth_cagewildcomb_summary_mov%>%filter(!movement_class=="Unknown"),
             aes(x=Inc_no_max, y=FL_mm, shape=movement_class), size=2)+
  # geom_text(data=growth_cagewildcomb_summary_mov%>%filter(!movement_class=="Unknown"),
  #            aes(x=Inc_no_max+5, y=FL_mm, label=Fish_ID), size=2)+
  scale_x_continuous(name="Age (days)", limits=c(0,100))+
  scale_y_continuous(name="Fork length (mm)")+
  # scale_fill_manual(values=c("steelblue","red"))+
  # scale_color_manual(values=c("steelblue","red"))+
  # scale_shape_manual(values=c(22,21,23))+
  theme_bw()+
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6))
p_FL_age_org
ggsave(plot=p_FL_age_org, "output/FL_age_org.png", width = 6, height = 4)



#Plot 87Sr86Sr spot data from the chemistry file for just wild fish
p_sr87sr86_spots_wild<- ggplot(data=data_clean%>%filter(Type=="Wild")) +
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.705610, ymax=0.70620, 
           alpha = 0.4, fill = 'steelblue')+ #Sac River
  annotate('ribbon', x = c(-Inf, Inf), ymin = 0.70420, ymax =0.705610, 
           alpha = 0.4, fill = 'darkblue')+ #Butte Creek
  annotate('text', x = 100,  y =0.7058100,alpha = 0.4, label = 'Lower Sacramento River')+ #Sac River
  annotate('text', x = 100, y= 0.7047566,alpha = 0.4, label = 'Butte Creek watershed')+ #Butte Creek
  geom_pointrange(data=strontium_join%>% filter(grepl("SBP", Sample_ID))%>%filter(Sr8786_norm>0), aes(x=Distance_um,y=Sr8786_norm, group=Sample_ID, 
                      ymin=Sr8786_norm-SE2, ymax=Sr8786_norm+SE2), 
                  color="grey")+
  # geom_line(aes(x=Distance_um,y=Sr87Sr86_assigned, group=Sample_ID))+
  geom_pointrange(aes(x=Distance_um,y=Sr87Sr86_assigned, group=Sample_ID, 
                      ymin=Sr87Sr86_assigned-SE2, ymax=Sr87Sr86_assigned+SE2),
                  color="black")+
  scale_y_continuous(name=expression(""^"87"*"Sr/"^"86"*"Sr"))+
  scale_x_continuous(name=expression("Distance ["*mu*"m]"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(.~Sample_ID, ncol=4, scales="free_x")
p_sr87sr86_spots_wild
ggsave(plot=p_sr87sr86_spots_wild, "output/sr87sr86_spots_wild.png", width = 10, height = 7)


data_clean_check <- data_clean %>%
  distinct(Sample_ID, .keep_all=T)



#Lifetime growth reconstructions and place reconstructions for wild fish linefacts
p_growth_place_lines <-ggplot (data=data_clean%>%filter(Type=="Wild"))+
  geom_line(aes(x=Inc_no, y=growth_bimf, group=Fish_ID, color=movement_class), alpha=0.5)+
  geom_smooth(aes(x=Inc_no, y=growth_bimf, color=movement_class), span=0.8)+
  scale_x_continuous(name="Increment number")+
  scale_y_continuous(name="Growth rate (mm/d)")+
  # geom_vline(aes(xintercept=transition), linetype="dashed")+
  # annotate('ribbon',xmin=5, xmax=20, y = c(-Inf, Inf),
  #          alpha = 0.4, fill = 'steelblue')+ #Sac River
  # annotate('ribbon', xmin=25, xmax=60, y = c(-Inf, Inf),
  #          alpha = 0.4, fill = 'darkblue')+ #Butte Creek
  # annotate('text', x = 10,  y =0.5, label = 'Lower Sacramento River')+ #Sac River
  # annotate('text', x = 3, y= 0.5, label = 'Butte Creek watershed')+ #Butte Creek
  theme_bw()+
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6))
  # facet_wrap(~Fish_ID)
p_growth_place_lines 
ggsave(plot=p_growth_place_lines , "output/growth_place_lines .png", width = 8, height = 5)


# Growth comparison among habitats wild fish ------------------------------
growthcomb_m <- data_clean%>%filter(Type=="Wild")%>%
  drop_na(Inc_no,growth_bimf)%>%
  mutate(Fish_ID=as.factor(Fish_ID))%>%
  mutate(movement_class=factor(movement_class))

#Fit gam model
m1 <- gam(growth_bimf ~ 
            s(Inc_no, k=10) +
            s(Fish_ID, bs = "re") + 
            movement_class,
          data = growthcomb_m, discrete = FALSE, method = "ML", select = TRUE)

#Explore results
summary(m1)
plot(m1, all.terms=T, pages=1, rug=T)

#Compare movement class
mc_comb <- emmeans(m1, pairwise~movement_class, data=growthcomb_m)
mc_comb
mc_comb_plot <-plot(mc_comb, 
     ylab="Movement history", 
     xlab="Estimated marginal mean growth rates (mm/d)",
     color="grey8")+
  theme_classic()
mc_comb_plot
ggsave(plot=mc_comb_plot, "output/mc_comb_plot.png", width = 6, height = 5)
