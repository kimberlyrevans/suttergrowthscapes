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
library(MASS)# classification analyses
library(viridis)#Color scale
library(ggforce) #Options for ggplot
library(ggpubr)#Options for ggplot
library(lme4) # Mixed models
library(modelsummary) # Mixed models evaluation
library(lmerTest) # Mixed models evaluation
library(RFishBC) # Different otolith models
library(FSA) #Fish stock analysis
library(mgcv)
library(emmeans)


#Clear the workspace
rm(list = ls())

# Load data ---------------------------------------------------------------
#Load otolith increment data from Sutter
increments_readin <- read.xlsx(here::here("data", "Inc_count_sutter.xlsx"),
                         sheet = "Inc_Counts")

#Load and clean otolith and fish size data for building the models
exogint_readin <- read.xlsx(here::here("data", "exog_intercept_20210430.xlsx"),
                               sheet = "sample_data")

## Load and clean Sutter fish data
fish_readin <- read.csv(here::here("data", "fish.csv"))
fish_data <- fish_readin%>%
  mutate(Date=as.Date(Date))%>%
  mutate(Year=lubridate::year(Date))%>%
  mutate(Fish_ID= case_when(Experiment=="Wild" ~ ID,
                          Experiment=="Cage"~str_sub(ID, 1, -10)))


## Load and clean Sutter cage locations
cage_location_readin <- read.csv(here::here("data", "Sutter_cage_locations.csv"))
cage_location <- cage_location_readin

#Load and clean Sutter genetic data
genetics_readin <- read.csv(here::here("data", "Sutter_wild_fish_with_genetics.csv"))
genetics <- genetics_readin

# Clean increment data ----------------------------------------------------
#Remove increment reads that are not accurate (QAQC: vaterite or poor increment expression)
fish_to_remove <- c("4283_20190402", 
                    "4466_20190401", 
                    "4362_20190403", 
                    "4348_20190401",
                    "4028_20190401",
                    "SBP19-149",
                    "SBSNWR3939_20200317",#Lots of vaterite on edge
                    "SBP22A010", #Remove the 1 sample from 2022,
                    "SBP20-002"
                    ) 

increments <- increments_readin %>%
  #Classify samples as either wild or cage based on their sample ID
  mutate (Experiment=ifelse(substr(Sample_ID, 1,3) == "SBP", "Wild", "Cage"))%>%
  #Filter data to remove the SBWS fish
  filter(!substr(Sample_ID,1,4)=="SBWS") %>%
  #Filter data to only keep the first analysis
  filter(Analysis_no==1)%>%
  #Assign an increment number to the edge
  group_by(Sample_ID)%>%
  mutate(Inc_no=case_when(Tag=="EDGE" ~ max(Inc_no, na.rm=T)+1,
                          Tag=="EXOG" ~ 0,
         TRUE~Inc_no))%>%
  ungroup()%>%
  filter(!Sample_ID %in% fish_to_remove) %>%
  left_join(dplyr::select(genetics, Sample_ID, Genetics_specific), by="Sample_ID")%>%
  dplyr::select(-Row_no)


#How many fish with inc counts?
n_distinct(increments$Sample_ID)


# QAQC Increment data -----------------------------------------------------
#Histogram of increment counts
ggplot()+geom_histogram(aes(increments$Inc_no))
#Plot forklengths
ggplot()+geom_histogram(aes(increments$FL_mm))
#Are the inbetween distance increments reasonable?
ggplot()+geom_histogram(aes(increments$Between_inc_dist))
#Any correlation between inc number and between inc dist
plot(increments$Inc_no, increments$Between_inc_dist)

# Plot increment profiles for each fish
p_inc<- ggplot(data=increments) +
  geom_line(aes(x=Inc_no,y=Between_inc_dist, group=Sample_ID, color=Experiment))+
  scale_y_continuous(name=expression("Between increment distances ["*mu*"m]"))+
  scale_x_continuous(name="Increment number")+
  theme_bw() +
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  facet_wrap(~Sample_ID)
p_inc
ggsave(plot=p_inc, "output/Between_inc_distances.png", width = 10, height = 7)

#Plot increment confidence profiles for each fish
p_inc_conf<- ggplot(data=increments) +
  geom_line(aes(x= Inc_no,y=Inc_confidence, group=Sample_ID), color="black")+
  theme_bw() +
  theme(legend.position="right")+
  facet_wrap(~Sample_ID)
p_inc_conf
ggsave(plot=p_inc_conf, "output/inc_conf.png", width = 10, height = 7)

#Plot FL vs oto size as a QAQC
p_oto_fl_check <- ggplot(increments%>%filter(Tag=="EDGE"), aes(x = Inc_distance, y = FL_mm)) + 
  geom_point() +
  geom_smooth(method = "lm", col = "black",formula=y~x-1)+
  geom_text(aes(label=Sample_ID))+
  # coord_cartesian(xlim=c(0, NA), ylim=c(0,NA))+
  scale_x_continuous(name=expression("Increment distance ["*mu*"m]"))+
  scale_y_continuous(name="FL measured (mm)")+
  theme_bw() +
  theme(legend.position="right", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p_oto_fl_check
ggsave(plot=p_oto_fl_check, "output/oto_fl_check.png", width = 5, height = 5)

# Parameterize the Fraslee and Bioint models ------------------------------
#Find intercepts for Fraslee and Bioint model using our caged fish and larger dataset

#Prepare otolith data for sutter cages
oto_model_data <-increments %>%
  filter(Tag=="EDGE")%>%
  filter(Experiment=="Cage")

exogint<- exogint_readin %>%
  drop_na(exog_dist, edge_dist, fork_length)%>%
  #Restrict data
  filter(fork_length<150)%>%
  filter(edge_dist< 600)%>%
  #Remove flagged fish
  filter(!lab_id%in% c("TT19363",
                       "DM15024",
                       "DM16003",
                       "DM16366",
                       "TT19389",
                       "TT19325",
                       "TT19312",
                       "DM16038",
                       "DM15026",
                       "DM17145",
                       "DM16271",
                       "DM16339",
                       "DM16148",
                       "DM17281",
                       "TT10377",
                       "DM17562",
                       "DM16435",
                       "DM16436",
                       "TT19377"
  ))

#Model for the exog int dataset
model_fraslee_lm <- lm(fork_length~edge_dist, data=exogint)
summary(model_fraslee_lm)

#Model just for the sutter dataset
model_fraslee_lm_sutter <- lm(FL_mm~Inc_distance, data=oto_model_data)
summary(model_fraslee_lm_sutter)

# plot(model_fraslee_lm)
p_fraslee_intercept <- ggplot() + 
  geom_point(data=exogint, aes(x = edge_dist, y = fork_length), color="gray46", alpha=0.65) +
  geom_point(data=oto_model_data, aes(x = Inc_distance, y = FL_mm), color="orange") +
  geom_smooth(data=exogint, aes(x = edge_dist, y = fork_length),method = "lm", col = "black",formula=y~x-1)+
  # geom_text(data=exogint, aes(x = edge_dist, y = fork_length, label=lab_id))+
  # coord_cartesian(xlim=c(0, NA), ylim=c(0,NA))+
  scale_x_continuous(name=expression("Otolith radius ["*mu*"m]"))+
  scale_y_continuous(name="Fork Length (cm)")+
  theme_bw() +
  theme(legend.position="right", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
p_fraslee_intercept
ggsave(plot=p_fraslee_intercept, "output/fraslee_intercept.png", width = 5, height = 5)

#from the lm of the entire dataset
Fraslee_int =-3.14

# #from average of Titus et al., 2004
# Biont_flsize=33.4
# Biont_otosize=220

#Old G. database version 
# Biont_flsize=29
# Biont_otosize=220

#New G. database version
# # Fraslee_int =-2.8335
Biont_flsize=33.5
Biont_otosize=220
#MF
#a is fish size at otolith formation
Biont_a=18

# Prepare dataframe for analyses ---------------------------------------------------

#Merge increments and sutter fish data and apply otolith models
otofish_prep <- increments %>%
  left_join(fish_data, by=c("Sample_ID"="ID", "Experiment"))%>%
  dplyr::select(-FL, -Inc_confidence, -Analysis_no,-Reader, -Inc_counts_Notes)%>%
  #Calculate date at otolith increments
  group_by(Sample_ID)%>%
  mutate(Date=as.Date(Date, format = "%y/%m/%d/"))%>%
  mutate(backinc=max(Inc_no, na.rm=T)-Inc_no,
         backdate=as.Date(as.Date(Date, format = "%y/%m/%d/") - backinc))%>%
  mutate(Edge=max(Inc_distance, na.rm=T))%>%
  ungroup()%>%
  #Join cage location info
  left_join(dplyr::select(cage_location,-Type,-Region), by=c("Location"="Site"))%>%
  ungroup()%>%
  mutate(Fish_ID= case_when(Experiment=="Wild" ~ Sample_ID,
                          Experiment=="Cage"~str_sub(Sample_ID, 1, -10)))%>%
  left_join(dplyr::select(fish_data, Fish_ID, Date, FL_interim=FL), by=c("Fish_ID", "backdate"="Date"))

#Prep data for backcalculation
backcomb<-otofish_prep %>%
  dplyr::select(id=Fish_ID, ann=Inc_no, rad=Inc_distance, lencap=FL_mm)%>%
  group_by(id)%>%
  filter(!is.na(ann))%>%
  mutate(agecap=max(ann),
         radcap=max(rad))%>%
  ungroup()

#Backcalc equations
#Fraslee
backcomb_fraslee <- backCalc(backcomb,lencap,BCM="FRALE",
                             inFormat="long",outFormat="long",a=Fraslee_int, digits=2)%>%
  mutate(Method="fraslee")

#Bioint
backcomb_bioint <- backCalc(backcomb,lencap,BCM="BI",L0p = Biont_flsize, R0p=Biont_otosize,
                            inFormat="long",outFormat="long",digits=2)%>%
  mutate(Method="bioint")

#MF
backcomb_MF <- backCalc(backcomb,lencap,BCM="MF",L0p = Biont_flsize, R0p=Biont_otosize,a=Biont_a,
                        inFormat="long",outFormat="long",digits=2)%>%
  mutate(Method="MF")

#Join dataframes back together
backcomb_results <- rbind(backcomb_fraslee, backcomb_bioint,backcomb_MF)%>%
  pivot_wider(names_from="Method", values_from="bclen")

#Final dataframe with otolith models applied

otofish <-otofish_prep %>%
  left_join(backcomb_results, by=c("Fish_ID"="id", "Inc_no"="ann"))%>%
  mutate(FL_fraslee=fraslee, FL_bim=bioint, FL_bimf=MF)%>%
  dplyr::select(-lencap, -agecap, -fraslee, -bioint, -MF)%>%
  group_by(Sample_ID)%>%
  mutate(Inc_no_max=max(Inc_no, na.rm=T))%>%
  mutate(Inc_distance_max=max(Inc_distance, na.rm=T))%>%
  ungroup()%>%
  #Calculate wateryear
  mutate(wateryear=case_when (lubridate::month(Date) %in% c(10,11,12) ~ Year +1,
                              TRUE ~Year))

#Export dataframe
write.csv(otofish, "output/otofish.csv", row.names=F)

#Sample overview table
sample_overview <- otofish%>%
  group_by(Sample_ID)%>%
  distinct(Sample_ID, Type, Region,Year, Date, wateryear, FL_mm, Wt, 
                Inc_no_max,Inc_distance_max)%>%
  arrange(Type, Region, Year)

write.csv (sample_overview, "output/sample_overview.csv", row.names=F)


# Objective 1: Evaluate the otolith models --------------------------------

#Just using caged fish
p_fl_recon_check<- ggplot(data=otofish%>% filter (Experiment=="Cage")) +
  geom_line(aes(x= backdate,y=FL_fraslee, group=Sample_ID), color="purple")+
  geom_line(aes(x= backdate,y=FL_bim, group=Sample_ID), color="orange")+
  geom_line(aes(x= backdate,y=FL_bimf, group=Sample_ID), color="darkgreen")+
  geom_point(aes(x= backdate,y=FL_interim, group=Sample_ID), color="black")+
  scale_y_continuous(name="Fork length")+
  scale_x_date(name="Date")+
  # coord_cartesian(xlim=as.Date(c('2019-02-10', '2020-04-01'), format="%Y-%m-%d"), ylim=c(35,75))+
  theme_bw() +
  theme(legend.position="right")+
  facet_wrap(~Fish_ID, scales="free_x", ncol=4)
p_fl_recon_check
ggsave(plot=p_fl_recon_check, "output/fl_recon_check.png", width = 8, height =12)

samples_to_exclude_2020 <- c("BCGR3599_20200316", "BCGR4191_20200316", "BCGR4479_20200316", "FRLF3543_20200316", "RGSAC4221_20200317", "SACWL4192_20200318", "SBSNWR3553_20200317", "SBSNWR3619_20200317")
p_fl_recon_check_2019 <- ggplot(data = otofish %>% filter(Experiment == "Cage" & !Sample_ID %in% samples_to_exclude_2020)) +
  geom_line(aes(x = backdate, y = FL_fraslee, group = Sample_ID), color = "purple") +
  geom_line(aes(x = backdate, y = FL_bim, group = Sample_ID), color = "orange") +
  geom_line(aes(x = backdate, y = FL_bimf, group = Sample_ID), color = "darkgreen") +
  geom_point(aes(x = backdate, y = FL_interim, group = Sample_ID), color = "black") +
  scale_y_continuous(name = "Fork length") +
  scale_x_date(name = "Date", date_labels = "%b") +
  coord_cartesian(xlim = as.Date(c('2019-02-10', '2019-04-01'), format = "%Y-%m-%d"), ylim = c(35, 75)) +
  theme_bw() +
  theme(legend.position = "right") +
  facet_wrap(~ Fish_ID, scales = "free_x", ncol = 5)

p_fl_recon_check_2019
ggsave(plot=p_fl_recon_check_2019, "output/fl_recon_check_2019.png", width = 8, height =12)

## Work In Progress
## comment: "depict the 'goodness' of fits in a single plot"
p_fl_recon_check_comb_2019 <- ggplot(data = otofish %>% filter(Experiment == "Cage" & !Sample_ID %in% samples_to_exclude_2020)) +
  geom_line(aes(x = backdate, y = FL_fraslee, group = Sample_ID), color = "#D3D3D3") +
  geom_line(aes(x = backdate, y = FL_bim, group = Sample_ID), color = "#D3D3D3") +
  geom_line(aes(x = backdate, y = FL_bimf, group = Sample_ID), color = "#D3D3D3") +
  geom_point(aes(x = backdate, y = FL_interim, group = Sample_ID), color = "black") +
  scale_y_continuous(name = "Fork length") +
  scale_x_date(name = "Date", date_labels = "%b") +
  coord_cartesian(xlim = as.Date(c('2019-02-10', '2019-04-01'), format = "%Y-%m-%d"), ylim = c(35, 75)) +
  theme_bw() +
  theme(legend.position = "right")+
  geom_smooth(aes(x=backdate, y=FL_bimf, color= "FL_bimf"), span=0.8)+
  geom_smooth(aes(x=backdate, y=FL_bim, color= "FL_bim"), span=0.8)+
  geom_smooth(aes(x=backdate, y=FL_fraslee, color= "fraslee"), span=0.8)+
  geom_smooth(aes(x=backdate, y=FL_interim, color= "FL_interim"), span=0.8)

p_fl_recon_check_comb_2019

oto_flcomp <- otofish %>%
  distinct(Fish_ID,FL_interim, backinc, .keep_all=T) %>%
  filter (Experiment=="Cage")%>%
  drop_na(FL_interim)%>%
  #Determine dates between field measurements
  group_by(Fish_ID)%>%
  arrange(Fish_ID,desc(backdate))%>%
  mutate(growth_days=lag(backdate, n=1L)-backdate)%>%
  mutate(growth_days=as.numeric(growth_days))%>%
  #Determine FL growth between field measurements
  mutate(FL_growth_field=lag(FL_interim, n=1L)-FL_interim)%>%
  #Determine FL growth fraslee
  mutate(FL_growth_fraslee=lag(FL_fraslee, n=1L)-FL_fraslee)%>%
  #Determine FL growth bim
  mutate(FL_growth_bim=lag(FL_bim, n=1L)-FL_bim)%>%
  #Determine FL growth mf bim
  mutate(FL_growth_bimf=lag(FL_bimf, n=1L)-FL_bimf)%>%
  ungroup()%>%
  #Remove the last value
  filter(backinc>0)%>%
  #Determine Growthrates
  group_by(Fish_ID)%>%
  #field observed growth
  mutate(observed_growth=FL_growth_field/growth_days)%>%
  mutate(observed_growth_mean=(max(FL_interim)-min(FL_interim))/as.numeric(max(backdate)-min(backdate)))%>%
  #fraslee growth rates
  mutate(fraslee_growth=FL_growth_fraslee/growth_days)%>%
  mutate(fraslee_growth_mean=(max(FL_fraslee)-min(FL_fraslee))/as.numeric(max(backdate)-min(backdate)))%>%
  #bim growth rates
  mutate(bim_growth=FL_growth_bim/growth_days)%>%
  mutate(bim_growth_mean=(max(FL_bim)-min(FL_bim))/as.numeric(max(backdate)-min(backdate)))%>%
  #mf bim growth rates
  mutate(bimf_growth=FL_growth_bimf/growth_days)%>%
  mutate(bimf_growth_mean=(max(FL_bimf)-min(FL_bimf))/as.numeric(max(backdate)-min(backdate)))%>%
  ungroup()%>%
  #Define growth categories
  mutate(Growth_cat=case_when (observed_growth_mean >0.5 ~ "High",
                               observed_growth_mean <=0.5 & observed_growth_mean >0.2 ~"Medium",
                               observed_growth_mean <=0.2 ~ "Low",
                               Type=="Wild"~"Wild"))%>%
  mutate(Growth_cat=factor(Growth_cat, levels=c("Low","Medium","High", "Wild")))

#Export otoflcomp dataframe 
write.csv(oto_flcomp, "output/oto_flcomp.csv", row.names=F)

#Mean growth rates by Method for caged fish
oto_flcomp_mean <-oto_flcomp %>%
 distinct(Fish_ID, Year, Type, Report_ID,Growth_cat,FL_mm, Inc_no_max,Inc_distance_max, observed_growth_mean, fraslee_growth_mean, bim_growth_mean, bimf_growth_mean)%>%
  pivot_longer(names_to="Method", values_to="Growth_rate_mean",9:12) %>%
  ungroup()%>%
  mutate(Method=factor(Method, levels=c ("fraslee_growth_mean",
                                         "bim_growth_mean", 
                                         "bimf_growth_mean",
                                         "observed_growth_mean")))

# code fix for line 380
oto_flcomp_mean <- oto_flcomp %>%
  distinct(Fish_ID, Year, Type, Report_ID, Growth_cat, FL_mm, Inc_no_max, Inc_distance_max, observed_growth_mean, fraslee_growth_mean, bim_growth_mean, bimf_growth_mean) %>%
  pivot_longer(cols = c(observed_growth_mean, fraslee_growth_mean, bim_growth_mean, bimf_growth_mean), names_to = "Method", values_to = "Growth_rate_mean") %>%
  ungroup() %>%
  mutate(Method = factor(Method, levels = c("fraslee_growth_mean", "bim_growth_mean", "bimf_growth_mean", "observed_growth_mean")))

#Export otoflcomp_mean dataframe 
write.csv(oto_flcomp_mean, "output/oto_flcomp_mean.csv", row.names=F)

#Define the growth categories using only the observed growth
p_growth_cat<- ggplot(data=oto_flcomp_mean%>% filter(Method =="observed_growth_mean")) +
  geom_boxplot(aes(x=reorder(Growth_cat,Growth_rate_mean), y=Growth_rate_mean,  fill=Growth_cat), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"))+
  geom_point(aes(x=reorder(Growth_cat,Growth_rate_mean),y=Growth_rate_mean, fill=Growth_cat), 
             color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_x_discrete(name="Growth category")+
  scale_y_continuous(name="Observed growth (mm/day)")+
  scale_fill_manual (name="Growth Category", values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  scale_color_manual (name="Growth Category",values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_cat
ggsave(plot=p_growth_cat, "output/growth_cat.png", width = 5, height = 4)

p_growth_cat_gray<- ggplot(data=oto_flcomp_mean%>% filter(Method =="observed_growth_mean")) +
  geom_boxplot(aes(x=reorder(Growth_cat,Growth_rate_mean), y=Growth_rate_mean,  fill=Growth_cat), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"))+
  geom_point(aes(x=reorder(Growth_cat,Growth_rate_mean),y=Growth_rate_mean, fill=Growth_cat), 
             color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_x_discrete(name="Growth category")+
  scale_y_continuous(name="Observed growth (mm/day)")+
  scale_fill_manual (name="Growth Category", values=c("#D3D3D3","#D3D3D3","#D3D3D3", "#228B22"))+
  scale_color_manual (name="Growth Category",values=c("#D3D3D3","#D3D3D3","#D3D3D3", "#228B22"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_cat_gray
ggsave(plot=p_growth_cat_gray, "output/growth_cat_gray.png", width = 5, height = 4)


#Boxplot growth rate comparison from cage experiment
p_growth_comp<- ggplot(data=oto_flcomp_mean) +
  geom_boxplot(aes(x=Method, y=Growth_rate_mean,  fill=Method), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5)+
  geom_point(aes(x=Method, y=Growth_rate_mean, fill=Method), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Method", labels=c("Fraslee","BI","MF","Observed"))+
  scale_fill_manual(label=c("Fraslee","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  scale_color_manual(label=c("Fraslee","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_comp
ggsave(plot=p_growth_comp, "output/growth_comp.png", width = 5, height = 5)

p_growth_comp_type<- ggplot(data=oto_flcomp_mean) +
  geom_boxplot(aes(x=Growth_cat, y=Growth_rate_mean, group=interaction(Method,Growth_cat), fill=Method), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5)+
  geom_point(aes(x=Growth_cat, y=Growth_rate_mean, group=interaction(Method,Growth_cat), fill=Method), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Growth potential")+
  scale_fill_manual(label=c("Fraslee","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  scale_color_manual(label=c("Fraslee","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_comp_type
ggsave(plot=p_growth_comp_type, "output/growth_comp_type.png", width = 6, height = 5)

#Combine growth boxplots from cage experiment
p_cagebox_growth <-plot_grid(p_growth_comp+theme(legend.position="none"),p_growth_comp_type+theme(legend.position="none"), labels = c("A", "B"), rel_widths = c(2, 3))
p_cagebox_growth
ggsave("output/cagebox_growth.png", p_cagebox_growth,  width = 8, height = 4)

#Test if there is a significant difference between predicted and observed growth rates
#Using a One-way ANOVA: https://www.datanovia.com/en/lessons/anova-in-r/
#Check assumptions:
#Assumption 1: Independence of the observations= Each subject should belong to only one group (no repeated measures)
#True because we are looking only at mean growth rate over the 6 week experiment 

#Assumption 2: No significant outliers 

oto_flcomp_mean %>% 
  group_by(Method, Growth_cat) %>%
  identify_outliers(Growth_rate_mean)
#We do have extreme outliers

#Assumption 3: Normality
# Build the linear model
model_fornorm  <- lm(Growth_rate_mean ~ Method, data = oto_flcomp_mean)

# Create a QQ plot of residuals
ggqqplot(residuals(model_fornorm))
#Create a QQ plot of residuals by Method
ggqqplot(oto_flcomp_mean, "Growth_rate_mean", facet.by = "Method")

# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model_fornorm))

# Compute Shapiro-Wilk test of normality by Method and Growth_cat
oto_flcomp_mean %>% 
  group_by(Method, Growth_cat) %>%
  shapiro_test(Growth_rate_mean)
#QQ plot looks ok and shapiro-wilk test =  we can't assume ~normality

#Assumption 4: Homogeneity of variances
#residuals versus fits plot
plot(model_fornorm, 1)
#Levene test
oto_flcomp_mean %>% 
  levene_test(Growth_rate_mean ~ as.factor(Method))

#Homogeneity of variances is ~ met. We should use the a one-way ANOVA test
#Compute the analysis of variance by Method
growth_comb_model <- welch_anova_test(Growth_rate_mean ~ Method, data = oto_flcomp_mean)
growth_comb_model
#Not significant but can check group wise comparisons anyways
# Pairwise comparisons
growth_comb_model_pwc <- oto_flcomp_mean %>%  games_howell_test(Growth_rate_mean ~ Method)
growth_comb_model_pwc

#Compute the analysis of variance by Method and Growth_cat two way ANOVA
growth_comb_model_mt <- welch_anova_test(Growth_rate_mean ~ Method+Growth_cat, data = oto_flcomp_mean)
growth_comb_model_mt
# Pairwise comparisons
growth_comb_model_mt_pwc <- oto_flcomp_mean %>% 
  group_by(Growth_cat)%>%
  games_howell_test(Growth_rate_mean ~ Method)
growth_comb_model_mt_pwc

# Pairwise comparisons with cat groups
pwc_plot_data <-growth_comb_model_pwc%>% add_xy_position(x = "Method")
pwc_plot <-ggboxplot(oto_flcomp_mean, x = "Method", y = "Growth_rate_mean",
                     fill='Method', alpha=0.5) +
  stat_pvalue_manual(pwc_plot_data, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(growth_comb_model, detailed = TRUE),
    # caption = get_pwc_label(pwc_plot_data)
  )+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Method", label=c("F-L","BI","MF","Observed"))+
  scale_fill_manual(label=c("F-L","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  scale_color_manual(label=c("F-L","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  theme(legend.position="top", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pwc_plot

 # Pairwise comparisons with cat groups
pwc_plot_data_group <-growth_comb_model_mt_pwc %>% add_xy_position(x = "Growth_cat")
pwc_plot_group <-ggboxplot(oto_flcomp_mean, x = "Growth_cat", y = "Growth_rate_mean",
                     fill='Method', group="Growth_cat", alpha=0.5) +
  stat_pvalue_manual(pwc_plot_data_group, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(growth_comb_model_mt, detailed = TRUE),
    caption = get_pwc_label(pwc_plot_data_group)
  )+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Growth Category")+
  scale_fill_manual(label=c("F-L","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  scale_color_manual(label=c("F-L","BI","MF","Observed"),values=c("purple","orange","darkgreen","grey"))+
  theme(legend.position="top", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pwc_plot_group

#Combine growth boxplots from cage experiment
p_growth_pwc <-plot_grid(pwc_plot+theme(legend.position="none"),
                         pwc_plot_group+theme(legend.position="none"), 
                         labels = c("A", "B"), ncol=1)
p_growth_pwc
ggsave("output/growth_pwc.png", p_growth_pwc,  width = 6, height = 10)

# Evaluate the otolith models using continuous approach --------------------

#For the average growth rate per fish from the cage
oto_con_data <-oto_flcomp_mean%>%
  pivot_wider(names_from="Method", values_from="Growth_rate_mean")%>%
  #Residual Observed-predicted
  mutate(d_BIM=observed_growth_mean-bim_growth_mean,
         d_BIMF=observed_growth_mean-bimf_growth_mean,
         d_Fraslee=observed_growth_mean-fraslee_growth_mean)%>%
  group_by(Fish_ID)%>%
  #Mean percent errors
  mutate(BIM_mpe=mpe_vec(observed_growth_mean, bim_growth_mean),
         BIMF_mpe=mpe_vec(observed_growth_mean, bimf_growth_mean),
         Fraslee_mpe=mpe_vec(observed_growth_mean, fraslee_growth_mean))%>%
  ungroup()%>%
  pivot_longer(names_to="Method", values_to="rate", 10:18)

#Mean abs error
mean_abs_error<-oto_con_data%>% 
  filter(Method %in% c("d_BIM","d_BIMF","d_Fraslee"))%>%
  group_by(Method)%>%
  mutate(abs_error=mean(abs(rate)))%>%
  distinct(Method, abs_error)

#Mean abs error filter to observed growth mean <0.2
mean_abs_error_filtered<-oto_con_data%>% 
  filter(Method %in% c("d_BIM","d_BIMF","d_Fraslee"))%>%
  filter(observed_growth_mean <0.2)%>%
  group_by(Method)%>%
  mutate(abs_error=mean(abs(rate)))%>%
  distinct(Method, abs_error)      
         
p_growth_con<- ggplot(data=oto_con_data%>%filter(Method %in% c("BIM_mpe","BIMF_mpe","Fraslee_mpe"))) +
  geom_hline(yintercept=0)+
  geom_point(aes(x=observed_growth_mean, y=rate, color=Method), alpha=0.5)+
  scale_y_continuous(name="Mean Percent Error")+
  scale_x_continuous(name=" Mean Observed Growth Rate (mm/day)")+
  scale_fill_manual(label=c("BI","MF","Fraslee"),values=c("orange","darkgreen","purple"))+
  scale_color_manual(label=c("BI","MF","Fraslee"),values=c("orange","darkgreen","purple"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_con
ggsave(plot=p_growth_con, "output/growth_con.png", width = 5, height = 5)

p_growth_con_error<- ggplot(data=oto_con_data%>%filter(Method %in% c("d_BIM","d_BIMF","d_Fraslee"))) +
  geom_hline(yintercept=0)+
  geom_point(aes(x=observed_growth_mean, y=rate, color=Method), alpha=0.5)+
  # geom_smooth(aes(x=observed_growth_mean, y=rate, color=Method), alpha=0.5, se=F, method="lm")+
  scale_y_continuous(name="Mean Error")+
  scale_x_continuous(name=" Mean Observed Growth Rate (mm/day)")+
  scale_fill_manual(label=c("BI","MF","Fraslee"),values=c("orange","darkgreen","purple"))+
  scale_color_manual(label=c("BI","MF","Fraslee"),values=c("orange","darkgreen","purple"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_con_error
ggsave(plot=p_growth_con_error, "output/growth_con_error.png", width = 5, height = 5)

lm_dBIM<-lm(rate~observed_growth_mean, data=oto_con_data%>%filter(Method=="d_BIM"))
summary(lm_dBIM)

lm_dBIMF<-lm(rate~observed_growth_mean, data=oto_con_data%>%filter(Method=="d_BIMF"))
summary(lm_dBIMF)

lm_dFraslee<-lm(rate~observed_growth_mean, data=oto_con_data%>%filter(Method=="d_Fraslee"))
summary(lm_dFraslee)

# Mixed models for repeated measures ------------------------------------------------------------
#https://m-clark.github.io/mixed-models-with-R/random_intercepts.html
rep_model_data <-oto_flcomp%>% 
  #Residual Observed-predicted
  mutate(d_BIM=observed_growth-bim_growth,
         d_BIMF=observed_growth-bimf_growth,
         d_Fraslee=observed_growth-fraslee_growth)%>%
  dplyr::select (Fish_ID, Inc_no, backdate, Growth_cat, observed_growth, bim_growth, bimf_growth, fraslee_growth, d_BIM, d_BIMF, d_Fraslee)%>%
  arrange(Fish_ID, backdate)%>%
  group_by(Fish_ID)%>%
  mutate(sample_point=row_number())%>%
  group_by(Fish_ID, sample_point)%>%
  #Mean percent errors
  mutate(BIM_mpe=mpe_vec(observed_growth, bim_growth),
         BIMF_mpe=mpe_vec(observed_growth, bimf_growth),
         Fraslee_mpe=mpe_vec(observed_growth, fraslee_growth))%>%
  ungroup()%>%
  pivot_longer(names_to="Method", values_to="rate", c(6:11,13:15))

#Plot of the growth rates for each time point
rep_model_data_facetplot <- rep_model_data%>%filter(Method%in%c("bim_growth","bimf_growth","fraslee_growth"))%>%
  mutate(Method=case_when(Method=="bim_growth" ~"BI",
                          Method=="bimf_growth" ~"MF",
                          Method=="fraslee_growth" ~"F-L"))%>%
  mutate(Method=factor(Method, levels=c("F-L", "BI","MF")))

p_growth_con_time<- ggplot(data=rep_model_data_facetplot) +
  geom_abline(intercept=0, slope=1, linetype="dashed")+
  geom_point(aes(x=observed_growth, y=rate, color=Method), alpha=0.5)+
  scale_y_continuous(name="Predicted Growth Rate (mm/day)", limits=c(-0.1,1.4))+
  scale_x_continuous(name="Observed Growth Rate (mm/day)", limits=c(-0.1,1.4))+
  scale_fill_manual(label=c("BI","MF","F-L"),values=c("purple","orange","darkgreen"))+
  scale_color_manual(label=c("BI","MF","F-L"),values=c("purple","orange","darkgreen"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  facet_wrap(~Method)
p_growth_con_time
ggsave(plot=p_growth_con_time, "output/growth_con_time.png", width = 6, height = 3)

#Plot of the residuals of the growth rates for each time point
p_growth_con_rep<- ggplot(data=rep_model_data%>%filter(Method%in%c("d_BIM","d_BIMF","d_Fraslee"))) +
  geom_hline(yintercept=0, linetype="dashed")+
  geom_point(aes(x=observed_growth, y=rate, color=Method), alpha=0.5)+
  scale_y_continuous(name=" Mean Observed - Predicted Growth Rate (mm/day)")+
  scale_x_continuous(name=" Mean Observed Growth Rate (mm/day)")+
  scale_fill_manual(label=c("BI","MF","F-L"),values=c("orange","darkgreen","purple"))+
  scale_color_manual(label=c("BI","MF","F-L"),values=c("orange","darkgreen","purple"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_con_rep
ggsave(plot=p_growth_con_rep, "output/growth_con_rep.png", width = 5, height = 5)

#Combine growth boxplots from cage experiment
p_growth_error <-plot_grid(p_growth_con_time+theme(legend.position="none"),
                           p_growth_con_rep+theme(legend.position="right"), 
                         labels = c("A", "B"), ncol=2, rel_widths=c(1,1.2))
p_growth_error
ggsave("output/growth_error.png", p_growth_error,  width = 11, height = 5)

#Assess accuracy and precision of the three different backcalculation models
#https://mspeekenbrink.github.io/sdam-r-companion/linear-mixed-effects-models.html

#prep data
rep_model_data_lmer <- rep_model_data %>%
  pivot_wider(names_from="Method", values_from="rate")

fraslee_lmer_m = lmer(observed_growth ~fraslee_growth +  (1 | Fish_ID),
             data =rep_model_data_lmer)
summary(fraslee_lmer_m)
get_gof(fraslee_lmer_m)

bim_lmer_m = lmer(observed_growth ~bim_growth + (1 | Fish_ID),
                  data =rep_model_data_lmer)
summary(bim_lmer_m)
get_gof(bim_lmer_m)

bimf_lmer_m = lmer(observed_growth ~bimf_growth + (1 | Fish_ID),
                   data =rep_model_data_lmer)
summary(bimf_lmer_m)
get_gof(bimf_lmer_m)
plot(bimf_lmer_m)

#compare the basic models
anova(fraslee_lmer_m, bim_lmer_m, bimf_lmer_m)

#Also create models with age included
fraslee_lmer_m_age= lmer(observed_growth ~fraslee_growth +  Inc_no + (1 | Fish_ID),
                      data =rep_model_data_lmer)
bim_lmer_m_age = lmer(observed_growth ~bim_growth + Inc_no + (1 | Fish_ID),
                  data =rep_model_data_lmer)
bimf_lmer_m_age = lmer(observed_growth ~bimf_growth + Inc_no + (1 | Fish_ID),
                   data =rep_model_data_lmer)
summary(bimf_lmer_m_age)
summary(bim_lmer_m_age)
summary(fraslee_lmer_m_age)

#compare all models 
anova(fraslee_lmer_m, bim_lmer_m, bimf_lmer_m, fraslee_lmer_m_age,bim_lmer_m_age,bimf_lmer_m_age )

#create output table
lmer_model_table <- rbind(get_gof(fraslee_lmer_m),get_gof(bim_lmer_m),get_gof(bimf_lmer_m),
                          get_gof(fraslee_lmer_m_age),get_gof(bim_lmer_m_age),get_gof(bimf_lmer_m_age))%>%
  mutate(Model=c("F-L","BI","MF", "F-L age","BI age","MF age"))
write.csv(lmer_model_table, "output/lmer_model_table.csv", row.names=F)

#How do the residuals change with fish growth rate for our best model

#Intercept only model
#This fits the mean intercept and also the mean intercept of the repeated measure for each fish
g_mm0=lmer(d_BIMF ~ 1 + (1 | Fish_ID),
           data =rep_model_data_lmer)
summary(g_mm0)
#The mean d_rate is the intercept
#Std. Error is how much the fish vary around that

#How does the residual change with observed growth rate, controlling for Fish ID
g_mm1 = lmer(d_BIMF  ~ observed_growth  + (1 | Fish_ID),
             data =rep_model_data_lmer)
summary(g_mm1)
get_gof(g_mm1)
#On average across the fish
confint(g_mm1)
Anova(g_mm1)
plot(g_mm1)
#Compare
anova(g_mm0, g_mm1)

# Two way mixed anova -----------------------------------------------------
#https://www.datanovia.com/en/lessons/mixed-anova-in-r/
res.aov2 <- anova_test(
  data = rep_model_data%>%filter(Method %in%c("d_BIM","d_BIMF","d_Fraslee")), dv = rate, wid = Fish_ID,
  between = Method, within = sample_point
)
get_anova_table(res.aov2)

#Average percent errors
mpe_table <- rep_model_data %>%
  group_by (Method)%>%
  mutate(mean_rate= mean(rate, na.rm=T))%>%
  distinct(Method, mean_rate)


# Application to Wild fish ------------------------------------------------
#This dataframe contains the wild fish and cage fish as comparisons
oto_growth_lifetime <-otofish %>%
  arrange(Fish_ID, Inc_no)%>%
  group_by(Fish_ID)%>%
  #Calculate growth rate for whole life of the fish
  #use lag1 for daily growth rates
  mutate(growth_fraslee=FL_fraslee-lag(FL_fraslee, n=1L))%>%
  mutate(growth_bim=FL_bim-lag(FL_bim, n=1L))%>%
  mutate(growth_bimf=FL_bimf-lag(FL_bimf, n=1L))%>%
  ungroup()%>%
  # filter(backinc>0)%>%
  filter(Inc_no>0)%>%
  #Calculate life time average growth rates for each fish
  group_by(Fish_ID)%>%
  mutate(growth_fraslee_lifetime=mean(growth_fraslee, na.rm=T),
         growth_bim_lifetime=mean(growth_bim, na.rm=T),
         growth_bimf_lifetime=mean(growth_bimf, na.rm=T))%>%
  ungroup()%>%
  mutate(Type=factor(Type, levels=c("Wild", "Wetland","Agriculture","Channel")))%>%
  filter(!is.na(Type))%>%
  ungroup()%>%
  mutate(region_label= case_when (Region =="Butte Sink" ~"Butte Sink Enclosures",
                                  Region =="Butte Sink Wetland" ~"Butte Sink Wild",
                                  Region =="Feather River" ~"Feather River Enclosures",
                                  Region =="Sacramento River" ~"Sacramento River Enclosures",
                                  Region =="Lower Bypass" ~"Lower Bypass Enclosures",
                                  Region =="Upper Bypass" ~"Upper Bypass Enclosures",
                                  Region =="Sutter Bypass" ~"Upper and Lower Bypass Wild"
  ))%>%
  mutate(region_label=factor(region_label, levels=c("Butte Sink Enclosures","Butte Sink Wild",
                                                    "Upper Bypass Enclosures","Lower Bypass Enclosures",
                                                    "Upper and Lower Bypass Wild",
                                                    "Feather River Enclosures","Sacramento River Enclosures")))%>%
  mutate(capture_location=case_when(Location=="Lundberg Farms phase I field" ~"Butte Creek watershed",
                                    Location=="Boat channel at Sanborn Slough" ~"Butte Creek watershed",
                                    Location=="North wetland at Mallard Ranch" ~"Butte Creek watershed",
                                    Location=="BSSS" ~"Butte Creek watershed",
                                    Location=="SBSNWR" ~"Butte Creek watershed",
                                    Location=="BCMR" ~"Butte Creek watershed",
                                    Location=="BCLR" ~"Butte Creek watershed",
                                    Location=="East margin of Sutter Bypass" ~"Butte Creek watershed",
                                    Location=="RIC" ~"Butte Creek watershed",
                                    Location=="SRCCL" ~"Sacramento River",
                                    Location=="SRSCCW" ~"Sacramento River",
                                    TRUE~Location))


#Summary table
sutter_summary_table <-oto_growth_lifetime %>%
  distinct(Fish_ID, Year, Experiment, Location, FL_mm, Inc_no_max, Inc_distance_max)%>%
  arrange(Experiment, Year)

write.csv(sutter_summary_table, "output/sutter_summary_table.csv", row.names=F)


#Lifetime growth reconstructions for wild fish 
#Using MF
p_growth_heat_date_loc <-ggplot (data=oto_growth_lifetime%>%filter(Type=="Wild")
                                )+
  geom_tile(aes(x=Inc_no, y=Fish_ID, fill=growth_bimf, group=Fish_ID))+
  # scale_x_date(name="Date")+
  scale_x_continuous(name="Increment number (days)", expand=c(0,0))+
  scale_y_discrete(name="Fish ID")+
  scale_fill_viridis(name="Growth rate (mm/d)")+
  theme_bw()+
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14))+
  # facet_wrap(~Year, scales="free")
  facet_col(vars(Year, capture_location), scales = "free_y", space = "free",
            strip.position = "right", labeller = label_wrap_gen(width=10),)
p_growth_heat_date_loc
ggsave(plot=p_growth_heat_date_loc, "output/growth_heat_date_loc.png", width = 9, height = 11)

#As a line plot
p_growth_comb_line <-ggplot (data=oto_growth_lifetime)+
  geom_line(aes(x=Inc_no, y=growth_bimf, group=Fish_ID,color=Region))+
  scale_x_continuous(name="Increment number (days)", expand=c(0,0))+
  scale_y_continuous(name="Growth rate (mm/d)")+
  theme_bw()+
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6))+
  facet_grid(Year~Experiment)

p_growth_comb_line 
ggsave(plot=p_growth_comb_line, "output/growth_comb_line.png", width = 8, height = 8)

#Growth rates by days
p_growth_age_summary<- ggplot(data=oto_growth_lifetime) +
  geom_line(aes(x=Inc_no,y=growth_bimf, group=Fish_ID), color="grey")+
  geom_smooth(aes(x=Inc_no,y=growth_bimf, group=Type, color=Type), method="loess", span=0.5)+
  scale_y_continuous(name="Growth rate (mm)")+
  scale_x_continuous(name="Increment number (days)")+
  theme_bw() +
  theme(legend.position="right", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_age_summary
ggsave(plot=p_growth_age_summary, "output/growth_age_summary.png", width = 6, height = 4)


# #Compare wild growth and observed growth categories OR type -------------
#Lifetime
wild_join_growth <-oto_growth_lifetime %>%
  filter(Type=="Wild")%>%
  group_by(Fish_ID)%>%
  mutate(Growth_rate_mean=mean(growth_bimf, na.rm=T))%>%
  ungroup()%>%
  distinct(Fish_ID,Year,Report_ID,Type,Growth_cat="Wild",
           Method="BIM_MF", Growth_rate_mean, Location, Inc_no_max, Inc_distance_max, FL_mm)%>%
  mutate(capture_location=case_when(Location=="Lundberg Farms phase I field" ~"Butte Creek watershed",
                                    Location=="Boat channel at Sanborn Slough" ~"Butte Creek watershed",
                                    Location=="North wetland at Mallard Ranch" ~"Butte Creek watershed",
                                    Location=="BSSS" ~"Butte Creek watershed",
                                    Location=="SBSNWR" ~"Butte Creek watershed",
                                    Location=="BCMR" ~"Butte Creek watershed",
                                    Location=="BCLR" ~"Butte Creek watershed",
                                    Location=="East margin of Sutter Bypass" ~"Butte Creek watershed",
                                    Location=="RIC" ~"Butte Creek watershed",
                                    Location=="SRCCL" ~"Sacramento River",
                                    Location=="SRSCCW" ~"Sacramento River",
                                    TRUE~Location))%>%
  dplyr::select(-Location)


growth_cagewildcomb<- oto_flcomp_mean%>% 
  filter(Method =="observed_growth_mean")%>%
  mutate(capture_location=Type)%>%
  rbind(wild_join_growth)%>%
  mutate(capture_location=factor(capture_location, levels=c( "Channel","Agriculture","Wetland",
                                                            "Sacramento River","Butte Creek watershed")))

p_growth_cagewildcomb<- ggplot(data=growth_cagewildcomb) +
  geom_boxplot(aes(x=Growth_cat, y=Growth_rate_mean,  fill=Growth_cat), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"))+
  geom_point(aes(x=Growth_cat,y=Growth_rate_mean, fill=Growth_cat), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_x_discrete(name="Growth category")+
  scale_y_continuous(name="Growth rate (mm/day)")+
  scale_fill_manual (name="Growth Category", values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  scale_color_manual (name="Growth Category",values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_cagewildcomb
ggsave(plot=p_growth_cagewildcomb, "output/growth_cagewildcomb.png", width = 5, height = 4)

## average growths by cage growth potential
average_growth_cat <- growth_cagewildcomb %>%
  group_by(Growth_cat) %>%
  summarize(average_growth_rate = mean(Growth_rate_mean))
average_growth_cat

p_growth_cagewildcomb_type<- ggplot(data=growth_cagewildcomb) +
  geom_boxplot(aes(x=capture_location, y=Growth_rate_mean,  fill=Type), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5, width=0.5)+
  geom_point(aes(x=capture_location,y=Growth_rate_mean, fill=Type), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_x_discrete(name="Location")+
  scale_y_continuous(name="Growth rate (mm/day)")+
  # scale_fill_manual (name="Growth Category", values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  # scale_color_manual (name="Growth Category",values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  scale_fill_manual (name="Growth Category", values=c("grey20","grey50","grey80", "firebrick"))+
  scale_color_manual (name="Growth Category", values=c("grey20","grey50","grey80", "firebrick"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_cagewildcomb_type
ggsave(plot=p_growth_cagewildcomb_type, "output/growth_cagewildcomb_type.png", width = 8, height = 5)
write.csv(growth_cagewildcomb,"output/growth_cagewildcomb.csv", row.names=F)

## stats for the location-based comparison w pwc / anova

growth_comb_model_location <- welch_anova_test(Growth_rate_mean ~ capture_location, data = growth_cagewildcomb)
growth_comb_model_location
#significant

# Pairwise comparisons
growth_comb_model_pwc_location <- growth_cagewildcomb %>%  games_howell_test(Growth_rate_mean ~ capture_location)
growth_comb_model_pwc_location

# pwc plot data group
pwc_plot_data_group_location <-growth_comb_model_pwc_location %>% add_xy_position(x = "capture_location")
pwc_plot_group_location <-ggboxplot(growth_cagewildcomb, x = "capture_location", y = "Growth_rate_mean",
                                    fill='Type', group="capture_location", alpha=0.5) +
  stat_pvalue_manual(pwc_plot_data_group_location, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(growth_comb_model_mt, detailed = TRUE),
    caption = get_pwc_label(pwc_plot_data_group)
  )+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Type")+
  scale_fill_manual(label=c("Channel","Agriculture","Wetland","Sacramento River", "Butte Creek watershed"),values=c("grey20","grey50","grey80", "firebrick"))+
  scale_color_manual(label=c("Channel","Agriculture","Wetland","Sacramento River", "Butte Creek watershed"),values=c("grey20","grey50","grey80", "firebrick"))+
  theme(legend.position="top", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.position = 'none')
pwc_plot_group_location

ggsave("output/pwc_plot_group_location.png", pwc_plot_group_location,  width = 10, height = 7)


### average mean for growth of each cage / wild capture location
average_growth <- growth_cagewildcomb %>%
  group_by(capture_location) %>%
  summarize(average_growth_rate = mean(Growth_rate_mean))
average_growth

growth_cagewildcomb_summary<- growth_cagewildcomb%>%
  group_by(Growth_cat, Year, capture_location)%>%
  mutate(N=n())%>%
  distinct(Growth_cat, Year, capture_location, N)%>%
  dplyr::select(Growth_cat, Year, capture_location, N)%>%
  arrange(Growth_cat, Year, capture_location)
write.csv(growth_cagewildcomb_summary,"output/growth_cagewildcomb_summary.csv", row.names=F)

##### all caged vs all wild growth rate comparison
growth_cagewildcomb_twocat <- growth_cagewildcomb %>%
  mutate(combined_type = ifelse(Type %in% c("Agriculture", "Channel", "Wetland"), "Caged", "Wild"))

p_growth_cagewildcomb_twocat<- ggplot(data=growth_cagewildcomb_twocat) +
  geom_boxplot(aes(x=combined_type, y=Growth_rate_mean,  fill=combined_type), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5, width=0.5)+
  geom_point(aes(x=combined_type,y=Growth_rate_mean, fill=combined_type), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_x_discrete(name="Type")+
  scale_y_continuous(name="Growth rate (mm/day)")+
  # scale_fill_manual (name="Growth Category", values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  # scale_color_manual (name="Growth Category",values=c("firebrick","#E69F01","#56B4E9", "#228B22"))+
  scale_fill_manual (name="Growth Category", values=c("grey20", "firebrick"))+
  scale_color_manual (name="Growth Category", values=c("grey20", "firebrick"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
p_growth_cagewildcomb_twocat

## standard deviation for wild fish
sd_wild <- growth_cagewildcomb_twocat %>%
  group_by(combined_type) %>%
  summarize(sd_growth_rate = sd(Growth_rate_mean))
sd_wild

p_cagewildcomb_grid <-plot_grid(p_growth_cagewildcomb_twocat+theme(legend.position="none"),
                                p_growth_cagewildcomb_type+theme(legend.position="none"),
                                labels= c("A", "B"), rel_widths = c(2,3))

p_cagewildcomb_grid

## Show what wild fish data would look like with the other models
oto_flcomp_wild <- otofish %>%
  distinct(Fish_ID, backinc, .keep_all=T) %>%
  filter (Experiment=="Wild")%>%
  #Determine dates between field measurements
  group_by(Fish_ID)%>%
  arrange(Fish_ID,desc(backdate))%>%
  mutate(growth_days=lag(backdate, n=1L)-backdate)%>%
  mutate(growth_days=as.numeric(growth_days))%>%
  #Determine FL growth fraslee
  mutate(FL_growth_fraslee=lag(FL_fraslee, n=1L)-FL_fraslee)%>%
  #Determine FL growth bim
  mutate(FL_growth_bim=lag(FL_bim, n=1L)-FL_bim)%>%
  #Determine FL growth mf bim
  mutate(FL_growth_bimf=lag(FL_bimf, n=1L)-FL_bimf)%>%
  ungroup()%>%
  #Remove the last value
  filter(backinc>0)%>%
  #Determine Growthrates
  group_by(Fish_ID)%>%
  #fraslee growth rates
  mutate(fraslee_growth=FL_growth_fraslee/growth_days)%>%
  mutate(fraslee_growth_mean=(max(FL_fraslee)-min(FL_fraslee))/as.numeric(max(backdate)-min(backdate)))%>%
  #bim growth rates
  mutate(bim_growth=FL_growth_bim/growth_days)%>%
  mutate(bim_growth_mean=(max(FL_bim)-min(FL_bim))/as.numeric(max(backdate)-min(backdate)))%>%
  #mf bim growth rates
  mutate(bimf_growth=FL_growth_bimf/growth_days)%>%
  mutate(bimf_growth_mean=(max(FL_bimf)-min(FL_bimf))/as.numeric(max(backdate)-min(backdate)))%>%
  ungroup()%>%
  mutate(Growth_cat=case_when (bimf_growth_mean >0.5 ~ "High",
                               bimf_growth_mean <=0.5 & bimf_growth_mean >0.2 ~"Medium",
                               bimf_growth_mean <=0.2 ~ "Low",
                               Type=="Wild"~"Wild"))%>%
  mutate(Growth_cat=factor(Growth_cat, levels=c("Low","Medium","High", "Wild")))
  
#### nothing in the low category-- changing above code to exclude it

oto_flcomp_wild <- otofish %>%
  distinct(Fish_ID, backinc, .keep_all=T) %>%
  filter (Experiment=="Wild")%>%
  #Determine dates between field measurements
  group_by(Fish_ID)%>%
  arrange(Fish_ID,desc(backdate))%>%
  mutate(growth_days=lag(backdate, n=1L)-backdate)%>%
  mutate(growth_days=as.numeric(growth_days))%>%
  #Determine FL growth fraslee
  mutate(FL_growth_fraslee=lag(FL_fraslee, n=1L)-FL_fraslee)%>%
  #Determine FL growth bim
  mutate(FL_growth_bim=lag(FL_bim, n=1L)-FL_bim)%>%
  #Determine FL growth mf bim
  mutate(FL_growth_bimf=lag(FL_bimf, n=1L)-FL_bimf)%>%
  ungroup()%>%
  #Remove the last value
  filter(backinc>0)%>%
  #Determine Growthrates
  group_by(Fish_ID)%>%
  #fraslee growth rates
  mutate(fraslee_growth=FL_growth_fraslee/growth_days)%>%
  mutate(fraslee_growth_mean=(max(FL_fraslee)-min(FL_fraslee))/as.numeric(max(backdate)-min(backdate)))%>%
  #bim growth rates
  mutate(bim_growth=FL_growth_bim/growth_days)%>%
  mutate(bim_growth_mean=(max(FL_bim)-min(FL_bim))/as.numeric(max(backdate)-min(backdate)))%>%
  #mf bim growth rates
  mutate(bimf_growth=FL_growth_bimf/growth_days)%>%
  mutate(bimf_growth_mean=(max(FL_bimf)-min(FL_bimf))/as.numeric(max(backdate)-min(backdate)))%>%
  ungroup()%>%
  mutate(Growth_cat=case_when (bimf_growth_mean >0.5 ~ "High",
                               bimf_growth_mean <=0.5  ~"Medium",
                               Type=="Wild"~"Wild"))%>%
  mutate(Growth_cat=factor(Growth_cat, levels=c("Medium","High", "Wild")))


#Export otoflcomp dataframe 
write.csv(oto_flcomp_wild, "output/oto_flcomp_wild.csv", row.names=F)

oto_flcomp_mean_wild <- oto_flcomp_wild %>%
  distinct(Fish_ID, Year, Type, Report_ID, FL_mm, Inc_no_max, Inc_distance_max, fraslee_growth_mean, bim_growth_mean, bimf_growth_mean, Growth_cat) %>%
  pivot_longer(cols = c(fraslee_growth_mean, bim_growth_mean, bimf_growth_mean), names_to = "Method", values_to = "Growth_rate_mean") %>%
  ungroup() %>%
  mutate(Method = factor(Method, levels = c("fraslee_growth_mean", "bim_growth_mean", "bimf_growth_mean")))

write.csv(oto_flcomp_mean_wild, "output/oto_flcomp_mean_wild.csv", row.names=F)


# graph of average growths for wild fish reconstructed from 3 models------------
p_growth_comp_wild<- ggplot(data=oto_flcomp_mean_wild) +
  geom_boxplot(aes(x=Method, y=Growth_rate_mean,  fill=Method), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5)+
  geom_point(aes(x=Method, y=Growth_rate_mean, fill=Method), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Method", labels=c("Fraslee","BI","MF"))+
  scale_fill_manual(label=c("Fraslee","BI","MF"),values=c("purple","orange","darkgreen"))+
  scale_color_manual(label=c("Fraslee","BI","MF"),values=c("purple","orange","darkgreen"))+
  theme_bw() +
  theme(legend.position="none", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_comp_wild
ggsave(plot=p_growth_comp_wild, "output/p_growth_comp_wild.png", width = 10, height = 7)

# graph of average growths delineated by growth category
# note: nothing fits low growth category?
p_growth_comp_type_wild<- ggplot(data=oto_flcomp_mean_wild) +
  geom_boxplot(aes(x=Growth_cat, y=Growth_rate_mean, group=interaction(Method,Growth_cat), fill=Method), outlier.shape=NA,
               position = position_dodge(1, preserve = "single"), alpha=0.5)+
  geom_point(aes(x=Growth_cat, y=Growth_rate_mean, group=interaction(Method,Growth_cat), fill=Method), color="black", shape=21,
             position = position_jitterdodge(dodge.width=1,jitter.width=0.1), size=1.5)+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Growth potential")+
  scale_fill_manual(label=c("Fraslee","BI","MF"),values=c("purple","orange","darkgreen"))+
  scale_color_manual(label=c("Fraslee","BI","MF"),values=c("purple","orange","darkgreen"))+
  theme_bw() +
  theme(legend.position="top", panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_growth_comp_type_wild
ggsave(plot=p_growth_comp_type_wild, "output/p_growth_comp_type_wild.png", width = 10, height = 7)

p_comp_comb_wild <- plot_grid(p_growth_comp_wild, p_growth_comp_type_wild, labels = c("A", "B"), ncol = 2)
ggsave(plot=p_comp_comb_wild, "output/comp_comb_wild.png", width = 15, height = 7)


# stats for comparing reconstructions in wild fish-----------------------
## testing assumptions for using a parametric vs nonparametric test:

# assumption 1: no repeat measures = met

# assumption 2: no significant outliers
oto_flcomp_mean_wild %>% 
  group_by(Method, Growth_cat) %>%
  identify_outliers(Growth_rate_mean)
# we have 6 outliers

#Assumption 3: Normality
# Build the linear model
model_fornorm_wild  <- lm(Growth_rate_mean ~ Method, data = oto_flcomp_mean_wild)
# Create a QQ plot of residuals
ggqqplot(residuals(model_fornorm_wild))
#Create a QQ plot of residuals by Method
ggqqplot(oto_flcomp_mean_wild, "Growth_rate_mean", facet.by = "Method")
# Compute Shapiro-Wilk test of normality
shapiro_test(residuals(model_fornorm))
# Compute Shapiro-Wilk test of normality by Method and Growth_cat
oto_flcomp_mean_wild %>% 
  group_by(Method, Growth_cat) %>%
  shapiro_test(Growth_rate_mean)
#shapiro-wilk test =  we can't assume normality

#Assumption 4: Homogeneity of variances
#residuals versus fits plot
plot(model_fornorm_wild, 1)
#Levene test
oto_flcomp_mean_wild %>% 
  levene_test(Growth_rate_mean ~ as.factor(Method))

#Homogeneity of variances is not met. Ran a Kruskal-Wallis test
##########come back here to return back
growth_comb_model_wild_kwt <- kruskal.test(Growth_rate_mean ~ Method, data = oto_flcomp_mean_wild)
growth_comb_model_wild_kwt
# not significant

growth_comb_model_wild <- welch_anova_test(Growth_rate_mean ~ Method, data = oto_flcomp_mean)
growth_comb_model_wild

growth_comb_model_pwc_wild <- oto_flcomp_mean_wild %>%  games_howell_test(Growth_rate_mean ~ Method)
growth_comb_model_pwc_wild

growth_comb_model_mt_wild <- welch_anova_test(Growth_rate_mean ~ Method+Growth_cat, data = oto_flcomp_mean)
growth_comb_model_mt_wild

growth_comb_model_mt_pwc_wild <- oto_flcomp_mean_wild %>% 
  group_by(Growth_cat)%>%
  games_howell_test(Growth_rate_mean ~ Method)
growth_comb_model_mt_pwc_wild

pwc_plot_data_wild <-growth_comb_model_pwc_wild%>% add_xy_position(x = "Method")
pwc_plot_wild <-ggboxplot(oto_flcomp_mean_wild, x = "Method", y = "Growth_rate_mean",
                          fill='Method', alpha=0.5) +
  stat_pvalue_manual(pwc_plot_data_wild, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(growth_comb_model_wild, detailed = TRUE),
    # caption = get_pwc_label(pwc_plot_data)
  )+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Method", label=c("F-L","BI","MF"))+
  scale_fill_manual(label=c("F-L","BI","MF"),values=c("purple","orange","darkgreen"))+
  scale_color_manual(label=c("F-L","BI","MF"),values=c("purple","orange","darkgreen"))+
  theme(legend.position="top", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pwc_plot_wild

pwc_plot_data_group_wild <-growth_comb_model_mt_pwc_wild %>% add_xy_position(x = "Growth_cat")

pwc_plot_group_wild <-ggboxplot(oto_flcomp_mean_wild, x = "Growth_cat", y = "Growth_rate_mean",
                                fill='Method', group="Growth_cat", alpha=0.5) +
  stat_pvalue_manual(pwc_plot_data_group_wild, hide.ns = TRUE, position = position_dodge(width = 0.8)) +
  labs(
    subtitle = get_test_label(growth_comb_model_mt_wild, detailed = TRUE),
    caption = get_pwc_label(pwc_plot_data_group_wild)
  )+
  scale_y_continuous(name="Growth (mm/day)")+
  scale_x_discrete(name="Growth Category")+
  scale_fill_manual(label=c("F-L","BI","MF"),values=c("purple","orange","darkgreen"))+
  scale_color_manual(label=c("F-L","BI","MF"),values=c("purple","orange","darkgreen"))+
  theme(legend.position="top", panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
pwc_plot_group_wild

p_growth_pwc_wild <-plot_grid(pwc_plot_wild+theme(legend.position="none"),
                         pwc_plot_group_wild+theme(legend.position="none"), 
                         labels = c("A", "B"), ncol=1)
p_growth_pwc_wild
ggsave("output/growth_pwc.png", p_growth_pwc,  width = 6, height = 10)


# Growth comparison among habitats all fish including age correction------------------------------
growthcomb_m_all <- oto_growth_lifetime%>%
  dplyr::select(Fish_ID, Inc_no,growth_bimf,Type, Location)%>%
  mutate(Fish_ID=as.factor(Fish_ID))%>%
  mutate(capture_location=case_when(
    Type=="Wild" & Location=="Lundberg Farms phase I field" ~"Butte Creek watershed",
    Type=="Wild" & Location=="Boat channel at Sanborn Slough" ~"Butte Creek watershed",
    Type=="Wild" &Location=="North wetland at Mallard Ranch" ~"Butte Creek watershed",
    Type=="Wild" &Location=="BSSS" ~"Butte Creek watershed",
    Type=="Wild" &Location=="SBSNWR" ~"Butte Creek watershed",
    Type=="Wild" &Location=="BCMR" ~"Butte Creek watershed",
    Type=="Wild" &Location=="BCLR" ~"Butte Creek watershed",
    Type=="Wild" &Location=="East margin of Sutter Bypass" ~"Butte Creek watershed",
    Type=="Wild" &Location=="RIC" ~"Butte Creek watershed",
    Type=="Wild" &Location=="SRCCL" ~"Sacramento River",
    Type=="Wild" &Location=="SRSCCW" ~"Sacramento River",
                             TRUE~Type))%>%
  mutate(capture_location=factor(capture_location, levels=c( "Channel","Agriculture","Wetland",
                                                             "Sacramento River","Butte Creek watershed")))

#Fit gam model
m2 <- bam(growth_bimf ~ 
            s(Inc_no, k=10) +
            s(Fish_ID, bs = "re") + 
            capture_location,
          data = growthcomb_m_all, discrete = FALSE, method = "ML", select = TRUE)

#Explore results
summary(m2)
plot(m2, all.terms=T, pages=1, rug=T)

#Compare movement class
mc_comb_all <- emmeans(m2, pairwise~capture_location, data=growthcomb_m_all)
mc_comb_all
mc_comb_all_plot <-plot(mc_comb_all, 
                        ylab="Habitat", 
                        xlab="Estimated marginal mean growth rates (mm/d)",
                        color="grey8")+
  theme_classic()
mc_comb_all_plot
ggsave(plot=mc_comb_all_plot, "output/mc_comb_all_plot.png", width = 6, height = 5)

#Emmeans estimates and pairwise table 
emmeans_pairwise_table <-as.data.frame(mc_comb_all$contrasts)
emmeans_estimates_table<-as.data.frame(mc_comb_all$emmeans)
write.csv(emmeans_pairwise_table,"output/emmeans_pairwise_table.csv",row.names=F)
write.csv(emmeans_estimates_table,"output/emmeans_estimates_table.csv",row.names=F)
