### Load packages
require(data.table)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggExtra)
library(BIEN)
library(smatr)
require(data.table)
library(hrbrthemes)

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                          PREPARE DATA FROM GHANA                         ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



rm(list = ls())
#setwd("F:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data")
setwd("~/GitHub/Ghana_rainfall_trait_variation_optimality_github/input_data")

##---------------------------------------------------------------
##                          phenology                           -
##---------------------------------------------------------------

#pheno<- openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Decideous_or_evergreen/leaftr_ghana_Jesus_051017_Phenology2.xlsx')%>%
pheno<- openxlsx::read.xlsx('leaftr_ghana_Jesus_051017_Phenology2.xlsx')%>%
  dplyr::rename(Species=!!as.name("Species_JA_trimfinalspace"))%>%
  replace_na(list(Phenology='unkown_in_Jesus'))%>%
  mutate(Phenology=stringi::stri_replace_all_charclass(Phenology, "\\p{WHITE_SPACE}", ""))%>%
  mutate(Phenology=recode(Phenology,"?ed"='evergreen',"ed"='evergreen',"de"='deciduous'))
unique(pheno$Phenology)

##----------------------------------------------------------------
##                              Ks                               -
##----------------------------------------------------------------
List_of_site<-c('KOG','BOB','ANK')

#WD_Ks<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Ks')%>%
WD_Ks<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  dplyr::rename(Ks=!!as.name("theoretical.conductivity.poorter2010.kg.m-1.s-1.Mpa-1"))%>%
  dplyr::rename(Vessel_diameter = !!as.name('hydraulically.weighted.diam.micron'))


#the left join will bring species to the twig density table
#Census_big_table<-read.csv("F:/Oxford/Chapter_one/From_imma/Traits/Census_data/Census_table_all_plots_Species_list.csv")%>%
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
WD_Ks<-left_join(WD_Ks,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
WD_Ks<-drop_na(WD_Ks,any_of('Species'))%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%left_join(pheno[,c('Species','Phenology')],by="Species")


##---------------------------------------------------------------
##                          Huber value                         -
##---------------------------------------------------------------


#Huber<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
Huber<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('SapwoodArea/LeafArea_mm2/cm2'))%>%
  mutate(Huber=Huber*10000) # get into mm2/m2 then in the following calculation, it will change into cm2/m2 when /100

##----------------------------------------------------------------
##                          Npercent                         -
##----------------------------------------------------------------

#NNNNN=openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Npercent')%>%
NNNNN=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Npercent')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE","Branch"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Plot_code=!!as.name("PLOT"))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")%>%
  left_join(pheno[,c('Species','Phenology')],by="Species")%>%
  mutate(sun_shade=str_detect(Branch, ".*SH"))%>%
  filter(!(sun_shade))%>%
  mutate(sun_shade='Sun')

NNNNN_lma<- NNNNN %>% drop_na(LMA)%>%
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  dplyr::select(Site,Species,LMA,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(sqrt(LMA),type="iqr"))%>%
  filter(Outlier_id<1.5)

NNNNN_leaf_area<- NNNNN %>% 
  drop_na(LeafArea.m2)%>%
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  dplyr::select(Site,Species,LeafArea.m2,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(sqrt(LeafArea.m2),type="iqr"))%>%
  filter(Outlier_id<1.5)

NNNNN_LeafDrymass<- NNNNN %>% 
  drop_na(LeafDrymass.g)%>%
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  dplyr::select(Site,Species,LeafDrymass.g,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(log2(LeafDrymass.g),type="iqr"))%>%
  filter(Outlier_id<1.5)


hist(log2(NNNNN_LeafDrymass$LeafDrymass.g))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                          THE DATASET BY MAURIZIO                         ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#setwd('F:/Oxford/Chapter_one/TRY_global')

#Maurizio_data<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/TRY_global/Huber_values_HuanYuan.xlsx')%>%
Maurizio_data<-openxlsx::read.xlsx('Huber_values_HuanYuan.xlsx')%>%
  rename(AccSpeciesName=Species)%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=(exp(LogHv)),
         Ks=(exp(LogKs)),
         leaf_area_m2=(exp(Log_wlsize)),
         SLA_m2_kg1=(exp(LogSLA)))
species_vector<-unique(Maurizio_data$AccSpeciesName)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Ks, log = "xy", data = Maurizio_data)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Sapwood_per_leaf_area_cm2_m2_ks.png", width = 400, height = 400)
plot(ftComSlope, pch=20, cex=0.5, col="black",
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     cex.lab = 1.5)
dev.off()



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              vessel diameter                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#median_wd <- BIEN_trait_mean(species = species_vector, trait = "vessel lumen area")
#save(median_wd,file='vessel_lumen_area_Maurizio_data.rda')
load('vessel_lumen_area_Maurizio_data.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_vessel<-median_wd%>%
  rename(AccSpeciesName=species)%>%
  mutate(Vessel_diameter_mm=sqrt(mean_value/pi))%>%
  filter(level_used=='Species' | level_used=='Genus')%>%
  drop_na()%>%
  dplyr::select(AccSpeciesName,Vessel_diameter_mm)

To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  dplyr::select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)

To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  dplyr::select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
To_draw=rbind(To_draw1,To_draw2)

pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology')],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  # filter(Ks < 250)%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
         Vessel_diameter_mm=Vessel_diameter*0.001
  )%>%
  dplyr::select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Vessel_diameter_mm)%>%
  mutate(data_source='Ghana')%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            Vessel_diameter_mm=mean(Vessel_diameter_mm,na.rm=T),
            data_source=first(data_source))


rbind_to_draw<-rbind(pairwise_table,To_draw)


ggplot(rbind_to_draw,aes(x=log10(Vessel_diameter_mm), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)


ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Vessel_diameter_mm*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Sapwood_per_leaf_area_cm2_m2_Vessel_diameter.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     xlab = "Vessel diameter (mm) [log scale]",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend(x="bottomleft",
       legend=c('Ghana rainfall gradient','(Mencuccini et al., 2019)','TROP_FOR in (Mencuccini et al., 2019)'),
       col=c("blue","orange","red"),
       lty=1,
       bty='n',
       cex=1.1)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              KS                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology')],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  # filter(Ks < 250)%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
         Vessel_diameter_mm=Vessel_diameter*0.001,
         Ks=Ks*0.1
  )%>%
  filter(Ks<20)%>%
  dplyr::select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Ks,SITE)%>%
  mutate(data_source='Ghana')%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            Ks=mean(Ks,na.rm=T),
            data_source=first(data_source))


To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='(Mencuccini et al., 2019)')%>%
  dplyr::select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)

To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR in \n (Mencuccini et al., 2019)')%>%
  dplyr::select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
To_draw=rbind(To_draw1,To_draw2)

rbind_to_draw<-rbind(pairwise_table,To_draw)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Ks*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Sapwood_per_leaf_area_cm2_m2_Ks2.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     xlab = "Ks [log scale]",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend(x="bottomleft",
       legend=c('Ghana rainfall gradient','(Mencuccini et al., 2019)','TROP_FOR in (Mencuccini et al., 2019)'),
       col=c("blue","orange","red"),
       lty=1,
       bty='n',
       cex=1.1)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Coefficient of Variation                        ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Coefficient_of_Variation<-function(data) {
  data<-as.numeric(data)
  sd(data) / mean(data) * 100
}

rbind_to_draw2<-rbind_to_draw[complete.cases(rbind_to_draw),]%>%
  mutate( Ks_timed_vH = Ks*Sapwood_per_leaf_area_cm2_m2)%>%
  pivot_longer(cols=c('Sapwood_per_leaf_area_cm2_m2','Ks','Ks_timed_vH'))%>%
  group_by(data_source,name)%>%
  summarise(CV=Coefficient_of_Variation(value))

library(ggpattern)

ggplot(rbind_to_draw2, aes(x = factor(data_source), y = CV,  fill = name)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+
  labs( x = 'Source of data', y = 'Coefficient of Variation(%)', fill = "Traits")+
  scale_fill_manual(values=c("blue","red","darkgreen"),labels = c("Ks", "Ks*AS/AL","AS/AL"))+
  theme(legend.position="bottom",text = element_text(size=12)) 

ggsave('../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/CoefficientofVariation_ks_huber.jpeg',width=4.22 , height=3.56,dpi=700)

inversed_rbind_to_draw<-rbind_to_draw%>%mutate(inversed_Ks=1/Ks)
ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ inversed_Ks*data_source, log = "xy", data = inversed_rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2) [log scale]",
     xlab = "inversed Ks [log scale]",
     cex.lab = 1.7,
     cex.axis = 1.1)



pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology')],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  # filter(Ks < 250)%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
         Vessel_diameter_mm=Vessel_diameter*0.001,
         Ks=Ks*0.1
  )%>%
  filter(Ks<20)%>%
  dplyr::select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Ks,SITE)%>%
  mutate(data_source=SITE)%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            Ks=mean(Ks,na.rm=T),
            data_source=first(data_source))


To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  dplyr::select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)

To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  dplyr::select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
To_draw=rbind(To_draw1,To_draw2)

rbind_to_draw<-rbind(pairwise_table,To_draw)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Ks*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Sapwood_per_leaf_area_cm2_m2_Ks2_plot_scale.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     xlab = "Ks [log scale]",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend(x="bottomleft",
       legend=c('KOG(dry)','BOB(mid)','ANK(wet)','(Mencuccini et al., 2019)','TROP_FOR in (Mencuccini et al., 2019)'),
       col=c("red","darkgoldenrod3","blue","chartreuse4","cyan2"),
       lty=1,
       bty='n',
       cex=1.1)
dev.off()

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              SLA                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


pairwise_table<-full_join(Huber,NNNNN_lma[,c('PLOT_TREE',"LMA")],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  mutate(SLA_m2_kg1=1/(LMA/1000))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  filter(log10(SLA_m2_kg1) < 1.7)%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100)%>%
  rename(AccSpeciesName=Species )%>%
  dplyr::select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,SLA_m2_kg1)%>%
  mutate(data_source='Ghana')%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
            data_source=first(data_source))

To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  dplyr::select(SLA_m2_kg1,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)

To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  dplyr::select(SLA_m2_kg1,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
To_draw=rbind(To_draw1,To_draw2)


rbind_to_draw<-rbind(pairwise_table,To_draw)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ SLA_m2_kg1*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Sapwood_per_leaf_area_cm2_m2_SLA_m2_kg1.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     xlab = "SLA (m^2/kg) [log scale]",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend(x="bottomleft",
       legend=c('Ghana rainfall gradient','(Mencuccini et al., 2019)','TROP_FOR in (Mencuccini et al., 2019)'),
       col=c("blue","orange","red"),
       lty=1,
       bty='n',
       cex=1.1)
dev.off()



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              KS vs Twig Density                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
setwd("~/GitHub/Ghana_rainfall_trait_variation_optimality_github/input_data")

# Read in data:
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Twig_density_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")
WD_Ks<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  dplyr::rename(Ks=!!as.name("theoretical.conductivity.poorter2010.kg.m-1.s-1.Mpa-1"))%>%
  dplyr::rename(Vessel_diameter = !!as.name('hydraulically.weighted.diam.micron'))

# Merge the data:
pairwise_table<-full_join(Twig_density_table,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter',"vessel.per.area.mm2")],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  filter(Ks < 250)%>%
  filter(TwigDensity.g.cm3 < 1)

# Plot:
Main_plot<-ggplot(pairwise_table, aes(x=TwigDensity.g.cm3, y=log10(Ks), color=SITE)) + 
  geom_point(size=2,alpha=0.5) +
  theme_bw()+
  scale_color_brewer(palette="Dark2")+
  geom_smooth(method=lm,se=F, fullrange=TRUE)+ stat_ellipse(type = "norm")
Main_plot

# Make new plot
# Calculate the mean for each species at each site:
pairwise_table2 <- pairwise_table%>%group_by(Species, SITE)%>%
  summarise(mean_Ks=mean(Ks), mean_TD=mean(TwigDensity.g.cm3))
# Model:
ftComSlope = sma(mean_Ks ~ mean_TD * SITE, log = "xy", data = pairwise_table2)
summary(ftComSlope)
# Make the plot:
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/Ks_TwigDensity.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Kp (kg m-1Mpa-1s-1)",
     xlab = "Twig density (g/cm^3)",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend(x="bottomleft",
       legend=c('ANK','BOB','KOG'),
       col=c("blue","orange","red"),
       lty=1,
       bty='n',
       cex=1.1)
dev.off()


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              XFT based Kp versus AS/AL                   ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#median_vessel<-read.csv('F:/Oxford/Chapter_one/TRY_global/XFT_full_database_download_20211101-182931.csv')%>%
#median_vessel<-read.csv('XFT_full_database_download_20211101-182931.csv')
median_vessel<-openxlsx::read.xlsx('XFT_full_database_download_20211101-182931_kp_calculated.xlsx')%>%
  rename(Vessel_diameter_mm="Av..conduit.diameter...m.")%>% #this is um micrometer
  rename(vessel_per_area_m2="vessel.per.area.m-2")%>%
  rename(Kp="theoretical.conductivity.poorter2010")%>%
  mutate(Vessel_diameter_mm=Vessel_diameter_mm/1000,
         Kp=Kp*10)%>%
  mutate(AccSpeciesName=paste0(Genus,' ',Species))%>%
  dplyr::select(Vessel_diameter_mm,AccSpeciesName,vessel_per_area_m2,Kp)



To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='(Mencuccini et al., 2019)')%>%
  dplyr::select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName,vessel_per_area_m2,Kp)

To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR in (Mencuccini et al., 2019)')%>%
  dplyr::select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName,vessel_per_area_m2,Kp)
To_draw=rbind(To_draw1,To_draw2)

pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology','vessel.per.area.mm2')],by='PLOT_TREE')%>%
  mutate(SITE=substring(PLOT_TREE,1,3))%>%
  filter(Ks < 250)%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
         Vessel_diameter_mm=Vessel_diameter*0.001,
         vessel_per_area_m2 = vessel.per.area.mm2/1000000,
  )%>%
  dplyr::select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Vessel_diameter_mm,Ks,vessel_per_area_m2)%>%
  mutate(data_source='Ghana')%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            Vessel_diameter_mm=mean(Vessel_diameter_mm,na.rm=T),
            Kp=mean(Ks,na.rm=T),
            vessel_per_area_m2=mean(vessel_per_area_m2,na.rm=T),
            data_source=first(data_source))


rbind_to_draw<-rbind(pairwise_table,To_draw)


ggplot(rbind_to_draw,aes(x=log10(Kp), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)


ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Kp*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
png("../Ghana_rainfall_trait_variation_optimality/Ghana_global_compare/XFTbasedKpversusASAL.png", width = 600, height = 500)
par(mar = c(6,6,1,1))
plot(ftComSlope, pch=20, lwd=2, pt.cex=1.5,
     ylab = "Sapwood per leaf area (cm^2/m^2)",
     xlab = "Kp (kg m-1Mpa-1s-1)",
     cex.lab = 1.7,
     cex.axis = 1.1)
legend('bottomleft',
       legend=c('(Mencuccini et al., 2019)','Ghana rainfall gradient','TROP_FOR in (Mencuccini et al., 2019)'),
       col=c("blue","orange","red"),
       lty=1,
       bty='n',
       cex=1.1, pch = c(19, 19, 19))
dev.off()


Coefficient_of_Variation<-function(data) {
  data<-as.numeric(data)
  sd(data) / mean(data) * 100
}

rbind_to_draw2<-rbind_to_draw[complete.cases(rbind_to_draw),]%>%
  mutate( Kp_timed_vH = Kp*Sapwood_per_leaf_area_cm2_m2)%>%
  pivot_longer(cols=c('Sapwood_per_leaf_area_cm2_m2','Kp','Kp_timed_vH'))%>%
  group_by(data_source,name)%>%
  summarise(CV=Coefficient_of_Variation(value))

library(ggpattern)

ggplot(rbind_to_draw2, aes(x = factor(data_source), y = CV,  fill = name)) + 
  geom_bar(stat = "identity", position = "dodge", alpha = 0.5)+
  labs( x = 'Source of data', y = 'Coefficient of Variation(%)', fill = "Traits")+
  scale_fill_manual(values=c("blue","red","darkgreen"),labels = c("Kp", "Kp*AS/AL","AS/AL"))+
  theme(legend.position="bottom",text = element_text(size=12)) 
