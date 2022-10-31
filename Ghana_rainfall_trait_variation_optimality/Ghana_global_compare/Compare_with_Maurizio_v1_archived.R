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
#setwd("F:/Oxford/Chapter_one/From_imma/Traits/plotV3/")
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


#Huber<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Huber')%>%
Huber<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('LeafArea/SapwoodArea_m2/mm2'))%>%
  mutate(Huber=1/Huber)


#Huber<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
Huber<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('SapwoodArea/LeafArea_mm2/cm2'))%>%
  mutate(Huber=Huber*10000) # get into mm2/m2 then in the following calculation, it will change into cm2/m2 when /100

##----------------------------------------------------------------
##                        tugor loss point                       -
##----------------------------------------------------------------

List_of_site<-c('KOG','BOB','ANK')

#TLPoint<-openxlsx::read.xlsx("F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx",sheet='turgor')%>%
TLPoint<-openxlsx::read.xlsx("Overall_trait_raw_data.xlsx",sheet='turgor')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM","MONTH","YEAR"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))

TLPoint[87,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake
TLPoint[96,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake
TLPoint<-TLPoint%>%
  rename(Turgor=!!as.name("Psi_TLP.(MPa)"))


##----------------------------------------------------------------
##                          Twig density                         -
##----------------------------------------------------------------
List_of_site<-c('KOG','BOB','ANK')

#read data
#Twig_density_table<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
Twig_density_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")

##----------------------------------------------------------------
##                          Vcmax                         -
##----------------------------------------------------------------

List_of_site<-c('KOG','BOB','ANK')

#Vcmax_table<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='VCmax')%>%
Vcmax_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='VCmax')%>%
  separate("File_Name", sep='-',into=c("PLOT","TREE","STEM","TIME"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))
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
  select(Site,Species,LMA,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(sqrt(LMA),type="iqr"))%>%
  filter(Outlier_id<1.5)

NNNNN_leaf_area<- NNNNN %>% 
  drop_na(LeafArea.m2)%>%
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  select(Site,Species,LeafArea.m2,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(sqrt(LeafArea.m2),type="iqr"))%>%
  filter(Outlier_id<1.5)

NNNNN_LeafDrymass<- NNNNN %>% 
  drop_na(LeafDrymass.g)%>%
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  select(Site,Species,LeafDrymass.g,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=outliers::scores(log2(LeafDrymass.g),type="iqr"))%>%
  filter(Outlier_id<1.5)


hist(log2(NNNNN_LeafDrymass$LeafDrymass.g))
##----------------------------------------------------------------
##                          13C                                  -
##----------------------------------------------------------------

get_chi_from_13C<-function(delta_13C){
  compensation=42
  Ca=404
  lower_f=0.016
  asss=0.0044
  bsss=0.03
  delta_13C_air = -8.7 #see reference here https://zenodo.org/record/569501#.YGGCe68zbLs
  Big_delta_13C = (delta_13C_air - delta_13C)/(1+delta_13C/1000)
  Chi_yunke = (Big_delta_13C/1000 - asss + lower_f * compensation / Ca)/(bsss - asss)
  Chi_wang = (Big_delta_13C/1000 - asss) / (0.027 - asss) 
  return(list(Big_delta_13C,Chi_wang,Chi_yunke))
} #note that delta_13C is in unit 1/1000


#C13_table=openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='C13')%>%
C13_table=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='C13')%>%
  rename(C13=!!as.name("13CVPDB"))%>%
  separate("ID", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  mutate(ci_ca_ratio=get_chi_from_13C(C13)[[3]])%>%
  left_join(Census_big_table[,c('PLOT_TREE','Species')],by="PLOT_TREE")%>%
  left_join(pheno[,c('Species','Phenology')],by="Species")


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



Try_vessel_diameter <- fread("17021.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
      filter(TraitName=="Branch vessel diameter") %>%
  mutate(OrigValue=as.numeric(OrigValueStr)*0.001)%>%
   group_by(AccSpeciesName)%>%
  summarise(Vessel_diameter_mm2= mean(OrigValue,na.rm=T) )
#TRY does not match up with Maurizio well, give up

# XFT https://xylemfunctionaltraits.org/

#median_vessel<-read.csv('F:/Oxford/Chapter_one/TRY_global/XFT_full_database_download_20211101-182931.csv')%>%
median_vessel<-read.csv('XFT_full_database_download_20211101-182931.csv')%>%
  rename(Vessel_diameter_mm="Av..conduit.diameter...m.")%>% #this is um micrometer
  mutate(Vessel_diameter_mm=Vessel_diameter_mm/1000)%>%
  mutate(AccSpeciesName=paste0(Genus,' ',Species))%>%
  select(Vessel_diameter_mm,AccSpeciesName)


#median_wd <- BIEN_trait_mean(species = species_vector, trait = "vessel lumen area")
#save(median_wd,file='vessel_lumen_area_Maurizio_data.rda')
load('vessel_lumen_area_Maurizio_data.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_vessel<-median_wd%>%
  rename(AccSpeciesName=species)%>%
  mutate(Vessel_diameter_mm=sqrt(mean_value/pi))%>%
  filter(level_used=='Species' | level_used=='Genus')%>%
  drop_na()%>%
  select(AccSpeciesName,Vessel_diameter_mm)
  
  To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  
  To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  select(Vessel_diameter_mm,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  To_draw=rbind(To_draw1,To_draw2)
  
pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology')],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                # filter(Ks < 250)%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
                   rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
                         Vessel_diameter_mm=Vessel_diameter*0.001
         )%>%
  select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Vessel_diameter_mm)%>%
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
  select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Ks)%>%
  mutate(data_source='Ghana')%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
                           Ks=mean(Ks,na.rm=T),
                           data_source=first(data_source))


 To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  
  To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_RF')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  select(Ks,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
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
##                              SLA                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


pairwise_table<-full_join(Huber,NNNNN_lma[,c('PLOT_TREE',"LMA")],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 mutate(SLA_m2_kg1=1/(LMA/1000))%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
                 filter(log10(SLA_m2_kg1) < 1.7)%>%
                 mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100)%>%
                 rename(AccSpeciesName=Species )%>%
                 select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,SLA_m2_kg1)%>%
                 mutate(data_source='Ghana')%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
                           SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
                           data_source=first(data_source))

 To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  select(SLA_m2_kg1,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  
  To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_FOR')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  select(SLA_m2_kg1,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
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
##                              leaf area (no mechanistic reason)           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pairwise_table<-full_join(Huber,NNNNN_leaf_area[,c('PLOT_TREE',"LeafArea.m2")],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
                         leaf_area_m2=LeafArea.m2
         )%>%
                   rename(AccSpeciesName=Species )%>%
  select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,leaf_area_m2)%>%
  mutate(data_source='Ghana')  %>%
                 group_by(AccSpeciesName)%>%
                 summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
                           leaf_area_m2=mean(leaf_area_m2,na.rm=T),
                           data_source=first(data_source))
  
 To_draw1<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  mutate(data_source='Maurizio')%>%
  select(leaf_area_m2,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  
  To_draw2<-Maurizio_data%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  #filter(biome=='TROP_RF'|biome=='TROP_FOR')%>%
  filter(biome=='TROP_RF')%>%
  mutate(data_source='TROP_FOR_Maurizio')%>%
  select(leaf_area_m2,Sapwood_per_leaf_area_cm2_m2,data_source,AccSpeciesName)
  To_draw=rbind(To_draw1,To_draw2)
  
  
  rbind_to_draw<-rbind(pairwise_table,To_draw)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ leaf_area_m2*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is global,\n red is tropical rainforest species in global dataset, blue is Ghana data')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                             TRY TRAIT DATABASE                           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 Huber value                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scores_na <- function(x, ...) {
  not_na <- !is.na(x)
  scores <- rep(NA, length(x))
  scores[not_na] <- outliers::scores(na.omit(x), ...)
  scores
}


setwd('F:/Oxford/Chapter_one/TRY_global')

global_huber <- fread("16387.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
filter(TraitName=="Stem sapwood cross-sectional area per supported leaf surface area (Huber value)")%>%  
  select(LastName,AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,OrigUnitStr,StdValue,UnitName,TraitName,DataName)%>%
  mutate(OrigValue=as.numeric(OrigValueStr))

##~~~~~~~~~~~~~~~~~~~~~~~~~
##  ~ unit conversion  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~

global_huber%>%
ggplot(aes(log10(OrigValue))) + geom_histogram() + facet_wrap(~LastName)
# this will expose the necessary unit conversion

#........step one correct it to sapwood area per leaf area.......

global_huber<-global_huber%>%
   mutate(
         OrigValue=case_when(
           LastName   =="Cornwell" ~ 1/OrigValue,
           LastName   =="Mencuccini" ~ 1/OrigValue,
           LastName   =="He" ~ 1/OrigValue,
           TRUE ~ OrigValue
         ))

global_huber%>%
ggplot(aes(log10(OrigValue))) + geom_histogram() + facet_wrap(~LastName)

#....................step two use cm2 m-2 unit...................

global_huber<-global_huber%>%
  mutate(UnitName='cm2 m-2',
         StdValue=case_when(
           
           LastName   =="Cornwell" ~ OrigValue*100, #they said mm2/mm2, probably write the unit name wrong?
           LastName   =="Jansen" ~ OrigValue*10000,
           TRUE ~ OrigValue
         ))

global_huber%>%
ggplot(aes((StdValue))) + geom_histogram() + facet_wrap(~LastName)

global_huber%>%
  group_by(LastName)%>%
  summarise(medianhuber=median(StdValue))

Specie_huber_individual<-global_huber%>%
  group_by(ObservationID)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(StdValue,na.rm=T) )

Specie_huber<-global_huber%>%
  group_by(AccSpeciesName)%>%
  summarise(Sapwood_per_leaf_area_cm2_m2 = mean(StdValue,na.rm=T) )




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                specific leaf area                        ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



global_SLA <- fread("16787.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
  filter(TraitName=="Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded")%>%  
  select(LastName,AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,OrigUnitStr,StdValue,UnitName,TraitName)


Specie_SLA_individual<-global_SLA%>%
  group_by(ObservationID)%>%
  summarise(SLA_m2_kg1 = mean(StdValue,na.rm=T)*0.000001/0.000001 ) # from mm2 mg-1 to m2 kg-1


Specie_SLA<-global_SLA%>%
  group_by(AccSpeciesName)%>%
  summarise(SLA_m2_kg1 = mean(StdValue,na.rm=T)*0.000001/0.000001 ) # from mm2 mg-1 to m2 kg-1



To_draw<-Specie_huber%>%
  left_join(Specie_SLA,by='AccSpeciesName')%>%
  drop_na()


  ggplot(To_draw,aes(x=log10(SLA_m2_kg1), y=log10(Sapwood_per_leaf_area_cm2_m2))) + 
  geom_point( color="#69b3a2") +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))
  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               Nmass                                      ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

global_Nmass <- fread("17240.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
  filter(TraitID==14)%>%
    filter(LastName=='Wright')%>%
  select(LastName,AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,OrigUnitStr,StdValue,UnitName,TraitName)
  
Specie_Nmass<-global_Nmass%>%
  group_by(AccSpeciesName)%>%
  summarise(Ngkg = mean(StdValue,na.rm=T)) 

global_Pmass <- fread("17240.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
  filter(TraitID==15)%>%
    filter(LastName=='Wright')%>%
  select(LastName,AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,OrigUnitStr,StdValue,UnitName,TraitName)
  
Specie_Pmass<-global_Pmass%>%
  group_by(AccSpeciesName)%>%
  summarise(Pgkg = mean(StdValue,na.rm=T)) 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                wood density                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Strange enough, wood density could not be combined with huber value on species level with TRY
  
global_WD <- fread("16787.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)%>%
  filter(TraitName=="Stem dry mass per stem fresh volume (stem specific density, SSD, Wood density) Twigs and branches smaller 2cm sapwood")%>%  
  select(LastName,AccSpeciesID,AccSpeciesName,ObservationID,OrigValueStr,OrigUnitStr,StdValue,UnitName,TraitName)

  Specie_WD<-global_WD%>%
  group_by(AccSpeciesName)%>%
  summarise(WD_g_cm3 = mean(StdValue,na.rm=T))

      #....................the above is not working....................
species_vector<-unique(Specie_huber$AccSpeciesName)
median_wd <- BIEN_trait_mean(species = species_vector, trait = "stem wood density")
save(median_wd,file='wood_density.rda')
load('wood_density.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_wd<-median_wd%>%
  rename(AccSpeciesName=species,
         WD_g_cm3=mean_value)%>%
 # filter(level_used=='Species')%>%
  drop_na()
  
  To_draw<-Specie_huber%>%
  left_join(median_wd,by='AccSpeciesName')%>%
  drop_na()
 
ggplot(To_draw,aes(x=log10(WD_g_cm3), y=log10(Sapwood_per_leaf_area_cm2_m2))) + 
  geom_point( color="#69b3a2") +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)


ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ WD_g_cm3, log = "xy", data = To_draw)
summary(ftComSlope)
plot(ftComSlope)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                          Vessel diameter                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#median_wd <- BIEN_trait_mean(species = species_vector, trait = "vessel lumen area")
#save(median_wd,file='vessel_lumen_area.rda')
load('vessel_lumen_area.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_vessel<-median_wd%>%
  rename(AccSpeciesName=species)%>%
  mutate(Vessel_diameter_mm=sqrt(mean_value/pi))%>%
  #filter(level_used=='Species')%>%
  drop_na()%>%
  select(AccSpeciesName,Vessel_diameter_mm)
  
  To_draw<-Specie_huber%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  drop_na()

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Vessel_diameter_mm, log = "xy", data = To_draw)
summary(ftComSlope)
plot(ftComSlope)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  leaf area                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#median_wd <- BIEN_trait_mean(species = species_vector, trait = "leaf area")
#save(median_wd,file='leaf_area_area.rda')
load('leaf_area_area.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_leaf_area<-median_wd%>%
  rename(AccSpeciesName=species)%>%
  mutate(leaf_area_m2=(mean_value/1000000))%>%
  #filter(level_used=='Species')%>%
  drop_na()%>%
  select(AccSpeciesName,leaf_area_m2)

  To_draw<-Specie_huber%>%
  left_join(median_leaf_area,by='AccSpeciesName')%>%
  drop_na()
 ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ leaf_area_m2, log = "xy", data = To_draw)
summary(ftComSlope)
plot(ftComSlope) 


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  try BIEN SLA                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(BIEN)
species_vector<-unique(Specie_huber$AccSpeciesName)

#median_wd <- BIEN_trait_mean(species = species_vector, trait = "leaf area per leaf dry mass")
#save(median_wd,file='leaf area per leaf dry mass.rda')
load('leaf area per leaf dry mass.rda')
median_wd$mean_value<-as.numeric(median_wd$mean_value)
median_wd<-median_wd%>%
  rename(AccSpeciesName=species,
         SLA_m2_kg1=mean_value)%>%
  filter(level_used=='Species')%>%
  drop_na()


To_draw<-Specie_huber%>%
  left_join(median_wd,by='AccSpeciesName')%>%
  mutate( Outlier_id = scores_na((SLA_m2_kg1), type = "iqr"),
          SLA_m2_kg1 = ifelse(abs(Outlier_id) < 0.3, SLA_m2_kg1, NA))%>% 
    drop_na()

  ggplot(To_draw,aes(x=log10(SLA_m2_kg1), y=log10(Sapwood_per_leaf_area_cm2_m2))) + 
  geom_point( color="#69b3a2") +
  theme_ipsum()+
  ggtitle(paste0('BEIN number of species: ', nrow(To_draw)))
  
# Seems worse than TRY
  
  

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                            OVERLAP TWO STUDIES                           ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             huber value and SLA                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



  # Unit for my LMA is g/m2
  # My huber value is mm2/m2 sapwood/leaves
pairwise_table<-full_join(Huber,NNNNN_lma[,c('PLOT_TREE',"LMA")],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 mutate(SLA_m2_kg1=1/(LMA/1000))%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
                 filter(log10(SLA_m2_kg1) < 1.7)%>%
                 mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100)%>%
                 rename(AccSpeciesName=Species )%>%
                 select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,SLA_m2_kg1)%>%
                 mutate(data_source='Ghana')%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
                           SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
                           data_source=first(data_source))
                 
                 

To_draw<-Specie_huber%>%
  left_join(Specie_SLA,by='AccSpeciesName')%>%
  drop_na()%>%
  mutate(data_source='TRY')


rbind_to_draw<-rbind(pairwise_table,To_draw)%>%
  filter(SLA_m2_kg1>2)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ SLA_m2_kg1*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is BIEN vessel, blue is Ghana data')

  ggplot(rbind_to_draw,aes(x=log10(SLA_m2_kg1), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)
summary(lm(data=pairwise_table, log10(Sapwood_per_leaf_area_cm2_m2)~log10(SLA_m2_kg1)))
  
pairwise_table<-full_join(Huber,NNNNN_lma[,c('PLOT_TREE',"LMA")],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 mutate(SLA_m2_kg1=1/(LMA/1000))%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
                 filter(log10(SLA_m2_kg1) < 1.7)%>%
                 mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100)%>%
                 rename(AccSpeciesName=Species )%>%
                 select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,SLA_m2_kg1,SITE)%>%
                 rename(data_source=SITE)%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
                           SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
                           data_source=first(data_source))

rbind_to_draw<-rbind(pairwise_table,To_draw)%>%
  filter(SLA_m2_kg1>2)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ SLA_m2_kg1*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is BIEN vessel, blue is Ghana data')

  ggplot(rbind_to_draw,aes(x=log10(SLA_m2_kg1), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                       huber value and vessel diameter                    ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To_draw<-Specie_huber%>%
  left_join(median_vessel,by='AccSpeciesName')%>%
  drop_na()%>%
  mutate(data_source='TRY')
  
pairwise_table<-full_join(Huber,WD_Ks[,c('PLOT_TREE','Ks','Vessel_diameter','WD','Phenology')],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                # filter(Ks < 250)%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
                   rename(AccSpeciesName=Species )%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
                         Vessel_diameter_mm=Vessel_diameter*0.001
         )%>%
  select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,Vessel_diameter_mm)%>%
  mutate(data_source='Ghana')

            #     group_by(AccSpeciesName)%>%
            #     summarise(Sapwood_per_leaf_area_cm2_m2 = mean(Sapwood_per_leaf_area_cm2_m2,na.rm=T),
            #               Vessel_diameter_mm=mean(Vessel_diameter_mm,na.rm=T),
            #               data_source=first(data_source))


rbind_to_draw<-rbind(pairwise_table,To_draw)


  ggplot(rbind_to_draw,aes(x=log10(Vessel_diameter_mm), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)


ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ Vessel_diameter_mm*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is BIEN vessel, blue is Ghana data')

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                       huber value and leaf area                          ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pairwise_table<-full_join(Huber,NNNNN_leaf_area[,c('PLOT_TREE',"LeafArea.m2")],by='PLOT_TREE')%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  mutate(Sapwood_per_leaf_area_cm2_m2=Huber/100,
                         leaf_area_m2=LeafArea.m2
         )%>%
                   rename(AccSpeciesName=Species )%>%
  select(AccSpeciesName,Sapwood_per_leaf_area_cm2_m2,leaf_area_m2)%>%
  mutate(data_source='Ghana')  

  To_draw<-Specie_huber%>%
  left_join(median_leaf_area,by='AccSpeciesName')%>%
  drop_na()%>%
  mutate(data_source='TRY')

rbind_to_draw<-rbind(pairwise_table,To_draw)

ftComSlope = sma(Sapwood_per_leaf_area_cm2_m2 ~ leaf_area_m2*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is BIEN vessel, blue is Ghana data') 

  ggplot(rbind_to_draw,aes(x=log10(leaf_area_m2), y=log10(Sapwood_per_leaf_area_cm2_m2),color=data_source)) + 
  geom_point( ) +
  theme_ipsum()+
  ggtitle(paste0('TRY number of species: ', nrow(To_draw)))+
  geom_smooth(method=lm,se=F, fullrange=TRUE)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                       SLA and NgkG                                       ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

pairwise_table<- NNNNN%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 mutate(SLA=1/(LMA/1000))%>%
                 filter(log10(SLA) < 1.7)%>%
               
                 rename(AccSpeciesName=Species,
                        SLA_m2_kg1=SLA)%>%
                 select(AccSpeciesName,Ngkg,SLA_m2_kg1,SITE)%>%
                 rename(data_source=SITE)%>%
                 mutate(data_source='Ghana')%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Ngkg = mean(Ngkg,na.rm=T),
                           SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
                           data_source=first(data_source))

 To_draw<-Specie_Nmass%>%
  left_join(Specie_SLA,by='AccSpeciesName')%>%
  drop_na()%>%
  mutate(data_source='TRY')

rbind_to_draw<-rbind(pairwise_table,To_draw)


ftComSlope = sma(SLA_m2_kg1 ~ Ngkg*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is TRY global, blue is Ghana data') 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                       SLA and PgkG                                       ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 

pairwise_table<- NNNNN%>%
                 mutate(SITE=substring(PLOT_TREE,1,3))%>%
                 mutate(SLA=1/(LMA/1000))%>%
                 filter(log10(SLA) < 1.7)%>%
                 filter((SLA) > 2.5)%>%
                 rename(AccSpeciesName=Species,
                        SLA_m2_kg1=SLA)%>%
                 select(AccSpeciesName,Pgkg,SLA_m2_kg1,SITE)%>%
                 rename(data_source=SITE)%>%
                 mutate(data_source='Ghana')%>%
                 group_by(AccSpeciesName)%>%
                 summarise(Pgkg = mean(Pgkg,na.rm=T),
                           SLA_m2_kg1=mean(SLA_m2_kg1,na.rm=T),
                           data_source=first(data_source))

 To_draw<-Specie_Pmass%>%
  left_join(Specie_SLA,by='AccSpeciesName')%>%
  drop_na()%>%
  mutate(data_source='TRY')

rbind_to_draw<-rbind(To_draw,pairwise_table)


ftComSlope = sma(SLA_m2_kg1 ~ Pgkg*data_source, log = "xy", data = rbind_to_draw)
summary(ftComSlope)
plot(ftComSlope,main='yellow is TRY global, blue is Ghana data') 
