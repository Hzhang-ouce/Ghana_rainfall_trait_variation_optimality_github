
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##............................READ.IN.TRAIT.VALUE...........................----
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())
library(tidyverse)
library(Hmisc)
standard_error_calc <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

weighted.var.se <- function(x, weights, na.rm = TRUE){
  var <- Hmisc::wtd.var(x, weights, na.rm)
  weights <- sum( (weights / sum(weights))^2 )
  sqrt(var*weights)
}


scores_na <- function(x, ...) {
  not_na <- !is.na(x)
  scores <- rep(NA, length(x))
  scores[not_na] <- outliers::scores(na.omit(x), ...)
  scores
} # This functions is to add outlier id 
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            ~~
##                          PREPARE DATA FROM GHANA                         ----
##                                                                            ~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




  setwd("H:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data")



##---------------------------------------------------------------
##                          phenology                           -
##---------------------------------------------------------------

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

WD_Ks<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  dplyr::rename(Ks=!!as.name("theoretical.conductivity.poorter2010.kg.m-1.s-1.Mpa-1"))%>%
  dplyr::rename(Vessel_diameter = !!as.name('hydraulically.weighted.diam.micron'))


#the left join will bring species to the twig density table
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


Huber<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('LeafArea/SapwoodArea_m2/mm2'))%>%
  mutate(Huber=1/Huber)


Huber<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('SapwoodArea/LeafArea_mm2/cm2'))%>%
  mutate(Sapwood_to_LeafArea_mm2_m2=Huber*10000)%>%# get into mm2/m2 then in the following calculation, it will change into cm2/m2 when /100
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  mutate(Outlier_id=scores_na((Sapwood_to_LeafArea_mm2_m2),type="iqr"))%>%
  mutate(Sapwood_to_LeafArea_mm2_m2=ifelse(abs(Outlier_id) < 1.2, Sapwood_to_LeafArea_mm2_m2, NA)) 
##----------------------------------------------------------------
##                        tugor loss point                       -
##----------------------------------------------------------------

List_of_site<-c('KOG','BOB','ANK')

TLPoint<-openxlsx::read.xlsx("Overall_trait_raw_data.xlsx",sheet='turgor')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM","MONTH","YEAR"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))

TLPoint[87,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake
TLPoint[96,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake
TLPoint<-TLPoint%>%
  rename(Turgor=!!as.name("Psi_TLP.(MPa)"))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")


##----------------------------------------------------------------
##                          Twig density                         -
##----------------------------------------------------------------
List_of_site<-c('KOG','BOB','ANK')

#read data
Twig_density_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")

##----------------------------------------------------------------
##                          Vcmax                         -
##----------------------------------------------------------------

List_of_site<-c('KOG','BOB','ANK')

Vcmax_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='VCmax')%>%
  separate("File_Name", sep='-',into=c("PLOT","TREE","STEM","TIME"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  left_join(Census_big_table[c('Species','PLOT_TREE')],by="PLOT_TREE")%>%
  mutate(Outlier_id=scores_na(log2(Vcmax_ACi_25C),type="iqr"))%>%
  mutate(Vcmax_ACi_25C=ifelse(abs(Outlier_id) < 1.5, Vcmax_ACi_25C, NA))
##----------------------------------------------------------------
##                          Npercent                         -
##----------------------------------------------------------------


NNNNN=openxlsx::read.xlsx('F:/Oxford/Chapter_one/From_imma/Traits/Overall_trait_raw_data.xlsx',sheet='Npercent')%>%
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
  mutate(LeafDrymass.g=ifelse(abs(Outlier_id) < 1.5, LeafDrymass.g, NA))

hist(log2(NNNNN_LeafDrymass$LeafDrymass.g))

NNNNN_ngkg<- NNNNN %>% 
  drop_na(Species)%>%
  mutate(Site=substr(Plot_code,0,3))%>%
  select(Site,Species,Pgkg,Ngkg,Plot_code,Leafcode,PLOT_TREE)%>%
  unique()%>%
  mutate(Outlier_id=scores_na((Ngkg),type="iqr"))%>%
  mutate(Ngkg=ifelse(abs(Outlier_id) < 1.5, Ngkg, NA))%>%
   mutate(Outlier_id=scores_na((Pgkg),type="iqr"))%>%
  mutate(Pgkg=ifelse(abs(Outlier_id) < 1.5, Pgkg, NA))
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


C13_table=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='C13')%>%
  rename(C13=!!as.name("13CVPDB"))%>%
  separate("ID", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  mutate(ci_ca_ratio=get_chi_from_13C(C13)[[3]])%>%
  left_join(Census_big_table[,c('PLOT_TREE','Species')],by="PLOT_TREE")%>%
  left_join(pheno[,c('Species','Phenology')],by="Species")

##----------------------------------------------------------------
##                             Amax                             -
##----------------------------------------------------------------


Amax_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Amax')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  filter(Light=="Sun")%>%
  left_join(Census_big_table[,c('PLOT_TREE','Species')],by="PLOT_TREE")


##----------------------------------------------------------------
##                             Tree height table
##----------------------------------------------------------------


Files_to_read_census<-list.files(path = "/Census_data", pattern = '*cleaned*', full.names = T) #Input folder

kog_04_tree_height<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data/Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='Savanna 1')[,c(22,9)]%>%
  mutate(plot_code='KOG_04')
kog_01_tree_height<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data/Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='remnant forest 1 ')[,c(20,9)]%>%
  mutate(plot_code='KOG_01')
kog_03_tree_height<-openxlsx::read.xlsx('F:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data/Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='transition')[,c(13,4)]%>%
  mutate(plot_code='KOG_03')

colnames(kog_01_tree_height)<-c('Height','Species','plot_code')
 colnames(kog_03_tree_height)<-c('Height','Species','plot_code')
colnames(kog_04_tree_height)<-c('Height','Species','plot_code')


kog_tree_height<-rbind(kog_04_tree_height,kog_03_tree_height)%>%
  mutate(Height=as.numeric(Height))%>%
    drop_na()

height_table <- lapply(Files_to_read_census, read.csv) %>%
  data.table::rbindlist()%>%
  filter(plot_code=='ANK_01' |plot_code=='ANK_03' |plot_code=='BOB_01' |plot_code=='BOB_02'|plot_code=='KOG_02')%>%
  mutate(year=substr(census_date,1,4))%>%
  filter(year==2012)%>%
  mutate(census_year_plot=paste0(str_sub(census_date,1,4),plot_code))%>%
  mutate(Height=ifelse(plot_code=='BOB_01',exp( 0.6152+ 0.6362*log(DBH4/10)+0.0109*25813371/1000000),Height),
         Height=ifelse(plot_code=='BOB_02',exp( 0.6152+ 0.6362*log(DBH4/10)+0.0109*32814313/1000000),Height))%>%
  mutate(Height=ifelse(plot_code=='KOG_02',exp( 0.4619+ 0.6362*log(DBH4/10)+0.0109*21446610/1000000),Height))%>%
  filter(Height>0)%>%
 #   mutate(Height=ifelse(plot_code=='BOB_01',exp( 0.8946+ 0.6365*log(DBH4/10)),Height),
 #        Height=ifelse(plot_code=='BOB_02',exp( 0.8946+ 0.6365*log(DBH4/10)),Height))%>%filter(Height>0)%>%
  drop_na(Height)%>%
  mutate(plot_code = str_replace_all(plot_code,'KOG_02','KOG_01'))
  
height_table<-rbind(height_table,kog_tree_height,fill=TRUE)%>%
  dplyr::select(Species,plot_code,Height)%>%
  rename(PLOT=plot_code)
  


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##................................FILL.IN.TRAIT.............................----
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    huber                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

full_census<-openxlsx::read.xlsx('Census_data/Full_census_Ghana.xlsx')%>%
  group_by(TreeID)%>%
  summarise(PLOT=first(Plot.Code),
            TREE=first(Tag.Number),
            Species=first(Species),
            DBH=mean(D2,na.rm=T))%>%
  mutate(basal= pi*(DBH/2)^2)

Huber<-Huber%>%
  group_by(Species)%>%
  summarise(Sapwood_to_LeafArea_mm2_m2=mean(Sapwood_to_LeafArea_mm2_m2,na.rm=T))

full_census_Huber<-full_census%>%
  left_join(Huber,by='Species')%>%
  group_by(PLOT)%>%
  summarise(Sapwood_to_LeafArea_mm2_m2_ste=weighted.var.se(Sapwood_to_LeafArea_mm2_m2,basal,na.rm=T),
            AS_AL_cm2_m2_mean=wtd.mean(Sapwood_to_LeafArea_mm2_m2,basal,na.rm=T))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                             Ks vessel diameter                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





WD_Ks_species<-WD_Ks%>%
  group_by(Species)%>%
  summarise(Ks_ste=standard_error_calc(Ks,na.rm=T),
            Ks=mean(Ks,na.rm=T),
            Vessel_diameter_ste=standard_error_calc(Vessel_diameter,na.rm=T),
            Vessel_diameter=mean(Vessel_diameter,na.rm=T),
           WD_ste=standard_error_calc(WD,na.rm=T),
            WD=mean(WD,na.rm=T),
           numberoftree=n())
full_census_WD_KS<-full_census%>%
  left_join(WD_Ks_species,by='Species')%>%
  group_by(PLOT)%>%
  summarise(Kp_ste=weighted.var.se(Ks,basal,na.rm=T),
            Kp_mean=wtd.mean(Ks,basal,na.rm=T),
            Vessel_diameter_ste=weighted.var.se(Vessel_diameter,basal,na.rm=T),
            Vessel_diameter_mean=wtd.mean(Vessel_diameter,basal,na.rm=T),
            WD_ste=weighted.var.se(WD,basal,na.rm=T),
            Wood_density_mean=wtd.mean(WD,basal,na.rm=T))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Amax                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Amax_table<-Amax_table%>%
  group_by(Species)%>%
  summarise(Asat=mean(Asat,na.rm=T),
            Rd=mean(Rd,na.rm=T),
            Amax=mean(Amax,na.rm=T))
full_census_Amax<-full_census%>%
  left_join(Amax_table,by='Species')%>%
  group_by(PLOT)%>%
  summarise(Asat_ste=weighted.var.se(Asat,basal,na.rm=T),
            Asat_mean=wtd.mean(Asat,basal,na.rm=T),
            Rd_ste=weighted.var.se(Rd,basal,na.rm=T),
            Rd_mean=wtd.mean(Rd,basal,na.rm=T),
            Amax_ste=weighted.var.se(Amax,basal,na.rm=T),
            Amax_mean=wtd.mean(Amax,basal,na.rm=T))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              turgor loss point                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TLPoint<-TLPoint%>%
  group_by(Species)%>%
  summarise(Turgor=mean(Turgor,na.rm=T))

full_census_Turgor<-full_census%>%
  left_join(TLPoint[,c('Species','Turgor')],by='Species')%>%
  group_by(PLOT)%>%
  summarise(TLP_ste=weighted.var.se(Turgor,basal,na.rm=T),
            TLP_mean=wtd.mean(-Turgor,basal,na.rm=T))
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 leaf trait                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NNNNN_LeafDrymass<-NNNNN_LeafDrymass%>%
  group_by(Species)%>%
  summarise(LeafDrymass.g=mean(LeafDrymass.g,na.rm=T))

NNNNN_leaf_area<-NNNNN_leaf_area%>%
  group_by(Species)%>%
  summarise(LeafArea.m2=mean(LeafArea.m2,na.rm=T))

NNNNN_lma<-NNNNN_lma%>%
  group_by(Species)%>%
  summarise(LMA=mean(LMA,na.rm=T))

NNNNN_ngkg<-NNNNN_ngkg%>%
  group_by(Species)%>%
  summarise(Ngkg=mean(Ngkg,na.rm=T),
            Pgkg=mean(Pgkg,na.rm=T))

full_census_drymass<-full_census%>%
  left_join(NNNNN_LeafDrymass[,c('Species','LeafDrymass.g')],by='Species')%>%
  group_by(PLOT)%>%
  summarise(LeafDrymass_ste=weighted.var.se(LeafDrymass.g,basal,na.rm=T),
            LeafDrymass_mean=wtd.mean(LeafDrymass.g,basal,na.rm=T))

full_census_area<-full_census%>%
  left_join(NNNNN_leaf_area[,c('Species','LeafArea.m2')],by='Species')%>%
  group_by(PLOT)%>%
  summarise(LeafArea_ste=weighted.var.se(LeafArea.m2,basal,na.rm=T),
            LeafArea_mean=wtd.mean(LeafArea.m2,basal,na.rm=T))

full_census_lma<-full_census%>%
  left_join(NNNNN_lma[,c('Species','LMA')],by='Species')%>%
  group_by(PLOT)%>%
  summarise(LMA_ste=weighted.var.se(LMA,basal,na.rm=T),
            LMA_mean=wtd.mean(LMA,basal,na.rm=T))

full_census_ngkg<-full_census%>%
  left_join(NNNNN_ngkg,by='Species')%>%
  group_by(PLOT)%>%
  summarise(Ngkg_ste=weighted.var.se(Ngkg,basal,na.rm=T),
            Ngkg_mean=wtd.mean(Ngkg,basal,na.rm=T),
            Pgkg_ste=weighted.var.se(Pgkg,basal,na.rm=T),
            Pgkg_mean=wtd.mean(Pgkg,basal,na.rm=T))


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Vcmax                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Vcmax_table<-Vcmax_table%>%
  group_by(Species)%>%
  summarise(Vcmax_ACi_25C=mean(Vcmax_ACi_25C,na.rm=T),
            Jmax_ACi_25C=mean(Jmax_ACi_25C,na.rm=T))

full_census_vcmax<-full_census%>%
  left_join(Vcmax_table,by='Species')%>%
  group_by(PLOT)%>%
  summarise(Vcmax25_ste=weighted.var.se(Vcmax_ACi_25C,basal,na.rm=T),
            Vcmax25_mean=wtd.mean(Vcmax_ACi_25C,basal,na.rm=T),
            Jmax25_ste=weighted.var.se(Jmax_ACi_25C,basal,na.rm=T),
            Jmax25_mean=wtd.mean(Jmax_ACi_25C,basal,na.rm=T)            )
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Cica ratio                            ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


C13_table<-C13_table%>%
  group_by(Species)%>%
  summarise(ci_ca_ratio=mean(ci_ca_ratio,na.rm=T))

full_census_ci_ca_ratio<-full_census%>%
  left_join(C13_table,by='Species')%>%
  group_by(PLOT)%>%
  summarise(ci_ca_ratio_ste=weighted.var.se(ci_ca_ratio,basal,na.rm=T),
            ci_ca_ratio_mean=wtd.mean(ci_ca_ratio,basal,na.rm=T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                               tree height                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


height_table<-height_table%>%
  group_by(Species)%>%
  summarise(Height=max(Height,na.rm=T))

full_census_height<-full_census%>%
  left_join(height_table,by='Species')%>%
  group_by(PLOT)%>%
  summarise(height_ste=weighted.var.se(Height,basal,na.rm=T),
            height_mean=wtd.mean(Height,basal,na.rm=T))

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                draw a plot                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

grab_moisture<-openxlsx::read.xlsx("soil_moisture_all_plot.xlsx")
colnames(grab_moisture)<-c('PLOT','volumetric_water_content_mean','volumetric_water_content_ste')

joined_table<-full_census_WD_KS%>%
 # left_join(full_census_drymass,by='PLOT')%>%
 # left_join(full_census_area,by='PLOT')%>%
    left_join(full_census_Huber,by='PLOT')%>%
  #left_join(full_census_lma,by='PLOT')%>%
  left_join(full_census_vcmax,by='PLOT')%>%
  left_join(full_census_Turgor,by='PLOT')%>%
  left_join(full_census_Amax,by='PLOT')%>%
  #left_join(full_census_ngkg,by='PLOT')%>%
   left_join(full_census_ci_ca_ratio,by='PLOT')%>%
  left_join(grab_moisture,by='PLOT')
  #left_join(full_census_height,by='PLOT')

final_table_mean<-joined_table%>%
  dplyr::select(PLOT,ends_with('mean'))%>%
  pivot_longer(cols = !PLOT,
               names_pattern = "(.*)_mean",
               values_to = 'mean',
               names_to = 'Traits' )

final_table_ste<-joined_table%>%
  dplyr::select(ends_with('ste'))%>%
  pivot_longer(cols = everything(),
               names_pattern = "(.*)_ste",
               values_to = 'ste',
               names_to = 'Traits' )
final_table<-bind_cols(final_table_mean,final_table_ste)

final_table$PLOT <- factor(final_table$PLOT, levels = c("ANK-03", "BOB-03", "BOB-04" ,"BOB-05" ,"BOB-06", "BOB-02" ,"ANK-02" ,"BOB-01",
                                                                      "ANK-01", "KOG-02", "KOG-03", "KOG-05", "KOG-04", "KOG-06"))
final_table$PLOT <- factor(final_table$PLOT, levels = c("ANK-03","ANK-02","ANK-01","BOB-03", "BOB-04" ,"BOB-05" ,"BOB-06", "BOB-02" ,"BOB-01",
                                                                       "KOG-02", "KOG-03", "KOG-04", "KOG-05", "KOG-06"))

ggplot(final_table) +
  geom_bar( aes(x=PLOT, y=mean), stat="identity", fill="grey", alpha=0.7,width = 0.8) +
  geom_errorbar(aes(x=PLOT, y=mean, ymin=mean-ste, ymax=mean+ste),alpha=0.9, size=1.3,width=0.4)+ theme_bw()+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  facet_wrap(~Traits...2, scales = "free",ncol = 3)
ggsave(file='H:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/Ghana_rainfall_trait_variation_optimality/trait_extrapolation/Trait_extrapolation.jpg', width=9.12, height=12.51)
