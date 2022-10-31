
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##---------- COMMUNITY WEIGHTED MEAN IN GHANA PLOT TO PLOT COMPARISON-----------
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#code following F:/Oxford/Chapter_one/From_imma/Traits/Community Weighted Mean Calculation.pdf
#code following http://rstudio-pubs-static.s3.amazonaws.com/502799_bf78c38d5f0a49a39207a0b0039fd5f5.html

# For ANOVA, following https://statdoe.com/step-by-step-barplots-for-one-factor-in-r/
# 
rm(list=ls())
library(multcompView)
library(tidyverse)
library(Hmisc)
library(multcomp)
library(vegan)
library(ARTofR)
setwd("H:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/input_data")
#setwd("~/GitHub/Ghana_rainfall_trait_variation_optimality_github/input_data")

scores_na <- function(x, ...) {
  not_na <- !is.na(x)
  scores <- rep(NA, length(x))
  scores[not_na] <- outliers::scores(na.omit(x), ...)
  scores
}

#The main reference is this paper, by Donald F. Gatz and Luther Smith, where 3
#formula based estimators are compared with bootstrap results. The best
#approximation to the bootstrap result comes from Cochran (1977):
# https://stats.stackexchange.com/questions/25895/computing-standard-error-in-weighted-mean-estimation
# the above method is wrong, now replaced with:
# https://githubmemory.com/repo/harrelfe/Hmisc/issues/138
# https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Statistical_properties

weighted.var.se <- function(x, weights, na.rm = TRUE){
  var <- Hmisc::wtd.var(x, weights, na.rm)
  weights <- sum( (weights / sum(weights))^2 )
  
  sqrt(var*weights)
}

standard_error_calc <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')
species_basal<-read.csv('Census_data/Census_table_all_plots_Species_list_basal_area.csv')


#...............................................................................
#                                                                              .
#  Common work flow:                                                           .
#  1. Read in trait information and match up with Cesus data (including        .
#     species and basal area)                                                  .
#  2. Check normal distribution, transformed if necessary                      .
#  3. Remove extreme outlier with IQR method                                   .
#  4. ANOVA on leaf (or branch) level and obtain letters for significance at   .
#     0.05                                                                     .
#  5. Community weighted mean and weighted standard error                      .
#                                                                              .
#...............................................................................




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                twig density                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#read data
Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))


Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
                  filter(grepl('^[0-9]{1,}$', Tag.No))%>%
                  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
                  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))


#the left join will bring species to the trait table
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  drop_na(.,any_of('TwigDensity.g.cm3'))%>%
  mutate(Outlier_id=outliers::scores((TwigDensity.g.cm3),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(TwigDensity.g.cm3)))

# Check normal distribution (twig density)
plot_histogram=qplot(TwigDensity.g.cm3, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="TwigDensity.g.cm3_hist.jpg",width = 3, height = 3)
#Twig density is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$TwigDensity.g.cm3, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/TwigDensity.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Twig density (g cm-3)")
dev.off()


Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>% 
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(TwigDensity.g.cm3 ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(TwigDensity.g.cm3,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(TwigDensity.g.cm3,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/TwigDensity.g.cm3_cwd.xlsx')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Wood density                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#read data
Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Twig_density')%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))
species_basal<-read.csv('Census_data/Census_table_all_plots_Species_list_basal_area.csv')


Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))


#the left join will bring species to the trait table
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  drop_na(.,any_of('WD'))%>%
  mutate(Outlier_id=outliers::scores((WD),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(WD)))

# Check normality (wood density)
plot_hist=qplot(WD, data=Trait_big_table, geom="histogram") 
plot_hist
ggsave(plot_hist,file="wood_density_hist.jpg",width = 3, height = 3)
# WD is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$WD, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/WD.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Wood density (g cm-3)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(WD ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(WD,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(WD,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/WD_cwd.xlsx')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        sapwood specific conductivity                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  dplyr::rename(Ks=!!as.name("theoretical.conductivity.poorter2010.kg.m-1.s-1.Mpa-1"))%>%
  dplyr::rename(Vessel_diameter = !!as.name('hydraulically.weighted.diam.micron'))%>%
  mutate(Outlier_id=scores_na(log2(Ks),type="iqr"))%>%
  mutate(Ks=ifelse(abs(Outlier_id) < 1.5, Ks, NA))
#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table%>%drop_na(any_of('Species'))%>%
  drop_na(any_of('Ks'))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Ks)))

# Check normal distribution
plot_histogram=qplot(Ks, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Ks_hist.jpg",width = 3, height = 3)
# The data is right-skewed, log transform it
plot_histogram=qplot(log10(Ks), data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Ks_log10_hist.jpg",width = 3, height = 3)
Trait_big_table$Ks_t <- log10(Trait_big_table$Ks)

#Variance partitioning
varp <- varpart (Trait_big_table$Ks_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Ks.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Kp (kg m-1 MPa-1 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Ks_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Ks_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Ks_t,weightsss,na.rm=T))%>%
  mutate(mean_original = 10^mean_value_trandformed,
         upper = 10^(mean_value_trandformed+ standard_d_trandformed),
         lower = 10^(mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Ks_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Narea                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



NNNNN=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Npercent')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE","Branch"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  mutate(sun_shade=str_detect(Branch, ".*SH"))%>%
  filter(!(sun_shade))%>%
  mutate(sun_shade='Sun')

List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-NNNNN

#the left join will bring species to the twig density table
Trait_big_table_species<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table_species%>%
  drop_na(.,any_of('Species'))%>%
  drop_na(.,any_of('Narea_g_m2'))%>%
  mutate(Outlier_id=outliers::scores((Narea_g_m2),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Narea_g_m2)))

# Check normal distribution
plot_histogram=qplot(Narea_g_m2, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Narea_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Narea_g_m2, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Narea.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Narea (g m-2)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Narea_g_m2 ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Narea_g_m2,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Narea_g_m2,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Narea_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Ngkg                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')
Trait_big_table<-NNNNN

#the left join will bring species to the twig density table
Trait_big_table_species<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table_species%>%
  drop_na(.,any_of('Species'))%>%
  drop_na(.,any_of('Ngkg'))%>%
  mutate(Outlier_id=outliers::scores(sqrt(Ngkg),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Ngkg)))

# Check normal distribution
plot_histogram=qplot(Ngkg, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Ngkg_hist.jpg",width = 3, height = 3)
# The data is right-skewed, square root transform it
plot_histogram=qplot(sqrt(Ngkg), data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Ngkg_sqrt_hist.jpg",width = 3, height = 3)
Trait_big_table$Ngkg_t <- sqrt(Trait_big_table$Ngkg)

#Variance partitioning
varp <- varpart (Trait_big_table$Ngkg_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Ngkg.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Nmass (g/kg)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Ngkg_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Ngkg_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Ngkg_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Ngkg_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Pgkg                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')
Trait_big_table<-NNNNN

#the left join will bring species to the twig density table
Trait_big_table_species<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table_species%>%
  drop_na(.,any_of('Species'))%>%
  drop_na(.,any_of('Pgkg'))%>%
  mutate(Outlier_id=outliers::scores((Pgkg),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Pgkg)))

# Check normal distribution
plot_histogram=qplot(Pgkg, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Pgkg_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Pgkg, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Pgkg.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Pmass (g/kg)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Pgkg ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Pgkg,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Pgkg,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Pgkg_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                     LMA                                  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')
Trait_big_table<-NNNNN

#the left join will bring species to the twig density table
Trait_big_table_species<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table_species%>%
  drop_na(.,any_of('Species'))%>%
  drop_na(.,any_of('LMA'))%>%
  mutate(Outlier_id=outliers::scores((LMA),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(LMA)))

# Check normal distribution
plot_histogram=qplot(LMA, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="LMA_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$LMA, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/LMA.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "LMA (g m-2)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(LMA ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(LMA,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(LMA,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/LMA_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Vcmax_ACi_25C                             ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='VCmax')%>%
  separate("File_Name", sep='-',into=c("PLOT","TREE","STEM","TIME"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))


#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  drop_na(.,any_of('Vcmax_ACi_25C'))%>%
  mutate(Outlier_id=outliers::scores((Vcmax_ACi_25C),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Vcmax_ACi_25C)))

# Check normal distribution
plot_histogram=qplot(Vcmax_ACi_25C, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Vcmax_ACi_25C_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Vcmax_ACi_25C, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Vcmax_ACi_25C.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Vcmax25(umol CO2 m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Vcmax_ACi_25C ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Vcmax_ACi_25C,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Vcmax_ACi_25C,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Vcmax_ACi_25C_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                Jmax_ACi_25C                              ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='VCmax')%>%
  separate("File_Name", sep='-',into=c("PLOT","TREE","STEM","TIME"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))
species_basal<-read.csv('Census_data/Census_table_all_plots_Species_list_basal_area.csv')

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  drop_na(.,any_of('Jmax_ACi_25C'))%>%
  mutate(Outlier_id=outliers::scores(sqrt(Jmax_ACi_25C),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Jmax_ACi_25C)))

# Check normal distribution
plot_histogram=qplot(Jmax_ACi_25C, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Jmax_ACi_25C_hist.jpg",width = 3, height = 3)
# The data is slightly right-skewed, square root it
Trait_big_table$Jmax_ACi_25C_t <- sqrt(Trait_big_table$Jmax_ACi_25C)
plot_histogram=qplot(Jmax_ACi_25C_t, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Jmax_ACi_25C_hist_sqrt.jpg",width = 3, height = 3)

#Variance partitioning
varp <- varpart (Trait_big_table$Jmax_ACi_25C_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Jmax_ACi_25C.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Jmax25 (umol CO2 m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Jmax_ACi_25C_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Jmax_ACi_25C_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Jmax_ACi_25C_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Jmax_ACi_25C_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Turgor loss point                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx("Overall_trait_raw_data.xlsx",sheet='turgor')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM","MONTH","YEAR"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))
colnames(Trait_big_table)[7]<-"Turgor"
Trait_big_table[87,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake
Trait_big_table[96,'PLOT']<-'BOB02' #It was BOB01 but looks like a typo mistake

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  mutate(Turgor=log10(-Turgor))%>%
  drop_na(Turgor)%>%
  mutate(Outlier_id=outliers::scores((Turgor),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Turgor)))

# Check normal distribution
plot_histogram=qplot(Turgor, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Turgor_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Turgor, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Turgor.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "TLP (MPa)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Turgor ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Turgor,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Turgor,weightsss,na.rm=T))%>%
  mutate(mean_original = -10^mean_value_trandformed,
         upper = -10^(mean_value_trandformed-standard_d_trandformed),
         lower = -10^(mean_value_trandformed+standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Turgor_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                              Huber value AS/AL                           ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('SapwoodArea/LeafArea_mm2/cm2'))%>%
  mutate(Huber=Huber*100) # change into cm2/m2 when *100



#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  mutate(Outlier_id=scores_na(log10(Huber),type="iqr"))%>%
  mutate(Huber=ifelse(abs(Outlier_id) < 1.5, Huber, NA))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))
Trait_big_table<-drop_na(Trait_big_table,any_of('Huber'))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Huber)))

# Check normal distribution
plot_histogram=qplot(Huber, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Huber_hist.jpg",width = 3, height = 3)
# The data is right-skewed, log transform it
Trait_big_table$Huber_t <- log10(Trait_big_table$Huber)
plot_histogram=qplot(Huber_t, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Huber_log10_hist.jpg",width = 3, height = 3)

#Variance partitioning
varp <- varpart (Trait_big_table$Huber_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Huber.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "AS/AL(cm-2 m-2)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Huber_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Huber_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Huber_t,weightsss,na.rm=T))%>%
  mutate(mean_original = 10^mean_value_trandformed,
         upper = 10^(mean_value_trandformed+ standard_d_trandformed),
         lower = 10^(mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Huber_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                 A2000 Amax                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Amax')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  filter(Light=="Sun")

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD','basal')],by="PLOT_TREE")%>%
  mutate(Outlier_id=outliers::scores((Amax),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Amax)))

# Check normal distribution
plot_histogram=qplot(Amax, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Amax_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Amax, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Amax.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Amax (umol CO2 m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Amax ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Amax,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Amax,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Amax_cwd.xlsx')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  A400 Asat                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Amax')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  filter(Light=="Sun")

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD','basal')],by="PLOT_TREE")%>%
  mutate(Outlier_id=outliers::scores(sqrt(Asat),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Asat)))

# Check normal distribution
plot_histogram=qplot(Asat, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Asat_hist.jpg",width = 3, height = 3)
# The data is right-skewed.
Trait_big_table$Asat_t <- sqrt(Trait_big_table$Asat)
plot_histogram=qplot(Asat_t, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Asat_sqrt_hist.jpg",width = 3, height = 3)

#Variance partitioning
varp <- varpart (Trait_big_table$Asat_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Asat.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Asat (umol CO2 m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Asat_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Asat_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Asat_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Asat_cwd.xlsx')




##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                  Dark resp                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Amax')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  filter(Light=="Sun")

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD','basal')],by="PLOT_TREE")%>%
  drop_na(any_of('DResp'))%>%
  mutate(Outlier_id=outliers::scores((DResp),type="iqr"),
         DResp=-DResp)%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(DResp)))

# Check normal distribution
plot_histogram=qplot(DResp, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="DResp_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$DResp, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/DResp.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Rd (umol CO2 m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(DResp ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(DResp,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(DResp,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/DResp_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    ci/ca                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')


#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))


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


Trait_big_table=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='C13')%>%
  dplyr::rename(C13=!!as.name(colnames(.)[2]))%>%
  separate("ID", sep='-',into=c("PLOT","TREE"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  mutate(Chi_chi=get_chi_from_13C(C13)[[3]])%>%
  left_join(Census_big_table[,c('PLOT_TREE','Species')],by="PLOT_TREE")

Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))%>%
  mutate(Outlier_id=outliers::scores((Chi_chi),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Chi_chi)))

# Check normal distribution
plot_histogram=qplot(Chi_chi, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="chi_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Chi_chi, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/chi.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "ci/ca (unitless)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Chi_chi ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Chi_chi,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Chi_chi,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/chi_cwd.xlsx')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        tree height                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read in datasheets
ANK01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Ankasa_Ank1')[,c(3,13)]%>%
  drop_na()%>%
  mutate(PLOT='ANK01')
ANK03<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Ankasa_Ank3')[,c(3,13)]%>%
  drop_na()%>%
  mutate(PLOT='ANK03')
BOB01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Bobiri_SNR')[,c(4,14)]%>%
  drop_na()%>%
  mutate(PLOT='BOB01')%>%
  rename(Tag.No=TreeTag)
BOB02<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Bobiri_Y53')[,c(5,15)]%>%
  drop_na()%>%
  mutate(PLOT='BOB02')%>%
  rename(Tag.No=TreeTag)
KOG01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_F1')[,c(2,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG01')%>%
  rename(Tag.No=Tag)
KOG03<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_FST')[,c(2,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG03')%>%
  rename(Tag.No=New.Tag)
KOG04<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_Sav1')[,c(2,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG04')%>%
  rename(Tag.No=Tag)

Tree_height_table <- rbind(ANK01, ANK03, BOB01, BOB02, KOG01, KOG03, KOG04)

#the left join will bring species to the table
Census_big_table2<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  #filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  rename(PLOT=plot_code)
Trait_big_table<-left_join(Tree_height_table,Census_big_table2[c('Tag.No','Species','PLOT')],by=c("PLOT","Tag.No"))%>%
  drop_na(any_of('Species'))%>%
  mutate(Outlier_id=outliers::scores((Htotal_m),type="iqr"))%>%
  filter(Outlier_id<2)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Htotal_m)))

# Check normal distribution
plot_histogram=qplot(Htotal_m, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="TreeHeight_hist.jpg",width = 3, height = 3)
# The data is normally distributed

#Variance partitioning
varp <- varpart (Trait_big_table$Htotal_m, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/TreeHeight.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green')) +
  title(main = "Hmax (m)")
dev.off()

Trait_big_table2_tree_height<-Trait_big_table%>%
  # group_by(PLOT)%>%
  # mutate(rank=ecdf(Htotal_m)(Htotal_m))%>%
  # filter(rank<1.1 & rank>0.8)
 group_by(PLOT,Species)%>%
  summarise(Htotal_m=max(Htotal_m,na.rm=T),
           Species_PLOT=first(Species_PLOT))

Trait_big_table2<-Trait_big_table2_tree_height%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table2_tree_height$Species_PLOT)==Species_PLOT)))%>% 
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))

# lm for significant test
mod <- lm(Htotal_m ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)
Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')
Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Htotal_m,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Htotal_m,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed,
         upper = (mean_value_trandformed+ standard_d_trandformed),
         lower = (mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/TreeHeight_cwd.xlsx')


#...............................................................................
#                                                                              .
#  important note: BOB02 has some extremely tall tree, if we calculate max     .
#  tree height, it is overwelming, we ignore those several trees, then we      .
#  find ANK is highest. But if we just simply average across canopy, then      .
#  KOG01 is very height, it does have lots of tall tree without small tree     .
#  layers.                                                                     .
#                                                                              .
#...............................................................................



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                    Parea                                 ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

NNNNN=openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Npercent')%>%
  separate("Leafcode", sep='-',into=c("PLOT","TREE","Branch"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  mutate(sun_shade=str_detect(Branch, ".*SH"))%>%
  filter(!(sun_shade))%>%
  mutate(sun_shade='Sun')

List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Trait_big_table<-NNNNN

#the left join will bring species to the twig density table
Trait_big_table_species<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-Trait_big_table_species%>%
  drop_na(.,any_of('Species'))%>%
  drop_na(.,any_of('Parea'))%>%
  mutate(Outlier_id=outliers::scores(sqrt(Parea),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Parea)))

# Check normal distribution
plot_histogram=qplot(Parea, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Parea_hist.jpg",width = 3, height = 3)
# The data is slightly right-skewed, square root it
Trait_big_table$Parea_t <- sqrt(Trait_big_table$Parea)
plot_histogram=qplot(Parea_t, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="Parea_t_hist.jpg",width = 3, height = 3)

#Variance partitioning
varp <- varpart (Trait_big_table$Parea_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/Parea.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Parea (g m-2)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>% 
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(Parea_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(Parea_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(Parea_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom)

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Parea_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Vessel diameter                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))
Trait_big_table<-drop_na(Trait_big_table,any_of('hydraulically.weighted.diam.micron'))%>%
  mutate(Outlier_id=outliers::scores(sqrt(hydraulically.weighted.diam.micron),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(hydraulically.weighted.diam.micron)))

# Check normal distribution
plot_histogram=qplot(hydraulically.weighted.diam.micron, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="VD_hist.jpg",width = 3, height = 3)
# The data is slightly right-skewed, square root it
plot_histogram=qplot(sqrt(hydraulically.weighted.diam.micron), data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="VD_sqrt_hist.jpg",width = 3, height = 3)
Trait_big_table$VD_t <- sqrt(Trait_big_table$hydraulically.weighted.diam.micron)

#Variance partitioning
varp <- varpart (Trait_big_table$VD_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/VD.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Vessel diameter (micron)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(VD_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(VD_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(VD_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom)

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/VesselDiameter_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Vessel density (vessel per area)                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Trait_big_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Ks')%>%
  separate("Sample/Tree.code", sep='-',into=c("PLOT","TREE", "stem"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))

#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Trait_big_table<-left_join(Trait_big_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Trait_big_table<-drop_na(Trait_big_table,any_of('Species'))
Trait_big_table<-drop_na(Trait_big_table,any_of('vessel.per.area.mm2'))%>%
  mutate(Outlier_id=outliers::scores(sqrt(vessel.per.area.mm2),type="iqr"))%>%
  filter(Outlier_id<1.5)%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(vessel.per.area.mm2)))

# Check normal distribution
plot_histogram=qplot(vessel.per.area.mm2, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="VesselDensity_hist.jpg",width = 3, height = 3)
# The data is right-skewed, square root it
plot_histogram=qplot(sqrt(vessel.per.area.mm2), data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="VesselDensity_sqrt_hist.jpg",width = 3, height = 3)
Trait_big_table$VD_t <- sqrt(Trait_big_table$vessel.per.area.mm2)

#Variance partitioning
varp <- varpart (Trait_big_table$VD_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/VesselDensity.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Vessel density (mm-2)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(VD_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(VD_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(VD_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom)

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/VesselDensity_cwd.xlsx')



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Stem respiration                     ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# this paragraph explain how to change in-situ stem respiration measured with chamber into a trait that indicates maintanace cost in xylem:


#  In Colin Prentice's Ecology Letter 2014 paper, you mentioned Rs as leaf-specific         
#  maintenance respiration cost of the sapwood.                                 
#                                                                               
#  Rs = vH * rs * h * ps / 2                                                    
#                                                                               
#  Cost of water transportation = a = Rs / Evapotranspiration                   
#                                                                               
#  Where, vH (Huber value); rs (sapwood specific respiration rate (s-1) and     
#  ps(mol m-3) sapwood density, h (tree height)                                 
                                                                          
#  I recently plot out stem respiration (measured with EGM-5 Portable CO 2      
#  Gas Analyzer attached to stem surface) per stem surface area in Ghana        
#  (attached in this email), my question is that, what would such a             
#  measurement represent? Rs or rs?                                             


#  Dear Huanyuan,                                                               
#                                                                               
#  These are in situ measurements of stem respiration, expressed per unit of    
#  stem surface area. That is not quite the same as either Rs or rs, but it     
#  is related to rs x ps.                                                      


#...............................................................................
#                                                                              .
#  It is about the conversion from unit stem surface area to unit stem mass    .
#  to unit leaf area, so, let in-situ stem respiration as Rarea, then          .
#                                                                              .
#  Rarea * stem surface area = rs * ps * stem mass                             .
#                                                                              .
#  Rarea * DBH * h * pi = rs * ps * pi * (DBH/2)*2 * h                         .
#                                                                              .
#  So, Rarea = rs * ps * DBH/4                                                 .
#                                                                              .
#  Rs = vH * rs * h * ps / 2 = vH * 2 * Rarea * h / DBH                        .
#                                                                              .
#  (DBH - diameter of stem at breast height)                                   .
#                                                                              .
#  Now Rs and Rarea have the same unit (umol m-2 s-1). I have all value on     .
#  hand so I have Rs calculated, and it looks like this (a nice and            .
#  significant increasing trend from wet to dry!):                             .
#                                                                              .
#...............................................................................



Trait_big_table <- read.csv("stem_respiration_Ghana_per_tree_20220101_HYZZ.csv")%>%
  separate("plot_code", sep='-',into=c("LETTERS","NUMBER"),remove = F)%>%
  mutate(PLOT=paste0(LETTERS,NUMBER))%>%
  drop_na(any_of('flux_umolm2sec_no_outlier'))%>%
  mutate(SR_final=ifelse(flux_umolm2sec_no_outlier>=0, flux_umolm2sec_no_outlier, 0))%>%
  mutate(Outlier_id_=outliers::scores(sqrt(SR_final),type="iqr"))%>%
  filter(Outlier_id_<1.5)%>%
  mutate(PLOT=str_replace(PLOT,"KOG02","KOG01"))%>%
  mutate(PLOT=str_replace(PLOT,"KOG04","KOG03"))%>%
  mutate(PLOT=str_replace(PLOT,"KOG05","KOG04"))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))



# Count number of measurements
nom <- Trait_big_table%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(SR_final)))

# Check normal distribution
plot_histogram=qplot(SR_final, data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="StemRespiration_hist.jpg",width = 3, height = 3)
# The data is right-skewed, square root it
plot_histogram=qplot(sqrt(SR_final), data=Trait_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="StemRespiration_sqrt_hist.jpg",width = 3, height = 3)
Trait_big_table$SR_t <- sqrt(Trait_big_table$SR_final)

#Variance partitioning
varp <- varpart (Trait_big_table$SR_t, ~ PLOT, ~Species, data = Trait_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/StemRespiration.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))+
  title(main = "Rstem_stem (umol m-2 s-1)")
dev.off()

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))


# lm for significant test
mod <- lm(SR_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(SR_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(SR_t,weightsss,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)%>%
  left_join(nom)

openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/StemRespiration_cwd.xlsx')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                        Specific stem respiration                         ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in datasheets
ANK01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Ankasa_Ank1')[,c(3,8,13)]%>%
  drop_na()%>%
  mutate(PLOT='ANK01')
ANK03<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Ankasa_Ank3')[,c(3,8,13)]%>%
  drop_na()%>%
  mutate(PLOT='ANK03')
BOB01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Bobiri_SNR')[,c(4,9,14)]%>%
  drop_na()%>%
  mutate(PLOT='BOB01')%>%
  rename(Tag.No=TreeTag)
BOB02<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Bobiri_Y53')[,c(5,10,15)]%>%
  drop_na()%>%
  mutate(PLOT='BOB02')%>%
  rename(Tag.No=TreeTag)
KOG01<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_F1')[,c(2,8,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG01')%>%
  rename(Tag.No=Tag)
KOG03<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_FST')[,c(2,8,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG03')%>%
  rename(Tag.No=New.Tag)
KOG04<-openxlsx::read.xlsx('Census_data/Ghana_Census2020_April2021_JA_tree_height.xlsx',sheet='Kogyae_Sav1')[,c(2,8,13)]%>%
  drop_na()%>%
  mutate(PLOT='KOG04')%>%
  rename(Tag.No=Tag)

Tree_height_table <- rbind(ANK01, ANK03, BOB01, BOB02, KOG01, KOG03, KOG04)
#the left join will bring species to the table
Census_big_table2<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  #filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  rename(PLOT=plot_code)
Trait_big_table<-left_join(Tree_height_table,Census_big_table2[c('Tag.No','Species','PLOT')],by=c("PLOT","Tag.No"))%>%
  drop_na(any_of('Species'))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))

Trait_big_table2_tree_height<-Trait_big_table%>%
 mutate(height_dbh=Htotal_m/DBH_cm*100)%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(PLOT,'_T',Tag_ID))



Stem_resp <- read.csv("stem_respiration_Ghana_per_tree_20220101_HYZZ.csv")%>%
  separate("plot_code", sep='-',into=c("LETTERS","NUMBER"),remove = F)%>%
  mutate(PLOT=paste0(LETTERS,NUMBER))%>%
  drop_na(any_of('flux_umolm2sec_no_outlier'))%>%
  mutate(SR_final=ifelse(flux_umolm2sec_no_outlier>=0, flux_umolm2sec_no_outlier, 0))%>%
  mutate(PLOT=str_replace(PLOT,"KOG02","KOG01"))%>%
  mutate(PLOT=str_replace(PLOT,"KOG04","KOG03"))%>%
  mutate(PLOT=str_replace(PLOT,"KOG05","KOG04"))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))%>%
                     filter(grepl('^[0-9]{1,3}$', tree_tag))%>%
                   mutate(Tag_ID=sprintf("%03d", as.numeric(tree_tag)))%>%
                   mutate(PLOT_TREE=paste0(PLOT,'_T',Tag_ID))


List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')

Huber_table<-openxlsx::read.xlsx('Overall_trait_raw_data.xlsx',sheet='Huber_new')%>%
  separate("Code", sep='-',into=c("PLOT","TREE","STEM"),remove = F)%>%
  mutate(PLOT_TREE=paste0(PLOT,'_',TREE))%>%
  rename(Huber=!!as.name('SapwoodArea/LeafArea_mm2/cm2'))%>%
  mutate(Huber=Huber/100) # change into m2/m2 when /100






#the left join will bring species to the twig density table
Census_big_table<-read.csv("Census_data/Census_table_all_plots_Species_list.csv")%>%
  filter(grepl('^[0-9]{1,3}$', Tag.No))%>%
  mutate(Tag_ID=sprintf("%03d", as.numeric(Tag.No)))%>%
  mutate(PLOT_TREE=paste0(plot_code,'_T',Tag_ID))
Huber_table<-left_join(Huber_table,Census_big_table[c('Species','PLOT_TREE','WD')],by="PLOT_TREE")
Huber_table<-drop_na(Huber_table,any_of('Species'))%>%

  mutate(Outlier_id=scores_na((Huber),type="iqr"))%>%
  mutate(Huber=ifelse(abs(Outlier_id) < 1.5, Huber, NA))%>%
  mutate(Species_PLOT = paste0(Species,PLOT))
Huber_table<-drop_na(Huber_table,any_of('Huber'))%>%
  group_by(Species_PLOT)%>%
  summarise(Huber=mean(Huber,na.rm=T))

 Trait_big_table<-Stem_resp%>%
   left_join(Trait_big_table2_tree_height[,c('PLOT_TREE','height_dbh')],by='PLOT_TREE')%>%

   left_join(Huber_table[,c('Huber','Species_PLOT')],by='Species_PLOT')%>%
      mutate(RS=SR_final*height_dbh*Huber * 2)%>%
   group_by(Species_PLOT)%>%
   summarise(RS=mean(RS,na.rm=T),
             Species=first(Species),
             PLOT=first(PLOT))%>%
   filter(RS>0)
#   left_join(Huber_table[,c('Huber','PLOT_TREE')],by='PLOT_TREE')%>%
#   mutate(RS = SR_final*Huber*Htotal_m)%>%
#   filter(RS>0)


#...............................................................................
#                                                                              .
#  I tried to link three datasets on individual level but this might miss      .
#  some key species (no common tree), so I have to link them on species level  .
#                                                                              .
#...............................................................................


plot_histogram=qplot(log10(RS), data=Trait_big_table, geom="histogram") 
plot_histogram
Trait_big_table$RS_t <- log10(Trait_big_table$RS)

Trait_big_table2<-Trait_big_table%>%
  left_join(species_basal[,c('Species','total_basal','PLOT')],by=c('Species','PLOT'))%>%
  rowwise()%>%
  mutate(number_record= (sum(as.vector(Trait_big_table$Species_PLOT)==Species_PLOT)))%>%
  # This is to count how many individuals are measured under one species, the basal area is the total basal area for a species in one plot, but in Trait_big_table, only some individuals are presented
  mutate(weightsss = total_basal/ number_record)%>%
  mutate(PLOT=as.factor(PLOT))

# lm for significant test
mod <- lm(RS_t ~ PLOT, data = Trait_big_table2, weights = weightsss)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod ,weights = weightsss)
plot(mod.cld)


Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

# Count number of measurements
nom <- Trait_big_table2%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(RS_t)))

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = wtd.mean(RS_t,weightsss,na.rm=T),
            standard_d_trandformed=weighted.var.se(RS_t,weightsss,na.rm=T))%>%
  mutate(mean_original = 10^mean_value_trandformed,
         upper = 10^(mean_value_trandformed+ standard_d_trandformed),
         lower = 10^(mean_value_trandformed- standard_d_trandformed))%>%
  left_join(Letter_significance)%>%
  left_join(nom, by="PLOT")
openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Specific_StemRespiration_cwd.xlsx')

# Try to just join on plot level:

Stem_respiration<-openxlsx::read.xlsx('../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/StemRespiration_cwd.xlsx')%>%
  dplyr::select(PLOT,mean_original)%>%
  rename(Stem_resp=mean_original)

Huber_table<-openxlsx::read.xlsx('../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/Huber_cwd.xlsx')%>%
  dplyr::select(PLOT,mean_original)%>%
  rename(Huber_value=mean_original)

overall_table<-openxlsx::read.xlsx('../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/TreeHeight_cwd.xlsx')%>%
  dplyr::select(PLOT,mean_original)%>%
  rename(Max_Height=mean_original)%>%
  left_join(Huber_table)%>%
  left_join(Stem_respiration)%>%
  mutate(Rs = Huber_value*Stem_resp)

# OK it is shit idea... worse results.

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                                                            --
##................................. DRAWING.................................----
##                                                                            --
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#setwd('F:/Oxford/Chapter_one/Ghana_rainfall_trait_variation_optimality_github/Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm')
setwd('../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm')
#List_of_f<-list.files(pattern = '*_cwd.xlsx$')
library(tidyverse)
library(ggplot2)
library(ggtext)
library(ggsci)
library(ggpubr)

#List_of_f
list1 <- c("Amax_cwd.xlsx","Asat_cwd.xlsx","chi_cwd.xlsx","DResp_cwd.xlsx","Jmax_ACi_25C_cwd.xlsx","Vcmax_ACi_25C_cwd.xlsx",'Specific_StemRespiration_cwd.xlsx','StemRespiration_cwd.xlsx')
list2 <- c("Huber_cwd.xlsx","Ks_cwd.xlsx","VesselDensity_cwd.xlsx","VesselDiameter_cwd.xlsx","TwigDensity.g.cm3_cwd.xlsx","WD_cwd.xlsx","Turgor_cwd.xlsx","TreeHeight_cwd.xlsx")
list3 <- c("LMA_cwd.xlsx","Narea_cwd.xlsx","Ngkg_cwd.xlsx","Pgkg_cwd.xlsx","Parea_cwd.xlsx")

# Read in standardized names
SN <- read.csv("StandadizedName.csv")

# Give colours to each plot
cls <- data.frame(plot_code=List_of_site, colour=c("#deebf7","#deebf7","#deebf7","#9ecae1","#9ecae1","#3182bd","#3182bd"))

plot_function <- function(List_of_f) {
  # Initialize
  points_plot <- list()
  
for (file_ind in 1:length(List_of_f)) {
cwm_table<-openxlsx::read.xlsx(List_of_f[file_ind])

ylabel_name<-filter(SN, FileName==List_of_f[file_ind])$SN
colnames(cwm_table)[c(1,4)]<-c("plot_code","value")
#cwm_table$standard_dev<-as.numeric(cwm_table$standard_dev)  #change variance into standard deviation
#cwm_table[,2]<-as.numeric(cwm_table[,2])
cwm_table <- left_join(cwm_table, cls, by="plot_code")
cwm_table$plot_code <- factor(cwm_table$plot_code, levels = c("ANK03", "ANK01", "BOB02", "BOB01", "KOG01", "KOG03", "KOG04"))

if (min(cwm_table$value,na.rm=T)>0) {
  points_plot[[file_ind]]<-ggplot(cwm_table) +
  geom_bar( aes(x=plot_code, y=value), stat="identity", fill=cwm_table$colour, alpha=0.7,width = 0.8) +
  geom_errorbar(aes(x=plot_code, y=value, ymin=lower, ymax=upper),alpha=0.9, size=1.3,width=0.4)+ theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian( ylim = c(min(cwm_table$value,na.rm=T)*0.8,max(cwm_table$upper,na.rm=T)*1.2)) +
  geom_text(y = max(cwm_table$upper)*1.18,aes(x=plot_code, y=value,label=Letters)) +
  labs(title = ylabel_name) +
  theme(title = element_text(size=9),
        plot.title = element_markdown()) +
  geom_text(y = max(cwm_table$upper)*1.07,aes(x=plot_code, y=value,label=NumberOfMeasurements), colour="#636363", size=3)+
  scale_x_discrete(labels=c('ANK-03','ANK-01','BOB-02','BOB-01','KOG-02','KOG-04','KOG-05'))
}else{
  points_plot[[file_ind]]<-ggplot(cwm_table) +
  geom_bar( aes(x=plot_code, y=value), stat="identity", fill=cwm_table$colour, alpha=0.7,width = 0.8) +
  geom_errorbar(aes(x=plot_code, y=value, ymin=lower, ymax=upper),alpha=0.9, size=1.3,width=0.4)+ theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian( ylim = c(min(cwm_table$lower,na.rm=T)*1.2,max(cwm_table$value,na.rm=T)*0.8)) +
  geom_text(y = min(cwm_table$lower)*1.18,aes(x=plot_code, y=value,label=Letters)) +
  labs(title = ylabel_name) +
  theme(title = element_text(size=9),
        plot.title = element_markdown()) +
  geom_text(y = min(cwm_table$lower)*1.09,aes(x=plot_code, y=value,label=NumberOfMeasurements), colour="#636363", size=3)+
  scale_x_discrete(labels=c('ANK-03','ANK-01','BOB-02','BOB-01','KOG-02','KOG-04','KOG-05'))
}


#ggsave(filename = paste0(List_of_f[file_ind],".jpg"),bar_plot_errorbar, width = 4, height = 4)

# points_plot<-ggplot(cwm_table, aes_string(x="PLOT", y=colnames(cwm_table)[2])) + 
#   geom_point(
#     color="orange",
#     fill="#69b3a2",
#     shape=21,
#     alpha=0.5,
#     size=6,
#     stroke = 2
#   ) +geom_text(aes(label=plot_code))
#ggsave(filename = paste0(List_of_f[file_ind],".jpg"),points_plot, width = 3, height = 3)
}
  return(points_plot)
}


points_plot <- plot_function(list1)
figure <- ggarrange(plotlist=points_plot,
                    ncol = 4, nrow = 2)
figure
ggsave(filename = "PLOT_trait_variation_overall_new1_new.pdf",figure,height=5, width = 10)

points_plot <- plot_function(list2)
figure <- ggarrange(plotlist=points_plot,
                    ncol = 4, nrow = 2)
figure
ggsave(filename = "PLOT_trait_variation_overall_new2_new.pdf",figure,height=5, width = 10)

points_plot <- plot_function(list3)
figure <- ggarrange(plotlist=points_plot,
                    ncol = 4, nrow = 2)
figure
ggsave(filename = "PLOT_trait_variation_overall_new3_new.pdf",figure,height=5, width = 10)



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                temperature and moisture                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

List_of_site<-c('KOG01','KOG03','KOG04','BOB01','BOB02' ,'ANK01','ANK03')
List_of_f<-c("mean_annual_temp.xlsx", "moisture.xlsx")

# Read in standardized names
SN <- read.csv("StandadizedName.csv")

### The temperature plot

# Give colours to each plot
cls <- data.frame(plot_code=c("KOG","BOB","ANK"), colour=c("#deebf7","#9ecae1","#3182bd"))

cwm_table<-openxlsx::read.xlsx("mean_annual_temp.xlsx")
ylabel_name<-filter(SN, FileName=="mean_annual_temp.xlsx")$SN
colnames(cwm_table)[c(1,4)]<-c("plot_code","value")
cwm_table <- left_join(cwm_table, cls, by="plot_code")
cwm_table$plot_code <- factor(cwm_table$plot_code, levels = c("ANK", "BOB", "KOG"))

temp<-ggplot(cwm_table) +
    geom_bar( aes(x=plot_code, y=value), stat="identity", fill=cwm_table$colour, alpha=0.7,width = 0.8) +
    theme_bw()+
    ylab("")+
    xlab("")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_cartesian( ylim = c(min(cwm_table$value,na.rm=T)*0.8,max(cwm_table$upper,na.rm=T)*1.2)) +
    geom_text(y = max(cwm_table$upper)*1.18,aes(x=plot_code, y=value,label=Letters)) +
    labs(title = ylabel_name) +
    theme(title = element_text(size=9))
temp  

### The moisture plot

# Give colours to each plot
cls <- data.frame(plot_code=List_of_site, colour=c("#deebf7","#deebf7","#deebf7","#9ecae1","#9ecae1","#3182bd","#3182bd"))

cwm_table<-openxlsx::read.xlsx("moisture.xlsx")
ylabel_name<-filter(SN, FileName=="moisture.xlsx")$SN
colnames(cwm_table)[c(1,4)]<-c("plot_code","value")
cwm_table <- left_join(cwm_table, cls, by="plot_code")
cwm_table$plot_code <- factor(cwm_table$plot_code, levels = c("ANK03", "ANK01", "BOB02", "BOB01", "KOG01", "KOG03", "KOG04"))

moi<-ggplot(cwm_table) +
  geom_bar( aes(x=plot_code, y=value), stat="identity", fill=cwm_table$colour, alpha=0.7,width = 0.8) +
  geom_errorbar(aes(x=plot_code, y=value, ymin=lower, ymax=upper),alpha=0.9, size=1.3,width=0.4)+ theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian( ylim = c(min(cwm_table$value,na.rm=T)*0.8,max(cwm_table$upper,na.rm=T)*1.2)) +
  geom_text(y = max(cwm_table$upper)*1.18,aes(x=plot_code, y=value,label=Letters)) +
  labs(title = ylabel_name) +
  theme(title = element_text(size=9))+
  scale_x_discrete(labels=c('ANK-03','ANK-01','BOB-02','BOB-01','KOG-02','KOG-04','KOG-05'))
moi

### Precipitation
cls <- data.frame(plot_code=c("KOG","BOB","ANK"), colour=c("#deebf7","#9ecae1","#3182bd"))

cwm_table<-openxlsx::read.xlsx("Mean_annul_precipitation.xlsx")
ylabel_name<-filter(SN, FileName=="Mean_annul_precipitation.xlsx")$SN
colnames(cwm_table)[c(1,4)]<-c("plot_code","value")
cwm_table <- left_join(cwm_table, cls, by="plot_code")
cwm_table$plot_code <- factor(cwm_table$plot_code, levels = c("ANK", "BOB", "KOG"))

precipitation_plot<-ggplot(cwm_table) +
  geom_bar( aes(x=plot_code, y=value), stat="identity", fill=cwm_table$colour, alpha=0.7,width = 0.8) +
  theme_bw()+
  ylab("")+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  coord_cartesian( ylim = c(min(cwm_table$value,na.rm=T)*0.8,max(cwm_table$upper,na.rm=T)*1.2)) +
  geom_text(y = max(cwm_table$upper)*1.18,aes(x=plot_code, y=value,label=Letters)) +
  labs(title = ylabel_name) +
  theme(title = element_text(size=9))
precipitation_plot  


points_plot <- list(temp, moi,precipitation_plot)
figure <- ggarrange(plotlist=points_plot,
                    ncol = 3, nrow = 1)
figure
ggsave(filename = "TempAndMoi.jpg",figure,height=2.3, width = 2.3*4)


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                census data                               ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#...............................................................................
#                                                                              .
#  this is to generate species basal area files, run only once                 .
#                                                                              .
#...............................................................................




# We first deploy F:/Side_project/GEM/Organize_census_data_for_Marc.R to tidy up data downloaded from forest plot
setwd("F:/Oxford/Chapter_one/From_imma/Traits/plotV3_log/")

Files_to_read_census<-list.files(path = "../Census_data", pattern = '*cleaned*', full.names = T) #Input folder

Census_big_table <- lapply(Files_to_read_census, read.csv) %>%
  data.table::rbindlist()%>%
  mutate(census_year_plot=paste0(str_sub(census_date,1,4),plot_code))


#Use this to see all the plot_year combination: unique(Census_big_table$census_year_plot)
#we only want year close to 2016 because the trait data were collected in 2016
Census_big_table<-Census_big_table%>%
  filter(
    (census_year_plot=="2016ANK_01")|
      (census_year_plot=="2016ANK_03")|
      (census_year_plot=="2016BOB_01")|
      (census_year_plot=="2016BOB_02")|
      (census_year_plot=="2015KOG_02")|
      (census_year_plot=="2015KOG_04")|
      (census_year_plot=="2015KOG_05")
  )%>%
  select(Tag.No,Family,Species,DBH1,Census.Notes,census_date,plot_code,DBH0,WD)%>%
  mutate(basal= pi*(DBH1/2)^2)%>%
mutate(plot_code=str_replace(plot_code,"KOG_02","KOG01"))%>%
  mutate(plot_code=str_replace(plot_code,"KOG_04","KOG03"))%>%
  mutate(plot_code=str_replace(plot_code,"BOB_02","BOB02"))%>%
  mutate(plot_code=str_replace(plot_code,"KOG_05","KOG04"))%>%
  mutate(plot_code=str_replace(plot_code,"BOB_01","BOB01"))%>%
  mutate(plot_code=str_replace(plot_code,"ANK_01","ANK01"))%>%
  mutate(plot_code=str_replace(plot_code,"ANK_03","ANK03"))%>%
  mutate(SITE=substring(plot_code,1,3))
  
List_of_site<-c('KOG','BOB','ANK')
#we change plot name becausse Census table does not match with Trait table
for (num_of_site in 1:length(List_of_site)) {
  Census_summary_table<-Census_big_table%>%
    filter(SITE==List_of_site[num_of_site])%>%
    group_by(Species) %>%
    dplyr::summarise(total_basal = sum(basal,na.rm = TRUE), count=n(), wood_density=mean(WD,na.rm = TRUE))%>%
    arrange(desc(total_basal))
   write.csv(Census_summary_table,file=paste0('../Census_data/',List_of_site[num_of_site],'Species_list.csv'))
}
write.csv(Census_big_table,file='../Census_data/Census_table_all_plots_Species_list.csv')

Census_big_table_by_species_basal<-Census_big_table%>%
    group_by(Species,plot_code) %>%
    dplyr::summarise(total_basal = sum(basal,na.rm = TRUE), count=n(), wood_density=mean(WD,na.rm = TRUE))%>%
    arrange(desc(total_basal))%>%
   mutate(PLOT=plot_code)
write.csv(Census_big_table_by_species_basal,file='../Census_data/Census_table_all_plots_Species_list_basal_area.csv')


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##                                old version of tree_height (archived)  ----
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Files_to_read_census<-list.files(path = "Census_data", pattern = '*cleaned*', full.names = T) #Input folder

kog_04_tree_height<-openxlsx::read.xlsx('Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='Savanna 1')[,c(22,9)]%>%
  mutate(plot_code='KOG_04')
kog_01_tree_height<-openxlsx::read.xlsx('Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='remnant forest 1 ')[,c(20,9)]%>%
  mutate(plot_code='KOG_01')
kog_03_tree_height<-openxlsx::read.xlsx('Census_data/Kogyae census main plots(WORK ON THIS).xlsx',sheet='transition')[,c(13,4)]%>%
  mutate(plot_code='KOG_03')

colnames(kog_01_tree_height)<-c('Height','Species','plot_code')
 colnames(kog_03_tree_height)<-c('Height','Species','plot_code')
colnames(kog_04_tree_height)<-c('Height','Species','plot_code')


kog_tree_height<-rbind(kog_04_tree_height,kog_03_tree_height,kog_01_tree_height)%>%
  mutate(Height=as.numeric(Height))%>%
    drop_na()

Census_big_table <- lapply(Files_to_read_census, read.csv) %>%
  data.table::rbindlist()%>%
  filter(plot_code=='ANK_01' |plot_code=='ANK_03' |plot_code=='BOB_01' |plot_code=='BOB_02'|plot_code=='KOG_02')%>%
  mutate(year=substr(census_date,1,4))%>%
  filter(year==2012)%>%
  mutate(census_year_plot=paste0(str_sub(census_date,1,4),plot_code))%>%
  mutate(Height=ifelse(plot_code=='BOB_01',exp( 0.6152+ 0.6362*log(DBH4/10)+0.0109*25813371/1000000),Height),
         Height=ifelse(plot_code=='BOB_02',exp( 0.6152+ 0.6362*log(DBH4/10)+0.0109*32814313/1000000),Height))%>%
 # mutate(Height=ifelse(plot_code=='KOG_02',exp( 0.4619+ 0.6362*log(DBH4/10)+0.0109*21446610/1000000),Height))%>%
  filter(Height>0)%>%
 #   mutate(Height=ifelse(plot_code=='BOB_01',exp( 0.8946+ 0.6365*log(DBH4/10)),Height),
 #        Height=ifelse(plot_code=='BOB_02',exp( 0.8946+ 0.6365*log(DBH4/10)),Height))%>%filter(Height>0)%>%
  drop_na(Height)%>%
  mutate(plot_code = str_replace_all(plot_code,'KOG_02','KOG_01'))
  
Census_big_table<-rbind(Census_big_table,kog_tree_height,fill=TRUE)


# Check normal distribution
plot_histogram=qplot(Height, data=Census_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="TreeHeight_hist.jpg",width = 3, height = 3)
# The data is right-skewed, sqrt it
plot_histogram=qplot(sqrt(Height), data=Census_big_table, geom="histogram") 
plot_histogram
ggsave(plot_histogram,file="TreeHeight_sqrt_hist.jpg",width = 3, height = 3)
Census_big_table$Height_t <- sqrt(Census_big_table$Height)

#Variance partitioning
varp <- varpart (Census_big_table$Height_t, ~ plot_code, ~Species, data = Census_big_table)
png("../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/varp_figures/TreeHeight.png", width = 400, height = 300)
plot (varp, digits = 2, Xnames = c('Plot', 'Species'), bg = c('navy', 'green'))
dev.off()

Trait_big_table2 <- Census_big_table%>%separate("plot_code", sep='_',into=c("LETTER","NUMBER"),remove = F)%>%
  mutate(PLOT=paste0(LETTER,NUMBER))%>%
  mutate(PLOT=as.factor(PLOT))%>%
  group_by(PLOT)%>%
  mutate(rank=ecdf(Height)(Height))%>%
  filter(rank<0.96 & rank>0.9)
 
# Count number of measurements
nom <- Trait_big_table2%>%group_by(PLOT)%>%
  summarise(NumberOfMeasurements=sum(!is.na(Height)))

# for tree height, we don't need community weighted mean because there might be lots of small tree which will not represent canopy height, we choose 85% to 90% of tree height


# lm for significant test
mod <- lm(Height_t ~ PLOT, data = Trait_big_table2)
summary(glht.mod <- glht(mod, mcp(PLOT = "Tukey")))
mod.cld <- cld(glht.mod)
plot(mod.cld)

Letter_significance <- as.data.frame(mod.cld$mcletters$Letters)%>%
  rename(Letters=!!as.name('mod.cld$mcletters$Letters'))%>%
  rownames_to_column(var='PLOT')

Trait_big_table_cwm<-Trait_big_table2%>%
  group_by(PLOT)%>%
  summarise(mean_value_trandformed = mean(Height_t,na.rm=T),
            standard_d_trandformed=standard_error_calc(Height_t,na.rm=T))%>%
  mutate(mean_original = mean_value_trandformed^2,
         upper = (mean_value_trandformed+ standard_d_trandformed)^2,
         lower = (mean_value_trandformed- standard_d_trandformed)^2)%>%
  left_join(Letter_significance)


# Re-readin the data to merge with NOM
Trait_big_table_cwm<-Trait_big_table_cwm%>%
  left_join(nom, by="PLOT")
openxlsx::write.xlsx(Trait_big_table_cwm,file='../Ghana_rainfall_trait_variation_optimality/Site_to_site_compare_and_variance_partitioning/processed_trait_cwm/TreeHeight_cwd.xlsx',overwrite = T)
