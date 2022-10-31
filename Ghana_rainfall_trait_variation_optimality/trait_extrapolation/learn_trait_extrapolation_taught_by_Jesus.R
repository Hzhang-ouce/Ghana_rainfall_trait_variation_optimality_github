####################################

library(BHPMF)

    data(trait.info)  # Read the matrix X
    data(hierarchy.info) # Read the hierarchy information
setwd('F:/Oxford/Chapter_one/From_imma/Traits/trait_extrapolation')

hierarchy.info<-hierarchy.info[,1:2]
trait.info<-trait.info[,1]
alista<-cbind(trait.info,hierarchy.info)%>%
  group_by(species)%>%
  summarise(trait1=mean(trait.info,na.rm=T))%>%
  filter(is.nan(trait1))
alista<-alista$species
full_table<-cbind(trait.info,hierarchy.info)%>%
  filter(!species %in% alista)

 trait.info <-as.matrix(full_table$trait.info)
 hierarchy.info <-as.matrix(full_table[,2:3])

    ###################
    # Usage 1: Figure 1
    ###################
    GapFilling(as.matrix(trait.info), hierarchy.info,
  	          mean.gap.filled.output.path = "mean_gap_filled.txt",
    	       	std.gap.filled.output.path="std_gap_filled.txt")
    
    ###################
    # Usage 2: Figure 2
    ###################
    GapFilling(trait.info, hierarchy.info,
              prediction.level = 2,
              used.num.hierarchy.levels = 1,
  	          mean.gap.filled.output.path = "mean_gap_filled.txt",
    	       	std.gap.filled.output.path="std_gap_filled.txt")
    
    ###################
    # Usage 3: Figure 3
    ###################
    GapFilling(trait.info, hierarchy.info,
              prediction.level = 4,
              used.num.hierarchy.levels = 3,
  	          mean.gap.filled.output.path = "mean_gap_filled.txt",
    	       	std.gap.filled.output.path="std_gap_filled.txt")

m<-read.csv("mean_gap_filled.txt", sep="")
m_filled<-cbind(hierarchy.info,m)   
m_ori<-cbind(hierarchy.info,trait.info) 

m_filled<-m_filled%>%
  group_by(species)%>%
  summarise(Trait7=mean(Trait7,na.rm=T))

m_ori<-m_ori%>%
  group_by(species)%>%
  summarise(Trait7=mean(Trait7,na.rm=T))

Df5_2<-Df5[-which(rowSums(is.na(Df5[,4:23]))==ncol(Df5[,4:23])),]#exclude  records with only NA's 

hierarchy<-select(Df5_2,c(1:3));hierarchy<-mutate(hierarchy,id=row.names(hierarchy))%>%select(4,1,2,3)

trait                   <-select(Df5_2,c(4:9,12,13,15:23))

trait_log10             <-log10(trait[,c(2:9,11:17)])%>%mutate(Psi_TLP._Mpa4.073=log10(trait$Psi_TLP._Mpa+4.073))#Only log10

trait_log10_psi_diam_Out<-select(trait, c(2:9,11:17))%>%log10                                                   #Only log10 

trait_log10_psi_log10_scale<-scale(trait_log10_psi_diam_Out)%>%as.data.table                     #Also Scaled


### when I want to convert it back use for the log10 transformed traits .....  x.orig = 10^(t(apply(trait1, 1, 
#function(r)r*attr(trait1,'scaled:scale') + attr(trait1, 'scaled:center'))))

GapFilling(as.matrix(trait_log10), hierarchy,
           prediction.level = 4,
           used.num.hierarchy.levels = ,
           mean.gap.filled.output.path = "D:/PostDoc/Oxford/DATA/Ghana_SophieFauset/2_Spp_Traits_Plot_Tables/GF_mean_ind_log10.csv",
           std.gap.filled.output.path="D:/PostDoc/Oxford/DATA/Ghana_SophieFauset/2_Spp_Traits_Plot_Tables/GF_std_ind_log10.csv")
### the resulting RMSE for the test data: 0.174


###############################################################################################
### Once the table to use is selected then put back the hierarchy infomration on the results
##########################
m<-read.csv("mean_gap_filled.txt", sep="")
m1<-cbind(hierarchy,m)            #Put the species, genus and family details back in the table
m2<-as.data.table(m1)
tomean <- c(names(m2[,c(5:20)]))  # Average to have 1 record per spp and be able to merge them to Df # or 5:20 or 5:19 for other files
m2[, (tomean) := lapply(tomean, function(x) {
  x <- get(x)
  x <- mean(x, na.rm = TRUE)
  x
}), by = Species] 
m2 <-unique(m2[,c(2:20)])         #Extract unique species

###############################################################################
########### Eliminate gap filled records with Coefficient of Variation above 1 
m<-read.csv("D:/PostDoc/Oxford/DATA/Ghana_SophieFauset/2_Spp_Traits_Plot_Tables/GF_mean_ind_log10.csv", sep="")
s<-read.csv("std_gap_filled.txt", sep="")

cv<-s/m; head(cv)

list_t<-list()
for (i in 1:16){      # 1:15 for GF_mean_ind_log10_PsiDiameterOut.csv
  c<-cv[i]
  me<-m[i]
  cme<-cbind(c,me)
  cme$final<-ifelse(cme[,1]>1,NA,cme[,2])
  cme<-cme[3]
  names(cme)<-paste(names(me))
  list_t[[i]]<-cme
} 
m_CV_above1<-do.call(cbind,list_t)
head(m_CV_above1)
nrow(m_CV_above1)

###Because I have values that do not make sense lets replace these outliers with NA's

outlierKD <- function(dt, var) {
  var_name <- eval(substitute(var),eval(dt))
  tot <- sum(!is.na(var_name))
  na1 <- sum(is.na(var_name))
  m1 <- mean(var_name, na.rm = T)
  par(mfrow=c(2, 2), oma=c(0,0,3,0))
  boxplot(var_name, main="With outliers")
  hist(var_name, main="With outliers", xlab=NA, ylab=NA)
  outlier <- boxplot.stats(var_name)$out
  mo <- mean(outlier)
  var_name <- ifelse(var_name %in% outlier, NA, var_name)
  boxplot(var_name, main="Without outliers")
  hist(var_name, main="Without outliers", xlab=NA, ylab=NA)
  title("Outlier Check", outer=TRUE)
  na2 <- sum(is.na(var_name))
  message("Outliers identified: ", na2 - na1, " from ", tot, " observations")
  message("Proportion (%) of outliers: ", (na2 - na1) / tot*100)
  message("Mean of the outliers: ", mo)
  m2 <- mean(var_name, na.rm = T)
  message("Mean without removing outliers: ", m1)
  message("Mean if we remove outliers: ", m2)
  
  dt[as.character(substitute(var))] <- invisible(var_name)
  assign(as.character(as.list(match.call())$dt), dt, envir = .GlobalEnv)
  message("Outliers successfully removed", "\n")
  return(invisible(dt))
} 

outlierKD(m_CV_above1,LMA)
outlierKD(m_CV_above1,SLA.g.m2)
outlierKD(m_CV_above1,Area.Cm2)
outlierKD(m_CV_above1,leafthickness.mm)
outlierKD(m_CV_above1,Ngkg)
outlierKD(m_CV_above1,Pgkg)
outlierKD(m_CV_above1,Asat)
outlierKD(m_CV_above1,Amax)
outlierKD(m_CV_above1,HuberValue_m2.mm2)
outlierKD(m_CV_above1,Vessel.lumen.fraction)
outlierKD(m_CV_above1,Diamenter.mean.micron)
outlierKD(m_CV_above1,Kh.mean)
outlierKD(m_CV_above1,Nr.vessels)
outlierKD(m_CV_above1,Vessel.per.area.mm2)
outlierKD(m_CV_above1,WD)
outlierKD(m_CV_above1,max.H)
outlierKD(m_CV_above1,Psi_TLP._Mpa4.073)

summary(m_CV_above1)

m_CV_above1<-cbind(hierarchy,m_CV_above1)           # Put the species, genus and family details back in the table

m_CV_above1<-as.data.table(m_CV_above1)             
tomean <- c(names(m_CV_above1[,c(5:20)]))           # Average to have 1 record per spp and be able to merge them to Df; 5:19 for other files
m_CV_above1[, (tomean) := lapply(tomean, function(x) {
  x <- get(x)
  x <- mean(x, na.rm = TRUE)
  x
}), by = Species] 

m_CV_above1 <-unique(m_CV_above1[,c(2:20)])      #Extract unique species
m_CV_above2<-merge(m_CV_above1,hpg[,c(1,6,8)],all.x=T,all.y=F) 
names(m_CV_above2)[19]<-"Psi_TLP._Mpa"           # Only for  GF_mean_ind_log10.csv
m_CV_above2<-merge(m_CV_above2,nitrogen,all.x=T,all.y=F)

write.csv(m_CV_above2,"D:/PostDoc/Oxford/DATA/Ghana_SophieFauset/2_Spp_Traits_Plot_Tables/Unique_Spp_Traits_ALLplots_GF_log10_NoOutliers_100117.csv", row.names=F)


####################################

