##get forcing and run splash Ghana plots GEM
library(climateR)
library(rgdal)
library(rsplash)
library(splashTools)
library(xts)
#################################################################################################
#01. load the data
#################################################################################################
#load(file="C:/dsval/Huanyuan_soil/sites_Ghana_forcing.RData")
#load(file="C:/dsval/Huanyuan_soil/splash_GEM_Ghana.RData")
Meta_data_plots_GEM_Ghana <- read.table(file = "C:/dsval/Huanyuan_soil/Meta_data_plots_GEM_Ghana.csv", header = TRUE, sep = ",", fileEncoding = "WINDOWS-1252", quote = "\"", stringsAsFactors = FALSE, comment.char = "", na.strings = c("",'NA'))

sites_centroids<-aggregate(Meta_data_plots_GEM_Ghana, list(Meta_data_plots_GEM_Ghana$Site.Name), mean,na.action = na.omit)

#load soil data
soil_property_plots_GEM_Ghana<- read.table(file = "C:/dsval/Huanyuan_soil/soil_property_from_Malhi_20220227_simplified.csv", header = TRUE, sep = ",", fileEncoding = "windows-1252", quote = "\"", stringsAsFactors = FALSE, comment.char = "", na.strings = "NA")
soil_property_plots_GEM_Ghana<-soil_property_plots_GEM_Ghana[1:11,]
soil_property_plots_GEM_Ghana<-soil_property_plots_GEM_Ghana[order(soil_property_plots_GEM_Ghana$GEM.Plot.Code),]

plots_GEM_Ghana<-subset(Meta_data_plots_GEM_Ghana,Meta_data_plots_GEM_Ghana$ï..GEM.plot.code %in% soil_property_plots_GEM_Ghana$GEM.Plot.Code)

plots_GEM_Ghana<-plots_GEM_Ghana[order(plots_GEM_Ghana$ï..GEM.plot.code),]
plots_GEM_Ghana$rdpth<-pmin(plots_GEM_Ghana$wtd,plots_GEM_Ghana$rdpth)
plots_GEM_Ghana$ups_area[1]<-0

###organize soil for splash
soil_data_GEM_Ghana<-list()
for(i in 1:11){soil_data_GEM_Ghana[[i]]<-c(sand=soil_property_plots_GEM_Ghana$Perc_Sand[i],clay=soil_property_plots_GEM_Ghana$Perc_Clay[i],som=soil_property_plots_GEM_Ghana$C_Perc[i]*1.72,fgravel=0,bd=soil_property_plots_GEM_Ghana$Bulk.density[i],depth=2)}
names(soil_data_GEM_Ghana)<-plots_GEM_Ghana$ï..GEM.plot.code
soil_data_GEM_Ghana$`BOB-01`[6]<-3.5
soil_data_GEM_Ghana$`BOB-02`[6]<-3.5
soil_data_GEM_Ghana$`BOB-04`[6]<-3.5

### stonines values
soil_data_GEM_Ghana$`ANK-01`[4]<-35.94658449
soil_data_GEM_Ghana$`ANK-02`[4]<-35.94658449
soil_data_GEM_Ghana$`ANK-03`[4]<-2.85804574

soil_data_GEM_Ghana$`BOB-02`[4]<-8.690775572
soil_data_GEM_Ghana$`BOB-04`[4]<-0.882750341
soil_data_GEM_Ghana$`BOB-01`[4]<-mean(c(soil_data_GEM_Ghana$`BOB-02`[4]),soil_data_GEM_Ghana$`BOB-04`[4])
soil_data_GEM_Ghana$`KOG-02`[4]<-0.030995778
soil_data_GEM_Ghana$`KOG-05`[4]<-0.122173226
soil_data_GEM_Ghana$`KOG-03`[4]<-mean(c(soil_data_GEM_Ghana$`KOG-02`[4]),soil_data_GEM_Ghana$`KOG-05`[4])
soil_data_GEM_Ghana$`KOG-04`[4]<-soil_data_GEM_Ghana$`KOG-03`[4]
soil_data_GEM_Ghana$`KOG-06`[4]<-soil_data_GEM_Ghana$`KOG-03`[4]

#################################################################################################
### 02. get forcing data
#################################################################################################
getForcing<-function(lat,lon,startDate='1980-01-01', endDate='2015-12-31'){
	###georeference
	site<-SpatialPoints(coords=cbind(x=lon,y=lat),proj4string =CRS("+init=epsg:4326"))
	site<-sf::st_as_sf(site)
	result<-mapply(getTerraClim,param=vars$terraclim$common.name,MoreArgs = list(AOI=site,startDate=startDate, endDate=endDate),SIMPLIFY = F)
	res <- lapply(result,function(df){df[,5]})
	res<-do.call(cbind,res)
	ind<-strptime(paste0(as.character(result$aet$date),'-01'),format='%Y-%m-%d')
	res<-xts(res,as.Date(ind))
	res$tmean<-(res$tmax+res$tmin)/2
	res
	
}

Ghana_forcing<-mapply(getForcing,lat=sites_centroids$Latitude,lon=sites_centroids$Longitude,MoreArgs = list(startDate='1980-01-01', endDate='2015-12-31'),SIMPLIFY = FALSE)
gc()
names(Ghana_forcing)<-sites_centroids$Group.1
save(Ghana_forcing,file="C:/dsval/Huanyuan_soil/sites_Ghana_forcing.RData")
#################################################################################################
### 03. run splash
#################################################################################################

runsplash<-function(df,lat,elev,slop,asp,au,soil_data){
	splash.point(
		sw_in=df$srad,	# shortwave radiation W/m2
		tc=df$tmean,		# air temperature C
		pn= df$prcp,		# precipitation mm
		lat=lat,		# latitude deg
		elev=elev,		# elevation masl
		slop=slop,	# slope deg
		asp=asp,		# aspect deg
		soil_data=soil_data, 		# soil data: sand,clay,som in w/w %. Gravel v/v %, bulk density g/cm3, and depth to the bedrock (m)**
		Au=au,		# upslope area m2
		resolution=250  			# resolution pixel dem used to get Au
	)
}
##################################################################################################
#####Run splash Ankasa
##################################################################################################
splash_Ankasa<-mapply(runsplash,lat=plots_GEM_Ghana$Latitude[plots_GEM_Ghana$Site.Name=='Ankasa'],
					elev=plots_GEM_Ghana$Elevation..m.[plots_GEM_Ghana$Site.Name=='Ankasa'],
					slop=plots_GEM_Ghana$slop[plots_GEM_Ghana$Site.Name=='Ankasa'],
					asp=plots_GEM_Ghana$asp[plots_GEM_Ghana$Site.Name=='Ankasa'],
					au=plots_GEM_Ghana$ups_area[plots_GEM_Ghana$Site.Name=='Ankasa'],
					soil_data=subset(soil_data_GEM_Ghana,names(soil_data_GEM_Ghana)%in% plots_GEM_Ghana$ï..GEM.plot.code[plots_GEM_Ghana$Site.Name=='Ankasa']),
					MoreArgs = list(df=Ghana_forcing$Ankasa),SIMPLIFY = FALSE)
##################################################################################################
#####Run splash Bobiri
##################################################################################################			
splash_Bobiri<-mapply(runsplash,lat=plots_GEM_Ghana$Latitude[plots_GEM_Ghana$Site.Name=='Bobiri'],
	elev=plots_GEM_Ghana$Elevation..m.[plots_GEM_Ghana$Site.Name=='Bobiri'],
	slop=plots_GEM_Ghana$slop[plots_GEM_Ghana$Site.Name=='Bobiri'],
	asp=plots_GEM_Ghana$asp[plots_GEM_Ghana$Site.Name=='Bobiri'],
	au=plots_GEM_Ghana$ups_area[plots_GEM_Ghana$Site.Name=='Bobiri'],
	soil_data=subset(soil_data_GEM_Ghana,names(soil_data_GEM_Ghana)%in% plots_GEM_Ghana$ï..GEM.plot.code[plots_GEM_Ghana$Site.Name=='Bobiri']),
	MoreArgs = list(df=Ghana_forcing$Ankasa),SIMPLIFY = FALSE)
	
##################################################################################################
#####Run splash Kogyae
##################################################################################################	
splash_Kogyae<-mapply(runsplash,lat=plots_GEM_Ghana$Latitude[plots_GEM_Ghana$Site.Name=='Kogyae  '],
	elev=plots_GEM_Ghana$Elevation..m.[plots_GEM_Ghana$Site.Name=='Kogyae  '],
	slop=plots_GEM_Ghana$slop[plots_GEM_Ghana$Site.Name=='Kogyae  '],
	asp=plots_GEM_Ghana$asp[plots_GEM_Ghana$Site.Name=='Kogyae  '],
	au=plots_GEM_Ghana$ups_area[plots_GEM_Ghana$Site.Name=='Kogyae  '],
	soil_data=subset(soil_data_GEM_Ghana,names(soil_data_GEM_Ghana)%in% plots_GEM_Ghana$ï..GEM.plot.code[plots_GEM_Ghana$Site.Name=='Kogyae  ']),
	MoreArgs = list(df=Ghana_forcing$Ankasa),SIMPLIFY = FALSE)
splash_GEM_Ghana<-c(splash_Ankasa,splash_Bobiri,splash_Kogyae)
names(splash_GEM_Ghana)<-plots_GEM_Ghana$ï..GEM.plot.code
# save(splash_GEM_Ghana,file="C:/dsval/Huanyuan_soil/splash_GEM_Ghana.RData")
##################################################################################################
##### 04. get drought splash
##################################################################################################	
getDrought<-function(df,soil_data){
	soil_info<-soil_hydro(sand=soil_data[1],clay=soil_data[2],OM=soil_data[3],fgravel = 0,bd=soil_data[5])
	#swc<-unSWC(soil_data,1,wn=df$wn,units = 'vol')
	swc<-df$wn/(soil_data[6]*1000)
	vwc_mean<-mean(swc[,1])
	vwc_sd<-sd(swc[,1])
	Se_mean<-mean(swc[,1]/soil_info$SAT,na.rm=T)
	Se_sd<-sd(swc[,1]/soil_info$SAT,na.rm=T)
	# fFC_mean<-mean(swc[,1]/soil_info$FC,na.rm=T)
	# fFC_sd<-sd(swc[,1]/soil_info$FC,na.rm=T)
	aet<-apply.yearly(df$aet,sum)
	pet<-apply.yearly(df$pet,sum)
	alpha<-aet/pet
	alpha_mean<-mean(alpha,na.rm=T)
	alpha_sd<-sd(alpha,na.rm=T)
	data.frame(vwc_mean=vwc_mean,vwc_sd=vwc_sd,Se_mean,Se_sd,alpha_mean,alpha_sd,
			aet=mean(aet),pet=mean(pet),ro=mean(apply.yearly(df$ro,sum)),cond=mean(apply.yearly(df$cond,sum)),
			perc=mean(apply.yearly(df$bf,sum)))
	
}

WB_GEM_Ghana<-mapply(getDrought,splash_GEM_Ghana,soil_data_GEM_Ghana,SIMPLIFY = F)
WB_GEM_Ghana<-do.call(rbind,WB_GEM_Ghana)
### calc annual means environmental vars
get_env<-function(df){
	MAP_mean<-mean(apply.yearly(df$prcp,sum))
	MAP_sd<-sd(apply.yearly(df$prcp,sum))
	MAT_mean<-mean(apply.yearly(df$tmean,mean))
	MAT_sd<-sd(apply.yearly(df$tmean,mean))
	VPD_mean<-mean(apply.yearly(df$vpd,mean))
	VPD_sd<-sd(apply.yearly(df$vpd,mean))
	data.frame(MAP_mean,MAP_sd,MAT_mean,MAT_sd,VPD_mean,VPD_sd)
}

env_sites<-lapply(Ghana_forcing,get_env)
env_sites<-do.call(rbind,env_sites)
env_sites$site<-row.names(env_sites)
##########################
#merge it all
plots_GEM_Ghana<-cbind(plots_GEM_Ghana,WB_GEM_Ghana)
plots_GEM_Ghana<-merge(plots_GEM_Ghana,env_sites,by.x='Site.Name',by.y='site')
plots_GEM_Ghana$dryness<-plots_GEM_Ghana$pet/plots_GEM_Ghana$MAP_mean


write.csv(plots_GEM_Ghana,"C:/dsval/Huanyuan_soil/plots_GEM_Ghana_sims_WB.csv")


