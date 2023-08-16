####factors
##most important predictor fapar hpc
#2023/6/21
library(raster)
rasterOptions(timer = TRUE, tolerance = 0.5,progress = 'text')
outf<-"Y:/Nepal_final/factors"
#tmpdir<-"/rds/general/user/ds6915/ephemeral/splash_test"
setwd(outf)

## "/rds/general/user/ds6915\projects\leverhulme_wildfires_life_sciences\live\data\Nepal_project\Nepal_data\simulations_terraclimate\2003_2022.gpp+sm_lim.nc"
###################################################################################################
# 1. load the data 
###################################################################################################
GPP_pmodel<-brick("Y:/Nepalbiome/2003_2022.gpp.nc")
GPP_pmodel<-GPP_pmodel[[1:216]]
GPP_pmodel<-setZ(GPP_pmodel,as.Date(ind.months))
##water balance splash
Nepal_wb<-list()
Nepal_wb$sm_lim<-brick("Y:/Nepalbiome/SPLASH_v2.0_sm_lim_monthly_2003-2022.nc")
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")
Nepal_wb$sm_lim<-Nepal_wb$sm_lim[[1:216]]
Nepal_wb$sm_lim<-setZ(Nepal_wb$sm_lim,as.Date(ind.months))
start.time<-Sys.time()
#################################### forcing ###################################################### 
##biomes
#Nepal_Biomes<-raster("C:/Users/Silwood/OneDrive - Imperial College London/Documents/papers/7_Nepal/Nepal_Biomes.tif")
##### precipitation cru
pn<-brick("Y:/Nepal_final/Nepal_pr_2003_2022.tif")
pn<-pn[[1:216]]
pn0<-setZ(pn,as.Date(ind.months))
#pn<-mask(pn,Nepal_Biomes)
##### air temperature 
tc<-brick("Y:/Nepal/r_scaled_tmean.tif")
tc<-tc[[1:216]]
tc0<-setZ(tc,as.Date(ind.months))
########...............
##### solar radiation
sw_in<-brick("Y:/Nepal_final/Nepal_srad_2003_2022.tif")
scale_0d1<-function(x){
  return((x*0.1))
}
sw_in<-calc(sw_in,fun=scale_0d1)
writeRaster(sw_in,"Y:/Nepal/r_scaled_sw_in.tif")
sw_in<-sw_in[[1:216]]
sw_in<-setZ(sw_in,as.Date(ind.months))
##### fapar
fapar<-brick("Y:/Nepal/r_scaled_fapar.tif")
fapar<-fapar[[1:216]]
fapar<-setZ(fapar,as.Date(ind.months))
####vpd
vpd<-brick("Y:/Nepal/r_scaled_vpd .tif")
vpd<-vpd[[1:216]]
vpd<-setZ(vpd,as.Date(ind.months))
##ppfd
ppfd<-brick("Y:/Nepal/ppfd.tif")
ppfd<-ppfd[[1:216]]
ppfd<-setZ(ppfd,as.Date(ind.months))


################## CO2 ppm
CO2_maunaLoa<-read.table('ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt')
CO2_maunaLoa[CO2_maunaLoa==-99.99]<-NA
names(CO2_maunaLoa)<-c('year','month','decimal_date','average', 'interpolated','trend','#days')
CO2_maunaLoa$yearmonth<-paste0(CO2_maunaLoa$year,'-',sprintf("%02d", CO2_maunaLoa$month))
monthly_in_situ_co2<-subset(CO2_maunaLoa,CO2_maunaLoa$year>=2003 & CO2_maunaLoa$year<2021)
CO2_yr<-as.numeric(tapply(monthly_in_situ_co2$average,monthly_in_situ_co2$year,FUN=mean))
cat("time taken reading")
end.time<-Sys.time()
end.time-start.time
##biomes
#Nepal_Biomes<-raster("C:/Users/Silwood/OneDrive - Imperial College London/Documents/papers/7_Nepal/Nepal_Biomes.tif")
#############################################################################
#### get annual values
#############################################################################
formatted_years <- format(ind.months, format = "%Y")
########______________
GPP_pmodel_yr<-zApply(GPP_pmodel,formatted_years,sum,na.rm=T)
sm_lim_yr<-zApply(Nepal_wb$sm_lim,formatted_years,mean,na.rm=T)
pn_yr<-zApply(pn,formatted_years,sum,na.rm=T)
tc_yr<-zApply(tc,formatted_years,mean,na.rm=T)
fapar_yr<-zApply(fapar,formatted_years,mean,na.rm=T)
vpd_yr<-zApply(vpd,formatted_years,mean,na.rm=T)
ppfd_yr<-zApply(ppfd,formatted_years,sum,na.rm=T);gc()
#############

writeRaster(,"Y:/Nepal_final/factors/GPP_pmodel_yr.tif")
GPP_pmodel_yr<-brick("Y:/Nepal_final/factors/GPP_pmodel_yr.tif")
sm_lim_yr<-brick("Y:/Nepal_final/factors/sm_lim_yr.tif")
pn_yr<-brick("Y:/Nepal_final/factors/pn_yr.tif")
tc_yr<-brick("Y:/Nepal_final/factors/tc_yr.tif")
fapar_yr<-brick("Y:/Nepal_final/factors/fapar_yr.tif")
vpd_yr<-brick("Y:/Nepal_final/factors/vpd_yr.tif")
ppfd_yr<-brick("Y:/Nepal_final/factors/ppfd_yr.tif")
#############################################################################
#### 2. create the functions
#############################################################################
categories<-data.frame(var=c("tc","ppfd","vpd","sm_lim","co2","fapar","pn" ,"Residuals"),value=1:8)
#################
## percentace explained
get_var_exp<-function(lmod){
  # aov
  lmod_aov<- anova(lmod)
  #percentage variance explained
  lmod_aov$PctExp<-lmod_aov$`Sum Sq`/sum(lmod_aov$`Sum Sq`)*100;lmod_aov$PctExp<-round(lmod_aov$PctExp,2)
  #lmod_aov<-lmod_aov[order(lmod_aov$PctExp,decreasing =T),]
  lmod_aov
}
#############################################################################
#### 3. create the empty rasters
#############################################################################

perc_expl_tc<-raster(tc_yr[[1]]*NA)
perc_expl_vpd<-perc_expl_tc
perc_expl_ppfd<-perc_expl_tc
perc_expl_sm_lim<-perc_expl_tc
perc_expl_co2<-perc_expl_tc
perc_expl_fapar<-perc_expl_tc
perc_expl_pn<-perc_expl_tc
most_important<-perc_expl_tc





#############################################################################
#### 4. loop through cells
# #############################################################################
# for(i in 1:ncell(perc_expl_tc)){
# 	##extract the values from cells
# 	#i=150000
# 	tc_temp<-extract(tc,i)
# 	ppfd_temp<-extract(ppfd,i)
# 	vpd_temp<-extract(vpd,i)
# 	fAPAR_temp<-extract(fAPAR,i)
# 	vcmax25_temp<-extract(vcmax25_star,i)
# 	#######################build the data frame
# 	df<-data.frame(tc=as.numeric(tc_temp),
# 		ppfd=as.numeric(ppfd_temp),
# 		vpd=as.numeric(vpd_temp),
# 		fapar=as.numeric(fAPAR_temp),
# 		vcmax=as.numeric(vcmax25_temp),
# 		co2=CO2_f_1983$CO2			
# 	)
# 	
# 	if(anyNA(df)==T){
# 		perc_expl_tc<-setValues(perc_expl_tc,values=NA,index=i)
# 		perc_expl_vpd<-setValues(perc_expl_vpd,values=NA,index=i)
# 		perc_expl_ppfd<-setValues(perc_expl_ppfd,values=NA,index=i)
# 		perc_expl_fapar<-setValues(perc_expl_fapar,values=NA,index=i)
# 		perc_expl_co2<-setValues(perc_expl_co2,values=NA,index=i)
# 		most_important<-setValues(most_important,values=NA,index=i)
# 	}else{
# 		##fit the model
# 		fit <- lm(vcmax~.,data = df)
# 		###get the variance explained
# 		varexp<-get_var_exp(fit)
# 		varexp$var<-row.names(varexp)
# 		##get the most important
# 		most_impor<-varexp[varexp$PctExp==max(varexp$PctExp),]
# 		most_impor<-merge(most_impor,categories,by='var')
# 		####populate the raster
# 		perc_expl_tc<-setValues(perc_expl_tc,values=varexp['tc','PctExp'],index=i)
# 		perc_expl_vpd<-setValues(perc_expl_vpd,values=varexp['vpd','PctExp'],index=i)
# 		perc_expl_ppfd<-setValues(perc_expl_ppfd,values=varexp['ppfd','PctExp'],index=i)
# 		perc_expl_fapar<-setValues(perc_expl_fapar,values=varexp['fapar','PctExp'],index=i)
# 		perc_expl_co2<-setValues(perc_expl_co2,values=varexp['co2','PctExp'],index=i)
# 		most_important<-setValues(most_important,values=most_impor$value,index=i)
# 	}	
# 	
# 
# 	
# }

#############################################################################
#### 4. parallel version
#############################################################################		


#ncores<-parallel::detectCores()-1
install.packages("doSNOW")
library(doSNOW)
cl <- makeCluster(10, type='SOCK')
cl
registerDoSNOW(cl)
niter <- ncell(perc_expl_tc)
pb <- txtProgressBar(max = niter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)	

#most_impor_vec<-foreach (i = 1:niter,.packages = "raster",
#                         .combine=c,.multicombine=TRUE, 
 #                        .inorder=TRUE,
 #                        .export=ls(envir=globalenv()),
 #                        .options.snow = opts) %dopar% {
most_impor_vec <- foreach (i = 1:niter, .packages = "raster",
                           .combine=c, .multicombine=TRUE, 
                           .inorder=TRUE,
                           .export=ls(envir=globalenv()),
                           .options.snow = opts) %dopar% {
  tc_temp<-extract(tc_yr,i)
  gpp_temp<-extract(GPP_pmodel_yr,i)
  pn_temp<-extract(pn_yr,i)
  ppfd_temp<-extract(ppfd_yr,i)
  vpd_temp<-extract(vpd_yr,i)
  fAPAR_temp<-extract(fapar_yr,i)
  sm_lim_temp<-extract(sm_lim_yr,i)
  #######################build the data frame
  df<-data.frame(tc=as.numeric(tc_temp),
                 gpp=as.numeric(gpp_temp),
                 pn=as.numeric(pn_temp),
                 ppfd=as.numeric(ppfd_temp),
                 vpd=as.numeric(vpd_temp),
                 fapar=as.numeric(fAPAR_temp),
                 sm_lim=as.numeric(sm_lim_temp),
                 co2=CO2_yr			
  )
  # if(anyNA(df)==T){
  # 	# perc_expl_tc<-setValues(perc_expl_tc,values=NA,index=i)
  # 	# perc_expl_vpd<-setValues(perc_expl_vpd,values=NA,index=i)
  # 	# perc_expl_ppfd<-setValues(perc_expl_ppfd,values=NA,index=i)
  # 	# perc_expl_fapar<-setValues(perc_expl_fapar,values=NA,index=i)
  # 	# perc_expl_co2<-setValues(perc_expl_co2,values=NA,index=i)
  # 	# most_important<-setValues(most_important,values=NA,index=i)
  # 	return(NA)
  # }else{most_impor_vec

  # 	##fit the model
  # 	fit <- lm(fapar~.,data = df)
  # 	###get the variance explained
  # 	varexp<-get_var_exp(fit)
  # 	varexp$var<-row.names(varexp)
  # 	##get the most important
  # 	most_impor<-varexp[varexp$PctExp==max(varexp$PctExp),]
  # 	most_impor<-merge(most_impor,categories,by='var')
  # 	####populate the raster
  # 	# perc_expl_tc<-setValues(perc_expl_tc,values=varexp['tc','PctExp'],index=i)
  # 	# perc_expl_vpd<-setValues(perc_expl_vpd,values=varexp['vpd','PctExp'],index=i)
  # 	# perc_expl_ppfd<-setValues(perc_expl_ppfd,values=varexp['ppfd','PctExp'],index=i)
  # 	# perc_expl_fapar<-setValues(perc_expl_fapar,values=varexp['fapar','PctExp'],index=i)
  # 	# perc_expl_co2<-setValues(perc_expl_co2,values=varexp['co2','PctExp'],index=i)
  # 	
  # 	return(most_impor$value[1])
  # }	
  if(all(is.na(df$gpp))){
    return(NA)
  }else{
    clcheck<-try(fit <- lm(gpp~.,data = df), silent=TRUE)
    if(class(clcheck)=="try-error"){
     return(list(var = NA, coef = NA))
    } else {
      varexp <- get_var_exp(fit)
      most_impor_var <- row.names(varexp)[which.max(varexp$PctExp)]
      return(list( coef = coef(fit)[most_impor_var]))
   
#  if(all(is.na(df$gpp))){
#    return(NA)var = most_impor_var,
#  } else {
#    clcheck <- try(fit <- lm(gpp ~ ., data = df), silent = TRUE)
 #   if(class(clcheck) == "try-error"){
  #    return(list(var = NA, coef = NA))
  #  } else {
  #    varexp <- get_var_exp(fit)
   #   most_impor_var <- row.names(varexp)[which.max(varexp$PctExp)]
    #  return(coef(fit)[most_impor_var]) # 直接返回最重要预测变量的系数
    }
  }
}
    #if(class(clcheck)=="try-error"){
      #return(NA)
    #}else{
      #fit <- lm(fapar~.,data = df)
      
      #varexp<-get_var_exp(fit)
      #varexp$var<-row.names(varexp)
      
     #most_impor<-varexp[varexp$PctExp==max(varexp$PctExp),]
      #most_impor<-merge(most_impor,categories,by='var')
      #return(most_impor$value[1])	
   # }
    

results <- foreach(i = 1:niter, ...) %dopar% {...}
most_important_var_raster <- setValues(perc_expl_tc, sapply(results, function(x) x$var))
coef_raster <- setValues(perc_expl_tc, most_impor_vec)
coef_raster <- setValues(perc_expl_tc, sapply(results, function(x) x$coef))
###############
most_important<-setValues(most_important,values=most_impor_vec)

stopCluster(cl)
gc()

most_important@data@isfactor<-TRUE
categories$ID<-categories$value
most_important@data@attributes[[1]] <- data.frame(ID=categories$ID,var=categories$var)
writeRaster(most_important,"Y:/Nepal_final/factors/most_important.tif")
most<-raster("Y:/Nepal_final/factors/most_important.tif")
Nepal_Biomes<-raster("Y:/Nepalbiome/Nepal_Biomes.tif")
most_important<-mask(most_important,Nepal_Biomes)

rasterVis::levelplot(most)
plot(most)
