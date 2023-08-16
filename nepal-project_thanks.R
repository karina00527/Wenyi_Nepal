#########Nepal Simulations vegetation
library(raster)
library(rpmodel)
library(rsplashtest)
#rasterOptions(todisk = TRUE, overwrite=TRUE, tolerance = 0.5,progress = 'text')
rasterOptions(timer = TRUE, tolerance = 0.5,progress = 'text')
##
###################################################################################################
# 1. load the data 
###################################################################################################
Nepal_wb<-list()
Nepal_wb$sm_lim<-brick("X:/projects/leverhulme_wildfires_life_sciences/live/data/Nepal_project/Nepal_data/SPLASH_v2.0_sm_lim_monthly_2003-2022.nc")
Nepal_wb$aet<-brick("X:/projects/leverhulme_wildfires_life_sciences/live/data/Nepal_project/Nepal_data/SPLASH_v2.0_aet_monthly_2003-2022.nc")
Nepal_wb$pet<-brick("X:/projects/leverhulme_wildfires_life_sciences/live/data/Nepal_project/Nepal_data/SPLASH_v2.0_pet_monthly_2003-2022.nc");gc()

Nepal_GPP_0<-brick("X:/projects/leverhulme_wildfires_life_sciences/live/data/Nepal_project/Nepal_data/pot_2003_2022.gpp.nc");gc()

ind.months<-seq(as.Date('2003-01-01'),as.Date('2022-12-31'), by="month")
start.time<-Sys.time()
##############ancillary##############
elev<-raster("/rds/general/user/ds6915/home/WORK/data_input/global_1km/elevation_1KMmd_GMTED_int.tif")
filenames.soil<- list.files(path="/rds/general/user/ds6915/home/WORK/data_input/global_1km/soil_data", pattern = ".*.tif$", full.names=TRUE)
soil<-stack(filenames.soil);gc()
#################################### forcing ###################################################### 
##### precipitation cru
pn<-brick("Nepal_pr_2003_2022.tif")
pn<-setZ(pn,as.Date(ind.months))
##### air temperature 
tc<-brick("tc_mo.grd")
tc<-setZ(tc,as.Date(ind.months))
##### solar radiation
sw_in<-brick("sw_in_mo.grd")
##### fapar
fapar<-brick("fapar_mo.grd")
####vpd
vpd<-brick("vpd_mo.grd")
################## CO2 ppm
CO2_maunaLoa<-read.table('ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt')
CO2_maunaLoa[CO2_maunaLoa==-99.99]<-NA
names(CO2_maunaLoa)<-c('year','month','decimal_date','average', 'interpolated','trend','#days')
CO2_maunaLoa$yearmonth<-paste0(CO2_maunaLoa$year,'-',sprintf("%02d", CO2_maunaLoa$month))
monthly_in_situ_co2<-subset(CO2_maunaLoa,CO2_maunaLoa$year>=2003 & CO2_maunaLoa$year<2023)
CO2_ave<-as.numeric(tapply(monthly_in_situ_co2$average,monthly_in_situ_co2$month,FUN=mean))
cat("time taken reading")
end.time<-Sys.time()
end.time-start.time
####################################################################################################
# 2. data pre-processing, calculate ppfd [mol/m2/month], vpd[Pa], add time info
####################################################################################################
#beginCluster(40,'SOCK')
start.time<-Sys.time()
#setwd(dirname(rasterTmpFile()))
# define the functions to paralellize
monthlyppfd<-function(x){
  ##factor from integers
  x<-x*0.1
  # w/m2 to mol/m2/day
  # sw_in*secsperday*kfFEC*1e-6*ndaysmonth
  return((x*86400*2.04*1e-6 *30))
}
scale_0d1<-function(x){
  return((x*0.1))
}
scale_0d01<-function(x){
  return((x*0.01))
}
scale_10<-function(x){
  return((x*10))
}
# ############ correct units
# ppfd<-calc(sw_in,fun=monthlyppfd,filename='ppfd_mol_mo.grd')
ppfd<-brick("ppfd_mol_mo.grd")
# fapar<-calc(fapar,fun=scale_0d01,filename='fapar_mo.grd')
# tc<-calc(tc,fun=scale_0d1,filename='tc_mo.grd')
# sw_in<-calc(sw_in,fun=scale_0d1,filename='sw_in_mo.grd')
# vpd<-calc(vpd,fun=scale_10,filename='vpd_mo.grd',overwrite=TRUE)
elev<-crop(elev,tc)
elev<-projectRaster(elev,tc)
soil<-crop(soil,tc)
soil<-projectRaster(soil,tc)
####################################################################################################
# 5. run splash
####################################################################################################
beginCluster(24,'SOCK')
start.time<-Sys.time()
Nepal_wb<-splash.grid(sw_in=sw_in,	# shortwave radiation W/m2
                      tc=tc,	# air temperature C
                      pn=pn,	# precipitation mm/mo
                      elev=elev,# elevation masl
                      soil=soil,# soil data: sand,clay,som in w/w %. Gravel v/v %, bulk density g/cm3, and depth
                      outdir=outf)
endCluster();gc()
#### calc climatological alpha
alpha<-sum(Nepal_wb$aet)/sum(Nepal_wb$pet)

cat("run splash")
end.time-start.time
####################################################################################################
# 5. run p model
####################################################################################################
cat("starting pmodel")
start.time<-Sys.time()
dummy<-elev>1
#gc()
beginCluster(12,'SOCK')
gpp_0<-rpmodel.grid(tc=tc,vpd=vpd,elev = elev,co2=CO2_ave,fapar = fapar,ppfd = ppfd,soilm = Nepal_wb$sm_lim,meanalpha = alpha,outdir=outf)
cat("time taken computing gpp")
end.time<-Sys.time()
end.time-start.time
endCluster();gc()

####################################################################################################
# Nepal biome distribution analysis season
# load libraries
library(raster)
library(tidyr)
library(dplyr)
library(plotly)
### NEPAL reduce NPP0
Nepal_GPP_0<-brick("Y:/Nepalbiome/pot_2003_2022.gpp.nc")
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")
index_agg<-format(ind.months,'%m')
gpp_zeng<-setZ(gpp_zeng, ind.months)
NEPAL_season_gpp0<-zApply(gpp_zeng,index_agg,mean)
NEPAL_season_gpp0_sd<-zApply(gpp_zeng,index_agg,sd)
names(NEPAL_season_gpp0)<-paste0('mean',month.abb)
names(NEPAL_season_gpp0_sd)<-paste0('sd',month.abb)
##### add biomes
Nepal_Biomes<-raster("Y:/Nepalbiome/Nepal_Biomes.tif")
###
Nepal_veg_season<-sampleRandom(stack(Nepal_Biomes,NEPAL_season_gpp0,NEPAL_season_gpp0_sd), ncell(Nepal_Biomes)*.3,xy=TRUE)
Nepal_veg_season<-as.data.frame(Nepal_veg_season)
Nepal_veg_season <- rename(Nepal_veg_season, Nepal_Biomes = category)
#creat attribute data
values <- c(1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)
attribute <- c("ENF", "EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","PWE","CRO",
               "URB","CVM","SNO","BSV","WAT")
igbp_md<- data.frame(value = values, category = attribute)
Nepal_veg_season<-merge(Nepal_veg_season,igbp_md,by.x='Nepal_Biomes',by.y='value')
write.csv(Nepal_veg_season,"Y:/Nepalbiome/Nepal_veg_season.csv")
Nepal_veg_season<-read.csv("Y:/Nepalbiome/Nepal_veg_season.csv")
#########
df_long <- pivot_longer(Nepal_veg_season, 
                        cols =  c(starts_with("mean"), starts_with("sd")),
                        names_to = c(".value", "Month"), 
                        names_pattern = "(\\w+)(\\w{3})")
#############filter zero values
df <- df_long %>%filter(mean != 0)%>%filter(sd!=0)
df<-df %>% filter(category!= "CRO")
df<-df %>% filter(category!= "URB")
df<-df %>% filter(category!= "SNO")
df<-df %>% filter(category!= "BSV")
df<-df %>% filter(category!= "WAT")
df<-df %>% filter(category!= "CVM")
#creat new df
df2 <- df %>%
  group_by(category, Month) %>%
  summarize(mean = mean(mean),
            sd = mean(sd),.groups = "drop")
#cal upper and lower
df3 <- df2 %>%
  mutate(upper = mean + sd,
         lower = mean - sd)
df3 <- df3 %>%
  arrange(match(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
df3$Month <- factor(df3$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), ordered = TRUE)
#########   
fig <-df3%>%
  ggplot(aes(x = Month, y = mean, group = category, color = category))+
  geom_ribbon(aes(x = Month, ymin = lower, ymax = upper, fill =category),linetype = 0,alpha=0.5) +
  geom_line(alpha=1,linewidth = 1)+
  geom_point()+
  labs(x = "Month", y = "GPP", title = "GPP by Biome Type") +
  theme_light()+
  scale_fill_manual(values = c(ENF = "#05450a", DBF = "#78d203", EBF = "#086a10",
                               GRA = "#b6ff05", MF = "#009900",
                               SAV = "#fbff13",OSH="#dcd159",
                               WSA = "#dade48")) +
  scale_color_manual(values = c(ENF = "#05450a", DBF = "#78d203", EBF = "#086a10",
                                GRA = "#b6ff05", MF = "#009900",
                                SAV = "#fbff13",OSH="#dcd159",
                                WSA = "#dade48"))
fig <- ggplotly()
fig
### NEPAL reduce NPP0
Nepal_GPP_2<-brick("Y:/Nepalbiome/2003_2022.gpp.nc")
ind.months<-seq(as.Date('2003-01-01'),as.Date('2022-12-31'), by="month")
index_agg<-format(ind.months,'%m')
NEPAL_season_gpp1<-zApply(Nepal_GPP_1,index_agg,mean)
NEPAL_season_gpp1_sd<-zApply(Nepal_GPP_1,index_agg,sd)
names(NEPAL_season_gpp1)<-paste0('mean',month.abb)
names(NEPAL_season_gpp1_sd)<-paste0('sd',month.abb)
##### add biomes
Nepal_Biomes<-raster("Y:/Nepalbiome/Nepal_Biomes.tif")
###
Nepal_veg_season1<-sampleRandom(stack(Nepal_Biomes,NEPAL_season_gpp1,NEPAL_season_gpp1_sd), ncell(Nepal_Biomes)*.3,xy=TRUE)
Nepal_veg_season1<-as.data.frame(Nepal_veg_season1)
Nepal_veg_season1 <- rename(Nepal_veg_season1, Nepal_Biomes = category)
#creat attribute data
values <- c(1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17)
attribute <- c("ENF", "EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","PWE","CRO",
               "URB","CVM","SNO","BSV","WAT")
igbp_md<- data.frame(value = values, category = attribute)
Nepal_veg_season1<-merge(Nepal_veg_season1,igbp_md,by.x='Nepal_Biomes',by.y='value')
write.csv(Nepal_veg_season1,"Y:/Nepalbiome/Nepal_veg_season1.csv")
Nepal_veg_season<-read.csv("Y:/Nepalbiome/Nepal_veg_season.csv")
#########
df_long1 <- pivot_longer(Nepal_veg_season1, 
                         cols =  c(starts_with("mean"), starts_with("sd")),
                         names_to = c(".value", "Month"), 
                         names_pattern = "(\\w+)(\\w{3})")
#############filter zero values
df1 <- df_long1 %>%filter(mean != 0)%>%filter(sd!=0)
df1<-df1 %>% filter(category!= "URB")
df1<-df1 %>% filter(category!= "SNO")
df1<-df1 %>% filter(category!= "BSV")
df1<-df1 %>% filter(category!= "WAT")
#creat new df
df22 <- df1 %>%
  group_by(category, Month) %>%
  summarize(mean = mean(mean),
            sd = mean(sd),.groups = "drop")
#cal upper and lower
df33 <- df22 %>%
  mutate(upper = mean + sd,
         lower = mean - sd)
df33 <- df33 %>%
  arrange(match(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
df33$Month <- factor(df33$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), ordered = TRUE)
#filter 
write.csv(df33,"Y:/Nepalbiome/df33.csv")
#########   
df33%>%
  ggplot(aes(x = Month, y = mean, group = category, color = category))+
  geom_ribbon(aes(x = Month, ymin = lower, ymax = upper, fill =category),linetype = 0,alpha=0.5) +
  geom_line(alpha=1,linewidth = 1)+
  geom_point()+
  labs(x = "Month", y = "GPP", title = "Actual GPP by Biome Type") +
  theme_light()+
  scale_fill_manual(values = c(ENF = "#05450a", DBF = "#78d203", EBF = "#086a10",
                               GRA = "#b6ff05", MF = "#009900",
                               SAV = "#fbff13",OSH="#dcd159",
                               WSA = "#dade48")) +
  scale_color_manual(values = c(ENF = "#05450a", DBF = "#78d203", EBF = "#086a10",
                                GRA = "#b6ff05", MF = "#009900",
                                SAV = "#fbff13",OSH="#dcd159",
                                WSA = "#dade48"))


####################################################################################################
Nepal_GPP_1<-brick("Y:/Nepalbiome/2003_2022.gpp.nc")
ind.months<-seq(as.Date('2003-01-01'),as.Date('2022-12-31'), by="month")
index_agg<-format(ind.months,'%m')
NEPAL_season_gpp1<-zApply(Nepal_GPP_1,index_agg,mean)
NEPAL_season_gpp1_sd<-zApply(Nepal_GPP_1,index_agg,sd)
names(NEPAL_season_gpp1)<-paste0('mean',month.abb)
names(NEPAL_season_gpp1_sd)<-paste0('sd',month.abb)
##### add biomes
Nepal_Biomes<-raster("Y:/Nepalbiome/Nepal_Biomes.tif")
###
Nepal_veg_season1<-sampleRandom(stack(Nepal_Biomes,NEPAL_season_gpp1,NEPAL_season_gpp1_sd), ncell(Nepal_Biomes)*.3,xy=TRUE)
Nepal_veg_season1<-as.data.frame(Nepal_veg_season1)
Nepal_veg_season1 <- rename(Nepal_veg_season1, Nepal_Biomes = category)
#creat attribute data
values <- c(1, 2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17)
attribute <- c("ENF", "EBF","DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","PWE","CRO",
               "URB","CVM","SNO","BSV","WAT")
igbp_md<- data.frame(value = values, category = attribute)
Nepal_veg_season1<-merge(Nepal_veg_season1,igbp_md,by.x='Nepal_Biomes',by.y='value')
write.csv(Nepal_veg_season1,"Y:/Nepalbiome/Nepal_veg_season1.csv")

#########
df_long1 <- pivot_longer(Nepal_veg_season1, 
                         cols =  c(starts_with("mean"), starts_with("sd")),
                         names_to = c(".value", "Month"), 
                         names_pattern = "(\\w+)(\\w{3})")
#############filter zero values
df1 <- df_long1 %>%filter(mean != 0)%>%filter(sd!=0)
df1<-df1 %>% filter(category!= "URB")
df1<-df1 %>% filter(category!= "SNO")
df1<-df1 %>% filter(category!= "BSV")
df1<-df1 %>% filter(category!= "WAT")
df1<-df1 %>% filter(category!= "OSH")
#creat new df
df22 <- df1 %>%
  group_by(category, Month) %>%
  summarize(mean = mean(mean),
            sd = mean(sd),.groups = "drop")
#cal upper and lower
df33 <- df22 %>%
  mutate(upper = mean + sd,
         lower = mean - sd)
df33 <- df33 %>%
  arrange(match(Month, c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")))
df33$Month <- factor(df33$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), ordered = TRUE)
df33$category <- factor(df33$category, levels = c("ENF","EBF","DBF","MF","WSA","SAV","GRA","CRO","CVM"), ordered = TRUE)
#filter 
write.csv(df33,"Y:/Nepalbiome/df33.csv")
##############
#mk trend for gpp
library(raster)
library(tseries)
library(forecast)
library(rgdal)
library(RColorBrewer)
library(stats)
library(base)
library(ggplot2)
library(Kendall)
library(trend)
Nepal_GPP<-brick("Y:/Nepalbiome/2003_2022.gpp.nc")
nepalgpp<-Nepal_GPP[[1:216]]
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")
nepalgpp<-setZ(nepalgpp, ind.months)
#make month to year
nbands <- nlayers(tc)
nmonths <- nbands / 12 
annual_data <- list()
for (i in 1:nmonths) {
  start_band <- (i - 1) * 12 + 1
  end_band <- i * 12
  annual_data[[i]] <- sum(nepalgpp[[start_band:end_band]], na.rm = TRUE)
}
nepalgpp_annual<- do.call(brick, annual_data)
writeRaster(nepalgpp_annual,"Y:/Nepal_final/hurst/annual.tif")
#pre-monsoon
#monsoon
#postmonsoon
#winter
nbands <- nlayers(nepalgpp)
nmonths <- nbands / 12 
annual_data <- list()
for (i in 1:nmonths) {
  start_band <- (i - 1) * 12 + 1
  end_band <- i * 12
  annual_data[[i]] <- sum(nepalgpp[[start_band:end_band]][[c(1,2,12)]], na.rm = TRUE)
}
nepal_gpp_winter<- do.call(brick, annual_data)
writeRaster(nepal_gpp_winter,"Y:/Nepal_final/hurst/premonsoon.tif")
#mk trend
rast(nepal_gpp_winter)@ptr[['names']] |> 
  length() -> times #用于时间序列MK趋势检验的时间跨度
fun_sen <- function(x){
  
  if (length(na.omit(x)) < 18 || all(x == 0)) return(c(NA, NA, NA))   
  MK_estimate = trend::sens.slope(
    ts(na.omit(x),start=2003,end = 2020,frequency=1)
  ) 
  slope <- MK_estimate$estimates
  MK_test <- MK_estimate$p.value
  Zs <- MK_estimate$statistic
  
  return(c(slope,MK_test,Zs))
}
winter = app(rast(nepal_gpp_winter), fun_sen, cores=4)
writeRaster(winter[[1]],"Y:/Nepal_final/gpp_slope_winter.tif")
writeRaster(winter[[2]],"Y:/Nepal_final/gpp_p_winter.tif")
#cal different area
annual<-raster("Y:/Nepal_final/gpp_slope_new_sum.tif")
winter<-raster("Y:/Nepal_final/gpp_slope_winter.tif")
pre_monsoon<-raster("Y:/Nepal_final/gpp_slope_premonsoon.tif")
monsoon<-raster("Y:/Nepal_final/gpp_slope_monsoon.tif")
post_monsoon<-raster("Y:/Nepal_final/gpp_slope_postmonsoon.tif")
#biome
ENF<-Nepal_Biomes==1
EBF<-Nepal_Biomes==2
DBF<-Nepal_Biomes==4
MF<-Nepal_Biomes==5
WSA<-Nepal_Biomes==8
SAV<-Nepal_Biomes==9
GRA<-Nepal_Biomes==10
CRO<-Nepal_Biomes==12
CVM<-Nepal_Biomes==14
plot(CRO)
# 计算ENF
#annual
annual_enf<-annual*ENF
annual_enf[annual_enf == 0] <- NA
mean_annual_enf <- cellStats(annual_enf, mean)
sd_annual_enf<-cellStats(annual_enf,sd)
#winter
winter_enf<-winter*ENF
winter_enf[winter_enf == 0] <- NA
mean_winter_enf <- cellStats(winter_enf, mean)
sd_winter_enf<-cellStats(winter_enf,sd)
# pre_monsoon
pre_monsoon_enf <- pre_monsoon * ENF
pre_monsoon_enf[pre_monsoon_enf == 0] <- NA
mean_pre_monsoon_enf <- cellStats(pre_monsoon_enf, mean)
sd_pre_monsoon_enf <- cellStats(pre_monsoon_enf, sd)
# monsoon
monsoon_enf <- monsoon * ENF
monsoon_enf[monsoon_enf == 0] <- NA
mean_monsoon_enf <- cellStats(monsoon_enf, mean)
sd_monsoon_enf <- cellStats(monsoon_enf, sd)
# post_monsoon
post_monsoon_enf <- post_monsoon * ENF
post_monsoon_enf[post_monsoon_enf == 0] <- NA
mean_post_monsoon_enf <- cellStats(post_monsoon_enf, mean)
sd_post_monsoon_enf <- cellStats(post_monsoon_enf, sd)

# 创建一个空的 data.frame
result_df <- data.frame(Group = character(), Season = character(), Mean = numeric(), SD = numeric(), stringsAsFactors = FALSE)
result_df <- rbind(result_df, data.frame(Group = "ENF", Season = "Annual", Mean = mean_annual_enf, SD = sd_annual_enf))# 添加 annual_enf 的结果
result_df <- rbind(result_df, data.frame(Group = "ENF", Season = "Winter", Mean = mean_winter_enf, SD = sd_winter_enf))# 添加 winter_enf 的结果
result_df <- rbind(result_df, data.frame(Group = "ENF", Season = "Pre Monsoon", Mean = mean_pre_monsoon_enf, SD = sd_pre_monsoon_enf))
result_df <- rbind(result_df, data.frame(Group = "ENF", Season = "Monsoon", Mean = mean_monsoon_enf, SD = sd_monsoon_enf))# 添加 monsoon_enf 的结果
result_df <- rbind(result_df, data.frame(Group = "ENF", Season = "Post Monsoon", Mean = mean_post_monsoon_enf, SD = sd_post_monsoon_enf))# 添加 post_monsoon_enf 的结果
#figure
library(ggplot2)
df<-result_df
write.csv(df,"Y:/Nepal_final/mk.csv")
df<-read.csv("Y:/Nepal_final/mk.csv")
# 创建 ggplot2 图表，添加误差棒
my_guide <- guide_legend(
  title = element_blank(),
  label = element_text(size = 12, family = "Arial", face = "bold"),
  keywidth = 1.2,
  keyheight = 0.8,
  label.position = "top"
)

#library(ggplot2)
df$Group <- factor(df$Group, levels = c("ENF", "EBF", "DBF", "MF", "WSA", "SAV", "GRA", "CRO", "CVM",ordered = TRUE))
df$Season<-factor(df$Season,levels=c("Annual","Winter","Pre Monsoon","Monsoon","Post Monsoon"))
# Create a vector of colors for the fill scale
colors <- c("#A52A2A", "#00FFFF","#008001", "#2AE41A", "#FFFF00")
# Create the ggplot2 chart with updated styles and parameters
ggplot(df, aes(x = Group, y = Mean, fill = Season)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD),
                width = 0.2, position = position_dodge(0.7),
                color = "black", size = 0.5, show.legend = FALSE) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha=0.9) +
  labs(x="Biome Type",y = "Trend (g C/m2/year)",
       fill = "Season") +
  scale_y_continuous(limits=c(-14,18),breaks = c(-10,-6,-3,0,3,6,10,15))+
  theme_minimal(base_family = "Arial") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 16,family = "Arial",face = "bold"),
        axis.title = element_text(size = 16, face = "bold", family = "Arial"),
        axis.text = element_text(size = 16, family = "Arial", face = "bold"),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = colors[1:length(unique(df$Season))])
#刻度线----------------------------------
gpp<-brick("Y:/Nepal_final/factors/GPP_pmodel_yr.tif")
tc_yr<-brick("Y:/Nepal_final/factors/tc_yr.tif")
ppfd_yr<-brick("Y:/Nepal_final/factors/ppfd_yr.tif")
fapar_yr<-brick("Y:/Nepal_final/factors/fapar_yr.tif")
#########cal


############
#CV
SDpremonsoon<- app(nepal_gpp_premonsoon, sd, cores=4)  #计算标准差
meanpremonsoon<- app(nepal_gpp_premonsoon, mean) 
SDmonsoon<- app(nepal_gpp_monsoon, sd, cores=4)  #计算标准差
meanmonsoon<- app(nepal_gpp_monsoon, mean) 
SDpostmonsoon<- app(nepal_gpp_postmonsoon, sd, cores=4)  #计算标准差
meanpostmonsoon<- app(nepal_gpp_postmonsoon, mean) 
SDwinter<- app(nepal_gpp_winter, sd, cores=4)  #计算标准差
meanwinter<- app(nepal_gpp_winter, mean) 
winter_cv <- SDwinter/meanwinter
monsoon_cv <- SDmonsoon/meanmonsoon
postmonsoon_cv <- SDpostmonsoon/meanpostmonsoon
premonsoon_cv <- SDpremonsoon/meanpremonsoon
##########
#################cal hurst
gpp_ly<-brick("Y:/Nepal_final/hurst/monsoon.tif")
ann<-brick("Y:/Nepal_final/hurst/annual.tif")
library(xts)
library(pracma)
#yearly
rasterOptions(timer = TRUE, tolerance = 0.5,progress = 'text')
######
get_hurst<-function(x){
  x<-as.numeric(x)
  if(anyNA(x)|sum(x==0)==length(x)){
    hrst<-NA
  }else{
    hrst<-as.numeric(pracma::hurstexp(x,display=FALSE)$Hs)
  }
  return(hrst)
}
hurst_premonsoon=calc(x=nepal_gpp_winter,fun=get_hurst)
writeRaster(hurst_premonsoon,"Y:/Nepal_final/hurst/winter2.tif")
###############
##########################
#comparison for zhang
gpp_zeng<-brick("Y:/Nepalbiome/Nepal_GPP_2003_2022.tif")
zeng<-brick("Y:/Nepalbiome/Nepal_GPP_2003_2022.tif")
zeng<-raster::reclassify(zeng, cbind(0, NA))
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")

index_agg<-format(ind.months,'%m')
gpp_zeng<-setZ(zeng, ind.months)
p<-setZ(p,ind.months)
zhang<-zApply(gpp_zeng,index_agg,mean)
p<-zApply(p,index_agg,mean)
names(zhang)<-paste0('mean',month.abb)
names(p)<-paste0("meanp",month.abb)
Nepal_veg_season<-sampleRandom(stack(zhang,p), ncell(p)*.5,xy=TRUE)
Nepal_veg_season<-as.data.frame(Nepal_veg_season)

p<-brick("Y:/Nepalbiome/2003_2022.gpp.nc")
p<-p[[1:216]]
p<-raster::reclassify(p, cbind(0, NA))
#set time
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")
gpp_zeng<-setZ(zeng,as.Date(ind.months))
formatted_years <- format(ind.months, format = "%Y")
gpp_zeng<-zApply(gpp_zeng,formatted_years,sum,na.rm=T)
gpp_zeng<-raster::reclassify(gpp_zeng, cbind(0, NA))
#data for p model
gpp_p<-brick("Y:/Nepal_final/factors/GPP_pmodel_yr.tif")
gpp_p<-raster::reclassify(gpp_p, cbind(0, NA))
##########_______________
p<-mean(gpp_p)
z<-mean(gpp_zeng)
diff_stack <- overlay(p, z, fun=function(t, g) {
  return((t - g) / t * 100)
})
result <- calc(diff_stack, fun=mean, na.rm=TRUE)
result3<-raster("Y:/Nepalbiome/pzre.tif")
plot(result3, breaks=c(-100, seq(-90, 90, by=20), 100), zlim=c(-100, 100), col=terrain.colors(10))

#pick
validation<-sampleRandom(stack(mean(gpp_p),mean(gpp_zeng)),
                         ncell(gpp_zeng)*.3,xy=TRUE)
validation<-as.data.frame(validation)
validation <- rename(validation,zhang = layer.2)
df<-sampleRandom(stack(mean(p),mean(zeng)),
                 ncell(p)*.2,xy=TRUE)
df<-as.data.frame(df)
df<-rename(df,p_model=layer.1)
#make a validation graph
dev.new(width = 6, height = 6)
heatscatter(validation$p_model,validation$zhang,cor=TRUE,ggplot = TRUE,xlim = c(0,3000),ylim=c(0,3000))
heatscatter(df$p_model,df$zeng,cor=TRUE,ggplot=TRUE,xlim = c(0,300),ylim=c(0,300),xlab ="P-model GPP (g C/m2/month)",ylab="Zang-model GPP (g C/m2/month)")
abline(a=6.812397,b=0.985218,col="red",lwd=2)
abline(a=0,b=1,lwd =2,lty = 4)
annotate(geom="text", x=3, y=30, label="Scatter plot",
         color="red")
# cal for remse
rmse<-sqrt(mean((validation$p_model - validation$zhang)^2))
#############
#tc and pn
pn<-raster::reclassify(pn, cbind(0, NA))
tc<-raster::reclassify(tc, cbind(0, NA))
##########
ind.months<-seq(as.Date('2003-01-01'),as.Date('2020-12-31'), by="month")
index_agg<-format(ind.months,'%m')
##############
tc0<-zApply(tc0,index_agg,mean)
pn0<-zApply(pn0,index_agg,mean)

names(tc0)<-paste0('meantc',month.abb)
names(pn0)<-paste0('meanpn',month.abb)

pntc<-sampleRandom(stack(pn0,tc0), ncell(pn0)*.3,xy=TRUE)
pntc<-as.data.frame(pntc)

df_long1 <- pivot_longer(pntc, 
                         cols =  c(starts_with("mean")),
                         names_to = c(".value", "Month"), 
                         names_pattern = "(\\w+)(\\w{3})")

df222 <- df_long1 %>%
  group_by(Month) %>%
  summarize(mean_value = mean(meantc),
            .groups = "drop")

df222<-rename(df222,tcmean=mean_value)
pntc<-merge(df22,df222,by="Month")
pntc$Month <- factor(pntc$Month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"), ordered = TRUE)
############plot
write.csv(pntc,"Y:/Nepal_final/tcpn/pntc.csv")
################
#try for arima(but did used)
#arima?
library(tseries)
library(forecast)
library(raster)
library(rgdal)
library(RColorBrewer)
library(stats)
library(base)
library(ggplot2)

###########
twsa_st <- stack("Y:/Nepalbiome/2003_2022.gpp.nc")
twsa_br <- brick(twsa_st)
nr <-457          ###number of row, requiring modifying
nc <-907          ###number of col, requiring modifying
twsa_data <- as.array(twsa_br)
# 
twsa_raster_s <- twsa_br[[1:204]]    #### number of bands, requiring modifying
twsa_simulated_data <- as.array(twsa_raster_s)
twsa_raster_p <- twsa_br[[1:12]]
twsa_predict_data <- as.array(twsa_raster_p)
# 
# 
for (i in 1:nr) for(j in 1:nc) {
  twsa_ARIMA_data <- as.vector(twsa_data[i,j,])
  twsa_ARIMA_data[twsa_ARIMA_data < 0 ] <- NA
  twsa_ARIMA_data[twsa_ARIMA_data > 100 ] <- NA
  
  if (sum(is.na(twsa_ARIMA_data))<200) {            
    Y <- twsa_ARIMA_data
    ts_data = ts(Y, frequency = 12, start = c(2003,1))        ###the start year, requiring modifying
    data_train = window(ts_data, start = c(2003,1), end = c(2019,12))           ###start and end year, requiring modifying
    arima1 = auto.arima(data_train, trace = TRUE, test = 'kpss', ic = 'bic')
    
    #Use auto.arima to fit the model
    arima1_res = arima1$residual
    arima1_simulated = arima1$fitted
    arima1_predict <- forecast(arima1, h = 12)
    #       
    index = 204-length(arima1_simulated)
    if (index>1) { 
      twsa_simulated_data[i,j,1:index] <- NA
    } else {
      twsa_simulated_data[i,j,1] <- NA  
    }
    twsa_simulated_data[i,j,(index+1):204] <- arima1_simulated     #### number of bands, requiring modifying
    twsa_predict_data[i,j,1:12] <- arima1_predict$mean
    
  } else {
    twsa_simulated_data[i,j,] = NA
    twsa_predict_data[i,j,] = NA
  }
  
  print(i)
  print(j)
  
}

# ############
# 
twsa_raster_s <- setValues(twsa_raster_s,twsa_simulated_data)
writeRaster(x = twsa_raster_s,
            filename="Y:/AVHRR_simulated.tif",
            format = "GTiff", # save as a tif
            datatype="FLT4S", # save as a INTEGER rather than a float
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.

twsa_raster_p <- setValues(twsa_raster_p,twsa_predict_data)
writeRaster(x = twsa_raster_p,
            filename="D:/AVHRR_predict.tif",
            format = "GTiff", # save as a tif
            datatype="FLT4S", # save as a INTEGER rather than a float
            overwrite = TRUE)  # OPTIONAL - be careful. This will OVERWRITE previous files.

# 
# #------------------------------------------------------------------
#partial correlation
GPP_pmodel_yr<-rast("Y:/Nepal_final/factors/GPP_pmodel_yr.tif")
GPP_pmodel_yr[GPP_pmodel_yr==0]<-NA
sm_lim_yr<-rast("Y:/Nepal_final/factors/sm_lim_yr.tif")
CO2_ly<-rast("Y:/Nepal_final/factors/CO2_ly.tif")
pn_yr[pn_yr==0]<-NA
tc_yr[tc_yr==0]<-NA
ppfd_yr[ppfd_yr==0]<-NA
fapar_yr[fapar_yr==0]<-NA
sm_lim_yr[sm_lim_yr==0]<-NA
pn_yr<-rast("Y:/Nepal_final/factors/pn_yr.tif")
tc_yr<-rast("Y:/Nepal_final/factors/tc_yr.tif")
fapar_yr<-rast("Y:/Nepal_final/factors/fapar_yr.tif")
vpd_yr<-rast("Y:/Nepal_final/factors/vpd_yr.tif")
ppfd_yr<-rast("Y:/Nepal_final/factors/ppfd_yr.tif")
library(terra)
library(Hmisc)
z=c(GPP_pmodel_yr,tc_yr,fapar_yr,ppfd_yr,pn_yr,sm_lim_yr,vpd_yr,CO2_ly)
fun_cor=function(x){
  if (any(is.na(x))) {
    return(rep(NA, 7))
  }
  gtc=ppcor::pcor.test(x[1:18],x[19:36],list(x[37:54],x[55:72],x[73:90],x[91:108],x[109:126],x[127:144]),method="pearson")
  gfa=ppcor::pcor.test(x[1:18],x[37:54],list(x[19:36],x[55:72],x[73:90],x[91:108],x[109:126],x[127:144]),method="pearson")
  gppfd=ppcor::pcor.test(x[1:18],x[55:72],list(x[19:36],x[37:54],x[73:90],x[91:108],x[109:126],x[127:144]),method="pearson")
  gpn=ppcor::pcor.test(x[1:18],x[73:90],list(x[19:36],x[37:54],x[55:72],x[91:108],x[109:126],x[127:144]),method="pearson")
  gsm=ppcor::pcor.test(x[1:18],x[91:108],list(x[19:36],x[37:54],x[55:72],x[73:90],x[109:126],x[127:144]),method="pearson")
  gvpd=ppcor::pcor.test(x[1:18],x[109:126],list(x[19:36],x[37:54],x[55:72],x[73:90],x[91:108],x[127:144]),method="pearson")
  gco2=ppcor::pcor.test(x[1:18],x[127:144],list(x[19:36],x[37:54],x[55:72],x[73:90],x[91:108],x[109:126]),method="pearson")
  gtcr=gtc$estimate
  #tct=(sqrt(18-2-2)*gtcr)/sprt(1-gtcr^2)
  gfar=gfa$estimate
  #gfat=(sqrt(18-2-2)*gfar)/sprt(1-gfar^2)
  gppfdr=gppfd$estimate
  #gppfdt=(sqrt(18-2-2)*gppfdr)/sprt(1-gppfdr^2)
  gpnr=gpn$estimate
  #gpnt=(sqrt(18-2-2)*gpnr)/sprt(1-gpnr^2)
  gsmr=gsm$estimate
  #gsmt=(sqrt(18-2-2)*gsmr)/sprt(1-gsmr^2)
  gvpdr=gvpd$estimate
  gco2r=gco2$estimate
  #gvpdt=(sqrt(18-2-2)*gvpdr)/sprt(1-gvpdr^2)
  return(c(gtcr,gfar,gppfdr,gpnr,gsmr,gvpdr,gco2r))
}

r_gpp2=app(z,fun_cor,cores=4)
writeRaster(r_gpp$lyr.1,"Y:/Nepal_final/factor2/tccor.tif")
writeRaster(r_gpp$lyr.2,"Y:/Nepal_final/factor2/fapcor.tif")
writeRaster(r_gpp$lyr.3,"Y:/Nepal_final/factor2/ppfdcor2.tif")
writeRaster(r_gpp$lyr.4,"Y:/Nepal_final/factor2/pncor.tif")
writeRaster(r_gpp$lyr.5,"Y:/Nepal_final/factor2/smcor.tif")
writeRaster(r_gpp$lyr.6,"Y:/Nepal_final/factor2/vpdcor.tif")
writeRaster(r_gpp2$lyr.7,"Y:/Nepal_final/factor2/co2cor.tif")








