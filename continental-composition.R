## Zircon project 

## Authors: Clement Bataille and Amy Willis

## This code accompanies our paper:
## "Continental igneous rock composition: A major control of past global chemical weathering"
## Clément P. Bataille, Amy Willis, Xiao Yang and Xiao-Ming Liu
## Science Advances 2017

## Please feel free to contact either Clement or
## Amy with any questions

#set-up working directory
#setwd("C:/Users/clement/Dropbox/Amy")
#setwd("/Users/adw96/Dropbox/Zircon function CB AW")

#set-up library packages
library(RcmdrMisc);library(data.table);library(zoo);library(dplyr);
library(mgcv);library(pracma);library(DataCombine);library(rootSolve)
library(BB);library(TTR);library(dlnm);library(sampling);
library(ggplot2);library(fields);library(geosphere);library(Hmisc);library(zoo)
library(bigsplines); library(plyr)

#######################STEP 1################################
# Debiasing the detrital zircon dataset

# Convert degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1);   delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = R * c
  return(d) # Distance in km
}

###function for geographic debiasing as in Keller
proximity<-function(long1, lat1, long2, lat2, age1, age2) {
  delta.long <- (long2 - long1);   delta.lat <- (lat2 - lat1); 
  k=(1/((acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(delta.long)))^2+1))
  return (k)
}

##input datasets
options( stringsAsFactors=TRUE)
require(readxl)

# Different zircon dataset inputs
# Zircon_input_71516_cm.xlsx=only Cenozoic and Mesozoic data
#Zircon_input_63016_cmep.xlsx=Cenozoic, Mesozoic, early Paleozoic data
#Zircon_input_63016_cmep.xls=All Phanerozoic age data
#Zircon_input_63016_all.xlsx= All depositional age data

zircon <- readxl::read_excel("New_dataset_6_2_2016/Database S1_1312017.xlsx",col_names=TRUE, na="NaN", sheet="Screening_results_P+NP")
#zircon <- readxl::read_excel("New_dataset_6_2_2016/Zircon_input_63016_cmep.xlsx",col_names=TRUE, na="NaN", sheet="Sheet1")

#classifying RT age by U_Pb age and calculating RT median for each 1myr U-Pb
zircon<-zircon[order(zircon$U_Pb),]
head(zircon)

# turn number column names to usual "X" convention
names(zircon)[which(!is.na(as.numeric(substr(names(zircon),1,1))))] <- paste("X", names(zircon)[which(!is.na(as.numeric(substr(names(zircon),1,1))))], sep="")

#input main measured zircon data age=crystallization age and Hf176/Hf178 is the hafnium isotope ratio measured
Hf176_Hf177<-zircon$X176Hf_177Hf_normalized
age<-zircon$U_Pb

#Caluclate the change in initial Hf176/Hf178 using the radiogenic equation lambda 2 is the decay constant
lambda2=1.867*10^-11 #range of value in the litterature 1.867*10^-11 and 1.93*10^-11 (Blichert-Toft and Albar?de 1997) 
Lu176_Hf177<-zircon$X176Lu_177Hf
Hf176_Hf177i<-Hf176_Hf177-Lu176_Hf177*(exp(lambda2*age*10^6)-1)

#Caluclate eHf a normalized value relative to chondrite
I<-Hf176_Hf177
Ii<-Hf176_Hf177i
#caluclating EHf
ICHUR=0.282772 #?0.00011 range of value in the litterature for this parameter 0.282785 (Bouvier et al. 2008) and 0.28277 (Blichert-Toft and Albar?de 1997) 
RCHUR=0.0332 #?0.0001 range of value some uncertinaty to consider +1 other value in the litterature 0.0336 Bouvier et al. 2008 0.0332 Blichert-Toft and Albar?de (1997) 
ICHUR_t=ICHUR-RCHUR*(exp(lambda2*age*10^6)-1)
eHf_t=((Ii/ICHUR_t)-1)*10000
zircon$eHf_t<-eHf_t

## Debiasing procedure
#########################################################
step<-seq(0,4500,1)
step1<-seq(0,600,10)





complete <- !is.na(zircon$latitude) & !is.na(zircon$longitude) 
zircon_complete <- zircon[complete, ]

###Sampling debiasing
#count the number of zircon grains per sediment
count <-by(zircon_complete$Sediment_ID, as.factor(zircon_complete$Sediment_ID), {function (x) sum(x)/(mean(x))})
#count <-by(zircon_complete$Sediment_ID, as.factor(zircon_complete$Sediment_ID), length)
zircon_complete$count <- count[as.factor(zircon_complete$Sediment_ID)]
zircon_complete<-subset(zircon_complete,count>10)

count2 <-by(zircon_complete$Sediment_ID, as.factor(zircon_complete$Sediment_ID), length)
zircon_complete$count2 <- count2[as.factor(zircon_complete$Sediment_ID)]
zircon_complete<-subset(zircon_complete,count2>10)
sampling <-by(zircon_complete$`count2`, as.factor(zircon_complete$Sediment_ID), min)

#Sampling weight can be change to 1/n or 1/n2 or any other possibility
weight_sampling<-1/(sampling)

###watershed size and reprentativity debiasing (not included in the science advance paper)
### first consider within sand variability
#x11()
sand_variances <- by(zircon_complete$`U_Pb`, as.factor(zircon_complete$Sediment_ID), var, na.rm=TRUE)
sand_variances[which(is.na(sand_variances))] <- 0 ## variance of groups with 1 measurement is "NA". change to 0
min(sand_variances[sand_variances>0], na.rm=T)
sand_variances[sand_variances==0] <- min(sand_variances[sand_variances>0], na.rm=T)
head(sand_variances, 20)
hist(sand_variances, step1*3*10^3)
#sand_variances <- age_variance[as.factor(zircon_complete$Sediment_ID)]
weight_sand <- sand_variances


###Geographic debiasing
#convert latitude and longitude in radians
zircon_complete$longitude_rad<-deg2rad(zircon_complete$longitude)
zircon_complete$latitude_rad<-deg2rad(zircon_complete$latitude)

latitude <- by(zircon_complete$`latitude_rad`, as.factor(zircon_complete$Sediment_ID), median, na.rm=TRUE)
longitude <- by(zircon_complete$`longitude_rad`, as.factor(zircon_complete$Sediment_ID), median, na.rm=TRUE)

#### below here

# next consider geographic variability
## find calibration weights
## find calibration weights
#calibrate_weight_distance <- function(my_longitude, my_latitude) {
#  spatial_difference <- mapply(gcd.hf, 
#                               long1=longitude, lat1=latitude,
#                               long2=my_longitude, lat2=my_latitude)
#  my_means <- mean(spatial_difference)
#my_pie <- 1/var(spatial_difference)
#  my_pie <- sum(1/((spatial_difference)^2+1))
#  if (is.na(my_pie)) print("Uh-oh, missing value?")
#  return(my_pie)
#  return(mean(spatial_difference^2))
#  return(1/var(spatial_difference))
#}

calibrate_weight <- function(my_longitude, my_latitude) {
  difference <- mapply(proximity, 
                       long1=longitude, lat1=latitude,
                       long2=my_longitude, lat2=my_latitude)
  my_pie<-sum(difference)
  return(my_pie)  
}

prox <- mapply(calibrate_weight, my_longitude=longitude, my_latitude=latitude)
prox[is.na(prox)] <- max(prox, na.rm=TRUE)
write.csv(x = prox, file = "prox_loc.csv")
prox <- c(read.csv("prox_loc.csv", row.names = NULL)[,2])

zz <- function(x) {
  z_scores <- scale(x) # subtracts the mean and divides by the standard deviation
  z_scores*50 + 100 # rescales to have mean 100 and standard deviation 10
}

#weight_location <- mapply(calibrate_weight_distance, my_longitude=longitude, my_latitude=latitude)
weight_location=zz(1/prox)
hist(weight_location,breaks=step1*5*10^-10)
plot(longitude,latitude, col=rbPal(10)[as.numeric(cut(weight_location,breaks = 10))], pch=16)
weight_location <- weight_location[as.factor(zircon_complete$Sediment_ID)]

plot(zircon_complete$longitude,zircon_complete$latitude, col=rbPal(10)[as.numeric(cut(weight_location,breaks = 10))], pch=16)
zircon_complete$weight_location<-weight_location

weight_location <- by(zircon_complete$weight_location, as.factor(zircon_complete$Sediment_ID), mean, na.rm=TRUE)


###Calculate the total weight for each sample
###weights_equal assume equal weights for each bias
###weights_more_loc assume more weight for location
###weights_more_sampling assume more weight for sampling

weights_equal <- zz(weight_location) + zz(weight_sand)+zz(weight_sampling)
weights_more_loc <- zz(weight_location)+zz(weight_sampling)
weights_more_sampling <- zz(weight_sampling)

dev.off()

###Generate some maps to vizualize weights distribution
###Map 1 to vizualize weight distribution after location + sampling debiasing
par(mfrow=c(1,1))
rbPal <- colorRampPalette(c('blue','white','red'))
map("world", fill=TRUE, col="white", bg="lightblue", xlim=c(-180,180), ylim=c(-90, 90), mar=c(0,0,0,0))
points(longitude*180/pi,latitude*180/pi, col=rbPal(10)[as.numeric(cut(weights_more_loc,breaks = 10))], pch=16)
#points(longitude2*180/pi,latitude2*180/pi, col="black", pch=10)
#points(longitude*180/pi,latitude*180/pi, col="red", pch=10)


id<-as.numeric(rownames(weights_equal))
weights_final<-weights_equal[,1]
weights_equal<-data.frame(id, weights_final)
##### above here

###Some issues with dimension of the vector in that part of the script
#zircon_complete<-zircon_complete[order(zircon_complete$Sediment_ID),]
#zircon_complete <- merge(zircon_complete, weights_equal, by.x= 'Sediment_ID', by.y='id')
#plot(zircon_complete$longitude,zircon_complete$latitude, col=rbPal(10)[as.numeric(cut(zircon_complete$weights_final,breaks = 10))], pch=16)

#par(mfrow=c(2,3))
#plot(weights_more_loc[as.factor(zircon_complete$Sediment_ID)], weights_equal, xlab="Geographic Weights", ylab="Combined Weights") ## more central => less weight, good
#plot(weight_sampling[as.factor(zircon_complete$Sediment_ID)], weights_equal, xlab="Representativity Weights", ylab="Combined Weights") ## looks good
#plot(zz(weight_sampling)[as.factor(zircon_complete$Sediment_ID)], weights_equal, xlab="Sampling Weights", ylab="Combined Weights") ## looks good

#### 
zircon_complete<-zircon_complete[order(zircon_complete$U_Pb),]
## output
#source("/Users/adw96/Dropbox/Zircon function CB AW/amy_stats_bin_weighted.R")
source("C:/Users/clement/Dropbox/Amy/clem_stats_bin_weighted.R")
#pdf("E:/clement/Other projects/Strontium/R_work/R_code/figures_3_cm.pdf", height=12, width=12)
#pdf("figures_17_june_mc.pdf", height=12, width=12)
set.seed(3)
my_size <- 1e6
my_sample_indices_unw <- sample(dim(zircon_complete)[1], size = my_size, replace = T)
my_sample_indices_equal <- sample(dim(zircon_complete)[1], prob = zircon_complete$weights_final, size = my_size, replace = T)
my_sample_indices_loc <- sample(dim(zircon_complete)[1], prob = zircon_complete$weight_location, size = my_size, replace = T)
#my_sample_indices_sampling <- sample(dim(zircon_complete)[1], prob = weights_more_sampling, size = my_size, replace = T)
eHf_resampled_equal<-extract_resampled_dataset(zircon_complete$U_Pb, zircon_complete$eHf_t, zircon_complete$latitude, zircon_complete$longitude, zircon_complete$Sediment_ID, my_sample_indices_equal)
eHf_resampled_loc<-extract_resampled_dataset(zircon_complete$U_Pb, zircon_complete$eHf_t, zircon_complete$latitude, zircon_complete$longitude, zircon_complete$Sediment_ID, my_sample_indices_loc)
par(mfrow=c(1,1))
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 1, my_sample_indices_unw, my_main="Unweighted")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 1, my_sample_indices_loc, my_main="Location and sand weighted equally")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 1, my_sample_indices_sampling, my_main="More weight on sampling")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 1, my_sample_indices_equal, my_main="More weight on location")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 2, my_sample_indices_unw, my_main="Unweighted")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 2, my_sample_indices_loc, my_main="Location and sand weighted equally")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 2, my_sample_indices_sampling, my_main="More weight on sampling")
#draw_figure(zircon_complete$U_Pb, zircon_complete$eHf_t, 2, my_sample_indices_equal, my_main="More weight on location")
dev.off()

###eHf_resampled_equal accoutn for sampling, representativity and geographic debiasing
#write.csv(x = eHf_resampled_equal, file = "C:/Users/clement/Dropbox/Amy/EHf_resampled_mc_final.csv")
###eHf_resampled_loc accoutn for sampling, and geographic debiasing
write.csv(x = eHf_resampled_loc, file = "C:/Users/clement/Dropbox/Amy/EHf_resampled_mc_final_loc.csv")

###Histograms showing the results of resampling on the final dataset
par(mfrow=c(2,2))
step1<-seq(-90,90,5)
hist(zircon_complete$latitude, breaks=step1, xlab="latitude (raw)")
hist(eHf_resampled_loc$latitude, breaks=step1, xlab="latitude (resampled)")

step2<-seq(-180,180,10)
hist(zircon_complete$longitude, breaks=step2, xlab="latitude (raw)")
hist(eHf_resampled_loc$longitude, breaks=step2, xlab="latitude (resampled)")





########STEP 2################
# Sensitivity analysis of eHf secular trends to uncertainty in U-Pb and Hf isotopes analysis using Gaussian simulation

###eHf input can change to test the sensitivity to data inputs and debiasing choices for instance:
# "EHf_resampled_mc_final_loc.csv"=sampling (1/n) and geographic bias accounted for and only Cenozoic and Mesozoic data
# "EHf_resampled_mc_final_all.csv"=sampling (1/n) and geographic bias accounted for and all data
# "EHf_resampled_mcep_final_loc.csv"= sampling (1/n) and geographic bias and only mesozoic cenozoic and early paleozoic data
# "EHf_resampled_mcep_final_loc_n2.csv"= sampling (1/n2) and geographic bias and only mesozoic cenozoic and early paleozoic data



####### Now permute the age values


## set-up plot
xlim = c(0, 1000)
ylim=c(-35,20)
mat = matrix(1:2, ncol = 1)

#pdf(file = "yx_file/Fig1.pdf", width = 8, height = 9)
layout(mat, heights = c(.7, 1, .5))

#Funtion to draw figure for Gaussian simulation to test the sensitivity to eHf variations
draw_pic <- function(permuteHf = 0, permuteUPb = 0, typeHf = "absolute", typeUPb = "relative", write = FALSE) {
  
  eHf = read.csv("EHf_resampled_mc_final_loc.csv")
  eHf = eHf[, 2:3]
  
  if (typeUPb == "relative") {
    eHf[,1] <- eHf[,1] + rnorm(n = length(eHf[,1]), mean = 0, sd = permuteUPb*eHf[,1])
  } else if (typeUPb == "absolute")  {
    eHf[,1] <- eHf[,1] + rnorm(n = length(eHf[,1]), mean = 0, sd = permuteUPb)
  }
  
  if (typeHf == "relative") {
    eHf[,2] <- eHf[,2] + rnorm(n = length(eHf[,2]), mean = 0, sd = permuteHf*eHf[,1])
  } else if (typeHf == "absolute") {
    eHf[,2] <- eHf[,2] + rnorm(n = length(eHf[,2]), mean = 0, sd = permuteHf)
  } 
  
  
  ind_in = order(eHf[, 1])
  eHf = eHf[ind_in, ]
  ##
  tr = c(0, 1000)
  ind = which(eHf[, 1] >= tr[1] & eHf[, 1] <= tr[2])
  eHf1000 = eHf[ind, ]
  #plot(eHf[,1],eHf[,2], xlim=c(0,1000))
  
  te = eHf1000[, 1]
  wd = 1
  
  t.incre = wd / 1
  start.time = seq(min(te), max(te) - t.incre, by = t.incre)
  ma = matrix(ncol = 6, nrow = length(start.time))
  
  for (i in 1:length(start.time)) {
    indtemp = which(te > start.time[i] & te < (start.time[i] + wd))
    ma[i, 1] = mean(te[indtemp])
    ma[i, 2:5] = summary(eHf1000[indtemp, 2])[c(3, 4, 2, 5)]
    ma[i, 6] = length(indtemp)
  }
  colnames(ma) <- c("age", "Median", "mean", "firstQ", "thirdQ", "sampleNo")
  
  ## draw
  layout(mat, heights = c(1, .5))
  op = par(mar = c(5, 4, 2, 2), bty = "n")
  par(op)
  plot(eHf1000[, 1], eHf1000[, 2], type = "n", xlim = xlim, ylim = ylim, ylab = "??Hf", xlab = "Age (Ma)", sub = paste("wd = ", wd, " Ma; step = ", t.incre, " Ma", sep = ""))
  #points(eHf[,1],eHf[,2], col = "black", pch=19, cex = .1)
  lines(ma[, 1], ma[, 2], col = "red", lwd = 2)
  lines(ma[, 1], ma[, 4], lty = 3, col = "green", lwd = 2)
  lines(ma[, 1], ma[, 5], lty = 3, col = "green", lwd = 2)
  plot(ma[, 1], ma[, 6], type = "l", ylab = "Number of samples", xlab = "Age (Ma)")
  
  if (write == TRUE) {
    write.csv(x = eHf1000, file = "Gaussian_noise_test/eHf1000_file_Hf_06_UPb_2_loc_mc.csv")
    write.csv(x = ma, file = "Gaussian_noise_test/ma_file_Hf_06_UPb_2_loc_mc.csv")
  }
}

set.seed(20161010)

pdf("Gaussian_noise_test/fig1_with_gaussian_devs_Hf_1_UPb_2_loc.pdf")
## permute eHf uncertainty range= median 0.62 (Q1=0.45 and Q3=1).
## permute U-Pb uncertainty range from median age discordance of 2% (Q1=0.7% and Q3=4.1%)
## for the Science advance paper we choose the median uncertainty values for our reference curve eHf?0.3; U-Pb age?1% 
draw_pic(permuteHf = 0.3, typeHf = "absolute", permuteUPb = 0.01, typeUPb = "relative", write = TRUE)
dev.off()





###############STEP 3###################### Normalized Sr/Sr in seawater curve
#Sensitivity analysis of normalized Sr_Sr in seawater to different parametization choices

dSr.raw = readxl::read_excel("New_dataset_6_2_2016/Strontium_input.xlsx")
Sr_Sr_seawater = data.frame(dSr.raw[which(!is.na(dSr.raw[, 3])), c(1, 3)])

#1st correction remove the radiogenic decay
#Several possible input for Rb_Sr_crust testing the sensitivity of this variable Rb/Sr can increase, decrease or be constant
#In the Science Advances paper Rb_Sr_crust=0.098
#Rb_Sr_crust=(0.0000000754*(1000-Sr_Sr_seawater[, 1])^2+7.54E-6*(1000-Sr_Sr_seawater[, 1])+0.098)/2
#Rb_Sr_crust=0.098*2-0.000098*Sr_Sr_seawater[, 1]
Rb_Sr_crust=0.098
Sr_Sr_decay<-0.70537+2.89*Rb_Sr_crust*(exp(1.42*10^-11*(1000-Sr_Sr_seawater[, 1])*10^6)-1)
dSr<-Sr_Sr_seawater[, 2] - Sr_Sr_decay

### linear interpolation to 1 ma interval
age_Sr_i = seq(min(Sr_Sr_seawater[, 1], na.rm = T), 998, by = 1)
dSr_i = (approx(Sr_Sr_seawater[, 1], dSr, xout = age_Sr_i, rule = 2)[[2]]) ## now the age_Sr_i and d_Sr_i are both sampled at 1 ma year interval

#2nd correction remove the effect of carbonates assuming carbonate represent 67% of the Sr flux to the ocean Mokadem et al. 2015
# this correction was not included in the Science Advances paper but should be in the future
#dSr_100<-rollapply(dSr_i, 100, mean, align = "left")
#age_Sr_100<-rollapply(age_Sr_i, 100, min, align = "left")

#age_Sr_1 = seq(min(age_Sr_100), max(age_Sr_100), by = 1)
#dSr_1 = approx(Sr_Sr_seawater[, 1], dSr, xout = age_Sr_1)[[2]] ## now the age_Sr_i and d_Sr_i are both sampled at 1 ma year interval

#dSr_1<-(dSr_1-0.65*dSr_100)/0.35
#dSr<-data.frame(age_Sr_1,dSr_1)
x11()
plot(age_Sr_i, dSr_i, xlim=c(0,1000), col="black", bg="blue", pch=21, xlab="Age (Ma)", ylab="87Sr/86Srseawater")
dSr<-data.frame(age_Sr_i,dSr_i)



################STEP4################
#Export data for time-series analysis
eHf = read.csv("Gaussian_noise_test/eHf1000_file_Hf_06_UPb_2_loc_mc.csv")
eHf = eHf[, c(2, 3)]
## plot
ind_in = order(eHf[, 1])
eHf = eHf[ind_in, ]
##
tr = c(0, 1000)
ind = which(eHf[, 1] >= tr[1] & eHf[, 1] <= tr[2])
eHf1000 = eHf[ind, ]
##
## plot
xlim = c(0, 1000)
x11()
op = par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))
plot(dSr, type = "l", xlim = xlim)
plot(eHf1000, type = "l", xlim = xlim)
par(op)
#plot(eHf[,1],eHf[,2], xlim=c(0,1000))

te = eHf1000[, 1]
wd = 1

t.incre = wd / 1
start.time = seq(min(te), max(te) - t.incre, by = t.incre)
ma = matrix(ncol = 6, nrow = length(start.time))

for (i in 1:length(start.time)) {
  
  # indtemp = start.ind[i]:(which(te - (te[start.ind[i]] + wd) >= 0)[1] - 1)
  indtemp = which(te > start.time[i] & te < (start.time[i] + wd))
  ma[i, 1] = mean(te[indtemp])
  ma[i, 2:5] = summary(eHf1000[indtemp, 2])[c(3, 4, 2, 5)]
  ma[i, 6] = length(indtemp)
}


colnames(ma) <- c("age", "Median", "mean", "firstQ", "thirdQ", "sampleNo")



## plot
xlim = c(0, 1000)
mat = matrix(1:3, ncol = 1)

pdf(file = "yx_file/data_import_for_comparison_v4.pdf", width = 8, height = 9)
layout(mat, heights = c(.7, 1, .5))

op = par(mar = c(5, 4, 2, 2), bty = "n")
plot(dSr[, 1], dSr[, 2], type = "l", col = "black", lwd = 1, ylab = " dSr", xlab = "Age (Ma)")
plot(eHf1000[, 1], eHf1000[, 2], type = "n", xlim = xlim, ylim = range(ma[, 4:5]), ylab = "??Hf", xlab = "Age (Ma)", sub = paste("wd = ", wd, " Ma; step = ", t.incre, " Ma", sep = ""))
lines(ma[, 1], ma[, 2], col = "black", lwd = 1)
lines(ma[, 1], ma[, 4], lty = 3, col = "grey")
lines(ma[, 1], ma[, 5], lty = 3, col = "grey")
plot(ma[, 1], ma[, 6], type = "h", ylab = "Number of samples", xlab = "Age (Ma)")
# lines(ma[, 1], - ma[, 3], col = "blue", lwd = 2)
par(op)

dev.off()

save(dSr, ma, file = "yx_file/finalV1/imported_data_for_compare_v4_1Ma.RData")


####STEP5###########
###Time-series analysis

###Function for the time-series analysis

##stdize function
stdize <- function(x) {
  stopifnot(is.vector(x, mode = "numeric"))
  if (any(is.na(x))) stop("NA value(s) exist in the input data")
  
  x_dm = dm(x)
  x_stdize = x_dm / sd(x, na.rm = FALSE)
  return(x_stdize)
}


###elementary functions
dm <- function(x) {
  stopifnot(is.vector(x, mode = "numeric"))
  if (any(is.na(x))) stop("NA value(s) exist in the input data")
  
  x_dm = x - mean(x, na.rm = FALSE)
  return(x_dm)
}

stdize <- function(x) {
  stopifnot(is.vector(x, mode = "numeric"))
  if (any(is.na(x))) stop("NA value(s) exist in the input data")
  
  x_dm = dm(x)
  x_stdize = x_dm / sd(x, na.rm = FALSE)
  return(x_stdize)
}

process_let <- function(type) {
  stopifnot(is.vector(type, mode = "character"))
  
  if (type == "dm") out_fun = dm
  if (type == "stdize") out_fun = stdize
  if (type == "taper") out_fun = spec.taper
  return(out_fun)
}


process <- function(x, type) {
  stopifnot(is.vector(x, mode = "numeric"))
  if (any(is.na(x))) stop("NA value(s) exist in the input data")
  
  x_temp = x
  
  N = length(type)
  for (i in 1:N) {
    x_temp = process_let(type[i])(x_temp)
  }
  
  x_out = x_temp
  
  return(x_out)
}

#bwfilter function
#' Filter the input signal using butterworth filter
#'
#' This function is a wrapper function for the butter and filtfilt functions from package signal that implement the butterworth filter more straightforward.
#' 
#' @param signal One-dimensional numeric vector representing the input signal to be filtered
#' @param cf Corner frequency(ies) for the filter
#' @param delta_t Sampling interval for the input signal
#' @param type The type of filter to be used with the option in "low" (lowpass), "high" (highpass), or "pass" (bandpass).
#' @param PLOT Whether to plot the filtered signal along with its original signal
#' @keywords butterworth filter, filter
#' @export
#' @return A vector of the input length.
#' @examples
#' y = rnorm(100)
#' cf_low = 0.3
#' bwfilter(y, cf_low, type = "low", PLOT = T)
#'
#' cf_band = c(0.1, 0.4)
#' bwfilter(y, cf_band, type = "pass", PLOT = T)

bwfilter <- function(signal, cf, delta_t = 1, type, PLOT=FALSE) {
  ### butterworth filter function
  
  stopifnot(
    is.vector(signal, mode = "numeric"),
    is.numeric(delta_t),
    is.numeric(cf),
    any(type == c("low", "high", "pass") & length(type) == 1)
  )
  
  library(signal)
  
  y = signal
  fNy = .5 / delta_t ## Nyquist frequency
  y_out = vector(mode = 'numeric', length = length(y))
  mean_y = mean(y)
  ym = y - mean_y
  
  if (type == 'low' | type == 'high') {
    stopifnot(length(cf) == 1)
    fil.but = butter(4, cf / fNy, type = type, plane = 'z')
    y_temp = filtfilt(fil.but, ym)
    y_out = y_temp + mean_y 
    if (type == 'low') filter.title = paste('lowpass filter:', cf)
    if (type == 'high') filter.title = paste('highpass filter:', cf)
  } else if (type == 'pass') {
    stopifnot(length(cf) == 2)
    fil.but1 = butter(4, W = cf[1] / fNy, type = "high", plane = "z")
    fil.but2 = butter(4, W = cf[2] / fNy, type = "low", plane = "z")
    y_temp1 = filtfilt(fil.but1, ym)
    y_temp = filtfilt(fil.but2, y_temp1)
    y_out = y_temp + mean_y
    filter.title = paste('bandpass filter:', cf[1], 'to', cf[2])
  }	
  
  if (PLOT) {
    dev.new()
    plot(y, type = 'l', col = 'grey', main = filter.title, xlab = 'Index')
    lines(y_out, col='red')
    legend('topright', legend = c('Original signal', 'filtered signal'), lty = 1, col = c('grey', 'red'))
  }
  
  invisible(y_out)
}








load("yx_file/finalV1/imported_data_for_compare_v4_1Ma.RData", verbose = T)
library(devtools)
library(signal)
load_all("~/Dropbox/apts/")

t1 = dSr[, 1]
t2 = ma[, 1]
#t2 = Sr_Sri[, 1]
#Convert eHf into Sr_Sri
Sr_Sri_median<-((-5.0741*ma[, 2]+49.433)/10000+1)*(0.699+0.0824*(exp(1.42*10^-11*(4.55*10^9-ma[, 1]*10^6))-1))
Sr_Sri_mean<-((-5.0741*ma[, 3]+49.433)/10000+1)*(0.699+0.0824*(exp(1.42*10^-11*(4.55*10^9-ma[, 1]*10^6))-1))
#Sr_Sri_median<-ma[, 2]
#Sr_Sri_mean<-ma[, 3]
Sr_Sri_median<-cbind(ma[, 1], Sr_Sri_median)
Sr_Sri_mean<-cbind(ma[, 1], Sr_Sri_mean)

y = list(dSr = dSr, Sr_Sri_median, Sr_Sri_mean)
#y = list(dSr = dSr, Sr_Sri_median = Sr_Sri[, 1:2], Sr_Sri_mean = Sr_Sri[, 1:2])

tr = c(min(t2), max(t1))
dt = .25  ## scale interested

ti = seq(tr[1], tr[2], by = dt)
yi = lapply(y, function(x, tout) approx(x[, 1], x[, 2], xout = tout)[[2]], tout = ti)

op = par(mfrow = c(2, 1), bty = "n", mar = c(4, 4, 2, 2))
plot(ti, yi[[1]], type = "l", ylab = "dSr")
plot(ti, yi[[2]], type = "l", ylab = "Sr_Sri")
lines(ti, yi[[3]], lty = 2)
par(op)

yi[[2]] = yi[[2]]
yi[[3]] = yi[[3]]

cf = c(1 / 700, 1 / 30)

yi_low = lapply(yi, function(x, cf, dt) bwfilter(signal = x, cf = cf[1], delta_t = dt, type = "low", PLOT = F), cf, dt)
yi_high = lapply(yi, function(x, cf, dt) bwfilter(signal = x, cf = cf[2], delta_t = dt, type = "high", PLOT = F), cf, dt)
yi_pass = lapply(yi, function(x, cf, dt) bwfilter(signal = x, cf = cf, delta_t = dt, type = "pass", PLOT = F), cf, dt)

pdf(file = "yx_file/finalV1/compare_v7.pdf", width = 14, height = 7, colormodel = "cmyk")
op = par(mfcol = c(3, 4), bty = "n", mar = c(4, 4, 2, 2))
plot(ti, yi[[1]], type = "l", ylab = "dSr", xlab = "Age /Ma", col = "blue")
plot(ti, yi[[2]], type = "l", ylab = "Sr_Sri", xlab = "Age /Ma", col = "red")
legend("topleft", legend = c("Median", "mean"), lty = 1, col = c("red", "green"))
lines(ti, yi[[3]], col = "green")
yi_ccf_median = ccf(x = yi[[1]], y = yi[[2]], lag.max = 1000, plot = F)
yi_ccf_mean = ccf(x = yi[[1]], y = yi[[3]], lag.max = 1000, plot = F)
plot(yi_ccf_median[[4]] * dt, yi_ccf_median[[1]], type = "l", col = "red", xlab = "Lag /Ma (positive lag = dSr leads)", ylab = "Correlation coefficient", xlim=c(-200, 200), ylim=c(-0.2,0.6))
lines(yi_ccf_mean[[4]] * dt, yi_ccf_mean[[1]], col = "green")
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)



plot(ti, yi_low[[1]], type = "l", ylab = "dSr (low freq)", xlab = "Age /Ma", col = "blue")
plot(ti, yi_low[[2]], type = "l", ylab = "Sr_Sri (low freq)", xlab = "Age /Ma", col = "red", ylim = range(c(yi_low[[2]], yi_low[[3]])))
lines(ti, yi_low[[3]], col = "green")
yi_low_ccf_median = ccf(x = yi_low[[1]], y = yi_low[[2]], lag.max = 1000, plot = F)
yi_low_ccf_mean = ccf(x = yi_low[[1]], y = yi_low[[3]], lag.max = 1000, plot = F)
plot(yi_low_ccf_median[[4]] * dt, yi_low_ccf_median[[1]], type = "l", col = "red", xlab = "Lag /Ma (positive lag = dSr leads)", ylab = "Correlation coefficient", xlim=c(-200, 200), ylim=c(-0.2,1))
lines(yi_low_ccf_mean[[4]] * dt, yi_low_ccf_mean[[1]], col = "green")
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

plot(ti, yi_pass[[1]], type = "l", ylab = "dSr (mid freq)", xlab = "Age /Ma", col = "blue")
plot(ti, yi_pass[[2]], type = "l", ylab = "Sr_Sri (mid freq)", xlab = "Age /Ma", col = "red", ylim = range(c(yi_pass[[2]], yi_pass[[3]])))
lines(ti, yi_pass[[3]], col = "green")
yi_pass_ccf_median = ccf(x = yi_pass[[1]], y = yi_pass[[2]], lag.max = 1000, plot = F)
yi_pass_ccf_mean = ccf(x = yi_pass[[1]], y = yi_pass[[3]], lag.max = 1000, plot = F)
plot(yi_pass_ccf_median[[4]] * dt, yi_pass_ccf_median[[1]], type = "l", col = "red", xlab = "Lag /Ma (positive lag = dSr leads)", ylab = "Correlation coefficient", xlim=c(-200, 200), ylim=c(-0.2,0.6))
lines(yi_pass_ccf_mean[[4]] * dt, yi_pass_ccf_mean[[1]], col = "green")
abline(v = 0, lty = 2)
abline(h = 0, lty = 2)

plot(ti, yi_high[[1]], type = "l", ylab = "dSr (high freq)", xlab = "Age /Ma", col = "blue")
plot(ti, yi_high[[2]], type = "l", ylab = "Sr_Sri (high freq)", xlab = "Age /Ma", col = "red")
lines(ti, yi_high[[3]], col = "green")
yi_high_hist1 = hist(stdize(yi_high[[1]]), plot = F)
yi_high_hist2 = hist(stdize(yi_high[[2]]), plot = F)
yi_high_hist3 = hist(stdize(yi_high[[3]]), plot = F)
plot(yi_high_hist1[[4]], yi_high_hist1[[3]], type = "n", ylab = "Probability", xlab = "Normalized yi_high", ylim = range(c(yi_high_hist1[[3]], yi_high_hist2[[3]], yi_high_hist3[[3]])))
x.norm = seq(-5, 5, length = 100)
y.norm = dnorm(x.norm)
lines(x.norm, y.norm, lty = 1, col = "grey")
lines(yi_high_hist2[[4]], yi_high_hist2[[3]], col = "red", lty = 1)
lines(yi_high_hist1[[4]], yi_high_hist1[[3]], col = "blue", lty = 1)
lines(yi_high_hist3[[4]], yi_high_hist3[[3]], col = "green", lty = 1)

par(op)
dev.off()

wd = 150
incre = 10
t.start = seq(tr[1], tr[2] -wd, by = incre)
tm = cbind(t.start, t.start + wd)
tm = cbind(tm, apply(tm, 1, mean))

yi_comb<-lapply(seq_along(yi_pass),function(i) (unlist(yi_pass[i])+unlist(yi_low[i])))
dsr = yi_comb[[1]]-mean(yi_comb[[1]]/2)
ehf_mean = yi_comb[[3]]-mean(yi_comb[[3]]/2)
ehf_median=yi_comb[[2]]-mean(yi_comb[[2]]/2)
ehf=list(ehf_median,ehf_mean)

ccfm = list()

N = nrow(tm)
for (i in 1:N) {
  tind = which(ti >= tm[i, 1] & ti <= tm[i, 2])
  ccfm[[i]] = lapply(ehf, function(x, tind, y, wd) ccf(x = x[tind], y = y, lag.max = floor(wd / 8 / dt), plot = F), tind = tind, y = dsr[tind], wd = wd)
}

ind.lagneg = which(ccfm[[1]][[1]][[4]] <= 10)

maxpos_median = lapply(ccfm, function(x) {
  temp = which.max(x[[1]][[1]][ind.lagneg])
  cbind(x[[1]][[4]][temp], x[[1]][[1]][temp])
})

maxpos_mean = lapply(ccfm, function(x) {
  temp = which.max(x[[2]][[1]][ind.lagneg])
  cbind(x[[1]][[4]][temp], x[[1]][[1]][temp])
})

maxpos_median_m = matrix(unlist(maxpos_median), ncol = 2, byrow = T)
maxpos_mean_m = matrix(unlist(maxpos_mean), ncol = 2, byrow = T)

par(mfrow=c(1,1))
step<-seq(0,length(dsr)-1,1)/4
plot(step,dsr, typ="l", col='blue')
plot(step,ehf[[2]], typ="l", col='green')
pdf(file = "yx_file/finalV1/compare_v6_lagmw.pdf", width = 12, height = 7, colormodel = "cmyk")
op = par(mfrow = c(2, 1), bty = "n", mar = c(4, 4, 2, 2))
plot(tm[, 1], maxpos_median_m[, 1] * dt, type = "o"
     , ylab = "Lag at maximum pos correlation /Ma", xlab = "Age /Ma"
     , col = "red", ylim = wd / 2 * c(-1, 0))
lines(tm[, 1], maxpos_mean_m[, 1] * dt, type = "o", col = "green")


plot(tm[, 1], maxpos_median_m[, 2], col = "red", type = "o"
     , ylab = "Maximum postive correlation coefficients", xlab = "Age /Ma"
     , ylim = range(c(maxpos_median_m[, 2], maxpos_mean_m[, 2])))
abline(h = .5, col = "grey")
lines(tm[, 1], maxpos_mean_m[, 2], col = "green", type = "o")
par(op)
dev.off()

# maxneg_median = lapply(ccfm, function(x) {
# 	temp = which.min(x[[1]][[1]])
# 	cbind(x[[1]][[4]][temp], x[[1]][[1]][temp])
# })
# 
# maxneg_mean = lapply(ccfm, function(x) {
# 	temp = which.min(x[[2]][[1]])
# 	cbind(x[[1]][[4]][temp], x[[1]][[1]][temp])
# })