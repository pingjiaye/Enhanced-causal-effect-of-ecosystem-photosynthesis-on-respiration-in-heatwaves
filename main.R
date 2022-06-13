rm(list = ls())
library(dplyr)
library(plyr)
library(rEDM)
library(Kendall)

# Load data (heatwave 2018)
file <- Sys.glob("~/NutstoreCloudBridge/Data/Europe/FLX*", dirmark = FALSE)
EuroDat <- NULL
for (i in 1:length(file)) {
  flux.list <- lapply(file[i], read.csv)
  flux.dataframe <- ldply(flux.list,data.frame) %>%
    select(TIMESTAMP,GPP_NT_VUT_MEAN,RECO_NT_VUT_MEAN) %>%
    mutate(SITE = rep(substr(file[i],32,37)))
  YEAR <- substr(as.character(flux.dataframe$TIMESTAMP),1,4)
  MONTH <- substr(as.character(flux.dataframe$TIMESTAMP),5,6)
  DAY <- substr(as.character(flux.dataframe$TIMESTAMP),7,8)
  flux.dataframe <- flux.dataframe %>%
    select(-TIMESTAMP) %>%
    mutate(YEAR=YEAR,MONTH=MONTH,DAY=DAY,NO=rep(i)) %>%
    filter(YEAR == 2017 | YEAR == 2018)
  EuroDat <- rbind(EuroDat, flux.dataframe)
}
DatCCM <- EuroDat %>%
  filter(as.numeric(EuroDat$MONTH)==06|as.numeric(EuroDat$MONTH)==07|as.numeric(EuroDat$MONTH)==08) %>%
  select(GPP_NT_VUT_MEAN,RECO_NT_VUT_MEAN,SITE,YEAR,NO)
colnames(DatCCM)=c("GPP","ER","Site","Year","SiteID")

# Function for generating lag time series
laf=function(x,y,lagf){
  n <- NROW(x)
  x.t=x;y.t=y
  if(lagf<=0){x.t=x.t[(1-lagf):n];y.t=y.t[1:(n+lagf)]} # if lagf<0, y is leading
  if(lagf>0){x.t=x.t[1:(n-lagf)];y.t=y.t[(1+lagf):n]}  # if lagf>0, x is leading           
  return(cbind(x.t,y.t))
}

# Calculate the CCM for each site and combined the results
GEccm.output <- NULL
star_time <- Sys.time()
for (k in unique(DatCCM$SiteID)) {
  DHdat <- DatCCM %>%
    filter(SiteID == k)
  SiteName <- as.character(unique(DHdat$Site))
  DHdat <- DHdat %>%
    select(Year,GPP,ER)
  n = nrow(DHdat)
  seed = 209059
  
  # Set causal links
  indcauz <- matrix(0,2,2)
  indcauz[,1] <- c("ER","GPP")
  indcauz[,2] <- c("GPP","ER")
  colnames(indcauz) = c("Effect","Cause")
  
  # Determine the embedding dimension E
  Emax = 10
  Ed = NULL
  for (i in 1:nrow(indcauz)) {
    E.test = NULL
    for (E.t in 2:Emax) {
      ccmE.t <- ccm(DHdat, E = E.t, lib_column = indcauz[i,1], target_column = indcauz[i,2],
                    lib_sizes = n, tp=-1,random_libs = F)
      E.test = c(E.test, mean(ccmE.t$rho))
    }
    Ed = c(Ed, which.max(E.test)+1)
  }
  
  # Calculate CCM
  LibSize = sort(c(5,10,20,30,seq(40,n,50),n))
  ccmda = NULL
  for (i in 1:nrow(indcauz)) {
    ccmda.t = NULL
    for (j in 0:-3) {
      ccm.t=laf(DHdat[,indcauz[i,1]],DHdat[,indcauz[i,2]],lagf=j)
      colnames(ccm.t) = indcauz[i,]
      Efect_Xmap_Cauz <- ccm(ccm.t, E = Ed[i], lib_column = indcauz[i,"Effect"],
                             target_column = indcauz[i,"Cause"], lib_sizes = LibSize,
                             tp = 0, RNGseed = seed, num_samples = 100,replace = F)
      # Average each libraty size's cross-map skill
      aveg=cbind(unique(Efect_Xmap_Cauz$lib_size),
                 aggregate(Efect_Xmap_Cauz[,c('rho')], by=list(as.factor(Efect_Xmap_Cauz$lib_size)), mean)[,'x'],
                 aggregate(Efect_Xmap_Cauz[,c('mae')], by=list(as.factor(Efect_Xmap_Cauz$lib_size)), mean)[,'x'],
                 aggregate(Efect_Xmap_Cauz[,c('rmse')],by=list(as.factor(Efect_Xmap_Cauz$lib_size)), mean)[,'x'])
      ccm_mean = data.frame(lag = rep(j, nrow(aveg)), Efect_Xmap_Cauz[1:nrow(aveg),])
      ccm_mean[,c('lib_size','rho','mae','rmse')]=aveg
      ccm_mean[ccm_mean[,'rho']<0, 'rho'] = 0
      # CCM Convergence test
      # Kendall's tau test
      if(length(ccm_mean$rho)>3){
        kend=MannKendall(ccm_mean$rho)
        kend.tau=kend$tau[1]
        kend.p=kend$sl[[1]]
      }else{
        kend.tau=NA
        kend.p=NA
      }
      # Fisher's delta rho Z test
      rho.Lmax=ccm_mean$rho[which.max(ccm_mean$lib_size)]
      rho.Lmin=ccm_mean$rho[1]
      ns=min(sum(!is.na(DHdat[,indcauz[i,1]])),sum(!is.na(DHdat[,indcauz[i,2]])))
      delta_rho=rho.Lmax-rho.Lmin
      z=abs(0.5*(log((1+rho.Lmax)/(1-rho.Lmax))-log((1+rho.Lmin)/(1-rho.Lmin)))*(2/(ns-3))^-0.5)
      z.p=(1-pnorm(z))
      # Compile all the testing results
      ccmda.t = rbind(ccmda.t,
                      unlist(c(ccm_mean[1,1:5],rho_Lmax=rho.Lmax,rho_Lmin=rho.Lmin,
                               Ztest=z,p_Ztest=z.p,Kendall_tau=kend.tau,Kendall_p=kend.p)))
    }
    # Select the CCM results based on the maximum cros-map skill
    ccmda = rbind(ccmda,ccmda.t[which.max(ccmda.t[,'rho_Lmax']),])
  }
  
  ccmda <- data.frame(indcauz,ccmda)
  convergence = ccmda$p_Ztest < 0.05 & ccmda$Kendall_p <= 0.05 & ccmda$Kendall_tau>0
  ccmda <- data.frame(ccmda, convergence)
  
  istd=data.frame(system=rep(SiteName),year=unique(DHdat$Year), ccmda[,c('Cause','Effect','rho_Lmax','p_Ztest','Kendall_p','convergence')])
  
  # Save output
  outputdat <-left_join(istd,ccmda)
  GEccm.output <- rbind(GEccm.output,outputdat)
}
end_time <- Sys.time()
runT <- end_time-star_time

write.csv(GEccm.output, "~/NutstoreCloudBridge/Data/Europe/CCM_DH18.csv")
