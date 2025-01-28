################################################################################
## Code for simulating survey data for monitoring ruffed grouse (Bonasa umbellus)
##  across a 858,000 ha study area in North Georgia (GA), USA, from Clayton D. Delancey,
##  William B. Lewis, Gregory T. Wann, Richard B. Chandler, Emily Rushton, and James A.
##  Martin. Assessing the utility of autonomous recording units and spring drumming surveys
##  for monitoring abundance of ruffed grouse (Bonasa umbellus).
## Corresponding author: Will Lewis (wblewis7@gmail.com), University of Georgia.
################################################################################



# Simulating monitoring data for ruffed grouse (RUGR) via spring roadside drumming 
#   point counts with or without deployment of Autonomous Recording Units (ARUs).
#   Grouse abundance in this area has been monitored using spring roadside drumming 
#   counts at 59 roadside survey routes, each consisting of 8 - 15 survey points
#   (633 survey point total).
# Local abundance is simulated as a function of 5 landscape covariates, while the
#   detection process is simulated as a function of availability and conditional
#   detectability. Point counts and ARUs use a similar state process but have
#   different observation processes, though both jointly estimate the seasonal 
#   change in availability with drumming availability peaking on April 19th.
# Simulating data for monitoring strategies varying in the average density across
#   the study area (0.003, 0.006, 0.011, 0.022), the number of survey routes (10, 20,
#   30, 40, 50), the number of replicate point count surveys (1 - 3), and the number
#   of survey points at which ARUs are deployed (0, 40, 80). ARU deployment is 
#   simulated in two ways. First, ARUs are randomly deployed at 40/80 point count 
#   locations to assess occupancy and availability (Method OcAv). Estimating availability
#   is conditional on detections, so the second method (Method Av) preferentially deploys
#   ARUs at survey locations with the highest likelihood of detections (i.e., highest 
#   expected abundance). Since ARUs are not randomly deployed in this method, the data
#   shouldn't be used to assess occupancy and landscape covariate effects.
# Simulation parameters are taken from a previous analysis of point count data (Lewis et al. 
#   2022. Abundance and distribution of ruffed grouse Bonasa umbellus at the southern 
#   periphery of the range. Wildlife Biology) and a preliminary analysis
#   of ARU data from the study region.





# Simulation Data --------------------------------------------------------------
################################################################################

load("RUGR.NGA.simulation.data.gzip")

# surveydata is a matrix giving the survey-level landscape covariates for each survey location
#   in northern Georgia. Suvey points (Point) are organized into 59 roadside survey
#   routes (Route). Survey-level landscape covariates are mean elevation (elev), Latitude,
#   aspect (decomposed to eastness and northness), and canopy height (canheight). All 
#   covariates are standardized.
# lambda_params_intercepts is a vector giving the intercept for the abundance modeling. Data
#   were simulated under four different levels of abundance, corresponding to average 
#   grouse densities of 0.003, 0.006, 0.011, and 0.022 males/ha, respectively.
# lambda_params_cov is a vector giving the effects of elevation, Latitude, eastness, northness
#   and canopy height (respectively) on grouse abundance.
# phi_params_PCintercept gives the intercept for the availability portion of the point count
#   observation process.
# phi_params_ARUintercept gives the intercept for the availability portion of the ARU
#   observation process. To compare with phi_params_PCintercept and allow for similar modeling
#   of the effects of day of year on availability (phi_params_doy), this parameter represents
#   the availability intercept for a 5-minute segment of an ARU recording.
# phi_params_doy is a vector giving the linear and quadratic effect of day of year, respectively,
#   on the availability portion of the observation process (shared between point counts and ARUs).
# p_params_PC gives the sigma parameter informing the scale of the conditional detectability
#   portion of the point count observation process (modeled with a Half-Normal distribution).
# pdet_ARU gives the probability of detecting at least one grouse drum in a 3-hour ARU recording
#   for informing the conditional detectability portion of the ARU observation process.
# ARU_dates is a vector giving the days of year (standardized) for which ARUs were simulated as being
#   deployed over. This range corresponds to day of year 60 (March 1) to 151 (May 31st).
# PC_period_dates is a 3-element list giving the range of simulated survey dates for point counts
#   repeated either 1 (first element), 2 (second element), or 3 times (3rd element). Within each element,
#   rows represent sampling periods. Each point count sampling period lasts 2 weeks, with periods centered
#   on the date of peak drumming activity in the region (April 19th). Surveys occurred from
#   from April 12th – April 25th for scenarios with 1 survey, from April 5th – April 18th (1st survey)
#   and April 19th – May 2nd (2nd survey) for scenarios with 2 repeat surveys, and from March 28th –
#   April 11th (1st survey), April 12th – April 25th (2nd survey), and April 26th – May 9th (3rd survey)
#   for scenarios with 3 repeat surveys. Day of year values are standardized.
# surveyarea gives the area in hectares surveyed on each point count (200m fixed-radius survey).
# PC_det_bands gives the boundaries of distance bands (out to 200m) in which detected grouse are binned
#   into during point count surveys.


# Calculating relative area and detection probability (based on Half-Normal detection function)
#   for each 20-m distance band for point counts
Nbands <- length(RUGR.NGA.simulation.data$PC_det_bands) - 1 # Number of bands
band_area <- rep(NA,times=Nbands) # Area of bands
for (i in 1:Nbands){
  band_area[i] <- (RUGR.NGA.simulation.data$PC_det_bands[i+1]^2*pi) - (RUGR.NGA.simulation.data$PC_det_bands[i]^2*pi)
}
rel_band_area <- band_area/sum(band_area) # Relative area of bands
band_det_prob <- rep(NA,times=Nbands) # Detection probability in each band
for (i in 1:Nbands){
  band_det_prob[i] <- (RUGR.NGA.simulation.data$p_params_PC^2*(1-exp(-RUGR.NGA.simulation.data$PC_det_bands[i+1]^2/(2*RUGR.NGA.simulation.data$p_params_PC^2)))-RUGR.NGA.simulation.data$p_params_PC^2*(1-exp(-RUGR.NGA.simulation.data$PC_det_bands[i]^2/(2*RUGR.NGA.simulation.data$p_params_PC^2))))*2*3.1416/(RUGR.NGA.simulation.data$surveyarea*10000*rel_band_area[i])
}




# Simulation Parameters --------------------------------------------------------
################################################################################

Density <- RUGR.NGA.simulation.data$lambda_params_intercepts # Intercepts for abundance modeling varying background grouse abundance
Density.val <- c(0.003,0.006,0.011,0.022) # Average density (males/ha) across the North GA study area corresponding to each intercept in Density
Abundance.val <- c(2675, 5116, 9785, 18714) # Total population abundance across the North GA study area corresponding to each intercept in Density
Periods <- 1:3 # Number of repeat point count surveys
Routes <- seq(10,50,by=10) # Number of survey routes
ARUs <- c(0,40,80) # Number of survey sites with ARU deployment
sims <- 500 # Number of simulation per scenario







# Looping over parameters and simulating datasets ------------------------------
################################################################################


for(z in 1:length(Density)){
  
  # Calculating expected abundance at each survey location. Accounting for offset
  #   of survey area to account for different area between surveys (12.6 ha) and size
  #   of grid cells in North Georgia study area (0.81 ha) for prediction.
  surveys <- RUGR.NGA.simulation.data$surveydata
  surveys$EN <- exp(as.matrix(cbind(1,surveys[,3:7],log(RUGR.NGA.simulation.data$surveyarea))) %*% c(RUGR.NGA.simulation.data$lambda_params_intercepts[z],RUGR.NGA.simulation.data$lambda_params_cov,1))
  
  for(r in 1:length(Routes)){
    
    for(p in 1:length(Periods)){
      
      sim_covs <- sim_pc_y <- sim_pc_ydb <- sim_pc_doy <- sim_ARU_Av <- sim_ARU_OcAv <- vector(mode="list", length=sims)
      
      for(x in 1:sims){
        
        surveydata <- surveys
        
        # Simulating abundance at each point based on expected abundance
        surveydata$N <- 0
        for (i in 1:nrow(surveydata)){
          surveydata$N[i] <- rpois(1,surveydata$EN[i])
        } #i
        
        # Simulating which distance band each bird is in. Assuming not moving 
        #   between repeat surveys
        point_dist <- matrix(NA, ncol=Nbands, nrow=nrow(surveydata))
        for (i in 1:nrow(point_dist)){
          point_dist[i,] <- rmultinom(1,size=surveydata$N[i],prob=rel_band_area)
        } #i
        point_dist <- cbind(surveydata[1:2],point_dist)
        
        # Selecting r routes
        goodroutes <- sample(unique(surveydata$Route), Routes[r], replace=F)
        surveydata_good <- surveydata[surveydata$Route %in% goodroutes,]
        point_dist_good <- point_dist[point_dist$Route %in% goodroutes,]
        
        
        ## Simulating point count data -----------------------------------------
        
        # Simulating day of year and availability for each survey at each point during 1-3 periods
        # Assuming all surveys from same route done on same day
        PC.avail <- doy_sim <- matrix(NA, ncol=Periods[p], nrow=length(goodroutes))
        for(j in 1:Periods[p]){
          doy_sim[,j] <- sample(RUGR.NGA.simulation.data$PC_period_dates[[Periods[p]]][j,],nrow(PC.avail),replace=T)
          doy_sim_mat <- cbind(1,doy_sim[,j],doy_sim[,j]^2)
          PC.avail[,j] <- plogis(doy_sim_mat %*% c(RUGR.NGA.simulation.data$phi_params_PCintercept,RUGR.NGA.simulation.data$phi_params_doy))
        } #j
        PC.avail <- merge(surveydata_good[,1:2],cbind(data.frame(Route=goodroutes),PC.avail),all.x=T,by="Route")
        doy_sim <- merge(surveydata_good[,1:2],cbind(data.frame(Route=goodroutes),doy_sim),all.x=T,by="Route")
        
        # Simulating observation process for point counts
        ydb_RUGR <- array(NA, dim=c(nrow=nrow(point_dist_good), ncol=ncol(point_dist_good)-2,Periods[p]))
        for (i in 1:Periods[p]){
          for (j in 1:nrow(point_dist_good)){
            for (k in 1:(ncol(point_dist_good)-2)){
              RUGR_per_avail <- rbinom(1,point_dist_good[j,k+2],PC.avail[j,i+2]) # Number available
              ydb_RUGR[j,k,i] <- rbinom(1,RUGR_per_avail,band_det_prob[k]) # Number detected
            } #j
          } #k
        } #i
        # Total detections per survey
        y_RUGR <- apply(ydb_RUGR,3,rowSums)
        
        
        ## Simulating ARU data -------------------------------------------------
        
        # Can handle ARUs in two ways. First, simulating ARUs at randomly-selected
        #   locations to inform occupancy and effect of day of year on availability
        #   (Method OcAv). Estimating availability requires sites to be occupied,
        #   so could also preferentially deploy at areas with high expected abundance
        #   to assess seasonal availability (Method Av).
        
        # Availability on ARUs
        ARU.avail.5min <- plogis(cbind(1,RUGR.NGA.simulation.data$ARU_dates,RUGR.NGA.simulation.data$ARU_dates^2) %*% c(RUGR.NGA.simulation.data$phi_params_ARUintercept,RUGR.NGA.simulation.data$phi_params_doy))
        # This is availability for 5 minutes, changing to 3-hour recording
        ARU.avail <- 1-(1-ARU.avail.5min)^36
        
        # Simulating data for both methods
        ARU.det.OcAv <- ARU.det.Av <- vector(mode="list",length=length(ARUs))
        
        for(a in 1:length(ARUs)){
          
          if(ARUs[a]==0){
            
            # If no ARUs, setting as NA
            ARU.det.Av[[1]] <- matrix(NA, nrow=1, ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            ARU.det.OcAv[[1]] <- matrix(NA, nrow=nrow(surveydata_good), ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            
          } else{
            
            ### Method Av ------------------------------------------------------
            
            ### Deploying at areas of highest expected abundance, not necessarily
            ###   same points doing point counts at.
            Occupancy.Av <- ifelse(surveydata$N>0, 1, 0)
            ARU.sites.numbers.Av <- sample(1:nrow(surveydata),ARUs[a],prob=surveydata$EN,replace=F)
            ARU.Av <- matrix(0, nrow=length(ARU.sites.numbers.Av),ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            for(i in 1:nrow(ARU.Av)){
              for(j in 1:ncol(ARU.Av)){
                ARU.Av[i,j] <- rbinom(1,1,Occupancy.Av[ARU.sites.numbers.Av[i]] * ARU.avail[j,1] * RUGR.NGA.simulation.data$pdet_ARU)
              } #j
            } #i
            # Removing ARUs without detections
            ARU.det.Av[[a]] <- ARU.Av[rowSums(ARU.Av)>0,]
            
            ### Method OcAv ----------------------------------------------------
            
            ### Deploying randomly at sites at which point counts are performed.
            Occupancy.OcAv <- ifelse(surveydata_good$N>0, 1, 0)
            ARU.sites.numbers.OcAv <- sample(1:nrow(surveydata_good),ARUs[a],replace=F)
            ARU.OcAv <- matrix(NA, nrow=nrow(surveydata_good),ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            for(i in 1:length(ARU.sites.numbers.OcAv)){
              for(j in 1:ncol(ARU.OcAv)){
                ARU.OcAv[ARU.sites.numbers.OcAv[i],j] <- rbinom(1,1,Occupancy.OcAv[ARU.sites.numbers.OcAv[i]] * ARU.avail[j,1] * RUGR.NGA.simulation.data$pdet_ARU)
              } #i
            } #j
            ARU.det.OcAv[[a]] <- ARU.OcAv
          } #else
        } #a
        
        sim_covs[[x]] <- surveydata_good
        sim_pc_y[[x]] <- y_RUGR
        sim_pc_ydb[[x]] <- ydb_RUGR
        sim_pc_doy[[x]] <- doy_sim
        sim_ARU_Av[[x]] <- ARU.det.Av
        sim_ARU_OcAv[[x]] <- ARU.det.OcAv
        
      } #x
      
      # Writing out
      sim_meta <- list(AvgDensity = Density.val[z],
                       TotAbund = Abundance.val[z],
                       Routes = Routes[r],
                       Periods = Periods[p],
                       ARUs = ARUs,
                       Nsims = sims,
                       PC_detbands = list(Nbands = Nbands,
                                          det_bands = RUGR.NGA.simulation.data$PC_det_bands,
                                          rel_band_area = rel_band_area,
                                          survey.area = RUGR.NGA.simulation.data$surveyarea),
                       lambda_params = c(RUGR.NGA.simulation.data$lambda_params_intercepts[z],RUGR.NGA.simulation.data$lambda_params_cov),
                       phi.ARUintercept.params = RUGR.NGA.simulation.data$phi_params_ARUintercept,
                       phi.PCintercept.params = RUGR.NGA.simulation.data$phi_params_PCintercept,
                       phi.doy.params = RUGR.NGA.simulation.data$phi_params_doy,
                       p.PC.sigma = RUGR.NGA.simulation.data$p_params_PC,
                       p.ARU = RUGR.NGA.simulation.data$pdet_ARU,
                       ARU.dates = RUGR.NGA.simulation.data$ARU_dates)
      sim_data <- list(sim_covs = sim_covs,
                       sim_pc_y = sim_pc_y,
                       sim_pc_ydb = sim_pc_ydb,
                       sim_pc_doy = sim_pc_doy,
                       sim_ARU_Av = sim_ARU_Av,
                       sim_ARU_OcAv = sim_ARU_OcAv)
      
      sim_out <- list(sim_meta = sim_meta,
                      sim_data = sim_data)
      
      filename <- paste("SimulationDataRUGR_",Routes[r],"Routes_",Periods[p],"Periods_",Density.val[z],"avgdensity.gzip",sep="")
      save(sim_out,file=filename)
    } #p
  } #r
} #z
