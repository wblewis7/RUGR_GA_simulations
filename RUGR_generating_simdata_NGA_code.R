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
#   point counts along 59 roadside survey routes, each consisting of 8 - 15 survey points
#   (633 survey point total).
# Local abundance is simulated as a function of 5 landscape covariates while the
#   detection process is simulated as a function of availability and conditional
#   detectability. Point counts and ARUs use a similar state process but have
#   different observation processes, though both jointly estimate the seasonal 
#   change in availability with drumming availability peaking on April 19th.
# Simulating data for monitoring strategies varying in the average density across
#   the study area (0.003, 0.006, 0.011, 0.022), the number of survey routes (20,
#   30, 40, 50), the number of replicate point count surveys (1 - 3), and the number
#   of survey points at which ARUs are deployed (0, 40, 80). The peak of drumming
#   activity can vary regionally and is often not known with certainty, so we model
#   surveys as being performed over a 12-week period (March 9 - May 31).
# We use the hierarchical distance-sampling models of Royle et al. (2004), modified
#   to account for repeat surveys and to separate availability from detectability,
#   to model abundance from point counts.
# ARUs are simulated as being randomly deployed at 0/40/80 point count locations. The
#   identity of grouse detected on ARU recordings is not known with certainty, so we
#   treat ARU data as a measure of occupancy (i.e., a 1 if at least one drum is detected).
#   We use the methods of Royle & Nichols (2003) for relating occupancy data to abundance,
#   using each ARU recording as a replicate survey. We decompose the per-individual
#   detection rate into availability (based on day of year) and the probability of detecting
#   at least one drum per individual (which should be near 1 based on preliminary analysis).
# Simulation parameters are taken from a previous analysis of point count data (Lewis et al. 
#   2022) and a preliminary analysis of ARU data from the study region.





# Simulation Data --------------------------------------------------------------
################################################################################

load("RUGR.NGA.simulation.data.gzip")

# surveydata is a matrix giving the survey-level landscape covariates for each survey location
#   in northern Georgia. Survey points (Point) are organized into 59 roadside survey
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
# pdet_ARU gives the probability of detecting at least one drum from an available grouse during a 3-hour ARU recording
#   for informing the conditional detectability portion of the ARU observation process.
# ARU_dates is a vector giving the days of year (standardized) for which ARUs were simulated as being
#   deployed over. This range corresponds to day of year 68 (March 9) to 151 (May 31st).
# PC_dates is a 12 x 7 matrix giving the day of year (standardized) for each day of the 12-week season.
#   Rows correspond to weeks of the season (1-12) while columns correspond to days of the week (1-7). 
#   This matrix is used to simulate the day of year for point counts.
# PC_periods is a 12 x 3 matrix giving the weeks (rows) in which to perform surveys when using either 1
#   (1st column), 2 (2nd column), or 3 (3rd column) repeat point counts. With only 1 repeat point count 
#   survey, survey dates were roughly evenly divided across the entire 12-week period to ensure adequate
#   temporal coverage. With 2 repeat point count surveys, survey dates for the first survey were roughly
#   evenly divided across the first 6 weeks and the dates for the second survey were roughly evenly
#   divided across the second 6 weeks. For 3 repeat point count surveys, survey dates for the first survey
#   were roughly evenly divided across the first 4 weeks, the dates for the second survey were roughly 
#   evenly divided across weeks 5 - 8, and the dates for the third survey were roughly evenly divided
#   across weeks 9 - 12.
# surveyarea gives the area in hectares surveyed on each point count (200m fixed-radius survey).
# PC_det_bands gives the boundaries of distance bands (out to 200m) in which detected grouse are binned
#   into during point count surveys.
# N_GA_covariates_map gives the standardized landscape covariates (elev, Latitude, eastness, northness, 
#   and canheight) for each of the 57070 390-390m grid cells across the North Georgia study area. This 
#   matrix is used to predict RUGR abundance across the study area via simulation posterior samples in
#   RUGR_analyzing_simdata_NGA_code.R. The area of each grid cell (ha) is included as an intercept to 
#   adjust between the area of surveys (12.6 ha) and grid cells (15.21 ha).
# abund_N_GA_actual gives the actual grouse abundance across the study area under the four different
#   avrage RUGR densities (0.003, 0.006, 0.011, 0.022 males/ha).



# Calculating the relative area and detection probability (based on Half-Normal detection function)
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
Periods <- 1:3 # Number of repeat point count surveys
Routes <- seq(20,50,by=10) # Number of survey routes
ARUs <- c(0,40,80) # Number of survey sites with ARU deployment
sims <- 200 # Number of simulations per scenario







# Looping over parameters and simulating datasets ------------------------------
################################################################################


for(z in 1:length(Density)){
  
  # Calculating expected abundance at each survey location. Including offset
  #   of survey area to account for different area between surveys (12.6 ha) and size
  #   of grid cells in North Georgia study area (15.21 ha) for prediction.
  surveys <- RUGR.NGA.simulation.data$surveydata
  surveys$EN <- exp(as.matrix(cbind(1,surveys[,3:7],log(RUGR.NGA.simulation.data$surveyarea))) %*% c(RUGR.NGA.simulation.data$lambda_params_intercepts[z],RUGR.NGA.simulation.data$lambda_params_cov,1))
  
  for(r in 1:length(Routes)){
    
    for(p in 1:length(Periods)){
      
      sim_covs <- sim_pc_y <- sim_pc_ydb <- sim_pc_doy <- sim_ARU <- vector(mode="list", length=sims)
      
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
        PC.avail <- doy_sim <- matrix(NA, ncol=Periods[p], nrow=nrow(surveydata_good))
        for(j in 1:Periods[p]){
          # Distributing surveys across weeks
          weeks_use <- which(RUGR.NGA.simulation.data$PC_periods[,Periods[p]]==j)
          random <- sample(1:length(goodroutes), replace = FALSE, length(goodroutes))
          weeks_route <- data.frame(Route=goodroutes,
                                    Week=weeks_use[ceiling(random/(length(goodroutes)/length(weeks_use)))])
          # Randomly picking dates within weeks for surveys
          weeks_route$DOY <- NA
          for(q in 1:nrow(weeks_route)){
            weeks_route$DOY[q] <- RUGR.NGA.simulation.data$PC_dates[weeks_route$Week[q],sample(1:ncol(RUGR.NGA.simulation.data$PC_dates),1)]
          } #q
          
          survey_data_merge <- merge(surveydata_good, weeks_route, by="Route", all.x=T)
          doy_sim[,j] <- survey_data_merge$DOY
          doy_sim_mat <- cbind(1,doy_sim[,j],doy_sim[,j]^2)
          PC.avail[,j] <- plogis(doy_sim_mat %*% c(RUGR.NGA.simulation.data$phi_params_PCintercept,RUGR.NGA.simulation.data$phi_params_doy))
        } #j
        
        # Simulating observation process for point counts
        ydb_RUGR <- array(NA, dim=c(nrow=nrow(point_dist_good), ncol=ncol(point_dist_good)-2,Periods[p]))
        for (i in 1:Periods[p]){
          for (j in 1:nrow(point_dist_good)){
            for (k in 1:(ncol(point_dist_good)-2)){
              RUGR_per_avail <- rbinom(1,point_dist_good[j,k+2],PC.avail[j,i]) # Number available
              ydb_RUGR[j,k,i] <- rbinom(1,RUGR_per_avail,band_det_prob[k]) # Number detected
            } #j
          } #k
        } #i
        # Total detections per survey
        y_RUGR <- apply(ydb_RUGR,3,rowSums)
        
        
        
        
        ## Simulating ARU data -------------------------------------------------
        
        
        # Availability on ARUs
        ARU.avail.5min <- plogis(cbind(1,RUGR.NGA.simulation.data$ARU_dates,RUGR.NGA.simulation.data$ARU_dates^2) %*% c(RUGR.NGA.simulation.data$phi_params_ARUintercept,RUGR.NGA.simulation.data$phi_params_doy))
        # This is availability for 5 minutes, changing to 3-hour recording
        ARU.avail <- 1-(1-ARU.avail.5min)^36
        
        # Simulating data for both methods
        ARU.det <- vector(mode="list",length=length(ARUs))
        
        for(a in 1:length(ARUs)){
          
          if(ARUs[a]==0){
            
            # If no ARUs, setting as NA
            ARU.det[[a]] <- matrix(NA, nrow=1, ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            
          } else{
            
            ### Deploying randomly at sites at which point counts are performed.
            ARU.sites.numbers <- sample(1:nrow(surveydata_good),ARUs[a],replace=F)
            
            # Using the methods of Royle & Nichols (2003), which models occupancy data
            #   as the probability of not detecting any of the N grouse at a location.
            # Calculating per bird detection probability as a function of availability and conditional detection
            ARU.dat <- matrix(NA, nrow=nrow(surveydata_good), ncol=length(RUGR.NGA.simulation.data$ARU_dates))
            for(i in 1:length(ARU.sites.numbers)){
              for(j in 1:ncol(ARU.dat)){
                pdetARU <- ARU.avail[j] * RUGR.NGA.simulation.data$pdet_ARU
                ARU.dat[ARU.sites.numbers[i],j] <- rbinom(1, 1, 1 - (1-pdetARU)^surveydata_good$N[ARU.sites.numbers[i]])
              } #i
            } #j
            ARU.det[[a]] <- ARU.dat
          } #else
        } #a
        
        sim_covs[[x]] <- surveydata_good
        sim_pc_y[[x]] <- y_RUGR
        sim_pc_ydb[[x]] <- ydb_RUGR
        sim_pc_doy[[x]] <- doy_sim
        sim_ARU[[x]] <- ARU.det
        
      } #x
      
      # Writing out
      sim_meta <- list(AvgDensity = Density.val[z],
                       TotAbund = RUGR.NGA.simulation.data$abund_N_GA_actual[z],
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
                       ARU.dates = RUGR.NGA.simulation.data$ARU_dates,
                       N_GA_covariates = RUGR.NGA.simulation.data$N_GA_covariates)
      sim_data <- list(sim_covs = sim_covs,
                       sim_pc_y = sim_pc_y,
                       sim_pc_ydb = sim_pc_ydb,
                       sim_pc_doy = sim_pc_doy,
                       sim_ARU = sim_ARU)
      
      sim_out <- list(sim_meta = sim_meta,
                      sim_data = sim_data)
      
      filename <- paste("SimulationDataRUGR_",Routes[r],"Routes_",Periods[p],"Periods_",Density.val[z],"avgdensity.gzip",sep="")
      save(sim_out,file=filename)
    } #p
  } #r
} #z
