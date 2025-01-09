################################################################################
## Code for analyzing simulated survey data for monitoring ruffed grouse (Bonasa umbellus)
##  across a 858,000 ha study area in North Georgia (GA), USA, from Clayton D. Delancey,
##  William B. Lewis, Gregory T. Wann, Richard B. Chandler, Mark D. McConnell, Emily
##  Rushton, and James A. Martin. Assessing the utility of Automated Recording Units 
##  and spring drumming surveys for monitoring abundance of Ruffed Grouse (Bonasa umbellus).
## Corresponding author: Will Lewis (wblewis7@gmail.com), University of Georgia.
################################################################################



# Monitoring data simulated for ruffed grouse (RUGR) via spring roadside drumming 
#   point counts with or without deployment of Automated Recording Units (ARUs).
#   Grouse abundance in this area has been monitored using spring roadside drumming 
#   counts at 59 roadside survey routes, each consisting of 8 - 15 survey points
#   (633 survey point total).
# Local abundance was simulated as a function of 5 landscape covariates, while the
#   detection process was simulated as a function of availability and conditional
#   detectability. Point counts and ARUs use a similar state process but have
#   different observation processes, though both jointly estimate the seasonal 
#   change in availability with drumming availability peaking on April 19th.
# Simulated data for monitoring strategies varying in the average density across
#   the study area (0.003, 0.006, 0.011, 0.022), the number of survey routes (10, 20,
#   30, 40, 50), the number of replicate point count surveys (1 - 3), and the number
#   of survey points at which ARUs were deployed (0, 40, 80). ARU deployment was 
#   simulated in two ways. First, ARUs were randomly deployed at 40/80 point count 
#   locations to assess occupancy and availability (Method OcAv). Estimated availability
#   requires detections, so the second method (Method Av) preferentially deployed
#   ARUs at survey locations with the highest likelihood of detections (i.e., highest 
#   expected abundance). Since ARUs were not randomly deployed in this method, the data
#   shouldn't be used to assess occupancy and landscape covariate effects.
# Simulation parameters were taken from a previous analysis of point count data (Lewis et al. 
#   2022. Abundance and distribution of ruffed grouse Bonasa umbellus at the southern 
#   periphery of the range) and a preliminary analysis of ARU data from the study region.


# Data was simulated and saved out in the script RUGR_generating_simdata_NGA_code.R
# Analyzing using Bayesian hierarchical models in NIMBLE. The basis of the model is
#   the distance-sampling N-mixture model of Royle et al (2004, Modeling abundance
#   effects in distance sampling) modified to decompose the detection process into
#   availability and conditional detectability and to incorporate ARUs.
# Setting an informative Uniform prior on sigma, the parameter informing the scale 
#   of the Half-Normal point-count detection process, from 90 - 130 (true value: 112.4).
#   In preliminary analyses, the model has trouble estimating this parameter if only 
#   a few detections and will sometimes gives extreme values (e.g., 200+)
# In the ARU submodel, there is not enough data to separate availability from conditional
#   detectability. Placing an informative Uniform prior on the detection process from 0.9
#   to 0.999999 (true value: 0.95).
# Going to save a subset of parameter estimates (to space) and assessing convergence based
#   on the Gelman-Rubin diagnostic (Gelman and Rubin. 1992. Inference from iterative
#   simulation using multiple sequences).


require(nimble)
require(coda)


scenario.files <- list.files(pattern="SimulationDataRUGR")


for(q in 1:length(scenario.files)){
  
  load(scenario.files[q])
  
  monitor <-  c("a0","a1","a2","a3","a4","a5","b0.ARU","b0.PC","b1","b2","sigma")
  
  GR.noARU <- GR.ARU.Av <- GR.ARU.OcAv <- array(NA, dim=c(sim_out$sim_meta$Nsims, length(monitor), length(sim_out$sim_meta$ARUs)))
  mcmc.samps.noARU <- mcmc.samps.ARU.Av <- mcmc.samps.ARU.OcAv <- vector(mode="list", length=sim_out$sim_meta$Nsims)
  
  for(x in 1:sim_out$sim_meta$Nsims){
    
    mcmc.samps.noARU.temp <- mcmc.samps.ARU.Av.temp <- mcmc.samps.ARU.OcAv.temp <- vector(mode="list", length=length(sim_out$sim_meta$ARUs))
    
    for(a in 1:length(sim_out$sim_meta$ARUs)){
      
      if(sim_out$sim_meta$ARUs[a] == 0){
        
        # Model without ARUs ---------------------------------------------------
        RUGR_nimble <- nimbleCode({
          
          # Priors
          a0 ~ dnorm(0,sd=10)
          a1 ~ dnorm(0,sd=10)
          a2 ~ dnorm(0,sd=10)
          a3 ~ dnorm(0,sd=10)
          a4 ~ dnorm(0,sd=10)
          a5 ~ dnorm(0,sd=10)
          b0.PC ~ dnorm(0,sd=10)
          b0.ARU ~ dnorm(0,sd=10)
          b1 ~ dnorm(0,sd=10)
          b2 ~ dnorm(0,sd=10)
          sigma ~ dunif(90,130)
          
          for(z in 1:nSites){
            
            # State process
            EN[z] <- exp(a0 + a1 * elev[z] + a2 * lat[z] + a3 * east[z] + a4 * north[z] + a5 * canheight[z] + log(area))
            N[z] ~ dpois(EN[z])
            
            # Detection process for point counts
            for (v in 1:nVisits){
              
              # Availability parameter based on day of year
              logit(avail[z,v]) <- b0.PC + b1 * PC.doy[z,v] + b2 * pow(PC.doy[z,v],2)
              # Number of detections on a survey are conditional on N, avail, and pd
              y[z,v] ~ dbin(avail[z,v] * pd, N[z])
              # Bands of detections of observed birds are conditional on detection probability for each band and number of detected individuals
              ydb[z,1:nB,v] ~ dmulti(pd_bins_adj[1:nB], y[z,v])
              
            } #v
          } # z
          
          # Detection probability for point counts estimated in each distance band based on Half-Normal detection function
          for(b in 1:nB){
            
            pd_bins[b] <- (sigma^2*(1-exp(-dB[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-dB[b]^2/(2*sigma^2))))*2*3.1416/(area*10000*pix[b])
            pd_bins_adj[b] <- pd_bins[b]*pix[b]
            
          } #b
          
          # Overall detection probability in sampling area (constant across sites)
          pd <- sum(pd_bins_adj[1:nB])
          
        })
        
        # Setting initial values
        # Making sure initial values for N are higher than number of birds detected
        N.init <- rep(max(sim_out$sim_data$sim_pc_y[[x]])+1,times=nrow(sim_out$sim_data$sim_pc_y[[x]]))
        init.function <- function() list(a0=runif(1,-12,-2),
                                         a1=runif(1,-1,1),
                                         a2=runif(1,-1,1),
                                         a3=-runif(1,-1,1),
                                         a4=runif(1,-1,1),
                                         a5=runif(1,-1,1),
                                         b0.PC=runif(1,-1,1),
                                         b0.ARU=runif(1,-1,1),
                                         b1=runif(1,-1,1),
                                         b2=runif(1,-1,1),
                                         sigma=runif(1,90,130),
                                         N=N.init)
        
        mod.data <- list(y = sim_out$sim_data$sim_pc_y[[x]],
                         ydb = sim_out$sim_data$sim_pc_ydb[[x]])
        
        mod.constants <- list(nSites = nrow(sim_out$sim_data$sim_pc_y[[x]]),
                             elev = sim_out$sim_data$sim_covs[[x]]$elev,
                             lat = sim_out$sim_data$sim_covs[[x]]$Latitude,
                             east = sim_out$sim_data$sim_covs[[x]]$eastness,
                             north = sim_out$sim_data$sim_covs[[x]]$northness,
                             canheight = sim_out$sim_data$sim_covs[[x]]$canheight,
                             area = sim_out$sim_meta$PC_detbands$survey.area,
                             nVisits = sim_out$sim_meta$Periods,
                             PC.doy = as.matrix(sim_out$sim_data$sim_pc_doy[[x]][,3:ncol(sim_out$sim_data$sim_pc_doy[[x]])]),
                             nB = sim_out$sim_meta$PC_detbands$Nbands,
                             dB = sim_out$sim_meta$PC_detbands$det_bands,
                             pix = sim_out$sim_meta$PC_detbands$rel_band_area)
        
        RUGR.mod <- nimbleModel(code = RUGR_nimble,
                                data = mod.data,
                                constants = mod.constants,
                                inits = init.function())
        RUGR.mcmc.out  <- nimbleMCMC(model = RUGR.mod, 
                                    niter = 20000, nchains = 3, nburnin = 1000,
                                    monitor=monitor, thin=5, samplesAsCodaMCMC=TRUE)
        
        # Subsetting posterior samples and convergence
        GR.noARU[x,1:length(monitor),a] <- gelman.diag(RUGR.mcmc.out)$psrf[,1]
        mc.noARU <- rbind(RUGR.mcmc.out[[1]], RUGR.mcmc.out[[2]], RUGR.mcmc.out[[3]])
        mcmc.samps.noARU.temp[[a]] <- mc.noARU[sample(1:nrow(mc.noARU), 1000, replace=FALSE),]
        
      } else{
        
        # ARU model Av -------------------------------------------------------
        RUGR_nimble_Av <- nimbleCode({
          
          # Priors
          a0 ~ dnorm(0,sd=10)
          a1 ~ dnorm(0,sd=10)
          a2 ~ dnorm(0,sd=10)
          a3 ~ dnorm(0,sd=10)
          a4 ~ dnorm(0,sd=10)
          a5 ~ dnorm(0,sd=10)
          b0.PC ~ dnorm(0,sd=10)
          b0.ARU ~ dnorm(0,sd=10)
          b1 ~ dnorm(0,sd=10)
          b2 ~ dnorm(0,sd=10)
          sigma ~ dunif(90,130)
          p.det.ARU ~ dunif(0.9,0.99999)
          
          for(z in 1:nSites){
            
            # State process
            EN[z] <- exp(a0 + a1 * elev[z] + a2 * lat[z] + a3 * east[z] + a4 * north[z] + a5 * canheight[z] + log(area))
            N[z] ~ dpois(EN[z])
            
            # Detection process for point counts
            for (v in 1:nVisits){
              
              # Availability parameter based on day of year
              logit(avail[z,v]) <- b0.PC + b1 * PC.doy[z,v] + b2 * pow(PC.doy[z,v],2)
              # Number of detections on a survey are conditional on N, avail, and pd
              y[z,v] ~ dbin(avail[z,v] * pd, N[z])
              # Bands of detections of observed birds are conditional on detection probability for each band and number of detected individuals
              ydb[z,1:nB,v] ~ dmulti(pd_bins_adj[1:nB], y[z,v])
              
            } #v
          } # z
          
          # ARU availability submodel
          for (d in 1:deployment.dates.n){
            
            # Availability parameter based on day of year, have to convert between 5 minute count and 3 hour recording
            logit(avail.ARU.5min[d]) <- b0.ARU + b1 * ARU.doy[d] + b2 * pow(ARU.doy[d],2)
            avail.ARU[d] <- 1 - pow(1-avail.ARU.5min[d], timeadj)
            
            for (a in 1:NARUs){
              
              # Detection on a day is function of availability and probability that at least one drum is detected, which should be high
              # Already subset to sites with detections. There are no simulated false positives, so know that this means all sites in likelihood are occupied
              ARUdetects[a,d] ~ dbern(avail.ARU[d] * p.det.ARU)
              
            } #a        
          } #d
          
          # Detection probability for point counts estimated in each distance band based on Half-Normal detection function
          for(b in 1:nB){
            
            pd_bins[b] <- (sigma^2*(1-exp(-dB[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-dB[b]^2/(2*sigma^2))))*2*3.1416/(area*10000*pix[b])
            pd_bins_adj[b] <- pd_bins[b]*pix[b]
            
          } #b
          
          # Overall detection probability in sampling area (constant across sites)
          pd <- sum(pd_bins_adj[1:nB])
          
        })
        
        # Setting initial values
        # Making sure initial values for N are higher than number of birds detected
        N.init_Av <- rep(max(sim_out$sim_data$sim_pc_y[[x]])+1,times=nrow(sim_out$sim_data$sim_pc_y[[x]]))
        init.function_Av <- function() list(a0=runif(1,-12,-2),
                                         a1=runif(1,-1,1),
                                         a2=runif(1,-1,1),
                                         a3=-runif(1,-1,1),
                                         a4=runif(1,-1,1),
                                         a5=runif(1,-1,1),
                                         b0.PC=runif(1,-1,1),
                                         b0.ARU=runif(1,-1,1),
                                         b1=runif(1,-1,1),
                                         b2=runif(1,-1,1),
                                         sigma=runif(1,90,130),
                                         p.det.ARU=runif(1,0.9,0.999),
                                         N=N.init)
        
        mod.data_Av <- list(y = sim_out$sim_data$sim_pc_y[[x]],
                         ydb = sim_out$sim_data$sim_pc_ydb[[x]],
                         ARUdetects = matrix(sim_out$sim_data$sim_ARU_Av[[x]][[a]], ncol=length(sim_out$sim_meta$ARU.dates)))
        
        mod.constants_Av <- list(nSites = nrow(sim_out$sim_data$sim_pc_y[[x]]),
                              elev = sim_out$sim_data$sim_covs[[x]]$elev,
                              lat = sim_out$sim_data$sim_covs[[x]]$Latitude,
                              east = sim_out$sim_data$sim_covs[[x]]$eastness,
                              north = sim_out$sim_data$sim_covs[[x]]$northness,
                              canheight = sim_out$sim_data$sim_covs[[x]]$canheight,
                              area = sim_out$sim_meta$PC_detbands$survey.area,
                              nVisits = sim_out$sim_meta$Periods,
                              PC.doy = as.matrix(sim_out$sim_data$sim_pc_doy[[x]][,3:ncol(sim_out$sim_data$sim_pc_doy[[x]])]),
                              deployment.dates.n = length(sim_out$sim_meta$ARU.dates),
                              ARU.doy = sim_out$sim_meta$ARU.dates,
                              timeadj = 3*60/5,
                              NARUs = nrow(matrix(sim_out$sim_data$sim_ARU_Av[[x]][[a]], ncol=length(sim_out$sim_meta$ARU.dates))),
                              nB = sim_out$sim_meta$PC_detbands$Nbands,
                              dB = sim_out$sim_meta$PC_detbands$det_bands,
                              pix = sim_out$sim_meta$PC_detbands$rel_band_area)
        
        RUGR.mod.Av <- nimbleModel(code = RUGR_nimble_Av,
                                data = mod.data_Av,
                                constants = mod.constants_Av,
                                inits = init.function_Av())
        RUGR.mcmc.out.Av  <- nimbleMCMC(model = RUGR.mod.Av, 
                                     niter = 20000, nchains = 3, nburnin = 1000,
                                     monitor=monitor, thin=5, samplesAsCodaMCMC=TRUE)
        
        # Subsetting posterior samples and convergence
        GR.ARU.Av[x,1:length(monitor),a] <- gelman.diag(RUGR.mcmc.out.Av)$psrf[,1]
        mc.ARU.Av <- rbind(RUGR.mcmc.out.Av[[1]], RUGR.mcmc.out.Av[[2]], RUGR.mcmc.out.Av[[3]])
        mcmc.samps.ARU.Av.temp[[a]] <- mc.ARU.Av[sample(1:nrow(mc.ARU.Av), 1000, replace=FALSE),]
        
        
        
        # ARU model OcAv -------------------------------------------------------
        RUGR_nimble_OcAv <- nimbleCode({
          
          # Priors
          a0 ~ dnorm(0,sd=10)
          a1 ~ dnorm(0,sd=10)
          a2 ~ dnorm(0,sd=10)
          a3 ~ dnorm(0,sd=10)
          a4 ~ dnorm(0,sd=10)
          a5 ~ dnorm(0,sd=10)
          b0.PC ~ dnorm(0,sd=10)
          b0.ARU ~ dnorm(0,sd=10)
          b1 ~ dnorm(0,sd=10)
          b2 ~ dnorm(0,sd=10)
          sigma ~ dunif(90,130)
          p.det.ARU ~ dunif(0.9,0.99999)
          
          for(z in 1:nSites){
            
            # State process
            EN[z] <- exp(a0 + a1 * elev[z] + a2 * lat[z] + a3 * east[z] + a4 * north[z] + a5 * canheight[z] + log(area))
            N[z] ~ dpois(EN[z])
            Occ[z] <- N[z] > 0
            
            # Detection process for point counts
            for (v in 1:nVisits){
              
              # Availability parameter based on day of year
              logit(avail[z,v]) <- b0.PC + b1 * PC.doy[z,v] + b2 * pow(PC.doy[z,v],2)
              # Number of detections on a survey are conditional on N, avail, and pd
              y[z,v] ~ dbin(avail[z,v] * pd, N[z])
              # Bands of detections of observed birds are conditional on detection probability for each band and number of detected individuals
              ydb[z,1:nB,v] ~ dmulti(pd_bins_adj[1:nB], y[z,v])
              
            } #v
          } # z
          
          # ARU availability and occupancy submodel
          for (d in 1:deployment.dates.n){
            
            # Availability parameter based on day of year, have to convert between 5 minute count and 3 hour recording
            logit(avail.ARU.5min[d]) <- b0.ARU + b1 * ARU.doy[d] + b2 * pow(ARU.doy[d],2)
            avail.ARU[d] <- 1 - pow(1-avail.ARU.5min[d], timeadj)
            
            for (a in 1:n.whichARU){
              
              # Detection on a day is function of occupancy, availability, and probability that at least one drum is detected, which should be high
              ARUdetects[whichARU[a],d] ~ dbern(Occ[whichARU[a]] * avail.ARU[d] * p.det.ARU)
              
            } #a        
          } #d
          
          # Detection probability for point counts estimated in each distance band based on Half-Normal detection function
          for(b in 1:nB){
            
            pd_bins[b] <- (sigma^2*(1-exp(-dB[b+1]^2/(2*sigma^2)))-sigma^2*(1-exp(-dB[b]^2/(2*sigma^2))))*2*3.1416/(area*10000*pix[b])
            pd_bins_adj[b] <- pd_bins[b]*pix[b]
            
          } #b
          
          # Overall detection probability in sampling area (constant across sites)
          pd <- sum(pd_bins_adj[1:nB])
          
        })
        
        # Setting initial values
        # Making sure initial values for N are higher than number of birds detected
        N.init_OcAv <- rep(max(sim_out$sim_data$sim_pc_y[[x]])+1,times=nrow(sim_out$sim_data$sim_pc_y[[x]]))
        init.function_OcAv <- function() list(a0=runif(1,-12,-2),
                                            a1=runif(1,-1,1),
                                            a2=runif(1,-1,1),
                                            a3=-runif(1,-1,1),
                                            a4=runif(1,-1,1),
                                            a5=runif(1,-1,1),
                                            b0.PC=runif(1,-1,1),
                                            b0.ARU=runif(1,-1,1),
                                            b1=runif(1,-1,1),
                                            b2=runif(1,-1,1),
                                            sigma=runif(1,90,130),
                                            p.det.ARU=runif(1,0.9,0.999),
                                            N=N.init)
        
        mod.data_OcAv <- list(y = sim_out$sim_data$sim_pc_y[[x]],
                            ydb = sim_out$sim_data$sim_pc_ydb[[x]],
                            ARUdetects = sim_out$sim_data$sim_ARU_OcAv[[x]][[a]])
        
        mod.constants_OcAv <- list(nSites = nrow(sim_out$sim_data$sim_pc_y[[x]]),
                                 elev = sim_out$sim_data$sim_covs[[x]]$elev,
                                 lat = sim_out$sim_data$sim_covs[[x]]$Latitude,
                                 east = sim_out$sim_data$sim_covs[[x]]$eastness,
                                 north = sim_out$sim_data$sim_covs[[x]]$northness,
                                 canheight = sim_out$sim_data$sim_covs[[x]]$canheight,
                                 area = sim_out$sim_meta$PC_detbands$survey.area,
                                 nVisits = sim_out$sim_meta$Periods,
                                 PC.doy = as.matrix(sim_out$sim_data$sim_pc_doy[[x]][,3:ncol(sim_out$sim_data$sim_pc_doy[[x]])]),
                                 deployment.dates.n = length(sim_out$sim_meta$ARU.dates),
                                 ARU.doy = sim_out$sim_meta$ARU.dates,
                                 timeadj = 3*60/5,
                                 whichARU = which(!is.na(sim_out$sim_data$sim_ARU_OcAv[[x]][[a]][,1])),
                                 n.whichARU = length(which(!is.na(sim_out$sim_data$sim_ARU_OcAv[[x]][[a]][,1]))),
                                 nB = sim_out$sim_meta$PC_detbands$Nbands,
                                 dB = sim_out$sim_meta$PC_detbands$det_bands,
                                 pix = sim_out$sim_meta$PC_detbands$rel_band_area)
        
        RUGR.mod.OcAv <- nimbleModel(code = RUGR_nimble_OcAv,
                                   data = mod.data_OcAv,
                                   constants = mod.constants_OcAv,
                                   inits = init.function_OcAv())
        RUGR.mcmc.out.OcAv  <- nimbleMCMC(model = RUGR.mod.OcAv, 
                                        niter = 20000, nchains = 3, nburnin = 1000,
                                        monitor=monitor, thin=5, samplesAsCodaMCMC=TRUE)
        
        # Subsetting posterior samples and convergence
        GR.ARU.OcAv[x,1:length(monitor),a] <- gelman.diag(RUGR.mcmc.out.OcAv)$psrf[,1]
        mc.ARU.OcAv <- rbind(RUGR.mcmc.out.OcAv[[1]], RUGR.mcmc.out.OcAv[[2]], RUGR.mcmc.out.OcAv[[3]])
        mcmc.samps.ARU.OcAv.temp[[a]] <- mc.ARU.OcAv[sample(1:nrow(mc.ARU.OcAv), 1000, replace=FALSE),]
        
      } # else
    } #a
    
    mcmc.samps.noARU[[x]] <- mcmc.samps.noARU.temp
    mcmc.samps.ARU.Av[[x]] <- mcmc.samps.ARU.Av.temp
    mcmc.samps.ARU.OcAv[[x]] <- mcmc.samps.ARU.OcAv.temp
    
  } #x
  
  mcmc_out <- list(mcmc.samps = list(noARU = mcmc.samps.noARU,
                                     ARU.Av = mcmc.samps.ARU.Av,
                                     ARU.OcAv = mcmc.samps.ARU.OcAv),
                   GR = list(noARU = GR.noARU,
                             ARU.Av = GR.ARU.Av,
                             ARU.OcAv = GR.ARU.OcAv))
  sim_mcmc_out <- list(sim_meta = sim_out$sim_meta,
                       mcmc_out = mcmc_out)
  mcmc.name.file.sim <- paste("RUGRmcmcsamps_",sim_out$sim_meta$Routes,"Routes_",sim_out$sim_meta$Periods,"Periods_",sim_out$sim_meta$AvgDensity,"Density.gzip", sep="")
  save(sim_mcmc_out, file=mcmc.name.file.sim)
  
} #z
