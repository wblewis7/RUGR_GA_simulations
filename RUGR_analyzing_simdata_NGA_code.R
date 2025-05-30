################################################################################
## Code for analyzing simulated survey data for monitoring ruffed grouse (Bonasa umbellus)
##  across a 858,000 ha study area in North Georgia (GA), USA, from Clayton D. Delancey,
##  William B. Lewis, Gregory T. Wann, Richard B. Chandler, Emily Rushton, and James A. 
##  Martin. Assessing the utility of autonomous recording units and spring drumming 
##  surveys for monitoring abundance of ruffed grouse (Bonasa umbellus).
## Corresponding author: Will Lewis (wblewis7@gmail.com), University of Georgia.
################################################################################



# Monitoring data were simulated for ruffed grouse (RUGR) via spring roadside drumming 
#   point counts with or without deployment of Autonomous Recording Units (ARUs).
#   Grouse abundance in this area has been monitored using spring roadside drumming 
#   point counts along 59 roadside survey routes, each consisting of 8 - 15 survey points
#   (633 survey point total).
# Local abundance was simulated as a function of 5 landscape covariates, while the
#   detection process was simulated as a function of availability and conditional
#   detectability. Point counts and ARUs use a similar state process but have
#   different observation processes, though both jointly estimate the seasonal 
#   change in availability with drumming availability peaking on April 19th.
# Data simulated for monitoring strategies varying in the average density across
#   the study area (0.003, 0.006, 0.011, 0.022), the number of survey routes (20,
#   30, 40, 50), the number of replicate point count surveys (1 - 3), and the number
#   of survey points at which ARUs are deployed (0, 40, 80). The peak of drumming
#   activity can vary regionally and is often not known with certainty, so we modeled
#   surveys as being performed over a 12-week period (March 9 - May 31).
# Simulation parameters are taken from a previous analysis of point count data (Lewis et al. 
#   2022) and a preliminary analysis of ARU data from the study region.


# Data were simulated and saved out in the script RUGR_generating_simdata_NGA_code.R
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


# Setting an informative Uniform prior on sigma, the parameter informing the scale 
#   of the Half-Normal point-count detection process, from 90 - 130 (true value: 112.4).
#   In preliminary analyses, the model has trouble estimating this parameter with only 
#   a few detections and will sometimes gives extreme values (e.g., 200+)
# In the ARU submodel, there are not enough data sources to separate availability from conditional
#   detectability. Placing an informative Uniform prior on the detection process from 0.9
#   to 0.999999 (true value: 0.95).



require(nimble)
require(coda)


scenario.files <- list.files(pattern="SimulationDataRUGR")


for(q in 1:length(scenario.files)){
  
  load(scenario.files[q])
  
  monitor <-  c("a0","a1","a2","a3","a4","a5","b0.ARU","b0.PC","b1","b2","sigma")
  
  GR <- mcmc.samps <- GA.pred <- param.quants <- vector(mode="list", length=sim_out$sim_meta$Nsims)
  
  for(x in 1:sim_out$sim_meta$Nsims){
    
    GR.temp <- matrix(NA, nrow=length(sim_out$sim_meta$ARUs), ncol=length(sim_out$sim_meta$lambda_params)+length(sim_out$sim_meta$phi.ARUintercept.params)+length(sim_out$sim_meta$phi.PCintercept.params)+length(sim_out$sim_meta$phi.doy.params)+length(sim_out$sim_meta$p.PC.sigma))
    mcmc.samps.temp <- param.quants.temp <- vector(mode="list", length=length(sim_out$sim_meta$ARUs))
    GA_mat <- sim_out$sim_meta$N_GA_covariates[,c("elev","Latitude","eastness","northness","canheight","area")]
    GA_mat[,ncol(GA_mat)] <- log(GA_mat[,ncol(GA_mat)])
    GA_mat <- as.matrix(cbind(matrix(rep(1, times=nrow(GA_mat), ncol=1)),
                    GA_mat))
    GA.pred.temp <- matrix(NA, nrow=length(sim_out$sim_meta$ARUs), ncol=11400)
    
    
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
                             PC.doy = as.matrix(sim_out$sim_data$sim_pc_doy[[x]]),
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
        
        # Convergence
        GR.temp[a,] <- gelman.diag(RUGR.mcmc.out)$psrf[,1]
        
        # Posterior summaries
        mc.noARU <- rbind(RUGR.mcmc.out[[1]], RUGR.mcmc.out[[2]], RUGR.mcmc.out[[3]])
        param.quants.temp[[a]] <- apply(mc.noARU, 2, quantile, probs=c(0.025,0.5,0.975))
        
        # Predicting abundance across the study area
        samps.pred <- cbind(mc.noARU[,1:6],
                            matrix(rep(1, times=nrow(mc.noARU), ncol=1)))
        GA.pred.temp[a,] <- colSums(exp(GA_mat %*% t(samps.pred)))
        
        # Saving subset of posterior samples
        mcmc.samps.temp[[a]] <- mc.noARU[sample(1:nrow(mc.noARU), 1000, replace=FALSE),]
        
        
        
      } else{
        
        # ARU model -------------------------------------------------------
        RUGR_nimble_ARU <- nimbleCode({
          
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
          
          # ARU availability and occupancy submodel
          for (d in 1:deployment.dates.n){
            
            # Availability parameter based on day of year, have to convert between 5 minute count and 3 hour recording
            logit(avail.ARU.5min[d]) <- b0.ARU + b1 * ARU.doy[d] + b2 * pow(ARU.doy[d],2)
            avail.ARU[d] <- 1 - pow(1-avail.ARU.5min[d], timeadj)
            
            for (a in 1:n.whichARU){
              
              # Detection on a day is function of abundance, availability, and probability that at least one drum is detected, which should be high
              ARUdetects[whichARU[a],d] ~ dbern(1 - ((1 - (avail.ARU[d] * p.det.ARU))^N[whichARU[a]]))
              
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
        N.init_ARU <- rep(max(sim_out$sim_data$sim_pc_y[[x]])+1,times=nrow(sim_out$sim_data$sim_pc_y[[x]]))
        init.function_ARU <- function() list(a0=runif(1,-12,-2),
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
        
        mod.data_ARU <- list(y = sim_out$sim_data$sim_pc_y[[x]],
                            ydb = sim_out$sim_data$sim_pc_ydb[[x]],
                            ARUdetects = sim_out$sim_data$sim_ARU[[x]][[a]])
        
        mod.constants_ARU <- list(nSites = nrow(sim_out$sim_data$sim_pc_y[[x]]),
                                 elev = sim_out$sim_data$sim_covs[[x]]$elev,
                                 lat = sim_out$sim_data$sim_covs[[x]]$Latitude,
                                 east = sim_out$sim_data$sim_covs[[x]]$eastness,
                                 north = sim_out$sim_data$sim_covs[[x]]$northness,
                                 canheight = sim_out$sim_data$sim_covs[[x]]$canheight,
                                 area = sim_out$sim_meta$PC_detbands$survey.area,
                                 nVisits = sim_out$sim_meta$Periods,
                                 PC.doy = as.matrix(sim_out$sim_data$sim_pc_doy[[x]]),
                                 deployment.dates.n = length(sim_out$sim_meta$ARU.dates),
                                 ARU.doy = sim_out$sim_meta$ARU.dates,
                                 timeadj = 3*60/5,
                                 whichARU = which(!is.na(sim_out$sim_data$sim_ARU[[x]][[a]][,1])),
                                 n.whichARU = length(which(!is.na(sim_out$sim_data$sim_ARU[[x]][[a]][,1]))),
                                 nB = sim_out$sim_meta$PC_detbands$Nbands,
                                 dB = sim_out$sim_meta$PC_detbands$det_bands,
                                 pix = sim_out$sim_meta$PC_detbands$rel_band_area)
        
        RUGR.mod.ARU <- nimbleModel(code = RUGR_nimble_ARU,
                                   data = mod.data_ARU,
                                   constants = mod.constants_ARU,
                                   inits = init.function_ARU())
        RUGR.mcmc.out.ARU  <- nimbleMCMC(model = RUGR.mod.ARU, 
                                        niter = 20000, nchains = 3, nburnin = 1000,
                                        monitor=monitor, thin=5, samplesAsCodaMCMC=TRUE)
        
        # Convergence
        GR.temp[a,] <- gelman.diag(RUGR.mcmc.out.ARU)$psrf[,1]
        
        # Posterior summaries
        mc.ARU <- rbind(RUGR.mcmc.out.ARU[[1]], RUGR.mcmc.out.ARU[[2]], RUGR.mcmc.out.ARU[[3]])
        param.quants.temp[[a]] <- apply(mc.ARU, 2, quantile, probs=c(0.025,0.5,0.975))
        
        # Predicting abundance across the study area
        samps.pred <- cbind(mc.ARU[,1:6],
                            matrix(rep(1, times=nrow(mc.ARU), ncol=1)))
        GA.pred.temp[a,] <- colSums(exp(GA_mat %*% t(samps.pred)))
        
        # Saving subset of posterior samples
        mcmc.samps.temp[[a]] <- mc.ARU[sample(1:nrow(mc.ARU), 1000, replace=FALSE),]
        
      } # else
    } #a
    
    mcmc.samps[[x]] <- mcmc.samps.temp
    GA.pred[[x]] <- GA.pred.temp
    GR[[x]] <- GR.temp
    param.quants[[x]] <- param.quants.temp
    
  } #x
  
  mcmc_out <- list(mcmc.samps = mcmc.samps,
                   GR = GR,
                   GA.pred = GA.pred,
                   param.quants = param.quants)
  sim_mcmc_out <- list(sim_meta = sim_out$sim_meta,
                       mcmc_out = mcmc_out)
  mcmc.name.file.sim <- paste("RUGRmcmcsamps_",sim_out$sim_meta$Routes,"Routes_",sim_out$sim_meta$Periods,"Periods_",sim_out$sim_meta$AvgDensity,"Density.gzip", sep="")
  save(sim_mcmc_out, file=mcmc.name.file.sim)
  
} #z
