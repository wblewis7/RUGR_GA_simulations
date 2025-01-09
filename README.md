# RUGR_GA_simulations
R and NIMBLE code for simulating and analyzing data to asses the efficacy of integrating Automated Recording Units (ARUs) and spring drumming point counts for monitoring populations of ruffed grouse (Bonasa umbellus).
----
authors: William B. Lewis, Gregory T. Wann, Richard B. Chandler, Mark D. McConnell, Emily Rushton, and James A. Martin


---

# Metadata

# RUGR_generating_simdata_NGA_code.R

R code for generating simulated monitoring data for ruffed grouse (Bonasa umbellus, RUGR) via traditional spring roadside drumming point count surveys, with or without deployment of Automated Recording Units (ARUs) across a 858,000 ha study area in Northeastern Georgia, USA. RUGR abundance in this area has been monitored using spring roadside drumming counts at 59 roadside survey routes, each consisting of 8 - 15 survey points (633 survey point total). Local abundance within 200 m of point count centers is simulated as a function of 5 landscape covariates: mean elevation, Lattitude, mean aspect (northness and eastness), and mean canopy height. The detection process is decomposed into two components (Amundson et al. 2014. A hierarchical model combining distance sampling and time removal to estimate detection probability during avian point counts): availability (the probability of a RUGR being within the survey area and drumming during a survey), and condtional detectability (the probability of detecting an available grouse). Availability is simulated based on a quadratic effect of day of year, peaking on April 19th. Both ARUs and point counts jointly estimate the linear and quadratic effects of day of year on availability but have separate intercepts. Point count detections within 20 m distance bands are simulated based on a Half-Normal detection function. ARUs are simulated as recording for 3 hours each morning from March 1 - May 31st, with the probability of detecting at least one drum on a 3 hour recording (conditional on availability) being 95%. ARU data is simulated via two methods. First, ARUs are randomly deployed at a series of point count survey locations to estimate occupancy and availability (Method OcAv). In the second method (Method Av), ARUs are preferentially deployed at survey locations with the highest expected grouse abundance (based on landscape covariates) to have the highest probability of being deployed at occupied sites in order to estimate availability. Since ARUs were not randomly deployed in this method, the data
shouldn't be used to assess occupancy and landscape covariate effects.

# RUGR.NGA.simulation.data.gzip

Parameters for generating simulated monitoring data for RUGR are stored in the 'RUGR.NGA.simulation.data.gzip' file. The parameters are for simulating montiroing data for RUGR  The study design and simulation parameters are based on a previous analysis of RUGR point count data (Lewis et al. 2022. Abundance and distribution of ruffed grouse Bonasa umbellus at the southern periphery of the range. Wildlife Biology) and a preliminary analysis of ARU data from the study region.
