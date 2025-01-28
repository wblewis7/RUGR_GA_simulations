# RUGR_GA_simulations
R and NIMBLE code for: Assessing the utility of autonomous recording units and spring drumming surveys for monitoring abundance of ruffed grouse (Bonasa umbellus).
----
authors: Clayton D. Delancey, William B. Lewis*, Gregory T. Wann, Richard B. Chandler, Emily Rushton, and James A. Martin
* corresonding author: William B. Lewis, University of Georgia, wblewis7@gmail.com

---

# Metadata

# RUGR_generating_simdata_NGA_code.R

R code for generating simulated monitoring data for ruffed grouse (Bonasa umbellus, RUGR) via traditional spring roadside drumming point count surveys, with or without deployment of Automated Recording Units (ARUs) across a 858,000 ha study area in North Georgia, USA. RUGR abundance in this area has been monitored using spring roadside drumming counts at 59 roadside survey routes, each consisting of 8 - 15 survey points (633 survey point total). Local abundance within 200 m of point count centers is simulated as a function of 5 landscape covariates: mean elevation, Lattitude, mean aspect (northness and eastness), and mean canopy height. The detection process is decomposed into two components (Amundson et al. 2014. A hierarchical model combining distance sampling and time removal to estimate detection probability during avian point counts): availability (the probability of a RUGR being within the survey area and drumming during a survey), and condtional detectability (the probability of detecting an available grouse). Availability is simulated based on a quadratic effect of day of year, peaking on April 19th. Both ARUs and point counts jointly estimate the linear and quadratic effects of day of year on availability but have separate intercepts. Point count detections within 20 m distance bands are simulated based on a Half-Normal detection function. ARUs are simulated as recording for 3 hours each morning from March 1 - May 31st, with the data containing a binary response of whether or not at least one drum is detected on a recording. Preliminary analyses suggest that the probability of detecting at least one drum should be > 95% if available grouse drum at least 15 times total during a 3-hour recording. ARU data is simulated via two methods. First, ARUs are randomly deployed at a series of point count survey locations to estimate occupancy and availability (Method OcAv). In the second method (Method Av), ARUs are preferentially deployed at survey locations with the highest expected grouse abundance (based on landscape covariates) to have the highest probability of being deployed at occupied sites in order to estimate availability. Since ARUs are not randomly deployed in this method, the data shouldn't be used to assess occupancy and landscape covariate effects. Data are simulated under four different levels of grouse densities, implemented by modifiying the intercept of the abundance submodel. Predicting across the entire North Georgia study area, these intercepts correspond to total male population abundance of 2675, 5116, 9785, and 18714, corresponding to average male grouse densities of 0.003, 0.006, 0.011, and 0.022 males/ha, respectively. Data are also simulated under different monitoring scenarios varying the number of survey routes (10, 20, 30, 40, or 50), the number of replicate point count surveys (1 - 3), and the number of survey points at which ARUs are deployed (0, 40, 80).

# RUGR_analyzing_simdata_NGA_code.R

R and NIMBLE code for analyzing simulated RUGR survey data generated in RUGR_generating_simdata_NGA_code.R. Point count and ARU data are analyzed using Bayesian hierarchical models. The basis of the model is the distance-sampling N-mixture model of Royle et al (2004, Modeling abundance effects in distance sampling) modified to decompose the detection process into availability and conditional detectability and to incorporate ARUs.

# RUGR.NGA.simulation.data.gzip

Parameters for generating simulated monitoring data for RUGR are stored in the 'RUGR.NGA.simulation.data.gzip' file. The study design and simulation parameters are based on a previous analysis of RUGR point count data (Lewis et al. 2022. Abundance and distribution of ruffed grouse Bonasa umbellus at the southern periphery of the range) and a preliminary analysis of ARU data from the study region.
### surveydata
A matrix giving the landscape covariates for each of 633 survey locations in North Georgia. Suvey locations (Point) are organized into 59 roadside survey routes (Route). Survey-level landscape covariates are mean elevation (elev), Latitude, aspect (decomposed to eastness and northness), and canopy height (canheight) within 200 m of the point count center. All landscape covariates are standardized.
### lambda_params_intercepts
A vector giving the intercepts for the abundance modeling. Data are simulated under four different levels of abundance, corresponding to average RUGR densities of 0.003, 0.006, 0.011, and 0.022 males/ha, respectively, across the study area.
### lambda_params_cov
A vector giving the covariate effects of elevation, Latitude, eastness, northness, and canopy height (respectively) on RUGR abundance.
### phi_params_PCintercept
The intercept for modeling the availability portion of the point count observation process.
### phi_params_ARUintercept
The intercept for modeling the availability portion of the ARU observation process. This parameter represents the availability for a 5-minute segment of an ARU recording to be comparable with phi_params_PCintercept and allow for similar modeling of the effects of day of year on availability (phi_params_doy). The estimate of 5-minute availability will be converted to availability for a 3-hour ARU recording.
### phi_params_doy
A vector giving the linear and quadratic effects of day of year, respectively, on the availability portion of the observation process. These parameters are shared between the point count and ARU observation processes.
### p_params_PC
The sigma parameter for the Half-Normal detection function for the conditional detectability portion of the point count observation process. Sigma informs the scale of the relationship between distance and detection probability.
### pdet_ARU
The probability of detecting at least one RUGR drum during a 3-hour recording, conditional on availability (0.95).
### ARU_dates
A vector giving the day of year (standardized) of ARU recordings.
### PC_period_dates
A 3-element list giving the range of survey dates for point counts replicated either 1 (first element), 2 (second element), or 3 (third element) times. Within each element, rows represent sampling periods while columns represent the first and last days of the sampling period (respectively). Each point count sampling period lasts 2 weeks, with periods centered on the date of peak drumming activity in the region (April 19th). Surveys occurred from from April 12th – April 25th for scenarios with 1 survey, from April 5th – April 18th (1st survey) and April 19th – May 2nd (2nd survey) for scenarios with 2 replicate surveys, and from March 28th – April 11th (1st survey), April 12th – April 25th (2nd survey), and April 26th – May 9th (3rd survey) for scenarios with 3 replicate surveys. Day of year values are standardized.
### surveyarea
The area in hectares surveyed on each point count (200 m fixed-radius surveys).
### PC_det_bands
A vector giving the boundaries of 20 m distance bands (out to 200 m) in which detected RUGR are binned into during point count surveys. 
