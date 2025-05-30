# RUGR_GA_simulations
R and NIMBLE code for: Assessing the utility of autonomous recording units and spring drumming surveys for monitoring abundance of ruffed grouse (Bonasa umbellus).
----
authors: Clayton D. Delancey, William B. Lewis*, Gregory T. Wann, Richard B. Chandler, Emily Rushton, and James A. Martin
* corresonding author: William B. Lewis, University of Georgia, wblewis7@gmail.com

---

# Metadata

# RUGR_generating_simdata_NGA_code.R

R code for generating simulated monitoring data for ruffed grouse (Bonasa umbellus, RUGR) via traditional spring roadside drumming point count surveys, with or without deployment of Autonomous Recording Units (ARUs) across a 858,000 ha study area in North Georgia, USA. The study design is based off ongoing RUGR abundance monitoring via spring roadside drumming counts at 59 roadside survey routes in this area. Simulation parameters are based on a previous analysis of RUGR point count data (Lewis et al. 2022) and a preliminary analysis of ARU data from the study region. Local abundance within 200 m of point count centers is simulated as a function of 5 landscape covariates: mean elevation, Lattitude, mean aspect (decomposed into additive effects of northness and eastness), and mean canopy height. The detection process is decomposed into two components: availability (the probability of a RUGR being within the survey area and drumming during a survey), and condtional detectability (the probability of detecting an available grouse). Point count and ARU surveys are simulated as occurring over a 12-week period from March 9 - May 31st. Availability is simulated based on a quadratic effect of day of year, peaking on April 19th. Both ARUs and point counts jointly estimate the linear and quadratic effects of day of year on availability but have separate intercepts to account for differences in the detection process between methods. Point count detections within 20 m distance bands are simulated based on a Half-Normal detection function. ARUs are simulated as recording for 3 hours each morning over the 12-week period, with the data containing a binary response of whether or not at least one drum is detected on a recording. This binary response is related to local abundance via the per-individual detection probability (ARU availability * ARU detectability) through the methods of Royle & Nichols (2003). Preliminary analyses suggest that the probability of detecting at least one drum on an ARU recordingshould be > 95% if an available grouse drums at least 15 times total during a 3-hour recording. ARUs are simulated as being randomly deployed at a series of point count survey locations. Data are simulated under four different levels of grouse densities, implemented by modifiying the intercept of the abundance submodel. Predicting across the entire North Georgia study area, these intercepts correspond to total male population abundance of 2643, 5056, 9669, and 18494, corresponding to average male grouse densities of 0.003, 0.006, 0.011, and 0.022 males/ha, respectively. Data are also simulated under different monitoring scenarios varying the number of point count survey routes (20, 30, 40, or 50), the number of replicate point count surveys (1 - 3), and the number of survey points at which ARUs are deployed (0, 40, 80).

# RUGR_analyzing_simdata_NGA_code.R

R and NIMBLE code for analyzing simulated RUGR survey data generated in RUGR_generating_simdata_NGA_code.R. Point count and ARU data are analyzed using Bayesian hierarchical models, with datasets jointly estimating the latent abundance state but each utilizing a separate detection process. Abundance is modeled as a function of mean elevation, Latitude, northness, eastness, mean canopy height, and an offset for survey area. The point count detection process is based on the hierarchical distance-sampling models of Royle et al. (2004), modified to account for repeat surveys, while the ARU detection process is modeled based on the methods of Royle & Nichols (2003). For both datasets, the detection process is decomposed into availability and conditional detectability.

# RUGR.NGA.simulation.data.gzip

Parameters for generating simulated monitoring data for RUGR are stored in the 'RUGR.NGA.simulation.data.gzip' file.
### surveydata
A matrix giving the landscape covariates for each of 633 potential survey locations in North Georgia. Suvey locations (Point) are organized into 59 roadside survey routes (Route). Survey-level landscape covariates are mean elevation (elev), Latitude, aspect (decomposed to eastness and northness), and canopy height (canheight) within 200 m of the survey location centers. All landscape covariates are standardized.
### N_GA_covariates
A matrix giving the standardized landscape covariates (elev, Latitude, eastness, northness, and canheight) for each of the 57070 390-390m grid cells across the North Georgia study area. This matrix is used to predict RUGR abundance across the study area via simulation posterior samples in RUGR_analyzing_simdata_NGA_code.R. The area of each grid cell (ha) is included as an intercept to adjust between the area of surveys (12.6 ha) and grid cells (15.21 ha).
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
The probability of detecting at least one drum from an available RUGR during a 3-hour recording.
### ARU_dates
A vector giving the day of year (standardized) of ARU recordings.
### PC_dates
A 12 x 7 matrix giving the day of year (standardized) for each day of the 12-week season. Rows correspond to weeks of the season (1-12) while columns correspond to days of the week (1-7). This matrix is used to simulate the day of year for point counts.
### PC_periods
A 12 x 3 matrix giving the weeks (rows) in which to perform surveys when using either 1 (1st column), 2 (2nd column), or 3 (3rd column) repeat point counts. With only 1 repeat point count survey, survey dates were roughly evenly divided across the entire 12-week period to ensure adequate temporal coverage. With 2 repeat point count surveys, survey dates for the first survey were roughly evenly divided across the first 6 weeks and the dates for the second survey were roughly evenly divided across the second 6 weeks. For 3 repeat point count surveys, survey dates for the first survey were roughly evenly divided across the first 4 weeks, the dates for the second survey were roughly evenly divided across weeks 5 - 8, and the dates for the third survey were roughly evenly divided across weeks 9 - 12.
### surveyarea
The area in hectares surveyed on each point count (200 m fixed-radius surveys).
### PC_det_bands
A vector giving the boundaries of 20 m distance bands (out to 200 m) in which detected RUGR are binned into during point count surveys. 
### abund_N_GA_actual
The actual abundance predicted across the study area under the four different average RUGR densities (0.003 - 0.022 males/ha).
