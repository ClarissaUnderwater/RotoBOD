#######################################################################################################
###### RotoBOD data processing script ###### 

# this script reads any .log files in the current directory, estimates oxygen concentrations and calculates respiration
# the output is roto.dt, containing all the data and roto.dt_unique, containing the respiration rates per sample
# samples are identified by filename, RotoBOD wheel and RotoBOD position, format your metadata accordingly


###### important output variables within roto.dt ###### 

# mc_pairs number of data point pairs that were used in the Monte-Carlo simulation - if this is too low, the standard error will
# be high, consider decreasing timedist to allow more pairings, or record moredata

# mc_slope oxygen drawdown in your sample in µmol O2 L-1 day-1
# mc_SE standard error of mc_slope

# resp_mc oxygen respiration per organism or particle in your bottle (standardized to bottle volume) in µmol O2 individual-1 day-1
# resp_SE standard error of resp_mc
#######################################################################################################


#######################################################################################################
### adapt these parameters to your data set
#######################################################################################################
setwd("/my_directory/") # the directory where your data files are stored
salinity <- 35 # the salinity of your samples in psu, average for marine samples would be 35
# salinity affects absolute oxygen concentrations, won't affect respiration slopes
timedist <- 0.7 # minimum time in hours that the data points for the Monte-Carlo-simulation should be apart from each other
# this is to avoid including random slopes between pulse replicate data points
# necessary distance depends on number of wheels in use
bottle_volume <- 10 # volume in mL


#######################################################################################################
# necessary packages to run this script
#######################################################################################################
library(data.table)
library(dplyr)
library(broom)

#######################################################################################################
# reading the files into a data table
#######################################################################################################
# Initialize data.table
roto.dt <- data.table()

# Read data from text files
files <- dir(pattern = ".log") # lists the log files (text files ending in .log)
roto.data <- vector("list", length = length(files))
names(roto.data) <- gsub(".log", "", files, fixed = TRUE)

for (i in 1:length(files)) {
  print(files[i])
  roto.data[[i]] <- read.table(files[i], header = FALSE, sep = "", quote = "", skip = 3, stringsAsFactors = FALSE)
  colnames(roto.data[[i]]) <- c(
    "amplitude", "phase", "oxygen", "error", "position", "bottle", "measurement",
    "Date", "Time", "IRdetT", "IRBotT", "locating_spot", "wheel"
  )
}

# Combine into one data.table
roto.dt <- rbindlist(roto.data, idcol = "ID")
rm(roto.data)
roto.dt <- as.data.table(roto.dt)

# create unique identifier for each sample based on file of origin and position in the RotoBOD
roto.dt[, sample := paste(ID, wheel, bottle)]


#######################################################################################################
# calculating accurate oxygen concentrations
#######################################################################################################

# remove values with low amplitudes, these can occur if the sensor is not centered on the spot
roto.dt <- subset(roto.dt, roto.dt$amplitude > 10000)

# Calculate new_phase
roto.dt[, new_phase := phase / 100]

# Temperature Correction based on temperature detected on the spot
roto.dt[, airsat := calc_air_sat(new_phase, IRBotT)]
roto.dt[, o2conc := calc_o2_conc(airsat, IRBotT)]
roto.dt[, o2conc_umol_L := o2conc * 31.25]
roto.dt[, datetime := as.POSIXct(lubridate::mdy_hms(paste(Date, Time, sep = " ")))]############
roto.dt[, unclassdatetime := as.POSIXct(datetime)]

# Salinity Correction based on sample properties
roto.dt[, salinity := ifelse(bottle %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"), 33, NA)]
roto.dt[, scaledT := log((298.15 - IRBotT) / (273.15 + IRBotT))]
roto.dt[, oxy_corrected_salinity := calc_o2_conc_salinity(salinity, scaledT, o2conc_umol_L)]


#######################################################################################################
# calculating the duration of the run in hours
#######################################################################################################

# Convert time into seconds
roto.dt[, timenum := as.numeric(hms(Time))]

# Estimate day
roto.dt[, day := case_when(
  Date %in% c("01/01/23") ~ 1,
  Date %in% c("01/02/23") ~ 2,
  Date %in% c("01/03/23") ~ 3,
  Date %in% c("01/04/23") ~ 4,
  Date %in% c("01/05/23") ~ 5
)]

# Calculate timenum and sample
roto.dt[, timenum := timenum + (day - 1) * 24 * 60 * 60]

# Calculate hours
roto.dt[, hours := timenum / 60 / 60]


#######################################################################################################
# calculating linear slopes 
#######################################################################################################

# Estimate slopes by bottle and perform statistical tests
roto.dt[, Slopes := coef(lm(o2conc_umol_L ~ hours, data = .SD))[[2]], by = sample]
roto.dt[, Shapiro := shapiro.test(o2conc_umol_L)$p.value, by = sample]
roto.dt[, testchoice := ifelse(Shapiro > 0.05, "Pearson", "Spearman"), by = sample]
roto.dt[, coefficient := {
  if (testchoice[1] == "Pearson") {
    cor.test(hours, o2conc_umol_L, method = "pearson")$estimate
  } else {
    cor.test(hours, o2conc_umol_L, method = "spearman")$estimate
  }
}, by = sample]


#######################################################################################################
#### function for the Monte Carlo simulation
#######################################################################################################

Monte_Carlo_Slope_Sim <- function(subset_data, sample_id){ 
  
  num_rows <- length(subset_data$ID)  # store the number of observations within the time range
  
  empty_slopes <- data.frame(Iteration = c(1:(num_rows%/%2)),  #                    Then set up an empty data frame with half the rows (with no remainder)
                             Row_1 = as.numeric(rep(NA, (num_rows%/%2))),  #       Set up an empty dataframe as detailed in the intro
                             Row_2 = as.numeric(rep(NA, (num_rows%/%2))),
                             Time_1 = as.numeric(rep(NA, (num_rows%/%2))),
                             Time_2 = as.numeric(rep(NA, (num_rows%/%2))),
                             O2_1 = as.numeric(rep(NA, (num_rows%/%2))),
                             O2_2 = as.numeric(rep(NA, (num_rows%/%2))),
                             dT = as.numeric(rep(NA, (num_rows%/%2))),
                             dO2 = as.numeric(rep(NA, (num_rows%/%2))),
                             dO2_dT_per_hour = as.numeric(rep(NA, (num_rows%/%2))))

  pair_no <- 1 # increase number of pairs by multiplying 
  samples <- replicate(pair_no * num_rows%/%2, sample(c(1:num_rows), 2, replace = FALSE)) # generate the rows in the bottle data to sample
 
  
  slopes <- empty_slopes %>% mutate(Row_1 = samples[1, Iteration], #              Now we're basically just pulling all of our info from the Bottle df
                                    Row_2 = samples[2, Iteration], #              And populating the new df with calcs, then we can filter it so
                                    Time_1 = subset_data$hours[Row_1], #         the minimum dT is an hour
                                    Time_2 = subset_data$hours[Row_2],
                                    O2_1 = subset_data$o2conc_umol_L[Row_1],
                                    O2_2 = subset_data$o2conc_umol_L[Row_2],
                                    dT = Time_1 - Time_2,
                                    dO2 = O2_1 - O2_2,
                                    dO2_dT_per_hour = (dO2/dT))
  
  slopes_with_ok_range <- slopes %>% filter(abs(dT) > timedist) # only include slopes where the time difference was greater than an hour (in seconds)
  
  slope_statistics <- data.frame(Number_of_Pairs = length(slopes_with_ok_range$Iteration),
                                 Avg_dO2_umol_hour = mean(slopes_with_ok_range$dO2_dT_per_hour),
                                 Slope_SE = sd(slopes_with_ok_range$dO2_dT_per_hour)/sqrt(length(slopes_with_ok_range$Iteration)))
  

  #############
  # Create a histogram for the dO2_dT_per_hour values
  n <- n+1
  histogram <- ggplot(data=slopes_with_ok_range, aes(x = dO2_dT_per_hour)) +
    geom_histogram(binwidth = 0.1, fill = "blue", color = "black") +
    labs(title = paste("Histogram for Sample", sample_id[n], length(slopes_with_ok_range$Iteration)),
         x = "dO2_dT_per_hour",
         y = "Frequency")

  # Print the histogram to a PDF
  print(histogram)
  #############
  
  return(slope_statistics)
}

#######################################################################################################
#### make histogram file and run simulation for each sample in roto.dt
#######################################################################################################

pdf("histograms.pdf", width = 8, height = 6)
n <- 0

results_list <- lapply(unique_samples, function(sample_id) {
  subset_data <- roto.dt %>% filter(sample == sample_id)
  subset_result <- Monte_Carlo_Slope_Sim(subset_data, sample_id)
  return(subset_result)
})

dev.off() 

##### write results into roto.dt

# Naming results_list items by unique_samples names
names(results_list) <- unique_samples

# Initialize an empty data frame for slope_table
slope_table <- data.frame(sample = character(),
                          mc_slope = numeric(),
                          mc_pairs = numeric(),
                          mc_SE = numeric(),
                          stringsAsFactors = FALSE)

# Loop through unique samples
for (sample_name in unique_samples) {
  # Find the corresponding result from 'results_list'
  sample_result <- results_list[[sample_name]]
  
  # Extract the three items from the sample result
  mc_pairs <- sample_result$Number_of_Pairs
  mc_slope <- sample_result$Avg_dO2_umol_hour * 24
  mc_SE <- sample_result$Slope_SE * 24
  
  # Create a data frame for the current sample
  sample_table <- data.frame(sample = sample_name,
                             mc_slope = mc_slope,
                             mc_pairs = mc_pairs,
                             mc_SE = mc_SE)
  
  # Append the sample_table to the slope_table
  slope_table <- rbind(slope_table, sample_table)
}

##### merge 
roto.dt <-merge(roto.dt,slope_table, by.x="sample", by.y="sample") # merge two data frames by overlapping variables
roto.dt$resp_mc <- - roto.dt$mc_slope /1000 * bottle_volume  # calculate respiration per organism
roto.dt$resp_SE <- - roto.dt$mc_SE /1000 * bottle_volume  # calculate respiration per organism
roto.dt_unique <- roto.dt[!duplicated(roto.dt$sample), ]

#######################################################################################################
#### test if the simulation was successful by comparing to linear estimates (slopes should be close to identical)
#######################################################################################################
# test for linear correlation
reg_model <- lm(Slopes ~ mc_slope, data = roto.dt_unique)

# Get the regression coefficients
m_sample_b <- coef(reg_model)[2]  # Slope
m_sample_a <- coef(reg_model)[1]  # Intercept

# Get the R-squared value
slopes_vs_mc_r2 <- summary(reg_model)$r.squared


#######################################################################################################
#### plot oxygen over time
#######################################################################################################

ggplot(data = roto.dt) +
  geom_point(mapping = aes(y=o2conc_umol_L, x=hours), shape = 1, size = 0.00001, color = "darkred") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border= element_rect(fill="transparent",color="black", size = 1) , legend.position = "none") +
  labs(x = "time (h)", y = "Oxygen (µmol L-1)") +
  facet_wrap(facets = vars(sample), ncol = 6) 

  ggsave("oxygen_raw.pdf", width=6, height=9, dpi=600)

#######################################################################################################
#### plot respiration per sample 
#######################################################################################################

ggplot(data = roto.dt_unique) +
  geom_point(mapping = aes(y = resp_mc, x = hours), shape = 21, size = 2) +
  geom_errorbar(aes(y = resp_mc, x = hours, ymin = resp_mc - resp_SE, ymax = resp_mc + resp_SE), size=0.2, width = 0.2,)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        panel.border = element_rect(fill = "transparent", color = "black", size = 1), legend.position = "none") +
  labs(x = "sample", y = "respiration (µmol ind-1 day-1)") +
  facet_wrap(facets = vars(wheel), ncol = 1) 
  
  ggsave("respiration.pdf", width = 8, height = 4, dpi = 300)

