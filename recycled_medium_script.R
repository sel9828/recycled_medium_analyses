# Project: Microalgae Recycled Growth Medium Experiments
# Author: Sarah Loftus, sarah.e.loftus@gmail.com

# This code organizes, analyses, and visualizes data from mircoalgae recycled medium studies. 
# Results of analyses are published in Loftus & Johnson, 2019, "   ", Journal. 

# Data is available at Figshare: DOI

#######################################################################################

### Install packages ###

  # Data management
    library(tidyverse) # includes ggplot2
    library(lubridate)
    library(broom)
  
  # Analyses
    library(nlme)
    library(car)
  
  # Plotting
    library(gtable)
    library(grid)

### Read Data (see metadata for CSV file column descriptions)
  growth_df <- read.csv("input_data/DataTable1_Growth.csv")        # Algae growth and experimental data for each culture replicate
  daily_df <- read.csv("input_data/DataTable2_Daily.csv")          # Daily sampling data 
  filtrate_df <- read.csv("input_data/DataTable3_Filtrate.csv")    # Post-experiment saved filtrate data
    
#### Parse dates and times, and calculate Days elapsed ####

    daily_df %>%
      mutate(Date = mdy(Date),   
             Tstart = hm(Tstart),
             Tend = hm(Tend),
             Tmid = Tstart + as.period(0.5*as.duration(Tend - Tstart))) %>%   # Tmid = midpoint between start and end of sampling time
      group_by(Algae) %>% 
      do(mutate(., DaysElapsed = ((Date + Tmid) - (Date[Round == 0 & Day == 0] + Tmid[Round == 0 & Day == 0]) )/ 86400) ) %>%   # DaysElapsed =  days elapsed since start of each experiment
      select(Algae, Round, Day, Chl_medium, OD750_medium, DaysElapsed) -> daily_df2 
      # Note: warning messages do not pose problem to further analysis
  
  # Find average and standard deviation of OD of ASW growth medium across all experiments
    avg_ASW_OD <- mean(daily_df$OD750_medium)
    sd_ASW_OD <- sd(daily_df$OD750_medium)

#### Define and calculate additional variables ####

  C_mol_wt <- 12.0107 # g/mol
  
  growth_df %>%    
    inner_join(daily_df2, by = c("Algae", "Round", "Day") ) %>%   # Add daily measurement data to growth data frame
    mutate(uM_PC = PC/(Vol/1000),                      # uM_PC = Blank-corrected PC concentration of biomass in micromolar
           uM_PN = PN/(Vol/1000),                      # uM PN = Blank-corrected PN concentration of biomass in micromolar
           CN_ratio = PC/PN,                           # CN_raio = C/N ratio of algae biomass (unitless)
           C_per_cell = (C_mol_wt*uM_PC/(AlgaeConc*10^9))*10^6,   # C_per_cell = Carbon content of cells in picogram C/cell
           AlgaeLipids = BulkLipids - ExtLipids,       # AlgaeLipids = Lipid concentration of algae biomass in relative fluorescence units (RFU)
           AlgaeLipidsPerCell = AlgaeLipids/AlgaeConc, # AlgaeLipidsPerCell = Lipid content of cells in RFU/(10^6 cells/mL)
           DOC_rate = TOC_rate - POC_rate,             # DOC_rate = DOC release rate in uM/day. Note: For Chlorella sp. D046 experiment, this is not the calculated rate but is the raw, blank-corrected, mean DPM (disintegrations per minute) divided by 14C DPM added. Therefore absolute values cannot be compared among replicates or treatments (however, comparing fractions of POC/TOC, etc. is still allowable).
           fraction_DOC_rate = DOC_rate/TOC_rate,      # fraction_DOC_rate = DOC release rate as a fraction of TOC production rate. Note: For Chlorella sp. D046 experiment, this is the fraction of the raw, blank-corrected, mean DPM (disintegrations per minute).
           biomass_OD = OD750 -  OD750_medium,         # biomass_OD = Blank-corrected OD750 of the culture, in arbitrary units
           biomass_chl = Chl - Chl_medium              # biomass_chl = Blank-corrected chlorophyll-a concentration of the culture, in relative fluorescence units.
           ) -> growth_df2    

  # Calculate daily averages and standard deviations of all variables in growth data frame. For use in plotting. ###

    growth_df2 %>%
      select(-Replicate) %>% 
      group_by(Algae, Round, Treatment, Day) %>% 
      summarize_all(funs(mean, sd), na.rm = TRUE) %>%
      replace(., is.na(.), NA)-> growth_df_avgs  
    
    # Find % difference in average bacteria concentrations between fresh and recycled treatments
    
    growth_df_avgs %>%  
      select(Algae, Round, Treatment, Day, BacteriaConc_mean) %>% 
      filter(Round != 0) %>% # Remove Round 0, when all treatments were in Fresh Medium
      spread(Treatment, BacteriaConc_mean) %>%  # Converts Fresh and Recycled to individual columns, with the values in the columns being the mean bacteria concentration
      mutate(percent_more_bacteria = 100*(R-F)/F) -> bacteria_percent_diff  # calcualtes the % increase in bacteria in Recycled treatment compared to Fresh control

        # Find average percent_more_bacteria on Day 5, removing C323 Round 3
        avg_percent_more_bacteria <- mean(bacteria_percent_diff$percent_more_bacteria[bacteria_percent_diff$Algae != "C323" & bacteria_percent_diff$Round != 3 & bacteria_percent_diff$Day == 5])
        sd_percent_more_bacteria <- sd(bacteria_percent_diff$percent_more_bacteria[bacteria_percent_diff$Algae != "C323" & bacteria_percent_diff$Round != 3 & bacteria_percent_diff$Day == 5])
        
#### Fill in missing data for multivariate analyses of growth parameters ####
  
  # Create regression equations of missing data with algae cell concentration.   
  # Missing PC data in Round 4 of "F" treatment in Navicula experiment for Replicates A, B, C. Replace with estimated PC based on algae cell concentration linear regression.
      PC_vs_cell_Navicula <- lm(uM_PC ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "Navicula", ]) 
      # R^2 = 0.967, p-val = 2.2e-16; uM_PC = 9113.0*AlgaeConc - 218.8
  
  # Missing lipid data. Replace with estimated lipids based on algae cell concentration linear regression with bulk lipid data from Day 5 only.
    
    # Missing lipid data in Round 2, Day 5, "R" Treatment, Replicate F in Navicula dataset. 
      lipid_vs_cell_Navicula <- lm(AlgaeLipids ~ AlgaeConc, data = growth_df2[growth_df2$Algae == "Navicula" & growth_df2$Day == 5, ]) # only use the Day 5 lipids for regression 
      # R^2 = 0.404, p-val = 0.00066; AlgaeLipids = 165039*AlgaeConc - 78695
    
    # Missing lipid data in Round 1, Day 5, "F" Treatment, Replicate A of C323 data.
      lipid_vs_cell_C323 <- lm(AlgaeLipids ~ AlgaeConc, data = growth_df2[growth_df2$Algae == 'C323' & growth_df2$Day == 5, ]) # only use the Day 5 lipids for regression 
      # R^2 = 0.8502, p-val = 1.427e-07; AlgaeLipids = 23395*AlgaeConc + 984
  
  # Use regression to fill in missing PC & lipid values for Navicula  (missing PC values: Treatment F, Day 5, Round 4, all Replicates; missing lipid value: Treatment R, Day 5, Round 2, Replicate F)
    growth_df2$uM_PC[growth_df2$Treatment == "F" & growth_df2$Round == 4 & growth_df2$Day == 5 & growth_df2$Algae == "Navicula"] <- 
          PC_vs_cell_Navicula$coefficients[1] + growth_df2$AlgaeConc[growth_df2$Treatment == "F" & growth_df2$Round == 4 & growth_df2$Day == 5 & growth_df2$Algae == "Navicula"]*PC_vs_cell_Navicula$coefficients[2]
    
    growth_df2$AlgaeLipids[growth_df2$Treatment == "R" & growth_df2$Round == 2 & growth_df2$Day == 5 & growth_df2$Replicate == "F" & growth_df2$Algae == "Navicula"] <- 
         lipid_vs_cell_Navicula$coefficients[1] + growth_df2$AlgaeConc[growth_df2$Treatment == "R" & growth_df2$Round == 2 & growth_df2$Day == 5 & growth_df2$Replicate == "F" & growth_df2$Algae == "Navicula"]*lipid_vs_cell_Navicula$coefficients[2]
   
     # Re-calculate ALgaeLipidsPerCell since using in further analyses. 
      growth_df2$AlgaeLipidsPerCell[growth_df2$Treatment == "R" & growth_df2$Round == 2 & growth_df2$Day == 5 & growth_df2$Replicate == "F" & growth_df2$Algae == "Navicula"] <-
           growth_df2$AlgaeLipids[growth_df2$Treatment == "R" & growth_df2$Round == 2 & growth_df2$Day == 5 & growth_df2$Replicate == "F" & growth_df2$Algae == "Navicula"] / growth_df2$AlgaeConc[growth_df2$Treatment == "R" & growth_df2$Round == 2 & growth_df2$Day == 5 & growth_df2$Replicate == "F" & growth_df2$Algae == "Navicula"]
      
  # Use regression to fill in missing lipid value for C323 (missing lipid value: Round 1, Day 5, Treatment F, Replicate A)
    growth_df2$AlgaeLipids[growth_df2$Algae == 'C323' & growth_df2$Treatment == "F" & growth_df2$Round == 1 & growth_df2$Day == 5 & growth_df2$Replicate == "A"] <- 
          lipid_vs_cell_C323$coefficients[1] + growth_df2$AlgaeConc[growth_df2$Algae == 'C323' & growth_df2$Treatment == "F" & growth_df2$Round == 1 & growth_df2$Day == 5 & growth_df2$Replicate == "A"]*lipid_vs_cell_C323$coefficients[2]
  
    # Re-calculate ALgaeLipidsPerCell since using in further analyses:
    growth_df2$AlgaeLipidsPerCell[growth_df2$Algae == 'C323' & growth_df2$Treatment == "F" & growth_df2$Round == 1 & growth_df2$Day == 5 & growth_df2$Replicate == "A"] <- 
           growth_df2$AlgaeLipids[growth_df2$Algae == 'C323' & growth_df2$Treatment == "F" & growth_df2$Round == 1 & growth_df2$Day == 5 & growth_df2$Replicate == "A"] / growth_df2$AlgaeConc[growth_df2$Algae == 'C323' & growth_df2$Treatment == "F" & growth_df2$Round == 1 & growth_df2$Day == 5 & growth_df2$Replicate == "A"]

#### Calculate variables for algae growth statistics ####

    day_start <- 0          # day_start = Starting day of each round
    day_end <- 5            # day_end = Last day of each round
    mu_period_D046 <- 0:2   # mu_period_D046 = Days when Chlorella sp. D046 was in exponential pahse, to use for calculating specific growth rate
    mu_period_diatom <- 0:3 # mu_period_diatom = Days when C323 and Navicula sp. were in exponential phase, to use for calculating specific growth rate
  
  # Specific growth rate 
    growth_df2 %>%
      filter( (Algae == "D046" & Day %in% mu_period_D046) | 
                (Algae %in% c("C323","Navicula") & Day %in% mu_period_diatom) ) %>%    # Based on the algae, select only the days used to calculate specific growth rate
      group_by(Algae, Round, Treatment, Replicate) %>% 
      do(tidy(lm(log(AlgaeConc) ~ DaysElapsed, data = .))) %>%   # linear regression of natural log algae cell concentration versus time
      filter(term == "DaysElapsed") %>%                          # where "DaysElapsed" is the row term in the linear regression output to indicate the slope, i.e., the specific growth rate
      rename(mu = estimate) %>%                                  # rename the slope of the linear regression ("estimate") as mu, mu = the specific growth rate (units = 1/day)
      select(Algae, Round, Treatment, Replicate, mu) -> mu_data        
  
  # Additional variables to summarize: Final algae lipid content, Change in algae carbon, and Final Fv/Fm
    growth_df2 %>%
      group_by(Algae, Round, Treatment, Replicate) %>% 
      summarize(AlgaeLipidsPerCell_final = AlgaeLipidsPerCell[Day == day_end],        # AlgaeLipidsPerCell_final = Final (Day 5) algae lipid content 
                PC_net = (uM_PC[Day == day_end] - uM_PC[Day == day_start])/1000,      # PC_net = Change in algae carbon concentration over the round, in mM
                FvFm_final = FvFm[Day == day_end]) %>%                                # FvFm_final = Final (Day 5) Fv/Fm 
      inner_join(mu_data, by = c("Algae", "Round", "Treatment", "Replicate") ) -> MultivarGrowthData  #join the mu data
  
  # Create summary table of averages and standard deviations 
    MultivarGrowthData %>%
      select(-Replicate) %>% # Remove 'Replicate' column
      filter(Round != 0) %>% # Remove Round 0, when all treatments were in Fresh Medium
      group_by(Algae, Round, Treatment) %>%
      summarize_all(funs(mean, sd), na.rm = TRUE) -> MultivarGrowthData_Summary # Take average and standard deviation of the different Reuse + Treatment combinations, save the table

  
#### Compare cumulative and net DOC produced, and DOC consumed by bacteria ####
  
  # Estimate production rates for missing days (Day 1 & 3) based on rates measured on following and previous days
  
  # Make function to calculate missing data
    missing_fun <- function(x){
      x <- ifelse( is.na(x) == TRUE, (lead(x, 1) + lag(x, 1))/2, x) # take average of value before and after the missing value
      return(x)
    }
  
  # Use function to replace missing values & calculate additional DOC fractions
  growth_df2 %>%
    filter(Round != 0 & Algae != "D046") %>%   # do not include Round 0; do not include Chlorella sp. D046 because there is no absolute rate data available
    select(Algae, Round, Treatment, Replicate, Day, DaysElapsed, TOC_rate, POC_rate, DOC_rate) %>%
    group_by(Algae, Round, Treatment, Replicate) %>%  
    do(mutate_at(., vars(ends_with("rate")), funs(missing_fun)) ) %>% # Apply the function to replace NA values
    mutate(fraction_DOC_rate = DOC_rate/TOC_rate) -> carbon_rate_df  # Re-calculate fraction_DOC_rate for all days, now that missing data has been replaced
  
  ### Use rates to estimate cumulative production over a round.
  
  # Light hours per photoperiod (with photoperiod being 12 hours) (carbon production rates are already in units of days, with 12 hrs light/day as defined by the experiment protocol)
  light_time <- 12/12      # light_time = light hours per hours of photoperiod (unitless). We assume carbon production is negligible in the dark, so rates are only applicable during the photoperiod
  day0_light_time <- 7/12  # day0_light_time = hours of light (12:00 to 19:00) the cultures experienced on Day 0 after sampling
  day5_light_time <- 5/12  # day5_light_time = hours of light (7:00 to 12:00) the cultures experienced on Day 5 until the end of sampling
  light_per_round <- c(day0_light_time, rep(light_time, 4), day5_light_time) # light_per_round = vector of light hours experienced by the culture from Days 0 through 5, respectively
  
  carbon_rate_df %>%
    select(-fraction_DOC_rate) %>% 
    group_by(Algae, Round, Treatment, Replicate) %>%
    do(summarize_at(., vars(ends_with("rate")), funs(sum(.*light_per_round/1000)))) %>%         # Multiply the carbon production rates by the light_per_round vector. Divide by 1000 to convert to mM
    rename(TOC_cumulative = TOC_rate, POC_cumulative = POC_rate, DOC_cumulative = DOC_rate) %>% # Rename the multiplication products (which still have "_rate" suffixes, as "_cumulative" suffixes to indicate cumulative carbon produced over the round)
    mutate(fraction_DOC_cumulative = DOC_cumulative/TOC_cumulative) -> carbon_cumulative_df     # fraction_DOC_cumulative = cumulative DOC released over the round as a fraction of cumulative carbon produced over the round
  
  # Calculate change in measured (net) PC (GF/F filter, 0.8 µm) and DOC (0.2 µm) over growth period and add to carbon_cumulative_df, then calculate means and st devs
  
  growth_df2 %>% 
    filter(Round != 0) %>%   # do not include Round 0
    group_by(Algae, Round, Treatment, Replicate) %>% 
    summarize(PC_net = ( (uM_PC[Day == day_end] - uM_PC[Day == day_start])/1000 ),    # Same as PC_net in MultivarGrowthData; convert to mM
              DOC_net = (DOC[Day == day_end] - DOC[Day == day_start])/1000 ) %>%      # DOC_net = change in DOC concentration over the round; convert to mM
    mutate(fraction_DOC_net = DOC_net/(PC_net + DOC_net)) %>%                         # fraction_DOC_net = DOC_net as a fraction of change in total carbon over the round (with total carbon calcualted as the sum of net_PC and net_DOC)                         
    full_join(carbon_cumulative_df, by = c("Algae", "Round", "Treatment", "Replicate")) -> carbon_cum_net_reps_df   # data with individual replicates still intact
    
  carbon_cum_net_reps_df %>%   
    filter(!is.na(DOC_net)) %>%      # For Chlorella sp. D046, RB R3 D5 measured DOC is missing, so remove NAs before calculating means 
    mutate(DOC_percent_lost = 100*(DOC_cumulative - DOC_net)/DOC_cumulative) %>%   # DOC_percent_lost = percent of released (cumulative) DOC missing per round after accounting for DOC measured (net) in the growth medium
    group_by(Algae, Round, Treatment) %>% 
    select(-Replicate) %>% # Remove 'Replicate' column
    summarize_all(funs(mean, sd)) -> carbon_cumulative_net_df  # Means and standard deviations summary table for both estimated (cumulative) and measured (net) carbon production
    
  ### Estiamte carbon consumed by bacteria and add it to the carbon_cumulative_net_df data frame
  
  # Define bacteria growth efficiency and bacteria carbon content 
  bge <- 0.15 # bge = bacteria growth efficiency, unitless; fraction of carbon consumed that becomes biomass
  C_content <- 1*10^-13 # C_content = bacteria carbon content in grams C/bacteria cell
   
  growth_df2 %>%
    filter(Round != 0) %>%   # do not include Round 0
    select(Algae, Round, Treatment, Replicate, Day, BacteriaConc) %>%
    group_by(Algae, Round, Treatment, Replicate) %>%
    summarize(BacteriaConc_change = (BacteriaConc[Day == day_end] - BacteriaConc[Day == day_start]) *10^9 ) %>%  # BacteriaConc_change = change in bacteria concentration from Day 0 to Day 5. Convert from 10^6 cells/mL to cells/Liter.
    mutate(bact_C_consumed = (BacteriaConc_change*C_content/bge)* 1000/C_mol_wt) -> bact_C_consumed_df # bact_C_consumed = carbon consumed by bacteria over a round. Multiply bacteria concentration change by C content, and divide by growth efficiency. Convert to mM C.
  
  # Join the bacteria_C_consumed data frame to the measured (net) & estimated (cumulative) carbon data frame
    
  bact_C_consumed_df %>%    
    select(-BacteriaConc_change) %>%              
    full_join(carbon_cum_net_reps_df, by = c("Algae", "Round", "Treatment", "Replicate")) %>% 
    mutate(fraction_C_bact_consumed = bact_C_consumed/(DOC_cumulative - DOC_net)) %>%  # fraction_C_bact_consumed = fraction of the difference between cumulative and net DOC that the bacteria C consumption accounts for
    select(-Replicate) %>% 
    summarize_all(funs(mean, sd)) -> carbon_all_df  # Means and standard deviations for replicates


#### Comparison of normalized biomass yield vs initial DOC in the recycled medium treatments ####

  growth_df2 %>% 
    group_by(Algae, Round, Treatment, Replicate) %>% 
    summarize(initial_DOC = DOC[Day == day_start]) %>%  # initial_DOC = biologically-derived DOC concentration on Day 0 of each round, units: µM
    inner_join(MultivarGrowthData, by = c("Algae", "Round", "Treatment", "Replicate")) %>% 
    select(-mu, -FvFm_final, -AlgaeLipidsPerCell_final) -> initialDOC_PC_df
   
  # Find the mean change in biomass for fresh medium (control) replicates within each round of the experiment
    initialDOC_PC_df %>% 
      group_by(Algae, Round, Treatment) %>% 
      filter(Treatment == "F") %>% 
      summarize(mean_Fresh_PCnet = mean(PC_net)) -> initialDOC_PC_Fresh_df  #mean_Fresh_PCnet = average PC_net of Fresh medium replicates within a round
  
    # Add the mean PC_net of Fresh replicates as a column in the original initialDOC_PC_df data frame, 
      # and calculate normalized biomass by dividing each recycled medium replicate PC_net 
      # by mean_Fresh_PCnet. Only retain rows of the recycled medium treatments. 
     
      initialDOC_PC_df %>% 
        filter(Treatment == "R") %>% 
        inner_join(initialDOC_PC_Fresh_df, by = c("Algae", "Round")) %>% 
        ungroup() %>% 
        select(-Treatment.x, -Treatment.y) %>% 
        mutate(PC_R_normalized = PC_net/mean_Fresh_PCnet)-> initialDOC_PC_final  #PC_R_normalized = PC_net of Recycled medium treatments divided by the mean PC_net of Fresh medium control replicates from the same round
      
        # Check out trends in data
        plot(PC_net ~ initial_DOC, data = initialDOC_PC_final)
        plot(PC_R_normalized ~ initial_DOC, data = initialDOC_PC_final)
        abline(lm(initialDOC_PC_final$PC_R_normalized ~ initialDOC_PC_final$initial_DOC))

        # Remove outliers (the 2nd and 3rd Round for Algae C323), which heavily influence the trendline
          initialDOC_PC_final%>% 
            filter(! (Algae == "C323" & Round > 1)) -> initialDOC_PC_no_outliers
        
          # Check out resulting data with outliers removed
          plot(PC_R_normalized ~ initial_DOC, data = initialDOC_PC_no_outliers)
          abline(lm(initialDOC_PC_no_outliers$PC_R_normalized ~ initialDOC_PC_no_outliers$initial_DOC))

#### Summarize data from post-experiment saved filtrate ####

  filtrate_df %>% 
    group_by(Algae, ElapsedDays, Treatment) %>% 
    select(-Replicate, -Date) %>% 
    summarize_all(funs(mean, sd)) -> filtrate_df_avgs
  
  # Calculate C consumed by bacteria in the filtrate
  
    filtrate_df %>% 
      select(-Date, -DOC, -TDN) %>% 
      group_by(Algae, Treatment, Replicate) %>% 
      summarize(BacteriaConc_change = (BacteriaConc[ElapsedDays == max(ElapsedDays)] - BacteriaConc[ElapsedDays == min(ElapsedDays)]) *10^9 ) %>%  # BacteriaConc_change = change in bacteria concentration over the experiment. Convert from 10^6 cells/mL to cells/Liter.
      mutate(bact_C_consumed = (BacteriaConc_change*C_content/bge)* 1000/C_mol_wt) %>%   # bact_C_consumed = carbon consumed by bacteria. Multiply bacteria concentration change by C content, and divide by growth efficiency. Convert to mM C.
      select(-Replicate) %>% 
      group_by(Algae, Treatment) %>% 
      summarize_all(funs(mean, sd)) -> filtrate_C_consumed_df


##### Statistics #####

#### 1. Do growth-related variables differ in fresh versus recycled medium across multiple medium reuses? ####

    # # Use MANOVA to determine differences between treatments with 4 dependent variables across rounds 1-4
    # 
    # # Check for colinearity of response variables
    # pairs(MultivarGrowthData[,5:8])
    # cor(MultivarGrowthData[MultivarGrowthData$Round != 0,5:8]) # Pearson's 
    # 
    # # MANOVA with error term for Replicate to account for dependency of recycled medium replicate bottles across rounds
    # 
    # manova.mod.D046 <- manova(cbind(AlgaeLipidsPerCell_final, PC_net, FvFm_final, mu) ~ Treatment*Round + Error(Replicate), 
    #                           data = MultivarGrowthData[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0, ]) 
    # 
    # manova.mod.C323 <- manova(cbind(AlgaeLipidsPerCell_final, PC_net, FvFm_final, mu) ~ Treatment*Round + Error(Replicate), 
    #                           data = MultivarGrowthData[MultivarGrowthData$Algae == "C323" & MultivarGrowthData$Round != 0, ]) 
    # 
    # manova.mod.Navi <- manova(cbind(AlgaeLipidsPerCell_final, PC_net, FvFm_final, mu) ~ Treatment*Round + Error(Replicate), 
    #                           data = MultivarGrowthData[MultivarGrowthData$Algae == "Navicula" & MultivarGrowthData$Round != 0, ]) 
    # 
    # # Check significance of model terms
    # summary(manova.mod.D046)
    # summary(manova.mod.C323)
    # summary(manova.mod.Navi)

  # UPDATED METHOD - ANALYZE EACH VARIABLE SEPARATELY WITH NLME MODEL AS DONE IN QUESTION 2 #
    
    # C323
    # 1. Biomass yield
    PC_net.lme.C323 <- lme(PC_net ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "C323" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(PC_net.lme.C323)
    
    # 2. Specific growth rate
    mu.lme.C323 <- lme(mu ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "C323" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(mu.lme.C323)
    
    # 3. Final lipid content
    lipids.lme.C323 <- lme(AlgaeLipidsPerCell_final ~ Treatment*Round,
                       random = ~1 | Replicate,
                       correlation = corAR1(form = ~ Round | Replicate), 
                       data = MultivarGrowthData[MultivarGrowthData$Algae == "C323" & MultivarGrowthData$Round != 0, ],
                       method = "REML")
    summary(lipids.lme.C323)
    
    # 4. Final Fv/Fm
    FvFm.lme.C323 <- lme(FvFm_final ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "C323" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(FvFm.lme.C323)
    
    # Navi
    # 1. Biomass yield
    PC_net.lme.Navi <- lme(PC_net ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "Navicula" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(PC_net.lme.Navi)
    
    # 2. Specific growth rate
    mu.lme.Navi <- lme(mu ~ Treatment*Round,
                       random = ~1 | Replicate,
                       correlation = corAR1(form = ~ Round | Replicate), 
                       data = MultivarGrowthData[MultivarGrowthData$Algae == "Navicula" & MultivarGrowthData$Round != 0, ],
                       method = "REML")
    summary(mu.lme.Navi)
    
    # 3. Final lipid content
    lipids.lme.Navi <- lme(AlgaeLipidsPerCell_final ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "Navicula" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(lipids.lme.Navi)
    
    # 4. Final Fv/Fm
    FvFm.lme.Navi <- lme(FvFm_final ~ Treatment*Round,
                         random = ~1 | Replicate,
                         correlation = corAR1(form = ~ Round | Replicate), 
                         data = MultivarGrowthData[MultivarGrowthData$Algae == "Navicula" & MultivarGrowthData$Round != 0, ],
                         method = "REML")
    summary(FvFm.lme.Navi)
    
    # D046
    # 1. Biomass yield
    PC_net.lme.D046 <- lme(PC_net ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(PC_net.lme.D046)
    
    # 2. Specific growth rate
    mu.lme.D046 <- lme(mu ~ Treatment*Round,
                       random = ~1 | Replicate,
                       correlation = corAR1(form = ~ Round | Replicate), 
                       data = MultivarGrowthData[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0, ],
                       method = "REML")
    summary(mu.lme.D046)
    
    # 3. Final lipid content
    lipids.lme.D046 <- lme(AlgaeLipidsPerCell_final ~ Treatment*Round,
                           random = ~1 | Replicate,
                           correlation = corAR1(form = ~ Round | Replicate), 
                           data = MultivarGrowthData[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0, ],
                           method = "REML")
    summary(lipids.lme.D046)
    
    # 4. Final Fv/Fm
    FvFm.lme.D046 <- lme(FvFm_final ~ Treatment*Round,
                         random = ~1 | Replicate,
                         correlation = corAR1(form = ~ Round | Replicate), 
                         data = MultivarGrowthData[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0, ],
                         method = "REML")
    summary(FvFm.lme.D046)
    
    
    # Check significance of model terms
    PC_net.lme.C323.Anova <- Anova(PC_net.lme.C323)
    PC_net.lme.Navi.Anova <- Anova(PC_net.lme.Navi)
    PC_net.lme.D046.Anova <- Anova(PC_net.lme.D046)
    mu.lme.C323.Anova <- Anova(mu.lme.C323)
    mu.lme.Navi.Anova <- Anova(mu.lme.Navi)
    mu.lme.D046.Anova <- Anova(mu.lme.D046)
    lipids.lme.C323.Anova <- Anova(lipids.lme.C323)
    lipids.lme.Navi.Anova <- Anova(lipids.lme.Navi)
    lipids.lme.D046.Anova <- Anova(lipids.lme.D046)
    FvFm.lme.C323.Anova <- Anova(FvFm.lme.C323)
    FvFm.lme.Navi.Anova <- Anova(FvFm.lme.Navi)
    FvFm.lme.D046.Anova <- Anova(FvFm.lme.D046)
    
    # interaction plot for selected variables/algae
    interaction.plot(MultivarGrowthData$Round[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0], 
                     MultivarGrowthData$Treatment[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0],
                     MultivarGrowthData$FvFm_final[MultivarGrowthData$Algae == "D046" & MultivarGrowthData$Round != 0])
    

#### 2.	Is DOC production rate in recycled medium different from that in fresh medium? Does DOC production rate change across multiple reuses of the medium? ####

    # Use the variable fraction_DOC_cumulative since this normalizes DOC release rate by TOC production rate over the round
    # Keep Round as an integer variable using lme with Replicate bottles as a random factor
      
      # C323
      DOC.lme.C323 <- lme(fraction_DOC_cumulative ~ Treatment*Round, 
                            random = ~1 | Replicate, 
                            correlation = corAR1(form = ~ Round | Replicate), 
                            data = carbon_cumulative_df[carbon_cumulative_df$Algae == "C323", ], 
                            method = "REML")
      summary(DOC.lme.C323)
        
      # Navicula
      DOC.lme.Navicula <- lme(fraction_DOC_cumulative ~ Treatment*Round, 
                            random = ~1 | Replicate, 
                            correlation = corAR1(form = ~ Round | Replicate), 
                            data = carbon_cumulative_df[carbon_cumulative_df$Algae == "Navicula", ], 
                            method = "REML")
      summary(DOC.lme.Navicula)  
      
      
      # Check significance of model terms
      DOC.lme.C323.Anova <- Anova(DOC.lme.C323)
      DOC.lme.Navi.Anova <- Anova(DOC.lme.Navicula)
      
      # interaction plots
      interaction.plot(carbon_cumulative_df$Round[carbon_cumulative_df$Algae == "Navicula"], 
                       carbon_cumulative_df$Treatment[carbon_cumulative_df$Algae == "Navicula"],
                       carbon_cumulative_df$fraction_DOC_cumulative[carbon_cumulative_df$Algae == "Navicula"])
      
      interaction.plot(carbon_cumulative_df$Round[carbon_cumulative_df$Algae == "C323"], 
                       carbon_cumulative_df$Treatment[carbon_cumulative_df$Algae == "C323"],
                       carbon_cumulative_df$fraction_DOC_cumulative[carbon_cumulative_df$Algae == "C323"])

#### 3. Relationship between normalized biomass change and initial DOC in the recycled medium ####

    # Model with lme including autocorrelation, with Replicates as random effects
    ## Without "outliers" (C323 Rounds 2 and 3) removed
    PC.lme.outliers <- lme(PC_R_normalized ~ initial_DOC, 
                              random = ~1 | Algae/Replicate, 
                              correlation = corAR1(form = ~ Round | Algae/Replicate), 
                              data = initialDOC_PC_final, 
                              method = "REML")
    
      summary(PC.lme.outliers)
      
      # Check significance of model
      PC.lme.outliers.Anova <- Anova(PC.lme.outliers)
      
      # Check residuals of model
      PC.lme.resid.outliers <- residuals(PC.lme.outliers)
      hist(PC.lme.resid.outliers) 
      plot(fitted(PC.lme.outliers), PC.lme.resid.outliers)
    
    ## With "outliers" (C323 Rounds 2 and 3) removed
    PC.lme <- lme(PC_R_normalized ~ initial_DOC, 
                     random = ~1 | Algae/Replicate, 
                     correlation = corAR1(form = ~ Round | Algae/Replicate), 
                     data = initialDOC_PC_no_outliers, 
                     method = "REML")
    
      summary(PC.lme)
      
      # Check significance of model
      PC.lme.Anova <- Anova(PC.lme)
      
      # Check residuals of model
      PC.lme.resid <- residuals(PC.lme)
      hist(PC.lme.resid) 
      plot(fitted(PC.lme), PC.lme.resid)

#### Figures ####

#### Figure 1: Daily Algae Concentrations ####

  # C323 #### 
    C323_algaeConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$AlgaeConc_mean), ], 
                                  aes(x = as.numeric(DaysElapsed_mean), y = AlgaeConc_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 2.5) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=AlgaeConc_mean - AlgaeConc_sd, ymax=AlgaeConc_mean + AlgaeConc_sd), width= 0.4) +
      labs(x = "Days", y = bquote('Algae'~(10^6~'cells/mL'))) +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(limits = c(0, 4.5), breaks = seq(0, 4.5, 2), expand = expand_scale(mult = c(0.05,0.1))) +  
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 32, 5)) +
      annotate("text", x = 2.4, y = 4.4, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             legend.position = c(0.93,0.23),
             axis.title.y = element_blank(),  # moves axis title away from axis label
             axis.title.x = element_text(margin = margin(r = 15), size = 14),   
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA), 
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12)) 
            # (note: missing several error bars in R0 because of missing data)
    
  # D046 ####    
    D046_algaeConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046", ], 
                                  aes(x = as.numeric(DaysElapsed_mean), y = AlgaeConc_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 2.5) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=AlgaeConc_mean - AlgaeConc_sd, ymax=AlgaeConc_mean + AlgaeConc_sd), width= 0.4) +
      labs(y = " ") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(breaks = seq(0, 10, 2), expand = expand_scale(mult = c(0.05,0.15))) +  
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 30, 5)) +
      annotate("text", x = 2.3, y = 9, label = expression(paste(bold("A"), italic("  Chlorella"), " sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 9, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),
             axis.title.x = element_blank(),   
             legend.position = c(0.94,0.23),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA), 
             panel.background = element_rect(fill = "white"),
             axis.ticks = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 12)) 
    
  # Navicula ####  
    Navicula_algaeConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula", ], 
                                      aes(x = as.numeric(DaysElapsed_mean), y = AlgaeConc_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 2.5) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=AlgaeConc_mean - AlgaeConc_sd, ymax=AlgaeConc_mean + AlgaeConc_sd), width= 0.4) +
      labs(x = "Days", y = bquote('Algae'~(10^6~'cells/mL'))) +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(breaks = seq(0, 1.2, 0.4), expand = expand_scale(mult = c(0.05,0.12))) +  # Change scale depending on values
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 32, 5)) +
      annotate("text", x = 1.5, y = 1.05, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),   
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),            
             panel.background = element_rect(fill = "white"),
             axis.ticks = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 12)) 
    
  # Combine Daily Algae Concentration plots in a grid & Save as pdf ####
  
    # Convert ggplots to grobs
    D046_algaeConc_grob <- ggplotGrob(D046_algaeConc_plot)
    C323_algaeConc_grob <- ggplotGrob(C323_algaeConc_plot)
    Navi_algaeConc_grob <- ggplotGrob(Navicula_algaeConc_plot)
    
    # Combine in grid
    algaeConc_gridplot <- rbind(D046_algaeConc_grob, Navi_algaeConc_grob, C323_algaeConc_grob, size = "first")
    
    grid.newpage()
    grid.draw(algaeConc_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/Fig1_AlgaeConcGrid.pdf", 
        width = 9.5, height = 5.5)
    grid.draw(algaeConc_gridplot)
    dev.off()
    
#### Figure S1: Daily OD ####

  # D046 ####    
    D046_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046", ], 
                           aes(x = as.numeric(DaysElapsed_mean), y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.4) +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(breaks = seq(0, 0.15, 0.05), expand = expand_scale(mult = c(0.05,0.2))) +  
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 32, 5)) +
      labs( y = " ") +
      annotate("text", x = 2.5, y = 0.17, label = expression(paste(bold("A"), italic("  Chlorella"), " sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 0.17, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),
             axis.title.x = element_blank(),
             legend.position = c(0.94,0.23),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA), 
             panel.background = element_rect(fill = "white"),
             axis.ticks = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 12))   
  
  # Navicula #### 
     Navicula_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula", ], 
                                aes(x = as.numeric(DaysElapsed_mean), y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.4) +
      labs(y = expression("OD"["750"])) +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(breaks = seq(0, 0.16, 0.05), expand = expand_scale(mult = c(0.05, 0.12))) +  
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 32, 5)) +
      annotate("text", x = 1.5, y = 0.16, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             legend.position = c(0.94,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  
             axis.title.x = element_blank(),   
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),            
             panel.background = element_rect(fill = "white"),
             axis.ticks = element_blank(),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 12)) 
   
  # C323 ####  
    C323_OD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323", ], 
                           aes(x = as.numeric(DaysElapsed_mean), y = biomass_OD_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=biomass_OD_mean - biomass_OD_sd, ymax=biomass_OD_mean + biomass_OD_sd), width= 0.4) +
      labs(x = "Days") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
      scale_y_continuous(breaks = seq(0, 0.2, 0.05), expand = expand_scale(mult = c(0.05,0.2))) +  
      scale_x_continuous(limits = c(-0.2,30), breaks = seq(0, 32, 5)) +
      annotate("text", x = 2.4, y = 0.2, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9),
             legend.position = c(0.94,0.23),
             axis.title.y = element_blank(),  
             axis.title.x = element_text(margin = margin(r = 35), size = 14),   
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA), 
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12))   
    
  # Combine Daily OD plots in a grid & Save as pdf ####
    
    # Convert ggplots to grobs
    D046_OD_grob <- ggplotGrob(D046_OD_plot)
    C323_OD_grob <- ggplotGrob(C323_OD_plot)
    Navi_OD_grob <- ggplotGrob(Navicula_OD_plot)
    
    # Combine in grid
    OD_gridplot <- rbind(D046_OD_grob, Navi_OD_grob, C323_OD_grob, size = "first")
    
    grid.newpage()
    grid.draw(OD_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/FigS1_ODGrid.pdf", 
        width = 9.5, height = 5.5)
    grid.draw(OD_gridplot)
    dev.off()

#### Figure S2: Daily Bacteria Concentrations ####  

  # C323 ####
  C323_bacteriaConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323", ], 
                               aes(x = as.numeric(DaysElapsed_mean), y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 0.4) +
    labs(x = "Days") +
    scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
    scale_y_continuous(breaks = seq(0, 25, 5), expand = expand_scale(mult = c(0.1,0.2))) +  
    annotate("text", x = 2.4, y = 23, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.text = element_text(size = 10),
           legend.position = c(0.07,0.55),
           legend.title = element_blank(),
           axis.title.y = element_blank(),  
           axis.title.x = element_text(margin = margin(r = 15), size = 14),
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),
           axis.ticks = element_blank(),
           axis.text = element_text(size = 12)) 
  
  # D046 ####  
  D046_bacteriaConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046", ], 
                                   aes(x = as.numeric(DaysElapsed_mean), y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 0.4) +
    labs(y = " ") +
    scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
    scale_y_continuous(breaks = seq(0, 20, 5), expand = expand_scale(mult = c(0.1,0.2))) +  # Change scale depending on values
    annotate("text", x = 2.5, y = 20, label = expression(paste(bold("A"), italic("  Chlorella"), " sp. D046")), size = 5) + 
    annotate("text", x = c(9, 15, 21, 27), y = 20, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.text = element_text(size = 10),
           legend.position = c(0.07,0.55),
           legend.title = element_blank(),
           axis.title.y = element_text(margin = margin(r = 15), size = 14),
           axis.title.x = element_blank(),
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),
           axis.ticks = element_blank(),
           axis.text.y = element_text(size = 12),
           axis.text.x = element_blank()) 
  
  # Navicula ####  
  Navi_bacteriaConc_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula", ], 
                                   aes(x = as.numeric(DaysElapsed_mean), y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 0.4) +
    labs(y = bquote('Bacteria'~(10^6~'cells/mL'))) +
    scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +
    scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
    scale_y_continuous(breaks = seq(0, 55, 10), expand = expand_scale(mult = c(0.05,0.15))) +  # Change scale depending on values
    annotate("text", x = 1.5, y = 52, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.text = element_text(size = 10),
           legend.position = c(0.07,0.55),
           legend.title = element_blank(),
           axis.title.y = element_text(margin = margin(r = 15), size = 14),  
           axis.title.x = element_blank(),
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),
           axis.ticks = element_blank(),
           axis.text.y = element_text(size = 12),
           axis.text.x = element_blank()) 
  
  # Combine bacteria Conc plots in a grid and Save as pdf ####
  
    # Convert ggplots to grobs
    D046_bacteriaConc_grob <- ggplotGrob(D046_bacteriaConc_plot)
    C323_bacteriaConc_grob <- ggplotGrob(C323_bacteriaConc_plot)
    Navi_bacteriaConc_grob <- ggplotGrob(Navi_bacteriaConc_plot)
    
    # Combine in grid
    bacteriaConc_gridplot <- rbind(D046_bacteriaConc_grob, Navi_bacteriaConc_grob, C323_bacteriaConc_grob, size = "first")
    
    grid.newpage()
    grid.draw(bacteriaConc_gridplot)
  
    # Save plot as pdf file
    pdf("Figures/FigS2_bacteriaConcGrid.pdf", 
        width = 9.5, height = 5.5)
    grid.draw(bacteriaConc_gridplot)
    dev.off()
  
#### Figure 2: Multivariate Algae Growth Parameters ####
  
  # C323 ####
    # PC net ####
    PC_plot_C323 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "C323", ], 
                           aes(x = Round, y = PC_net_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=PC_net_mean - PC_net_sd, ymax=PC_net_mean + PC_net_sd), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Algae Biomass (mM C)", title = expression(paste(italic("S. sourniae"), " C323"))) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0,8, by = 2)) +
      annotate("text", x = 0.75, y = 8, label = expression(paste(bold("C"))), size = 4) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9, margin = margin(r = 8)),
             legend.position = "top",
             legend.key.size = unit(3, "mm"),
             legend.justification = "center",
             legend.margin=margin(t = 0, b = 0, l = 10),
             legend.box.margin=margin(t = 0, b = 0),
             legend.spacing.x = unit(1, "mm"),
             plot.title = element_text(hjust = 0.5, size = 10),
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.text.x = element_blank()) 
    
    # Lipids ####
    lipids_plot_C323 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "C323", ],  
                               aes(x = Round, y = AlgaeLipidsPerCell_final_mean/10^6, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(AlgaeLipidsPerCell_final_mean - AlgaeLipidsPerCell_final_sd)/10^6, ymax=(AlgaeLipidsPerCell_final_mean + AlgaeLipidsPerCell_final_sd)/10^6), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Lipid Content (RFU/cells/mL) ") +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.15, by = 0.05)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43')) +   
      annotate("text", x = 0.75, y = 0.12, label = expression(paste(bold("I"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) 
    
    # mu ####
    mu_plot_C323 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "C323", ], 
                           aes(x = Round, y = mu_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(mu_mean - mu_sd), ymax=(mu_mean + mu_sd)), width= 0.4, position = position_dodge(0.9)) +
      scale_y_continuous(expand = expand_scale(mult = c(0.05, 0.1)), breaks = seq(-0.5, 1.5, by = 0.5)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43')) +   
      annotate("text", x = 0.75, y = 1.4, label = expression(paste(bold("F"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) 
  
    # FvFm final ####
    FvFm_plot_C323 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "C323", ], 
                             aes(x = Round, y = FvFm_final_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(FvFm_final_mean - FvFm_final_sd), ymax=(FvFm_final_mean + FvFm_final_sd)), width= 0.4, position = position_dodge(0.9)) +
      labs(y = 'Final Fv/Fm', x = "Medium Reuses") +
      scale_fill_manual(labels = c("Fresh", "Recycled"), 
                        values = c('gray63','gray43')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.6, by = 0.2)) +
      annotate("text", x = 0.75, y = 0.6, label = expression(paste(bold("L"))), size = 4) + 
      theme( axis.title.x = element_blank(),
             legend.position = "none",                              
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.text = element_text(size = 9),
             axis.title.y = element_blank()) 
    
  # D046 ####
    # PC net ####
    PC_plot_D046 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "D046", ], 
                           aes(x = Round, y = PC_net_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=PC_net_mean - PC_net_sd, ymax=PC_net_mean + PC_net_sd), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Algae Biomass (mM C)", title = expression(paste(italic("Chlorella"), " sp. D046"))) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0,8, by = 2)) +
      annotate("text", x = 0.75, y = 6.3, label = expression(paste(bold("A"))), size = 4) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9, margin = margin(r = 8)),
             legend.position = "top",
             legend.key.size = unit(3, "mm"),
             legend.justification = "center",
             legend.margin=margin(t = 0, b = 0, l = 10),
             legend.box.margin=margin(t = 0, b = 0),
             legend.spacing.x = unit(1, "mm"),
             plot.title = element_text(hjust = 0.5, size = 10),
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(margin = margin(r = 8),size = 10),
             axis.title.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.text.x = element_blank()) 
    
    # Lipids ####
    lipids_plot_D046 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "D046", ],  
                               aes(x = Round, y = AlgaeLipidsPerCell_final_mean/10^6, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(AlgaeLipidsPerCell_final_mean - AlgaeLipidsPerCell_final_sd)/10^6, ymax=(AlgaeLipidsPerCell_final_mean + AlgaeLipidsPerCell_final_sd)/10^6), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Lipid Content (RFU/cells/mL) ") +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.01, by = 0.004)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4')) +   
      annotate("text", x = 0.75, y = 0.0095, label = expression(paste(bold("G"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(margin = margin(r = 8),size = 10),
             axis.title.x = element_blank()) 
    
    # mu ####
    mu_plot_D046 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "D046", ], 
                           aes(x = Round, y = mu_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(mu_mean - mu_sd), ymax=(mu_mean + mu_sd)), width= 0.4, position = position_dodge(0.9)) +
      labs(y = expression(paste("Growth rate ", ("day"^-1)))) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 1.5, by = 0.5)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4')) +   
      annotate("text", x = 0.75, y = 1.7, label = expression(paste(bold("D"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_text(margin = margin(r = 8),size = 10),
             axis.title.x = element_blank()) 
    
    
    # FvFm final ####
    FvFm_plot_D046 <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "D046", ], 
                             aes(x = Round, y = FvFm_final_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(FvFm_final_mean - FvFm_final_sd), ymax=(FvFm_final_mean + FvFm_final_sd)), width= 0.4, position = position_dodge(0.9)) +
      labs(y = expression(paste("Final ", "F"["v"]*"/F"["m"])), x = " ") +  # Placeholder x-axis title to keep space for grid plot
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.6, by = 0.2)) +
      annotate("text", x = 0.75, y = 0.55, label = expression(paste(bold("J"))), size = 4) + 
      theme( axis.title.x = element_text(margin = margin(r = 25), size = 9),
             legend.position = "none",                              
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.text = element_text(size = 9),
             axis.title.y = element_text(margin = margin(r = 8), size = 10)) 
    
    
  # Navicula #####
    # PC net ####
    PC_plot_Navicula <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "Navicula", ], 
                               aes(x = Round, y = PC_net_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=PC_net_mean - PC_net_sd, ymax=PC_net_mean + PC_net_sd), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Algae Biomass (mM C)", title = expression(paste(italic("Navicula"), " sp."))) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.05)), breaks = seq(0,10, by = 2)) +
      annotate("text", x = 0.75, y = 10, label = expression(paste(bold("B"))), size = 4) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.title = element_blank(),
             legend.text = element_text(size = 9, margin = margin(r = 8)),
             legend.position = "top",
             legend.key.size = unit(3, "mm"),
             legend.justification = "center",
             legend.margin=margin(t = 0, b = 0, l = 10),
             legend.box.margin=margin(t = 0, b = 0),
             legend.spacing.x = unit(1, "mm"),
             plot.title = element_text(hjust = 0.5, size = 10),
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.title = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.text.x = element_blank()) 
    
    # Lipids ####
    lipids_plot_Navicula <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "Navicula", ],  
                                   aes(x = Round, y = AlgaeLipidsPerCell_final_mean/10^6, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(AlgaeLipidsPerCell_final_mean - AlgaeLipidsPerCell_final_sd)/10^6, ymax=(AlgaeLipidsPerCell_final_mean + AlgaeLipidsPerCell_final_sd)/10^6), width= 0.4, position = position_dodge(0.9)) +
      labs(y = "Lipid Content (RFU/cells/mL) ") +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.15, by = 0.05)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4')) +   
      annotate("text", x = 0.75, y = 0.13, label = expression(paste(bold("H"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) 
    
    # mu ####
    mu_plot_Navicula <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "Navicula", ], 
                               aes(x = Round, y = mu_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(mu_mean - mu_sd), ymax=(mu_mean + mu_sd)), width= 0.4, position = position_dodge(0.9)) +
      labs(y = 'Specific growth rate (1/day)') +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 1.5, by = 0.5)) +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4')) +   
      annotate("text", x = 0.75, y = 1.4, label = expression(paste(bold("E"))), size = 4) + 
      theme( legend.position = "none",
             panel.background = element_rect(fill = NA),
             #panel.border = element_rect(color = "gray60", fil = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.text.x = element_blank(),
             axis.text.y = element_text(size = 9),
             axis.ticks.x = element_blank(),
             axis.title.y = element_blank(),
             axis.title.x = element_blank()) 
    
    
    # FvFm final ####
    FvFm_plot_Navicula <- ggplot(data = MultivarGrowthData_Summary[MultivarGrowthData_Summary$Algae == "Navicula", ], 
                                 aes(x = Round, y = FvFm_final_mean, fill = Treatment))  +
      geom_col(position = position_dodge(0.9), color = "gray15") +
      geom_errorbar(aes(ymin=(FvFm_final_mean - FvFm_final_sd), ymax=(FvFm_final_mean + FvFm_final_sd)), width= 0.4, position = position_dodge(0.9)) +
      labs(y = 'Final Fv/Fm', x = "Medium Reuses") +
      scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4')) +
      scale_y_continuous(expand = expand_scale(mult = c(0, 0.1)), breaks = seq(0, 0.6, by = 0.2)) +
      annotate("text", x = 0.75, y = 0.6, label = expression(paste(bold("K"))), size = 4) + 
      theme( axis.title.x = element_text(margin = margin(t = 5), size = 11),  
             legend.position = "none",                              
             panel.background = element_rect(fill = NA),
             axis.line = element_line(color = "black", size = 0.4),
             axis.ticks.x = element_blank(),
             axis.text = element_text(size = 9),
             axis.title.y = element_blank()) 

  # Arrange multivariate plots in a grid and Save as pdf ####
    
    # Convert ggplots to grobs
    g1 <- ggplotGrob(PC_plot_D046)
    g2 <- ggplotGrob(PC_plot_C323)
    g3 <- ggplotGrob(PC_plot_Navicula)
    g4 <- ggplotGrob(mu_plot_D046)
    g5 <- ggplotGrob(mu_plot_C323)
    g6 <- ggplotGrob(mu_plot_Navicula)
    g7 <- ggplotGrob(lipids_plot_D046)
    g8 <- ggplotGrob(lipids_plot_C323)
    g9 <- ggplotGrob(lipids_plot_Navicula)
    g10 <- ggplotGrob(FvFm_plot_D046)
    g11 <- ggplotGrob(FvFm_plot_C323)
    g12 <- ggplotGrob(FvFm_plot_Navicula)
    
    # Bind the plots by column
    col1 <- rbind(g1, g4, g7, g10, size = "first")
    col3 <- rbind(g2, g5, g8, g11, size = "first")
    col2 <- rbind(g3, g6, g9, g12, size = "first")
    # Fix the column widths to the max plot width
    col1$widths <- unit.pmax(g1$widths, g4$widths, g7$widths, g10$widths)
    col2$widths <- unit.pmax(g2$widths, g5$widths, g8$widths, g11$widths)
    col3$widths <- unit.pmax(g3$widths, g6$widths, g9$widths, g12$widths)
    # Combine all plots
    g_final <- cbind(col1, col2, col3, size = "first")
    
    grid.newpage()
    grid.draw(g_final)
    
    # Save plot as pdf file
    pdf("Figures/Fig2_Multivar.pdf", 
        width = 6, height = 7.5)
    grid.draw(g_final)
    dev.off()
  
   
#### Figure S3: Daily pH of cultures ####

  # C323 ####    
    C323_pH_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323", ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = pH_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=pH_mean - pH_sd, ymax=pH_mean + pH_sd), width= 0.4) +
      labs(x = "Days") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(6,10),breaks = seq(6,10, 1), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 2.4, y = 9.7, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_blank(),
             axis.title.x = element_text(margin = margin(r = 15), size = 14),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12)) 
  
    # D046 ####  
    D046_pH_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046", ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = pH_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=pH_mean - pH_sd, ymax=pH_mean + pH_sd), width= 0.4) +
      labs(y = " ") +  # x axis label as space holder for space
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(6, 10),breaks = seq(0, 10, 1), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 2.5, y = 9.7, label = expression(paste(bold("A"), italic("  Chlorella"), " sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 9.7, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
  
    # Navicula ####  
    Navi_pH_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula", ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = pH_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=pH_mean - pH_sd, ymax=pH_mean + pH_sd), width= 0.4) +
      labs(y = "pH") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(6,10),breaks = seq(6, 10, 1), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 1.5, y = 9.7, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.grid.minor = element_line(colour = "gray95", size = 0.1),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    
  # Combine pH plots in a grid and Save as pdf ####
  
    # Convert ggplots to grobs
    D046_pH_grob <- ggplotGrob(D046_pH_plot)
    C323_pH_grob <- ggplotGrob(C323_pH_plot)
    Navi_pH_grob <- ggplotGrob(Navi_pH_plot)
    
    # Combine in grid
    pH_gridplot <- rbind(D046_pH_grob, Navi_pH_grob, C323_pH_grob, size = "first")
    
    grid.newpage()
    grid.draw(pH_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/FigS3_pHGrid.pdf", 
        width = 9, height = 5.5)
    grid.draw(pH_gridplot)
    dev.off()  

#### Figure S4: DIC over time ####

  # C323 ####    
    C323_DIC_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$DIC_mean), ],   
                            aes(x = as.numeric(DaysElapsed_mean), y = DIC_mean/1000, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=(DIC_mean - DIC_sd)/1000, ymax=(DIC_mean + DIC_sd)/1000), width= 0.4) +
      labs(x = "Days") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(0,12),breaks = seq(0,12, 4), expand = expand_scale(mult = c(0,0))) +
      annotate("text", x = 2.4, y = 9.7, label = expression(paste(bold("B"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_blank(),
             axis.title.x = element_text(margin = margin(r = 15), size = 14),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12)) 
  # Navicula ####  
    Navi_DIC_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$DIC_mean),],   
                            aes(x = as.numeric(DaysElapsed_mean), y = DIC_mean/1000, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=(DIC_mean - DIC_sd)/1000, ymax=(DIC_mean + DIC_sd)/1000), width= 0.4) +
      labs(x = "Days", y = "DIC (mM)") +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=DIC_mean - DIC_sd, ymax=DIC_mean + DIC_sd), width= 0.4) +
      labs(y = "DIC (mM)") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(0,5),breaks = seq(0, 5, 2), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 1.5, y = 4.5, label = expression(paste( bold("A"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.8),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    
  # Combine DIC plots in a grid ####
    
    # Convert ggplots to grobs
    C323_DIC_grob <- ggplotGrob(C323_DIC_plot)
    Navi_DIC_grob <- ggplotGrob(Navi_DIC_plot)
    
    # Combine in grid
    DIC_gridplot <- rbind(Navi_DIC_grob, C323_DIC_grob, size = "first")
    
    grid.newpage()
    grid.draw(DIC_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/FigS4_DICGrid.pdf", 
        width = 9, height = 4)
    grid.draw(DIC_gridplot)
    dev.off() 
  
#### Figure S5: Salinity ####
    
    # C323 ####    
    C323_salt_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$Salinity_mean), ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = Salinity_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=Salinity_mean - Salinity_sd, ymax=Salinity_mean + Salinity_sd), width= 0.4) +
      labs(x = "Days") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(25,39),breaks = seq(25,39, 5), expand = expand_scale(mult = c(0,0))) +
      annotate("text", x = 2.4, y = 37, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_blank(),
             axis.title.x = element_text(margin = margin(r = 15), size = 14),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12)) 
    # D046 ####  
    D046_salt_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !is.na(growth_df_avgs$Salinity_mean), ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = Salinity_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=Salinity_mean - Salinity_sd, ymax=Salinity_mean + Salinity_sd), width= 0.4) +
      labs(y = " ") +  # x axis label as space holder for space
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(25, 39),breaks = seq(25, 39, 5), expand = expand_scale(mult = c(0,0))) +
      annotate("text", x = 2.5, y = 37, label = expression(paste(bold("A"), italic("  Chlorella "), "sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 37, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    # Navicula ####  
    Navi_salt_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$Salinity_mean), ],   
                           aes(x = as.numeric(DaysElapsed_mean), y = Salinity_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=Salinity_mean - Salinity_sd, ymax=Salinity_mean + Salinity_sd), width= 0.4) +
      labs(y = "Salinity (ppt)") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(25,39),breaks = seq(25, 39, 5), expand = expand_scale(mult = c(0,0))) +
      annotate("text", x = 1.5, y = 37, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.grid.minor = element_line(colour = "gray95", size = 0.1),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    
    # Combine salinity plots in a grid and Save as pdf ####
    
    # Convert ggplots to grobs
    D046_salt_grob <- ggplotGrob(D046_salt_plot)
    C323_salt_grob <- ggplotGrob(C323_salt_plot)
    Navi_salt_grob <- ggplotGrob(Navi_salt_plot)
    
    # Combine in grid
    salt_gridplot <- rbind(D046_salt_grob, Navi_salt_grob, C323_salt_grob, size = "first")
    
    grid.newpage()
    grid.draw(salt_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/FigS5_saltGrid.pdf", 
        width = 9, height = 5.5)
    grid.draw(salt_gridplot)
    dev.off()  
    
#### Figure S6: OD filtrate ####  
    
    # C323 ####    
    C323_filtrateOD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$OD750_filt_mean), ],   
                             aes(x = as.numeric(DaysElapsed_mean), y = OD750_filt_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=OD750_filt_mean - OD750_filt_sd, ymax=OD750_filt_mean + OD750_filt_sd), width= 0.4) +
      labs(x = "Days") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(0,0.06),breaks = seq(0, 0.06, 0.02), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 2.4, y = 0.056, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_blank(),
             axis.title.x = element_text(margin = margin(r = 15), size = 14),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text = element_text(size = 12)) 
    
    # D046 ####  
    D046_filtrateOD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !is.na(growth_df_avgs$OD750_filt_mean), ],   
                             aes(x = as.numeric(DaysElapsed_mean), y = OD750_filt_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=OD750_filt_mean - OD750_filt_sd, ymax=OD750_filt_mean + OD750_filt_sd), width= 0.4) +
      labs(y = " ") +  # x axis label as space holder for space
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(0,0.06),breaks = seq(0,0.06, 0.02), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 2.5, y = 0.055, label = expression(paste(bold("A"), italic("  Chlorella "), "sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 37, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
   
     # Navicula ####  
    Navi_filtrateOD_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$OD750_filt_mean), ],   
                             aes(x = as.numeric(DaysElapsed_mean), y = OD750_filt_mean, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=OD750_filt_mean - OD750_filt_sd, ymax=OD750_filt_mean + OD750_filt_sd), width= 0.4) +
      labs(y = expression("OD"["750"]*" of culture filtrate")) +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(limits = c(0,0.08),breaks = seq(0,0.08, 0.02), expand = expand_scale(mult = c(0,0.1))) +
      annotate("text", x = 1.4, y = 0.076, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.93,0.23),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.grid.minor = element_line(colour = "gray95", size = 0.1),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    
    # Combine OD filtrate plots in a grid and Save as pdf ####
    
    # Convert ggplots to grobs
    D046_filtrateOD_grob <- ggplotGrob(D046_filtrateOD_plot)
    C323_filtrateOD_grob <- ggplotGrob(C323_filtrateOD_plot)
    Navi_filtrateOD_grob <- ggplotGrob(Navi_filtrateOD_plot)
    
    # Combine in grid
    filtrateOD_gridplot <- rbind(D046_filtrateOD_grob, Navi_filtrateOD_grob, C323_filtrateOD_grob, size = "first")
    
    grid.newpage()
    grid.draw(filtrateOD_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/FigS6_filtrateODGrid.pdf", 
        width = 9, height = 5.5)
    grid.draw(filtrateOD_gridplot)
    dev.off()  

#### Figure 3: Measured/Accumulated/Net DOC concentration over time ####

  # C323 ####
  C323_DOCmeas_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$DOC_mean),],   # Remove NAs
                          aes(x = as.numeric(DaysElapsed_mean), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 0.4) +
    labs(x = "Days") +
    scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
    scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
    scale_y_continuous(limits = c(0, 0.305), breaks = seq(0,0.3, 0.1), expand = expand_scale(mult = c(0.05,0.2))) +
    annotate("text", x = 2.4, y = 0.3, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.text = element_text(size = 10),
           legend.title = element_blank(),
           legend.position = c(0.07,0.55),
           axis.title.y = element_blank(),
           axis.title.x = element_text(margin = margin(r = 15), size = 14),
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),
           axis.ticks = element_blank(),
           axis.text = element_text(size = 12)) 
  
    # D046 ####
    D046_DOCmeas_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !is.na(growth_df_avgs$DOC_mean),],   # Remove NAs
                                aes(x = as.numeric(DaysElapsed_mean), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
      geom_point(size = 3) +
      geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
      geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 0.4) +
      labs(y = " ") +
      scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
      scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
      scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
      scale_y_continuous(breaks = seq(0, 0.7, 0.2), expand = expand_scale(mult = c(0.05,0.2))) +
      annotate("text", x = 2.6, y = 0.68, label = expression(paste( bold("A"), italic(" Chlorella"), " sp. D046")), size = 5) + 
      annotate("text", x = c(9, 15, 21, 27), y = 0.683, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
      theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
             legend.text = element_text(size = 10),
             legend.title = element_blank(),
             legend.position = c(0.07,0.55),
             axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
             axis.title.x = element_blank(),
             panel.background = element_rect(fill = "white"),
             panel.grid.major = element_line(colour = "gray87", size = 0.2),
             panel.border = element_rect(color = "gray60", fil = NA),
             axis.ticks = element_blank(),
             axis.text.y = element_text(size = 12),
             axis.text.x = element_blank()) 
    
  # Navicula ####
  Navi_DOCmeas_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$DOC_mean),],   # Remove NAs
                              aes(x = as.numeric(DaysElapsed_mean), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
    geom_point(size = 3) +
    geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
    geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 0.4) +
    labs(x = "Days", y = "Biologically-derived DOC (mM)") +
    scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
    scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
    scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
    scale_y_continuous(breaks = seq(0, 2.2, 0.5), expand = expand_scale(mult = c(0.05,0.2))) +
    annotate("text", x = 1.4, y = 2.1, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
    theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
           legend.text = element_text(size = 10),
           legend.title = element_blank(),
           legend.position = c(0.07,0.55),
           axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
           axis.title.x = element_blank(),
           panel.background = element_rect(fill = "white"),
           panel.grid.major = element_line(colour = "gray87", size = 0.2),
           panel.border = element_rect(color = "gray60", fil = NA),
           axis.ticks = element_blank(),
           axis.text.y = element_text(size = 12),
           axis.text.x = element_blank()) 
  
  # Combine DOC Conc plots in a grid and Save as pdf ####
  
    # Convert ggplots to grobs
    D046_DOCmeas_grob <- ggplotGrob(D046_DOCmeas_plot)
    C323_DOCmeas_grob <- ggplotGrob(C323_DOCmeas_plot)
    Navi_DOCmeas_grob <- ggplotGrob(Navi_DOCmeas_plot)
    
    # Combine in grid
    DOCmeas_gridplot <- rbind(D046_DOCmeas_grob, Navi_DOCmeas_grob, C323_DOCmeas_grob, size = "first")
    
    grid.newpage()
    grid.draw(DOCmeas_gridplot)
    
    # Save plot as pdf file
    pdf("Figures/Fig3_DOCmeasGrid.pdf", 
        width = 9.5, height = 6.5)
    grid.draw(DOCmeas_gridplot)
    dev.off()

#### Figure 4: Released and Accumulated DOC per Round (Percent) #### (Note: "released" = "cumulative", "accumulated" = "net")
      
   # D046 ####
      # Accumulated DOC ####
      netDOCpercent_plot_D046 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "D046", ], 
                                 aes(x = Round, y = fraction_DOC_net_mean*100, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=(fraction_DOC_net_mean - fraction_DOC_net_sd)*100, ymax=(fraction_DOC_net_mean + fraction_DOC_net_sd)*100), width= 0.4, position = position_dodge(0.9)) +
        labs(x = "Medium Reuses", title = expression( paste( italic("Chlorella"), " sp. D046"))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0.05)), breaks = seq(0,5, by = 1)) +
        annotate("text", x = 0.75, y = 5.2, label = expression(paste(bold("C"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 4),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title.x = element_text(margin = margin(r = 22),size = 11),
               axis.title.y = element_blank(),
               axis.text = element_text(size = 11))   
      # Empty plot ####
      empty_plot <- ggplot(data.frame()) +
        labs(x = " ") +
        theme(panel.background = element_rect(fill = NA))
      
   # C323 ####
      # Accumulated DOC ####  
      netDOCpercent_plot_C323 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "C323", ], 
                                 aes(x = Round, y = fraction_DOC_net_mean*100, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=(fraction_DOC_net_mean - fraction_DOC_net_sd)*100, ymax=(fraction_DOC_net_mean + fraction_DOC_net_sd)*100), width= 0.4, position = position_dodge(0.9)) +
        labs(title = expression(paste(italic("S. sourniae"), " C323"))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('gray63','gray43')) +
        scale_y_continuous(expand = expand_scale(mult = c(0.05,0)), breaks = seq(0,30, by = 10)) +
        annotate("text", x = 0.7, y = 32, label = expression(paste(bold("B"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 1),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_blank(),
               axis.text = element_text(size = 10))
      # Released DOC ####  
      cumDOCpercent_plot_C323 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "C323", ], 
                                 aes(x = Round, y = fraction_DOC_cumulative_mean*100, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=(fraction_DOC_cumulative_mean - fraction_DOC_cumulative_sd)*100, ymax=(fraction_DOC_cumulative_mean + fraction_DOC_cumulative_sd)*100), width= 0.4, position = position_dodge(0.9)) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('gray63','gray43')) +
        labs(x = "Medium Reuses") +
        scale_y_continuous(expand = expand_scale(mult = c(0,0.05)), breaks = seq(0,40, by = 10)) +
        annotate("text", x = 0.7, y = 32, label = expression(paste(bold("E"))), size = 5) + 
        theme( legend.position = "none",
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_text(margin = margin(r = 22),size = 10),
               axis.text = element_text(size = 10)) 
      
   # Navicula ####  
      # Accumulated DOC ####
      netDOCpercent_plot_Navicula <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "Navicula", ], 
                                     aes(x = Round, y = fraction_DOC_net_mean*100, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=(fraction_DOC_net_mean - fraction_DOC_net_sd)*100, ymax=(fraction_DOC_net_mean + fraction_DOC_net_sd)*100), width= 0.4, position = position_dodge(0.9)) +
        labs(y = "Accumulated DOC\n(% of TOC)", x = " ", title = expression(paste(italic("Navicula"), " sp."))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('dodgerblue2', 'dodgerblue4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0)), breaks = seq(0,10, by = 2)) +
        annotate("text", x = 0.75, y = 8, label = expression(paste(bold("A"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 1),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_text(margin = margin(r = 75),size = 10),
               axis.text = element_text(size = 10))
      
    # Released DOC ####
      cumDOCpercent_plot_Navicula <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "Navicula", ], 
                                     aes(x = Round, y = fraction_DOC_cumulative_mean*100, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=(fraction_DOC_cumulative_mean - fraction_DOC_cumulative_sd)*100, ymax=(fraction_DOC_cumulative_mean + fraction_DOC_cumulative_sd)*100), width= 0.4, position = position_dodge(0.9)) +
        labs(x = "Medium Reuses", y = "Released DOC\n(% of TOC)") +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('dodgerblue2', 'dodgerblue4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0.05)), breaks = seq(0,80, by = 20)) +
        annotate("text", x = 0.75, y = 80, label = expression(paste(bold("D"))), size = 5) + 
        theme( legend.position = "none",
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_text(margin = margin(r = 15),size = 10),
               axis.text = element_text(size = 10))
   # Combine net and cumulative DOC plots in grid and save as pdf ####
      
      # Convert ggplots to grobs
      D046_netDOCpercent_grob <- ggplotGrob(netDOCpercent_plot_D046)
      C323_netDOCpercent_grob <- ggplotGrob(netDOCpercent_plot_C323)
      Navi_netDOCpercent_grob <- ggplotGrob(netDOCpercent_plot_Navicula)
      C323_cumDOCpercent_grob <- ggplotGrob(cumDOCpercent_plot_C323)
      Navi_cumDOCpercent_grob <- ggplotGrob(cumDOCpercent_plot_Navicula)
      empty_grob <- ggplotGrob(empty_plot)
      
      DOC_col1 <- rbind(Navi_netDOCpercent_grob, Navi_cumDOCpercent_grob, size = "first")
      DOC_col2 <- rbind(C323_netDOCpercent_grob, C323_cumDOCpercent_grob, size = "first")
      DOC_col3 <- rbind(D046_netDOCpercent_grob, empty_grob, size = "first")  # need a placeholder plot
      
      DOC_col1$widths <- unit.pmax(Navi_netDOCpercent_grob$widths, Navi_cumDOCpercent_grob$widths)
      
      netcumDOCpercent_final <- cbind(DOC_col1, DOC_col2, DOC_col3, size = "first")
      
      grid.newpage()
      grid.draw(netcumDOCpercent_final)
      
      # Save plot as pdf file
      pdf("Figures/Fig4_netcumDOCpercent_grid.pdf", 
          width = 6.5, height = 4.5)
      grid.draw(netcumDOCpercent_final)
      dev.off()    


#### Figure S7: Absolute DOC release rates (mM C/day) over time (C323 & Navicula only) ####
      
      # C323 ####    
      C323_DOCrate_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$DOC_rate_mean), ],   
                              aes(x = as.numeric(DaysElapsed_mean), y = DOC_rate_mean/1000, color = Treatment, shape = Treatment))  +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        geom_errorbar(aes(ymin=(DOC_rate_mean - DOC_rate_sd)/1000, ymax=(DOC_rate_mean + DOC_rate_sd)/1000), width= 0.4) +
        labs(x = "Days") +
        scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
        scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
        scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
        scale_y_continuous(limits = c(-0.02,2), breaks = seq(0, 2, 0.5), expand = expand_scale(mult = c(0.05,0))) +
        annotate("text", x = 2.4, y = 1.7, label = expression(paste(bold("B"), italic("  S. sourniae"), " C323")), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.position = c(0.92,0.23),
               axis.title.y = element_blank(),
               axis.title.x = element_text(margin = margin(r = 15), size = 14),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "gray87", size = 0.2),
               panel.border = element_rect(color = "gray60", fil = NA),
               axis.ticks = element_blank(),
               axis.text = element_text(size = 12)) 
      
      # Navicula ####  
      Navi_DOCrate_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$DIC_mean),],   
                              aes(x = as.numeric(DaysElapsed_mean), y = DOC_rate_mean/1000, color = Treatment, shape = Treatment))  +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        geom_errorbar(aes(ymin=(DOC_rate_mean - DOC_rate_sd)/1000, ymax=(DOC_rate_mean + DOC_rate_sd)/1000), width= 0.4) +
        labs(y = "DOC Release Rate (mM/day)") +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
        scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
        scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
        scale_y_continuous(limits = c(0,7),breaks = seq(0, 7, 2), expand = expand_scale(mult = c(0.1,0.1))) +
        annotate("text", x = 1.5, y = 6.5, label = expression(paste( bold("A"), italic("  Navicula"), " sp.")), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               legend.position = c(0.93,0.7),
               axis.title.y = element_text(margin = margin(r = 15), size = 12), 
               axis.title.x = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "gray87", size = 0.2),
               panel.border = element_rect(color = "gray60", fil = NA),
               axis.ticks = element_blank(),
               axis.text.y = element_text(size = 12),
               axis.text.x = element_blank()) 
      
      # Combine DOC rate plots in a grid ####
      
      # Convert ggplots to grobs
      C323_DOCrate_grob <- ggplotGrob(C323_DOCrate_plot)
      Navi_DOCrate_grob <- ggplotGrob(Navi_DOCrate_plot)
      
      # Combine in grid
      DOCrate_gridplot <- rbind(Navi_DOCrate_grob, C323_DOCrate_grob, size = "first")
      
      grid.newpage()
      grid.draw(DOCrate_gridplot)
      
      # Save plot as pdf file
      pdf("Figures/FigS7_DOCrateGrid.pdf", 
          width = 9, height = 4)
      grid.draw(DOCrate_gridplot)
      dev.off() 
     
#### Figure S8: DOC rate (% of TOC) over time ####

      # C323 ####    
      C323_DOCratepercent_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "C323" & !is.na(growth_df_avgs$fraction_DOC_rate_mean), ],   
                                  aes(x = as.numeric(DaysElapsed_mean), y = fraction_DOC_rate_mean*100, color = Treatment, shape = Treatment))  +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        geom_errorbar(aes(ymin=(fraction_DOC_rate_mean - fraction_DOC_rate_sd)*100, ymax=(fraction_DOC_rate_mean + fraction_DOC_rate_sd)*100), width= 0.4) +
        labs(x = "Days") +
        scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
        scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
        scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
        scale_y_continuous(limits = c(-50,60), breaks = seq(-50, 60, 25), expand = expand_scale(mult = c(0.1,0.1))) +
        annotate("text", x = 2.7, y = 54, label = expression(paste(bold("C"), italic("  S. sourniae"), " C323")), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_text(margin = margin(r = 15), size = 14),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "gray87", size = 0.2),
               panel.border = element_rect(color = "gray60", fil = NA),
               axis.ticks = element_blank(),
               axis.text = element_text(size = 12)) 
      
      # Navicula ####  
      Navi_DOCratepercent_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "Navicula" & !is.na(growth_df_avgs$fraction_DOC_rate_mean),],   
                                  aes(x = as.numeric(DaysElapsed_mean), y = fraction_DOC_rate_mean*100, color = Treatment, shape = Treatment))  +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        geom_errorbar(aes(ymin=(fraction_DOC_rate_mean - fraction_DOC_rate_sd)*100, ymax=(fraction_DOC_rate_mean + fraction_DOC_rate_sd)*100), width= 0.4) +
        labs(y = "DOC Release Rate (% of TOC)") +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
        scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
        scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
        scale_y_continuous(limits = c(0,90),breaks = seq(0, 90, 20), expand = expand_scale(mult = c(0,0.1))) +
        annotate("text", x = 1.6, y = 84, label = expression(paste( bold("B"), italic("  Navicula"), " sp.")), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               axis.title.y = element_text(margin = margin(r = 15), size = 14),  # moves axis title away from axis label
               axis.title.x = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "gray87", size = 0.2),
               panel.border = element_rect(color = "gray60", fil = NA),
               axis.ticks = element_blank(),
               axis.text.y = element_text(size = 12),
               axis.text.x = element_blank()) 
      
      # D046 ####
      D046_DOCratepercent_plot <- ggplot(data = growth_df_avgs[growth_df_avgs$Algae == "D046" & !is.na(growth_df_avgs$fraction_DOC_rate_mean), ],   
                               aes(x = as.numeric(DaysElapsed_mean), y = fraction_DOC_rate_mean*100, color = Treatment, shape = Treatment))  +
        geom_point(size = 3) +
        geom_line(aes(group = interaction(Treatment, Round)), size = 1) +
        geom_errorbar(aes(ymin=(fraction_DOC_rate_mean - fraction_DOC_rate_sd)*100, ymax=(fraction_DOC_rate_mean + fraction_DOC_rate_sd)*100), width= 0.4) +
        labs(y = " ") +  # x axis label as space holder for space
        scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
        scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
        scale_x_continuous(limits = c(-0.2,30),breaks = seq(0, 30, 5)) +
        scale_y_continuous(limits = c(0, 43),breaks = seq(0, 43, 10), expand = expand_scale(mult = c(0.05,0.2))) +
        annotate("text", x = 2.5, y = 42, label = expression(paste(bold("A"), italic(" Chlorella "),"sp. D046")), size = 5) + 
        annotate("text", x = c(9, 15, 21, 27), y = 42, label = c("1st reuse", "2nd reuse", "3rd reuse", "4th reuse"), color = "gray40", size = 5) +
        theme( legend.key = element_rect(fill = NA),  
               legend.text = element_text(size = 10),
               legend.title = element_blank(),
               axis.title.y = element_text(margin = margin(r = 15), size = 14),
               axis.title.x = element_blank(),
               panel.background = element_rect(fill = "white"),
               panel.grid.major = element_line(colour = "gray87", size = 0.2),
               panel.border = element_rect(color = "gray60", fil = NA),
               axis.ticks = element_blank(),
               axis.text.y = element_text(size = 12),
               axis.text.x = element_blank()) 
      
      # Combine DOC rate plots in a grid ####
      
      # Convert ggplots to grobs
      C323_DOCratepercent_grob <- ggplotGrob(C323_DOCratepercent_plot)
      Navi_DOCratepercent_grob <- ggplotGrob(Navi_DOCratepercent_plot)
      D046_DOCratepercent_grob <- ggplotGrob(D046_DOCratepercent_plot)
      
      # Combine in grid
      DOCratepercent_gridplot <- rbind(D046_DOCratepercent_grob, Navi_DOCratepercent_grob, C323_DOCratepercent_grob, size = "first")
      
      grid.newpage()
      grid.draw(DOCratepercent_gridplot)
      
      # Save plot as pdf file
      pdf("Figures/FigS8_DOCratepercentGrid.pdf", 
          width = 9.5, height = 5)
      grid.draw(DOCratepercent_gridplot)
      dev.off() 

 #### Figure S9: Released and accumulated DOC produced per Round (Absolute values) #### (Note: "released" = "cumulative", "accumulated" = "net")
      
      # D046 ####
      # Released DOC ####
      netDOC_plot_D046 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "D046", ], 
                                 aes(x = Round, y = DOC_net_mean, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=DOC_net_mean - DOC_net_sd, ymax=DOC_net_mean + DOC_net_sd), width= 0.4, position = position_dodge(0.9)) +
        labs(x = "Medium Reuses", title = expression(paste( italic("Chlorella "), "sp. D046"))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('green3', 'green4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0)), breaks = seq(0,0.3, by = 0.1)) +
        annotate("text", x = 0.75, y = 0.23, label = expression(paste(bold("C"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 2),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title.x = element_text(margin = margin(r = 22),size = 10),
               axis.title.y = element_blank(),
               axis.text = element_text(size = 10))   
     
       # Empty plot ####
      empty_plot <- ggplot(data.frame()) +
        labs(x = " ") +
        theme(panel.background = element_rect(fill = NA))
      
      # C323 ####
      # Accumulated DOC ####  
      netDOC_plot_C323 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "C323", ], 
                                 aes(x = Round, y = DOC_net_mean, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=DOC_net_mean - DOC_net_sd, ymax=DOC_net_mean + DOC_net_sd), width= 0.4, position = position_dodge(0.9)) +
        labs(title = expression(paste(italic("S. sourniae "), " C323"))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('gray63','gray43')) +
        scale_y_continuous(expand = expand_scale(mult = c(0.05,0)), breaks = seq(0,0.3, by = 0.05)) +
        annotate("text", x = 0.7, y = 0.1, label = expression(paste(bold("B"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 1),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_blank(),
               axis.text = element_text(size = 10))
      # Released DOC ####  
      cumDOC_plot_C323 <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "C323", ], 
                                 aes(x = Round, y = DOC_cumulative_mean, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=DOC_cumulative_mean - DOC_cumulative_sd, ymax=DOC_cumulative_mean + DOC_cumulative_sd), width= 0.4, position = position_dodge(0.9)) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('gray63','gray43')) +
        labs(x = "Medium Reuses") +
        scale_y_continuous(expand = expand_scale(mult = c(0,0.05)), breaks = seq(0,3, by = 1)) +        
        annotate("text", x = 0.7, y = 2.8, label = expression(paste(bold("E"))), size = 5) + 
        theme( legend.position = "none",
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.title.x = element_text(margin = margin(r = 22),size = 10),
               axis.text = element_text(size = 10)) 
      
      # Navicula ####  
      # Accumulated DOC ####
      netDOC_plot_Navicula <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "Navicula", ], 
                                     aes(x = Round, y = DOC_net_mean, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=DOC_net_mean - DOC_net_sd, ymax=DOC_net_mean + DOC_net_sd), width= 0.4, position = position_dodge(0.9)) +
        labs(y = "Accumulated DOC (mM)", x = " ", title = expression(paste(italic("Navicula"), " sp."))) +
        scale_fill_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0)), breaks = seq(0,0.6, by = 0.2)) +
        annotate("text", x = 0.75, y = 0.58, label = expression(paste(bold("A"))), size = 5) + 
        theme( legend.key = element_rect(fill = NA),  # removes color from behind legned points/lines
               legend.title = element_blank(),
               legend.text = element_text(size = 10),
               legend.position = "top",
               legend.key.size = unit(3, "mm"),
               legend.justification = "center",
               legend.margin=margin(t = 0, b = 0, l = 1),
               legend.box.margin=margin(t = 0, b = 0),
               legend.spacing.x = unit(1, "mm"),
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_text(margin = margin(r = 15),size = 10),
               axis.text = element_text(size = 10))
      
      # Released DOC ####
      cumDOC_plot_Navicula <- ggplot(data = carbon_cumulative_net_df[carbon_cumulative_net_df$Algae == "Navicula", ], 
                                     aes(x = Round, y = DOC_cumulative_mean, fill = Treatment))  +
        geom_col(position = position_dodge(0.9), color = "gray15") +
        geom_errorbar(aes(ymin=DOC_cumulative_mean - DOC_cumulative_sd, ymax=DOC_cumulative_mean + DOC_cumulative_sd), width= 0.4, position = position_dodge(0.9)) +
        labs(x = "Medium Reuses", y = "Released DOC (mM)") +
        scale_fill_manual(labels = c("Fresh", "Recycled"), 
                          values = c('dodgerblue2', 'dodgerblue4')) +
        scale_y_continuous(expand = expand_scale(mult = c(0,0.05)), breaks = seq(0,15, by = 5)) +
        annotate("text", x = 0.75, y = 16, label = expression(paste(bold("D"))), size = 5) + 
        theme( legend.position = "none",
               plot.title = element_text(hjust = 0.5, size = 11),
               panel.background = element_rect(fill = NA),
               axis.line = element_line(color = "black", size = 0.4),
               axis.ticks.x = element_blank(),
               axis.title = element_text(margin = margin(r = 15),size = 10),
               axis.text = element_text(size = 10))
      
      # Combine net and cumulative DOC plots in grid and save as pdf ####
      
      # Convert ggplots to grobs
      D046_netDOC_grob <- ggplotGrob(netDOC_plot_D046)
      C323_netDOC_grob <- ggplotGrob(netDOC_plot_C323)
      Navi_netDOC_grob <- ggplotGrob(netDOC_plot_Navicula)
      C323_cumDOC_grob <- ggplotGrob(cumDOC_plot_C323)
      Navi_cumDOC_grob <- ggplotGrob(cumDOC_plot_Navicula)
      empty_grob <- ggplotGrob(empty_plot)
      
      DOC_col1 <- rbind(Navi_netDOC_grob, Navi_cumDOC_grob, size = "first")
      DOC_col2 <- rbind(C323_netDOC_grob, C323_cumDOC_grob, size = "first")
      DOC_col3 <- rbind(D046_netDOC_grob, empty_grob, size = "first")  # need a placeholder plot
      
      netcumDOC_final <- cbind(DOC_col1, DOC_col2, DOC_col3, size = "first")
      
      grid.newpage()
      grid.draw(netcumDOC_final)
      
      # Save plot as pdf file
      pdf("Figures/FigS9_netcumDOC_grid.pdf", 
          width = 6.5, height = 4.5)
      grid.draw(netcumDOC_final)
      dev.off()    
     
#### Figure 5: DOC and bacteria concentrations in saved filtrate experiment ####
    # D046 ####
        # Bacteria over time in filtrate ####
        bacteriafiltrate_D046_plot <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == 'D046', ], 
                                         aes(x = ElapsedDays, y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
          geom_point(size = 2.5) +
          geom_line(size = 1) +
          geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 10) +
          labs(y = bquote('Bacteria'~(10^6~'cells/mL')), x = " ") +
          scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
          scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) + 
          scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260)) +
          scale_y_continuous(limits = c(0,8), expand = expand_scale(mult = c(0,0.05))) +
          annotate("text", x = 10, y = 7.8, label = expression(paste(bold("D"))), size = 4.5) + 
          theme( legend.position = "none",
                 axis.title = element_text(margin = margin(r = 15)),  # moves axis title away from axis label
                 panel.grid.major = element_line(colour = "gray87", size = 0.2),
                 panel.background = element_rect(fill = NA),
                 axis.line = element_line(color = "black", size = 0.4),
                 axis.text = element_text(size = 12),
                 axis.ticks = element_blank()) 
        
        # DOC over time in filtrate ####
            filtrateDOC_plot_D046 <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == "D046", ], 
                                            aes(x = as.numeric(ElapsedDays), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
              geom_point(size = 2.5) +
              geom_line(size = 1) +
              geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 10) +
              scale_color_manual(labels = c("Fresh", "Recycled"), values = c('green3', 'green4'), guide = guide_legend(reverse=TRUE)) +
              scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
              labs(title = expression(paste(italic("Chlorella"), " sp. D046")), y = "DOC (mM)") +
              scale_y_continuous(expand = expand_scale(mult = c(0,0)), limits = c(0,1), breaks = seq(0,1, by = 0.2)) +
              scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260), expand = expand_scale(mult = c(0.05,0))) +
              annotate("text", x = 10, y = 0.94, label = expression(paste(bold("A"))), size = 4.5) + 
              theme( legend.key = element_rect(fill = NA),  
                     legend.title = element_blank(),
                     legend.text = element_text(size = 10),
                     legend.position = "top",
                     legend.key.size = unit(4, "mm"),
                     legend.justification = "center",
                     legend.margin=margin(t = 0, b = 0, l = 2),
                     legend.box.margin=margin(t = 0, b = 0),
                     legend.spacing.x = unit(1, "mm"),
                     plot.title = element_text(hjust = 0.5, size = 11),
                     panel.background = element_rect(fill = NA),
                     axis.line = element_line(color = "black", size = 0.4),
                     panel.grid.major = element_line(colour = "gray87", size = 0.2),
                     axis.ticks = element_blank(),
                     axis.title.y = element_text(margin = margin(r = 15),size = 11),
                     axis.title.x = element_blank(),                            
                     axis.text = element_text(size = 11))   
      
    # C323 ####    
        # Bacteria over time in filtrate ####
          bacteriafiltrate_C323_plot <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == 'C323', ], 
                                           aes(x = ElapsedDays, y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
            geom_point(size = 2.5) +
            geom_line(size = 1) +
            geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 10) +
            labs(y = bquote('Bacteria'~(10^6~'cells/mL'))) +
            scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
            scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) + 
            scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260), expand = expand_scale(mult = c(0.05,0))) +
            scale_y_continuous(limits = c(0,12), expand = expand_scale(mult = c(0.05,0))) +
            annotate("text", x = 10, y = 11, label = expression(paste(bold("F"))), size = 4.5) + 
            theme( legend.position = "none",
                   axis.title = element_blank(),
                   panel.grid.major = element_line(colour = "gray87", size = 0.2),
                   panel.background = element_rect(fill = NA),
                   axis.line = element_line(color = "black", size = 0.4),
                   axis.text = element_text(size = 12),
                   axis.ticks = element_blank()) 
          
        # DOC over time in filtrate ####
          filtrateDOC_plot_C323 <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == "C323", ], 
                                          aes(x = as.numeric(ElapsedDays), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
            geom_point(size = 2.5) +
            geom_line(size = 1) +
            geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 10) +
            scale_color_manual(labels = c("Fresh", "Recycled"), values = c('gray63','gray43'), guide = guide_legend(reverse=TRUE)) +
            scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
            labs(x = "Medium Reuses", title = expression(paste(italic("  S. sourniae"), " C323")), y = "DOC (µM)") +
            scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260), expand = expand_scale(mult = c(0.05,0))) +
            scale_y_continuous(expand = expand_scale(mult = c(0,0)), limits = c(0, 0.2), breaks = seq(0,0.2, by = 0.05)) +
            annotate("text", x = 10, y = 0.19, label = expression(paste(bold("C"))), size = 4.5) + 
            theme( legend.key = element_rect(fill = NA),  
                   legend.title = element_blank(),
                   legend.text = element_text(size = 10),
                   legend.position = "top",
                   legend.key.size = unit(4, "mm"),
                   legend.justification = "center",
                   legend.margin=margin(t = 0, b = 0, l = 4),
                   legend.box.margin=margin(t = 0, b = 0),
                   legend.spacing.x = unit(1, "mm"),
                   plot.title = element_text(hjust = 0.5, size = 11),
                   panel.background = element_rect(fill = NA),
                   axis.line = element_line(color = "black", size = 0.4),
                   panel.grid.major = element_line(colour = "gray87", size = 0.2),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),                            
                   axis.text = element_text(size = 11))   
        
    # Navi ####
        # Bacteria over time in filtrate ####
          bacteriafiltrate_Navi_plot <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == 'Navicula', ], 
                                               aes(x = ElapsedDays, y = BacteriaConc_mean, color = Treatment, shape = Treatment))  +
            geom_point(size = 2.5) +
            geom_line(size = 1) +
            geom_errorbar(aes(ymin=BacteriaConc_mean - BacteriaConc_sd, ymax=BacteriaConc_mean + BacteriaConc_sd), width= 10) +
            labs(x = "Days") +
            scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
            scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) + 
            scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260), expand = expand_scale(mult = c(0.05,0))) +
            scale_y_continuous(limits = c(0,8), breaks = seq(0, 8, by = 2), expand = expand_scale(mult = c(0.05,0))) +
            annotate("text", x = 10, y = 7.5, label = expression(paste(bold("E"))), size = 4.5) + 
            theme( legend.position = "none",
                   axis.title.x = element_text(margin = margin(r = 15), size = 12),  
                   axis.title.y = element_blank(),
                   panel.grid.major = element_line(colour = "gray87", size = 0.2),
                   panel.background = element_rect(fill = NA),
                   axis.line = element_line(color = "black", size = 0.4),
                   axis.text = element_text(size = 12),
                   axis.ticks = element_blank()) 
          
        # DOC over time in filtrate ####
          filtrateDOC_plot_Navi <- ggplot(data = filtrate_df_avgs[filtrate_df_avgs$Algae == "Navicula", ], 
                                          aes(x = as.numeric(ElapsedDays), y = DOC_mean/1000, color = Treatment, shape = Treatment))  +
            geom_point(size = 2.5) +
            geom_line(size = 1) +
            geom_errorbar(aes(ymin=(DOC_mean - DOC_sd)/1000, ymax=(DOC_mean + DOC_sd)/1000), width= 10) +
            scale_color_manual(labels = c("Fresh", "Recycled"), values = c('dodgerblue2', 'dodgerblue4'), guide = guide_legend(reverse=TRUE)) +
            scale_shape_manual(labels = c("Fresh", "Recycled"), values = c(16, 15), guide = guide_legend(reverse=TRUE)) +  
            labs(title = expression(paste(italic("  Navicula"), " sp."))) +
            scale_x_continuous(breaks = seq(0, 250, 50), limits = c(-5,260), expand = expand_scale(mult = c(0.05,0))) +
            scale_y_continuous(expand = expand_scale(mult = c(0,0)), limits = c(0, 2.5), breaks = seq(0,2.5, by = 0.5)) +
            annotate("text", x = 10, y = 2.35, label = expression(paste(bold("B"))), size = 4.5) + 
            theme( legend.key = element_rect(fill = NA),  
                   legend.title = element_blank(),
                   legend.text = element_text(size = 10),
                   legend.position = "top",
                   legend.key.size = unit(4, "mm"),
                   legend.justification = "center",
                   legend.margin=margin(t = 0, b = 0, l = 4),
                   legend.box.margin=margin(t = 0, b = 0),
                   legend.spacing.x = unit(1, "mm"),
                   plot.title = element_text(hjust = 0.5, size = 11),
                   panel.background = element_rect(fill = NA),
                   axis.line = element_line(color = "black", size = 0.4),
                   panel.grid.major = element_line(colour = "gray87", size = 0.2),
                   axis.ticks = element_blank(),
                   axis.title = element_blank(),                            
                   axis.text = element_text(size = 11))         
        
    # Combine filtrate plots and save plot as pdf ####
          # Convert plots to grobs
          bacteriafiltrate_D046_grob <- ggplotGrob(bacteriafiltrate_D046_plot)
          bacteriafiltrate_C323_grob <- ggplotGrob(bacteriafiltrate_C323_plot)
          bacteriafiltrate_Navi_grob <- ggplotGrob(bacteriafiltrate_Navi_plot)
          filtrateDOC_D046_grob <- ggplotGrob(filtrateDOC_plot_D046)
          filtrateDOC_C323_grob <- ggplotGrob(filtrateDOC_plot_C323)
          filtrateDOC_Navi_grob <- ggplotGrob(filtrateDOC_plot_Navi)
          
          # Make plot columns
          filtrate_col1 <- rbind(filtrateDOC_D046_grob, bacteriafiltrate_D046_grob, size = "first")
          filtrate_col2 <- rbind(filtrateDOC_Navi_grob, bacteriafiltrate_Navi_grob, size = "first")
          filtrate_col3 <- rbind(filtrateDOC_C323_grob, bacteriafiltrate_C323_grob, size = "first")
          
          # Set widths
          filtrate_col1$widths <- unit.pmax(filtrateDOC_D046_grob$widths, bacteriafiltrate_D046_grob$widths)
          filtrate_col2$widths <- unit.pmax(filtrateDOC_Navi_grob$widths, bacteriafiltrate_Navi_grob$widths)
          filtrate_col3$widths <- unit.pmax(filtrateDOC_C323_grob$widths, bacteriafiltrate_C323_grob$widths)
          
          # Combine plot columns
          saved_filtrate_plot <- cbind(filtrate_col1, filtrate_col2, filtrate_col3, size = "first")
          
          grid.newpage()
          grid.draw(saved_filtrate_plot)
          
          # Save plot as pdf
          pdf("Figures/Fig5_saved_filtrate.pdf", 
              width = 8, height = 4.5)
          grid.draw(saved_filtrate_plot)
          dev.off() 

#### Figure S10: Biomass response (PC_net normalized by mean PC_net in Fresh medium) versus initial DOC in the recycled medium ####

DOC_biomass_plot <- ggplot(data = initialDOC_PC_final, 
                           aes(x = initial_DOC, y = PC_R_normalized, color = Algae, shape = Algae))  +
  geom_point(alpha = 0.6, size = 3) +
  labs(x = "Initial DOC (µM) in Recycled Medium", y = "Normalized Biomass Yield in Recycled Medium") +
  scale_color_manual(labels = c(expression(paste(italic("S. sourniae"), " C323")), expression(paste(italic("Chlorella"), " sp. D046")), expression(paste(italic("Navicula"), " sp."))),
                     values = c('gray43','green4','dodgerblue4')) +
  scale_shape_manual(labels = c(expression(paste(italic("S. sourniae"), " C323")), expression(paste(italic("Chlorella"), " sp. D046")), expression(paste(italic("Navicula"), " sp."))), 
                     values = c(16, 17, 15)) +  
  geom_hline(yintercept = 1, color = "gray55", size = 0.5) +   # reference line for no effect
  xlim(0,1750) +
  scale_y_continuous(breaks = seq(0, 1.5, 0.5), limits = c(0,1.5)) +
  theme( legend.background = element_rect(fill = "white"),  # removes color from behind legned points/lines
         legend.box.background = element_rect(color = "gray70"),
         legend.key = element_rect(fill = "white"),
         legend.title = element_blank(),
         legend.text = element_text(size = 9),
         legend.justification = c(1,0),
         legend.position = c(0.95,0.05),
         legend.text.align = 0,
         axis.title.y = element_text(margin = margin(r = 5), size = 12),  
         axis.title.x = element_text(margin = margin(t = 5), size = 12),  
         panel.background = element_rect(fill = NA),
         panel.grid.major = element_line(colour = "gray87", size = 0.2),
         panel.grid.minor = element_line(colour = "gray90", size = 0.2),
         panel.border = element_rect(color = "gray60", fil = NA), 
         axis.ticks = element_blank(),
         axis.text = element_text(size = 12)) 

    # Save plot as pdf ####
        pdf("Figures/FigS10_BiomassvsDOC.pdf", 
            width = 4, height = 4)
        print(DOC_biomass_plot)
        dev.off()
