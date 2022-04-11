
############################################################################################################

######### Change estimates (COVID vs no COVID)
######### Guinevere Oliver
######### April 11, 2022

############################################################################################################


#################################### Compare COVID to no COVID ####################################

################### 40%

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(screen_0.4$Year, screen_0.4$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 2178586
HSIL_noCov_num_5y <- 2119390
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_abs_change_5y
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(screen_0.4$Year, screen_0.4$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 12603908
HSIL_noCov_num_30y <- 12534868
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_abs_change_30y
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 892440.8
HSIL_det_noCov_num_5y <- 913000.8
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_abs_change_5y
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 5369703
HSIL_det_noCov_num_30y <- 5386641
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_abs_change_30y
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(screen_0.4$Year, screen_0.4$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 206802.7
cancer_all_noCov_num_5y <- 198016.2
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_abs_change_5y
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(screen_0.4$Year, screen_0.4$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 1164670
cancer_all_noCov_num_30y <- 1151903
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_abs_change_30y
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 50561.41
cancer_det_noCov_num_5y <- 48678.25
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_abs_change_5y
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 285261.7
cancer_det_noCov_num_30y <- 282207.5
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_abs_change_30y
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 17727.02
cancer_death_noCov_num_5y <- 17161.86
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_abs_change_5y
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(screen_0.4$Year, screen_0.4$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 100097.5
cancer_death_noCov_num_30y <- 99020.17
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_abs_change_30y
cancer_death_perc_change_30y




################### 60%

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(screen_0.6$Year, screen_0.6$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 2208189
HSIL_noCov_num_5y <- 2119390
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_abs_change_5y
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(screen_0.6$Year, screen_0.6$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 12638417
HSIL_noCov_num_30y <- 12534868
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_abs_change_30y
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 882162.4
HSIL_det_noCov_num_5y <- 913000.8
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_abs_change_5y
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 5361228
HSIL_det_noCov_num_30y <- 5386641
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_abs_change_30y
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(screen_0.6$Year, screen_0.6$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 211199.8
cancer_all_noCov_num_5y <- 198016.2
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_abs_change_5y
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(screen_0.6$Year, screen_0.6$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 1171051
cancer_all_noCov_num_30y <- 1151903
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_abs_change_30y
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 51505.43
cancer_det_noCov_num_5y <- 48678.25
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_abs_change_5y
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 286789.1
cancer_det_noCov_num_30y <- 282207.5
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_abs_change_30y
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 18010.58
cancer_death_noCov_num_5y <- 17161.86
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_abs_change_5y
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(screen_0.6$Year, screen_0.6$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 100643.2
cancer_death_noCov_num_30y <- 99020.17
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_abs_change_30y
cancer_death_perc_change_30y



################### 80%

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(screen_0.8$Year, screen_0.8$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 2237785
HSIL_noCov_num_5y <- 2119390
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_abs_change_5y
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(screen_0.8$Year, screen_0.8$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 12672936
HSIL_noCov_num_30y <- 12534868
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_abs_change_30y
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 871886.3
HSIL_det_noCov_num_5y <- 913000.8
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_abs_change_5y
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 5352770
HSIL_det_noCov_num_30y <- 5386641
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_abs_change_30y
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(screen_0.8$Year, screen_0.8$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 215585.6
cancer_all_noCov_num_5y <- 198016.2
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_abs_change_5y
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(screen_0.8$Year, screen_0.8$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 1177415
cancer_all_noCov_num_30y <- 1151903
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_abs_change_30y
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 52443.8
cancer_det_noCov_num_5y <- 48678.25
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_abs_change_5y
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 288309.4
cancer_det_noCov_num_30y <- 282207.5
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_abs_change_30y
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 18291.55
cancer_death_noCov_num_5y <- 17161.86
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_abs_change_5y
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(screen_0.8$Year, screen_0.8$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(result_tot_noCov$Year, result_tot_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 101175
cancer_death_noCov_num_30y <- 99020.17
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_abs_change_30y
cancer_death_perc_change_30y



#################################### Analysis by Race ####################################


################### White non-Hispanic

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 1306490
HSIL_noCov_num_5y <- 1259015
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 7471033
HSIL_noCov_num_30y <- 7414902
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 513763.7
HSIL_det_noCov_num_5y <- 531583.6
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 3108160
HSIL_det_noCov_num_30y <- 3122973
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 144896
cancer_all_noCov_num_5y <- 139926.1
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 814863.1
cancer_all_noCov_num_30y <- 807024.2
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 34532.63
cancer_det_noCov_num_5y <- 33497.47
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 194153
cancer_det_noCov_num_30y <- 192350.5
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 12154.84
cancer_death_noCov_num_5y <- 11833.18
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 68177.08
cancer_death_noCov_num_30y <- 67544.09
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_perc_change_30y


################### Black non-Hispanic

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 268244.8
HSIL_noCov_num_5y <- 258971.8
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 1545440
HSIL_noCov_num_30y <- 1534519
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 108258.8
HSIL_det_noCov_num_5y <- 111911.7
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 658555.9
HSIL_det_noCov_num_30y <- 661579.6
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 23438.97
cancer_all_noCov_num_5y <- 22595.76
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 133261.9
cancer_all_noCov_num_30y <- 131957
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 5784
cancer_det_noCov_num_5y <- 5602.502
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 32934.77
cancer_det_noCov_num_30y <- 32624.6
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 2034.402
cancer_death_noCov_num_5y <- 1975.222
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 11562.79
cancer_death_noCov_num_30y <- 11457.61
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_perc_change_30y



################### Hispanic

### HSIL All

# Until 2025
HSIL_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_HSIL), 2020,2025)
HSIL_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_HSIL), 2020,2025)
HSIL_Cov_5y
HSIL_noCov_5y
HSIL_Cov_num_5y <- 387069.9
HSIL_noCov_num_5y <- 373114.4
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 2222872
HSIL_noCov_num_30y <- 2206350
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 148392.1
HSIL_det_noCov_num_5y <- 153738
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 902523.5
HSIL_det_noCov_num_30y <- 906977.5
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 41617.17
cancer_all_noCov_num_5y <- 40213.71
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 235961.6
cancer_all_noCov_num_30y <- 233727.2
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 9905.473
cancer_det_noCov_num_5y <- 9614.13
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 56196.87
cancer_det_noCov_num_30y <- 55678.93
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 3479.997
cancer_death_noCov_num_5y <- 3387.07
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 19712.59
cancer_death_noCov_num_30y <- 19523.92
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_perc_change_30y






