

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
HSIL_Cov_num_5y <- 1293102
HSIL_noCov_num_5y <- 1234999
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 7349669
HSIL_noCov_num_30y <- 7291124
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 656242.8
HSIL_det_noCov_num_5y <- 655356.5
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 3860504
HSIL_det_noCov_num_30y <- 3860720
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 243859.9
cancer_all_noCov_num_5y <- 245093.1
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 1434001
cancer_all_noCov_num_30y <- 1433383
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 118679.9
cancer_det_noCov_num_5y <- 118938.2
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 694567.6
cancer_det_noCov_num_30y <- 693934.1
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_white$Year, Results_white$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 40922.94
cancer_death_noCov_num_5y <- 41755.94
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_white$Year, Results_white$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_white_noCov$Year, Results_white_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 243259.1
cancer_death_noCov_num_30y <- 243331
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
HSIL_Cov_num_5y <- 257327.1
HSIL_noCov_num_5y <- 246932.6
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 1475416
HSIL_noCov_num_30y <- 1465166
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 137250.2
HSIL_det_noCov_num_5y <- 137357.4
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 813234.1
HSIL_det_noCov_num_30y <- 813790.4
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 40049.19
cancer_all_noCov_num_5y <- 40186.09
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 236832
cancer_all_noCov_num_30y <- 236765.3
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 19937.03
cancer_det_noCov_num_5y <- 19959.32
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 117549.3
cancer_det_noCov_num_30y <- 117476.2
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_black$Year, Results_black$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 6879.559
cancer_death_noCov_num_5y <- 7010.79
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_black$Year, Results_black$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_black_noCov$Year, Results_black_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 41176.08
cancer_death_noCov_num_30y <- 41190.32
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
HSIL_Cov_num_5y <- 380463.6
HSIL_noCov_num_5y <- 363648.4
HSIL_abs_change_5y <- HSIL_Cov_num_5y - HSIL_noCov_num_5y
HSIL_perc_change_5y <- HSIL_abs_change_5y/HSIL_noCov_num_5y*100
HSIL_perc_change_5y

# Until 2050
HSIL_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_HSIL), 2020,2050)
HSIL_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_HSIL), 2020,2050)
HSIL_Cov_30y
HSIL_noCov_30y
HSIL_Cov_num_30y <- 2171262
HSIL_noCov_num_30y <- 2154297
HSIL_abs_change_30y <- HSIL_Cov_num_30y - HSIL_noCov_num_30y
HSIL_perc_change_30y <- HSIL_abs_change_30y/HSIL_noCov_num_30y*100
HSIL_perc_change_30y


### HSIL Detected

# Until 2025
HSIL_det_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_HSIL), 2020,2025)
HSIL_det_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_HSIL), 2020,2025)
HSIL_det_Cov_5y
HSIL_det_noCov_5y
HSIL_det_Cov_num_5y <- 190938.2
HSIL_det_noCov_num_5y <- 190820.6
HSIL_det_abs_change_5y <- HSIL_det_Cov_num_5y - HSIL_det_noCov_num_5y
HSIL_det_perc_change_5y <- HSIL_det_abs_change_5y/HSIL_det_noCov_num_5y*100
HSIL_det_perc_change_5y

# Until 2050
HSIL_det_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_HSIL), 2020,2050)
HSIL_det_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_HSIL), 2020,2050)
HSIL_det_Cov_30y
HSIL_det_noCov_30y
HSIL_det_Cov_num_30y <- 1128172
HSIL_det_noCov_num_30y <- 1128346
HSIL_det_abs_change_30y <- HSIL_det_Cov_num_30y - HSIL_det_noCov_num_30y
HSIL_det_perc_change_30y <- HSIL_det_abs_change_30y/HSIL_det_noCov_num_30y*100
HSIL_det_perc_change_30y


### Cancer All

# Until 2025
cancer_all_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_Cancer), 2020,2025)
cancer_all_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_Cancer), 2020,2025)
cancer_all_Cov_5y
cancer_all_noCov_5y
cancer_all_Cov_num_5y <- 70587.71
cancer_all_noCov_num_5y <- 70959.67
cancer_all_abs_change_5y <- cancer_all_Cov_num_5y - cancer_all_noCov_num_5y
cancer_all_perc_change_5y <- cancer_all_abs_change_5y/cancer_all_noCov_num_5y*100
cancer_all_perc_change_5y

# Until 2050
cancer_all_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$All_Cancer), 2020,2050)
cancer_all_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$All_Cancer), 2020,2050)
cancer_all_Cov_30y
cancer_all_noCov_30y
cancer_all_Cov_num_30y <- 417279
cancer_all_noCov_num_30y <- 417127.3
cancer_all_abs_change_30y <- cancer_all_Cov_num_30y - cancer_all_noCov_num_30y
cancer_all_perc_change_30y <- cancer_all_abs_change_30y/cancer_all_noCov_num_30y*100
cancer_all_perc_change_30y


### Cancer Detected

# Until 2025
cancer_det_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_Cancer), 2020,2025)
cancer_det_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_Cancer), 2020,2025)
cancer_det_Cov_5y
cancer_det_noCov_5y
cancer_det_Cov_num_5y <- 34327.49
cancer_det_noCov_num_5y <- 34412.58
cancer_det_abs_change_5y <- cancer_det_Cov_num_5y - cancer_det_noCov_num_5y
cancer_det_perc_change_5y <- cancer_det_abs_change_5y/cancer_det_noCov_num_5y*100
cancer_det_perc_change_5y

# Until 2050
cancer_det_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Det_Cancer), 2020,2050)
cancer_det_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Det_Cancer), 2020,2050)
cancer_det_Cov_30y
cancer_det_noCov_30y
cancer_det_Cov_num_30y <- 202075.2
cancer_det_noCov_num_30y <- 201892.4
cancer_det_abs_change_30y <- cancer_det_Cov_num_30y - cancer_det_noCov_num_30y
cancer_det_perc_change_30y <- cancer_det_abs_change_30y/cancer_det_noCov_num_30y*100
cancer_det_perc_change_30y


### Cancer Death

# Until 2025
cancer_death_Cov_5y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Cancer_Death), 2020,2025)
cancer_death_noCov_5y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Cancer_Death), 2020,2025)
cancer_death_Cov_5y
cancer_death_noCov_5y
cancer_death_Cov_num_5y <- 11824.44
cancer_death_noCov_num_5y <- 12068.35
cancer_death_abs_change_5y <- cancer_death_Cov_num_5y - cancer_death_noCov_num_5y
cancer_death_perc_change_5y <- cancer_death_abs_change_5y/cancer_death_noCov_num_5y*100
cancer_death_perc_change_5y

# Until 2050
cancer_death_Cov_30y <- integrate(splinefun(Results_hisp$Year, Results_hisp$Cancer_Death), 2020,2050)
cancer_death_noCov_30y <- integrate(splinefun(Results_hisp_noCov$Year, Results_hisp_noCov$Cancer_Death), 2020,2050)
cancer_death_Cov_30y
cancer_death_noCov_30y
cancer_death_Cov_num_30y <- 70740.68
cancer_death_noCov_num_30y <- 70745.21
cancer_death_abs_change_30y <- cancer_death_Cov_num_30y - cancer_death_noCov_num_30y
cancer_death_perc_change_30y <- cancer_death_abs_change_30y/cancer_death_noCov_num_30y*100
cancer_death_perc_change_30y






