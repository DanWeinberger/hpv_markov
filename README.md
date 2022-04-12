# hpv_markov

### Guinevere Oliver
### Last updated: April 12, 2022

This repository contains the code associated with the thesis titled "Potential impact of COVID-19-related disruptions to cervical cancer screening on future disease burden: A modeling study." The thesis was written to satisfy graduation requirements for the Yale School of Public Health, Epidemiology of Microbial Disease Department.

### Objective
The objective of this analysis was to offer insight into possible future trends in HSIL cases, cervical cancer cases, and cervical cancer deaths.

### R scripts

##### analysis_no COVID.R
The function in this script creates a matrix with the prevalence of cervical lesions, cancer, and cancer death in a scenario with no COVID (baseline).

##### analysis.R
The function in this script creates a matrix with the prevalence of cervical lesions, cancer, and cancer death in a scenario with different levels of disruption to screening due to COVID.

##### analysis by race_no COVID.R
The function in this script creates a matrix with the prevalence of cervical lesions, cancer, and cancer death in a scenario with no COVID (baseline) for a single racial/ethnic category based on the percentage of people identifying as that race/ethnicity in the U.S., the screening rates for people identifying as that race/ethnicity, and vaccination rate for people identifying as that race/ethnicity.

##### analysis by race.R
The function in this script creates a matrix with the prevalence of cervical lesions, cancer, and cancer death in a scenario with different levels of disruption to screening due to COVID for a single racial/ethnic category based on the percentage of people identifying as that race/ethnicity in the U.S., the screening rates for people identifying as that race/ethnicity, and vaccination rate for people identifying as that race/ethnicity.

##### profile_likelihood_estimation.R
This script provides an example of how profile likelihood estimation was used to identify the best parameters for the model. It includes fitting for an earlier, more simplified version of the model.

##### change_estimates.R
This script contains code to generate the absolute change and percentage change estimates for HSIL cases, cervical cancer cases, and cervical cancer deaths between the "no disruptions due to COVID" and "disruptions due to COVID" lines. It also contains code to generate the estimates for the cervical cancer:HSIL indicator.

##### plots.R
This script contains code to generate plots of the counts of HSIL cases, cervical cancer cases, and cervical cancer deaths for each analysis. It also contains code to generate the plot for the cervical cancer:HSIL indicator.
