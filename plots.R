
### Project Plots
### Guinevere Oliver
### March 12

#################################### Compare COVID to no COVID ####################################

# Normal
p_Norm<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Normal, alpha="Normal no disruption"), color="red") + 
  geom_line(data=result_tot, aes(y=Normal, alpha="Normal"), color="red") + 
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(10000000,20000000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  scale_alpha_manual(values=c("Normal no disruption"=0.3, "Normal"=1), name="") +
  theme_classic()
p_Norm


# LSIL
p_LSIL<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_LSIL, alpha="No disruption"), color="red", linetype="solid") + 
  geom_line(data=result_tot, aes(y=Undet_LSIL, alpha="Disruption", linetype="Undetected"), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_LSIL, alpha="No disruption", linetype="Detected"), colour="red") +
  geom_line(data=result_tot, aes(y=Det_LSIL, alpha="Disruption",linetype="Detected"), colour="red") +
  geom_line(data=result_tot_noCov, aes(y=All_LSIL, alpha="No disruption", color="Combined Undetected and Detected")) +
  geom_line(data=result_tot, aes(y=All_LSIL, alpha="Disruption", color="Combined Undetected and Detected")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,1500000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  labs(y="LSIL Cases") +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_color_manual(values=c("Combined Undetected and Detected"="black"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="red", alpha=c(0.1,1)))) +
  theme_classic()
p_LSIL


# HSIL
p_HSIL<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_HSIL, alpha="No disruption"), color="red", linetype="solid") + 
  geom_line(data=result_tot, aes(y=Undet_HSIL, alpha="Disruption", linetype="Undetected"), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_HSIL, alpha="No disruption", linetype="Detected"), colour="red") +
  geom_line(data=result_tot, aes(y=Det_HSIL, alpha="Disruption",linetype="Detected"), colour="red") +
  geom_line(data=result_tot_noCov, aes(y=All_HSIL, alpha="No disruption", color="Combined Undetected and Detected")) +
  geom_line(data=result_tot, aes(y=All_HSIL, alpha="Disruption", color="Combined Undetected and Detected")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,750000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  labs(y="HSIL Cases") +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_color_manual(values=c("Combined Undetected and Detected"="black"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="red", alpha=c(0.1,1)))) +
  theme_classic()
p_HSIL

# Cancer Cases and Deaths
p_Cancer<-ggplot(data=result_tot, aes(x=Year)) +
  geom_line(data=result_tot_noCov, aes(y=Undet_Cancer,alpha="No disruption"), colour="red") + 
  geom_line(data=result_tot, aes(y=Undet_Cancer,alpha="Disruption",linetype="Undetected"), colour="red") + 
  geom_line(data=result_tot_noCov, aes(y=Det_Cancer,alpha="No disruption"), colour="red") +
  geom_line(data=result_tot, aes(y=Det_Cancer,alpha="Disruption"), colour="red") +
  geom_line(data=result_tot_noCov, aes(y=Cancer_Death,alpha="No disruption"), colour="red") +
  geom_line(data=result_tot, aes(y=Cancer_Death,alpha="Disruption",linetype="Cervical Cancer Death"), colour="red") +
  geom_line(data=result_tot_noCov, aes(y=All_Cancer,alpha="No disruption",color="Combined Undetected and Detected")) +
  geom_line(data=result_tot, aes(y=All_Cancer,alpha="Disruption",color="Combined Undetected and Detected")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,55000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) + 
  labs(y="Cancer Cases/Deaths") +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid","Cervical Cancer Death"="dotted"), name="") +
  scale_color_manual(values=c("Combined Undetected and Detected"="black"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="red", alpha=c(0.1,1)))) +
  theme_classic()
p_Cancer



#################################### Scenario Analysis: % Decrease ####################################

# LSIL
p_LSIL<-ggplot(screen_0.6, aes(x=Year)) +
  geom_line(data=screen_0.6, aes(y=Undet_LSIL, linetype="Undetected"), colour="red",) + 
  geom_line(data=screen_0.6, aes(y=Det_LSIL, linetype="Detected"), colour="red") +
  geom_line(data=screen_0.4, aes(y=Undet_LSIL), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.4, aes(y=Det_LSIL), colour="red", alpha=0.1) +
  geom_line(data=screen_0.8, aes(y=Undet_LSIL), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.8, aes(y=Det_LSIL), colour="red", alpha=0.1) +
  geom_ribbon(aes(ymin=screen_0.4$Undet_LSIL,ymax=screen_0.8$Undet_LSIL, fill="40%-80% range"), alpha=0.3) +
  geom_ribbon(aes(ymin=screen_0.8$Det_LSIL,ymax=screen_0.4$Det_LSIL, fill="40%-80% range"), alpha=0.3) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(250000,1000000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) + 
  labs(y="LSIL Cases") + 
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_fill_manual(values=c("40%-80% range"="red"),name="",guide=guide_legend(override.aes = list(alpha=c(0.3)))) + 
  theme_classic()
p_LSIL


# HSIL
p_HSIL<-ggplot(screen_0.6, aes(x=Year)) +
  geom_line(data=screen_0.6, aes(y=Undet_HSIL, linetype="Undetected"), colour="red",) + 
  geom_line(data=screen_0.6, aes(y=Det_HSIL, linetype="Detected"), colour="red") +
  geom_line(data=screen_0.4, aes(y=Undet_HSIL), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.4, aes(y=Det_HSIL), colour="red", alpha=0.1) +
  geom_line(data=screen_0.8, aes(y=Undet_HSIL), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.8, aes(y=Det_HSIL), colour="red", alpha=0.1) +
  geom_ribbon(aes(ymin=screen_0.4$Undet_HSIL,ymax=screen_0.8$Undet_HSIL, fill="40%-80% range"), alpha=0.3) +
  geom_ribbon(aes(ymin=screen_0.8$Det_HSIL,ymax=screen_0.4$Det_HSIL, fill="40%-80% range"), alpha=0.3) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,400000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) + 
  labs(y="HSIL Cases") + 
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_fill_manual(values=c("40%-80% range"="red"),name="",guide=guide_legend(override.aes = list(alpha=c(0.3)))) + 
  theme_classic()
p_HSIL


p_Cancer<-ggplot(screen_0.6, aes(x=Year)) +
  geom_line(aes(y=Undet_Cancer, linetype="Undetected"), colour="red") + 
  geom_line(aes(y=Det_Cancer, linetype="Detected"), colour="red") +
  geom_line(aes(y=Cancer_Death, linetype="Cervical Cancer Death"), colour="red") +
  geom_line(data=screen_0.4, aes(y=Undet_Cancer), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.4, aes(y=Det_Cancer), colour="red", alpha=0.1) +
  geom_line(data=screen_0.4, aes(y=Cancer_Death), colour="red", alpha=0.1) +
  geom_line(data=screen_0.8, aes(y=Undet_Cancer), colour="red", alpha=0.1) + 
  geom_line(data=screen_0.8, aes(y=Det_Cancer), colour="red", alpha=0.1) +
  geom_line(data=screen_0.8, aes(y=Cancer_Death), colour="red", alpha=0.1) +
  geom_ribbon(aes(ymin=screen_0.4$Undet_Cancer,ymax=screen_0.8$Undet_Cancer, fill="40%-80% range"), alpha=0.3) +
  geom_ribbon(aes(ymin=screen_0.8$Det_Cancer,ymax=screen_0.4$Det_Cancer, fill="40%-80% range"), alpha=0.3) +
  geom_ribbon(aes(ymin=screen_0.8$Cancer_Death,ymax=screen_0.4$Cancer_Death, fill="40%-80% range"), alpha=0.3) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,50000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  labs(y="Cancer Cases/Deaths") +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid", "Cervical Cancer Death"="dotted"), name="") +
  scale_fill_manual(values=c("40%-80% range"="red"),name="",guide=guide_legend(override.aes = list(alpha=c(0.3)))) + 
  theme_classic()
p_Cancer



#################################### Analysis by Race ####################################

# LSIL
p_LSIL_race<-ggplot(Results_white, aes(x=Year)) +
  geom_line(data=Results_white, aes(y=Undet_LSIL,linetype="Undetected",color="non-Hispanic White",alpha="Disruption")) + 
  geom_line(data=Results_white, aes(y=Det_LSIL,linetype="Detected",color="non-Hispanic White",alpha="Disruption")) +
  geom_line(data=Results_black, aes(y=Undet_LSIL,linetype="Undetected",color="non-Hispanic Black",alpha="Disruption")) + 
  geom_line(data=Results_black, aes(y=Det_LSIL,linetype="Detected",color="non-Hispanic Black",alpha="Disruption")) +
  geom_line(data=Results_hisp, aes(y=Undet_LSIL,linetype="Undetected",color="Hispanic",alpha="Disruption")) + 
  geom_line(data=Results_hisp, aes(y=Det_LSIL,linetype="Detected",color="Hispanic",alpha="Disruption")) +
  geom_line(data=Results_white_noCov, aes(y=Undet_LSIL,color="non-Hispanic White",alpha="No disruption")) + 
  geom_line(data=Results_white_noCov, aes(y=Det_LSIL,color="non-Hispanic White",alpha="No disruption")) +
  geom_line(data=Results_black_noCov, aes(y=Undet_LSIL,color="non-Hispanic Black",alpha="No disruption")) + 
  geom_line(data=Results_black_noCov, aes(y=Det_LSIL,color="non-Hispanic Black",alpha="No disruption")) +
  geom_line(data=Results_hisp_noCov, aes(y=Undet_LSIL,color="Hispanic",alpha="No disruption")) + 
  geom_line(data=Results_hisp_noCov, aes(y=Det_LSIL,color="Hispanic",alpha="No disruption")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,600000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_color_manual(values=c("non-Hispanic White"="red","non-Hispanic Black"="blue","Hispanic"="green"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="black", alpha=c(0.05,1)))) +
  labs(y="LSIL Cases") +
  theme_classic()
p_LSIL_race


# HSIL
p_HSIL_race<-ggplot(Results_white, aes(x=Year)) +
  geom_line(data=Results_white, aes(y=Undet_HSIL,linetype="Undetected",color="non-Hispanic White",alpha="Disruption")) + 
  geom_line(data=Results_white, aes(y=Det_HSIL,linetype="Detected",color="non-Hispanic White",alpha="Disruption")) +
  geom_line(data=Results_black, aes(y=Undet_HSIL,linetype="Undetected",color="non-Hispanic Black",alpha="Disruption")) + 
  geom_line(data=Results_black, aes(y=Det_HSIL,linetype="Detected",color="non-Hispanic Black",alpha="Disruption")) +
  geom_line(data=Results_hisp, aes(y=Undet_HSIL,linetype="Undetected",color="Hispanic",alpha="Disruption")) + 
  geom_line(data=Results_hisp, aes(y=Det_HSIL,linetype="Detected",color="Hispanic",alpha="Disruption")) +
  geom_line(data=Results_white_noCov, aes(y=Undet_HSIL,color="non-Hispanic White",alpha="No disruption")) + 
  geom_line(data=Results_white_noCov, aes(y=Det_HSIL,color="non-Hispanic White",alpha="No disruption")) +
  geom_line(data=Results_black_noCov, aes(y=Undet_HSIL,color="non-Hispanic Black",alpha="No disruption")) + 
  geom_line(data=Results_black_noCov, aes(y=Det_HSIL,color="non-Hispanic Black",alpha="No disruption")) +
  geom_line(data=Results_hisp_noCov, aes(y=Undet_HSIL,color="Hispanic",alpha="No disruption")) + 
  geom_line(data=Results_hisp_noCov, aes(y=Det_HSIL,color="Hispanic",alpha="No disruption")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,200000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid"), name="") +
  scale_color_manual(values=c("non-Hispanic White"="red","non-Hispanic Black"="blue","Hispanic"="green"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="black", alpha=c(0.05,1)))) +
  labs(y="HSIL Cases") +
  theme_classic()
p_HSIL_race


p_Cancer_race<-ggplot(result_tot, aes(x=Year)) +
  geom_line(data=Results_white, aes(y=Undet_Cancer,color="non-Hispanic White",linetype="Undetected",alpha="Disruption")) + 
  geom_line(data=Results_white, aes(y=Det_Cancer,color="non-Hispanic White",linetype="Detected",alpha="Disruption")) +
  geom_line(data=Results_white, aes(y=Cancer_Death,color="non-Hispanic White",linetype="Cervical Cancer Death",alpha="Disruption")) +
  geom_line(data=Results_black, aes(y=Undet_Cancer,color="non-Hispanic Black",linetype="Undetected",alpha="Disruption")) + 
  geom_line(data=Results_black, aes(y=Det_Cancer,color="non-Hispanic Black",linetype="Detected",alpha="Disruption")) +
  geom_line(data=Results_black, aes(y=Cancer_Death,color="non-Hispanic Black",linetype="Cervical Cancer Death",alpha="Disruption")) +
  geom_line(data=Results_hisp, aes(y=Undet_Cancer,color="Hispanic",linetype="Undetected",alpha="Disruption")) + 
  geom_line(data=Results_hisp, aes(y=Det_Cancer,color="Hispanic",linetype="Detected",alpha="Disruption")) +
  geom_line(data=Results_hisp, aes(y=Cancer_Death,color="Hispanic",linetype="Cervical Cancer Death",alpha="Disruption")) +
  geom_line(data=Results_white_noCov, aes(y=Undet_Cancer,color="non-Hispanic White",alpha="No disruption")) + 
  geom_line(data=Results_white_noCov, aes(y=Det_Cancer,color="non-Hispanic White",alpha="No disruption")) +
  geom_line(data=Results_white_noCov, aes(y=Cancer_Death,color="non-Hispanic White",alpha="No disruption")) +
  geom_line(data=Results_black_noCov, aes(y=Undet_Cancer,color="non-Hispanic Black",alpha="No disruption")) + 
  geom_line(data=Results_black_noCov, aes(y=Det_Cancer,color="non-Hispanic Black",alpha="No disruption")) +
  geom_line(data=Results_black_noCov, aes(y=Cancer_Death,color="non-Hispanic Black",alpha="No disruption")) +
  geom_line(data=Results_hisp_noCov, aes(y=Undet_Cancer,color="Hispanic",alpha="No disruption")) + 
  geom_line(data=Results_hisp_noCov, aes(y=Det_Cancer,color="Hispanic",alpha="No disruption")) +
  geom_line(data=Results_hisp_noCov, aes(y=Cancer_Death,color="Hispanic",alpha="No disruption")) +
  coord_cartesian(
    xlim = c(2000,2050),
    ylim = c(0,35000)) +
  geom_vline(xintercept = 2008,linetype="dashed",alpha=0.3) +
  geom_vline(xintercept = 2020,linetype="dashed",alpha=0.3) +
  scale_linetype_manual(values=c("Undetected"="dashed", "Detected"="solid","Cervical Cancer Death"="dotted"), name="") +
  scale_color_manual(values=c("non-Hispanic White"="red","non-Hispanic Black"="blue","Hispanic"="green"), name="") +
  scale_alpha_manual(values=c("No disruption"=0.5, "Disruption"=1), name="", guide=guide_legend(override.aes = list(color="black", alpha=c(0.03,1)))) +
  labs(y="Cancer Cases/Deaths") +
  theme_classic()
p_Cancer_race




