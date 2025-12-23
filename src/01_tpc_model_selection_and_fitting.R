#--------------------------------------------------
#--------------- Clean up workspace----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Load packages ----------------
#--------------------------------------------------
library(plyr)
library(dplyr)
library(nlstools)
library(minpack.lm)
library(readxl)
library(ggplot2)
library(patchwork)
library(zen4R)
library(readr)
library(AICcmodavg)
library(patchwork)
library(tidyr)


#--------------------------------------------------
#-------------Source helper functions -------------
#--------------------------------------------------
source(file.path("src", "helpers", "TPC_model_selection_helpers.R"))


#--------------------------------------------------
#--------------- Create folders  ----------------
#--------------------------------------------------
input_folder<-file.path("data", "Input", "Physiology")
output_folder<- file.path("results", "model_selection")
tpc_output_folder<- file.path("results", "tpc_metrics")

for(folder in c(input_folder, output_folder, tpc_output_folder)){
if(!dir.exists(folder)) dir.create(folder, recursive=TRUE)
}


#--------------------------------------------------
#--------------- Load Files ----------------
#--------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.17865936", 
                                    path=input_folder,
                                    files=c("Dictyota_growth.csv",
                                            "Dictyota_FvFm.csv",
                                            "Dictyota_germination.csv"),
                                    quiet=FALSE)

Growth <- read.csv2(file.path(input_folder,"Dictyota_growth.csv"))
FvFm <- read.csv2(file.path(input_folder,"Dictyota_FvFm.csv"))
Germination <- read.csv2(file.path(input_folder,"Dictyota_germination.csv"))


#--------------------------------------------------
#--------------- Combine files ----------------
#--------------------------------------------------
FvFm <- FvFm %>%
  dplyr::rename(Value = FvFm,
         Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "FvFm")%>%
  dplyr::select(Value, Trait, Temperature)

Growth <- Growth %>%
  dplyr::rename(Value = Growth_rate,
         Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "Growth")%>%
  dplyr::select(Value, Trait, Temperature)

Germination <- Germination %>%
  dplyr::rename(Value = Germination_rate_percent,
         Temperature = Temperature_C) %>%
  dplyr::mutate(Trait = "Germination")%>%
  dplyr::select(Value, Trait, Temperature)
  
Alldata<-bind_rows(FvFm, Growth, Germination)%>%
  dplyr::mutate(TemperatureK = Temperature + 273.15)

  
#------------------------------------------------------------
#------------------- Define TPC models-----------------------
#------------------------------------------------------------
#define formulas
formula1 <- as.formula(Value ~ Gmax * ((Tmax-Temperature)/(Tmax-Topt))^Betaa * exp(-Betaa*((Tmax-Temperature)/(Tmax-Topt)-1))) # Formula Blanchard et al

formula2 <- as.formula(Value ~ Gmax * exp(-1*((Temperature-Topt)/theta)^2)) # function of Drake et al

formula3 <- as.formula(Value ~ a * exp (-b/(0.001987*(Temperature + 273.15)))-c * exp(-d/(0.001987*(Temperature + 273.15))))# based on thermodynamics of chemical reactions

formula4 <- as.formula(Value~ (a*((Temperature + 273.15)/298.15)*exp((b/1.987)*(1/298.15-1/(Temperature + 273.15))))/(1+exp((c/1.987)*(1/d-1/(Temperature + 273.15))))) #Schoolfield et al., 1981, based on thermodynamics of chemical reactions

formula5 <- as.formula( Value ~ ifelse(Temperature <= Topt,
                                       Gmax * exp(-((Temperature - Topt) / (2 * a))^2),
                                       Gmax - Gmax * (((Temperature - Topt) / (Topt - Tmax))^2)))#Deutsch

formula6 <- as.formula(Value~ a*exp(-0.5*((Temperature-tref)/b)^2)) #modification of a gaussian function, but allows direct estimation of ecologically relevant parameters 

formula7 <- as.formula(Value~ a*exp(-0.5*(abs(Temperature-tref)/b)^c)) # modification of a gaussian function, but allows direct estimation of ecologically relevant parameters

formula8 <- as.formula(Value~ a*exp(c*Temperature)*(1-((Temperature-tref)/b)^2))

formula9<- as.formula(Value~ (a*(Temperature-Tmin))^2*(1-exp(b*(Temperature-Tmax)))^2)

formula10 <- as.formula(Value~ (a*(1-exp(-b*(Temperature-Tmin)))*(1-exp(-c*(Tmax-Temperature)))))

formula11 <- as.formula(Value~ Gmax*(sin(3.141593*((Temperature-Tmin)/(Tmax-Tmin))^a))^b)

formula12 <- as.formula(Value~ a/((1+exp(-b*(Temperature-c)))*(1+exp(d*(Temperature-e))))) #Paper Rodriguez et al


#------------------------------------------------------------
#------- Estimate start values of model parameters ----------
#------------------------------------------------------------    
test<-dplyr::filter(Alldata, Trait=="Growth")
preview(formula1,test, list(Gmax = 8, Tmax = 30, Topt = 22, Betaa=4), variable=3)
preview(formula2,test, list(Gmax = 8, Topt = 22,  theta=4),variable=3)
preview(formula3,test, list(a = 1.477e+23, b = 2.881e+01, c= 1.040e+26, d= 3.275e+01),variable=3)
preview(formula4,test, list(a = 10, b = 9000, c=300000, d=301),variable=3)
preview(formula5,test, list(Gmax = 8, Tmax = 30, Topt = 22, a=4),variable=3)
preview(formula6, test, list(a = 10, b = 10, tref=20),variable=3)
preview(formula7,test, list(a = 10, b = 10, tref=20, c=10),variable=3)
preview(formula8,test, list(a = 8, b = 15, c=0, tref=22),variable=3)
preview(formula9,test, list(a = 0.15, b = 1, Tmin=4, Tmax=30),variable=3)
preview(formula10,test, list(a = 10, b = 0.17, c=10, Tmin=7, Tmax=28.1),variable=3)
preview(formula11,test, list(a = 1, b = 3, Gmax=13, Tmin=4, Tmax=32),variable=3)
preview(formula12,test, list(a = 10, b = 0.5, c=15, d=0.7, e=28),variable=3)
rm(test)


#------------------------------------------------------------
#---------- Perform model fitting in loop -------------------
#------------------------------------------------------------ 

#Create a dataset with AICc and RSS values for every trait and formula
Trait <- rep(NA, 12)
Formula <- rep(c("formula1","formula2","formula3","formula4",
                 "formula5","formula6","formula7","formula8",
                 "formula9","formula10","formula11","formula12"), 3)
AICc_values <- rep(NA, 36)
Root_mean_squared_error <- rep(NA, 36)
Final_AIC <- data.frame(Trait,Formula,AICc_values,Root_mean_squared_error)

i = 1

for(trait in unique(Alldata$Trait)){
  
  trait_data <- filter(Alldata, Trait == trait) %>%
    select(Temperature, Value)
  
  #---------------------------
  #------- Fit models --------
  #---------------------------
  #Formula1 
  fit_1 <- nlsLM(formula1, 
                 data=trait_data, 
                 start = list(Gmax = 8, Tmax = 32, Topt = 22, Betaa=4), 
                 control = list(maxiter = 100))
  
  #Formula 2 
  fit_2<- nlsLM(formula2, 
                data=trait_data, 
                start = list(Gmax = 8, Topt = 22,  theta=4))
  
  # Formula 3 
  fit_3<- nlsLM(formula3, 
                data=trait_data, 
                start =list(a = 1.477e+23, b = 2.881e+01, c= 1.040e+26, d= 3.275e+01), 
                control = list(maxiter = 100)) 
  
  # Formula 4
  if(!trait=="Germination"){
  fit_4<- nlsLM(formula4, 
                data=trait_data, 
                start = list(a = 10, b = 9000, c=300000, d=301))
  }else{
    fit_4<- nlsLM(formula4, 
                  data=trait_data, 
                  start = list(a = 100, b = 9000, c=300000, d=301))
  }
  
  #Formula 5
  fit_5<- nlsLM(formula5, 
                data=trait_data, 
                start = list(Gmax = 8, Tmax = 32, Topt = 22, a=4), 
                control = list(maxiter = 100))
  
  #Formula 6
  fit_6<- nlsLM(formula6, 
                data=trait_data, 
                start = list(a = 10, b = 10, tref=20))
  
  #Formula 7
  if(!trait=="Germination"){
  fit_7<- nlsLM(formula7, 
                data=trait_data, 
                start = list(a = 10, b = 10, tref=20, c=10))
  }else{
  fit_7<- nlsLM(formula7, 
                data=trait_data, 
                start = list(a = 100, b = 10, tref=20, c=10))
  }
  
  #Formula 8
  fit_8<- nlsLM(formula8, 
                data=trait_data, 
                start = list(a = 8, b = 15, c=0, tref=22))
  
  #Formula 9
  fit_9<- nlsLM(formula9, 
                data=trait_data, 
                start = list(a = 0.15, b = 1, Tmin=4, Tmax=30), 
                control = list(maxiter = 100))
  
  #Formula 10
  fit_10<- nlsLM(formula10, 
                 data=trait_data, 
                 start =  list(a = 110, b = 0.007, c=10, Tmin=5, Tmax=28.1),
                 control = list(maxiter = 100))
  
  #Formula 11
  if(!trait=="Germination"){
  fit_11<- nlsLM(formula11, 
                 data=trait_data, 
                 start = list(a = 1, b = 3, Gmax=13, Tmin=4, Tmax=32))
  }else{
  fit_11<- nlsLM(formula11, data=trait_data, 
                 start = list(a = 1, b = 1, Gmax=100, Tmin=8, Tmax=30),
                 lower = c(a=0.1, b=0.1, Gmax=0, Tmin = 0, Tmax = 20),
                 upper = c(a=10, b=10, Gmax = 120, Tmin = 6, Tmax = 30)) 
  }
  
  #Formula 12
  if(!trait=="Germination"){
  fit_12<- nlsLM(formula12, 
                 data=trait_data, start = list(a = 10, b = 0.5, c=15, d=0.7, e=28))
  }else{
  fit_12<- nlsLM(formula12, 
                 data=trait_data, start = list(a = 100, b = 0.5, c=15, d=1, e=20), 
                 control=list(maxiter=100))
  }
  
  
  #---------------------------
  #----Store trait -----
  #---------------------------
  Final_AIC[(i:(i + 11)),1]<-trait
  
  #---------------------------
  #----Store AICc values -----
  #---------------------------
  Final_AIC[i,3]<-AICc(fit_1)
  Final_AIC[i+1,3]<-AICc(fit_2)
  Final_AIC[i+2,3]<-AICc(fit_3)
  Final_AIC[i+3,3]<-AICc(fit_4)
  Final_AIC[i+4,3]<-AICc(fit_5)
  Final_AIC[i+5,3]<-AICc(fit_6)
  Final_AIC[i+6,3]<-AICc(fit_7)
  Final_AIC[i+7,3]<-AICc(fit_8)
  Final_AIC[i+8,3]<-AICc(fit_9)
  Final_AIC[i+9,3]<-AICc(fit_10)
  Final_AIC[i+10,3]<-AICc(fit_11)
  Final_AIC[i+11,3]<-AICc(fit_12)
  
  
  #---------------------------
  #----Store RMSQ values -----
  #---------------------------
  #RMSQ: sqrt(RSS/n)
  Final_AIC[i,4]<-sqrt(deviance(fit_1)/length(trait_data$Value))
  Final_AIC[i+1,4]<-sqrt(deviance(fit_2)/length(trait_data$Value))
  Final_AIC[i+2,4]<-sqrt(deviance(fit_3)/length(trait_data$Value))
  Final_AIC[i+3,4]<-sqrt(deviance(fit_4)/length(trait_data$Value))
  Final_AIC[i+4,4]<-sqrt(deviance(fit_5)/length(trait_data$Value))
  Final_AIC[i+5,4]<-sqrt(deviance(fit_6)/length(trait_data$Value))
  Final_AIC[i+6,4]<-sqrt(deviance(fit_7)/length(trait_data$Value))
  Final_AIC[i+7,4]<-sqrt(deviance(fit_8)/length(trait_data$Value))
  Final_AIC[i+8,4]<-sqrt(deviance(fit_9)/length(trait_data$Value))
  Final_AIC[i+9,4]<-sqrt(deviance(fit_10)/length(trait_data$Value))
  Final_AIC[i+10,4]<-sqrt(deviance(fit_11)/length(trait_data$Value))
  Final_AIC[i+11,4]<-sqrt(deviance(fit_12)/length(trait_data$Value))
  
  #---------------------------
  #------   Clean up    ------
  #--------------------------- 
  rm(trait_data, fit_1, fit_2, fit_3, fit_4, fit_5, fit_6, fit_7, fit_8, fit_9, fit_10, fit_11, fit_12)
  
  i=i+12
  
}


#------------------------------------------------------------
#-------- Calculate delta AICc and Akaiki weights------------
#------------------------------------------------------------ 

# Calculate delta AICc per trait
Final_AIC <- Final_AIC %>% 
  group_by(Trait) %>% 
  dplyr::mutate(AICDifference = AICc_values - min(AICc_values, na.rm = TRUE),
                AICweights = exp(-0.5 * AICDifference) / sum(exp(-0.5 * AICDifference), na.rm = TRUE)) %>%
  ungroup()

#Export file
write.csv2(Final_AIC, file.path(output_folder,"Model_selection_overview.csv"), row.names=F)


#----------------------------------------------
#- Prepare df with fitted values for plotting -
#----------------------------------------------
#Initialize dataframe
modelfits<-data.frame(Temperature=  seq(0,50 , length.out = 1500))

#Initialize trait dataframes
growth_data<-filter(Alldata, Trait=="Growth")
pam_data<-filter(Alldata, Trait=="FvFm")
germination_data<-filter(Alldata, Trait=="Germination")

#Initialize trait dataframes with mean values per temperature
meangrowth<-plyr::ddply(growth_data, .(Temperature),summarize,
                                  Growth = round(mean(Value), 6))
meanfvfm<-plyr::ddply(pam_data, .(Temperature),summarize,
                        FvFm = round(mean(Value), 6))
meangermination<-plyr::ddply(germination_data, .(Temperature),summarize,
                        Germination = round(mean(Value), 6))

#Create fits for each trait using the best model
growth_fit<- nlsLM(formula5, 
                  data=growth_data, 
                  start = list(Gmax = 8, Tmax = 32, Topt = 22, a=4), 
                  control = list(maxiter = 100))

pam_fit<- nlsLM(formula7, 
                data=pam_data, 
                start = list(a = 10, b = 10, tref=20, c=10))

germination_fit<- nlsLM(formula7, 
                        data=germination_data, 
                        start = list(a = 100, b = 10, tref=20, c=10))

#Add fits to dataframe
modelfits$Growth<- predict(growth_fit,modelfits)
modelfits$FvFm<- predict(pam_fit,modelfits)
modelfits$Germination<- predict(germination_fit,modelfits)

#Convert negative growth rates to NA; FvFm and Germination are never negative
modelfits <- modelfits %>%
  mutate(Growth = ifelse(Growth < 0, NA_real_, Growth))


#----------------------------------------------
#-- Add bootstrap confidence intervals to df --
#----------------------------------------------
# Run bootstraps for each trait
set.seed(123)
boot_growth <- boot::boot(data = growth_data,
                          statistic = bootstrap_tpc,
                          R = 1000, #Number of bootstraps
                          newdata = modelfits,
                          formula = formula5,
                          start_coef = list(Gmax = 8, Tmax = 32, Topt = 22, a = 4)
                          )

set.seed(456)
boot_fvfm <- boot::boot(data = pam_data,
                        statistic =   bootstrap_tpc,
                        R = 1000,
                        newdata = modelfits,
                        formula = formula7,
                        start_coef = list(a = 10, b = 10, tref=20, c=10)
                        )

set.seed(30) 
boot_germination <- boot::boot(data = germination_data,
                               statistic =   bootstrap_tpc,
                               R = 1000,
                               newdata = modelfits,
                               formula = formula7,
                               start_coef = list(a = 100, b = 10, tref=20, c=10)
                               )

#Get confidence intervals: In boot_trait$t rows = replicate predictions and cols = temperatures; apply 2 stands for cols
ci_growth <- apply(boot_growth$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
ci_fvfm <- apply(boot_fvfm$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
ci_germination <- apply(boot_germination$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))

# Add confidence intervals to modelfits 
modelfits$lower_growth <- ci_growth[1, ]
modelfits$upper_growth <- ci_growth[2, ]

modelfits$lower_fvfm <- ci_fvfm[1, ]
modelfits$upper_fvfm<- ci_fvfm[2, ]

modelfits$lower_germination <- ci_germination[1, ]
modelfits$upper_germination<- ci_germination[2, ]

# Set lower growth and upper growth to NA for rows where Growth is NA
modelfits <- modelfits %>%
  mutate(lower_growth = ifelse(is.na(Growth), NA_real_, lower_growth),
         upper_growth = ifelse(is.na(Growth), NA_real_, upper_growth))


#----------------------------------------------
#------ Generate plots -------
#----------------------------------------------
growth<-ggplot() +
  geom_point(data=growth_data, color="#163f6bff", alpha=0.3, aes(x = Temperature, y = Value))+
  #geom_line(data = Predict_growth, color="#ef42f5", size=1) +
  geom_line(data = modelfits, color="#163f6bff", size=1)+
  geom_ribbon(data = modelfits, aes( ymin = lower_growth, ymax = upper_growth), 
              fill = "#163f6b33", alpha = 0.2) +
  geom_point(data = meangrowth, shape = 21, fill = "#163f6bff",size=4 )+
  aes(x = Temperature, y = Growth) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24,28,32,36) , limits= c(0, 36) ) +
  scale_y_continuous(breaks=c(0,5,10,15) , limits= c(-1.4, 16) ) +
  labs (x=  expression(paste("Temperature (",degree,"C)")),
        y= expression(paste('Relative growth rate (% day '^{-1},')')))  +
  theme_test()+
theme(panel.background = element_blank(),
                        legend.position = "none",
                        axis.text = element_text(colour = "black", size = rel(1.2)),
                        axis.title= element_text(colour = "black", size = rel(1.4)),
                        legend.text= element_text(colour = "black", size = rel(1.1)),
                        legend.title= element_text(colour = "black", size = rel(1.2)),
                        axis.title.y = element_text(vjust = 0.5),
                        axis.title.x = element_text(vjust = -1),
                        legend.key = element_blank(),
                        strip.background = element_rect(fill="white", colour="white", size=1.2),
                        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(c)", x = 1, y = 15.53, size = 6, colour = "black", fontface="bold")
# annotate( "text", label = "Growth", x = 27, y = 17.9, size = 6, colour = "black", fontface="bold") 


pam<-ggplot() +
  geom_point(data=pam_data, color="#f0a144ff", alpha=0.3, aes(x = Temperature, y = Value) )+
  #geom_line(data = Predict_pam, color="#ef42f5", size=1) +
  geom_line(data = modelfits, color="#f0a144ff", size=1) +
  geom_ribbon(data = modelfits, aes(ymin = lower_fvfm, ymax = upper_fvfm), 
              fill="#f0a14450", alpha = 0.2) +
  geom_point(data = meanfvfm,shape = 21, color="black",fill = "#f0a144ff",size=4 )+
  aes(x = Temperature, y = FvFm) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24,28,32,36) , limits= c(0, 36) ) +
  scale_y_continuous(breaks=c(0.0,0.2,0.4,0.6,0.8) , limits= c(-0.074, 0.85) ) +
  labs ( x=  expression(paste("Temperature (",degree,"C)")),
         y= expression(paste(F[v]/F[m])))  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = rel(1.2)),
        axis.title= element_text(colour = "black", size = rel(1.4)),
        legend.text= element_text(colour = "black", size = rel(1.1)),
        legend.title= element_text(colour = "black", size = rel(1.2)),
        axis.title.y = element_text(vjust = 0.5),
        axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(a)", x =1, y = 0.825, size = 6, colour = "black", fontface="bold")

germ<-ggplot() +
  geom_point(data=Germination, color="#48a8b3ff", alpha=0.3,  aes(x = Temperature, y = Value))+
  #geom_line(data = Predict_germ, color="#ef42f5", size=1) +
  geom_line(data = modelfits, color="#48a8b3ff", size=1) +
  geom_ribbon(data = modelfits, aes(ymin = lower_germination, ymax = upper_germination), 
              fill = "#48a8b350", alpha = 0.2) +
  geom_point(data = meangermination, shape = 21, color="black",fill = "#48a8b3ff",size=4)+
  aes(x = Temperature, y = Germination) +
  scale_x_continuous(breaks=c(0,4,8,12,16,20,24,28,32,36) , limits= c(0, 36) ) +
  scale_y_continuous(breaks=c(0,25,50,75,100) , limits= c(-10.06, 115) ) +
  labs (x=  expression(paste("Temperature (",degree,"C)")),
        y= "Germination rate (%)")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = rel(1.2)),
        axis.title= element_text(colour = "black", size = rel(1.4)),
        legend.text= element_text(colour = "black", size = rel(1.1)),
        legend.title= element_text(colour = "black", size = rel(1.2)),
        axis.title.y = element_text(vjust = 0.5),
        axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(b)", x =1, y = 111.62, size = 6, colour = "black", fontface="bold")

pam + germ + growth  #1025-400


#------------------------------------------------------------------------
#----------------------- Create a combined plot--------------------------
#------------------------------------------------------------------------
# compute maxima used to scale each trait
max_growth <- max(modelfits$Growth, na.rm = TRUE)
max_fvfm <- max(modelfits$FvFm, na.rm = TRUE)
max_germination <- max(modelfits$Germination, na.rm = TRUE)

# create relative trait values and scaled confidence limits
Relative_rates <- modelfits %>%
  transmute(
    Temperature = Temperature,
    Growth       = (Growth / max_growth) * 100,
    lower_growth = (lower_growth / max_growth) * 100,
    upper_growth = (upper_growth / max_growth) * 100,
    FvFm         = (FvFm / max_fvfm) * 100,
    lower_fvfm   = (lower_fvfm / max_fvfm) * 100,
    upper_fvfm   = (upper_fvfm / max_fvfm) * 100,
    Germination  = (Germination / max_germination) * 100,
    lower_germination = (lower_germination / max_germination) * 100,
    upper_germination = (upper_germination / max_germination) * 100
  )

Relative_fits <- Relative_rates%>%
  pivot_longer(cols = c(Growth, FvFm,Germination),  
               names_to = "Trait",          
               values_to = "Relative_value")%>%
  select(Temperature,Trait, Relative_value)

Relative_rates<-Relative_rates%>%
  select(Temperature, lower_growth, upper_growth, lower_fvfm, upper_fvfm, lower_germination, upper_germination)
    
    
all<-ggplot() +
  geom_line(data = Relative_fits, aes(group=Trait, color=Trait, x = Temperature, y = Relative_value), size=0.8) +
  geom_ribbon(data= Relative_rates, aes(ymin=lower_growth, ymax=upper_growth, x = Temperature), fill="#163f6b33", alpha = 0.2) +
  geom_ribbon(data= Relative_rates, aes(ymin=lower_fvfm, ymax=upper_fvfm, x = Temperature), fill="#f0a14450", alpha = 0.2) +
  geom_ribbon(data= Relative_rates, aes(ymin=lower_germination, ymax=upper_germination, x = Temperature), fill="#48a8b350", alpha = 0.2)  +
  scale_x_continuous(breaks=c(0, 4,8,12,16,20,24,28,32, 36) , limits= c(0, 36) ) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100), limits= c(-8.75, 110) ) +
  scale_color_manual(values=c("#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  labs (x=  expression(paste("Temperature (",degree,"C)")),
        y= "% of maximum rate")  +
  theme_test()+
  theme(panel.background = element_blank(),
                  legend.position = "top",
                  axis.text = element_text(colour = "black", size = rel(1.2)),
                  axis.title= element_text(colour = "black", size = rel(1.4)),
                  legend.text= element_text(colour = "black", size = rel(1.1)),
                  legend.title= element_blank(),
                  axis.title.y = element_text(vjust = 0.5),
                  axis.title.x = element_text(vjust = -1),
                  legend.key = element_blank(),
                  strip.background = element_rect(fill="white", colour="white", size=1.2),
                  panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(d)", x =1, y = 106.79, size = 6, colour = "black", fontface="bold")

all

(pam+germ)/(growth+all) # 1000 800; Inkscape: Figure to 1750 width (300dpi), white edge to 1961 width


#----------------------------------------------------------------------
# 95% CI for performance breadth and upper and lower performance limits
#----------------------------------------------------------------------
# Apply function to each row (bootstrap replicate; 1)
growth_metrics <- t(apply(boot_growth$t, 1, calc_thermal_metrics, temperatures = modelfits$Temperature))
fvfm_metrics <- t(apply(boot_fvfm$t, 1, calc_thermal_metrics, temperatures = modelfits$Temperature))
germination_metrics <- t(apply(boot_germination$t, 1, calc_thermal_metrics, temperatures = modelfits$Temperature))

# Calculate 95% confidence intervals for each metric
ci_lower_growth <- quantile(growth_metrics[, "lower"], c(0.025, 0.975), na.rm = TRUE)
ci_upper_growth <- quantile(growth_metrics[, "upper"], c(0.025, 0.975), na.rm = TRUE)
ci_breadth_growth <- quantile(growth_metrics[, "breadth"], c(0.025, 0.975), na.rm = TRUE)

ci_lower_fvfm <- quantile(fvfm_metrics[, "lower"], c(0.025, 0.975), na.rm = TRUE)
ci_upper_fvfm <- quantile(fvfm_metrics[, "upper"], c(0.025, 0.975), na.rm = TRUE)
ci_breadth_fvfm <- quantile(fvfm_metrics[, "breadth"], c(0.025, 0.975), na.rm = TRUE)

ci_lower_germination <- quantile(germination_metrics[, "lower"], c(0.025, 0.975), na.rm = TRUE)
ci_upper_germination <- quantile(germination_metrics[, "upper"], c(0.025, 0.975), na.rm = TRUE)
ci_breadth_germination <- quantile(germination_metrics[, "breadth"], c(0.025, 0.975), na.rm = TRUE)

# Calculate values of actual fit 
growth_idx <- which(modelfits$Growth >= (0.8 * max_growth))
fvfm_idx <- which(modelfits$FvFm>= (0.8 * max_fvfm))
germination_idx <- which(modelfits$Germination >= (0.8 * max_germination))

lower_growth <- modelfits$Temperature[min(growth_idx)]
upper_growth <- modelfits$Temperature[max(growth_idx)]
breadth_growth <- upper_growth - lower_growth

lower_fvfm <- modelfits$Temperature[min(fvfm_idx)]
upper_fvfm <- modelfits$Temperature[max(fvfm_idx)]
breadth_fvfm <- upper_fvfm - lower_fvfm

lower_germination <- modelfits$Temperature[min(germination_idx)]
upper_germination <- modelfits$Temperature[max(germination_idx)]
breadth_germination <- upper_germination - lower_germination

#Store in dataframes per metric
lower_performance_temp<- data.frame(
  Trait=c("Growth", "FvFm", "Germination"),
  lower_CI= c(ci_lower_growth[1], ci_lower_fvfm[1], ci_lower_germination[1]),
  upper_CI=c(ci_lower_growth[2], ci_lower_fvfm[2], ci_lower_germination[2]),
  Value = c(lower_growth, lower_fvfm, lower_germination)
)

upper_performance_temp<- data.frame(
  Trait=c("Growth", "FvFm", "Germination"),
  lower_CI= c(ci_upper_growth[1], ci_upper_fvfm[1], ci_upper_germination[1]),
  upper_CI=c(ci_upper_growth[2], ci_upper_fvfm[2], ci_upper_germination[2]),
  Value = c(upper_growth, upper_fvfm, upper_germination)
)

performance_breadth<- data.frame(
  Trait=c("Growth", "FvFm", "Germination"),
  lower_CI= c(ci_breadth_growth[1], ci_breadth_fvfm[1], ci_breadth_germination[1]),
  upper_CI=c(ci_breadth_growth[2], ci_breadth_fvfm[2], ci_breadth_germination[2]),
  Value = c(breadth_growth, breadth_fvfm, breadth_germination)
)


#--------------------------------------------------------------------------
# --------------------- Visualize thermal metrics -------------------------
#--------------------------------------------------------------------------
#Lower performance limit
lp<-ggplot() +
  geom_errorbar(data=lower_performance_temp,  aes(x=Trait, ymin=lower_CI,ymax=upper_CI), color="black", position=position_dodge(0.5), width=0.2) + 
  geom_point(data=lower_performance_temp,  aes(x=Trait, y=Value, fill=Trait), size = 3.4,shape=21, color="black", position=position_dodge(0.5)) + 
  scale_fill_manual(values=c("#f0a144ff","#48a8b3ff","#163f6bff"))+
  scale_y_continuous( limits= c(5, 20) ) +
  labs ( x=  "Trait",
         y=expression("Lower performance limit" * " "* (degree~C)))  +
  theme_test()+
  theme(panel.background = element_blank(),
                legend.position = "none",
                axis.text = element_text(colour = "black", size = rel(1.2)),
                axis.title= element_text(colour = "black", size = rel(1.4)),
                legend.text= element_text(colour = "black", size = rel(1.1)),
                legend.title= element_text(colour = "black", size = rel(1.2)),
                axis.title.y = element_text(vjust = 2),
                legend.key = element_blank(),
                axis.title.x = element_blank(),
                strip.background = element_rect(fill="white", colour="white", size=1.2),
                panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(b)", x =0.7, y = 19.8, size = 6, colour = "black", fontface="bold")+
  # annotate( "text", label = "Growth", x = 27, y = 17.9, size = 6, colour = "black", fontface="bold") 
  scale_x_discrete(labels=c(expression(F[v]/F[m]), "Germination","Growth"))

#Upper performance limit
up<-ggplot() +
  geom_errorbar(data=upper_performance_temp,  aes(x=Trait, ymin=lower_CI,ymax=upper_CI), color="black", position=position_dodge(0.5), width=0.2) + 
  geom_point(data=upper_performance_temp,  aes(x=Trait, y=Value, fill=Trait), size = 3.4,shape=21, color="black", position=position_dodge(0.5)) + 
  scale_fill_manual(values=c("#f0a144ff","#48a8b3ff","#163f6bff"))+
  scale_y_continuous( limits= c(21, 30) ) +
  labs ( x=  "Trait",
         y=expression("Upper performance limit" * " "* (degree~C)))  +
  theme_test()+
  theme(panel.background = element_blank(),
                legend.position = "none",
                axis.text = element_text(colour = "black", size = rel(1.2)),
                axis.title= element_text(colour = "black", size = rel(1.4)),
                legend.text= element_text(colour = "black", size = rel(1.1)),
                legend.title= element_text(colour = "black", size = rel(1.2)),
                axis.title.y = element_text(vjust = 2),
                legend.key = element_blank(),
                axis.title.x = element_blank(),
                strip.background = element_rect(fill="white", colour="white", size=1.2),
                panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(c)", x =0.7, y = 29.8, size = 6, colour = "black", fontface="bold")+
  # annotate( "text", label = "Growth", x = 27, y = 17.9, size = 6, colour = "black", fontface="bold") +
  scale_x_discrete(labels=c(expression(F[v]/F[m]), "Germination","Growth"))

#Performance breadth
pb<-ggplot() +
  geom_errorbar(data=performance_breadth,  aes(x=Trait, ymin=lower_CI, ymax=upper_CI), color="black", width=0.2) + 
  geom_point(data=performance_breadth,  aes(x=Trait, y=Value, fill=Trait), size = 3.4,shape=21, color="black") + 
  scale_fill_manual(values=c("#f0a144ff","#48a8b3ff","#163f6bff"))+
  scale_y_continuous( limits= c(5, 23) ) +
  labs (  x=  "Trait",
          y=expression("Performance breadth" * " "* (degree~C)))  +
  theme_test()+
  theme(panel.background = element_blank(),
               plot.title= element_text(hjust=0.5, size= rel(1.4)),
               legend.position = "none",
               axis.text = element_text(colour = "black", size = rel(1.2)),
               axis.title= element_text(colour = "black", size = rel(1.4)),
               legend.text= element_text(colour = "black", size = rel(1.1)),
               legend.title= element_text(colour = "black", size = rel(1.2)),
               axis.title.y = element_text(vjust = 2),
               legend.key = element_blank(),
               axis.title.x = element_blank(),
               strip.background = element_rect(fill="white", colour="white", size=1.2),
               panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(a)", x =0.7, y = 22.8, size = 6, colour = "black", fontface="bold")+
  scale_x_discrete(labels=c(expression(F[v]/F[m]), "Germination","Growth"))

#Combine plots
(pb+ lp + up) #1030-350


#--------------------------------------------------------------------------
# --------------------- Export datasets -------------------------
#--------------------------------------------------------------------------
write.csv2(lower_performance_temp, file.path(tpc_output_folder,"Lower_performance_limit_CI.csv"), row.names=F)
write.csv2(upper_performance_temp, file.path(tpc_output_folder,"Upper_performance_limit_CI.csv"), row.names=F)
write.csv2(performance_breadth, file.path(tpc_output_folder,"Performance_breadth_CI.csv"), row.names=F)
