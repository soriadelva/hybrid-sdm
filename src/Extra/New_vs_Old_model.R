#--------------------------------------------------
# Compare fits of the AICc best model with our original model
#--------------------------------------------------

#Initialize dataframe
modelfits<-data.frame(Temperature=  seq(0,50 , length.out = 1500))

#Initialize trait dataframes
growth_data<-filter(Alldata, Trait=="Growth")
pam_data<-filter(Alldata, Trait=="FvFm")
germination_data<-filter(Alldata, Trait=="Germination")


#--------------------------------------------------
# Create fits using the best model
#--------------------------------------------------
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

#--------------------------------------------------
# Create fits using the old model
#--------------------------------------------------
growth_oldfit<- nlsLM(formula12, 
                   data=growth_data, 
                   start = list(a = 10, b = 0.5, c=15, d=0.7, e=28), 
                   control = list(maxiter = 100))

pam_oldfit<- nlsLM(formula12, 
                data=pam_data, 
                start = list(a = 10, b = 0.5, c=15, d=0.7, e=28))

germination_oldfit<- nlsLM(formula12, 
                        data=germination_data, 
                        start = list(a = 100, b = 0.5, c=15, d=1, e=20))


#Add fits to dataframe
modelfits$Growth<- predict(growth_fit,modelfits)
modelfits$FvFm<- predict(pam_fit,modelfits)
modelfits$Germination<- predict(germination_fit,modelfits)
modelfits$oldGrowth<- predict(growth_oldfit,modelfits)
modelfits$oldFvFm<- predict(pam_oldfit,modelfits)
modelfits$oldGermination<- predict(germination_oldfit,modelfits)

#Convert negative growth rates to NA; FvFm and Germination are never negative
modelfits <- modelfits %>%
  mutate(Growth = ifelse(Growth < 0, NA_real_, Growth),
         oldGrowth = ifelse(oldGrowth < 0, NA_real_, oldGrowth))

#----------------------------------------------
#------ Generate plots -------
#----------------------------------------------
growth<-ggplot() +
  geom_point(data=growth_data, color="#163f6bff", alpha=0.3, aes(x = Temperature, y = Value))+
  #geom_line(data = Predict_growth, color="#ef42f5", size=1) +
  geom_line(data = modelfits, color="#163f6bff", size=1,   aes(x = Temperature, y = Growth))+
  geom_line(data = modelfits, color="red", size=1,  aes(x = Temperature, y = oldGrowth))+
  geom_point(data = meangrowth, shape = 21, fill = "#163f6bff",size=4,   aes(x = Temperature, y = Growth) )+
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
  geom_line(data = modelfits, color="red", size=1,  aes(x = Temperature, y = oldFvFm))+
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
  geom_line(data = modelfits, color="red", size=1,  aes(x = Temperature, y = oldGermination))+
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
