#--------------------------------------------------
#--------------- Start clean  ----------------
#--------------------------------------------------
rm(list=(ls()[ls()!="v"]))
gc(reset=TRUE)


#--------------------------------------------------
#--------------- Load packages ----------------
#--------------------------------------------------
library(dplyr)
library(ggplot2)

#---------------------------------------------
#-------Define folders and model types--------
#---------------------------------------------
effect_folder<-file.path("results", "final_brt_predictions", "response_curves")
Model_types<-c("Correlative", "FvFm", "Germination", "Growth")

#---------------------------------------------
#-------Load data and merge in one df---------
#---------------------------------------------
effects<-data.frame()
for(Model_type in Model_types){
  
  effect_real<-read.csv2(file.path(effect_folder, paste0(Model_type, "_response_curves.csv")))%>%
    mutate(Model_type = Model_type)
  
  if(nrow(effects)==0){
    effects<-effect_real
  }else{
    effects<-bind_rows(effects, effect_real) 
  }
}

  
#---------------------------------------------
#-----------------Create plots----------------
#---------------------------------------------

light<-filter(effects, Predictor_name=="Mean_light")%>%
  ggplot() +
  geom_line(aes(group=Model_type,color=Model_type, x = Variable_for_plotting, y = Effect), size=0.8) +
  #scale_y_continuous(breaks=c(0,20,40,60,80,100), limits= c(-8.75, 110) ) +
  scale_color_manual(values=c( "grey70","#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("Correlative","FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  labs (x=  "Mean light",
        y= "Effect")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = rel(1)),
        axis.title= element_text(colour = "black", size = rel(1.2)),
        legend.text= element_text(colour = "black", size = rel(1.2)),
        legend.title= element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        #axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(a)",    
            x = -Inf,  # right edge
            y = Inf,  # top edge
            hjust = -0.5,  # slightly inside the plot
            vjust = 1.5, 
            size = 5, 
            colour = "black", 
            fontface="bold")



salinity<-filter(effects, Predictor_name=="Min_Salinity")%>%
  ggplot() +
  geom_line(aes(group=Model_type,color=Model_type, x = Variable_for_plotting, y = Effect), size=0.8) +
  #scale_y_continuous(breaks=c(0,20,40,60,80,100), limits= c(-8.75, 110) ) +
  scale_color_manual(values=c( "grey70","#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("Correlative","FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  labs (x=  "Minimum salinity (PSS)",
        y= "Effect")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = rel(1)),
        axis.title= element_text(colour = "black", size = rel(1.2)),
        legend.text= element_text(colour = "black", size = rel(1.2)),
        legend.title= element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        #axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(b)",    
            x = -Inf,  # right edge
            y = Inf,  # top edge
            hjust = -0.5,  # slightly inside the plot
            vjust = 1.5, 
            size = 5, 
            colour = "black", 
            fontface="bold")


maxtemp<-filter(effects, Predictor_name=="Max_temperature")%>%
  ggplot() +
  geom_line(aes(group=Model_type,color=Model_type, x = Variable_for_plotting, y = Effect), size=0.8) +
  scale_y_continuous(breaks=c(0,0.3,0.6,0.9), limits= c(0, 1) ) +
  scale_color_manual(values=c( "grey70","#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("Correlative","FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  labs (x=  "Maximum temperature (°C)",
        y= "Effect")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = rel(1)),
        axis.title= element_text(colour = "black", size = rel(1.2)),
        legend.text= element_text(colour = "black", size = rel(1.2)),
        legend.title= element_blank(),
        axis.title.y = element_text(vjust = 0.5),
       # axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(d)",    
            x = -Inf,  # right edge
            y = Inf,  # top edge
            hjust = -0.5,  # slightly inside the plot
            vjust = 1.5, 
            size = 5, 
            colour = "black", 
            fontface="bold")


mintemp<-filter(effects, Predictor_name=="Min_temperature")%>%
  ggplot() +
  geom_line(aes(group=Model_type,color=Model_type, x = Variable_for_plotting, y = Effect), size=0.8) +
  #scale_y_continuous(breaks=c(0,20,40,60,80,100), limits= c(-8.75, 110) ) +
  scale_color_manual(values=c( "grey70","#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("Correlative","FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  labs (x=  "Minimum temperature (°C)",
        y= "Effect")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "top",
        axis.text = element_text(colour = "black", size = rel(1)),
        axis.title= element_text(colour = "black", size = rel(1.2)),
        legend.text= element_text(colour = "black", size = rel(1.2)),
        legend.title= element_blank(),
        axis.title.y = element_text(vjust = 0.5),
       # axis.title.x = element_text(vjust = -0.5),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white", colour="white", size=1.2),
        panel.border = element_rect(colour="black")) +
  annotate( "text", label = "(c)",    
            x = -Inf,  # right edge
            y = Inf,  # top edge
            hjust = -0.5,  # slightly inside the plot
            vjust = 1.5, 
            size = 5, 
            colour = "black", 
            fontface="bold")

mintemp
#-------------------------------
#Create combined graphs
#------------------------------

#--------------------------Create legends------------------------
legend_horizontal_comb <- ggpubr::get_legend(
  mintemp +
    guides(
      fill = guide_legend(
        title = NULL,          # position of the title
        # horizontal justification
        direction = "horizontal",   # make legend horizontal
        keywidth = unit(0.4 ,"cm"),   # width of each key
        keyheight = unit(0.4, "cm") # height of each key
      )
    ) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 12)
    )
)

legend_horizontal_comb <- ggpubr::as_ggplot(legend_horizontal_comb)


#--------------------------Create plots------------------------
final_plot_combined<-legend_horizontal_comb/wrap_plots(
  light,
  salinity, 
  mintemp,
  maxtemp,
  ncol = 2
)&
  theme(
    legend.position = "none")

final_plot_combined<- final_plot_combined+ plot_layout(heights = c( 0.1, 1))
final_plot_combined #width 800, height 631


#---------------------------------------------
#---------Plot real effect plots---------------
#---------------------------------------------
ggplot() +
  geom_line(data=effects,aes(color=Model_type, x = Variable, y = Effect), size=0.8) +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1), limits= c(0, 1) ) +
  scale_color_manual(values=c( "grey70","#f0a144ff","#48a8b3ff","#163f6bff"),
                     labels=c("Correlative","FvFm" = expression(F[v]/F[m]), "Germination", "Growth"))+
  facet_wrap(~ Model_type + Predictor_name, scales = "free") +
  labs (y= "effect")  +
  theme_test()+
  theme(panel.background = element_blank(),
        legend.position = "none",
        axis.text = element_text(colour = "black", size = rel(1)),
        axis.title= element_text(colour = "black", size = rel(1.2)),
        legend.text= element_text(colour = "black", size = rel(1.1)),
        legend.title= element_blank(),
        axis.title.y = element_text(vjust = 0.5),
        axis.title.x = element_text(vjust = -1),
        legend.key = element_blank(),
        strip.background = element_rect(fill="white",color="black", size=1.2),
        panel.border = element_rect(colour="black")) 
