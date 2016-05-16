library(vegan)
library(MASS)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)

load("Workspace/Multistress.RData")
Output_means$Response<-factor(Output_means$Response,levels=c("Species richness","Biomass", "Composition"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,levels=c("A_only","B_only","Actual","Additive","Multiplicative","Species_specific"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,labels=c("A","B","AB","Additive","Multiplicative","Compositional"))

ColV<-c("grey50","grey80","black",brewer.pal(3,"Set1"))
LineV<-c(4,2,1,2,3,4)
ggplot(filter(Output_means,Interactions=="No interactions",Response=="Species richness"),aes(x=Stress,y=Change_mean,color=Null_model, fill=Null_model,group=Null_model, linetype=Null_model))+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Stress_type~CoTolerance,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype_manual(values=LineV,name="")+
  theme_bw()+
  removeGrid()+
  ylab("Species richness change from control")
ggsave("./Figures/Null model - Richness.pdf",width = 11, height=8.5)


ggplot(filter(Output_means,Interactions=="No interactions",Response=="Biomass"),aes(x=Stress,y=Change_mean,color=Null_model, fill=Null_model,group=Null_model, linetype=Null_model))+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Stress_type~CoTolerance,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype_manual(values=LineV,name="")+
  theme_bw()+
  removeGrid()+
  ylab("Biomass change from control")
ggsave("./Figures/Null model - Biomass.pdf",width = 11, height=8.5)


ggplot(filter(Output_means,Interactions=="No interactions",Response=="Composition"),aes(x=Stress,y=Change_mean,color=Null_model, fill=Null_model,group=Null_model, linetype=Null_model))+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Stress_type~CoTolerance,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype_manual(values=LineV,name="")+
  theme_bw()+
  removeGrid()+
  ylab("Compositional change from control")
ggsave("./Figures/Null model - Composition.pdf",width = 11, height=8.5)
