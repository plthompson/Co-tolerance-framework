library(vegan)
library(MASS)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(data.table)

load("Workspace/Multistress.RData")
Output_means<-Output_means%>%
  ungroup()%>%
  mutate(CoTolerance=replace(CoTolerance,CoTolerance=="Negative" & Stress_type =="-,+","hold"),
         CoTolerance=replace(CoTolerance,CoTolerance=="Positive" & Stress_type =="-,+","Negative"),
         CoTolerance=replace(CoTolerance,is.na(CoTolerance) & Stress_type =="-,+","Positive"))

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
  xlab(expression(paste("Stressor magnitude (",italic("m"),")",sep="")))+
  ylab("Species richness relative to control community")
ggsave("./Figures/Null model - Figure 2.pdf",width = 8, height=6)


ggplot(filter(Output_means,Interactions=="No interactions",Response=="Biomass"),aes(x=Stress,y=Change_mean,color=Null_model, fill=Null_model,group=Null_model, linetype=Null_model))+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Stress_type~CoTolerance,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype_manual(values=LineV,name="")+
  theme_bw()+
  removeGrid()+
  xlab(expression(paste("Stressor magnitude (",italic("m"),")",sep="")))+
  ylab("Community biomass relative to control community")
ggsave("./Figures/Null model - Figure 4.pdf",width = 8, height=6)

