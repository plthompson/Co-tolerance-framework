require(tidyr)
require(dplyr)
require(ggplot2)
require(vegan)

null_compare<-function(Control,A,B,AB){
  output<-data.frame(Response=c("Biomass","Species richness","Composition"),Actual=NA,Compositional=NA,Additive=NA,Multiplicative=NA)
  
  output$Actual[1]<-(sum(AB)-sum(Control))/sum(Control)
  output$Actual[2]<-(sum(AB>0)-sum(Control>0))/sum(Control>0)
  output$Actual[3]<-vegdist(rbind(Control,AB))
  
  AB_predict<-Control+(A-Control)+(B-Control)
  AB_predict[AB_predict<0]<-0
  output$Compositional[1]<-(sum(AB_predict)-sum(Control))/sum(Control)
  output$Compositional[2]<-(sum(AB_predict>0)-sum(Control>0))/sum(Control>0)
  output$Compositional[3]<-vegdist(rbind(Control,AB_predict))
  
  pA<-(sum(A)-sum(Control))/sum(Control)
  pB<-(sum(B)-sum(Control))/sum(Control)
  output$Additive[1]<-pA+pB
  output$Multiplicative[1]<-pA+pB+pA*pB
  
  pA<-(sum(A>0)-sum(Control>0))/sum(Control>0)
  pB<-(sum(B>0)-sum(Control>0))/sum(Control>0)
  output$Additive[2]<-pA+pB
  output$Multiplicative[2]<-pA+pB+pA*pB
  
  pA<-vegdist(rbind(A,Control))
  pB<-vegdist(rbind(B,Control))
  output$Additive[3]<-pA+pB
  output$Multiplicative[3]<-pA+pB-pA*pB
  
  output<-gather(output,key = Null_model,value=Change,Actual:Multiplicative)
  output<-output %>%
    group_by(Response) %>%
    mutate(Difference=abs(Change[Null_model=="Actual"])-abs(Change))%>%
    mutate(Reversal = Change>0 & Change[Null_model=="Actual"]<0 | Change<0 & Change[Null_model=="Actual"]>0)
  
  return(output)
}

#MacClennan Experiment####
MM_data<-read.csv("./Experimental data/McClennan experiment.csv")
head(MM_data)

MM_wide<-spread(MM_data[,-c(5:6)],Species,Biomass,fill = 0)

for(lake in unique(MM_wide$Lake)){
  Lake_data<-filter(MM_wide,Lake==lake)
  
  Control<-Lake_data[1,-c(1:3)]
  A<-Lake_data[2,-c(1:3)]
  B<-Lake_data[3,-c(1:3)]
  AB<-Lake_data[4,-c(1:3)]
  
  output<-null_compare(Control,A,B,AB)
  output$Difference[is.nan(output$Difference)]<-100
  output$Lake<-lake
  if(lake == "CO"){
    output_full<-output
  } else {output_full<-rbind(output_full,output)}
}

#Thompson and Shurin 2012####
T_S_data<-read.csv("./Experimental data/Thompson and Shurin 2012.csv")


Treatment_means<-T_S_data%>%
  filter(Dispersal==0)%>%
  group_by(Temperature,Salt)%>%
  summarise_each(funs(mean))

Control<-Treatment_means[1,-c(1:3)]
A<-Treatment_means[2,-c(1:3)]
B<-Treatment_means[3,-c(1:3)]
AB<-Treatment_means[4,-c(1:3)]

output<-null_compare(Control,A,B,AB)
output$Difference[is.nan(output$Difference)]<-100
output$Lake<-"T_S_no_disp"
output_full<-rbind(output_full,output)

#with dispersal
Treatment_means<-T_S_data%>%
  filter(Dispersal==1)%>%
  group_by(Temperature,Salt)%>%
  summarise_each(funs(mean))

Control<-Treatment_means[1,-c(1:3)]
A<-Treatment_means[2,-c(1:3)]
B<-Treatment_means[3,-c(1:3)]
AB<-Treatment_means[4,-c(1:3)]

output<-null_compare(Control,A,B,AB)
output$Difference[is.nan(output$Difference)]<-100
output$Lake<-"T_S_disp"
output_full<-rbind(output_full,output)

#Shurin et al 2012####
Shurin_data<-read.csv("./Experimental data/Shurin et al. 2012.csv")

Treatment_means<-Shurin_data%>%
  filter(Nutrients==-1)%>%
  group_by(Warming,Fish)%>%
  summarise_each(funs(mean))

Control<-Treatment_means[1,-c(1:3)]
A<-Treatment_means[2,-c(1:3)]
B<-Treatment_means[3,-c(1:3)]
AB<-Treatment_means[4,-c(1:3)]

output<-null_compare(Control,A,B,AB)
output$Difference[is.nan(output$Difference)]<-100
output$Lake<-"Shurin_no_nut"
output_full<-rbind(output_full,output)

#with nutrients
Treatment_means<-Shurin_data%>%
  filter(Nutrients==1)%>%
  group_by(Warming,Fish)%>%
  summarise_each(funs(mean))

Control<-Treatment_means[1,-c(1:3)]
A<-Treatment_means[2,-c(1:3)]
B<-Treatment_means[3,-c(1:3)]
AB<-Treatment_means[4,-c(1:3)]

output<-null_compare(Control,A,B,AB)
output$Difference[is.nan(output$Difference)]<-100
output$Lake<-"Shurin_nut"
output_full<-rbind(output_full,output)

#clasify net response####
output_full$Syn_vs_Ant <- "NA"
output_full$Syn_vs_Ant[output_full$Difference>0 & output_full$Reversal == F] <- "Synergy"
output_full$Syn_vs_Ant[output_full$Difference<0 & output_full$Reversal == F] <- "Antagonism"
output_full$Syn_vs_Ant[output_full$Reversal == T] <- "Reversal"
output_full$Syn_vs_Ant[output_full$Difference>-0.05 & output_full$Difference<0.05 & output_full$Reversal == F] <- "Match"
output_full$Syn_vs_Ant<-factor(output_full$Syn_vs_Ant,levels = rev(c("Reversal","Antagonism","Match","Synergy")),ordered = T)

pdf("./Figures/Lake null model contrast actual change.pdf",width = 11,height = 8.5)
ggplot(output_full,aes(x=Lake,y=Change, group=interaction(lake,Null_model),color=Null_model,pch=Reversal))+
  geom_point(size=5)+
  facet_grid(Response~.,scale='free')+
  scale_color_brewer(type="qual",palette = 6)+
  theme_bw(base_size = 16)+
  geom_hline(yintercept = 0,linetype=2)
dev.off()

pdf("./Figures/Lake null model contrast.pdf",width = 11,height = 8.5)
ggplot(output_full,aes(x=Lake,y=Difference, group=interaction(lake,Null_model),color=Null_model,pch=Reversal))+
  #geom_point(size=3)+
  geom_jitter(size=5,width=0.25)+
  facet_grid(Response~.)+
  theme_bw(base_size = 16)+
  scale_color_brewer(type="qual",palette = 6)+
  geom_hline(yintercept = 100,linetype=2)+
  ylab("% of predicted")+
  scale_y_log10(breaks=c(100,1000),labels=c(100,1000))
dev.off()

Syn_counts<-output_full%>%
  group_by(Response,Null_model,Syn_vs_Ant) %>%
  summarise(Counts = n())

Syn_counts$Response<-factor(Syn_counts$Response,levels=c("Species richness","Biomass","Composition"),ordered = T)
Syn_counts$Null_model<-factor(Syn_counts$Null_model,levels=c("Additive","Multiplicative","Compositional"),ordered = T)

pdf("./Figures/Synergy vs. antagonism.pdf",width = 11,height = 8.5)
ggplot(filter(Syn_counts,Null_model!="Actual"),aes(x=Null_model,y=Counts,fill=Syn_vs_Ant))+
  facet_wrap(~Response)+
  #scale_fill_brewer(type='div',palette=8,name="",direction = -1)+
  scale_fill_manual(values = c("#d7191c","#1a9641","#2c7bb6","grey","lightgreen","black"),name="")+
  geom_bar(stat="identity",color="black")+
  theme_bw(base_size = 12)+
  xlab("Null model")+
  scale_y_continuous(breaks=seq(1:9))+
  xlab("Null model")
dev.off()

ggplot(filter(output_full,Null_model!="Actual"),aes(x=Null_model,y=Lake,fill=Syn_vs_Ant))+
  facet_wrap(~Response)+
  #scale_fill_brewer(type='div',palette=8,name="",direction = -1)+
  scale_fill_manual(values = c("#d7191c","#1a9641","#2c7bb6","grey","lightgreen","black"),name="")+
  geom_tile()+
  theme_bw(base_size = 16)+
  xlab("Null model")

is.same<-function(x){x[1]==x}

output_full%>%
  filter(Null_model=="Additive" | Null_model=="Compositional", Response=="Species richness")%>%
  group_by(Lake)%>%
  summarise(Same=is.same(Syn_vs_Ant)[2])


