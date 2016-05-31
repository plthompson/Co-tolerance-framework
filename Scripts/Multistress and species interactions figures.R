#Figure 1####
Stress_type<-c("-,-","-,+","+,+","mixed")
species<-20

Co_sensitivityV<-c("Random","Positive","Negative")


for(st in 1:length(Stress_type)){
  for(ct in 1:3){
    #environmental response####
    #start with sensitivity to stress
    sensitivity_scaler<-0.005
    min_sensitivity<-0.3
    if(ct==1){
      #sensitivity<-mvrnorm(n=species,c(0,0),Sigma = matrix(0,2,2)+diag(2)*0.05)
      sensitivity<--matrix(runif(species*2,min = min_sensitivity*sensitivity_scaler,1*sensitivity_scaler),species,2)
    } else {
      #sensitivity<-mvrnorm(n=species,c(0,0),Sigma = matrix(0.9,2,2)+diag(2)*0.05)
      sensitivity1<-mvrnorm(n=species,mu = c(0,0),Sigma = matrix(c(1, 0.9, 0.9, 1), nrow = 2))
      sensitivity<--(pnorm(sensitivity1)*(1-min_sensitivity)+min_sensitivity)*sensitivity_scaler
    }
    if(ct==3){
      sensitivity1[,2]<-sensitivity1[,2]*-1
      sensitivity <- -(pnorm(sensitivity1)*(1-min_sensitivity)+min_sensitivity)*sensitivity_scaler
    }
    if(st==2){
      sensitivity[,2]<-sensitivity[,2]*-1
    }
    if(st==3){
      sensitivity<-sensitivity*-1
    }
    if(st==4){
      if(ct==1){
        sensitivity<-matrix(runif(species*2,min = -sensitivity_scaler,max=sensitivity_scaler),species,2)
      } else {
        sensitivity1<-mvrnorm(n=species,mu = c(0,0),Sigma = matrix(c(1, 0.9, 0.9, 1), nrow = 2))
        sensitivity<-((pnorm(sensitivity1)*2)-1)*sensitivity_scaler
      }
      if(ct==3){
        sensitivity1[,2]<-sensitivity1[,2]*-1
        sensitivity<-((pnorm(sensitivity1)*2)-1)*sensitivity_scaler
      }    
    }
    sens.hold<-data.frame(A=sensitivity[,1],B=sensitivity[,2]) 
    sens.hold$ct=Co_sensitivityV[ct]
    sens.hold$st=Stress_type[st]
    if(ct==1 & st==1){
      sens.df<-sens.hold
    } else {
      sens.df<-rbind(sens.df,sens.hold)
    }
  }
}

sens.df$ct<-factor(sens.df$ct,levels=Co_sensitivityV,ordered = T)


line.df<-data.frame(x=seq(from = -0.005,to = 0.005,length=1000),y=-0.0055)

tri_fill <- data.frame(x= c(-0.005,0.005,0.005,-0.005), 
                       y = c(-0.005,-0.005,0.005,0.005)) 

n_polys <- 1000
#identify changepoint (very data/shape dependent)
change_point <-  max(tri_fill$x[which(tri_fill$y==max(tri_fill$y))])

#calculate slope and intercept
slope <- (max(tri_fill$y)-min(tri_fill$y))/ (change_point - max(tri_fill$x))
intercept <- max(tri_fill$y)
#calculate sequence of borders: x, and accompanying lower and upper y coordinates
x_seq <- seq(min(tri_fill$x),max(tri_fill$x),length.out=n_polys+1)
y_max_seq <- ifelse(x_seq<=change_point, max(tri_fill$y), intercept + (x_seq - change_point)*slope)
y_min_seq <- rep(min(tri_fill$y), n_polys+1)

#create polygons/rectangles
poly_list <- lapply(1:n_polys, function(p){
  res <- data.frame(x=rep(c(x_seq[p],x_seq[p+1]),each=2),
                    y = c(y_min_seq[p], y_max_seq[p:(p+1)], y_min_seq[p+1]))
  res$fill_id <- x_seq[p]
  res
}
)

poly_data <- do.call(rbind, poly_list)

ggplot(sens.df,aes(x=A,y=B))+
  facet_grid(ct~st)+
  geom_polygon(data=poly_data, aes(x=x,y=y, group=fill_id,fill=fill_id),alpha=1) +
  geom_polygon(data=poly_data, aes(y=x,x=y, group=fill_id,fill=fill_id),alpha=0.5)+
  scale_fill_gradient2(high = "dodgerblue1",mid="white",low="goldenrod1",name="Stress\neffect")+
  #scale_fill_viridis(option = "B")+
  theme_bw(base_size = 10)+
  theme(panel.margin = unit(0.5, "cm"))+
  removeGrid()+
  geom_hline(yintercept = 0, linetype=2)+
  geom_vline(xintercept = 0, linetype=2)+
  geom_point(data = sens.df,aes(x=A,y=B))+
  ylab("Response to stressor B")+
  xlab("Response to stressor A")
ggsave("./Figures/Species interactions - Fig 1.png",height=8.5,width=11)

#Figure 2####
library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(dplyr)
load("Workspace/Multistress.RData")
Output_means$Response<-factor(Output_means$Response,levels=c("Species richness","Biomass", "Composition"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,levels=c("A_only","B_only","Actual","Additive","Multiplicative","Species_specific"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,labels=c("A_only","B_only","Actual","Additive","Multiplicative","Compositional"))

lineV2<-c(1,2,1,1)

#Figure 2####
ColV<-c("grey30",brewer.pal(3,"Set1"))
ggplot(filter(Output_means,CoTolerance=="Random",Null_model=="Compositional" | Null_model=="Actual",Response != "Composition"),aes(x=Stress,y=Change_mean,color=Interactions, fill=Interactions,group=interaction(Interactions,Null_model), linetype=Null_model))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype(name="",labels=c("Realized","Predicted"))+
  theme_bw()+
  removeGrid()+
  ylab("Change from control")
ggsave("./Figures/Species interactions - Fig 2.pdf",width = 11, height=8.5)


#Figure 3####
ggplot(filter(Output_means,CoTolerance=="Random",Null_model=="Compositional",Response != "Composition"),aes(x=Stress,y=Difference_mean,color=Interactions, fill=Interactions,group=Interactions))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Difference_lower,ymax=Difference_upper),alpha=0.2,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  theme_bw()+
  removeGrid()+
  ylab("Difference from null model")
ggsave("./Figures/Species interactions - Fig 3.pdf",width = 11, height=8.5)

ggplot(filter(Output_means,CoTolerance=="Negative",Null_model=="Compositional",Response != "Composition"),aes(x=Stress,y=Difference_mean,color=Interactions, fill=Interactions,group=Interactions))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Difference_lower,ymax=Difference_upper),alpha=0.2,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  theme_bw()+
  removeGrid()+
  ylab("Difference from null model")
ggsave("./Figures/Species interactions - Fig 3a.pdf",width = 11, height=8.5)

ggplot(filter(Output_means,CoTolerance=="Positive",Null_model=="Compositional",Response != "Composition"),aes(x=Stress,y=Difference_mean,color=Interactions, fill=Interactions,group=Interactions))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Difference_lower,ymax=Difference_upper),alpha=0.2,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  theme_bw()+
  removeGrid()+
  ylab("Difference from null model")
ggsave("./Figures/Species interactions - Fig 3b.pdf",width = 11, height=8.5)


#Figure 4####
Output_means_trophic$Response<-factor(Output_means_trophic$Response,levels=c("Species richness","Biomass", "Composition"), ordered=T)

ggplot(filter(Output_means_trophic,CoTolerance=="Random",Null_model=="Species_specific" | Null_model=="Actual",Response != "Composition"),aes(x=Stress,y=Change_mean,color=Trophic_level, fill=Trophic_level,group=interaction(Trophic_level,Null_model), linetype=Null_model))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_brewer(palette = "Set2",name="")+
  scale_fill_brewer(palette = "Set2",name="")+
  scale_linetype(name="",labels=c("Realized","Predicted"))+
  theme_bw()+
  removeGrid()+
  ylab("Change from control")
ggsave("./Figures/Species interactions - Fig 4a.pdf",width = 11, height=8.5)


ggplot(filter(Output_means_trophic,CoTolerance=="Random",Null_model=="Species_specific",Response != "Composition"),aes(x=Stress,y=Difference_mean,color=Trophic_level, fill=Trophic_level,group=Trophic_level))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Difference_lower,ymax=Difference_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_brewer(palette = "Set2",name="")+
  scale_fill_brewer(palette = "Set2",name="")+
  theme_bw()+
  removeGrid()+
  ylab("Difference from null model")
ggsave("./Figures/Species interactions - Fig 4.pdf",width = 11, height=8.5)


#Figure 5####
ggplot(filter(Output_means,CoTolerance=="Random",Null_model=="Compositional", Response=="Biomass"),aes(x=Stress,y=Reversal_mean*1,color=Interactions, fill=Interactions,group=Interactions))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Reversal_lower,ymax=Reversal_upper),alpha=0.2,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  theme_bw()+
  removeGrid()+
  ylab("Frequency of reversals")
ggsave("./Figures/Species interactions - Fig. 5.pdf",width = 11, height=3.5)

