library(vegan)
library(MASS)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)

null_compare<-function(tseries){
  A<-tseries[,,"A"]
  B<-tseries[,,"B"]
  AB<-tseries[,,"AB"]
  Control<-tseries[1,,"A"]
  output<-data.frame(Stress=stressV[samp_stress],Response=rep(c("Biomass","Species richness","Composition"),each=length(samp_stress)),Actual=NA,Species_specific=NA,Additive=NA,Multiplicative=NA,A_only=NA,B_only=NA)
  
  output<-mutate(output,Actual=replace(Actual,Response=="Biomass",(rowSums(AB)-sum(Control))/sum(Control)))
  output<-mutate(output,Actual=replace(Actual,Response=="Species richness",(rowSums(AB>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Actual=replace(Actual,Response=="Composition",  c(0,vegdist(AB,na.rm=T)[1:(nrow(AB)-1)])))
  
  output<-mutate(output,A_only=replace(A_only,Response=="Biomass",(rowSums(A)-sum(Control))/sum(Control)))
  output<-mutate(output,A_only=replace(A_only,Response=="Species richness",(rowSums(A>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,A_only=replace(A_only,Response=="Composition",  c(0,vegdist(A,na.rm=T)[1:(nrow(A)-1)])))
  
  output<-mutate(output,B_only=replace(B_only,Response=="Biomass",(rowSums(B)-sum(Control))/sum(Control)))
  output<-mutate(output,B_only=replace(B_only,Response=="Species richness",(rowSums(B>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,B_only=replace(B_only,Response=="Composition",  c(0,vegdist(B,na.rm=T)[1:(nrow(B)-1)])))
  
  AB_predict<-rep(Control,each=nrow(A))+(A-rep(Control,each=nrow(A)))+(B-rep(Control,each=nrow(A)))
  AB_predict[AB_predict<0]<-0
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Biomass",(rowSums(AB_predict)-sum(Control))/sum(Control)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Species richness",(rowSums(AB_predict>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Composition",  c(0,vegdist(AB_predict)[1:(nrow(AB_predict)-1)])))
  
  pA<-(rowSums(A)-sum(Control))/sum(Control)
  pB<-(rowSums(B)-sum(Control))/sum(Control)
  addM<-pA+pB
  addM[addM<(-1)]<--1
  output<-mutate(output,Additive=replace(Additive,Response=="Biomass",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Biomass",pA+pB+(pA*pB)))
  
  pA<-(rowSums(A>0)-sum(Control>0))/sum(Control>0)
  pB<-(rowSums(B>0)-sum(Control>0))/sum(Control>0)
  addM<-pA+pB
  addM[addM<(-1)]<--1
  output<-mutate(output,Additive=replace(Additive,Response=="Species richness",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Species richness",pA+pB+(pA*pB)))
  
  pA<-c(0,vegdist(A)[1:(nrow(A)-1)])
  pB<-c(0,vegdist(B)[1:(nrow(B)-1)])
  addM<-pA+pB
  addM[addM>1]<-1
  output<-mutate(output,Additive=replace(Additive,Response=="Composition",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Composition",pA+pB-(pA*pB)))
  
  
  output<-gather(output,key = Null_model,value=Change,Actual:B_only)
  output<-output %>%
    group_by(Response,Stress) %>%
    mutate(Difference=abs(Change[Null_model=="Actual"])-abs(Change))%>%
    mutate(Reversal = Change>0 & Change[Null_model=="Actual"]<0 | Change<0 & Change[Null_model=="Actual"]>0)
  return(output)
}

nStressors<-2
reps<-20
Co_sensitivityV<-c("Random","Positive","Negative")
Com_types<-c("No interactions","Competitive","Facilitative","Trophic","Mixed")
Stress_type<-c("-,-","-,+","+,+","mixed")
int_strength<-seq(0,1,by=0.25)

#environmental change####
stress_increase<-6000
stress<-65
stressV<-c(seq(0,stress,length=stress*stress_increase+1))
Tmax<-length(stressV)

species<-20
nprey<-species*0.5
nherb<-species*0.3
npred<-species*0.2
trophV<-factor(c(rep("Plants",nprey),rep("Herbivores",nherb),rep("Predators",npred)),levels = c("Plants","Herbivores","Predators"),ordered = T)

samp_stress<-seq(1,Tmax,by=stress_increase)
recolV<-seq(100,Tmax,by=100)

#weight by species number
weight=1/(species*4)

#growth rate
C<-rep(0.1,species)
C_troph<-c(rep(0.1,nprey),rep(0,species-nprey))

for(r in 1:reps){
  #make communities####
  BB<-array(0,dim=c(species, species,length(Com_types)),dimnames=list(1:species,1:species,Com_types))
  
  intra<--.2
  diag(BB[,,"No interactions"])<-intra*weight
  
  
  #competition
  repeat{
    b11<--.15
    BI<-b11*matrix(runif(species*species),species,species)
    diag(BI)<-intra
    BI<-BI*weight
    
    if(sum(solve(-BI,C)>0)==species) break
  }
  BB[,,"Competitive"]<-BI
  
  #facilitation
  repeat{
    b11<--.15
    BI<-b11*matrix(runif(species*species),species,species)*-0.1
    diag(BI)<-intra
    BI<-BI*weight
    if(sum(solve(-BI,C)>0)==species) break
  }
  BB[,,"Facilitative"]<-BI
  
  #trophic
  repeat{
    b11<-0#-.15
    b12<--0.15#-0.3
    b21<-0.1
    b23<--.1
    b32<- 0.07#.08
    bdiag1<--.2
    bdiag2<--.2
    
    #tritrophic BB Matrix####
    B11<-b11*matrix(runif(nprey*nprey),nprey,nprey)
    B12<-b12*matrix(runif(nprey*nherb),nprey,nherb)
    B13<-matrix(0,nprey,npred)
    B21<-b21*matrix(runif(nherb*nprey),nherb,nprey)
    B22<-matrix(0,nherb,nherb)
    B23<-b23*matrix(runif(nherb*npred),nherb,npred)
    B31<-matrix(0,npred,nprey)
    B32<-b32*matrix(runif(npred*nherb),npred,nherb)
    B33<-matrix(0,npred,npred)
    BI<-rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
    diag(BI)<-bdiag1
    diag(BI[(nprey+nherb+1):species,(nprey+nherb+1):species])<-bdiag2
    if(sum(solve(-BI,C_troph)>0)==species) break
  }
  BB[,,"Trophic"]<-BI
  
  #mixed
  repeat{
    b11<--0.10
    b11_pos<-0.15
    b12<--0.15#-0.3
    b21<-0.1
    b23<--.1
    b32<- 0.07#.08
    bdiag1<--.2
    bdiag2<--.2
    
    #tritrophic BB Matrix####
    B11<-b11*matrix(runif(nprey*nprey),nprey,nprey)
    mutuals<-rbinom(n = ((nprey*nprey-nprey)/2),size = 1,prob = 0.2)
    B11[upper.tri(B11)][mutuals]<-B11[upper.tri(B11)][mutuals]*-b11_pos
    B11[lower.tri(B11)][mutuals]<-B11[lower.tri(B11)][mutuals]*-b11_pos
    B12<-b12*matrix(runif(nprey*nherb),nprey,nherb)
    B13<-matrix(0,nprey,npred)
    B21<-b21*matrix(runif(nherb*nprey),nherb,nprey)
    B22<-matrix(0,nherb,nherb)
    B23<-b23*matrix(runif(nherb*npred),nherb,npred)
    B31<-matrix(0,npred,nprey)
    B32<-b32*matrix(runif(npred*nherb),npred,nherb)
    B33<-matrix(0,npred,npred)
    BI<-rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
    diag(BI)<-bdiag1
    diag(BI[(nprey+nherb+1):species,(nprey+nherb+1):species])<-bdiag2
    if(sum(solve(-BI,C_troph)>0)==species) break
  }
  BB[,,"Mixed"]<-BI
  
  for(st in 1:length(Stress_type)){
    print(paste("Rep",r, "-",st))
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
      
      #starting abundances
      X<-array(NA,dim=c(species,nStressors+1,length(Com_types)),dimnames = list(1:species,c("A","B","AB"),Com_types))
      X[,,"No interactions"]<-solve(-BB[,,"No interactions"],C)
      X[,,"Competitive"]<-solve(-BB[,,"Competitive"],C)
      X[,,"Facilitative"]<-solve(-BB[,,"Facilitative"],C)
      X[,,"Trophic"]<-solve(-BB[,,"Trophic"],C_troph)
      X[,,"Mixed"]<-solve(-BB[,,"Mixed"],C_troph)
      
      Xsave<-array(NA,dim=c(length(samp_stress),species,nStressors+1,length(Com_types)),dimnames = list(samp_stress,1:species,c("A","B","AB"),Com_types))
      Xsave[1,,,]<-X
      Env_resp<-matrix(NA,species,3,dimnames = list(1:species,c("A","B","AB")))
      
      for(l in 1:Tmax){
        Env_resp[,"A"]<-sensitivity[,1]*stressV[l]
        Env_resp[,"B"]<-sensitivity[,2]*stressV[l]
        Env_resp[,"AB"]<-sensitivity[,1]*stressV[l]+sensitivity[,2]*stressV[l]
        
        for(com in 1:length(Com_types)){
          if(com<4){
            Xt<-X[,,com]*exp(C+BB[,,com]%*%X[,,com]+Env_resp)
          } else {
            Xt<-X[,,com]*exp(C_troph+BB[,,com]%*%X[,,com]+Env_resp)
          }
          Xt[Xt<0.01]<-0
          X[,,com]<-Xt
        }
        if(l>1 & sum(l==recolV)==1){
          recol<-rbinom(n = species,size = 1,prob = 0.01)
          X[recol==1,,]<-X[recol==1,,]+0.02
        }
        if(l>1 & sum(l==samp_stress)==1){
          Xsave[which(l==samp_stress),,,]<-X
        }
      }
      
      for(com in 1:length(Com_types)){
        hold<-null_compare(Xsave[,,,com])
        hold$Interactions<-Com_types[com]
        hold$CoTolerance<-c("Random","Positive","Negative")[ct]
        hold$Stress_type<-Stress_type[st]
        hold$Rep<-r
        if(com==4){ 
          for(tlevel in 1:3){
            hold2<-null_compare(Xsave[,trophV==unique(trophV)[tlevel],,com])
            hold2$Trophic_level<-unique(trophV)[tlevel]
            hold2$CoTolerance<-c("Random","Positive","Negative")[ct]
            hold2$Stress_type<-Stress_type[st]
            hold2$Rep<-r
            if(r==1 & ct==1 & st==1 & tlevel==1){
              Output_trophic<-hold2
            } else{
              Output_trophic<-rbind(Output_trophic,hold2)}
          }
        }
        if(r==1 & ct==1 & st==1 & com==1 ){
          Output<-hold
        } else{
          Output<-rbind(Output,hold)}
      }
    }
  }
}

Output$Interactions<-factor(Output$Interactions,levels=Com_types,ordered = T)
Output$CoTolerance<-factor(Output$CoTolerance,levels=c("Positive","Random","Negative"),ordered = T)
Output$Reversal<-as.numeric(Output$Reversal)

Output_means<-Output%>%
  group_by(Stress,Response,Null_model,Interactions,CoTolerance,Stress_type)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)),-Rep)

Output_trophic$Trophic_level<-factor(Output_trophic$Trophic_level,levels=unique(trophV),ordered = T)
Output_trophic$CoTolerance<-factor(Output_trophic$CoTolerance,levels=c("Positive","Random","Negative"),ordered = T)
Output_trophic$Reversal<-as.numeric(Output_trophic$Reversal)

Output_means_trophic<-Output_trophic%>%
  group_by(Stress,Response,Null_model,CoTolerance,Stress_type,Trophic_level)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)),-Rep)

save(Output_means,Output_means_trophic,file="./Workspace/Multistress - recol.RData")

library(ggplot2)
library(RColorBrewer)
library(ggExtra)
library(dplyr)
load("Workspace/Multistress - recol.RData")
Output_means<-Output_means%>%
  ungroup()%>%
  mutate(CoTolerance=replace(CoTolerance,CoTolerance=="Negative" & Stress_type =="-,+","hold"),
         CoTolerance=replace(CoTolerance,CoTolerance=="Positive" & Stress_type =="-,+","Negative"),
         CoTolerance=replace(CoTolerance,is.na(CoTolerance) & Stress_type =="-,+","Positive"))

Output_means$Response<-factor(Output_means$Response,levels=c("Species richness","Biomass", "Composition"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,levels=c("A_only","B_only","Actual","Additive","Multiplicative","Species_specific"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,labels=c("A_only","B_only","Actual","Additive","Multiplicative","Compositional"))

lineV2<-c(1,2,1,1)

#Figure 1####
ColV<-c("grey30",brewer.pal(4,"Set1"))
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
  ylab("Difference relative to control community")+
  xlab(expression(paste("Stressor magnitude (",italic("m"),")",sep="")))
#ggsave("./Figures/Species interactions - Fig S4.pdf",width = 8, height=5)

#Figure 2####
ggplot(filter(Output_means,CoTolerance=="Random",Null_model=="Compositional",Response != "Composition"),aes(x=Stress,y=Difference_mean,color=Interactions, fill=Interactions,group=Interactions))+
  geom_hline(yintercept = 0,linetype=2,col=1)+
  geom_ribbon(aes(ymin=Difference_lower,ymax=Difference_upper),alpha=0.2,color=NA)+
  geom_line(size=1)+
  facet_grid(Response~Stress_type,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  theme_bw()+
  removeGrid()+
  ylab(expression(paste("Difference from null model (",italic(D["AB"]),")",sep="")))+
  xlab(expression(paste("Stressor magnitude (",italic("m"),")",sep="")))
ggsave("./Figures/Species interactions - Fig S4.png",width = 8, height=5)

