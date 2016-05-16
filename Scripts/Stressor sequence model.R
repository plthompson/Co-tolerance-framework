library(vegan)
library(MASS)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)

null_compare2<-function(tseries){
  A<-tseries[3,,"A"]
  B<-tseries[3,,"B"]
  AB<-tseries[3,,"AB"]
  A_AB<-tseries[3,,"A->AB"]
  B_AB<-tseries[3,,"B->AB"]
  Control<-tseries[1,,"A"]
  output<-data.frame(Response=rep(c("Biomass","Species richness","Composition")),Actual=NA,Species_specific=NA,Additive=NA,Multiplicative=NA,A_only=NA,B_only=NA,A_first=NA,B_first=NA)
  
  output<-mutate(output,Actual=replace(Actual,Response=="Biomass",(sum(AB)-sum(Control))/sum(Control)))
  output<-mutate(output,Actual=replace(Actual,Response=="Species richness",(sum(AB>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Actual=replace(Actual,Response=="Composition",vegdist(rbind(Control,AB),na.rm=T)))
  
  output<-mutate(output,A_only=replace(A_only,Response=="Biomass",(sum(A)-sum(Control))/sum(Control)))
  output<-mutate(output,A_only=replace(A_only,Response=="Species richness",(sum(A>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,A_only=replace(A_only,Response=="Composition",vegdist(rbind(Control,A),na.rm=T)))
  
  output<-mutate(output,B_only=replace(B_only,Response=="Biomass",(sum(B)-sum(Control))/sum(Control)))
  output<-mutate(output,B_only=replace(B_only,Response=="Species richness",(sum(B>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,B_only=replace(B_only,Response=="Composition",vegdist(rbind(Control,B),na.rm=T)))
  
  output<-mutate(output,A_first=replace(A_first,Response=="Biomass",(sum(A_AB)-sum(Control))/sum(Control)))
  output<-mutate(output,A_first=replace(A_first,Response=="Species richness",(sum(A_AB>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,A_first=replace(A_first,Response=="Composition",vegdist(rbind(Control,A_AB),na.rm=T)))
  
  output<-mutate(output,B_first=replace(B_first,Response=="Biomass",(sum(B_AB)-sum(Control))/sum(Control)))
  output<-mutate(output,B_first=replace(B_first,Response=="Species richness",(sum(B_AB>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,B_first=replace(B_first,Response=="Composition",vegdist(rbind(Control,B_AB),na.rm=T)))
  
  AB_predict<-Control+(A-Control)+(B-Control)
  AB_predict[AB_predict<0]<-0
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Biomass",(sum(AB_predict)-sum(Control))/sum(Control)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Species richness",(sum(AB_predict>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Composition", vegdist(rbind(Control,AB_predict),na.rm=T)))
  
  pA<-(sum(A)-sum(Control))/sum(Control)
  pB<-(sum(B)-sum(Control))/sum(Control)
  addM<-pA+pB
  addM[addM<(-1)]<--1
  output<-mutate(output,Additive=replace(Additive,Response=="Biomass",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Biomass",pA+pB+(pA*pB)))
  
  pA<-(sum(A>0)-sum(Control>0))/sum(Control>0)
  pB<-(sum(B>0)-sum(Control>0))/sum(Control>0)
  addM<-pA+pB
  addM[addM<(-1)]<--1
  output<-mutate(output,Additive=replace(Additive,Response=="Species richness",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Species richness",pA+pB+(pA*pB)))
  

  pA<-vegdist(rbind(Control,A),na.rm=T)
  pB<-vegdist(rbind(Control,B),na.rm=T)
  addM<-pA+pB
  addM[addM>1]<-1
  output<-mutate(output,Additive=replace(Additive,Response=="Composition",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Composition",pA+pB-(pA*pB)))
  
  
  output<-gather(output,key = Null_model,value=Change,Actual:B_first)
  output<-output %>%
    group_by(Response) %>%
    mutate(Difference=abs(Change[Null_model=="Actual"])-abs(Change))%>%
    mutate(Reversal = Change>0 & Change[Null_model=="Actual"]<0 | Change<0 & Change[Null_model=="Actual"]>0)
  return(output)
}

nStressors<-2
reps<-100
Co_sensitivityV<-c("Random","Positive","Negative")
Com_types<-c("No interactions","Competitive","Facilitative","Trophic")
Stress_type<-c("-,-","-,+","+,+","mixed")

#environmental change####
stress_period<-5000
stress<-15
stressV<-rep(c(1,2),each=stress_period)
Tmax<-length(stressV)

species<-20
nprey<-species*0.5
nherb<-species*0.3
npred<-species*0.2
trophV<-factor(c(rep("Plants",nprey),rep("Herbivores",nherb),rep("Predators",npred)),levels = c("Plants","Herbivores","Predators"),ordered = T)

samp_stress<-c(0,stress_period,Tmax)

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
  
  for(st in c(1:4)){
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
      X<-array(NA,dim=c(species,nStressors+3,length(Com_types)),dimnames = list(1:species,c("A","B","AB","A->AB","B->AB"),Com_types))
      X[,,"No interactions"]<-solve(-BB[,,"No interactions"],C)
      X[,,"Competitive"]<-solve(-BB[,,"Competitive"],C)
      X[,,"Facilitative"]<-solve(-BB[,,"Facilitative"],C)
      X[,,"Trophic"]<-solve(-BB[,,"Trophic"],C_troph)
      
      Xsave<-array(NA,dim=c(length(samp_stress),species,nStressors+3,length(Com_types)),dimnames = list(samp_stress,1:species,c("A","B","AB","A->AB","B->AB"),Com_types))
      Xsave[1,,,]<-X
      Env_resp<-matrix(NA,species,5,dimnames = list(1:species,c("A","B","AB","A->AB","B->AB")))
      
      for(l in 1:Tmax){
        Env_resp[,"A"]<-sensitivity[,1]*stress
        Env_resp[,"B"]<-sensitivity[,2]*stress
        Env_resp[,"AB"]<-sensitivity[,1]*stress+sensitivity[,2]*stress
        Env_resp[,"A->AB"]<-sensitivity[,1]*stress+sensitivity[,2]*stress*(stressV[l]-1)
        Env_resp[,"B->AB"]<-sensitivity[,2]*stress+sensitivity[,1]*stress*(stressV[l]-1)
        
        for(com in 1:length(Com_types)){
          if(com<4){
            Xt<-X[,,com]*exp(C+BB[,,com]%*%X[,,com]+Env_resp)
          } else {
            Xt<-X[,,com]*exp(C_troph+BB[,,com]%*%X[,,com]+Env_resp)
          }
          Xt[Xt<0.01]<-0
          X[,,com]<-Xt
        }
        if(l>1 & sum(l==samp_stress)==1){
          Xsave[which(l==samp_stress),,,]<-X
        }
      }
      
      for(com in 1:length(Com_types)){
        hold<-null_compare2(Xsave[,,,com])
        hold$Interactions<-Com_types[com]
        hold$CoTolerance<-c("Random","Positive","Negative")[ct]
        hold$Stress_type<-Stress_type[st]
        hold$Rep<-r
        if(com==4){ 
          for(tlevel in 1:3){
            hold2<-null_compare2(Xsave[,trophV==unique(trophV)[tlevel],,com])
            hold2$Trophic_level<-unique(trophV)[tlevel]
            hold2$CoTolerance<-c("Random","Positive","Negative")[ct]
            hold2$Stress_type<-Stress_type[st]
            hold2$Rep<-r
            if(r==1 & ct==1 & tlevel==1){
              Output_trophic<-hold2
            } else{
              Output_trophic<-rbind(Output_trophic,hold2)}
          }
        }
        if(r==1 & ct==1 & com==1 ){
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
Output$Response<-factor(Output$Response,levels = c("Species richness", "Biomass","Composition"), ordered=T)

save(Output,Output_trophic,file = "Multistress_sequence.RData")

AB_diff<-Output%>%
  filter(Null_model=="A_first" |Null_model == "B_first")%>%
  group_by(Response,Interactions, CoTolerance,Rep,Stress_type)%>%
  summarise(order_diff=abs(Change[Null_model=="A_first"]-Change[Null_model=="B_first"]))

ggplot(filter(AB_diff, Response=="Species richness"),aes(x=CoTolerance,y=order_diff, fill=CoTolerance))+
  geom_boxplot()+
  #geom_violin()+
  theme_bw()+
  removeGrid()+
  scale_fill_manual(values = c("#d7191c","grey50","#2c7bb6"),guide=F)+
  facet_grid(Stress_type~Interactions, scales = "free")+
  ylab("Sequence effect (species richness)")+
  xlab("Co-tolerance")
ggsave("Sequence - species richness.pdf",width = 11, height=8.5)

ggplot(filter(AB_diff, Response=="Biomass"),aes(x=CoTolerance,y=order_diff, fill=CoTolerance))+
  geom_boxplot()+
  #geom_violin()+
  theme_bw()+
  removeGrid()+
  scale_fill_manual(values = c("#d7191c","grey50","#2c7bb6"),guide=F)+
  facet_grid(Stress_type~Interactions, scales = "free")+
  ylab("Sequence effect (community biomass)")+
  xlab("Co-tolerance")
ggsave("Sequence - biomass.pdf",width = 11, height=8.5)

ggplot(filter(AB_diff, Response=="Composition"),aes(x=CoTolerance,y=order_diff, fill=CoTolerance))+
  geom_boxplot()+
  #geom_violin()+
  theme_bw()+
  removeGrid()+
  scale_fill_manual(values = c("#d7191c","grey50","#2c7bb6"),guide=F)+
  facet_grid(Stress_type~Interactions, scales = "free")+
  ylab("Sequence effect (community composition)")+
  xlab("Co-tolerance")
ggsave("Sequence - composition.pdf",width = 11, height=8.5)
