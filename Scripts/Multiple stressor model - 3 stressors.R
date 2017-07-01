library(vegan)
library(MASS)
library(ggplot2)
library(ggExtra)
library(RColorBrewer)
library(dplyr)
library(tidyr)

null_compare3<-function(tseries){
  A<-tseries[,,"A"]
  B<-tseries[,,"B"]
  C<-tseries[,,"C"]
  ABC<-tseries[,,"ABC"]
  Control<-tseries[1,,"A"]
  output<-data.frame(Stress=stressV[samp_stress],Response=rep(c("Biomass","Species richness","Composition"),each=length(samp_stress)),Actual=NA,Species_specific=NA, SS_mod<-NA ,Additive=NA,Multiplicative=NA,A_only=NA,B_only=NA, C_only=NA)
  
  output<-mutate(output,Actual=replace(Actual,Response=="Biomass",(rowSums(ABC)-sum(Control))/sum(Control)))
  output<-mutate(output,Actual=replace(Actual,Response=="Species richness",(rowSums(ABC>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Actual=replace(Actual,Response=="Composition",  c(0,vegdist(ABC,na.rm=T)[1:(nrow(ABC)-1)])))
  
  output<-mutate(output,A_only=replace(A_only,Response=="Biomass",(rowSums(A)-sum(Control))/sum(Control)))
  output<-mutate(output,A_only=replace(A_only,Response=="Species richness",(rowSums(A>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,A_only=replace(A_only,Response=="Composition",  c(0,vegdist(A,na.rm=T)[1:(nrow(A)-1)])))
  
  output<-mutate(output,B_only=replace(B_only,Response=="Biomass",(rowSums(B)-sum(Control))/sum(Control)))
  output<-mutate(output,B_only=replace(B_only,Response=="Species richness",(rowSums(B>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,B_only=replace(B_only,Response=="Composition",  c(0,vegdist(B,na.rm=T)[1:(nrow(B)-1)])))
  
  output<-mutate(output,C_only=replace(C_only,Response=="Biomass",(rowSums(C)-sum(Control))/sum(Control)))
  output<-mutate(output,C_only=replace(C_only,Response=="Species richness",(rowSums(C>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,C_only=replace(C_only,Response=="Composition",  c(0,vegdist(C,na.rm=T)[1:(nrow(C)-1)])))
  
  ABC_predict<-rep(Control,each=nrow(A))+(A-rep(Control,each=nrow(A)))+(B-rep(Control,each=nrow(A)))+(C-rep(Control,each=nrow(A)))
  ABC_predict[ABC_predict<0]<-0
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Biomass",(rowSums(ABC_predict)-sum(Control))/sum(Control)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Species richness",(rowSums(ABC_predict>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,Species_specific=replace(Species_specific,Response=="Composition",  c(0,vegdist(ABC_predict)[1:(nrow(ABC_predict)-1)])))
  
  ABC_predict<-rep(Control,each=nrow(A))+(A-rep(Control,each=nrow(A)))+(B-rep(Control,each=nrow(A)))+(C-rep(Control,each=nrow(A)))
  ABC_present<-(rep(Control,each=nrow(A))+(A-rep(Control,each=nrow(A)))>0)*(rep(Control,each=nrow(A))+(B-rep(Control,each=nrow(A)))>0)*(rep(Control,each=nrow(A))+(C-rep(Control,each=nrow(A)))>0)
  ABC_predict<-ABC_predict*ABC_present
  ABC_predict[ABC_predict<0]<-0
  output<-mutate(output,SS_mod=replace(SS_mod,Response=="Biomass",(rowSums(ABC_predict)-sum(Control))/sum(Control)))
  output<-mutate(output,SS_mod=replace(SS_mod,Response=="Species richness",(rowSums(ABC_predict>0)-sum(Control>0))/sum(Control>0)))
  output<-mutate(output,SS_mod=replace(SS_mod,Response=="Composition",  c(0,vegdist(ABC_predict)[1:(nrow(ABC_predict)-1)])))
  
  pA<-(rowSums(A)-sum(Control))/sum(Control)
  pB<-(rowSums(B)-sum(Control))/sum(Control)
  pC<-(rowSums(C)-sum(Control))/sum(Control)
  addM<-pA+pB+pC
  addM[addM<(-1)]<--1
  multM<-((1+pA)*(1+pB)*(1+pC))-1

  output<-mutate(output,Additive=replace(Additive,Response=="Biomass",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Biomass",multM))
  
  pA<-(rowSums(A>0)-sum(Control>0))/sum(Control>0)
  pB<-(rowSums(A>0)-sum(Control>0))/sum(Control>0)
  pC<-(rowSums(A>0)-sum(Control>0))/sum(Control>0)
  addM<-pA+pB+pC
  addM[addM<(-1)]<--1
  multM<-((1+pA)*(1+pB)*(1+pC))-1
  output<-mutate(output,Additive=replace(Additive,Response=="Species richness",addM))
  output<-mutate(output,Multiplicative=replace(Multiplicative,Response=="Species richness",multM))
  
  
  output<-gather(output,key = Null_model,value=Change,Actual:SS_mod)
  output<-output %>%
    group_by(Response,Stress) %>%
    mutate(Difference=abs(Change[Null_model=="Actual"])-abs(Change))%>%
    mutate(Reversal = Change>0 & Change[Null_model=="Actual"]<0 | Change<0 & Change[Null_model=="Actual"]>0)
  return(output)
}

nStressors<-3
reps<-100
Co_sensitivityV<-c("Random","Positive","Negative")
Com_types<-c("No interactions","Competitive","Facilitative","Trophic")
Stress_type<-c("-,-,-","-,+,-","+,+,+","mixed")
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

#weight by species number
weight=1/(species*4)

#growth rate
C<-rep(0.1,species)
repeat{
  C<-rnorm(n = species,mean = 0.1,sd = 0.001)
  if(sum(C<=0)==0) break
}

repeat{
  C_hold<-rnorm(n = nprey,mean = 0.1,sd = 0.001)
  C_troph<-c(C_hold,rep(0,species-nprey))
  if(sum(C_hold<=0)==0) break
}

for(r in 1:reps){
  #make communities####
  BB<-array(0,dim=c(species, species,length(Com_types)),dimnames=list(1:species,1:species,Com_types))
  
  intra<-rnorm(n = species,mean = -.2, sd = 0.05)
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
  
  for(st in 1:length(Stress_type)){
    print(paste("Rep",r, "-",st))
    for(ct in 1:1){
      #environmental response####
      #start with sensitivity to stress
      sensitivity_scaler<-0.005
      min_sensitivity<-0.3
      if(ct==1){
        #sensitivity<-mvrnorm(n=species,c(0,0),Sigma = matrix(0,2,2)+diag(2)*0.05)
        sensitivity<--matrix(runif(species*nStressors,min = min_sensitivity*sensitivity_scaler,1*sensitivity_scaler),species,nStressors)
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
          sensitivity<-matrix(runif(species*3,min = -sensitivity_scaler,max=sensitivity_scaler),species,3)
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
      X<-array(NA,dim=c(species,nStressors+1,length(Com_types)),dimnames = list(1:species,c("A","B","C","ABC"),Com_types))
      X[,,"No interactions"]<-solve(-BB[,,"No interactions"],C)
      X[,,"Competitive"]<-solve(-BB[,,"Competitive"],C)
      X[,,"Facilitative"]<-solve(-BB[,,"Facilitative"],C)
      X[,,"Trophic"]<-solve(-BB[,,"Trophic"],C_troph)
      
      Xsave<-array(NA,dim=c(length(samp_stress),species,nStressors+1,length(Com_types)),dimnames = list(samp_stress,1:species,c("A","B","C","ABC"),Com_types))
      Xsave[1,,,]<-X
      Env_resp<-matrix(NA,species,nStressors+1,dimnames = list(1:species,c("A","B","C","ABC")))
      
      for(l in 1:Tmax){
        Env_resp[,"A"]<-sensitivity[,1]*stressV[l]
        Env_resp[,"B"]<-sensitivity[,2]*stressV[l]
        Env_resp[,"C"]<-sensitivity[,3]*stressV[l]
        Env_resp[,"ABC"]<-sensitivity[,1]*stressV[l]+sensitivity[,2]*stressV[l]+sensitivity[,3]*stressV[l]
        
        for(com in 1:1){
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
      
      for(com in 1:1){
        hold<-null_compare3(Xsave[,,,com])
        hold$Interactions<-Com_types[com]
        hold$CoTolerance<-c("Random","Positive","Negative")[ct]
        hold$Stress_type<-Stress_type[st]
        hold$Rep<-r
        if(com==4){ 
          for(tlevel in 1:3){
            hold2<-null_compare3(Xsave[,trophV==unique(trophV)[tlevel],,com])
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
#Output$CoTolerance<-factor(Output$CoTolerance,levels=c("Positive","Random","Negative"),ordered = T)
Output$Reversal<-as.numeric(Output$Reversal)

Output_means<-Output%>%
  group_by(Stress,Response,Null_model,Interactions,CoTolerance,Stress_type)%>%
  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)),-Rep)

#Output_trophic$Trophic_level<-factor(Output_trophic$Trophic_level,levels=unique(trophV),ordered = T)
#Output_trophic$CoTolerance<-factor(Output_trophic$CoTolerance,levels=c("Positive","Random","Negative"),ordered = T)
#Output_trophic$Reversal<-as.numeric(Output_trophic$Reversal)

#Output_means_trophic<-Output_trophic%>%
#  group_by(Stress,Response,Null_model,CoTolerance,Stress_type,Trophic_level)%>%
#  summarise_each(funs(mean(.,na.rm=T),lower=quantile(.,probs=c(0.25),na.rm=T),upper=quantile(.,probs=c(0.75),na.rm=T)),-Rep)


Output_means$Response<-factor(Output_means$Response,levels=c("Species richness","Biomass", "Composition"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,levels=c("A_only","B_only","C_only","Actual","Additive","Multiplicative","Species_specific","SS_mod"), ordered=T)
Output_means$Null_model<-factor(Output_means$Null_model,labels=c("A","B","C","ABC","Additive","Multiplicative","Compositional","Mod_comp"))

ColV<-c("black",brewer.pal(3,"Set1"),brewer.pal(3,"Set1")[3])
LineV<-c(4,2,4,1,2,3,4,4)
ggplot(filter(Output_means,Interactions=="No interactions",Response!="Composition",Null_model!="A",Null_model!="B",Null_model!="C", Null_model!="Mod_comp"),aes(x=Stress,y=Change_mean,color=Null_model, fill=Null_model,group=Null_model, linetype=Null_model))+
  geom_ribbon(aes(ymin=Change_lower,ymax=Change_upper),alpha=0.15,color=NA)+
  geom_line(size=1)+
  facet_grid(Stress_type~Response,scale="free")+
  scale_color_manual(values = ColV,name="")+
  scale_fill_manual(values = ColV,name="")+
  scale_linetype_discrete(name="")+
  theme_bw()+
  removeGrid()+
  xlab(expression(paste("Stressor magnitude (",italic("m"),")",sep="")))+
  ylab("Community property relative to control")
ggsave("./Figures/Null model - Figure S1.png", height = 8, width=10)


save(Output_means,file="Multistress3.RData")