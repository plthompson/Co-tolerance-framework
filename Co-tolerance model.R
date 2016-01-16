require(vegan)
require(MASS)
require(ggplot2)
require(dplyr)

nStressors<-2
reps<-20
Co_sensitivityV<-c("Random","Positive","Negative")
Com_types<-c("Competitive","Facilitative","Trophic")
int_strength<-seq(0,1,by=0.25)

#environmental change####
stress_increase<-6000
stress<-65
stressV<-c(seq(0,stress,length=stress*stress_increase+1))
Tmax<-length(stressV)

samp_stress<-seq(1,Tmax,by=stress_increase)

#function for additive null model####
additive_null<-function(A,B){
  SR<-rep(NA,nrow(A))
  bmass_loss<-SR
  bray<-SR
  if(sum(A[1,])>0 | sum(B[1,])>0){
  hold1<-apply(A,2,function(x){(x[1]-x)/x[1]})+apply(B,2,function(x){(x[1]-x)/x[1]})
  hold1[hold1>1]<-1
  bmass_loss<-rowMeans(hold1*rep(A[1,],each=nrow(hold1)))/mean(A[1,])
  
  hold3<-rep(A[1,],each=nrow(hold1))*(1-hold1)
  bray<-c(0,vegdist(hold3)[1:(nrow(hold3)-1)])
  
  SR<-rowSums(A==0 | B == 0)/sum(A[1,]>0)} 
  return(data.frame(SR=SR,Bmass=bmass_loss,Bray=bray))
}


#data storage####
Vary_int.df<-data.frame(SR=NA,Biomass=NA,Bray.Curtis=NA,Stress=stressV[samp_stress], Reps=rep(1:reps,each=length(samp_stress)),Co_sensitivity=factor(rep(Co_sensitivityV,each=length(samp_stress)*reps),levels = c("Positive","Random","Negative"),ordered = T),Int_strength=rep(int_strength,each=length(samp_stress)*reps*length(Co_sensitivityV)),Com_type=rep(Com_types[1:2],each=length(samp_stress)*reps*length(Co_sensitivityV)*length(int_strength)))
Vary_int_troph.df<-data.frame(SR=NA,Biomass=NA,Bray.Curtis=NA,Stress=stressV[samp_stress], Reps=rep(1:reps,each=length(samp_stress)),Co_sensitivity=factor(rep(Co_sensitivityV,each=length(samp_stress)*reps),levels = c("Positive","Random","Negative"),ordered = T),Int_strength=rep(int_strength,each=length(samp_stress)*reps*length(Co_sensitivityV)),Trophic_level=rep(c("Plant","Herbivore","Carnivore"),each=length(samp_stress)*reps*length(Co_sensitivityV)*length(int_strength)))
Vary_int.df2<-Vary_int.df
Vary_int_troph.df2<-Vary_int_troph.df

Raw_effects.df<-data.frame(SR=NA,Biomass=NA,Bray.Curtis=NA,Stress=stressV[samp_stress],Stressors=rep(factor(c("A","B","AB","Net_effect"),levels=c("AB","A","B","Net_effect"),ordered = T),each=length(samp_stress)),Reps=rep(1:reps,each=length(samp_stress)*4),Co_sensitivity=factor(rep(Co_sensitivityV,each=length(samp_stress)*reps*4),levels = c("Positive","Random","Negative"),ordered = T),Int_strength=rep(int_strength,each=length(samp_stress)*reps*length(Co_sensitivityV)*4),Com_type=rep(Com_types[1:2],each=length(samp_stress)*reps*length(Co_sensitivityV)*length(int_strength)*4))
Raw_effects_troph.df<-data.frame(SR=NA,Biomass=NA,Bray.Curtis=NA,Stress=stressV[samp_stress],Stressors=rep(factor(c("A","B","AB"),levels=c("AB","A","B","Net_effect"),ordered = T),each=length(samp_stress)), Reps=rep(1:reps,each=length(samp_stress)*4),Co_sensitivity=factor(rep(Co_sensitivityV,each=length(samp_stress)*reps*4),levels = c("Positive","Random","Negative"),ordered = T),Int_strength=rep(int_strength,each=length(samp_stress)*reps*length(Co_sensitivityV)*4),Trophic_level=rep(c("Plant","Herbivore","Carnivore"),each=length(samp_stress)*reps*length(Co_sensitivityV)*length(int_strength)*4))

for(com in 1:length(Com_types)){
  print(Com_types[com])
  species<-20
  
  for(r in 1:reps){
    print(paste("Rep",r))
    #make communities####
    BB<-array(0,dim=c(species, species),dimnames=list(1:species,1:species))
    
    intra<--.2
    
    repeat{
      
      #competition
      
      b11=-.15
      BB=b11*matrix(runif(species*species),species,species)
      diag(BB)<-intra
      
      #mixed
      #       hold<-matrix(-1,species,species)
      #       int.n<-sum(hold[upper.tri(hold)])*-1
      #       hold[upper.tri(hold)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
      #       hold[lower.tri(hold)][t(hold)[lower.tri(hold)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
      #       BB[,,"Mixed"]<-hold*-BB[,,"Competitive"]
      
      #Facilitative
      if(com==2){
        BB<-BB*-0.1#-0.5
        diag(BB)<-intra
      }
      
      #trophic
      if(com==3){
        nprey<-species*0.5
        nherb<-species*0.3
        npred<-species*0.2
        
        trophV<-factor(c(rep("Plants",nprey),rep("Herbivores",nherb),rep("Predators",npred)),levels = c("Plants","Herbivores","Predators"),ordered = T)
        
        b11=0#-.15
        b12=-0.15#-0.3
        b21=0.1
        b23=-.1
        b32= 0.07#.08
        bdiag1=-.2
        bdiag2=-.2
        
        #tritrophic BB Matrix####
        B11=b11*matrix(runif(nprey*nprey),nprey,nprey)
        B12=b12*matrix(runif(nprey*nherb),nprey,nherb)
        B13=matrix(0,nprey,npred)
        B21=b21*matrix(runif(nherb*nprey),nherb,nprey)
        B22=matrix(0,nherb,nherb)
        B23=b23*matrix(runif(nherb*npred),nherb,npred)
        B31=matrix(0,npred,nprey)
        B32=b32*matrix(runif(npred*nherb),npred,nherb)
        B33=matrix(0,npred,npred)
        BB=rbind(cbind(B11 ,B12, B13),cbind(B21,B22, B23),cbind(B31, B32, B33))
        diag(BB)<-bdiag1
        diag(BB[(nprey+nherb+1):species,(nprey+nherb+1):species])<-bdiag2
      }
      
      #weight by species number
      weight=1/(species*4)
      BBhold<-BB*weight
      
      C<-rep(0.1,species)
      if(Com_types[com]=="Trophic"){
        C<-c(rep(0.1,nprey),rep(0,species-nprey))
      }
      btemp1<-BBhold*int_strength[1]
      diag(btemp1)<-diag(BBhold)
      btemp2<-BBhold*int_strength[2]
      diag(btemp2)<-diag(BBhold)
      btemp3<-BBhold*int_strength[3]
      diag(btemp3)<-diag(BBhold)
      btemp4<-BBhold*int_strength[4]
      diag(btemp4)<-diag(BBhold)
      
      if(sum(c(sum(solve(-BBhold,C)>0)==species,sum(solve(-btemp1,C)>=0)==species,sum(solve(-btemp2,C)>0)==species,sum(solve(-btemp3,C)>0)==species,sum(solve(-btemp4,C)>0)==species))==5) break
    }
    
    #tapply(solve(-BBhold[,,com],C),INDEX = trophV,FUN = sum)
    
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
      
      for(k in 1:length(int_strength)){
        BBrun<-BBhold*int_strength[k]
        diag(BBrun)<-diag(BBhold)
        
        #starting abundances
        X<-array(NA,dim=c(Tmax,species,nStressors+1,length(Com_types)))
        X[1,,,com]<-solve(-BBrun,C)
        
        
        for(l in 1:(Tmax-1)){
          Env_resp<-sensitivity[,1]*stressV[l]+sensitivity[,2]*stressV[l]
          X[l+1,,1,com]<-X[l,,1,com]*exp(C+BBrun%*%X[l,,1,com]+Env_resp)
          X[l+1,,1,com][(X[l+1,,1,com]<0.01)]<-0
          
          Env_resp1<-sensitivity[,1]*stressV[l]
          X[l+1,,2,com]<-X[l,,2,com]*exp(C+BBrun%*%X[l,,2,com]+Env_resp1)
          X[l+1,,2,com][(X[l+1,,2,com]<0.01)]<-0
          
          Env_resp2<-sensitivity[,2]*stressV[l]
          X[l+1,,3,com]<-X[l,,3,com]*exp(C+BBrun%*%X[l,,3,com]+Env_resp2)
          X[l+1,,3,com][(X[l+1,,3,com]<0.01)]<-0
        }
        stressL<-20
        extCol<-c(X[samp_stress[stressL],,1,com]==0)
        extCol[X[samp_stress[stressL],,2,com]==0 | X[samp_stress[stressL],,3,com]==0]<-2
        plot(sensitivity, pch=19, ylim=c(-sensitivity_scaler*1,-sensitivity_scaler*min_sensitivity),xlim=c(-sensitivity_scaler*1,-sensitivity_scaler*min_sensitivity),col=extCol+1, xlab="Tolerance A",ylab="Tolerance B")
        
        if(Com_types[com]=="Trophic"){
          AB_plants<-X[samp_stress,trophV=="Plants",1,com]
          A_plants<-X[samp_stress,trophV=="Plants",2,com]
          B_plants<-X[samp_stress,trophV=="Plants",3,com]
          
          AB_herb<-X[samp_stress,trophV=="Herbivores",1,com]
          A_herb<-X[samp_stress,trophV=="Herbivores",2,com]
          B_herb<-X[samp_stress,trophV=="Herbivores",3,com]
          
          AB_pred<-X[samp_stress,trophV=="Predators",1,com]
          A_pred<-X[samp_stress,trophV=="Predators",2,com]
          B_pred<-X[samp_stress,trophV=="Predators",3,com]
          
          #species richness
          srA_prey<-(nprey-specnumber(A_plants))/nprey
          srB_prey<-(nprey-specnumber(B_plants))/nprey
          srAB_prey<-(nprey-specnumber(AB_plants))/nprey
          
          srA_herb<-(nherb-specnumber(A_herb))/nherb
          srB_herb<-(nherb-specnumber(B_herb))/nherb
          srAB_herb<-(nherb-specnumber(AB_herb))/nherb
          
          srA_pred<-(npred-specnumber(A_pred))/npred
          srB_pred<-(npred-specnumber(B_pred))/npred
          srAB_pred<-(npred-specnumber(AB_pred))/npred
          
          #biomass
          bA_prey<-(sum(X[1,trophV=="Plants",2,com])-rowSums(A_plants))/sum(X[1,trophV=="Plants",2,com])
          bB_prey<-(sum(X[1,trophV=="Plants",3,com])-rowSums(B_plants))/sum(X[1,trophV=="Plants",3,com])
          bAB_prey<-(sum(X[1,trophV=="Plants",1,com])-rowSums(AB_plants))/sum(X[1,trophV=="Plants",1,com])
          
          bA_herb<-(sum(X[1,trophV=="Herbivores",2,com])-rowSums(A_herb))/sum(X[1,trophV=="Herbivores",2,com])
          bB_herb<-(sum(X[1,trophV=="Herbivores",3,com])-rowSums(B_herb))/sum(X[1,trophV=="Herbivores",3,com])
          bAB_herb<-(sum(X[1,trophV=="Herbivores",1,com])-rowSums(AB_herb))/sum(X[1,trophV=="Herbivores",1,com])
          
          bA_pred<-(sum(X[1,trophV=="Predators",2,com])-rowSums(A_pred))/sum(X[1,trophV=="Predators",2,com])
          bB_pred<-(sum(X[1,trophV=="Predators",3,com])-rowSums(B_pred))/sum(X[1,trophV=="Predators",3,com])
          bAB_pred<-(sum(X[1,trophV=="Predators",1,com])-rowSums(AB_pred))/sum(X[1,trophV=="Predators",1,com])
          
          #Bray-Curtis dissimilarity
          bcAB_prey<-c(0,vegdist(AB_plants)[1:(length(samp_stress)-1)])
          bcA_prey<-c(0,vegdist(A_plants)[1:(length(samp_stress)-1)])
          bcB_prey<-c(0,vegdist(B_plants)[1:(length(samp_stress)-1)])
          
          bcAB_herb<-c(0,vegdist(AB_herb)[1:(length(samp_stress)-1)])
          bcA_herb<-c(0,vegdist(A_herb)[1:(length(samp_stress)-1)])
          bcB_herb<-c(0,vegdist(B_herb)[1:(length(samp_stress)-1)])
          
          bcAB_pred<-c(0,vegdist(AB_pred)[1:(length(samp_stress)-1)])
          bcA_pred<-c(0,vegdist(A_pred)[1:(length(samp_stress)-1)])
          bcB_pred<-c(0,vegdist(B_pred)[1:(length(samp_stress)-1)])
          
          
          Vary_int_troph.df$SR[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Plant"]<-srAB_prey-(srA_prey+(1-srA_prey)*srB_prey)
          Vary_int_troph.df$Biomass[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Plant"]<-bAB_prey-(bA_prey+(1-bA_prey)*bB_prey)
          Vary_int_troph.df$Bray.Curtis[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Plant"]<-bcAB_prey-(bcA_prey+(1-bcA_prey)*bcB_prey)
          
          Vary_int_troph.df$SR[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Herbivore"]<-srAB_herb-(srA_herb+(1-srA_herb)*srB_herb)
          Vary_int_troph.df$Biomass[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Herbivore"]<-bAB_herb-(bA_herb+(1-bA_herb)*bB_herb)
          Vary_int_troph.df$Bray.Curtis[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Herbivore"]<-bcAB_herb-(bcA_herb+(1-bcA_herb)*bcB_herb)
          
          Vary_int_troph.df$SR[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Carnivore"]<-srAB_pred-(srA_pred+(1-srA_pred)*srB_pred)
          Vary_int_troph.df$Biomass[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Carnivore"]<-bAB_pred-(bA_pred+(1-bA_pred)*bB_pred)
          Vary_int_troph.df$Bray.Curtis[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Carnivore"]<-bcAB_pred-(bcA_pred+(1-bcA_pred)*bcB_pred)
          
          add.null<-additive_null(A_plants,B_plants)
          DF_plants<-data.frame(SR=srAB_prey-add.null$SR,Biomass=bAB_prey-add.null$Bmass,Bray.Curtis=bcAB_prey-add.null$Bray)
          Vary_int_troph.df2[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Plant",1:3]<-DF_plants
          
          add.null<-additive_null(A_herb,B_herb)
          DF_herb<-data.frame(SR=srAB_herb-add.null$SR,Biomass=bAB_herb-add.null$Bmass,Bray.Curtis=bcAB_herb-add.null$Bray)
          Vary_int_troph.df2[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Herbivore",1:3]<-DF_herb
          
          add.null<-additive_null(A_pred,B_pred)
          DF_pred<-data.frame(SR=srAB_pred-add.null$SR,Biomass=bAB_pred-add.null$Bmass,Bray.Curtis=bcAB_pred-add.null$Bray)
          Vary_int_troph.df2[Vary_int_troph.df$Reps==r & Vary_int_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int_troph.df$Int_strength == int_strength[k] & Vary_int_troph.df$Trophic_level == "Carnivore",1:3]<-DF_pred
          
          Raw_effects_troph.df$SR[Raw_effects_troph.df$Reps==r & Raw_effects_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects_troph.df$Int_strength == int_strength[k]]<-c(srA_prey,srB_prey,srAB_prey,srA_herb,srB_herb,srAB_herb,srA_pred,srB_pred,srAB_pred,DF_plants$SR,DF_herb$SR,DF_pred$SR)
          Raw_effects_troph.df$Biomass[Raw_effects_troph.df$Reps==r & Raw_effects_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects_troph.df$Int_strength == int_strength[k]]<-c(bA_prey,bB_prey,bAB_prey,bA_herb,bB_herb,bAB_herb,bA_pred,bB_pred,bAB_pred,DF_plants$Biomass,DF_herb$Biomass,DF_pred$Biomass)
          Raw_effects_troph.df$Bray.Curtis[Raw_effects_troph.df$Reps==r & Raw_effects_troph.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects_troph.df$Int_strength == int_strength[k]]<-c(bcA_prey,bcB_prey,bcAB_prey,bcA_herb,bcB_herb,bcAB_herb,bcA_pred,bcB_pred,bcAB_pred,DF_plants$Bray.Curtis,DF_herb$Bray.Curtis,DF_pred$Bray.Curtis)
        } else{
          
          A<-X[samp_stress,,2,com]
          B<-X[samp_stress,,3,com]
          AB<-X[samp_stress,,1,com]
          
          #species richness
          srA<-(species-specnumber(A))/species
          srB<-(species-specnumber(B))/species
          srAB<-(species-specnumber(AB))/species
          
          #biomass
          bA<-(sum(A[1,])-rowSums(A))/sum(A[1,])
          bB<-(sum(B[1,])-rowSums(B))/sum(B[1,])
          bAB<-(sum(AB[1,])-rowSums(AB))/sum(AB[1,])
          
          #Bray-Curtis dissimilarity
          bcA<-c(0,vegdist(A)[1:(length(samp_stress)-1)])
          bcB<-c(0,vegdist(B)[1:(length(samp_stress)-1)])
          bcAB<-c(0,vegdist(AB)[1:(length(samp_stress)-1)])
          
          
          Vary_int.df$SR[Vary_int.df$Reps==r & Vary_int.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int.df$Int_strength == int_strength[k] & Vary_int.df$Com_type==Com_types[com]]<-srAB-(srA+(1-srA)*srB)
          Vary_int.df$Biomass[Vary_int.df$Reps==r & Vary_int.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int.df$Int_strength == int_strength[k] & Vary_int.df$Com_type==Com_types[com]]<-bAB-(bA+(1-bA)*bB)
          Vary_int.df$Bray.Curtis[Vary_int.df$Reps==r & Vary_int.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int.df$Int_strength == int_strength[k] & Vary_int.df$Com_type==Com_types[com]]<-bcAB-(bcA+(1-bcA)*bcB)
          
          add.null<-additive_null(A,B)
          DF<-data.frame(SR=srAB-add.null$SR,Biomass=bAB-add.null$Bmass,Bray.Curtis=bcAB-add.null$Bray)
          
          Vary_int.df2[Vary_int.df$Reps==r & Vary_int.df$Co_sensitivity==Co_sensitivityV[ct] & Vary_int.df$Int_strength == int_strength[k] & Vary_int.df$Com_type==Com_types[com],1:3]<-DF
        
          Raw_effects.df$SR[Raw_effects.df$Reps==r & Raw_effects.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects.df$Int_strength == int_strength[k] & Raw_effects.df$Com_type==Com_types[com]]<-c(srA,srB,srAB,srAB-add.null$SR)
          Raw_effects.df$Biomass[Raw_effects.df$Reps==r & Raw_effects.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects.df$Int_strength == int_strength[k] & Raw_effects.df$Com_type==Com_types[com]]<-c(bA,bB,bAB,bAB-add.null$Bmass)
          Raw_effects.df$Bray.Curtis[Raw_effects.df$Reps==r & Raw_effects.df$Co_sensitivity==Co_sensitivityV[ct] & Raw_effects.df$Int_strength == int_strength[k] & Raw_effects.df$Com_type==Com_types[com]]<-c(bcA,bcB,bcAB,bcAB-add.null$Bray)
          }
      }}}}

save(Vary_int.df,Vary_int_troph.df,Vary_int.df2,Vary_int_troph.df2,Raw_effects.df,Raw_effects_troph.df,file = "Co-tolerance.RData")


ColV<-c(2,"Grey30","Dodgerblue")

maxStress<-summarise(group_by(Vary_int.df,Co_sensitivity,Int_strength,Stress,Com_type),Mean=mean(SR,na.rm=T),SD=sd(SR,na.rm=T)) %>%
  group_by(Int_strength,Com_type,Co_sensitivity) %>%
  summarise(stressMax = Stress[which.max(abs(Mean))])

pdf("./Plots/Species richness.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int.df2,Stress,Co_sensitivity,Int_strength,Com_type),Mean=mean(SR,na.rm=T),SD=sd(SR,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity, fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))+
  geom_vline(data=maxStress,aes(xintercept=stressMax), linetype=2)
dev.off()

maxStress<-summarise(group_by(Vary_int.df,Co_sensitivity,Int_strength,Stress,Com_type),Mean=mean(SR,na.rm=T),SD=sd(SR,na.rm=T)) %>%
  filter(Co_sensitivity=="Random") %>%
  group_by(Int_strength,Com_type) %>%
  summarise(stressMax = Stress[which.max(Mean)])

ggplot(summarise(group_by(subset(Raw_effects.df,Co_sensitivity=="Random"),Stress,Int_strength,Com_type,Stressors),Mean=mean(SR,na.rm=T),SD=sd(SR,na.rm=T)),aes(x=Stress,y=Mean, group=Stressors,linetype=Stressors))+
  geom_line(size=1.2)+
  #geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))+
  geom_vline(data=maxStress,aes(xintercept=stressMax), linetype=2)

pdf("./Plots/Composition.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int.df2,Stress,Co_sensitivity,Int_strength,Com_type),Mean=mean(Bray.Curtis,na.rm=T),SD=sd(Bray.Curtis,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity, fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Compositional change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))
dev.off()

maxStress<-summarise(group_by(Vary_int.df,Co_sensitivity,Int_strength,Stress,Com_type),Mean=mean(Bray.Curtis,na.rm=T),SD=sd(Bray.Curtis,na.rm=T)) %>%
  filter(Co_sensitivity=="Random") %>%
  group_by(Int_strength,Com_type) %>%
  summarise(stressMax = Stress[which.max(Mean)])

ggplot(summarise(group_by(subset(Raw_effects.df,Co_sensitivity=="Random"),Stress,Int_strength,Com_type,Stressors),Mean=mean(Bray.Curtis,na.rm=T),SD=sd(Bray.Curtis,na.rm=T)),aes(x=Stress,y=Mean, group=Stressors,linetype=Stressors))+
  geom_line(size=1.2)+
  #geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))+
  geom_vline(data=maxStress,aes(xintercept=stressMax), linetype=2)

pdf("./Plots/Biomass.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int.df2,Stress,Co_sensitivity,Int_strength,Com_type),Mean=mean(Biomass,na.rm=T),SD=sd(Biomass,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity, fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Biomass change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))
dev.off()

maxStress<-summarise(group_by(Vary_int.df,Co_sensitivity,Int_strength,Stress,Com_type),Mean=mean(Biomass,na.rm=T),SD=sd(Biomass,na.rm=T)) %>%
  filter(Co_sensitivity=="Random") %>%
  group_by(Int_strength,Com_type) %>%
  summarise(stressMax = Stress[which.max(Mean)])

ggplot(summarise(group_by(subset(Raw_effects.df,Co_sensitivity=="Random"),Stress,Int_strength,Com_type,Stressors),Mean=mean(Biomass,na.rm=T),SD=sd(Biomass,na.rm=T)),aes(x=Stress,y=Mean, group=Stressors,linetype=Stressors))+
  geom_line(size=1.2)+
  #geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  facet_grid(Com_type~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))+
  geom_vline(data=maxStress,aes(xintercept=stressMax), linetype=2)

pdf("./Plots/Trophic species richness.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int_troph.df2,Stress,Co_sensitivity,Int_strength,Trophic_level),Mean=mean(SR,na.rm=T),SD=sd(SR,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity,fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Trophic_level~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))
dev.off()

pdf("./Plots/Trophic composition.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int_troph.df2,Stress,Co_sensitivity,Int_strength,Trophic_level),Mean=mean(Bray.Curtis,na.rm=T),SD=sd(Bray.Curtis,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity,fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Trophic_level~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Compositional change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))
dev.off()

pdf("./Plots/Trophic biomass.pdf",width = 11,height=8.5)
ggplot(summarise(group_by(Vary_int_troph.df,Stress,Co_sensitivity,Int_strength,Trophic_level),Mean=mean(Biomass,na.rm=T),SD=sd(Biomass,na.rm=T)),aes(x=Stress,y=Mean, group=Co_sensitivity, color=Co_sensitivity,fill=Co_sensitivity))+
  geom_line(size=1.2)+
  geom_ribbon(aes(ymin=Mean-SD,ymax=Mean+SD),alpha=0.3)+
  scale_color_manual(values = ColV)+
  scale_fill_manual(values = ColV)+
  facet_grid(Trophic_level~Int_strength)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Biomass change vs. predicted")+
  scale_y_continuous(breaks=seq(-1,1,by=0.1))
dev.off()
