require(vegan)
require(MASS)

species<-40
nStressors<-2
reps<-100
co_toleranceV<-c("Random","Positive","Negative")
Com_types<-c("No interactions","Competitive","Mixed","Foodweb")

#environmental change####
stress_increase<-8000
stress<-10
stressV<-seq(0,stress,length=stress*stress_increase)
Tmax<-length(stressV)

samp_stress<-seq(1,Tmax,length=100)

#data storage####
Co_tolerance_df<-data.frame(SR=NA,Biomass=NA,Bray.Curtis=NA,Stress=stressV[samp_stress], Reps=rep(1:reps,each=length(samp_stress)),Co_tolerance=rep(co_toleranceV,each=length(samp_stress)*reps),Com=rep(factor(Com_types[1:2],ordered = T,levels = Com_types[1:2]),each=length(samp_stress)*reps*length(co_toleranceV)))

for(r in 1:reps){
print(paste("Rep",r))
#make communities####
BB<-array(0,dim=c(species, species, length(Com_types)),dimnames=list(1:species,1:species,Com_types))

#no interactions
diag(BB[,,"No interactions"])<--.2

#competition
b11=-.15
BB[,,"Competitive"]=b11*matrix(runif(species*species),species,species)
diag(BB[,,"Competitive"])<--.2

#mixed
hold<-matrix(-1,species,species)
int.n<-sum(hold[upper.tri(hold)])*-1
hold[upper.tri(hold)][sample(int.n, replace=F)<=(0.35*int.n)]<-0.5
hold[lower.tri(hold)][t(hold)[lower.tri(hold)]>0][sample(0.35*int.n, replace=F)<(0.10*int.n)]<-0.5
BB[,,"Mixed"]<-hold*-BB[,,"Competitive"]

#no foodweb at the moment

#weight by species number
weight=1/(species*4)
BB<-BB*weight

#starting abundances
X<-array(rep(round(rnorm(n = species,mean = 5,sd = 2)),each=Tmax),dim=c(Tmax,species,nStressors+1,length(Com_types)))
X[X<=0]<-1

for(ct in 1:3){
#environmental response####
#start with sensitivity to stress
sensitivity_scaler<-0.004
if(ct==1){
  sensitivity<-mvrnorm(n=species,c(0,0),Sigma = matrix(0,2,2)+diag(2)*0.05)
} else {
  sensitivity<-mvrnorm(n=species,c(0,0),Sigma = matrix(0.9,2,2)+diag(2)*0.05)
}
if(ct==3){
  sensitivity[,2]<-sensitivity[,2]*-1
}
sensitivity<-(0.5+sensitivity/max(abs(sensitivity))/2)*-sensitivity_scaler
plot(sensitivity, pch=19, ylim=c(-sensitivity_scaler,0),xlim=c(-sensitivity_scaler,0))

#model####
for(com in 1:2){
  #intrinsic rate of increase
  C<-0.125
  C=-BB[,,com]%*%X[1,,1,com]
  for(l in 1:(Tmax-1)){
    Env_resp<-sensitivity[,1]*stressV[l]+sensitivity[,2]*stressV[l]
    X[l+1,,1,com]<-X[l,,1,com]*exp(C+BB[,,com]%*%X[l,,1,com]+Env_resp)
    X[l+1,,1,com][(X[l+1,,1,com]<0.01)]<-0
    
    Env_resp1<-sensitivity[,1]*stressV[l]
    X[l+1,,2,com]<-X[l,,2,com]*exp(C+BB[,,com]%*%X[l,,2,com]+Env_resp1)
    X[l+1,,2,com][(X[l+1,,2,com]<0.01)]<-0
    
    Env_resp2<-sensitivity[,2]*stressV[l]
    X[l+1,,3,com]<-X[l,,3,com]*exp(C+BB[,,com]%*%X[l,,3,com]+Env_resp2)
    X[l+1,,3,com][(X[l+1,,3,com]<0.01)]<-0
  }
#species richness
srA<-(species-specnumber(X[samp_stress,,2,com]))/species
srB<-(species-specnumber(X[samp_stress,,3,com]))/species
srAB<-(species-specnumber(X[samp_stress,,1,com]))/species

#biomass
bA<-(sum(X[1,,2,com])-rowSums(X[samp_stress,,2,com]))/sum(X[1,,2,com])
bB<-(sum(X[1,,3,com])-rowSums(X[samp_stress,,3,com]))/sum(X[1,,3,com])
bAB<-(sum(X[1,,1,com])-rowSums(X[samp_stress,,1,com]))/sum(X[1,,1,com])

#Bray-Curtis dissimilarity
bcAB<-c(0,vegdist(X[seq(1,Tmax,length=100),,1,com])[1:99])
bcA<-c(0,vegdist(X[seq(1,Tmax,length=100),,2,com])[1:99])
bcB<-c(0,vegdist(X[seq(1,Tmax,length=100),,3,com])[1:99])

Co_tolerance_df$SR[Co_tolerance_df$Reps==r & Co_tolerance_df$Co_tolerance==co_toleranceV[ct] & Co_tolerance_df$Com == Com_types[com]]<-srAB-(srA+(1-srA)*srB)
Co_tolerance_df$Biomass[Co_tolerance_df$Reps==r & Co_tolerance_df$Co_tolerance==co_toleranceV[ct] & Co_tolerance_df$Com == Com_types[com]]<-bAB-(bA+(1-bA)*bB)
Co_tolerance_df$Bray.Curtis[Co_tolerance_df$Reps==r & Co_tolerance_df$Co_tolerance==co_toleranceV[ct] & Co_tolerance_df$Com == Com_types[com]]<-bcAB-(bcA+(1-bcA)*bcB)
}}}

require(ggplot2)
require(dplyr)


pdf("Species richness change vs. predicted.pdf",width = 11,height = 8.5)
ggplot(summarise(group_by(Co_tolerance_df,Stress,Co_tolerance,Com),Mean=mean(SR),SD=sd(SR)),aes(x=Stress,y=Mean, group=Co_tolerance, color=Co_tolerance))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0.1)+
  facet_grid(.~Com)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Species richness change vs. predicted")
dev.off()

pdf("Biomass change vs. predicted.pdf",width = 11,height = 8.5)
ggplot(summarise(group_by(Co_tolerance_df,Stress,Co_tolerance,Com),Mean=mean(Biomass),SD=sd(Biomass)),aes(x=Stress,y=Mean, group=Co_tolerance, color=Co_tolerance))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0.1)+
  facet_grid(.~Com)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Biomass change vs. predicted")
dev.off()

pdf("Compositional change vs. predicted.pdf",width = 11,height = 8.5)
ggplot(summarise(group_by(Co_tolerance_df,Stress,Co_tolerance,Com),Mean=mean(Bray.Curtis),SD=sd(Bray.Curtis)),aes(x=Stress,y=Mean, group=Co_tolerance, color=Co_tolerance))+
  geom_line(size=2)+
  geom_errorbar(aes(ymin=Mean-SD,ymax=Mean+SD),width=0.1)+
  facet_grid(.~Com)+
  theme_bw(base_size = 15)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_abline(intercept=0,slope=0, linetype=2)+
  ylab("Compositional change vs. predicted")
dev.off()  
