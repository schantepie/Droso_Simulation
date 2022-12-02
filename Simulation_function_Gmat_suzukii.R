## Code for D. suzukii G-mat
library(MCMCglmm)
library(nadiv)
library(tidyverse)
setwd('/home/stephane/Dropbox/G_Droso/')
load("data_suzukii_Gmat.RData")
pops=split(pheno_suzukii,pheno_suzukii$pop)
ped.list={}
diag={}

p=26 # number of traits (PCs)

for (j in 1:length(pops)){
# make broken stick priors
  varexplained <- matrix(1/1:p, 2, p, byrow=T)
  for (i in 1:p){
    varexplained[2, i] <- sum(varexplained[1, i:p])/p
  }
  Vtotale <- sum(diag(cov(pops[[j]][,7:(p+6)])))
  variance <- varexplained[2, ] * (Vtotale/2)
  varmat<-as.matrix(variance,p,p)
  diag[[j]]<-diag(c(variance),p,p)
  Gprior <- diag
  Rprior <- diag
  
# get pedigrees  
  ped.list[[j]]=prepPed(pops[[j]][,1:3])
}
 


##############################################################################
############ simulate pedigree according to the Antoine's design
##############################################################################
nbgen=5
NbEgg=100
sexRatio=0.5
nbfounder=2
nbpop=6


pops[[1]]$line=paste0("L",pops[[1]]$dam) # WARNINGS pb with lines of the first pop

library( nadiv)
cpt_progressBar=0
pb <- txtProgressBar(min = 0,      
                     max = nlevels(pheno_suzukii$line),
                     style = 3,   
                     width = 50,  
                     char = "=") 
pedSimul=list()
Ainv=list()

for(P in 1:nbpop) {
       DummyparentIso=list()
              for(iso in levels(factor(pops[[P]]$line))){ #SIMULATE PEDIGREE BY ISOLINE
                        cpt_progressBar=cpt_progressBar+1
                        Dummyparent=list()
                        Dummyparent[[1]]=cbind(paste0(iso,"founder",1:nbfounder),NA,NA,paste0("founder"))
                        for (i in 1:nbgen){
                          nbReproducers=dim(Dummyparent[[i]])[1]
                          nfemale=round(nbReproducers*0.5)
                          nmale=nbReproducers-nfemale
                          if(i==1 & nbfounder==2){# special case to bypass undesired behaviour of the sample function
                            Dummyparent[[i+1]]=cbind(paste0(iso,"generation",i,"_offspring",1:NbEgg),
                                                     Dummyparent[[i]][rep(1,NbEgg),1],
                                                     Dummyparent[[i]][rep(2,NbEgg),1],paste0("gen",i))
                          }else{
                            Dummyparent[[i+1]]=cbind(paste0(iso,"generation",i,"_offspring",1:NbEgg),
                                                     Dummyparent[[i]][sample(x = 1:nfemale,size =  NbEgg,replace = TRUE),1],
                                                     Dummyparent[[i]][sample(x =((nmale+1):nbReproducers), size =  NbEgg,replace = TRUE),1],paste0("gen",i))
                          }
                          if(i==nbgen){ #crop to the number of female used after 5generations of repro and plug to existing pedigree
                            Dummyparent[[i+1]]=Dummyparent[[i+1]][1:2,,drop=FALSE]
                            Dummyparent[[i+1]][1:2,1]=unique(as.matrix(pops[[P]][pops[[P]]$line==iso,2:3]))
                          }
                        }
                        DummyparentIso[[iso]]=as.matrix(do.call(rbind,Dummyparent))
                        setTxtProgressBar(pb,cpt_progressBar) 
              }
    
       pedSimul[[P]]=prunePed(as.data.frame(rbind(as.matrix(do.call(rbind,DummyparentIso)),as.matrix(cbind(pops[[P]][,1:3], generation=paste0("gen",nbgen+1))))),ped.list[[P]]$animal)
       Ainv[[P]]=inverseA(pedSimul[[P]][,1:3], nodes=as.character(ped.list[[P]]$animal))$Ainv
}
# full simulated pedigree is in the list pedSimul[[P]], with P the pop number
# Ainv[[P]] is the inverse relatedness martix of simulated pedergee reduced which contain only animals present in data, with P the pop number

############



#############################################################################################
######### RUN MCMCglmm using the new inverse relatedness matrix (Ainv)  generated ###########
############################################################################################
#exemple for pop 1

pop=1
nbtraits=26

priorBS= list(R = list(V = Rprior[[pop]], nu = nbtraits+0.002),
              G = list(G1 = list(V = Gprior[[pop]], nu = nbtraits+0.002)))

mod <- MCMCglmm(cbind(Comp1,Comp2,Comp3,Comp4,
                      Comp5,Comp6,Comp7,Comp8,
                      Comp9,Comp10,Comp11,Comp12,
                      Comp13,Comp14,Comp15,Comp16,
                      Comp17,Comp18,Comp19,Comp20,
                      Comp21,Comp22,Comp23,Comp24,
                      Comp25,Comp26) ~ (trait-1) + trait:Csize,
                random = ~us(trait):animal,
                rcov=~ us(trait):units,
                family = rep("gaussian",nbtraits),
                ginverse = list(animal =Ainv[[pop]]), 
                data = pops[[pop]],pr=TRUE,
                prior = priorBS,
                burnin=500, nitt=1500, thin=1)

mod <- MCMCglmm(cbind(Comp1,Comp2,Comp3,Comp4,
                      Comp5,Comp6,Comp7,Comp8,
                      Comp9,Comp10,Comp11,Comp12,
                      Comp13,Comp14,Comp15,Comp16,
                      Comp17,Comp18,Comp19,Comp20,
                      Comp21,Comp22,Comp23,Comp24,
                      Comp25,Comp26) ~ (trait-1) + trait:Csize,
                random = ~us(trait):animal,
                rcov=~ us(trait):units,
                family = rep("gaussian",nbtraits),
                pedigree = pedSimul[[pop]][,1:3], 
                data = pops[[pop]],pr=TRUE,
                prior = priorBS,
                burnin=500, nitt=1500, thin=1)

##################################################################################################
############ Function that simulate Breeding values that correspond to null G matrix ##############
##################################################################################################

simul_Null_Gmatrix<-function(ped,G,animalToExtractByPop=NULL){
  
      require(pedigree)
      require(mvtnorm)
   #initialize container and parameter   
      if (length(G)!=length(ped)) stop("pb number pedigree vs G")
      nbpop=length(G)
      nbsample=dim(G[[1]])[1]
      nbTraits=sqrt(dim(G[[1]])[2])
      pb <- txtProgressBar(min = 0,      
                           max = nbsample, 
                           style = 3,   
                           width = 50,  
                           char = "=") 
      F=list()
      sizepedigree=NULL

      for(i in 1:nbpop) {
        if(!all(colnames(ped[[i]][1:3])==c("animal","sire","dam"))) stop("pb in pedigree name")
        ped[[i]]$animal=as.character(ped[[i]]$animal)
        ped[[i]]$sire=as.character(ped[[i]]$sire)
        ped[[i]]$dam=as.character(ped[[i]]$dam)
        nbInd=dim(ped[[i]])[1]
        F[[i]] <- calcInbreeding(ped[[i]][,1:3])
        sizepedigree= c(sizepedigree,dim(ped[[i]])[1])
      }
        #inialise the countainer of breeding values for all pop and all sample
        sizepedigree=sum(sizepedigree)
        fullpedigreeSample=array(NA,c(sizepedigree,nbTraits,nbsample))

    cpt_progressBar=0
    
     for (nsample in 1:nbsample){
       cpt_progressBar=cpt_progressBar+1
       fullpedigree=NULL
      
       for (popNum in 1:nbpop){
       ###### Transform NA to zero in pedigree 
          pedi= ped[[popNum]]
          pedi[is.na(ped[[popNum]])] <-0
       ###### bind all pedigree ,add columns for location of parents in pedigree, inbreeding coef and breeding values
          parent1=NA
          parent2=NA
          BreedingValue=matrix(NA,ncol=nbTraits,nrow=dim(ped[[popNum]])[1])
          colnames(BreedingValue)=paste("BreedingValue",(1:nbTraits),sep="")
          pedi=cbind(pedi[,1:4],popNum,parent1,parent2,F[[popNum]],BreedingValue)
          Founders=which(pedi$sire==0&pedi$dam==0)
          ## predict breeding values of founders
          pedi[Founders,grep("BreedingValue",colnames(pedi))]=rmvnorm(n = length(Founders), mean=rep(0,nbTraits), sigma=matrix(G[[popNum]][nsample,],nbTraits,nbTraits))
          fullpedigree=rbind(fullpedigree,pedi)
       }
      ##### locate parent 
      fullpedigree$parent1=match(fullpedigree$sire,fullpedigree$animal)
      fullpedigree$parent2=match(fullpedigree$dam,fullpedigree$animal)
      ##### store the rows corresponding to founders and predict their breeding values
      Founders_pos=which(fullpedigree$sire==0&fullpedigree$dam==0)
      ##### store the rows corresponding to non-founders
      nonFounders_pos=c(1:dim(fullpedigree)[1])[-Founders_pos]
      ##### Randomize founders among populations
      fullpedigree[Founders_pos,grep("BreedingValue",colnames(fullpedigree))]=fullpedigree[sample(Founders_pos,size = length(Founders_pos),replace = FALSE),grep("BreedingValue",colnames(fullpedigree))]
      ##### Estimate G based on the breeding values of randomized founders
      GAfterRand=matrix(NA,nrow = length(G),ncol=dim(G[[i]])[2])
      for (i in 1:nbpop) {
           if (dim(G[[i]])[2]==1){
             GAfterRand[i,]=var(fullpedigree[fullpedigree$popNum==i & fullpedigree$sire==0 & fullpedigree$dam==0 ,grep("BreedingValue",colnames(fullpedigree))])
           }else{
             GAfterRand[i,]=c(cov(fullpedigree[fullpedigree$popNum==i & fullpedigree$sire==0 & fullpedigree$dam==0 ,grep("BreedingValue",colnames(fullpedigree))]))
           }
      }
    
      ##### Predict breeding value of non-founders according to standard rules of polygenic inheritance 
      # Here special development for pedrigree with discrete generation that allow to speed up computation
       generation=as.character(levels(fullpedigree$generation)[grep("gen",levels(fullpedigree$generation))])
      # Estimate scaling factor of G using the inbreeding coeff of parents
       varscale=matrix(NA,dim(fullpedigree)[1],1)
       varscale[nonFounders_pos]=as.matrix(1-((fullpedigree[fullpedigree[nonFounders_pos,"parent1"],"F[[popNum]]"]+fullpedigree[fullpedigree[nonFounders_pos,"parent2"],"F[[popNum]]"])/2))
       
       for (gen in generation){
        
        ##### Estimate mean BV for all indiv at generation gen 
        fullpedigree[which(fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]=(fullpedigree[fullpedigree[which(fullpedigree$generation==gen),"parent1"],grep("BreedingValue",colnames(fullpedigree))]+
                                                                                          fullpedigree[fullpedigree[which(fullpedigree$generation==gen),"parent2"],grep("BreedingValue",colnames(fullpedigree))])/2
  

        ##### Add the deviation that depend on the GAfterRand of pop each pop and the scaling factor based on F of each individuals
           for (pop in 1:nbpop) {
             fullpedigree[which(fullpedigree$popNum==pop & fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]=
              fullpedigree[which(fullpedigree$popNum==pop & fullpedigree$generation==gen),grep("BreedingValue",colnames(fullpedigree))]+ #retrieve mean BV
                    rmvnorm(n = length(which(fullpedigree$popNum==pop & fullpedigree$generation==gen)), mean=rep(0,nbTraits), sigma=matrix(GAfterRand[pop,],nbTraits,nbTraits))* #draw MVT using GAfterRand to add to mean BV
                            matrix(rep(t(sqrt(as.matrix(varscale[which(fullpedigree$popNum==pop & fullpedigree$generation==gen)]))),nbTraits),ncol=nbTraits) # previously dranw MVT deviation is scaled using mean F of parents (varscale)
           }
    
        }
      fullpedigreeSample[,,nsample]=as.matrix(fullpedigree[,grep("BreedingValue",colnames(fullpedigree))])
      setTxtProgressBar(pb,cpt_progressBar) 
      }

      dimnames(fullpedigreeSample)=list(fullpedigree$animal,colnames(fullpedigree)[grep("BreedingValue",colnames(fullpedigree))])
      random_bv=list()
      for (pop in 1:nbpop) {
                 random_bv[[pop]]=fullpedigreeSample[match(animalToExtractByPop[[pop]],rownames(fullpedigreeSample)),,]
      }
 return(random_bv)
}



G<-list(mod_sok$VCV[,grep(pattern = "animal*",x = colnames(mod_sok$VCV))],
        mod_mon$VCV[,grep(pattern = "animal*",x = colnames(mod_mon$VCV))],
        mod_par$VCV[,grep(pattern = "animal*",x = colnames(mod_par$VCV))],
        mod_sap$VCV[,grep(pattern = "animal*",x = colnames(mod_sap$VCV))],
        mod_tok$VCV[,grep(pattern = "animal*",x = colnames(mod_tok$VCV))],
        mod_wat$VCV[,grep(pattern = "animal*",x = colnames(mod_wat$VCV))]
)
###### Necessary to extract the breeding value of last generation
animalToExtractByPop=list()
for(i in 1:length(pops)) animalToExtractByPop[[i]]=as.character(pops[[i]]$animal)

RandomBV_Droso=simul_Null_Gmatrix(ped = pedSimul,
                                  G=G,animalToExtractByPop =animalToExtractByPop )

# save(list = c("RandomBV_Droso"),file = "~/RandomBV_Droso.Rdata")



Ve<-list(mod_sok$VCV[,grep(pattern = "unit",x = colnames(mod_sok$VCV))],
        mod_mon$VCV[,grep(pattern = "unit",x = colnames(mod_mon$VCV))],
        mod_par$VCV[,grep(pattern = "unit",x = colnames(mod_par$VCV))],
        mod_sap$VCV[,grep(pattern = "unit",x = colnames(mod_sap$VCV))],
        mod_tok$VCV[,grep(pattern = "unit",x = colnames(mod_tok$VCV))],
        mod_wat$VCV[,grep(pattern = "unit",x = colnames(mod_wat$VCV))]
)


RandomPheno=RandomBV_Droso
nbT=sqrt(dim(G[[1]])[2])
for (i in 1:length(G)){
  for (j in 1:dim(G[[1]])[1]){
    RandomPheno[[i]][,,j]= RandomBV_Droso[[i]][,,j]+rmvnorm(n = dim(RandomBV_Droso[[i]][,,j])[1], mean=rep(0,nbT), sigma=matrix(Ve[[i]][j,],nbT,nbT))
  }
}

save(list = c("RandomPheno"),file = "~/NULLmodel_Pheno_Droso.Rdata")

