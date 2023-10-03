######"Rare and abundant bacterial communities in poplar rhizosphere soils respond differently to genetic effects"

setwd("C:/Users/maple/Desktop/code and data")

###########################
##Fig.1F:Differences in community structure between family groups as demonstrated by CPCoA. 
##Fig.1G: the contribution of genetic effects and soil environmental factors to bacterial community structure.

rm(list=ls())

library(ggthemes)
library(ggplot2)
library(factoextra)
library(ggsci)
library(vegan)
library(amplicon)
library(ggConvexHull)

#Fig.1F
otu<-read.table("ASV_table.txt",header = T)
group<-read.table("env.txt",header = T)
beta_cpcoa(otu,group,dis="bray",groupID = "family",ellipse = F)+geom_point(size=5)+
  geom_convexhull(alpha = 0.5, aes(fill = group))+theme_bw()+
  theme(text = element_text(size = 15,face="bold",family = "serif"))+
  scale_fill_npg()+scale_color_npg()

#Adonis
distance <- vegdist(t(otu), method = 'bray')
ub_dis_table <- as.dist(distance, diag = FALSE, upper = FALSE)
adonis(ub_dis_table~group$group, data=otu, permutations = 999) 
adonis(ub_dis_table~group$family, data=otu, permutations = 999) 

#Anosim 
summary(anosim(distance,group$group,permutations = 999))
summary(anosim(distance,group$family,permutations = 999))

#Fig.1G
library(rdacca.hp)
r<-rdacca.hp(t(otu),group[2:10],method="RDA",type="adjR2")
plot(r)



###########################
##Fig.2 Screening for rare and abundant taxa

#ASV table should first be converted to percentage (relative abundance)

otu_AT <- otu[apply(otu, 1, function(x) min(x)>0.001), ] ### abundant taxa,AT
otu_RT <- otu[apply(otu, 1, function(x) max(x)<0.0001 & sum(x)>0), ] ### rare taxa,RT
otu_MT <- otu[apply(otu, 1, function(x) max(x)<0.001 & min(x)>0.0001), ] ### moderate taxa,MT
otu_CRAT <- otu[apply(otu, 1, function(x) max(x)>0.001 & min(x)<0.0001), ] ### conditionally rare or abundant taxa, CRAT
otu_CRT <- otu[apply(otu, 1, function(x) max(x)<=0.001 & min(x)<0.0001), ]
otu_CRT<-otu_CRT[-which(row.names(otu_CRT)%in%row.names(otu_RT)),]   ### conditionally rare taxa,CRT
otu_CAT <- otu[apply(otu, 1, function(x) min(x)>=0.0001 & max(x)>0.001), ] 
otu_CAT<-otu_CAT[-which(row.names(otu_CAT)%in%row.names(otu_AAT)),] ### conditionally abundant taxa, CAT



###########################
##Fig. 3 A:Ternary plot representing the distribution of ASVs between progenies and their parents.

library(tidyverse)
library(ggtern)
library(reshape2)
library(doBy)

data_clean <- function(otu, design, type=c("relative", "absolute"), threshold=0.00001, times=100){
  
  if (type == "absolute"){
    otu_relative <- apply(otu, 2, function(x){x/sum(x)})
  }else {otu_relative <- otu}
  idx <- rowSums(otu_relative > threshold) >= 1
  otu_threshold <- as.data.frame(otu_relative[idx, ])
  otu_threshold$OTUs <- row.names(otu_threshold)
  otu_longer <- pivot_longer(data=otu_threshold,
                             cols=-OTUs,
                             names_to="SampleID",
                             values_to="value")
  merge_data <- merge(otu_longer, design, by ="SampleID")
  otu <- subset(merge_data, select=c("Group","OTUs","value"))
  otu_mean <-summaryBy(value~OTUs+Group,otu,FUN=mean)
  otu_tern <-dcast(otu_mean,OTUs~Group)
  otu_tern$size <- (apply(otu_tern[2:4], 1, mean))*times
  return(otu_tern)
}

otu<-read.table("ASV_table.txt",header = T)
design<-read.table("metadata.txt",header = T)

col<-c("#00468B","#EC0000","#41B43F","#0099B3","#925E9F","#FDAE91","#AC002A","#B2811E","#40A67D","#9E3B4F","#ACB6B6")
top<-c("Proteobacteria","Actinobacteria","Acidobacteria","Chloroflexi","Patescibacteria","Rokubacteria","Gemmatimonadetes","Bacteroidetes","Firmicutes","Cyanobacteria","Other")
taxonomy<-read.table("taxonomy2.txt",header = T)

otuA<-otu[,c(1:11)]
designA<-design[c(1:11),]
otu_tern <- data_clean(otuA, designA, type="absolute", threshold=0.000001, times=100)
head(otu_tern,n=3)
plotdata<-merge(x=otu_tern,y=taxonomy,by="OTUs",all=F)

for(i in 1:length(plotdata$phylum)){
  if (plotdata$phylum[i]%in%top)
  {plotdata$phylum[i]=plotdata$phylum[i]
  }else{
    plotdata$phylum[i]="Other"
  }
  
}
plotdata$phylum<-factor(plotdata$phylum,levels = top)

p1<-ggtern(data=plotdata, aes(x=AF,y=AM,z=AP1)) +
  geom_point(aes(size=size,col=plotdata$phylum), alpha=0.8, show.legend=T) +
  scale_size(range=c(1, 8)) + geom_mask() +scale_color_manual(values=col)+
  theme_bw()+theme(axis.text=element_blank(),
                   axis.ticks=element_blank())


otuA<-otu[,c(1:8,12:15)]
designA<-design[c(1:8,12:15),]
otu_tern <- data_clean(otuA, designA, type="absolute", threshold=0.000001, times=100)
head(otu_tern,n=3)
plotdata<-merge(x=otu_tern,y=taxonomy,by="OTUs",all=F)

for(i in 1:length(plotdata$phylum)){
  if (plotdata$phylum[i]%in%top)
  {plotdata$phylum[i]=plotdata$phylum[i]
  }else{
    plotdata$phylum[i]="Other"
  }
  
}
plotdata$phylum<-factor(plotdata$phylum,levels = top)

p2<-ggtern(data=plotdata, aes(x=AF,y=AM,z=AP2)) +
  geom_point(aes(size=size,col=plotdata$phylum), alpha=0.8, show.legend=T) +
  scale_size(range=c(1, 8)) + geom_mask() +scale_color_manual(values=col)+
  theme_bw()+theme(axis.text=element_blank(),
                   axis.ticks=element_blank())


otuA<-otu[,c(16:26)]
designA<-design[c(16:26),]
otu_tern <- data_clean(otuA, designA, type="absolute", threshold=0.000001, times=100)
head(otu_tern,n=3)
plotdata<-merge(x=otu_tern,y=taxonomy,by="OTUs",all=F)

for(i in 1:length(plotdata$phylum)){
  if (plotdata$phylum[i]%in%top)
  {plotdata$phylum[i]=plotdata$phylum[i]
  }else{
    plotdata$phylum[i]="Other"
  }
  
}
plotdata$phylum<-factor(plotdata$phylum,levels = top)

p3<-ggtern(data=plotdata, aes(x=BF,y=BM,z=BP)) +
  geom_point(aes(size=size,col=plotdata$phylum), alpha=0.8, show.legend=T) +
  scale_size(range=c(1, 8)) + geom_mask() +scale_color_manual(values=col)+
  theme_bw()+theme(axis.text=element_blank(),
                   axis.ticks=element_blank())


library(cowplot)
plot_grid(p1,p2,p3,ncol = 1,align = "hv")



###########################
####Fig.5 The co-occurrence networks,hub nodes and Network robustness
##Import the ASV table for each group separately!!!


##networks
library(WGCNA)
library(igraph)
library(ggplot2)

#Correlation matrix calculation
#It is necessary to remove the taxa with more 0 values in the group before analysis
occor<-corAndPvalue(t(otu),method="pearson",use="p")
occor.r = occor$cor 
occor.p = occor$p 
occor.r[occor.p>0.01|abs(occor.r)<0.8] = 0
diag(occor.r)<-0
occor.r[is.na(occor.r)]<-0
#edge
sum(abs(occor.r)>0)/2 
#node
sum(colSums(abs(occor.r))>0)
#Export the data and plot using Gephi
write.csv(occor.r,file="network.csv")


## hub nodes

occor.r[abs(occor.r)>0]=1
adjacency_unweight<-occor.r

# Obtain an undirected network with no weights
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph 

#computational node
V(igraph)$degree <- degree(igraph)

set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#integrated data
nodes_list <- data.frame(
  nodes_id = V(igraph)$name, 
  degree = V(igraph)$degree, 
  modularity = V(igraph)$modularity
)
head(nodes_list) 

#Calculate within-module connectivity (Zi) and among-module connectivity (Pi)
source('hub_node.r')

row.names(nodes_list)<-nodes_list$nodes_id
nodes_list<-nodes_list[-1]

zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

zi_pi <- na.omit(zi_pi)   # remove NA
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type,shape=type), alpha = 0.5, size = 2) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)


###Network robustness
#input network matrix, percentage of randomly removed species, and ra of all species
#return the proportion of species remained
rand.remov.once<-function(netRaw, rm.percent, sp.ra, abundance.weighted=T){
  id.rm<-sample(1:nrow(netRaw), round(nrow(netRaw)*rm.percent))
  net.Raw=netRaw #don't want change netRaw
  net.Raw[id.rm,]=0;  net.Raw[,id.rm]=0;   ##remove all the links to these species
  if (abundance.weighted){
    net.stength= net.Raw*sp.ra
  } else {
    net.stength= net.Raw
  }
  
  sp.meanInteration<-colMeans(net.stength)
  
  id.rm2<- which(sp.meanInteration<=0)  ##remove species have negative interaction or no interaction with others
  remain.percent<-(nrow(netRaw)-length(id.rm2))/nrow(netRaw)

  remain.percent
}

rm.p.list=seq(0.05,0.2,by=0.05)
rmsimu<-function(netRaw, rm.p.list, sp.ra, abundance.weighted=T,nperm=100){
  t(sapply(rm.p.list,function(x){
    remains=sapply(1:nperm,function(i){
      rand.remov.once(netRaw=netRaw, rm.percent=x, sp.ra=sp.ra, abundance.weighted=abundance.weighted)
    })
    remain.mean=mean(remains)
    remain.sd=sd(remains)
    remain.se=sd(remains)/(nperm^0.5)
    result<-c(remain.mean,remain.sd,remain.se)
    names(result)<-c("remain.mean","remain.sd","remain.se")
    result
  }))
}

sp.ra<-colMeans(otu) 
sp.ra2<-sp.ra[colSums(abs(occor.r))>0]

Weighted.simu<-rmsimu(netRaw=occor.r, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=T,nperm=100)

Unweighted.simu<-rmsimu(netRaw=occor.r.raw, rm.p.list=seq(0.05,1,by=0.05), sp.ra=sp.ra2, abundance.weighted=F,nperm=100)

dat1<-data.frame(Proportion.removed=rep(seq(0.05,1,by=0.05),2),rbind(Weighted.simu,Unweighted.simu),
                 weighted=rep(c("weighted","unweighted"),each=20))
currentdat<-dat1
ggplot(currentdat[currentdat$weighted=="unweighted",], aes(x=Proportion.removed, y=remain.mean)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()

ggplot(currentdat[currentdat$weighted=="weighted",], aes(x=Proportion.removed, y=remain.mean)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  xlab("Proportion of species removed")+
  ylab("Proportion of species remained")+
  theme_light()



###########################
#####Fig. 6 A: Dispersal-niche Continuum Index between groups 
library(DNCImper)
otu<-read.table("ASV_table.txt",header = T)
otu[otu!=0]=1
group<-c(1,1,1,1,2,2,2,2,3,3,3,4,4,4,4,5,5,5,5,6,6,6,7,7,7,7,8,8,8,8)

DNCI_multigroup (otu,group, id = “no_name”, Nperm = 1000, count = TRUE,
                  symmetrize = FALSE, plotSIMPER = TRUE)



#####Fig. 6 B:Community assembly analysis NTI, BNTI and RC_bray
library(picante)
## Import ASV tables and trees
comun<-read.table("ASV_table.txt",header = T)
phy<-read.tree("rooted_tree.tre")
comun<-t(comun)
prune_tree<-prune.sample(comun,phy)
phydist<-cophenetic(prune_tree)
##Calculate the NTI
mntd<-ses.mntd(comun,phydist,null.model="taxa.labels",abundance.weighted=T, runs=999)
##The observed βMNTD
comdist.resultS12<-comdistnt(comun,phydist)
##nullmodel
f<-function(){
  g<-randomizeMatrix(comun,null.model=c("frequency","richness","independentswap","trialswap"),iterations=1000)
  fc<-comdist.result<-comdistnt(g,phydist)
}
mt<-mean(replicate(999, f()))
mt.sd<-sd(replicate(999, f()))
BNTI=(comdist.resultS12-mt)/mt.sd


##RC_bray
raup_crick= function(comun, reps=999){
  require(ecodist) 
  
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(comun)
  gamma<-ncol(comun)
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(comun), row.names(comun)))
  ##make the comun matrix into a new, pres/abs. matrix:
  ceiling(comun/max(comun))->comun.inc
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(comun.inc, MARGIN=2, FUN=sum)
  ##create an abundance vector- used to give more weight to abundant species in the second step of the null model:
  abundance<-apply(comun, MARGIN=2, FUN=sum)
  ##make_null:
  ##looping over each pairwise community combination:
  for(null.one in 1:(nrow(comun)-1)){
    for(null.two in (null.one+1):nrow(comun)){
      null_bray_curtis<-NULL
      for(i in 1:reps){
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        ##add observed number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, sum(comun.inc[null.one,]), replace=FALSE, prob=occur)]<-1
        com1.samp.sp = sample(which(com1>0),(sum(comun[null.one,])-sum(com1)),replace=TRUE,prob=abundance[which(com1>0)]);
        com1.samp.sp = cbind(com1.samp.sp,1); # head(com1.samp.sp);
        com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum)); colnames(com1.sp.counts) = 'counts'; # head(com1.sp.counts);
        com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts)); # head(com1.sp.counts);
        com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts; # com1;
        #sum(com1) - sum(spXsite[null.one,]); ## this should be zero if everything work properly
        rm('com1.samp.sp','com1.sp.counts');      
        ##same for com2:
        com2[sample(1:gamma, sum(comun.inc[null.two,]), replace=FALSE, prob=occur)]<-1
        com2.samp.sp = sample(which(com2>0),(sum(comun[null.two,])-sum(com2)),replace=TRUE,prob=abundance[which(com2>0)]);
        com2.samp.sp = cbind(com2.samp.sp,1); # head(com2.samp.sp);
        com2.sp.counts = as.data.frame(tapply(com2.samp.sp[,2],com2.samp.sp[,1],FUN=sum)); colnames(com2.sp.counts) = 'counts'; # head(com2.sp.counts);
        com2.sp.counts$sp = as.numeric(rownames(com2.sp.counts)); # head(com2.sp.counts);
        com2[com2.sp.counts$sp] = com2[com2.sp.counts$sp] + com2.sp.counts$counts; # com2;
        # sum(com2) - sum(spXsite[null.two,]); ## this should be zero if everything work properly
        rm('com2.samp.sp','com2.sp.counts');
        null.comun = rbind(com1,com2); # null.comun;
        ##calculate null bray curtis
        null_bray_curtis[i] = distance(null.comun,method='bray-curtis');
      }; # end reps loop
      ## empirically observed bray curtis
      obs.bray = distance(comun[c(null.one,null.two),],method='bray-curtis');
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null = sum(null_bray_curtis==obs.bray);
      ##how many null values are smaller than the observed *dissimilarity*?
      num_less_than_in_null = sum(null_bray_curtis<obs.bray);
      rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
      ##modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
      rc = (rc-.5)*2
      results[null.two,null.one] = round(rc,digits=2); ##store the metric in the results matrix
      print(c(null.one,null.two,date()));
    }; ## end null.two loop
  }; ## end null.one loop
  
  results<-as.dist(results)
  return(results)
}
