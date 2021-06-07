# Model for ordination ####
# Markus Bauer



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# A Preparation ################################################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Packages ###
library(here)
library(tidyverse)
library(vegan)

### Start ###
rm(list = ls())
setwd(here("data/processed"))

### Load data ###
sites <- read_csv2("data_processed_sites.csv", col_names = T, na = "na", col_types = 
                     cols(
                       .default = col_double(),
                       id = col_factor(),
                       location = col_factor(),
                       block = col_factor(),
                       plot = col_factor(),
                       position = col_factor(),
                       dataset = col_factor(),
                       side = col_factor(),
                       exposition = col_factor(),
                       phosphorousClass = col_factor(),
                       biotopeType = col_factor(),
                       baykompv = col_factor(),
                       ffh = col_factor(),
                       min8 = col_factor(),
                       min9 = col_factor()
                     )        
)

species <- read_csv2("data_processed_species.csv", col_names = T, na = "na", col_types = 
                       cols(
                         .default = col_double(),
                         name = col_factor()
                       )) %>%  
  pivot_longer(-name, "site", "value") %>%
  pivot_wider(site, name) 
species6510 <- species %>%
  column_to_rownames("site")
boxplot(m)

traits <- read_csv2("data_processed_traits.csv", col_names = T, na = "na", col_types = 
                      cols(
                        .default = col_double(),
                        name = col_factor(),
                        abb = col_factor(),
                        family = col_factor(),
                        rlg = col_factor(),
                        rlb = col_factor(),
                        targetHerb = col_factor(),
                        targetGrass = col_factor(),
                        ffh6510 = col_factor(),
                        ffh6210 = col_factor(),
                        ruderalIndicator = col_factor(),
                        leanIndicator = col_factor(),
                        target = col_factor()
                      )        
)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# B Statistics ################################################################################################################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### 1 Model #####################################################################################

m <- raupcrick(species)
plot(m)

##################################################################################################
# Statistik #################################################################################
##################################################################################################
library(reshape); library(betapart);library(vegan);library(pastecs);library(devtools);library(DSUR.noof)

## Raup-Crick-Funktion nach Chase et al. 2011 Ecosphere

raup_crick=function(spXsite, plot_names_in_col1=TRUE, classic_metric=FALSE, split_ties=TRUE, reps=9999, set_all_species_equal=FALSE, as.distance.matrix=TRUE, report_similarity=FALSE){
  ##expects a species by site matrix for spXsite, with row names for plots, or optionally plots named in column 1.  By default calculates a modification of the Raup-Crick metric (standardizing the metric to range from -1 to 1 instead of 0 to 1). Specifying classic_metric=TRUE instead calculates the original Raup-Crick metric that ranges from 0 to 1. The option split_ties (defaults to TRUE) adds half of the number of null observations that are equal to the observed number of shared species to the calculation- this is highly recommended.  The argument report_similarity defaults to FALSE so the function reports a dissimilarity (which is appropriate as a measure of beta diversity).  Setting report_similarity=TRUE returns a measure of similarity, as Raup and Crick originally specified.  If ties are split (as we recommend) the dissimilarity (default) and similarity (set report_similarity=TRUE) calculations can be flipped by multiplying by -1 (for our modification, which ranges from -1 to 1) or by subtracting the metric from 1 (for the classic metric which ranges from 0 to 1). If ties are not split (and there are ties between the observed and expected shared number of species) this conversion will not work. The argument reps specifies the number of randomizations (a minimum of 999 is recommended- default is 9999).  set_all_species_equal weights all species equally in the null model instead of weighting species by frequency of occupancy.  
  
  ##Note that the choice of how many plots (rows) to include has a real impact on the metric, as species and their occurrence frequencies across the set of plots is used to determine gamma and the frequency with which each species is drawn from the null model	
  
  ##this section moves plot names in column 1 (if specified as being present) into the row names of the matrix and drops the column of names
  if(plot_names_in_col1){
    row.names(spXsite)<-spXsite[,1]
    spXsite<-spXsite[,-1]
  }
  ## count number of sites and total species richness across all plots (gamma)
  n_sites<-nrow(spXsite)
  gamma<-ncol(spXsite)
  
  ##make the spXsite matrix into a pres/abs. (overwrites initial spXsite matrix):
  ceiling(spXsite/max(spXsite))->spXsite
  
  ##create an occurrence vector- used to give more weight to widely distributed species in the null model:
  occur<-apply(spXsite, MARGIN=2, FUN=sum)
  
  ##NOT recommended- this is a non-trivial change to the metric:
  ##sets all species to occur with equal frequency in the null model
  ##e.g.- discards any occupancy frequency information
  if(set_all_species_equal){
    occur<-rep(1,gamma)
  }
  ## determine how many unique species richness values are in the dataset
  ##this is used to limit the number of null communities that have to be calculated
  alpha_levels<-sort(unique(apply(spXsite, MARGIN=1, FUN=sum)))
  
  ##make_null:
  
  ##alpha_table is used as a lookup to help identify which null distribution to use for the tests later.  It contains one row for each combination of alpha richness levels. 
  alpha_table<-data.frame(c(NA), c(NA))
  names(alpha_table)<-c("smaller_alpha", "bigger_alpha")
  col_count<-1
  
  ##null_array will hold the actual null distribution values.  Each element of the array corresponds to a null distribution for each combination of alpha values.  The alpha_table is used to point to the correct null distribution- the row numbers of alpha_table correspond to the [[x]] indices of the null_array.  Later the function will find the row of alpha_table with the right combination of alpha values.  That row number is used to identify the element of null_array that contains the correct null distribution for that combination of alpha levels. 
  null_array<-list()
  
  ##looping over each combination of alpha levels:
  for(a1 in 1:length(alpha_levels)){
  for(a2 in a1:length(alpha_levels)){
      
      ##build a null distribution of the number of shared species for a pair of alpha values:
      null_shared_spp<-NULL
      for(i in 1:reps){
        
        ##two empty null communities of size gamma:
        com1<-rep(0,gamma)
        com2<-rep(0,gamma)
        
        ##add alpha1 number of species to com1, weighting by species occurrence frequencies:
        com1[sample(1:gamma, alpha_levels[a1], replace=FALSE, prob=occur)]<-1
        
        ##same for com2:
        com2[sample(1:gamma, alpha_levels[a2], replace=FALSE, prob=occur)]<-1
        
        ##how many species are shared in common?
        null_shared_spp[i]<-sum((com1+com2)>1)
        
      }
      ##store null distribution, record values for alpha 1 and 2 in the alpha_table to help find the correct null distribution later:
      null_array[[col_count]]<-null_shared_spp
      
      alpha_table[col_count, which(names(alpha_table)=="smaller_alpha")]<-alpha_levels[a1]
      alpha_table[col_count, which(names(alpha_table)=="bigger_alpha")]<-alpha_levels[a2]
      
      #increment the counter for the columns of the alpha table/ elements of the null array
      col_count<-col_count+1
    }
}
  ##create a new column with both alpha levels to match on:
  alpha_table$matching<-paste(alpha_table[,1], alpha_table[,2], sep="_")
  
  #####################
  ##do the test:
  
  ##build a site by site matrix for the results, with the names of the sites in the row and col names:
  results<-matrix(data=NA, nrow=n_sites, ncol=n_sites, dimnames=list(row.names(spXsite), row.names(spXsite)))
  
  ##for each pair of sites (duplicates effort now to make a full matrix instead of a half one- but this part should be minimal time as compared to the null model building)
  for(i in 1:n_sites){
    for(j in 1:n_sites){
      
      ##how many species are shared between the two sites:
      n_shared_obs<-sum((spXsite[i,]+spXsite[j,])>1)
      
      ## what was the observed richness of each site?
      obs_a1<-sum(spXsite[i,])
      obs_a2<-sum(spXsite[j,])
      
      ##place these alphas into an object to match against alpha_table (sort so smaller alpha is first)
      obs_a_pair<-sort(c(obs_a1, obs_a2))
      
      ##match against the alpha table- row index identifies which element of the null array contains the correct null distribution for the observed combination of alpha values:
      null_index<-which(alpha_table$matching==paste(obs_a_pair[1], obs_a_pair[2], sep="_"))
      
      ##how many null observations is the observed value tied with?
      num_exact_matching_in_null<-sum(null_array[[null_index]]==n_shared_obs)
      
      ##how many null values are bigger than the observed value?
      num_greater_in_null<-sum(null_array[[null_index]]>n_shared_obs)
      
      rc<-(num_greater_in_null)/reps
      
      if(split_ties){
        
        rc<-((num_greater_in_null+(num_exact_matching_in_null)/2)/reps)
      }
      
      if(!classic_metric){
        
        ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
        rc<-(rc-.5)*2
      }
      ## at this point rc represents an index of dissimilarity- multiply by -1 to convert to a similarity as specified in the original 1979 Raup Crick paper
      if(report_similarity & !classic_metric){
        rc<- rc*-1
      }
      ## the switch to similarity is done differently if the original 0 to 1 range of the metric is used:
      if(report_similarity & classic_metric){
        rc<- 1-rc
      }
      ##store the metric in the results matrix:
      results[i,j]<-round(rc, digits=2)
      
      
    }
  }
  
  if(as.distance.matrix){
    results<-as.dist(results)
  }	
  
  return(results)
}

## Vorbereitung

mT <- raup_crick(vdataT, plot_names_in_col1 = F, classic_metric = T);mT84 <- raup_crick(vdataT84, plot_names_in_col1 = F);mT93 <- raup_crick(vdataT93, plot_names_in_col1 = F);mT18 <- raup_crick(vdataT18, plot_names_in_col1 = F);
mB <- raup_crick(vdataB, plot_names_in_col1 = F,classic_metric = T);mB03 <- raup_crick(vdataB03, plot_names_in_col1 = F);mB18 <- raup_crick(vdataB18, plot_names_in_col1 = F);
vT <- raupcrick(vdataT,nsimul=9999); vT84 <- raupcrick(vdataT84,nsimul=9999);vT93 <- raupcrick(vdataT93,nsimul=9999);vT18 <- raupcrick(vdataT18,nsimul=9999);
mTtest <- betadisper(mT, edataT$year);mBtest <- betadisper(mB, edataB$year);
vB <- raupcrick(vdataB,nsimul=9999); vB03 <- raupcrick(vdataB03,nsimul=9999); vB18 <- raupcrick(vdataB18,nsimul=9999);

## Statistik

boxplot(mT84,mT93,mT18)
#boxplot(vT84,vT93,vT18)
boxplot(mB03,mB18)
#boxplot(vB03,vB18)
permutest(mTtest);
permutest(mBtest);
plot(TukeyHSD(mTtest));TukeyHSD(mTtest)
plot(TukeyHSD(mBtest));TukeyHSD(mBtest)



##################################################################################################
# Plotten #################################################################################
##################################################################################################
library(ggplot2);library(cowplot)

## Vorbereitung

pdataT <- data.frame(matrix(unlist(mT84),byrow=T,dimnames = list(NULL,c("T84"))),stringsAsFactors=FALSE)
pdataT$T93 <- as.numeric(matrix(unlist(mT93),byrow=T,dimnames = list(NULL,"T93")))
pdataT$T18 <- as.numeric(matrix(unlist(mT18),byrow=T,dimnames = list(NULL,"T18")))
pdataB <- data.frame(matrix(unlist(mB03),byrow=T,dimnames = list(NULL,c("B03"))),stringsAsFactors=FALSE)
pdataB$B18 <- as.numeric(matrix(unlist(mB18),byrow=T,dimnames = list(NULL,"B18")))
pdataT$ID <- as.character(1:780);pdataB$ID <- as.character(1:861)
pdataT <- melt(pdataT, id="ID", measured=c("T84","T93","T18"));pdataB <- melt(pdataB, id="ID", measured=c("B03","B18"));

## Plot

###1984-1993-2018
TGraph <- ggplot(pdataT, aes(variable, value))+
  geom_pointrange(stat = "summary",fun.ymin = function(z) {quantile(z,0.25)},fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, size = 0.4)+
  geom_hline(yintercept = 0,linetype="dashed", size=.3,colour="grey40")+
  scale_x_discrete(labels=c("1984","1993","2018"))+
  scale_y_continuous(limits=c(-1,1), breaks=seq(-1, 1000, 0.5))+
  labs(x="",y="Raup-Crick dissimalirty")+
  theme(axis.text.x  = element_text(size=7), axis.title.y = element_text(size=7), axis.text.y = element_text(size=7), axis.line.y=element_line(linetype = 1), axis.ticks.x = element_blank(),axis.line.x = element_blank(),panel.background=element_rect(fill="white"))
###2003-2018
BGraph <- ggplot(pdataB, aes(variable,value))+
  geom_pointrange(stat = "summary",fun.ymin = function(z) {quantile(z,0.25)},fun.ymax = function(z) {quantile(z,0.75)},
                  fun.y = median, size = 0.4)+
  geom_hline(yintercept = 0,linetype="dashed", size=.3,colour="grey40")+
  scale_x_discrete(labels=c("2003","2018"))+
  scale_y_continuous(limits=c(-1,1), breaks=seq(-1, 1000, 0.5))+
  labs(x="",y="")+
  theme(axis.text.x  = element_text(size=7), axis.title.y = element_text(size=7), axis.text.y = element_text(size=7), axis.line.y=element_line(linetype = 1), axis.ticks.x = element_blank(),axis.line.x = element_blank(),panel.background=element_rect(fill="white"))
###Multiplot
p <- plot_grid(TGraph,BGraph, rel_widths=c(1.2,1), labels=c("(a)","(b)"),
               ncol=2, label_size=10, align="h")+
  theme(plot.margin = unit(c(0,0,-0.6,-0.1), "cm"));p
ggsave("raupcrickDiv_(800dpi_7.5x4cm).tiff", p, width=7.5, height=4, dpi=800, units = "cm",path = "C:/Users/HB/Documents/0_Uni/12_Semester/Masterthesis/3_Aufnahmen_und_Ergebnisse/Ergebnisse")
