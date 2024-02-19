##########################
#MEDIATION ANAL.##########
##########################

library(ggplot2)
library(phyloseq)
library(LDM)

#####Custom functions

phyloseq_2_GG<-function(which.taxa,which.phyloseq,manage.ussigned=NULL,unassign.string=NULL){  
  RESULT<-data.frame(stringsAsFactors = F)
  SSUMS<-sample_sums(which.phyloseq)
  for(i in 1: length(which.taxa))
  {
    actual.phylos<-prune_taxa(taxa_names(which.phyloseq)%in%which.taxa[i],which.phyloseq)
    Abundance<-as.numeric(otu_table(actual.phylos))  
    OTU<-rep(which.taxa[i],length(Abundance))
    Seq.tot<-rep(SSUMS[i],length(Abundance))
    Taxo.actual<-as.character(tax_table(actual.phylos))
    if(manage.ussigned==T){
      for(j in 2: length(Taxo.actual)){
        if(is.na(Taxo.actual[j])){Taxo.actual[j]<-ifelse(regexpr(unassign.string,Taxo.actual[j-1])>0,
                                                         Taxo.actual[j-1],paste(Taxo.actual[j-1],unassign.string,sep=""))}
      }
    }
    
    Taxo.actual.rbind<-do.call("rbind", replicate(length(sample_sums(actual.phylos)), Taxo.actual, simplify = FALSE))
    colnames(Taxo.actual.rbind)<-colnames(tax_table(which.phyloseq))
    Sd<-sample_data(actual.phylos)
    TEMP<-data.frame(Abundance,OTU,Seq.tot,Taxo.actual.rbind,Sd,stringsAsFactors = F)
    RESULT<-rbind(RESULT,TEMP)
  }
  RESULT
}

# adjusting parameters of ggplot object
small.fig<-function(x){x+theme_bw()+
    theme(axis.text = element_text(size=8),
          axis.title = element_text(size=8),
          strip.text = element_text(size=8),
          legend.text = element_text(size=8),
          legend.title = element_text(size=8))}


###############################
#DATA##########################
###############################


#####subsets ~ ASVs devected in > 5% samples
load("/media/kreising/DATA/data/Radka_Janet/DADA_PHYLOSEQ/PHYLOSEQ.final_New_met.R")
PHYLOSEQ.final<-PHYLOSEQ.final_New_met


PHYLOSEQ.final.sub<-prune_samples(!is.na(sample_data(PHYLOSEQ.final)$BMI),PHYLOSEQ.final)
PHYLOSEQ.final.sub.r<-rarefy_even_depth(PHYLOSEQ.final.sub)

PREV<-taxa_sums(transform_sample_counts(PHYLOSEQ.final,function(x) ifelse(x>0,1,0)))/
  nsamples(PHYLOSEQ.final)


PHYLOSEQ.final.sub.hp<-prune_taxa(PREV>0.05,PHYLOSEQ.final.sub)

##############
#ALL SAMPLES##
##############
SD<-data.frame(sample_data(PHYLOSEQ.final.sub.hp))
class(SD)<-"data.frame"
summary(as.factor(SD$Diagnosis))
summary(as.factor(SD$Diagnosis2))
SD$Diagnosis3<-ifelse(SD$Diagnosis=="CTRL","CTRL","narco")

OTU_T<-otu_table(PHYLOSEQ.final.sub.hp)
class(OTU_T)<-"matrix"


MED.all_b <- ldm(formula = OTU_T  ~ Diagnosis3+BMI,
                 data=SD,test.omni3=T,
                 test.mediation = TRUE,
                 seed=1234)

##############
#####NT1######
##############

SUBSET<-prune_samples(sample_data(PHYLOSEQ.final.sub.hp)$Diagnosis%in%c("CTRL","NT_1"),
                      PHYLOSEQ.final.sub.hp)
SUBSET<-prune_taxa(taxa_sums(SUBSET)>0,SUBSET)


SD<-data.frame(sample_data(SUBSET))
class(SD)<-"data.frame"

OTU_T<-otu_table(SUBSET)
class(OTU_T)<-"matrix"

MED.NT1_b <- ldm(formula = OTU_T  ~ Diagnosis+BMI,
                 data=SD,
                 test.omni3=T,
                 test.mediation = TRUE,
                 seed=1234)


##############
#####NT2######
##############

SUBSET<-prune_samples(sample_data(PHYLOSEQ.final.sub.hp)$Diagnosis%in%c("CTRL","NT_2"),
                      PHYLOSEQ.final.sub.hp)
SUBSET<-prune_taxa(taxa_sums(SUBSET)>0,SUBSET)

SD<-data.frame(sample_data(SUBSET))
class(SD)<-"data.frame"

OTU_T<-otu_table(SUBSET)
class(OTU_T)<-"matrix"

MED.NT2_b <- ldm(formula = OTU_T  ~ Diagnosis+BMI,
                 data=SD,
                 test.omni3=T,
                 test.mediation = TRUE,
                 seed=1234)


##############
#####IH######
##############

SUBSET<-prune_samples(sample_data(PHYLOSEQ.final.sub.hp)$Diagnosis%in%c("CTRL","IH"),
                      PHYLOSEQ.final.sub.hp)
SUBSET<-prune_taxa(taxa_sums(SUBSET)>0,SUBSET)

SD<-data.frame(sample_data(SUBSET))
class(SD)<-"data.frame"

OTU_T<-otu_table(SUBSET)
class(OTU_T)<-"matrix"

MED.IH_b <- ldm(formula = OTU_T  ~ Diagnosis+BMI,
                data=SD,
                test.omni3=T,
                test.mediation = TRUE,
                seed=1234)



# #significance of mediator effect
MED.all_b$med.p.global.omni3 
MED.NT1_b$med.p.global.omni3
MED.NT2_b$med.p.global.omni3
MED.IH_b$med.p.global.omni3


#mediator OTUs
MED.all_b$med.detected.otu.omni3
MED.NT1_b$med.detected.otu.omni3
MED.NT2_b$med.detected.otu.omni3
MED.IH_b$med.detected.otu.omni3

# values bivariate
PVALS<-t(MED.NT1$p.otu.omni3)
# MED.all_b$p.otu.omni3.com
unname(PVALS[rownames(PVALS)%in%MED.NT1_b$med.detected.otu.omni3,])
colnames(PVALS)


TT<-tax_table(PHYLOSEQ.final.sub.hp)
unname(TT[taxa_names(TT)==MED.all_b$med.detected.otu.omni3,])
unname(TT[taxa_names(TT)==MED.NT1_b$med.detected.otu.omni3,])


#######################
#plots for NT1 + CTRL##
#######################

# select subset and transform to proportions
SUBSET<-prune_samples(sample_data(PHYLOSEQ.final.sub.hp)$Diagnosis%in%c("CTRL","NT_1"),
                      PHYLOSEQ.final.sub.hp)
SUBSET<-prune_taxa(taxa_sums(SUBSET)>0,SUBSET)
SUBSET.tr<-transform_sample_counts(SUBSET, function(x) x/sum(x))

#data.frame for significant ASVs
FOR_PLOT<-
  phyloseq_2_GG(which.taxa=MED.NT1_b$med.detected.otu.omni3,
                which.phyloseq=SUBSET.tr,manage.ussigned=T,unassign.string="")

#plots for BMI
BMI.p<-
  ggplot(FOR_PLOT, aes(x=BMI,y=Abundance))+geom_point(alpha=0.3,size=0.5)+facet_grid(.~Genus)+
  geom_smooth(method = "lm",colour="black", size=0.7)

#plots for diagnosis (NT1 vs. CTRL)
DIAG.p<-
  ggplot(FOR_PLOT, aes(x=Diagnosis,y=Abundance))+geom_boxplot(outlier.shape = NA)+
  geom_jitter(alpha=0.3,size=0.7)+facet_grid(.~Genus,scales = "free",space="free")



#######################
BMI.p.s<-small.fig(BMI.p)
DIAG.p.s<-small.fig(DIAG.p)

GRID<-ggarrange(BMI.p.s+ggtitle("A)")+theme(plot.title = element_text(size=8,face = "bold")),
                DIAG.p.s+ggtitle("B)")+theme(plot.title = element_text(size=8,face = "bold")),
                nrow = 2)
GRID


ggsave(GRID, filename = "/media/kreising/DATA/data/Radka_Janet/RESULTS_NT1/Mediation_ASVs.pdf",
       width = 8, height = 10, units = "cm")
