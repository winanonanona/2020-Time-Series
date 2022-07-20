# Introduction
# Author: Cherlyn Chua
# Last updated: 19 Jul 2022
# FCM stuff


# Set/update working directory (make sure you have at least  
#one folder named as "fcsfiles_IV A and B", where the fcs files are saved 
#one folder named as "background", where the control samples and the samples you would like to check
#one folder named as "plots_IV A_B", where you can save the results)
#https://github.com/rprops/Phenoflow_package/wiki/1.-Phenotypic-diversity-analysis
#Rtools40 for R studio 4.1
setwd("/Users/cherlyn/Desktop/PhenoFlow-master/fcsfiles_IV A and B")


# Packages
#install.packages("changepoint")
#library("tidyverse") 
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("flowAI")
#install_github("rprops/Phenoflow_package")

library(tidyr)
library("Rcpp")
library("devtools")
library("sfsmisc")
library("zoo")
library("changepoint")
library("flowCore")
library("lattice")
library("flowViz")
library("ggplot2") 
library("flowAI")
library("cowplot")
library("scales") # for plotting time series with ggplot
library("openxlsx") # for working with excel files
library("dplyr")
library("flowClean")
library("foreach")
library("Phenoflow")
library("ggrepel")
library("ggpubr")
install.packages("ggplot2")


# Seed for reproducible analysis

#if see memory error
# memory.limit()
#[1] 8117
# memory.limit(size=56000)

set.seed(777)
#1 Import data without or with pattern (there are 2 ways, one is to get all files from one folder, and the other one is set conditions for the sample name)  
fcsfiles <- list.files(path = "/Users/cherlyn/Desktop/URECA/PhenoFlow-master/fcsfiles_IV A and B/total", recursive = TRUE, pattern = ".*STRT.*\\.fcs$",full.names = TRUE)
#process bacteria limit to strt 

#path="/Users/cherlyn/Desktop/URECA/PhenoFlow-master/fcsfiles_IV A and B/my"
flowData <- read.flowSet(files = fcsfiles, transformation = FALSE)
attributes(flowData)

#2. Denoise data
# Asinh transformation (Accuri: signal height)
#fitc is fl1, pe is fl2, vssc is fl4
flowData_transformed <- transform(flowData,`FL1-H`= asinh(`FL1-H`), 
                                  `SSC-H`= asinh(`SSC-H`),
                                  `FL4-H`= asinh(`FL4-H`),
                                  `FSC-H`= asinh(`FSC-H`))

# Select parameters of interest
param <- c("FL1-H","SSC-H","FSC-H","FL4-H")

# Remove original flowset to save space
remove(flowData)

# Extract metadata from sample names (remove first and last characters)
substring(flowCore::sampleNames(flowData_transformed),1)
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed)," ") , rbind)))
View(metadata)

#split station and sampling day
Station <- substr(metadata$X2,1,3)
Sampling_Day <- substr(metadata$X2,4,5)
Dilution <- substr(metadata$X1,1,2)
Replicate <- substr(metadata$X4,1,1)
Stain <- metadata$X3
metadata <- st_dt <- data.frame(Station = Station, Sampling_Day = Sampling_Day, Dilution= Dilution, Stain=Stain, Replicate= Replicate,stringsAsFactors = FALSE)


# Build gates using the background data
# Import background data for gating
bg <- list.files(path = "/Users/cherlyn/Desktop/URECA/PhenoFlow-master/background", recursive = TRUE, pattern = ".fcs$", full.names = TRUE)
BGData <- read.flowSet(files = bg, transformation = FALSE)
BGData_transformed <- transform(BGData,
                                `FL1-H`= asinh(`FL1-H`), 
                                `SSC-H`= asinh(`SSC-H`), 
                                `FL4-H`=asinh(`FL4-H`),
                                `FSC-H`= asinh(`FSC-H`))

# Remove original flowset to save space
remove(BGData)

#Update working directory (just to create a special folder for the results) to allow results generated into specified folder (not to mess up results from the original raw data folder)
setwd("/Users/cherlyn/Desktop/URECA/PhenoFlow-master/plots_IV A_B/total bac re-run")

#Apply the gating to current dataset
# bacteria
sqrcut3 <- matrix(c(13,13,16,16,8,13,15.5,8), ncol = 2, nrow = 4) 
colnames(sqrcut3) <- c("FL1-H","FL4-H")
polyGate1 <- polygonGate(.gate = sqrcut3, filterId = "Total Cells")

#virus
sqrcut3 <- matrix(c(10.5,9.8,13,13,6.8,11,13.5,5.8), ncol = 2, nrow = 4) 
colnames(sqrcut3) <- c("FL1-H","FL4-H")
polyGate1 <- polygonGate(.gate = sqrcut3, filterId = "Total Cells")

# Evaluation of the gates (#run xyplot and save the images manually without g1<- or just run the script below)
g1<-
  xyplot(`FL4-H` ~ `FL1-H`, 
         data = BGData_transformed[c(1,2,3,4)], # adjust sample(s)
         filter = polyGate1, # adjust filter
         scales = list(y = list(limits = c(2,16)), x = list(limits = c(0,16))),
         axis = axis.default, 
         nbin = 125, 
         par.strip.text = list(col = "white", font = 2, cex = 1), 
         smooth = FALSE)
trellis.device(device="tiff", filename="g1.tiff")
print(g1)
dev.off()

g2<-
  xyplot(`FL4-H` ~ `FL1-H`, 
         data = BGData_transformed[c(5,6,7,8)], # adjust sample(s)
         filter = polyGate1, # adjust filter
         scales = list(y = list(limits = c(2,16)), x = list(limits = c(2,16))),
         axis = axis.default, 
         nbin = 125, 
         par.strip.text = list(col = "white", font = 2, cex = 1), 
         smooth = FALSE)
trellis.device(device="tiff", filename="g2.tiff")
print(g2)
dev.off()

g3<-
  xyplot(`FL4-H` ~ `FL1-H`, 
         data = flowData_transformed[c(2,3,4,5)], # copy sample data into background folder before line 78 and then adjust sample(s)
         filter = polyGate1, # adjust filter
         scales = list(y = list(limits = c(2,16)), x = list(limits = c(2,16))),
         axis = axis.default, 
         nbin = 125, 
         par.strip.text = list(col = "white", font = 2, cex = 1), 
         smooth = FALSE)
trellis.device(device="tiff", filename="g3.tiff")
print(g3)
dev.off()

### Isolate only the cellular information based on the polyGate1
### Isolate only the cellular information based on the polyGate1
flowData_transformed <- Subset(flowData_transformed, polyGate1)

summary <- fsApply(x = flowData_transformed, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
maxval <- max(summary[,"FL1-H"]) #Replace with the column representing the green fluorescence channel (e.g. "FITC-H")
mytrans <- function(x) x/maxval
flowData_transformed <- transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
                                  `FL4-H`=mytrans(`FL4-H`), 
                                  `SSC-H`=mytrans(`SSC-H`),
                                  `FSC-H`=mytrans(`FSC-H`))



### Randomly resample to the lowest sample size
#flowData_transformed <- FCS_resample(flowData_transformed, replace=TRUE)
### Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed, param, nbin=128, 
                    bw=0.01,normalize=function(x) x)

#This function allows more accurate diversity and error estimation by bootstrapping from the initial FCS files, as well as denoising through flowAI: automatic and interactive anomaly discerning tools for flow cytometry data. 
### Calculate Diversity from normalized fingerprint (8 mins)
Diversity.fbasis <- Diversity(fbasis,d=3,plot=FALSE, R=999)

Diversity.fbasis2 <- cbind(Diversity.fbasis, metadata)

# Diversity assessment with cleaning (12:43 to 12:57, it takes about 15 mins for 3 bootstraps)
Diversity.clean <- Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3,
                                cleanFCS = TRUE)

#Making plots now
theme_set(theme_bw()+theme(panel.grid.major = element_blank(), 
                           panel.grid.minor = element_blank(), 
                           panel.background = element_rect(colour = "black", size=1, fill=NA), 
                           axis.line = element_line(colour = "black"), 
                           axis.title.x = element_text(size = 14), 
                           axis.text.x = element_text(size = 12, color = "black", margin = unit(c(t = 12, r = 0, b = 0, l = 0), "pt")), 
                           axis.text.y = element_text(size = 12, color = "black", margin = unit(c(t = 0, r = 12, b = 0, l = 0), "pt")), 
                           axis.title.y = element_text(size = 14), title = element_text(size = 14), axis.ticks.length.x = unit(-6, "pt"), 
                           axis.ticks.length.y = unit(-6, "pt")))

Diversity.fbasis2Experiment[Diversity.fbasis2Station=="STL"]<-"STL"
Diversity.fbasis2Experiment[Diversity.fbasis2Station=="SNB"]<-"SNB"

p1 <- ggplot(data = Diversity.fbasis2, aes(x = Station, y = D2, color = Station)) +
  geom_point(size = 8, alpha = 0.7) +
  geom_line() +
  theme_bw() +
  labs( y = "Phenotypic diversity (D2)", x = "Station")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")+ylim(0,4000)
trellis.device(device="tiff", filename="p1.tiff")
print(p1)
dev.off()

plog1 <- ggplot(data = Diversity.fbasis2, aes(x = Station, y = log10(D2), color = Station)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_line() +
  #stat_summary(fun=mean, geom="line", aes(group=Location)) +
  #facet_wrap(~Testbed + Phase, nrow = 2, scales = "free") +
  theme_bw() +
  labs( y = "Phenotypic diversity (log D2)", x = "Station")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")+ylim(2.8,3.6)

trellis.device(device="tiff", filename="plog1.tiff")
print(plog1)
dev.off()

#try box plots
library("readxl")
bacteria_all=read_excel("/Users/cherlyn/Desktop/URECA/PhenoFlow-master/plots_IV A_B/FCM diversity metrics/Bacteria/results.metrics_fbasis2_v2.xlsx")
virus_all =read_excel("/Users/cherlyn/Desktop/URECA/PhenoFlow-master/plots_IV A_B/FCM diversity metrics/Virus All/results.metrics_clean4_d_v1.xlsx")

#perform pairwise comparisons
compare_means(D2 ~ Station, data= virus_all)
my_comparisons <- list(c("SBW","WDL"),c("STL","SNB"), c("SBW", "SNB"), c("WDL","STL"))
barplot <- ggplot(data = bacteria_all, aes(x = Station, y = log10(D2), color = Station), palette="jco") + labs( y = "Phenotypic diversity (log D2)", x = "Station") +geom_boxplot(aes(colour=Station)) +geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")+ylim(2.8,3.6)

#barplot <- ggplot(data = virus_all, aes(x = Station, y = log10(D2), color = Station, add="mean_se" ,palette="jco")) + labs( y = "Phenotypic diversity (log D2)", x = "Station") +geom_boxplot(aes(colour=Station)) +stat_compare_means(comparisons= my_comparisons) + stat_compare_means(label.y=3.7)

trellis.device(device="tiff", filename="barplotvirus.tiff")
print(barplot)
dev.off()

### Calculate Diversity from clean fingerprint 

Diversity.fbasis_c <- Diversity.clean
Diversity.test1<-Diversity.clean
View(Diversity.test1)

#It is important to put as.character before the name to split.
d<- data.frame(do.call(rbind,lapply(strsplit(as.character(Diversity.test1$Sample_names)," "),rbind)))
Station <- substr(d$X2,1,3)
Sampling_Day <- substr(d$X2,4,5)
Dilution <- substr(d$X1,1,2)
Replicate <- substr(d$X4,1,1)
Stain <- d$X3
d <- data.frame(Station = Station, Sampling_Day = Sampling_Day, Dilution= Dilution, Stain=Stain, Replicate= Replicate,stringsAsFactors = FALSE)
View(d)

View(Diversity.test1)
Diversity.fbasis4 <- cbind(Diversity.fbasis_c, d)
#To make sure the split is correct
View(Diversity.fbasis4)
write.csv(Diversity.fbasis4, "results.metrics_clean4_d_v1.csv",row.names=T)
write.xlsx(Diversity.fbasis4, "results.metrics_clean4_d_v1.xlsx")


### Export diversity estimates to .csv file in the chosen directory
write.csv(Diversity.fbasis2, "results.metrics_fbasis2_v2.csv",row.names=T)
write.csv(Diversity.clean, "results.metrics_clean_orig.csv",row.names=T)
write.xlsx(Diversity.fbasis2, "results.metrics_fbasis2_v2.xlsx")
write.xlsx(Diversity.clean, "results.metrics_clean_orig.xlsx")

# Beta-diversity assessment of fingerprint
#Attention: this dimension reduction utilizes the density values as opposed to the count values which are used for beta-diversity analysis in community 16S data sets


# Beta diversity (now all samples can be included, even with relative low cell counts (but do not go below 100 cells))
  
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")
beta.div.co <- data.frame(beta.div[["points"]],beta.div$eig)
colnames(beta.div.co) <- c("Axis1","Axis2","eig")
var <-  round(100*beta.div$eig/(sum(beta.div$eig)),1)



Country <- metadata$Station
metadata<-cbind(metadata,Country)
metadata$Country[metadata$Station=="SBW"] <- "Singapore"
metadata$Country[metadata$Station=="WDL"] <- "Singapore"
metadata$Country[metadata$Station=="SNB"] <- "Malaysia"
metadata$Country[metadata$Station=="STL"] <- "Malaysia"
  


plotdata.beta <- cbind(metadata,beta.div.co)


View(plotdata.beta)
write.csv(plotdata.beta, "results_metrics_plotdatabeta-div_v1.csv",row.names=T)
write.xlsx(plotdata.beta, "results_metrics_plotdatabeta-div_v1.xlsx")
  
  # Plot PCoA
  #library(ggrepel)
  shape_values<-seq(1,23)
  
    pcoa2<-ggplot(data=plotdata.beta) +
    geom_point(aes(x=Axis1, y=Axis2, colour = Country, shape = Station), 
               alpha = 0.9, size = 3, show.legend = TRUE) +
    scale_shape_manual(name="Station", values = c(8,17,19,15,1,4,3,7,2,16))+
    scale_colour_manual(name="Country", values = c("#000000", "#cc6677", "#ddcc77", "#44aa99"))+ 
    labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),
         #title= "PCoA analysis", 
         colour = "Country", shape = "Station")+
    #coord_fixed(ratio = 1)+ 
    stat_ellipse(aes(x=Axis1, y=Axis2,color=Country),type = "norm", level=0.95,linetype = 2, show.legend=F, inherit.aes=T)
  
    
    pcoa2#with ellipse
 
  ggsave(paste0("Beta diversity PCoA plot_4 groups_country.tiff"), pcoa2, width = 6, height = 5, dpi=300) 
  
 
 
 pcoa6<-ggplot(data=plotdata.beta) +
   geom_point(aes(x=Axis1, y=Axis2, colour = Station), 
              alpha = 0.9, size = 3, show.legend = TRUE) +
   scale_colour_manual(values = c("green","#a50f15", "#ffa7e4","blue","#01ACB6","black","#fd8d3c","#4daf4a", "#826600"))+
   labs(x= paste0("Axis1 (",var[1], "%)"), y=paste0("Axis2 (",var[2], "%)"),
        #title= "PCoA analysis", 
        colour = "Station")+
   #coord_fixed(ratio = 1)+
   geom_text_repel(max.overlaps = 100, aes(x=Axis1, y=Axis2,label = Sampling_Day), size = 2, show.legend=F) #showsamplename
# stat_ellipse(aes(x=Axis1, y=Axis2,color=Testbed),type = "norm", level=0.95,linetype = 2, show.legend=F, inherit.aes=T)
 pcoa6#with ellipse
 
 ggsave(paste0("Beta diversity PCoA plot_with sample name_.tiff"), pcoa6, width = 6, height = 5, dpi=300)
 
 
 #diversity clean for alpha diversity diff between diversity clean and diversityfbasis is that clean use bootstrap
 #diversityfbasis for beta diversity 
 
 
 
 
 #Code for exporting the mean and individual D2 value (2nd order Hill number)
 setwd("/Users/cherlyn/Desktop/URECA/PhenoFlow-master/plots_IV A_B/TOTAL VIRUS")
 df<-read.csv("/Users/cherlyn/Desktop/URECA/PhenoFlow-master/plots_IV A_B/TOTAL VIRUS/results.metrics_clean4_d_v1.csv")

 
 #View(df)
 df1<-df
 df1$log<-log(df1$D2,base=10)
 df1$Sample_names<-gsub('.{5}$','',df1$Sample_names) #to remove the last 4 characters
 names(df1)[2] <- "ID"
 df2<-aggregate(df1[,c(5,14)],list(df1$ID),mean)
 colnames(df2)<-c("ID","D2","log_D2")
 
 df3<-df1[,c("ID","D2","log")]
 
#create another column called rep which label the replicates 1 or 2
 df3$rep<-with(df3,ave(ID,ID,FUN = seq_along))
 
 #merge the replicates
 df4<-reshape(df3,timevar="rep",idvar="ID",direction="wide")
 rownames(df4)<-df4$ID

 df4<-df4[,-1]
 df4<-df4[,order(names(df4))]
 df4$ID<-rownames(df4)
 df5<-merge(df4,df2,by="ID")
 
 write.csv(df5,"D2 Summary_virus_total.csv",row.names = F,na="")
 dev.set(dev.next())
 dev.off()
 
 
 
 
 #Bray-curtis matrix
 x=fbasis
 INDICES=NULL
 x <- x@basis/apply(x@basis, 1, max)
 x <- base::round(x, 4)
 if (!is.null(INDICES)) {
   x <- by(x, INDICES = INDICES, FUN = colMeans)
   x <- do.call(rbind, x)
 }


 library(vegan)
 #ADONIS
 adon.results<-adonis(x ~plotdata.beta$Country, method="bray",perm=999)
 results<- data.frame(results = adon.results$aov.tab)
 write.csv(results,file = "adonis_bac_country.csv")
 
print.adonis(adon.results)

 dev.set(dev.next())
 dev.off()
 View(adon.results)
 
 #ANOSIM
anosim.results <- anosim(x, plotdata.beta$Country)
summary(anosim.results)
 write.csv(anosim_results,file="anosim_bac.csv")
 
 #BRAY
 input.dist <- vegan::vegdist(x, method = "bray", binary = FALSE)
 matrix<-as.data.frame(as.matrix(input.dist))
 write.csv(matrix,file = "Bray_Curtis dissimilarity_bac.csv")
 dev.set(dev.next())
 dev.off()
 