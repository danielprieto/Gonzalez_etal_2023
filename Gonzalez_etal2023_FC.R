###############################################
# Flow cytometry analysis script  
# date:2023-05-30         
# author:Daniel Prieto dprieto(at)fcien.edu.uy
###############################################
#Load required libraries
library(flowCore)
library(flowDensity)
library(flowViz)
library(flowStats)
library(flowAI)
library(ggcyto)
#Open files
setwd <- ("~/Documentos/FCyt/Gonzalez_et_al2023/")
files.test1 = list.files("~/Documentos/FCyt/Gonzalez_et_al2023/data/", all.files = F, full.names = TRUE)#Read all the files in the directory
fs <- read.flowSet(files.test1[1:4], transformation = F, alter.names =T, emptyValue=FALSE)#Read FCS files and load them into a flowSet object
####################
#QC & data cleaning#
####################
flow_auto_qc(fs, html_report = T, mini_report = FALSE, fcs_QC = TRUE, fcs_highQ=TRUE, folder_results = "~/Documentos/FCyt/Gonzalez_et_al2023/dataHQ/")
#Load QC-filtered files to files.hq
files.hq = list.files("~/Documentos/FCyt/Gonzalez_et_al2023/dataHQ/", all.files = F, full.names = TRUE)#Lee todos los archivos del directorio
fs1 <- read.flowSet(files.hq, transformation = F, alter.names =T, emptyValue=FALSE)#[Read QC-filtered FCS files and load them into a flowSet object
#
plotDens (fs1[[1]], c("FSC.A" ,"SSC.A"))#Plot first flowFrame within the flowSet just to check
################
#Transform data#
################
bx <- logicleTransform()#Define data transformation
bxlist <- transformList(c("FSC.A", "SSC.A", "FSC.H", "VL1.A", "VL2.A", "BL1.A"), bx)#Select channels to transform
datostrans <- transform(fs1, bxlist)#Apply transformation
#####################################
#Clean dataset from debris (fcs/ssc)#
#####################################
clean.dt <- datostrans
for (i in 1:3) { # Loop over the length of the flowSet
  f <- datostrans[[i]]
  # First restrict the FSC - A values :
  fsc.indices <- intersect (which (exprs (f)[, "FSC.A"] < 5) , which (exprs (f) [, "FSC.A"] > 3))
  # Then restrict SSC - A values and intersect with FSC - A restriction above :
  ssc.indices <- intersect (which ( exprs ( f)[, "SSC.A"] > 3) ,
                            which ( exprs (f)[, "SSC.A"] < 5) )
  non.debris.indices <- intersect ( fsc.indices , ssc.indices )
  # Now only select the non debris cells and place the cleaned up flowFrame into the flowSet :
  f.clean <- f[non.debris.indices]
  clean.dt [[i]] <- f.clean
}
plotDens (clean.dt[[2]], c("FSC.A" ,"SSC.A"))#Plot to check data-cleaning results
#Define names
pData(clean.dt)$name <- c("+8Br-cGMP 1 min", "+8Br-cGMP 5 min", "-8Br-cGMP")
############################
#Set visual style for plots#
############################
library(effects)#Load library to set colorblind
lattice::trellis.par.set(effectsTheme(col="colorblind"))#Set pallette to colorblind
flowViz.par.set(theme = trellis.par.get())#Apply palette to flowViz graphics
#############
#Check plots#
#############
#VL1-CFP, VL2-FRET, BL1-YFP
plotDens (clean.dt[[1]], c("FSC.A" ,"SSC.A"))#
plotDens (clean.dt[[1]], c ("VL1.A", "FSC.A"))
autoplot(clean.dt[[1]]) + labs_cyto("marker")
##
############################
#Select cells elipsoid gate#
############################
cells <- lymphGate(clean.dt, channels=c("SSC.A", "FSC.A"), scale=1.6, bwFac = 1, plot = F)#Create data-driven gate
xyplot(`SSC.A` ~ `FSC.A`, clean.dt, filter = cells, xlim = c(3,6), ylim = c(3, 6), smooth = F)#Show gated data in plots
celulas <- Subset(clean.dt, cells)#Apply gate
#
#################################
#Singlet filtering-elipsoid gate#
#################################
singfilt1 <- lymphGate(celulas, channels = c("FSC.H", "FSC.A"), scale=4, bwFac = 2, 
                       filterId = "singGate", evaluate = T, plot =F)#Create data-driven singlet gate
xyplot(`FSC.A` ~ `FSC.H`, celulas, filter = singfilt1, xlim = c(3.5,5.5), ylim = c(3, 5), smooth = F)#Plot gate
sing <- Subset(celulas, singfilt1)#Apply singlet gate
#
######################################
##Normalize (limit) events displayed
#####################################
limit <- function(frame, limit = 10000) {
  frame@exprs <- frame@exprs[1:min(limit, nrow(frame@exprs)), ]
  frame
}#Create function limit
filtset1 <- fsApply(sing, limit, limit = 8000)#Filter dataset and srt limit to 8K
#################################
#Gate populations               #
#################################
CFPposFilt <- rangeGate(sing, "BL1.A", plot=T, refLine=0)#Define data-driven gate CFP+
YFPposFilt <- rangeGate(sing, "VL1.A", plot=T, refLine=0)#Define data-driven gate YFP+
CFPpos <- Subset(sing, CFPposFilt)#Apply CFP+ gate onto singlets (event number normalized)
Doublepos <- Subset(CFPpos, YFPposFilt)#Apply YFP+ gate onto CFP+-gated cells
##############################################
#Evaluate populations with data-driven filter#
##############################################
dposFilt <- curv2Filter("VL1.A", "BL1.A", filterId="data-driven filter", bwFac = 2)#Create a data-driven filter object that selects high-density regions in two dimensions.
dpos <- filter(sing, dposFilt)#Apply gate
dposdata <- split(sing, dpos)#Split area data
#
xyplot(`BL1.A` ~ `VL1.A`, data = filtset1, filter=dposFilt, smooth=F,  
       par.settings = list(panel.background=list(col="transparent"), axis.text=list(col="black"), 
                           strip.border=list(col="black"), axis.line=list(col="black"), 
                           strip.background=list(col="transparent")), xlim=c(0, 4), ylim=c(0, 3), 
       stats = T, xlab="CFP-A", ylab="YFP-A")#Plot gate
#
FRETposFilt <- rangeGate(Doublepos, "VL2.A", plot=T, refLine=0)#Define data-driven gate FRET+
########################################
#Evaluate double positives by quadrants#
########################################
qg <- quadrantGate(sing, c("VL1.A", "BL1.A"), plot=F, sd = c(-2.5, -0.5))#Create quadrant gate
qfs <- split(sing, qg)#split quadrant data
xyplot(`BL1.A` ~ `VL1.A`, data = filtset1, filter=qg, smooth=F, 
       par.settings = list(panel.background=list(col="transparent"), 
                           axis.text=list(col="black"), strip.border=list(col="black"), 
                           axis.line=list(col="black"), strip.background=list(col="transparent")), 
       xlim=c(0, 4), ylim=c(0, 3), stats = F)#Plot gate

#################################
#FRET-positive counting filter  #
#################################
FRETpos <- filter(Doublepos, FRETposFilt)#Apply FRET+ filter
FRETpos <- Subset(Doublepos, subset = FRETposFilt)
summary(FRETpos)#Retrieve descriptive stats from CFP+YFP+FRET+ cells
#
#####################################
##Normalize (limit) events displayed
#####################################
limit <- function(frame, limit = 10000) {
  frame@exprs <- frame@exprs[1:min(limit, nrow(frame@exprs)), ]
  frame
}#Create limit function
filtset <- fsApply(Doublepos, limit, limit = 15000)#Filter dataset and srt limit to 8K
##
#####################################
#Trellis plots for the whole flowSet#
#####################################
#evaluation by c2f
xyplot(`BL1.A` ~ `VL1.A`, data = sing, filter=dposFilt, smooth=F, 
       par.settings = list(panel.background=list(col="transparent"), axis.text=list(col="black"), 
                           strip.border=list(col="black"), axis.line=list(col="black"), strip.background=list(col="transparent")), 
       xlim=c(0, 4), ylim=c(0, 3))
#
###########################
#Comparative density plots#
###########################
#gated with c2f
densityplot (~ `VL2.A`,  dposdata$`area 2`, filter = NULL, layout = c(1,1), xlim = c(2.7,4), 
             xlab=("FRET-A"), stack = T, 
             refline=median(dposdata$`area 2`@frames$`FRET Y-C cGMP 2021-08-31_Group_sin tratar msimos ISTRUE.fcs`[,5]@exprs))
#EOF