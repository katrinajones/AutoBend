####Plots for manuscript validating in silico methods using cat and tegu data
####11.24.2020

library(dplyr)
library(tidyr)
library(plyr)
library(myfunctions)
library(abind)
library(ggplot2)
#library(ggtern)#only load when needed as messes up ggplot. clear environment after.
library(viridis)
library(RColorBrewer)
library(phytools)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(gridGraphics)
library(ggtree)
library(ggstance)
source('G:/My Drive/DOCUMENTS/WORK/RESEARCH/POSTDOC WORK/SYNAPSID PROJECT/MANUSCRIPTS/IN SILICO BENDING/VALIDATION/In silico modelling_R/R scripts/Functions.R')

#########################################################################################
#Load data from Maya
########################################################################################

###########################################################
#REAL BIOMECHANICAL DATA
###########################################################

#Load real biomechanical data
catbiom<-read.csv(file=c("./Data//Biomechanical/Cat Biomech Real.csv"))
catbiom$JOINT<-as.factor(gsub("27--28", "27--S", catbiom$JOINT))
catbiom$JOINT<-factor(catbiom$JOINT,levels(catbiom$JOINT)[c(19:25,1:18)])
#catbiom<-catbiom[-which(catbiom$JOINT=="3--4"|catbiom$JOINT=="4--5"|catbiom$JOINT=="5--6"|catbiom$JOINT=="6--7"|catbiom$JOINT=="7--8"),]
catbiom<-data.frame(file=rep("real", nrow(catbiom)), Species=rep("Cat", nrow(catbiom)),catbiom, Sag=catbiom$VENTRO_A+catbiom$DORSO_A, Latero2=catbiom$LATERO_A*2, 
                    Tor2=catbiom$TORSION_A*2, DV_Stiff=rowMeans(catbiom[,c("DORSO_B", "VENTRO_B")]), jointno=as.numeric(catbiom$JOINT))
tegubiom<-read.csv(file=c("./Data/Biomechanical/Tegu Biomech Real.csv"))
tegubiom$JOINT<-as.factor(gsub("-", "--", tegubiom$JOINT))
tegubiom$JOINT<-factor(tegubiom$JOINT,levels(tegubiom$JOINT)[c(17:23,1:16)])
#tegubiom<-tegubiom[-which(tegubiom$JOINT=="3--4"|tegubiom$JOINT=="4--5"|tegubiom$JOINT=="5--6"|tegubiom$JOINT=="6--7"|tegubiom$JOINT=="7--8"),]
tegubiom<-data.frame(file=rep("real", nrow(tegubiom)), Species=rep("Tegu", nrow(tegubiom)),tegubiom, Sag=tegubiom$VENTRO_A+tegubiom$DORSO_A, Latero2=tegubiom$LATERO_A*2, 
                     Tor2=rep(0, nrow(tegubiom)), DV_Stiff=rowMeans(tegubiom[,c("DORSO_B", "VENTRO_B")]),jointno=as.numeric(tegubiom$JOINT))
realdat<-list(Cat=catbiom, Tegu=tegubiom)
save(realdat, file="./Data/Reformatted data/Experimental data.Rdata")

###########################################################
#VALIDATION DATA OF CAT/TEGU
###########################################################

#Load bending files
#First combine all cat and tegu files into one file
#add an additional header for "conditions.zyg" and "conditions.cent"
val.dirs<-c("./Data/Validation results/Current/Final")
dir.missing<-c("./Data/Validation results/Missing")
dirnames<-list.files(val.dirs, ".txt")

#load validation data
ROM.val<-load.insilico(dir.bending=val.dirs, dir.missing)

# #Load each validation file
# ROM.val<-list()
# for(i in 2:length(val.dirs)){
#   ROM.vali<-load.insilico(dir.bending=val.dirs[i], dir.missing)
#   ROM.vali<-lapply(ROM.vali, function(x) data.frame(file=rep(dirnames[i-1], nrow(x)), x))
#   ROM.val<-append(ROM.val, ROM.vali)
# }

##############################################################
#Calculate stiffness

#Read morphological data
dir.morph<-c("./Data/Validation results/Morph vert")
morph.all<-load.morphdata(dir.morph, ROM=ROM.val)

#save joint files
for (i in 1:length(morph.all)){
  write.csv(morph.all[[i]], file=paste0("./Data/Validation results/Morph joint/", list.files(dir.morph, pattern=".csv")[i]))
}

#stiffness
func.val<-CalcStiffness(ROM=ROM.val, morphdat=morph.all)

#save data
func_long<-abind::abind(func.val, along=1)
write.csv(func_long, file=paste0("./Data/Reformatted data/Validation data.csv"))
save(func.val, file="./Data/Reformatted data/Validation data.Rdata")

########################################################################
#make neat
############################################################################

#load in silico data - func
load("./Data/Reformatted data/Validation data.Rdata")
#load experimental data - realdat
load("./Data/Reformatted data/Experimental data.Rdata")

#Combine
func.all<-do.call(rbind.data.frame, func.val)
#func.all<-func.all[which(func.all$Trans==levels(func.all$Trans)[1]),]#no translation

#Make neat dataframe
validation<-func.all[c("Trial","Species","Joint.name","Joint", 
                       "js" , "it", "cent.strain","zyg.strain", "Trans","Trans.fact","Zyg.const","Cent.const",
                       "Latero","DV_ROM", "Tor", "LATERO", "DV_Stiff")]
valreal<-lapply(realdat, function(x) x[c("file","Species","JOINT", "jointno",  "Latero2", "Sag", "Tor2", "LATERO_B",  "DV_Stiff")])
valreal<-do.call(rbind.data.frame, valreal)
valreal<-data.frame(valreal, js=rep(NA, nrow(valreal)), it=rep(NA, nrow(valreal)), cent.strain=rep(NA, nrow(valreal)),zyg.strain=rep(NA, nrow(valreal)),  
                    Trans=rep(NA, nrow(valreal)),Trans.fact=rep(NA, nrow(valreal)), Zyg.const=rep(NA, nrow(valreal)), Cent.const=rep(NA, nrow(valreal)))
colnames(validation)<- c("Trial", "Species","JOINT.name","JOINT",
                         "js" , "it", "cent.strain","zyg.strain","Trans","Trans.fact","Zyg.const","Cent.const",
                         "Latero", "Sag", "Tor", "LatStiff", "SagStiff")
colnames(valreal)<-c("Trial", "Species","JOINT.name","JOINT",  "Latero", "Sag", "Tor", "LatStiff", "SagStiff",
                     "js" , "it", "cent.strain","zyg.strain", "Trans","Trans.fact","Zyg.const","Cent.const")
valreal[,c("LatStiff", "SagStiff")]<-1/valreal[,c("LatStiff", "SagStiff")]#convert to stiffness
valreal<-valreal[,match(colnames(valreal), colnames(validation))]
validation<-rbind(validation, valreal)

#Size
#Cat - SEP38, 3.142kg, 0.52m, meanCL = 10.7396
#Tegu - SEP103 - 2.7kg, 0.4m, meanCL = 14.21082
scalefact<-list(Cat=10.7396, Tegu=14.21082)

#Scale by CL
validation[which(validation$Trial!="real"&validation$Species=="Cat"),c("LatStiff", "SagStiff")]<-validation[which(validation$Trial!="real"&validation$Species=="Cat"),c("LatStiff", "SagStiff")]/scalefact$Cat
validation[which(validation$Trial!="real"&validation$Species=="Tegu"),c("LatStiff", "SagStiff")]<-validation[which(validation$Trial!="real"&validation$Species=="Tegu"),c("LatStiff", "SagStiff")]/scalefact$Tegu

###################################################
#final data for suppl.
write.csv(validation, "./Results/Final ROM data.csv")
final.func<-func_long[,c("Trial", "Species" , "Joint.name", "Joint",
                         "js" , "it", "cent.strain","zyg.strain", "Trans","Zyg.const","Cent.const",
                         "Tor", "Latero", "DV_ROM", "TorL_type.", "TorR_type.", "LatL_type.",
                         "LatR_type.", "DFL_type." , "VFLtype..")]
colnames(final.func)<-c("Trial", "Species" , "Joint.name", "Joint", 
                        "Joint Spacing" , "Int. Thresh.", "cent.strain","zyg.strain", "Transl","Zyg.Constraint","Cent.Constraint",
                        "Axial", "Latero", "Sag", "AxL_type", "AxR_type", "LatL_type",
                        "LatR_type", "DFL_type" , "VFLtype")
write.csv(final.func, "./Results/Final constraint data.csv")
save(comp1, file="./Results/Example data.Rdata")

#########################################################################################
###Figures Testing different parameter assumptions
#########################################################################################
#6.7.2021
#Compare real with "best data"

#Set colors
linecols<-c("#425D78", "#E74C3C", "#FFB03B")
confcols<-c("#425D78", "#E74C3C", "#FFB03B")
confcols<-adjustcolor(confcols,alpha.f = 0.5)

#set up different comparisons
#comparison1 - real data vs "best" estimate with errors
comp1<-list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
            Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
            Catest=subset(validation, Trial=="Error"&Species=="Cat"),
            Tegutest=subset(validation, Trial=="Error"&Species=="Tegu"))

#Comparison 2 - Translation
#validation$Trans.fact<-as.factor(validation$Trans.fact)
comp2 <- list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
              Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
              CatNoTrans=subset(validation, Trial=="Translation"&Species=="Cat"&Trans.fact==0),
              TeguNoTrans=subset(validation, Trial=="Translation"&Species=="Tegu"&Trans.fact==0),
              CatTr0.005=subset(validation, Trial=="Translation"&Species=="Cat"&Trans.fact==0.005),
              TeguTr0.005=subset(validation, Trial=="Translation"&Species=="Tegu"&Trans.fact==0.005),
              CatTr0.01=subset(validation, Trial=="Translation"&Species=="Cat"&Trans.fact==0.01),
              TeguTr0.01=subset(validation, Trial=="Translation"&Species=="Tegu"&Trans.fact==0.01),
              CatTr0.05=subset(validation, Trial=="Translation"&Species=="Cat"&Trans.fact==0.05),
              TeguTr0.05=subset(validation, Trial=="Translation"&Species=="Tegu"&Trans.fact==0.05))


#Comparison 3 - Compare COR location
comp3<-list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
            Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
            CatDisc=subset(validation, Trial=="Translation"&Species=="Cat"&Trans.fact=="0"),
            TeguDisc=subset(validation, Trial=="Translation"&Species=="Tegu"&Trans.fact=="0"),
            CatZyg=subset(validation, Trial=="CORzyg"&Species=="Cat"),
            TeguZyg=subset(validation, Trial=="CORzyg"&Species=="Tegu"))

#comparison 4 - centrum displacement
comp4 <-list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
             Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
             CatCentNo=subset(validation, Trial=="CentDisart"&Species=="Cat"&cent.strain==0),
             TeguCentNo=subset(validation, Trial=="CentDisart"&Species=="Tegu"&cent.strain==0),
             CatCent0.5=subset(validation, Trial=="CentDisart"&Species=="Cat"&cent.strain==0.5),
             TeguCent0.5=subset(validation, Trial=="CentDisart"&Species=="Tegu"&cent.strain==0.5),
             CatCent0.25=subset(validation, Trial=="CentDisart"&Species=="Cat"&cent.strain==0.25),
             TeguCent0.25=subset(validation, Trial=="CentDisart"&Species=="Tegu"&cent.strain==0.25)
)


#Comparison 5 - Zygapophyseal displacement
comp5 <- list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
              Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
              CatdisNo=subset(validation, Trial=="ZygDisart"&Species=="Cat"&zyg.strain==0),
              TegudisNo=subset(validation, Trial=="ZygDisart"&Species=="Tegu"&zyg.strain==0),
              Catdis1=subset(validation, Trial=="ZygDisart"&Species=="Cat"&zyg.strain==1),
              Tegudis1=subset(validation, Trial=="ZygDisart"&Species=="Tegu"&zyg.strain==1),
              Catdis75=subset(validation, Trial=="ZygDisart"&Species=="Cat"&zyg.strain==0.75),
              Tegudis75=subset(validation, Trial=="ZygDisart"&Species=="Tegu"&zyg.strain==0.75),
              Catdis50=subset(validation, Trial=="ZygDisart"&Species=="Cat"&zyg.strain==0.5),
              Tegudis50=subset(validation, Trial=="ZygDisart"&Species=="Tegu"&zyg.strain==0.5)
)

#comparison 6 - Intersection threshold
comp6 <- list(Catreal=subset(validation, Trial=="real"&Species=="Cat"),
              Tegureal=subset(validation, Trial=="real"&Species=="Tegu"),
              CatInt0=subset(validation, Trial=="IntThresh"&Species=="Cat"&it==0),
              TeguInt0=subset(validation, Trial=="IntThresh"&Species=="Tegu"&it==0),
              CatInt0.25=subset(validation, Trial=="IntThresh"&Species=="Cat"&it==0.0025),
              TeguInt0.25=subset(validation, Trial=="IntThresh"&Species=="Tegu"&it==0.0025),
              CatInt0.5=subset(validation, Trial=="IntThresh"&Species=="Cat"&it==0.005),
              TeguInt0.5=subset(validation, Trial=="IntThresh"&Species=="Tegu"&it==0.005),
              CatInt1=subset(validation, Trial=="IntThresh"&Species=="Cat"&it==0.01),
              TeguInt1=subset(validation, Trial=="IntThresh"&Species=="Tegu"&it==0.01)
)

#############################################################################################
#plot

allbars<-list()
alldots<-list()

#Comparison1
factname<-rep(c("Experimental", "Autobend"), each=2)
for(i in 1:length(comp1)){
  comp1[[i]]$Trial<-factor(rep(factname[i], nrow(comp1[[i]])), levels=c("Experimental", "Autobend"))
}

pout<-compile.plots(plotdat=comp1, linecols, confcols, ylim_cc=c(0,50), ylim_bar=c(0,35), 
              pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), title=c("A. Experimental Data","","B. AutoBend",""),
              conf.int=rep(TRUE,4), FctOrder=c("Experimental", "Autobend"))


# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]]),
#                widths=c(2,2,1),
#                nrow = 2)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]]),
               widths=c(2,2),
               nrow = 2)

pdf(file=paste0("./Figures/New/Best model.pdf"), width=11, height=8)
grid.draw(g)
dev.off()

allbars<-c(allbars, "Exp dat"=list(print(pout$Bar[[1]]+ggtitle("A.Experimental Data"))))
alldots<-c(alldots, "Experimental data"=list(print(pout$dot+ggtitle("A.Autobend model"))))

#Comparison3
factname<-rep(c("Experimental", "CORdisc", "CORzyg"), each=2)
for(i in 1:length(comp3)){
  comp3[[i]]$Trial<-factor(rep(factname[i], nrow(comp3[[i]])), levels=c("Experimental", "CORdisc", "CORzyg"))
}

pout<-compile.plots(plotdat=comp3, linecols, confcols, ylim_cc=c(0,50), ylim_bar=c(0,35), 
                    pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), conf.int=c(rep(TRUE,2), rep(FALSE,4)),
                    title=c("A. Experimental Data","","B. CORdisc","","C. CORzyg",""), FctOrder=c("Experimental", "CORdisc", "CORzyg"))


# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]],
#                           pout$CCplot[[5]], pout$CCplot[[6]], pout$Bar[[5]]
# ),
# widths=c(2,2,1),
# nrow = 3)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]],
                          pout$CCplot[[5]], pout$CCplot[[6]]
),
widths=c(2,2),
nrow = 3)

pdf(file=paste0("./Figures/New/COR location.pdf"), width=11, height=12)
grid.draw(g)
dev.off()

allbars<-c(allbars, "CORdisc"=list(print(pout$Bar[[3]]+ggtitle("B.CORdisc"))),
           "CORzyg"=list(print(pout$Bar[[5]]+ggtitle("C.CORzyg"))))
alldots<-c(alldots, "COR location"=list(print(pout$dot+ggtitle("B.COR location"))))

#Comparison2
factname<-rep(c("Experimental", "No Translation", "0.5% Translation", "1% Translation", "5% Translation"), each=2)
for(i in 1:length(comp2)){
  comp2[[i]]$Trial<-factor(rep(factname[i], nrow(comp2[[i]])), levels=c("Experimental", "No Translation", "0.5% Translation", "1% Translation", "5% Translation"))
}


pout<-compile.plots(plotdat=comp2, linecols, confcols, ylim_cc=c(0,50), ylim_bar=c(0,35), 
                    pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), conf.int=c(rep(TRUE,2), rep(FALSE,8)),
                    title=c("A. Experimental Data","","B. No Translation","","C. 0.5% Translation","",
                            "D. 1% Translation","","E. 5% Translation",""),
                    FctOrder=c("Experimental", "No Translation", "0.5% Translation", "1% Translation", "5% Translation"))


# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]],
#                           pout$CCplot[[5]], pout$CCplot[[6]], pout$Bar[[5]],
#                           pout$CCplot[[7]], pout$CCplot[[8]], pout$Bar[[7]],
#                           pout$CCplot[[9]], pout$CCplot[[10]], pout$Bar[[9]]),
#                widths=c(2,2,1),
#                nrow = 5)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]],
                          pout$CCplot[[5]], pout$CCplot[[6]],
                          pout$CCplot[[7]], pout$CCplot[[8]],
                          pout$CCplot[[9]], pout$CCplot[[10]]),
               widths=c(2,2),
               nrow = 5)

pdf(file=paste0("./Figures/New/Translation.pdf"), width=11, height=20)
grid.draw(g)
dev.off()

allbars<-c(allbars, "No transl"=list(print(pout$Bar[[3]]+ggtitle("D.No Translation"))),
           "0.5% transl"=list(print(pout$Bar[[5]]+ggtitle("E.0.5% Translation"))),
           "1% transl"=list(print(pout$Bar[[7]]+ggtitle("F.1% Translation"))),
           "5% transl"=list(print(pout$Bar[[9]]+ggtitle("G.5% Translation"))))
alldots<-c(alldots, "Translation"=list(print(pout$dot+ggtitle("C.Translation"))))

#Comparison5
factname<-rep(c("Experimental", "No zyg constraint", "Total disart.", "75% disart.", "50% disart."), each=2)
for(i in 1:length(comp5)){
  comp5[[i]]$Trial<-factor(rep(factname[i], nrow(comp5[[i]])), levels=c("Experimental", "No zyg constraint", "Total disart.", "75% disart.", "50% disart."))
}


pout<-compile.plots(plotdat=comp5, linecols, confcols, ylim_cc=c(0,50), ylim_bar=c(0,35), 
                    pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), conf.int=c(rep(TRUE,2), rep(FALSE,6)),
                    title=c("A. Experimental Data","","B. No zyg constraint","",
                            "C. Total disart.","","D. 75% disart.","","E. 50% disart.",""),
                    FctOrder= c("Experimental", "No zyg constraint", "Total disart.", "75% disart.", "50% disart."))


# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]],
#                           pout$CCplot[[5]], pout$CCplot[[6]], pout$Bar[[5]],
#                           pout$CCplot[[7]], pout$CCplot[[8]], pout$Bar[[7]],
#                           pout$CCplot[[9]], pout$CCplot[[10]], pout$Bar[[9]]),
#                widths=c(2,2,1),
#                nrow = 5)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]],
                          pout$CCplot[[5]], pout$CCplot[[6]],
                          pout$CCplot[[7]], pout$CCplot[[8]],
                          pout$CCplot[[9]], pout$CCplot[[10]]),
               widths=c(2,2),
               nrow = 5)

pdf(file=paste0("./Figures/New/Zygapophyseal disarticulation.pdf"), width=11, height=20)
grid.draw(g)
dev.off()

allbars<-c(allbars, "No zyg"=list(print(pout$Bar[[3]]+ggtitle("H.No Zygapophyseal constraint"))),
           "100% zyg"=list(print(pout$Bar[[5]]+ggtitle("I.100% Zygapophyseal strain"))), 
           "75% zyg"=list(print(pout$Bar[[7]]+ggtitle("J.75% Zygapophyseal strain"))),
           "50% zyg"=list(print(pout$Bar[[9]]+ggtitle("K.50% Zygapophyseal strain"))))
alldots<-c(alldots, "Zyg Disarticulation"=list(print(pout$dot+ggtitle("D.Zyg. Disarticulation"))))

#Comparison4
factname<-rep(c("Experimental", "No cent constraint", "50% strain", "25% strain."), each=2)
for(i in 1:length(comp4)){
  comp4[[i]]$Trial<-factor(rep(factname[i], nrow(comp4[[i]])), levels=c("Experimental", "No cent constraint", "50% strain", "25% strain."))
}

pout<-compile.plots(plotdat=comp4, linecols, confcols, ylim_cc=c(0,65), ylim_bar=c(0,40), 
                    pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), conf.int=c(rep(TRUE,2), rep(FALSE,4)),
                    title=c("A. Experimental Data","","B. No cent constraint","",
                            "C. 50% strain","","D. 25% strain",""),
                    FctOrder= c("Experimental", "No cent constraint", "50% strain", "25% strain."))

# 
# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]],
#                           pout$CCplot[[5]], pout$CCplot[[6]], pout$Bar[[5]],
#                           pout$CCplot[[7]], pout$CCplot[[8]], pout$Bar[[7]]),
#                widths=c(2,2,1),
#                nrow = 4)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]],
                          pout$CCplot[[5]], pout$CCplot[[6]],
                          pout$CCplot[[7]], pout$CCplot[[8]]),
               widths=c(2,2),
               nrow = 4)

pdf(file=paste0("./Figures/New/Centrum strain.pdf"), width=11, height=16)
grid.draw(g)
dev.off()

allbars<-c(allbars, "No cent"=list(print(pout$Bar[[3]]+ggtitle("L.No centrum constraints"))),
           "50% cent"=list(print(pout$Bar[[5]]+ggtitle("M.50% Centrum strain"))),
           "25% cent"=list(print(pout$Bar[[7]]+ggtitle("N.25% Centrum strain"))))

alldots<-c(alldots, "Centrum strain"=list(print(pout$dot+ggtitle("E.Centrum strain"))))

#Comparison6
factname<-rep(c("Experimental", "No Intersection", "0.25% Intersection", "0.5% Intersection","1% Intersection"), each=2)
for(i in 1:length(comp6)){
  comp6[[i]]$Trial<-factor(rep(factname[i], nrow(comp6[[i]])), levels=c("Experimental", "No Intersection", "0.25% Intersection", "0.5% Intersection","1% Intersection"))
}


pout<-compile.plots(plotdat=comp6, linecols, confcols, ylim_cc=c(0,60), ylim_bar=c(0,35), 
                    pagedim=c(3,2), mar=c(7, 6, 4.1, 2.1), conf.int=c(rep(TRUE,2), rep(FALSE,6)),
                    title=c("A. Experimental Data","","B. No Intersection","",
                            "C. 0.25% Intersection","","D. 0.5% Intersection","","E. 1% Intersection",""),
                    FctOrder=c("Experimental", "No Intersection", "0.25% Intersection", "0.5% Intersection","1% Intersection"))


# g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]], pout$Bar[[1]],
#                           pout$CCplot[[3]], pout$CCplot[[4]], pout$Bar[[3]],
#                           pout$CCplot[[5]], pout$CCplot[[6]], pout$Bar[[5]],
#                           pout$CCplot[[7]], pout$CCplot[[8]], pout$Bar[[7]],
#                           pout$CCplot[[9]], pout$CCplot[[10]], pout$Bar[[9]]),
#                widths=c(2,2,1),
#                nrow = 5)
g<-arrangeGrob(grobs=list(pout$CCplot[[1]], pout$CCplot[[2]],
                          pout$CCplot[[3]], pout$CCplot[[4]],
                          pout$CCplot[[5]], pout$CCplot[[6]],
                          pout$CCplot[[7]], pout$CCplot[[8]],
                          pout$CCplot[[9]], pout$CCplot[[10]]),
               widths=c(2,2),
               nrow = 5)

pdf(file=paste0("./Figures/New/Intersection threshold.pdf"), width=11, height=20)
grid.draw(g)
dev.off()

allbars<-c(allbars, "No IT"=list(print(pout$Bar[[3]]+ggtitle("O.0% Intersection"))),
           "0.25% IT"=list(print(pout$Bar[[5]]+ggtitle("P.0.25% Intersection"))),
           "0.5% IT"=list(print(pout$Bar[[7]]+ggtitle("Q.0.5% Intersection"))),
           "1% IT"=list(print(pout$Bar[[9]]+ggtitle("R.1% Intersection"))))

alldots<-c(alldots, "Intersection threshold"=list(print(pout$dot+ggtitle("F.Intersection threshold"))))

#Comparison barchart

lay <- rbind(c(1,2,3,NA),
             c(4,5,6,7),
             c(8,9,10,11),
             c(12,13,14, NA),
             c(15,16,17,18))
g<-arrangeGrob(grobs=allbars, layout_matrix = lay)

pdf(file=paste0("./Figures/New/Barchart summary.pdf"), width=12, height=20)
grid.draw(g)
dev.off()

#Comparison dotplot

g<-arrangeGrob(grobs=alldots, nrow = 3)

pdf(file=paste0("./Figures/New/Dotplot summary.pdf"), width=14, height=14)
grid.draw(g)
dev.off()

############################################################################################
#Stiffness comparison figure
##########################################################################################

#Craniocaudal plots
ylab<-c(rep("Stiffness",2),rep("Est. Stiffness",2))
ylim<-list(c(0.1,0.35),c(0.1,0.35),c(0.4,1),c(0.4,1))
ROMplots<-list()
dev.new(width=3, height=2)
par(mar=c(7, 6, 4.1, 2.1))
for(i in 1:length(comp1)){
  
  plotdata<-comp1[[i]]
  plotdata$JOINT.name<-droplevels(plotdata$JOINT.name)
  plotdata<-tidyr::gather(plotdata,"Direction","Stiff", "LatStiff", "SagStiff")
  plotdata$Direction<-revalue(plotdata$Direction, c('SagStiff'="1Sag", "LatStiff"="2Latero"))
  plotdata$Direction<-as.factor(plotdata$Direction)
  craniocaudPlot(x=plotdata$JOINT.name,y=plotdata$Stiff,group=plotdata$Direction, 
                 linecols = linecols, confcols = confcols, legend=F, xlab="Joints",
                 ylab=ylab[i], ylim=ylim[[i]])
  grid.echo()
  ROMplots[[i]]<-grid.grab()
  
}

#Barplots
ylab<-c(rep("Stiffness",1),rep("Stiffness",1),rep("Est. Stiffness",1),rep("Est. Stiffness",1))
ylim<-list(c(0,0.3),c(0,0.3),c(0,1),c(0,1))
barplots<-list()
#ROM
for(i in seq(1,4, by=2)){
  Joints<-"JOINT"
  Vars<-c("SagStiff", "LatStiff")
  col.scale<-c("SagStiff"=linecols[1], "LatStiff"=linecols[2])
  
  valbar<-do.call(rbind.data.frame, list(comp1[[i]], comp1[[i+1]]))
  valbar<-split(valbar, valbar$Species)
  barplots[[i]]<-print(ROM.barplot(ROM=valbar, Joints, Vars,order=c("Cat", "Tegu"), col.scale, 
                             ylab=ylab[i], ylim=ylim[[i]])$plot)
}



#Arrange
g<-arrangeGrob(grobs=list(ROMplots[[1]], ROMplots[[2]], barplots[[1]],
                          ROMplots[[3]], ROMplots[[4]], barplots[[3]]),
               widths=c(2,2,1),
               nrow = 2)


pdf(file=paste0("./Figures/New/validation_Stiffness.pdf"), width=12, height=8)
grid.draw(g)
dev.off()


#######################################################################################
#Impact of different parameters
######################################################################################

res.plots<-list()
for(i in 1:2){
test<-func.val[[i]]
test<-subset(test, Trial=="Error")
colnames(test)[14]<-"Ax"

#test<-test[which(test$Trans==levels(test$Trans)[1]),]#no translation
test_long<-pivot_longer(test,cols=c(Latero , DV_ROM ,Ax), names_to="Direction", values_to = "Mobility")

#remove effect of joint
test_long$js<-as.factor(test_long$js)
test_long$it<-as.factor(test_long$it)
test_long$cent.strain<-as.factor(test_long$cent.strain)
lm.param<-lm(Mobility~Joint.name*Direction+js+it+cent.strain, data=test_long)
aov.param<-anova(lm.param)
write.csv(aov.param, file=paste0("./Figures/",names(func.val)[i],"_ANOVA.csv"))
lm.param.res<-lm(Mobility~Joint.name*Direction, data=test_long)
test_long<-data.frame(test_long, residuals=lm.param.res$residuals)

#revalue
test_long<-pivot_longer(test_long,cols=c(js ,it ,cent.strain), names_to="Parameter", values_to = "Value")
test_long$Value<-revalue(as.factor(test_long$Value),
                         c('0.1'="High", '-0.1'="Low",'0.0025'="Low", '0.005'="High", '0.45'="Low", '0.55'="High"))
test_long$Direction<-revalue(test_long$Direction, c(Latero="Lat", Ax="Ax", DV_ROM="Sag"))
test_long$Parameter<-revalue(test_long$Parameter, c(it="Intersection thresh.", js="Joint Spacing",
                                                    cent.strain="Joint strain"))

#Piechart
#pies<-c(aov.param$`Sum Sq`[c(1,2,6)], sum(aov.param$`Sum Sq`[3:5]), aov.param$`Sum Sq`[c(7)])
pdf(file=paste0("./Figures/New/ANOVApie_validation_", names(func.val)[i],".pdf"), width=6, height=6)
#plotrix::pie3D(aov.param$`Sum Sq`[c(1,2,6,3,4,5,7)],  col=c(brewer.pal(7,"Set1")), 
#                start=pi/2, explode = 0.1)
pie(aov.param$`Sum Sq`[c(1,2,6,3,4,5,7)],  col=c(brewer.pal(7,"Set1")), labels= rownames(aov.param)[c(1,2,6,3,4,5,7)])
dev.off()

#Boxplot
test_long<-group_by(data.frame(test_long), Parameter)
p <- ggplot(test_long, aes(x=Direction, y=residuals)) +
      geom_boxplot(aes(fill=Value))+
   labs(title=names(func.val)[i]) +
  scale_fill_manual(values=c("white", "grey")) +
  facet_grid(cols=vars(Parameter))
res.plots[[i]]<-p

}

g<-arrangeGrob(grobs=res.plots,
               nrow = 2)

pdf(file=paste0("./Figures/New/Parameter residuals barplot.pdf"), width=10, height=10)
grid.draw(g)
dev.off()

#####################################################################################
#Barchart of bony stops
######################################################################################
#Get bony stops data and reformat

#inter.val<-lapply(func.val[c(1:2)], function(x) x[which(x$Trans==levels(x$Trans)[1]),])#no translation
inter.val<-lapply(func.val, function(x) subset(x, Trial=="Error"))
inter_byspecies<-ConstonMotion(Funcdat=inter.val, splitby=c("Species"), rows=c("JOINT"))

#save results with most common intersections
func_long<-do.call(rbind.data.frame, inter.val)
func_long<-aggregate(func_long, by=list(spj=interaction(func_long$Species, func_long$Joint)), FUN=mean, na.rm=T)
out<-data.frame(func_long[colSums(!is.na(func_long)) > 0],abind::abind(inter_byspecies$comm_stops, along=1))
write.csv(out, file=paste0("./Results/Validation results_by joint.csv"))

#Arrange
g<-arrangeGrob(grobs=list(inter_byspecies$plots$Cat, inter_byspecies$plots$Tegu),
               nrow = 2)


pdf(file=paste0("./Figures/New/Intersections.pdf"), width=10, height=12)
grid.draw(g)
dev.off()


