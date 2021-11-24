### Script 1: Read data ###

setwd("/Users/katarinameramo/Desktop/PHD/Caatinga/data")
library(readxl)

da = read_excel("ASI_results.xlsx")
n = dim(da)[1]

plots = rep(NA,n)
dates = rep(NA,n)
for(i in 1:n){
  filename = da$Files[i]
  pos = gregexpr(pattern ='-',filename)[[1]]
  plots[i] = gsub(" ", "",substr(filename, 1, pos-1))
  pos = gregexpr(pattern ='_',filename)[[1]]
  dates[i] = substr(filename, pos[1]+1, pos[2]-1)
}
plots = as.factor(plots)
dates = as.factor(dates)

index = data.frame(file = da$Files, plot = plots, date = dates)

da$Files = NULL

da2 = matrix(0,nrow = dim(da)[1], ncol = dim(da)[2])
for(i in 1:n){
  da2[i,] = as.numeric(da[i,])
}
colnames(da2) = colnames(da)
da = da2

save(da,index,file = "processed_data/data.R")


### Script 2: Make matrix ###

setwd("/Users/katarinameramo/Desktop/PHD/Caatinga/data")

load("processed_data/data.R")
#da,index

n = dim(index)[1]
m = dim(da)[2]
plot.date = rep(NA,n)
for(i in 1:n){
  plot.date[i] = paste0(index$plot[i],"-",index$date[i])
}
index$plot.date = as.factor(plot.date)

plot.dates = levels(index$plot.date)
u = length(plot.dates)

abu.matrix = matrix(NA,nrow = u, ncol = m)
for(i in 1:u){
  abu.matrix[i,] = colSums(da[which(index$plot.date==plot.dates[i]),]>0.9)
}
rownames(abu.matrix)=plot.dates
colnames(abu.matrix)=colnames(da)

save(abu.matrix,file = "processed_data/abu.matrix.R")

#Combine with covariate data

setwd("/Users/katarinameramo/Desktop/PHD/Caatinga/data")
library(readxl)

da = read_excel("Vegetation_R_220421.xlsx")

load("processed_data/abu.matrix.R")
Plot = rownames(abu.matrix)
Plot = gsub("Plot", "P", Plot)  # changing cell names inside column Plot (PlotX -> PX)
Plot = gsub("P2-", "P02-", Plot)  # changing cell names inside column Plot (PlotX -> PX)
Plot = gsub("P8-", "P08-", Plot)  # changing cell names inside column Plot (PlotX -> PX)
rownames(abu.matrix) = Plot
abu.matrix = abu.matrix[-(1:5),] # remove extras
Plot = rownames(abu.matrix)
n = dim(abu.matrix)[1]
season = rep("autumn",n)
ii = rep(NA,n)
for(i in 1:n){
  if(substr(Plot[[i]],7,10)=="1218"){
    season[i]="winter"
  }
  ii[i] = match(substr(Plot[[i]],1,3),da$Plot)
}
all.data = cbind(abu.matrix,da[ii,],season)
all.data$Plot = as.factor(all.data$Plot)
all.data$season = as.factor(all.data$season)
all.data$Openness = as.factor(all.data$Openness)
save(all.data,file = "processed_data/all.data.R")


### Script 3: Fit models ###

setwd("C:/LocalData/meramo/HMSC")
library(Hmsc)

load("processed_data/all.data.R")
n=dim(all.data)[1] 
studyDesign = data.frame(plot = all.data$Plot, sample = as.factor(row.names(all.data)))
rL.plot = HmscRandomLevel(units = levels(studyDesign$plot))
rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
Y = as.matrix(all.data[,1:18])
prev = colSums(Y>0)
sel.sp = prev>4
Y = Y[,sel.sp]
XData = data.frame(GMDI = all.data$GMDI, season = all.data$season)
XFormula = ~ GMDI + season
#m.poisson = Hmsc(Y=Y, XData = XData, XFormula=XFormula,
#distr="lognormal poisson", studyDesign=studyDesign,
#ranLevels=list(plot=rL.plot,sample=rL.sample))

Y.pa = 1*(Y>0)

m.pa = Hmsc(Y=Y.pa, XData = XData, XFormula=XFormula,
            distr="probit", studyDesign=studyDesign,
            ranLevels=list(plot=rL.plot,sample=rL.sample))

Y.abu = Y
Y.abu[Y==0] = NA
Y.abu = scale(log(Y.abu))

m.abu = Hmsc(Y=Y.abu, XData = XData, XFormula=XFormula,
             distr="normal", studyDesign=studyDesign,
             ranLevels=list(plot=rL.plot,sample=rL.sample))

models = list(m.pa, m.abu)
modelnames = c("pa","abu")
save(models,modelnames,file = "models/unfitted_models")



load(file = "models/unfitted_models") #models, modelnames

samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nChains = 4
for(Lst in 1:length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  nm = length(models)
  for (model in 1:nm) {
    print(paste0("model = ",modelnames[model]))
    m = models[[model]]
    m = sampleMcmc(m, samples = samples, thin=thin,
                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr), 
                   transient = ceiling(0.5*samples*thin),
                   nChains = nChains) 
    models[[model]] = m
  }
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),
                   ".Rdata",sep = "")
  save(models,modelnames,file=filename)
}


### Script 4: Evaluate convergence ###

setwd("C:/LocalData/meramo/HMSC") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)
library(colorspace)
library(vioplot)

#include in samples_list and thin_list only those models that you have actually fitted!
samples_list = c(5,250,250,250,250,250)
thin_list = c(1,1,10,100,1000,10000)
nst = length(thin_list)
nChains = 4

ma = NULL
na = NULL
for (Lst in 1:nst) {
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = paste("models/models_thin_", as.character(thin),
                   "_samples_", as.character(samples),
                   "_chains_",as.character(nChains),".Rdata",sep = "")
  load(filename)
  nm = length(models)
  for(j in 1:nm){
    mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
    tmp = summary(psrf.beta)
    if(is.null(ma)){
      ma=psrf.beta[,1]
      na = paste0(as.character(thin),",",as.character(samples))
    } else {
      ma = cbind(ma,psrf.beta[,1])
      if(j==1){
        na = c(na,paste0(as.character(thin),",",as.character(samples)))
      } else {
        na = c(na,"")
      }
    }
  }
}

pdf(file=paste("panels/MCMC_convergence.pdf"))
par(mfrow=c(2,1))
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0,max(ma)),main="psrf(beta)")
vioplot(ma,col=rainbow_hcl(nm),names=na,ylim=c(0.9,1.1),main="psrf(beta)")
dev.off()


### Script 5: Compute model fit ###

setwd("C:/LocalData/meramo/HMSC") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)

#You may wish to loop over samples_list and thinning_list as done in Script S3 
nChains = 4
thin = 1000 #vaihdoin 1:stä tonniin
samples = 250 #vaihdoin 5:stä 250 #ajossa kesti 1,5 vrk
print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
filename_in = paste("models/models_thin_", as.character(thin),
                    "_samples_", as.character(samples),
                    "_chains_",as.character(nChains),
                    ".Rdata",sep = "")
load(file = filename_in) #models, modelnames
nm = length(models)

MF = list()
MFCV = list()
WAIC = list()

for(model in 1:nm){
  print(paste0("model = ",as.character(model)))
  m = models[[model]]
  preds = computePredictedValues(m)
  MF[[model]] = evaluateModelFit(hM=m, predY=preds)
  partition = createPartition(m, nfolds = 2)
  preds = computePredictedValues(m,partition=partition) #nParallel = nChains
  MFCV[[model]] = evaluateModelFit(hM=m, predY=preds)
  WAIC[[model]] = computeWAIC(m)
}

filename_out = paste("models/MF_thin_", as.character(thin),
                     "_samples_", as.character(samples),
                     "_chains_",as.character(nChains),
                     ".Rdata",sep = "")
save(MF,MFCV,WAIC,modelnames,file = filename_out)


### Script 6: Show model fit ###

setwd("C:/LocalData/meramo/HMSC") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)

#This script shows model fit for probit and linear models.
#For Poisson models, you may plot e.g. pseudoR2, see e.g. the book for examples

thin = 1000 #vaihdoin yhdestä tonniin taas
samples = 250 #vaihdoin vitosesta 250
nChains = 4

filename = paste("models/MF_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(MF)
filename = paste("panels/model_fit.pdf")
pdf(file = filename)
for(j in 1:nm){
  cMF = MF[[j]]
  cMFCV = MFCV[[j]]
  if(!is.null(cMF$TjurR2)){
    plot(cMF$TjurR2,cMFCV$TjurR2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": Tjur R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$R2)){
    plot(cMF$R2,cMFCV$R2,xlim=c(-1,1),ylim=c(-1,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": R2"))
    abline(0,1)
    abline(v=0)
    abline(h=0)
  }
  if(!is.null(cMF$AUC)){
    plot(cMF$AUC,cMFCV$AUC,xlim=c(0,1),ylim=c(0,1),
         xlab = "explanatory power",
         ylab = "predictive power",
         main=paste0(modelnames[[j]],", thin = ",as.character(thin),", samples = ",as.character(samples),": AUC"))
    abline(0,1)
    abline(v=0.5)
    abline(h=0.5)
  }
}
dev.off()


### Script 7: Show parameter estimates ###

setwd("C:/LocalData/meramo/HMSC") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)
library(Hmsc)
library(colorspace)
library(corrplot)
library(writexl)

nChains = 4 
samples = 250 
thin = 1000 

filename = paste("models/models_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(models)

filename = paste("panels/parameter_estimates.pdf")

pdf(file = filename)
for(j in 1:nm){
  m = models[[j]]
  
  VP = computeVariancePartitioning(m)
  vals = VP$vals
  mycols =  cm.colors(nrow(VP$vals)) #rainbow.. tsekkaa colorspace paketti
  plotVariancePartitioning(hM=m, VP=VP,cols = mycols, args.leg=list(bg="white",cex=0.7),
                           main = paste0("Proportion of explained variance, ",modelnames[[j]]),cex.main=0.8)
  preds = computePredictedValues(m)
  MF = evaluateModelFit(hM=m, predY=preds)
  
  R2 = NULL
  if(!is.null(MF$TjurR2)){
    TjurR2 = MF$TjurR2
    vals = rbind(vals,TjurR2)
    R2=TjurR2
  }
  if(!is.null(MF$R2)){
    R2=MF$R2
    vals = rbind(vals,R2)
  }
  
  filename =  paste0("panels/parameter_estimates_VP_",modelnames[[j]],".csv")
  write.csv(vals,file=filename)
  
  if(!is.null(R2)){
    VPr = VP
    for(k in 1:m$ns){
      VPr$vals[,k] = R2[k]*VPr$vals[,k]
    }
    
    VPr$vals = VPr$vals[,order(-R2)]
    plotVariancePartitioning(hM=m, VP=VPr,cols = mycols, args.leg=list(bg="white",cex=0.7),ylim=c(0,1),
                             main=paste0("Proportion of raw variance, ",modelnames[[j]]),cex.main=0.8)
  }
  
}

for(j in 1:nm){
  m = models[[j]]
  postBeta = getPostEstimate(m, parName="Beta")
  show.sp.names = (is.null(m$phyloTree) && m$ns<=20) 
  plotBeta(m, post=postBeta, supportLevel = 0.75,param="Sign", #muutin 0.95 0.75:ksi
           plotTree = !is.null(m$phyloTree),
           covNamesNumbers = c(TRUE,FALSE),
           spNamesNumbers=c(show.sp.names,FALSE),
           cex=c(0.6,0.6,0.8))
  mymain = paste0("BetaPlot, ",modelnames[[j]])
  if(!is.null(m$phyloTree)){
    mpost = convertToCodaObject(m)
    rhovals = unlist(poolMcmcChains(mpost$Rho))
    mymain = paste0(mymain,", E[rho] = ",round(mean(rhovals),2),", Pr[rho>0] = ",round(mean(rhovals>0),2))
  }
  title(main=mymain, line=2.5, cex.main=0.8)
  
  me = as.data.frame(t(postBeta$mean))
  me = cbind(m$spNames,me)
  colnames(me) = c("Species",m$covNames)
  po = as.data.frame(t(postBeta$support))
  po = cbind(m$spNames,po)
  colnames(po) = c("Species",m$covNames)
  ne = as.data.frame(t(postBeta$supportNeg))
  ne = cbind(m$spNames,ne)
  colnames(ne) = c("Species",m$covNames)
  vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
  filename = paste0("panels/parameter_estimates_Beta_",modelnames[j],".xlsx")
  writexl::write_xlsx(vals,path = filename)
}

for(j in 1:nm){
  if(m$nt>1){
    m = models[[j]]
    postGamma = getPostEstimate(m, parName="Gamma")
    plotGamma(m, post=postGamma, supportLevel = 0.9, param="Sign",
              covNamesNumbers = c(TRUE,FALSE),
              trNamesNumbers=c(m$nt<21,FALSE),
              cex=c(0.6,0.6,0.8))
    title(main=paste0("GammaPlot ",modelnames[[j]]), line=2.5,cex.main=0.8)
  }
}

for(j in 1:nm){
  m = models[[j]]
  OmegaCor = computeAssociations(m)
  supportLevel = 0.95
  for (r in 1:m$nr){
    plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
    toPlot = ((OmegaCor[[r]]$support>supportLevel) + (OmegaCor[[r]]$support<(1-supportLevel))>0)*sign(OmegaCor[[r]]$mean)
    if(m$ns>20){
      colnames(toPlot)=rep("",m$ns)
      rownames(toPlot)=rep("",m$ns)
    }
    mymain = paste0("Associations, ",modelnames[[j]], ": ",names(m$ranLevels)[[r]])
    if(m$ranLevels[[r]]$sDim>0){
      mpost = convertToCodaObject(m)
      alphavals = unlist(poolMcmcChains(mpost$Alpha[[1]][,1]))
      mymain = paste0(mymain,", E[alpha1] = ",round(mean(alphavals),2),", Pr[alpha1>0] = ",round(mean(alphavals>0),2))
    }
    corrplot(toPlot[plotOrder,plotOrder], method = "color",
             col=colorRampPalette(c("blue","white","red"))(3),
             mar=c(0,0,1,0),
             main=mymain,cex.main=0.8)
    
    me = as.data.frame(OmegaCor[[r]]$mean)
    me = cbind(m$spNames,me)
    colnames(me)[1] = ""
    po = as.data.frame(OmegaCor[[r]]$support)
    po = cbind(m$spNames,po)
    colnames(po)[1] = ""
    ne = as.data.frame(1-OmegaCor[[r]]$support)
    ne = cbind(m$spNames,ne)
    colnames(ne)[1] = ""
    vals = list("Posterior mean"=me,"Pr(x>0)"=po,"Pr(x<0)"=ne)
    filename = paste0("panels/parameter_estimates_Omega_",modelnames[[j]],"_",names(m$ranLevels)[[r]],".xlsx")
    writexl::write_xlsx(vals,path = filename)
  }
}
dev.off()


### Script 8: Make predictions ###

setwd("C:/LocalData/meramo/HMSC") # set directory to the folder where the folders "data", "models" and "panels" are
library(Hmsc)
library(ggplot2)

localDir = "."

nChains = 4
samples = 250 #
thin = 1000 #
filename = paste("models/models_thin_", as.character(thin),
                 "_samples_", as.character(samples),
                 "_chains_",as.character(nChains),
                 ".Rdata",sep = "")
load(filename)
nm = length(models)

filename =  paste("panels/predictions.pdf")
pdf(file = filename)
for(j in 1:nm){
  m = models[[j]]
  covariates = all.vars(m$XFormula)
  ex.sp = which.max(colMeans(m$Y,na.rm = TRUE)) #most common species as example species
  if(m$distr[1,1]==2){
    ex.sp = which.min(abs(colMeans(m$Y,na.rm = TRUE)-0.5)) #for probit models the species with prevalence closest to 0.5
  }
  for(k in 1:(length(covariates))){
    covariate = covariates[[k]]
    Gradient = constructGradient(m,focalVariable = covariate)
    Gradient2 = constructGradient(m,focalVariable = covariate,non.focalVariables = 1)
    predY = predict(m, Gradient=Gradient, expected = TRUE)  
    predY2 = predict(m, Gradient=Gradient2, expected = TRUE)  
    par(mfrow=c(2,1))
    pl = plotGradient(m, Gradient, pred=predY, yshow = 0, measure="S", showData = TRUE, 
                      main = paste0(modelnames[[j]],": summed response (total effect)"))
    if(inherits(pl, "ggplot")){
      print(pl + labs(title=paste0(modelnames[[j]],": summed response (total effect)")))
    }
    pl = plotGradient(m, Gradient2, pred=predY2, yshow = 0, measure="S", showData = TRUE, 
                      main = paste0(modelnames[[j]],": summed response (marginal effect)"))
    if(inherits(pl, "ggplot")){
      print(pl + labs(title=paste0(modelnames[[j]],": summed response (marginal effect)")))
    }
    par(mfrow=c(2,1))
    pl = plotGradient(m, Gradient, pred=predY, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp, showData = TRUE, 
                      main = paste0(modelnames[[j]],": example species (total effect)"))
    if(inherits(pl, "ggplot")){
      print(pl + labs(title=paste0(modelnames[[j]],": example species (total effect)")))
    }
    pl = plotGradient(m, Gradient2, pred=predY2, yshow = if(m$distr[1,1]==2){c(-0.1,1.1)}else{0}, measure="Y",index=ex.sp, showData = TRUE, 
                      main = paste0(modelnames[[j]],": example species (marginal effect)"))
    if(inherits(pl, "ggplot")){
      print(pl + labs(title=paste0(modelnames[[j]],": example species (marginal effect)")))
    }
    if(m$nt>1){
      for(l in 2:m$nt){
        par(mfrow=c(2,1))
        pl = plotGradient(m, Gradient, pred=predY, measure="T",index=l, showData = TRUE,yshow = 0,
                          main = paste0(modelnames[[j]],": community weighted mean trait (total effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[[j]],": community weighted mean trait (total effect)")))
        }
        pl = plotGradient(m, Gradient2, pred=predY2, measure="T",index=l, showData = TRUE, yshow = 0,
                          main = paste0(modelnames[[j]],": community weighted mean trait (marginal effect)"))
        if(inherits(pl, "ggplot")){
          print(pl + labs(title=paste0(modelnames[[j]],": community weighted mean trait (marginal effect)")))
        }
      }
    }
  }
}
dev.off()


#Making of stacked barplots

# library
library(ggplot2)

#working dir
setwd("C:/LocalData/meramo/HMSC/HMSC_thin_1000")

#Bring data, data is the FigB and datab is the figA in the ms
data <- read.csv("parameter_estimates_VP_ abu .csv",header = TRUE, sep = "," )
datab <- read.csv("parameter_estimates_VP_ pa .csv",header = TRUE, sep = "," )

data$X[2] <- "Survey time"
datab$X[2] <- "Survey time"
#data$X[data$X == "season"] <- "Survey time"

#Remove last row
data <- data[-nrow(data),]
datab <- datab[-nrow(datab),]



#Reformat data for plotting
#install.packages("reshape2")
library(reshape2)

data2 <- melt(data, id.vars = c("X"))
datab2 <- melt(datab, id.vars = c("X"))

#Plot, colors can be changed by changing the brewer palette
#to remove legend, change guides part to guides(fill = FALSE)
plot<- ggplot(data2, aes(fill=factor(X, levels=c("Random: sample","Random: plot", "Survey time", "GMDI" )), y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Species") +
  ylab("Variance proportion") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_brewer(palette = "PuBu")+
  scale_x_discrete(labels = abbreviate)
plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             axis.text.x = element_text(size=8, angle=45))

ggsave("PlotB.png")


plotb<- ggplot(datab2, aes(fill=factor(X, levels=c("Random: sample","Random: plot", "Survey time", "GMDI" )), y=value, x=variable)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("Species") +
  ylab("Variance proportion") +
  guides(fill=guide_legend(title=NULL))+
  scale_fill_brewer(palette = "PuBu")+
  scale_x_discrete(labels = abbreviate)
plotb + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text.x = element_text(size=8, angle=45))

ggsave("PlotA.png")

