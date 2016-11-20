args = commandArgs(trailingOnly=TRUE)
baseLocation = args[1]
fileLocation = args[2]
name = args[3]

library(amap, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(rgl, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(e1071, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(scales, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(bios2mds, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))

mytable = read.table(file=paste(fileLocation, "/Matrix.dist", sep = ""), sep=',')
headers = read.table(file=paste(fileLocation, "/Header.head", sep = ""), sep=',')
sampleDetails = read.table(file=paste(fileLocation, "/SampleInfo.info", sep=""), sep=',')
colnames(sampleDetails) = c('Sample', 'Taxon', 'Reef', 'Region', 'MajITS', 'Type')
colnames(mytable) = headers$V1
rownames(mytable) = headers$V1
myMatrix = as.matrix(mytable)
myMDS = mmds(myMatrix, pc=2)
twoDMDSPlot = myMDS$coord
plotDetails = cbind(twoDMDSPlot, sampleDetails)


svg(paste(name, ".svg", sep=""))
plot(plotDetails$PC1, plotDetails$PC2, col=c('darkblue','green3','red','deeppink', 'purple','darkorchid','darkgreen','darkorange3', 'dodgerblue3', 'darkred', 'darkolivegreen',  'cornflowerblue', 'chartreuse1', 'cyan1', 'brown4')[((as.numeric(plotDetails$MajITS))%%14)+1])
dev.off()



print("You made it to the end of the script with 0 problems! Good work dude!")
