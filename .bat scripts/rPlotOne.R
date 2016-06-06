args = commandArgs(trailingOnly=TRUE)

# Make a dataframe that has the observations as the information for the plots to be made
# so that the variables are baseLocation, fileLocation, name and isMaj
setOfDataToBePlotted = args[1]

for (argSets in 1:nrow(setOfDataToBePlotted)

baseLocation = args[1]
fileLocation = args[2]
name = args[3]
# plotType will be wither MAJ or CLADE
# If CLADE then we will make the Pareto chart and the Fst plot
# In doing this we will write out the colourVectorDictionary
# This will allow us to have the correct Maj colour when we come to do the MAJ plots
# If MAJ we will read in the colour 
isMaj = args[4]
setwd(fileLocation)

library(amap, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(rgl, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(e1071, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(scales, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(bios2mds, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))
library(qcc, lib.loc=paste(baseLocation, "/Rlibs", sep = ""))

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



if(isMaj == FALSE){# If this is a cladal plotting then generate this and write out the colourVectorDictionary
	print('IsMaj is FALSE')
	##### CREATION OF THE PARETO CHART
	# Get a table that is frequencies of Majs
	w = table(plotDetails$MajITS)
	# Turn this into dataframe so that it can be sorted
	t = as.data.frame(w)
	# Give the first column a name
	names(t)[1] = 'Maj'
	# Sort the table according to frequency
	t = t[with(t, order(-Freq)), ]
	#i <- sapply(newt, is.factor)
	#newt[i] <- lapply(newt[i], as.character)
	# Change the $Maj to be of type character rather than a factor
	t = transform(t, Maj = as.character(Maj))
	# Create a new data.frame that will become the dataframe we use for the Pareto plot
	newt = data.frame(Maj = c('Other'), Freq=c(0), Colour = c('black'), stringsAsFactors = FALSE)
	# Set up the counters that we will use to populate our new dataframe (newt)
	# The totoal number of samples found with Majs that are not the top 10 most abundant
	otherTotal = 0
	# The number of Majs not in the top 10
	otherCount = 0
	# Counter so we know when we have done the top 10 Majs
	counter = 1
	# Colours that we will asign to the top 10
	colours = c('darkblue', 'green3','red','deeppink', 'cornflowerblue', 'darkorchid','darkgreen','darkorange3', 'cyan1', 'dodgerblue3', 'darkred', 'darkolivegreen',  'purple', 'chartreuse1',  'brown4')
	# Loop that that cycles through all of the Majs of t
	# if they are the top 10 then we asign a coulour and add them to the new dataframe newt
	# if they are not top 10 then we add the number of samples they contained to the otherTotal
	# and add 1 to the other count
	for (index in 1:nrow(t)){
	  if (counter <= 10){#Only keep the first 10 values
		newt <- rbind(newt, c(as.character(t[index,][1]), as.numeric(t[index,][2]), colours[index]))
	  }else{
		otherTotal = otherTotal + as.numeric(t[index,][2])
		otherCount = otherCount + 1
	  }
	  counter = counter + 1
	}

	#If there is an others categorie then add it to then
	if(otherCount > 0){newt = rbind(newt, c('Other', otherTotal, 'black'))}
	# Delete the top row
	newt = newt[-1,]
	# Create the labels for the x axis
	lbls = ifelse(newt$Maj == 'Other', paste(newt$Maj, " (", newt$Freq, ")", sep=""), paste(newt$Maj, " (", newt$Freq, ")", sep="")) 
	# Add the lbls to the newt dataframe
	newt = data.frame(Maj = newt$Maj, Colour = newt$Colour, Freq = newt$Freq, Lbls = lbls, stringsAsFactors = FALSE)
	# Makes ure that the Freq variable is numeric
	newt = transform(newt, Freq = as.numeric(Freq))
	# Sort newt again so that the 'other' observation is in the right order and so so are the colours
	newt = newt[order(-newt[,3]),]


	# Make a vector (Freq) with names (lbls) from the newt dataframe to put into the pareto plot
	xData <- newt$Freq
	names(xData) <- newt$Lbls
	# Finally, create the pareto chart within svg() so that it is written to file to the specified size
	svg(paste(name, "_pareto.svg", sep=""), width = 5, height = 5)
	pareto.chart(xData, main='Samples with majority ITS2 seq', col = newt$Colour)
	dev.off()

	# Create and write out the colourVectorDictionary for use in making the MAJ plots
	colourVectorDictionary <- as.vector(newt$Maj)
	names(colourVectorDictionary) <- as.vector(newt$Colour)
	write.table(colourVectorDictionary, file = paste(getwd(), '/colourVectorDictionary.txt', sep=""), sep="\t")

	## PLOT THE FST PLOT
	# We need to have a way of knowing which colour the points should be
	# We will relate a colour according to the MAJ type to a new vector list of colours
	# Cycle through the plotDetials and at each row look search for the Maj in the as.vector(newt$Maj)
	# If true then assign this colour to the colourVector, else assign black
	# Then when making the plot assign colour to the colourVector
	counter = 1
	colourVector <- vector(mode='character', length=nrow(plotDetails))
	for (index in 1:nrow(plotDetails)){
	 if(as.character(plotDetails[index,]$MajITS) %in% as.vector(newt$Maj)){
	   colourVector[counter] <- as.vector(newt$Colour)[match(as.character(plotDetails[index,]$MajITS), as.vector(newt$Maj))]
	 }else{
	   colourVector[counter] <- 'black'
	 }
	  counter = counter + 1
	}
	# Check x axis
	# Only need to check either max or min as they should be symettrical around 0
	# If less thatn 1 then change xAxisChange to TRUE and enforce xlim=c(-1,1)
	if (max(plotDetails$PC1) < 1){xlimValue = c(-1,1)}else{xlimValue = c(min(plotDetails$PC1), max(plotDetails$PC1))}
	# Check y axis
	if (max(plotDetails$PC2) < 1){ylimValue = c(-1,1)}else{ylimValue = c(min(plotDetails$PC2), max(plotDetails$PC2))}

	svg(paste(name, "_FstPlot.svg", sep=""))
	#plot(plotDetails$PC1, plotDetails$PC2, col=colourVector)
	plot(plotDetails$PC1, plotDetails$PC2, col=colourVector, xlim=xlimValue, ylim=ylimValue)
	dev.off()

}else{#If this is a Maj plot
	print('IsMaj is True')
	#Read in the colourVectorDictionary and get it back to a useable vector with names
	# Set the working directory so that we are back in the cladal level directory where the colourVectorDictionary is saved
	setwd('../../')
	readInTable = read.table(file = paste(getwd(), '/colourVectorDictionary.txt', sep=""), sep="\t")
	colourVectorDictionary <- as.vector(readInTable$x)
	names(colourVectorDictionary) <- rownames(readInTable)
	# colourVectorDictionary
	# darkblue         green3            red       deeppink cornflowerblue          black 
    #      "C3"     "Otu15163"           "C1"     "Otu15723"          "C15"        "Other" 
    #darkorchid      darkgreen    darkorange3          cyan1    dodgerblue3 
    #"Otu18013"     "Otu17566"     "Otu23917"     "Otu24379"     "Otu29717" 
	
	# Change the working directory back to fileLocation so that the plot is put in the correct file
	setwd(fileLocation)

	
	# Use colourVectorDictionary as a dictionary and make a new vector colourVector that contains the colours required for each sample in the correct order
	counter = 1
	colourVector <- vector(mode='character', length=nrow(plotDetails))
	for (index in 1:nrow(plotDetails)){
	 if(as.character(plotDetails[index,]$MajITS) %in% as.vector(colourVectorDictionary)){
	   colourVector[counter] <- names(colourVectorDictionary)[match(as.character(plotDetails[index,]$MajITS), as.vector(colourVectorDictionary))]
	 }else{
	   colourVector[counter] <- 'black'
	 }
	  counter = counter + 1
	}
	
	# Check x axis
	# Only need to check either max or min as they should be symettrical around 0
	# If less thatn 1 then change xAxisChange to TRUE and enforce xlim=c(-1,1)
	if (max(plotDetails$PC1) < 1){
		xlimValue = c(-1,1)
		print(paste(plotDetails[1,]$MajITS, xlimValue, sep=""))
	}else{
		xlimValue = c(min(plotDetails$PC1), max(plotDetails$PC1))
		print(paste(plotDetails[1,]$MajITS, xlimValue, sep=""))
	}
	# Check y axis
	if (max(plotDetails$PC2) < 1){
		ylimValue = c(-1,1)
	}else{
		ylimValue = c(min(plotDetails$PC2), max(plotDetails$PC2))
	}
	#Plot the FstPlot using the identified colour
	svg(paste(name, "_FstPlot.svg", sep=""))
	#plot(plotDetails$PC1, plotDetails$PC2, col=colourVector)
	plot(plotDetails$PC1, plotDetails$PC2, col=colourVector, xlim=xlimValue, ylim=ylimValue)
	dev.off()
}



print("You made it to the end of the script with 0 problems! Good work dude!")

