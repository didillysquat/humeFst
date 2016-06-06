args = commandArgs(trailingOnly=TRUE)
baseLocation = args[1]
for (zipFile in list.files(paste(baseLocation, '/Rbins', sep=''), pattern="\\.zip$")){
  print(paste('Installing module: ', zipFile, sep=""))
  install.packages(paste(baseLocation, '/Rbins/', zipFile, sep=""), lib=paste(baseLocation,'/Rlibs/',sep=""), repos=NULL)
}
print('Groovy! All libraries were installed successfully. Continuing SymbiodiniumType')
