##########ENRICHR FUNCTIONS
#First, export gene list where we take top (avg_logFC) x.genes number of genes for each comparison
	#Function for generating python input
runEnrichrAPIFromR<-function(ifile="Enrichr.geneList.txt",ofile="Enrichr.output.txt",
	libraries="Enrichr.libraries.txt",minOverlap=5,minAdjPval=0.05,sortParam="adjPval",
	EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",pythonVersion="python2.7"){
	
	command<-paste(pythonVersion,EnrichrAPI_location,"--ifile",ifile,"--ofile",ofile,"--libraries",libraries,
		"--minOverlap",minOverlap,"--minAdjPval",minAdjPval,"--summarize","--sort",sortParam,sep=" ")
	system(command,intern=T)
}


#This function takes a DEGs list from seurat and formats it for EnrichrAPI and then runs Enrichr
#If ngenes is a number, that many of genes are taken from each cluster/module. If ngenes is NULL, then all genes in the table are used	
#cluster: /Seibold/home/jacksonna/enrichrAPI.py (must have loaded python/2.7.12 - might need to be on the fat node to load python 2 and R/3.5.1 simultaneously)
#laptop: /usr/local/bin/enrichrAPI.py
DEG_enrich_now<-function(DEG_table,dataset,ngenes=NULL,Enrichr_dir="Enrichr",
	enrichr.libraries=c("GO_Cellular_Component_2018","GO_Biological_Process_2018","GO_Molecular_Function_2018",
	"KEGG_2019_Human","Reactome_2016"),minOverlap=3,minAdjPval=0.2,EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",
	pythonVersion="python2.7"){

	if(!dir.exists(Enrichr_dir)){system(paste("mkdir ",Enrichr_dir,sep=""))}
	DEG_table<-read.table(DEG_table,header=T,stringsAsFactors=F)
	if(!is.null(ngenes)){
		d <- data.table(DEG_table, key="comparison")
		DEG_table_top<-d[,head(.SD, ngenes), by=comparison]
		DEG_table_top<-as.data.frame(DEG_table_top)[c(2,1)]
		write.table(DEG_table_top,sep="\t",quote=FALSE,col.names=c("gene","module"),row.names=FALSE,
			file=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep=""))
	}else{
		write.table(DEG_table[,c(1,grep("comparison",colnames(DEG_table)))],sep="\t",quote=FALSE,col.names=c("gene","module"),row.names=FALSE,
			file=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep=""))
	}

	#Print enrichr libraries file
	cat(paste(enrichr.libraries,collapse="\n"),file=paste(Enrichr_dir,"/Enrichr.libraries.txt",sep=""))

	#Run pythonAPI
	ifile=paste(Enrichr_dir,"/",dataset,"_Enrichr.geneList.txt",sep="")
	ofile=paste(Enrichr_dir,"/",dataset,"_Enrichr.output.xlsx",sep="")
	libraries=paste(Enrichr_dir,"/","Enrichr.libraries.txt",sep="")
	runEnrichrAPIFromR(ifile=ifile,ofile=ofile,libraries=libraries,minOverlap=minOverlap,minAdjPval=minAdjPval,
		EnrichrAPI_location=EnrichrAPI_location)
}



#Run Enrichr
#Submit tables one at a time to avoid time out errors, then stich the results together later
doEnrichOneAtATime<-function(DEG_table,dataset,ngenes=NULL,Enrichr_dir="Enrichr",
	enrichr.libraries=c("GO_Cellular_Component_2018","GO_Biological_Process_2018","GO_Molecular_Function_2018",
	"KEGG_2019_Human","Reactome_2016"),minOverlap=3,minAdjPval=0.2,EnrichrAPI_location="/Seibold/home/jacksonna/enrichrAPI.py",
	pythonVersion="python2.7"){

	for(i in 1:length(unique(DEG_table$comparison))){
		currData<-DEG_table[which(DEG_table$comparison == unique(DEG_table$comparison)[i]),]
		DEG_table_curr<-paste("Enrichr/",dataset,"_forEnrichr_",i,".txt",sep="")
		dataset_curr<-paste(dataset,"_",i,sep="")
		write.table(currData,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE,file=DEG_table_curr)
		DEG_enrich_now(DEG_table=DEG_table_curr,dataset=dataset_curr,ngenes=ngenes,Enrichr_dir=Enrichr_dir,EnrichrAPI_location=EnrichrAPI_location,
			enrichr.libraries=enrichr.libraries,minOverlap=minOverlap,minAdjPval=minAdjPval,pythonVersion=pythonVersion)
	}
	
	#Now read in all the outputs into a list and export them as a single file
	wb<-createWorkbook()
	for(i in 1:length(unique(DEG_table$comparison))){
		currPath<-paste("Enrichr/",dataset,"_",i,"_Enrichr.output.xlsx",sep="")
		a <- loadWorkbook(currPath)
		sheetNames <- names(a)
		#Bring in all sheets minus the summary
		for(j in 1:(length(sheetNames) - 1)){
	 		currData<-as.data.frame(readWorkbook(currPath,sheet = j))
	 	}
	 	#Now add the imported data to the workbook
	 	addWorksheet(wb,as.character(unique(DEG_table$comparison)[i]))
	 	writeData(wb, i, currData)
	 	widths<-c(25,80,8,8,8,12,15,100)
	 	setColWidths(wb, sheet = i, cols = 1:ncol(currData), widths = widths) #or use auto
	 	freezePane(wb,i,firstRow=T)
	 	saveWorkbook(wb, paste("Enrichr/",dataset,"_Enrichr.output.xlsx",sep=""), overwrite = TRUE)  	
	}
	
	#Remove all the single tables
	for(i in 1:length(unique(DEG_table$comparison))){
		file.remove(paste("Enrichr/",dataset,"_forEnrichr_",i,".txt",sep=""))
		file.remove(paste("Enrichr/",dataset,"_",i,"_Enrichr.output.xlsx",sep=""))
		file.remove(paste("Enrichr/",dataset,"_",i,"_Enrichr.geneList.txt",sep=""))
	}
}
