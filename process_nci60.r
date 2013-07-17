#--------------------------------EDIT---------------------------------------------
targetspro="targets.txt"
input_dir="/Users/rg353/Desktop/ICBI_NCI60"
premat="GSE_NCI60.txt"
promat="NCI60_Matrix"
image_name="NCI60_pre.to.pro_workspace"
#--------------------------------EDIT---------------------------------------------


Average.replicates<-function(targetspro,premat){
	targetspro<-read.table(file=targetspro,sep="\t",header=TRUE)
	datmat<-as.data.frame(read.table(file=premat,sep="\t",row.names=1,header=TRUE,check.names=FALSE))
	dfsubset<-datmat[,(colnames(datmat)%in%(targetspro$Filename))]
	names(dfsubset)<-targetspro[match(names(dfsubset),targetspro[,'Filename']),'Names']
	df<-data.matrix(dfsubset)
	df2<-do.call(cbind,(lapply(unique(colnames(df)), function(x) rowMeans(df[,colnames(df) == x,drop=FALSE]))))
	colnames(df2)<-unique(colnames(df))
	return(df2)
}

dataMatrix=Average.replicates(targetspro,premat)

write.table(dataMatrix,file=paste0(promat,paste0(".txt")),sep="\t",row.names=FALSE,quote=FALSE)
save(dataMatrix,file=paste0(promat,paste0(".Rda")))
save.image(file=paste0(image_name,paste0(".Rdata")))

