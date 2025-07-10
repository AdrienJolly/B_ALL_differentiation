library(Rsamtools)
library(DECIPHER)


what=c("qname","flag" ,"rname","strand","pos","qwidth","mapq" ,"cigar","mrnm" ,"mpos" ,"isize","seq","qual")  



# bam file already filtered for barcode reads
bam <- scanBam(file.choose(), param=ScanBamParam(tag = c("CB","CN","MA"),what=what ))




make_df=function(bam){

lst <- lapply(names(bam[[1]]), function(elt) {do.call(c, unname(lapply(bam, "[[", elt)))})

names(lst) <- names(bam[[1]])
df <- do.call("DataFrame", lst)
colnames(df) = names(bam[[1]])
nrow(df)

colnames(df)[14]="CB"
colnames(df)[15]="CN"
colnames(df)[16]="MA"


dfCells=df[which(df$CN%in%"T"),]
return(dfCells)
}


dfCells=make_df(bam)


cleanBarcodes=function(dfCells)
{
CBs=unique(dfCells$CB)
report=vector()
report["InitialCBNumber"]=length(CBs)
report["InitialMeanReadNbPerCB"]=length(dfCells$CB)/length(CBs)


library(DECIPHER)



#trim barcodes 
trimedreads=TrimDNA(dfCells$seq,leftPatterns="TGTACAACTAGGCGGCCGCTAGACTGACTGCAGTCTGAGTCTGACAG",rightPatterns="",type="both",minWidth= 26)


keptreadsIndex=as.data.frame(trimedreads[[1]])

keptreadsIndex=which(keptreadsIndex$width<=28 & keptreadsIndex$width>=26)


dfCellsKept=dfCells[keptreadsIndex,]

dfCellsKept$seq=TrimDNA(dfCellsKept$seq,leftPatterns="TGTACAACTAGGCGGCCGCTAGACTGACTGCAGTCTGAGTCTGACAG",rightPatterns="",type="sequences",minWidth= 26)



dfCellsKept$seq <- DNAStringSet(substr(dfCellsKept$seq, 1, 26))


CBclean=unique(dfCellsKept$CB)

# generate consensus barcode sequence per cell
consensusClean=DNAStringSet()

report["FullLengthCBNumber"]=length(CBclean)
report["MeanReadNbPerFullLengthCB"]=length(dfCellsKept$CB)/length(CBclean)


for(i in 1:length(CBclean))
{
  consensusClean=append(consensusClean,ConsensusSequence(dfCellsKept[which(dfCellsKept$CB%in%CBclean[i]),]$seq,threshold = 0.1))
}


# find non consensual barcodes (carrying bases ) and exclude them

seqAlphabet=alphabetFrequency(consensusClean, baseOnly=TRUE, as.prob=TRUE)

consensusCleaner=consensusClean[which(seqAlphabet[,5]==0)]

CBcleaner=CBclean[which(seqAlphabet[,5]==0)]

CBAmbiguous=CBclean[which(seqAlphabet[,5]!=0)]

consensusAmbiguous=consensusClean[which(seqAlphabet[,5]!=0)]



report["NonEmbiguousCBNumber"]=length(CBcleaner)
report["MeanReadNbPerFullLengthCB"]=length(which(dfCellsKept$CB%in%CBcleaner))/length(CBcleaner)




depths=vector(length=length(consensusCleaner))



for(i in 1:length(depths))
{
  depths[i]=length(which(dfCellsKept$seq%in%consensusCleaner[i]))
}



# filter barcodes based on number of reads per barcode

UMIdepths=vector(length=length(CBcleaner))

for(i in 1:length(UMIdepths))
{
  UMIdepths[i]=length(unique(dfCellsKept[which(dfCellsKept$CB%in%CBcleaner[i]),"MA"]))
}


threshold=4


#read count based cleaning

CBsGood=CBcleaner[which(depths>=threshold)]



consensusGood=consensusCleaner[which(depths>=threshold)]

UMIdepthsGood=UMIdepths[which(depths>=threshold)]

# identify barcodes with high similarity and merge clones of size 1 with larger clones 
#if small clone is supported by less than 3 UMI

LargerClones=names(table(consensusGood)[which(table(consensusGood)>1)])

DistMat=as.matrix(stringDist(consensusGood,method="hamming"))
CellsWithAfriend=vector()

for(i in 1:nrow(DistMat))
{
dist_1=which(DistMat[i,]%in%"1")
dist_1=c(dist_1,which(DistMat[i,]%in%"2"))
if(length(dist_1)!=0)
{
  CellsWithAfriend=c(CellsWithAfriend,i)
}
}

CellsWithAfriendSeq=consensusGood[CellsWithAfriend]

CellsWithAfriendRed=CellsWithAfriend[-which(CellsWithAfriendSeq%in%LargerClones)]


CellsWithAfriend_lowEvidence=CellsWithAfriendRed[which(UMIdepthsGood[CellsWithAfriendRed]<3)]



for(i in 1:length(CellsWithAfriend_lowEvidence))
{
  friends=unique(consensusGood[which(DistMat[CellsWithAfriend_lowEvidence[i],]%in%1)])
  LargeFriend=friends[which(friends%in%LargerClones)]
  if (length(LargeFriend)==1)
  {
    consensusGood[CellsWithAfriend_lowEvidence[i]]=LargeFriend
  }
}
result=list()

report["CBFinalNumber"]=length(CBsGood)
report["MeanReadNbPerCBFinal"]=length(which(dfCellsKept$CB%in%CBsGood))/length(CBsGood)

result[[1]]=CBsGood
result[[2]]=consensusGood
result[[3]]=report

result[[4]]=CBAmbiguous
result[[5]]=dfCellsKept
result[[6]]=UMIdepthsGood
return(result)
}

                     

barcodes=cleanBarcodes(dfCells)

# argument for the adding barcode function below
barcodedCells=cbind(barcodes[[1]],as.character(barcodes[[2]]),barcodes[[6]])



#adding barcode sequence and number of UMI per cell in coldata of a SingleCellExperiment object (sce) based on the colnames of the SingleCellExperiment
addingBarcodes=function(sce,barcodedCells)
{
uniqueBarcodes=names(table(barcodedCells[,2]))
cellBarcodes=vector(length=ncol(sce))
barcodedCellsIndexes=match(barcodedCells[,1],colnames(sce))

sce$LinBarcode=rep(NA,ncol(sce))
nas=is.na(barcodedCellsIndexes)
sce$LinBarcode[barcodedCellsIndexes[which(nas=="FALSE")]]=barcodedCells[which(nas=="FALSE"),2]
sce$LinBarcodeCounts=rep(NA,ncol(sce))
sce$LinBarcodeCounts[barcodedCellsIndexes[which(nas=="FALSE")]]=barcodedCells[which(nas=="FALSE"),3]
return(sce)
}

sce=addingBarcodes(sce,barcodedCells,barcodes[[2]])
