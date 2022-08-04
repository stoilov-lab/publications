library('Rsubread')
library('edgeR')
# library('org.Mm.eg.db')
library('biomaRt')


#these bams are produced by MATS and have the reads with undesirable flags (multimapped, wrong pair orientation) filtered out
# unfiltered data is deposited at NCBI SRA under accession PRJNA795137
files<-c('/projects/Retina/MSI/mapped/m2541/_m2541.sorted.bam',
'/projects/Retina/MSI/mapped/m2543/_m2543.sorted.bam',
'/projects/Retina/MSI/mapped/m2544/_m2544.sorted.bam',
'/projects/Retina/MSI/mapped/m2545/_m2545.sorted.bam',
'/projects/Retina/MSI/mapped/m2540/_m2540.sorted.bam',
'/projects/Retina/MSI/mapped/m2546/_m2546.sorted.bam',
'/projects/Retina/MSI/mapped/m2547/_m2547.sorted.bam',
'/projects/Retina/MSI/mapped/m2548/_m2548.sorted.bam')

#using custom modified GTF, but standard GTF from encode should work just as well
counts<-featureCounts(files,annot.ext='/projects/Spirou/annotation/Mus_musculus.GRCm38.86.spikes.ai9.gtf',strandSpecific=2, 
        isGTFAnnotationFile=TRUE,nthreads=12,allowMultiOverlap=FALSE,useMetaFeatures=TRUE,isPairedEnd=TRUE)
        
descriptions<-c(rep("MsiFL",4),rep("MsiKO",4))
names<-c('MsiFL-1','MsiFL-2','MsiFL-3','MsiFL-4','MsiKO-1','MsiKO-2','MsiKO-3','MsiKO-4')
group<-c(rep("FL",4),rep("KO",4))
# 
colnames(counts$counts)<-names

counts$targets<-cbind(names,group,descriptions)

# grab annotation from biomart
mart<- useMart(biomart='ENSEMBL_MART_ENSEMBL',dataset='mmusculus_gene_ensembl')
geneTR<-getBM(attributes=c("ensembl_gene_id","ensembl_transcript_id"),mart=mart)
geneName<-getBM(attributes=c("ensembl_gene_id","entrezgene_id","external_gene_name",'mgi_symbol','description'),mart=mart)
geneName$description<-gsub(' \\[.*\\]','',geneName$description)
m<-match(counts$annotation$GeneID,geneName$ensembl_gene_id)

# we will use this as featureData 
counts$annotation<-cbind(counts$annotation,geneName[m,])

# set up the DGElist for edgeR
# what follows is gene expression analysis per edgeR manual
raw_data<-DGEList(counts=counts$counts, genes=counts$annotation)

names(raw_data)


o<-order(rowSums(raw_data$counts[,1:4]))
head(o)
sorted_data<-raw_data[o,]

countsFilt<-apply(sorted_data$counts,1,function(x) !all(x==0))

d<-duplicated(sorted_data$genes$GeneID)
head(d)
head(d[d==TRUE])
filtered_data<-sorted_data[!d & countsFilt,]
nrow(filtered_data$counts)
nrow(raw_data$counts)

filtered_data$samples$lib.size<-colSums(filtered_data$counts)
rownames(filtered_data$counts)<-rownames(filtered_data$genes)<-filtered_data$genes$GeneID
filtered_data$genes$GeneID<-NULL
head(filtered_data$genes)
head(filtered_data$counts)


filtered_data<-calcNormFactors(filtered_data)
filtered_data$samples



pdf('MDS_plot.pdf')
plotMDS(filtered_data)
dev.off()



targets<-counts$targets
targets<-as.data.frame(targets)
rownames(targets)<-targets$names


colnames(targets)<-c('label','group','description')
targets

y<-DGEList(counts=filtered_data$counts,group=targets$group,genes=filtered_data$genes)
dim(y)


y<-estimateCommonDisp(y, verbose=TRUE)
y<-estimateTagwiseDisp(y, verbose=TRUE)

pdf('BCV_plot.pdf')
plotBCV(y)
dev.off()

et<-exactTest(y)

top<-topTags(et)
top
cpm(y)[rownames(top),]


detags<-rownames(y)[as.logical(de)]

summary(de<-decideTestsDGE(et))
detags<-rownames(y)[as.logical(de)]

pdf('Smear_plot.pdf')
plotSmear(et, de.tags=detags)
abline(h=c(-1,1),col='blue')
dev.off()

head(de)
head(et$comparison)
head(et$table)

et$table$FDR<-p.adjust(et$table$P, method='BH')
exportTable<-cbind(et$genes,et$table)[order(et$table$FDR),]


exportTable[,'Chr']<-apply(exportTable,1,function(x) strsplit(x['Chr'],';')[[1]][1])
exportTable[,'Strand']<-apply(exportTable,1,function(x) strsplit(x['Strand'],';')[[1]][1])
exportTable[,'Start']<-apply(exportTable,1,function(x) min(as.numeric(unlist(strsplit(x['Start'],';')))))
exportTable[,'End']<-apply(exportTable,1,function(x) max(as.numeric(unlist(strsplit(x['End'],';')))))
names(exportTable)[1]<-'Symbol'
write.csv(file='Floxed_vs_Msi1-2KO_Differential_expression.csv',exportTable)



## fpkm calculation for filtering mass-spec database
library(ensembldb)
library(EnsDb.Mmusculus.v79)

rpkm<-rpkm(filtered_data)
medianRPKM<-median(rpkm)
filter<-apply(rpkm, 1, function(x) median(x)>=medianRPKM)
length(filter)
sum(filter)

expressedGenes<-filtered_data$genes[filter,]
MMdb<-EnsDb.Mmusculus.v79
proteins<-proteins(MMdb, return.type="AAStringSet", filter=GeneIdFilter(rownames(expressedGenes)))
writeXStringSet(proteins,"refseq_expressed_retinal_proteins.fasta", format="fasta")



q()

