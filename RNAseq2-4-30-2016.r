# install.packages("gplots")
library(gplots)
# install bioconductor packages 
#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite(c( "GO.db","GSEABase","GOstats","biomaRt"))
library(biomaRt)
library("GO.db")
library("GSEABase")
library("GOstats")

setwd("C:/Users/Xijin.Ge/Google Drive/research/collaboration/Sen Subramanian")

hclust2 <- function(x, method="average", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))
library(gplots)
myheatmap <- function (x,n=-1) {
if(n == -1) n=dim(x)[1]
geneSD = apply(x,1,sd)
x = x[order(-geneSD),]
# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
labRow = rownames(x)
if(dim(x)[1] > 300) labRow = F 
hy <-  heatmap.2(x, distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
#,Colv=FALSE,
,labRow = labRow
,key=F
,margins = c(6, 8)
)}
myheatmap3 <- function (x,bar = bar) {
if(length(bar) == 0 | length(bar) != dim(x)[1]) {
cat( "length does not match!"); stop() }
# this will cutoff very large values, which could skew the color 
x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff

nclass = length(unique(bar))
set.seed(2)
ncolors = sample( rainbow(nclass) )

heatmap.2(x,#distfun = dist2,hclustfun=hclust2,
 col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
,dendrogram ="column"
,key=F, labRow = F
,Rowv = FALSE
,RowSideColors = ncolors[bar]
,margins = c(6, 8)
)
}
myheatmap2 <- function (x,bar = bar) {
if(length(bar) == 0 | length(bar) != dim(x)[1]) {
cat( "length does not match!"); stop() }
# this will cutoff very large values, which could skew the color 
x=as.matrix(x)-apply(x,1,mean)
cutoff = median(unlist(x)) + 3*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 3*sd (unlist(x)) 
x[x< cutoff] <- cutoff
x = x[length(bar):1,]; bar = bar[length(bar):1] # reverse order
nclass = length(unique(bar))
set.seed(2)
ncolors = sample( rainbow(nclass) )
# use heatmap as it is faster
hy <-  heatmap(x,distfun = dist2,hclustfun=hclust2
, col=greenred(75),labRow = F,Rowv = NA
,RowSideColors = ncolors[bar],margins = c(6, 8))
}


x = read.csv("expression.csv", header=T)
dim(x)
x = x[,-c(2:5)]

ix = grep("q_val.*",colnames(x))
sig = x[,ix ]
x <- x[, -ix]

if(1) {   # remove genes not significantly different from control
qmin = apply(sig, 1, min)    
x <- x[which(qmin <0.01),]   # get rid of genes with FDR > 0.01
dim(x) 
}
ix = grep("Fold_.*",colnames(x))
fold = x[,c(1,ix )]
x <- x[,-ix]

ix <- grep("sig_fold.*",colnames(x))
fold = x[,c(1,ix )]
x <- x[,-ix]

ix <- grep("FPKM.*",colnames(x))
FPKM <- x[,c(1,ix)]
x <- x[,-ix]

summary(x[,2])

Gmax = apply(x[,2:dim(x)[2]], 1, max)    # max by gene
x = x[order(x[,1],Gmax, decreasing=TRUE),]     # sort by gene symbol then max
x <- x[ !duplicated(x[,1]), ]   # remove duplicated gene symbol, keep highest
dim(x)

Gmax = apply(x[,2:dim(x)[2]], 1, max)    # max by gene
x = x[which(Gmax> 1),]   # remove lowly expressed
dim(x)   #

rownames(x) = x[,1]   # rownames ---> gene symbol 
x = x[,-1]

x= log2(x+1)
boxplot(x)

hist(x[,2])
 
plot( hclust2(dist2(t(x)) ), main="All samples" )

x= x[-1,] # remove novel gene
dim(x)
x_raw = x; # keep a copy of the whole data
write.csv(x,"genes12813.csv")
######################################

#myheatmap(x,12000)
bound(cor(x),2)
myheatmap(cor(x))

colnames(x)
myheatmap(x[,1:12],4000)  # just nodule
myheatmap(x[,13:24],4000)  # just LR
myheatmap(x[,1:12],100)  # just nodule

############################################################
##  investigating fold change values
x = FPKM
colnames(x) = gsub("FPKM_","",colnames(x))
x = x[-1,] # novel gene
Gmax = apply(abs(x[,2:dim(x)[2]]), 1, max)    # max by gene
x = x[order(x[,1],Gmax, decreasing=TRUE),]     # sort by gene symbol then max
x <- x[ !duplicated(x[,1]), ]   # remove duplicated gene symbol, keep highest
dim(x)

Gmax = apply(abs( x[,2:dim(x)[2]]), 1, max)    # max by gene
x = x[which(Gmax> 1),]   # remove lowly expressed
dim(x)   #


rownames(x) = x[,1]   # rownames ---> gene symbol 
x = x[,-1]
x = log2(x+1)   # add 1 to all FPKM
boxplot(x)

if(0) {
Gmax = apply(abs( x[,2:dim(x)[2]]), 1, max) - apply(abs( x[,2:dim(x)[2]]), 1, min)   # max by gene
x = x[which(Gmax> 2),]   # remove lowly expressed
dim(x)   #
}

x = x[order(-apply(x,1,sd)),]

barplot( apply(x,1,sd))

x_raw2 = x

#####################################################

x = x_raw2[1:5000,] # just the top 5000 genes

if(1){ # remove control samples
x = x_raw2[,c(2,4,6,8)] # just the top 5000 genes
#Gmax = apply(x, 1, max) - apply( x, 1, min)   # max by gene
#x = x[which(Gmax> 1.5),]   # remove lowly expressed

x = x[order(-apply(x,1,sd)),]
x = x[1:5000,]
dim(x) 

}

if(0){ # using original replicates
x = x_raw
ix = grep("AB",colnames(x))
x = x[,-ix]
x = x[order(-apply(x,1,sd)),]
x = x[1:5000,]  
}

###########################################
## just two class
if(0) {
x = read.csv("expression.csv", header=T)
dim(x)
x = x[,-c(2:5)]

ix = grep("q_val.*",colnames(x))
sig = x[,ix ]
x <- x[, -ix]
ix = grep("Fold_.*",colnames(x))
fold = x[,c(1,ix )]
x <- x[,-ix]

ix <- grep("sig_fold.*",colnames(x))
fold = x[,c(1,ix )]
x <- x[,-ix]

ix <- grep("FPKM.*",colnames(x))
FPKM <- x[,c(1,ix)]
x <- x[,-ix]
  
x <- x[which(sig[,1] <0.05),]   # get rid of genes with FDR > 0.01
dim(x) 

x = x[,1:7]

Gmax = apply(x[,2:dim(x)[2]], 1, max)    # max by gene
x = x[order(x[,1],Gmax, decreasing=TRUE),]     # sort by gene symbol then max
x <- x[ !duplicated(x[,1]), ]   # remove duplicated gene symbol, keep highest
dim(x)
rownames(x) = x[,1]
x=x[,-1]

fold = apply(x[,4:6],1,mean) -apply(x[,1:3],1,mean)
x <- x[which(fold >1),]   # at least 2 fold
}
#########################




#myheatmap(x)
##  kNN
#L1 norm normalization
x = 100* x / apply(x,1,sum)
#colnames(x) = gsub(".*_","",colnames(x))
# determining number of clusters
k = 25
wss <- (nrow(x)-1)*sum(apply(x,2,var))
  for (i in 2:k) wss[i] <- sum(kmeans(x,centers=i,iter.max = 30)$withinss)
plot(1:k, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")

# k=42  # inluding noncoding genes
k=10
set.seed(2)
cl = kmeans(x,k,iter.max = 20)
tem = cl$size
names(tem) = 1:k
barplot(sort(cl$size,decreasing=T))
sort(tem,decreasing=T)
myheatmap(cl$centers)

hc <- hclust2(dist2(cl$centers-apply(cl$centers,1,mean) )  )# perform cluster for the reordering of samples
plot(hc)

tem = match(cl$cluster,hc$order) #  new order 

x = x[order(tem),]
bar = sort(tem)
system.time( myheatmap2(x, bar) )
cl$cluster = bar # update clusters

write.csv(cbind(cl$cluster,x), "data_for_clustering.csv")
#barplot(table(bar),xlab="Gene cluster",ylab="#Genes in cluster")

### save gene list a GMT files
a=""
for( i in 1:k) {
   cat (i,"\n")
   a = paste(a,i,"\t\t",sep="")
    b = toupper( rownames(x)[which(cl$cluster == i)]  )
	b = paste(b, collapse ="\t",sep="")
	a = paste(a,b,"\n",sep="")
}
write(a,"Clusters.gmt")

results = overlap("Clusters.gmt","Soy_GO.txt" )

results$GO = Ontology( as.character( results[,3]) )
results$Term = Term(as.character (results[,3]) )
results[!duplicated(results[,1]),-6]
write.csv(results,"GO_enrichment_manual.csv",row.names=F)
#############################################################################
## retrieve Gene Onotology info from Ensembl and GO.db package
#############################################################################
myMart = useMart("plants_mart", dataset = "gmax_eg_gene",host="plants.ensembl.org")

 cat("Reading existing file for gene ontology...\n"); 
 
## retrieve GO annotation from ENSEMBL
 cat("Retieving info from ENSEMBL...\n")
 annot<-getBM(attributes = c("go_accession","go_linkage_type","ensembl_gene_id"),
               filters = 'biotype',values = "protein_coding",    # only focus on protein coding genes
               mart = myMart)
 goframeData = annot
 goframeData = goframeData[which(goframeData[,2]!=""),]  # if evidence is missing remove.
 goframeData = goframeData[which(goframeData[,3]!="NA"),]  # if evidence is missing remove.
 head(goframeData)
 goFrame=GOFrame(goframeData)
 goAllFrame=GOAllFrame(goFrame)
 #build gene set collection
 gsc <-  GeneSetCollection(goAllFrame, setType = GOCollection())
 toGmt(gsc,"Soy_GO.txt");  # save gene set to file  toGmt() is a GSEABase method
 
###########################################################################
### GO  Enrichment analysis 
###########################################################################

# define uinverse as all genes with GO annotation
universe =  toupper( rownames(x) )
genesWithGO = unique(goframeData[,3])
universe  = intersect(universe, genesWithGO) 
#universe = genesWithGO
length( universe)


first =TRUE
for( i in 1:k) {
   cat (i,"\n")
    b = toupper( rownames(x)[which(cl$cluster == i)]  )
	b = intersect(b,universe)
    if(length(b) <15) next; 
	#GOStats
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
   geneSetCollection=gsc,
   geneIds = b,
   universeGeneIds = universe,
   ontology = "BP",
   pvalueCutoff = 0.05,
   conditional = FALSE,
   testDirection = "over")
 
 #And finally we can call hyperGTest in the same way we always have before.
 Over <- hyperGTest(params)
 head(summary(Over))
 if( dim(summary(Over))[1] == 0 ) next;  # if no significant  
  r1 = data.frame( summary(Over) )
  r1[,8]= i
  if (first) {result = r1; first = FALSE }  else
  { result = rbind(result,r1) }
}

colnames(result)[8]="Cluster"
result$FDR = (p.adjust(result$Pvalue,method="fdr"))
result = result[order(result$Cluster, -result$FDR),]
head(result)
result[!duplicated(result$Cluster),]
write.csv(result,"GO_enrichment.csv")



###########################################################
# compute overlap 
################
# This program computes overlaps between two gene sets files in the GMT format. 
overlap <- function (file1, file2, total_elements = 35000, minFDR=0.05, minP=0.01 ) {
Pval_cutoff <- .01
#total_elements <- 30750 #24389 #31269   # This needs to be changed
#total_elements <- dim(genes)[1]
#total_elements <- length(universe)
#total_elements <- 35000
Min_overlap <- 1
Depeletion =FALSE    # Change to TRUE if want to identify depletion
minSetSize = 3; 
#file1 = "Clusters.gmt";
#file2 = "KEGG_Glycine_max_soybean_EnsemblID.gmt";
#file2 ="Soy_GO.txt"

# Read in the first file 
x <- scan(file1, what="", sep="\n")
#x <- gsub("\t\t.","",x)     # GMT files saved by Excel has a lot of empty cells "\t\t\t\t"   "\t." means one or more tab
#x <- gsub(" ","",x)  # remove white space
#x <- toupper(x)    # convert to upper case

#----Process the first file
# Separate elements by one or more whitespace
y <- strsplit(x, "\t")
# Extract the first vector element and set it as the list element name
names(y) <- sapply(y, `[[`, 1)
#names(y) <- sapply(y, function(x) x[[1]]) # same as above
# Remove the first vector element from each list element
y <- lapply(y, `[`, -c(1,2))
#y <- lapply(y, function(x) x[-1]) # same as above
# remove duplicated elements
for ( i in 1:length(y) )  y[[i]] <- unique(y[[i]])
# check the distribution of the size of gene lists sapply(y, length) hold a vector of sizes
if( max( sapply(y,length) ) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
y = y[which(sapply(y,length) > minSetSize)]  # gene sets smaller than 10 is ignored!!!

#---Process second file
x2 <- scan(file2, what="", sep="\n")
#x2 <- gsub("\t\t.","",x2)  
#x2 <- gsub(" ","",x2)
#x2 <- toupper(x2)
# Separate elements by one or more whitepace
y2 <- strsplit(x2, "\t")
# Extract the first vector element and set it as the list element name
names(y2) <- sapply(y2, `[[`, 1)
y2 <- lapply(y2, `[`, -c(1,2))
# remove duplicated elements
for ( i in 1:length(y2) )  y2[[i]] <- unique(y2[[i]])
if( max( sapply(y2,length)) <5) cat("Warning! Gene sets have very small number of genes!\n Please double check format.")
y2 = y2[which(sapply(y2,length) > minSetSize)]  # gene sets smaller than 10 is ignored!!!

#initialize a matrix to hold results
results <- matrix(ncol=8)
results_count <- matrix(ncol=length(y2), nrow=length(y))
Pval_table <- matrix(ncol=length(y2), nrow=length(y))
rownames(results_count) = names(y);   colnames(results_count) = names(y2)
rownames(Pval_table) = names(y);   colnames(Pval_table) = names(y2)
colnames(results) <- c("Cluster/Set1","n1","Set2","n2","#overlaped","Genes","Pval","FDR")
Pvalues = c()
Ntest=0;
#compute overlaps and P values
for ( i in 1:length(y) ) {
   cat(paste(i,"/",length(y),"\n"))
   if( length(y[[i]]) ==0    ) next; 
   for(j in 1:length(y2) ) {
       if( length(y2[[j]]) ==0    ) next; 
       ooo <- intersect(y[[i]], y2[[j]])
	   if (length(ooo) == 0) next; 
       results_count[i,j] = length(ooo)
       Ntest = Ntest+1;
       # if (length(ooo) <= Min_overlap) next  this may cause problem for detecting depletions
       xx <- length(ooo)
       mm <- length(y[[i]])
       nn <- total_elements - mm
       kk <- length(y2[[j]])
       if(nn<0 || kk> total_elements) { cat(paste(" Check total elements")); stop(); }
	   if (xx> 100) ooo <- ooo[1:100] # only keep the first 100 genes, if there is more
        # Note that it needs to be i-1 http://stats.stackexchange.com/questions/16247/calculating-the-probability-of-gene-list-overlap-between-an-rna-seq-and-a-chip-c
       Pval_deplete=phyper(xx,mm,nn,kk, lower.tail=TRUE );
       Pval_enrich=phyper(xx-1,mm,nn,kk, lower.tail=FALSE );
       #Pval_table[i,j] = (-log10(Pval_enrich));   # if enrichment a positive number otherwise negative
       #if(Pval_deplete < Pval_enrich) Pval_table[i,j] = -(-log10(Pval_deplete))^(1/2);
	 if(Depeletion)   Pval=Pval_deplete  else Pval=Pval_enrich;
       if( Pval < Pval_cutoff ) {
           newO =c(names(y)[i],mm,names(y2)[j],kk,xx,paste(ooo,collapse=";"),Pval,"NA")
           results <- rbind(results,newO) 
   	     Pvalues = c(Pvalues,Pval)
        }
   }
}

results <- results[-1,]   # remove the first row as it is empty
if(dim(results)[1] <1) { cat("\nNo significant overlap found! Please "); stop() }
results <- as.data.frame(results)  #convert to data frame
results$FDR = (p.adjust(Pvalues,method="fdr",n=Ntest))
results <- results[ order( as.numeric(results[,1]),results[,8])  ,]  # sort according to FDR

results <- results[ which(results[,8]<minFDR)  ,]  # filter by FDR

return(results)
#results = results[,-6] # remove genes as there are too many


}

########################################################


myPGSEA  <- function (exprs, cl, range = c(25, 500), ref = NULL, center = TRUE, 
    p.value = 0.005, weighted = TRUE, nPermutation=100, enforceRange = TRUE, ...) 
{
    if (is(exprs, "ExpressionSet")) 
        exprs <- exprs(exprs)
    if (!is.list(cl)) 
        stop("cl need to be a list")
    if (!is.null(ref)) {
        if (!is.numeric(ref)) 
            stop("column index's required")
    }
    if (!is.null(ref)) {
        if (options()$verbose) 
            cat("Creating ratios...", "\n")
        ref_mean <- apply(exprs[, ref], 1, mean, na.rm = TRUE)
        exprs <- sweep(exprs, 1, ref_mean, "-")
    }
    if (center) 
        exprs <- scale(exprs, scale = FALSE)         # column centering is done
    results <- matrix(NA, length(cl), ncol(exprs))
    rownames(results) <- names(cl)
    colnames(results) <- colnames(exprs)
    mode(results) <- "numeric"
	Setsize = c(rep(0,length(cl)))     # gene set size vector
	mean2 = c(rep(0,length(cl)))     # mean of the range of means 
	meanSD = c(rep(0,length(cl)))     # SD of the range of means	
    if (is.logical(p.value)) 
        { p.results <- results; mean.results <- results;}
    for (i in 1:length(cl)) {              # for each gene list
		cat("\nProcessing gene set",i);
        if (class(cl[[i]]) == "smc") {
            clids <- cl[[i]]@ids
        }
        else if (class(cl[[i]]) %in% c("GeneColorSet", "GeneSet")) {
            clids <- cl[[i]]@geneIds
        }
        else {
            clids <- cl[[i]]
        }
        if (options()$verbose) 
            cat("Testing region ", i, "\n")
        ix <- match(clids, rownames(exprs))
        ix <- unique(ix[!is.na(ix)])
        present <- sum(!is.na(ix))
		Setsize[i] <- present 
        if (present < range[1]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too small-", 
                  present, ",\n")
            next
        }
        if (present > range[2]) {
            if (options()$verbose) 
                cat("Skipping region ", i, " because too large-", 
                  present, "\n")
            next
        }
        texprs <- exprs[ix, ]           # expression matrix for genes in gene set
        if (any(is.na(texprs))) 
            cat("Warning - 'NA' values within expression data, enrichment scores are estimates only.\n")
        if (!is.matrix(texprs)) 
            texprs <- as.matrix(texprs)
                            
        stat <- try(apply(texprs, 2, t.test, ...))
		means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
		ps <- unlist(lapply(stat, function(x) x$p.value))
        stat <- unlist(lapply(stat, function(x) x$statistic))
        p.results[i, ] <- ps
		mean.results[i,] <- means
        results[i, ] <- as.numeric(stat)
		
		# permutation of gene sets of the same size
		if(nPermutation > 2 )  { # no permutation if <=2
			meansRanges = c(0, rep(nPermutation))
			for( k in 1:nPermutation ) {
				ix <- sample.int( dim(exprs)[1], length(ix) )
				texprs <- exprs[ix, ] 
				means <- try(apply(texprs, 2, mean,trim=0.1))   # trim mean
				meansRanges[k] = dynamicRange(means)
			}
			mean2[i] = mean(meansRanges)
			meanSD[i]= sd(meansRanges,na.rm=TRUE)   # NA are removed before calculating standard deviation
		}
    }
    return(list(results = results, p.results = p.results, means = mean.results, size=Setsize, mean2=mean2, meanSD=meanSD))
    
}

dynamicRange <- function( x ) {
y = sort(x)
   if(length(x)>=4)  k =2 else k =1;
   return( y[length(x)-k+1] - y[k]) 
}




#############################################
## prepare data
library(PGSEA)
#read GMT file.
gmt=readGmt("Soy_GO.txt" )   # read GMT file
gmt=readGmt("KEGG_Glycine_max_soybean_EnsemblID.gmt" )   # read GMT file

x=x_raw[,1:6] # just emerging nodule
subtype = colnames(x)  # sample groups
subtype = gsub("_.*","",subtype)

Pvalue = 0.01  # cut off to report in PGSEA. Otherwise NA
Pval_pathway = 0.001   # cut off for P value of ANOVA test  to writ to file 
top = 60  # number of pathways to show

if(0) {  # remove control samples
x_raw =x
ix = grep("AB",colnames(x))
x = x_raw[,-ix]
subtype = subtype[-ix]
}
#---------------------Pathways using PGSEA
rownames(x) = toupper(rownames(x))  # this is for consistency with pathway databases
x = x - rowMeans(x)   # centering by mean
pg = myPGSEA (x,cl=gmt,range=c(15,3000),p.value=TRUE, weighted=FALSE,nPermutation=100)
pg2 = pg$results;
pg2 = pg2[rowSums(is.na(pg2))<ncol(pg2) ,]  # remove se/wrts with all missing(non-signficant)
if (dim(pg2)[1] < 2 ) {  cat("Pathway matrix pg2 has less than 2 rows"); next;}
pg.pvalues = pg$p.results
pg.pvalues = pg.pvalues[rowSums(is.na(pg.pvalues))<ncol(pg2) ,]  # remove se/wrts with all missing(non-signficant)
pg.means = pg$means
pg.means = pg.means[rowSums(is.na(pg$results))<ncol(pg2) ,]  # remove se/wrts with all missing(non-signficant)
setSize =  pg$size[rowSums(is.na(pg$results))<ncol(pg2)] 
mean2 =   pg$mean2[rowSums(is.na(pg$results))<ncol(pg2)]
meanSD = pg$meanSD[rowSums(is.na(pg$results))<ncol(pg2)]
rownames(pg2) = gsub(".*_MM_","",rownames(pg2))
rownames(pg2) = gsub(" .*","",rownames(pg2))

max= c(rep(0,rep(dim(pg2)[1])))
nSample = dim(x)[2]
for( j in 1:dim(pg2)[1] )
{  # use the difference between the 2nd hights and 2nd lowest to indicate change.
   tem = pg.means[j,];  
   if(nSample>=4)  {  max[j] = sort(tem)[nSample-1] - sort(tem)[2]} 
   else  max[j] = max(tem)-min(tem)   # small sample.
}

Tmax= c(rep(0,rep(dim(pg2)[1])))
for( j in 1:dim(pg2)[1] )
{  # use the difference between the 2nd hights and 2nd lowest to indicate change.
   tem = pg2[j,];  
   if(nSample>=4)  {  Tmax[j] = sort(tem)[nSample-1] - sort(tem)[2]} 
   else  Tmax[j] = max(tem)-min(tem)   # small sample.
}

pg4=pg2
pg4[ is.na(pg4) ] <- 0   # replace NA with 0
pg3 = abs(pg4);  
best = max(pg3);
result <- pg4;
if(length(subtype) < 4 || length(unique(subtype)) <2) { 
 result = result[ order(-result[,1])   ,]
 write.csv(result,paste(wd,"/",GMTfile,"/",platform,"_",geoid,".csv",sep=""))
 cat("No P value calculated due to the lack of group definition.\n")
 next;    #skip to the next GSE
} 


cat("Computing P values using ANOVA\n");
Pvalues = c()   # compute P values using analysis of variance
for( k  in 1:dim(pg2)[1] ) {
  P=summary( aov(pg2[k,]~subtype) )[[1]][["Pr(>F)"]][1] 
  Pvalues = c(Pvalues,P)
  }
if(min(Pvalues) > Pval_pathway ) {  cat("nothing significant. Skip to the next GSE"); next; }  

NsigT = rowSums(pg.pvalues<Pvalue)
Zscore =(max-mean2)/meanSD

#result = cbind(NsigT,max,Tmax,setSize,mean2, meanSD,result)
result=cbind( as.matrix(Pvalues),Zscore,NsigT,result); 
colnames(result)[1]="Pval"
#rank2 = dim(pg2)[1]+1-rank(max) + sqrt(dim(pg2)[1]+1- rank(Tmax) )   # sum of two ranks, 1 means biggest
result = result[which(result[,1] < Pval_pathway)    ,]
result = result[which(result[,3] >2)    ,]

result = result[ order(-result[,2])   ,]
result = result[,-3]
write.csv(result,"tem.csv")

pg2 = result[,-2]

# when there is only 1 left in the matrix pg2 becomes a vector
if(sum( Pvalues<Pval_pathway) == 1) { pg3 = t( as.matrix(pg2));pg3 = rbind(pg3,pg3);} else
{ if(dim(pg2)[1] > top ) {  pg3 = pg2[1:top,]; } else {  pg3 = pg2[1:dim(pg2)[1],];  } }

a=sprintf("%-1.0e",pg3[,1])
rownames(pg3) = paste(a,rownames(pg3),sep=" ")
pg3 =pg3[,-1]

smcPlot(pg3,factor(subtype),scale = c(-best, best), show.grid = F, margins = c(3,1, 13, 13), col = .rwb,cex.lab=0.5)
# myheatmap(pg2[,-1])

hc <- hclust2(dist2(pg2[,-1]) ) # perform cluster for the reordering of samples

pg2 = as.data.frame(pg2)
pg2$GO = Ontology(gsub(" ", "", rownames(pg2)))
pg2$Term = Term(gsub(" ", "", rownames(pg2)))

pg2 = pg2[hc$order,]

write.csv(pg2,"PGSEA_result.csv")

