# glabrata_DGE.R
# Evan Lavelle
# Created: May 2019
library("rtracklayer")
library("DESeq2")
library("ggplot2")

source("annotation_search.R")

# Fixed parameters:

MIN.EXPRESSION.LEVEL <- 50

# Output files of figures:
file.pca <- "snail-parasite-pca.pdf"
file.top.genes <- "snail-parasite-top-genes.pdf"

# loads transcript counts
theCountData <- read.csv("gene_count_matrix.csv", header=TRUE)

# captures row names
names <- theCountData[,1]
# removes names column from file
cleanCounts <- theCountData[,-1]
librarySizes <- colSums(cleanCounts)

# ----------------------------------------------------------------
# Gene filtering: removing genes lowly expressed in all conditions
#   This step avoids unnecessary issues to deal with NAs in modeling.
#   It will also slightly increase the runtime with reduced number
#   of genes to work with.

keep <- rowSums(cleanCounts) >= MIN.EXPRESSION.LEVEL
cleanCounts <- cleanCounts[keep, ]
cat(sum(keep), "out of", length(keep), 
    "genes passed filtering at a minimum total counts of", MIN.EXPRESSION.LEVEL, "\n")

# sets row names as attribute:
row.names(cleanCounts) <- names[keep]

# -----------------------------------------------------------
# Define experimental design
#
# just to fill dataframe with
rVector <- rep(c("RR"), times = 18)
s1Vector <- rep(c("S1"), times = 9)
s2Vector <- rep(c("S2"), times = 9)
# because challenged are out of order compared to SRR no.
genotypeVector <- c(rVector,s1Vector,s2Vector,"RR","RR","RR","S1","S2","S2","RR","RR","RR","S1","S2","S2")

controlVector <- rep(c("control"), times = 36)
challengeVector <- rep(c("challenge"), times = 12)
treatmentVector <- c(controlVector, challengeVector)

zeroVector <- rep(c("0 hrs"), times = 36)
twoVector <- rep(c("2 hrs"), times = 6)
sixVector <- rep(c("6 hrs"), times = 6)
timeVector <- c(zeroVector,twoVector,sixVector)

controlTypeVector <- rep(c("resistant","susceptible"), each = 18)
challengeTypeVector <- rep(c("resistant","susceptible"), each = 3)
typeVector <- c(controlTypeVector,rep(challengeTypeVector, times = 2))

readVector <- rep(c("single"), times = 12)
readVector <- c(readVector,rep(c("paired"), times = 6))
readVector <- c(readVector,rep(c("single"), times = 6))
readVector <- c(readVector,c("paired","paired","paired"))
readVector <- c(readVector,rep(c("single"), times = 6))
readVector <- c(readVector,c("paired","paired","paired"))
readVector <- c(readVector,rep(c("single"), times = 12))

samples <- data.frame(Genotype = genotypeVector, 
                      Treatment = treatmentVector,
                      Timepoint = timeVector,
                      Phenotype = typeVector,
                      Sequencing_Type = readVector,
                      Lib_Size = librarySizes)

samples$Treatment <- relevel(samples$Treatment, "control")
samples$Phenotype <- relevel(samples$Phenotype, "susceptible")
samples$Sequencing_Type <- relevel(samples$Sequencing_Type, "single")

# ------------------------------------------------------
# DESeq2 analysis

dataSet <- DESeqDataSetFromMatrix(countData = cleanCounts, colData = samples, 
                                  design = ~Sequencing_Type + Phenotype + Treatment + Phenotype:Treatment, tidy = FALSE)

dds <- DESeq(dataSet, minReplicatesForReplace = 6) # , minReplicatesForReplace=Inf

# check for replacement
ddsreplaced <- which(rowData(dds)$replace == "TRUE")
# indeces of genes where 1 or more counts are replaced
ddsreplacedindeces <- match(names(ddsreplaced), rownames(counts(dds)))
# makes vector of maxCook values 
mcols(dds)$maxCooks <- apply(assays(dds)[["cooks"]],1,max)
# shows maxCook value for gene where replacement took place
ddsreplacedMaxes <- rowData(dds)$maxCooks[ddsreplacedindeces]

ddsResults <- results(dds, list("Phenotyperesistant.Treatmentchallenge"), cooksCutoff = .05)
summary(ddsResults)
testResults <- results(dds)
testResults <- testResults[which(testResults$padj < .05),]

ddsSigOrdered <- testResults[order(testResults$padj),]

# DEGs for susceptible strain after challenge
sus <- results(dds, contrast=c("Treatment","challenge","control"), cooksCutoff = .05)
sus <- sus[which(sus$padj < .05),]
summary(sus, alpha = .05)

# DEGs for resistant strain after challenge
res <- results(dds, contrast=list(c("Treatment_challenge_vs_control","Phenotyperesistant.Treatmentchallenge")), cooksCutoff = 0.05)
res <- res[which(res$padj < .05),]
summary(res, alpha = .05)

# intersection
exposureDEGs <- length(intersect(row.names(res),row.names(sus)))

#noFilter <- results(dds, list("Phenotyperesistant.Treatmentchallenge"), alpha = .05)
#summary(noFilter)
#noFilterOrdered <- noFilter[order(noFilter$padj),]
#grep(pattern = "MSTRG.46623", row.names(dds))
assays(dds[43703,])[["cooks"]]

# DEGs between strains in control groups
cont <- results(dds, contrast = c("Phenotype","resistant","susceptible"), alpha = .05)
cont <- cont[which(cont$padj < .05),]
#TEMPORARY REPLACEMENT FOR ddsSigOrdered
ddsSigOrdered <- cont[order(cont$padj),]
summary(cont)

# DEGs between strains in challenged groups
chal <- results(dds, contrast=list(c("Phenotype_resistant_vs_susceptible","Phenotyperesistant.Treatmentchallenge")), cooksCutoff = 0.05, alpha = .05)
chal <- chal[which(chal$padj < .05),]
summary(chal)

# intersection
strainDEGs <- length(intersect(row.names(cont),row.names(chal)))

ddsSig <- ddsResults[which(ddsResults$padj < 0.05), ]
ddsOrdered <- ddsResults[order(ddsResults$padj),]
ddsSigOrdered <- ddsSig[order(ddsSig$padj),]

# indeces in results matrix of DEGs
ddsSigOrderedIndeces <- match(row.names(ddsSigOrdered), rownames(counts(dds)))
# Cook's value maxes for DEGs
ddsSigOrderedCooks <- rowData(dds)$maxCooks[ddsSigOrderedIndeces]

# indeces in results matrix of selected DEGs
ddsSigOrderedIndeces <- match(row.names(ddsSigOrdered), rownames(counts(dds)))
# syntax to check all cookCounts of given row/gene
# assays(dds)[["cooks"]][*geneindex*,]

LOClist <- lapply(row.names(ddsSigOrdered), grep, pattern = "LOC")

Lindeces <- which((LOClist) == 1)
SigNames <- row.names(ddsSigOrdered)[Lindeces]
# takes out MSTRG id
SigNames <- sub("[^L]*", "", SigNames)
MSTRGNames <- sapply(strsplit(row.names(ddsSigOrdered),"|", fixed = TRUE),'[',1)

######################## output dataframe
productNames <- extract.annotation(SigNames)
geneNames <- c("")
# hardcoded for double-name
productNames[18] <- "uncharacterized LOC106076923|DNA endonuclease RBBP8 like"
rank <- seq(1:nrow(ddsSigOrdered))
blankVector <- vector("character", length=nrow(ddsSigOrdered)) 
Names <- c("AIG1","PARP14","PRRT1","MRPS2","HSPA1","MAP3K3","TNRC6C",
                "TUBA3","gacU","CLEC7A","STK","HSPA1","MIPEP","RBBP8")

regulation <- vector("character", length=nrow(ddsSigOrdered)) 
for(num in seq(1:nrow(ddsSigOrdered))){
  # if resistant phenotype regulation is less
  if(ddsSigOrdered$log2FoldChange[num] < 0){
    regulation[num] <- "-"
  }
  else{
    regulation[num] <- "+"
  }
}

DEGs <- data.frame(rank, blankVector, blankVector, blankVector, ddsSigOrdered$padj, regulation, stringsAsFactors = FALSE)
colnames(DEGs)[1] <- "Rank"
colnames(DEGs)[2] <- "Locus #"
colnames(DEGs)[3] <- "Annotation"
colnames(DEGs)[4] <- "Names"
colnames(DEGs)[5] <- "padj"
colnames(DEGs)[6] <- "Regulation"

rownames(DEGs) <- MSTRGNames

i <- 0
j <- 0
for (index in Lindeces){
  print(index)
  i <- i + 1
  DEGs$`Locus #`[index] <- SigNames[i] 
  DEGs$Annotation[index] <- productNames[i] 
  # skips uncharacterized genes
  if(!grepl("uncharacterized",productNames[i])){
    j <- j + 1
    DEGs$Names[index] <- Names[j]
  } 
}
# hardcode for DNA endonuclease @ 41st position
DEGs$Names[41] <- "|RBBP8"
write.csv(DEGs, "DEGs.csv", row.names = TRUE) # "f:/DEGs.csv"
write.csv(ddsSigOrdered, "/media/el239/MYDRIVE/baselineDEGs.csv", row.names = TRUE) # 
# ----------------------------------------------------------
# Perform PCA

pdf(file.pca, width=7, height=5)

pairedSamples <- which(readVector == "paired")

###before

lognormcounts <- DESeqTransform(dds)
assay(lognormcounts) <- log2(1 + counts(dds, normalized = FALSE))

data <- plotPCA(lognormcounts, ntop=nrow(dds), 
                intgroup=c("Phenotype", "Treatment", "Sequencing_Type",
                           "Lib_Size"), 
                returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, fill=Sequencing_Type, size=Lib_Size, 
                 color=Treatment, shape=Phenotype))+
  scale_shape_manual(values=c(22, 24)) + 
  scale_fill_manual(values=c("transparent", "gray")) + 
  geom_point(alpha=100/100, stroke=1.25) + 
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("principal component 1: ", percentVar[1], "% variance"))+
  ylab(paste0("principal component 2: ", percentVar[2], "% variance"))+
  labs(title = "Variance overview of log raw counts", subtitle="Before sequencing type effect removal")

lognormcounts <- DESeqTransform(dds)
assay(lognormcounts) <- log2(1 + counts(dds, normalized = TRUE))

data <- plotPCA(
  lognormcounts, ntop=nrow(dds), returnData=TRUE,
  intgroup=c("Phenotype", "Treatment", "Sequencing_Type", "Lib_Size")
)

percentVar <- round(100 * attr(data, "percentVar"))

ggplot(data, aes(PC1, PC2, fill=Sequencing_Type, size=Lib_Size,
                 color=Treatment, shape=Phenotype))+
  scale_shape_manual(values=c(22, 24)) + 
  scale_fill_manual(values=c("transparent", "gray")) + 
  geom_point(alpha=100/100, stroke=1.25) + 
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("principal component 1: ", percentVar[1], "% variance"))+
  ylab(paste0("principal component 2: ", percentVar[2], "% variance"))+
  labs(title = "Variance overview of log normalized counts", subtitle="Before sequencing type effect removal")

# -------------------------------------------------------
# sequencing type effect removal
adjustVector <- coef(dds)[,2]
adjustVector <- unname(adjustVector)
# changes NA to 0s
adjustVector <- replace(adjustVector, is.na(adjustVector), 0)

lognormcounts <- DESeqTransform(dds)
assay(lognormcounts) <- log2(1 + counts(dds, normalized = TRUE))

removed <- lognormcounts
for (index in which(samples$Sequencing_Type == "paired")){
  # print(index)
  assay(removed)[,index] <- assay(lognormcounts)[,index] - adjustVector
}

### after
data1 <- plotPCA(
  removed, ntop=nrow(removed), returnData=TRUE,
  intgroup=c("Phenotype", "Treatment", "Sequencing_Type", "Lib_Size"))

percentVar1 <- round(100 * attr(data1, "percentVar"))

ggplot(data1, aes(PC1, PC2, fill=Sequencing_Type, size=Lib_Size,
                  color=Treatment, shape=Phenotype))+
  scale_shape_manual(values=c(22, 24)) + 
  scale_fill_manual(values=c("transparent", "gray")) + 
  geom_point(alpha=100/100, stroke=1.25) + 
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21))) +
  xlab(paste0("principal component 1: ", percentVar1[1], "% variance")) +
  ylab(paste0("principal component 2: ", percentVar1[2], "% variance")) +
  labs(title = "Variance overview of log normalized counts", 
       subtitle="After sequencing type effect removal")

dev.off()

######################## to establish read effect values
typeremoved <- sapply(seq(ncol(lognormcounts)), function(j) {
  if(dds$Sequencing_Type[j] == "paired") {
    assay(lognormcounts)[, j] - coef(dds)[, "Sequencing_Type_paired_vs_single"]
  } else {
    assay(lognormcounts)[, j]
  } } )
removed2 <- lognormcounts
assay(removed2) <- typeremoved

# ------------------------------------------------------------------
# Visualize top genes of significant treatment:phenotype interaction

######################## top genes
# gene <- rownames((dds)[which.min(ddsResults$padj)])

#pdf(file.top.genes, width=5, height=4)

for(k in seq(nrow(ddsSigOrdered))) {  
  gene <- rownames(ddsSigOrdered)[k]
  pval <- format(ddsSigOrdered$padj[k], digits=2)
  
  d <- plotCounts(dds, gene=gene, # replaced = TRUE, 
                  intgroup = c("Treatment", "Phenotype", "Sequencing_Type"), 
                  returnData = TRUE)
  d$Treatment <- relevel(d$Treatment, ref="control")
  
  d$Reads <- dds$Sequencing_Type
  
  p <- ggplot(d, aes(x=Phenotype, y=count, fill=Reads,
                     color=Treatment, shape=Phenotype)) +
    scale_shape_manual(values=c(22, 24)) + 
    scale_fill_manual(values=c("transparent", "gray")) + 
    geom_point(position=position_jitter(w=0.2, h=0), 
               alpha=100/100, size=3, stroke=1.25) + 
    scale_y_log10() +
    guides(fill = guide_legend(override.aes=list(shape=21)),
           color = guide_legend(override.aes=list(shape=21))) +
    labs(title=paste0("Top ", k, ": ", gene, " (P=", pval, ")"), 
         subtitle = "Normalized counts")
  
  # ggsave(paste("F:/Top",k,"gene.pdf"))
  # print(p)
  
  d <- data.frame(LogCount = assay(removed)[gene, ], 
                  Treatment = dds$Treatment,
                  Phenotype = dds$Phenotype,
                  Reads = dds$Sequencing_Type)
  d$Treatment <- relevel(d$Treatment, ref="control")
  
  p <- ggplot(d, aes(x=Phenotype, y=LogCount, fill=Reads,
                     color=Treatment, shape=Phenotype)) +
    scale_shape_manual(values=c(22, 24)) + 
    scale_color_manual(values=c("orange", "blue")) + 
    scale_fill_manual(values=c("transparent", "gray")) + 
    geom_point(position=position_jitter(w=0.2, h=0), 
               alpha=100/100, size=3, stroke=1) + 
    guides(fill = guide_legend(override.aes=list(shape=21)),
           color = guide_legend(override.aes=list(shape=21))) +
    labs(title=paste0("Top ", k, ": ", gene, " (P=", pval, ")"), 
         subtitle = "Normalized; sequencing type effect removed")
  
  # ggsave(paste("F:/Top",k,"gene, effect removed.pdf"))
  print(p)
}

dev.off()
# paste("Top",k,"gene:",gene,"effect removed",".pdf",sep = " ")
# LOC106077580 (lectin gene)

