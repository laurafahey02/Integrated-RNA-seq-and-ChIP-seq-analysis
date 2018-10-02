
library(edgeR)
counts <-read.table(file = "gene_count_matrix.csv", header = TRUE, sep=",", row.names=1)
dgList <- DGEList(counts=counts, genes=rownames(counts))
# Filtering: Only keep genes which have at least 1 count per million in at least 2/3 samples. Wouldn’t have statistical power to detect difference in expression of these genes.
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 3)
dgList <- dgList[keep,]
summary(cpm(dgList))
# trimmed mean of M values (TMM) normalisation
dgList <- calcNormFactors(dgList, method="TMM")
#plotMDS(dgList, main = "MDS Plot for Count Data", labels = colnames(dgList$counts))
sampleType <- c("C", "T","C","T","C","T")
sampleReplicate <- c("S1", "S1", "S2", "S2", "S3", "S3")
designMat <- model.matrix(~sampleReplicate + sampleType)
dgList <- estimateGLMCommonDisp(dgList, design=designMat) #Estimate common (mean) dispersion
dgList <- estimateGLMTrendedDisp(dgList, design=designMat) #Estimate trended disperions
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat) #Estimate tag(gene)wise dispersion
#plotBCV(dgList, main="BVC PLOT")
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=4) # Performs likelihood ration test. sampleType is 3rd column of designMat
edgeR_result <- topTags(lrt, n = Inf, adjust.method = "BH") # inf=infinity # Benjamini–Hochberg (BH) - method which controls the false discovery rate is used to correct for multiple testing.
sig_results <- edgeR_result$table[edgeR_result$table$FDR<0.01,]
deGenes <- rownames(edgeR_result$table[edgeR_result$table$FDR < .01],)
#plotSmear(lrt, de.tags=deGenes)

### Convert ensembl IDs to gene symbols which are required by Beta ###
library("biomaRt")
library("biomaRt")
sig_results <- edgeR_result$table[edgeR_result$table$FDR<0.01,]

grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
deGenes_sym <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"), filters = 'ensembl_gene_id', values = rownames(sig_results), mart = grch37)
df <- deGenes_sym[!(is.na(deGenes_sym$hgnc_symbol) | deGenes_sym$hgnc_symbol==""), ]
write.table(df, "deGenes_sym.txt", row.names = F, quote = F)
