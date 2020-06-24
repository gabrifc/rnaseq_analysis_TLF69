## Import the count files and metadata (Sample Grouping) -----------------------
library("tximport")
txMapping <- read.csv("txMapping.txt", 
                      row.names=1, sep="\t", stringsAsFactors=FALSE)
tx2gene <- tx2gene <- read.delim("tx2gene.csv",
                                      stringsAsFactors = FALSE)
countDirectory <- "counts/"
samples <- list.files(path = countDirectory)
names(samples) <- sub(".txt", "", samples)
samples <- paste0(countDirectory, samples)
txi <- tximport(samples, 
                type = "salmon", 
                tx2gene = tx2gene)
sampleGrouping <- read.csv('sampleGrouping.csv', 
                           sep=";", 
                           row.names = 1)
all(rownames(sampleGrouping) == colnames(txi$counts))

## Specify the model and run DESeq2 --------------------------------------------
library('DESeq2')
ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = sampleGrouping,
                                   design = ~ mutation)
## Specify the base levels
ddsTxi$mutation <- relevel(ddsTxi$mutation, ref = "WT")


## Check the fitType -----------------------------------------------------------

# ddsPara <- DESeq(ddsTxi, fitType = "parametric")
# ddsLocal <- DESeq(ddsTxi, fitType = "local")
# #
# par(mfrow=c(1, 2))
# plotDispEsts(ddsPara, main = 'Parametric Fit')
# plotDispEsts(ddsLocal, main = 'Local Fit')
# par(mfrow=c(1, 1))
# #
# # ## Calculate the ratios of median absolute residue
# logParaDispGeneEst=log(mcols(ddsPara)$dispGeneEst)
# logParaDispFit=log(mcols(ddsPara)$dispFit)
# paraRatio <- abs(median(na.omit(logParaDispGeneEst))) -
#   abs(median(na.omit(logParaDispFit)))
# #
# logLocalDispGeneEst=log(mcols(ddsLocal)$dispGeneEst)
# logLocalDispFit=log(mcols(ddsLocal)$dispFit)
# localRatio <- abs(median(na.omit(logLocalDispGeneEst))) -
#    abs(median(na.omit(logLocalDispFit)))
# #
# localRatio < paraRatio

## Run DESeq2 ------------------------------------------------------------------
dds <- DESeq(ddsTxi)
resultsNames(dds)

## Mutation vs WT --------------------------------------------------------------

resMut <- results(dds,
                  name = "mutation_Mutant_vs_WT")
summary(resMut)

# ## Check that p-values are distributted as they should and we have no artifacts.
# hist(resMut$pvalue, breaks = 0:20/20,
#      col = "grey50", border = "white")
# ## Looks stepper to the right. This could be due to
# ## genes with low count reads. We can exlude them:
# hist(resMut$pvalue[resMut$baseMean > 1], breaks = 0:20/20,
#      col = "grey50", border = "white")

resMutMapped <- merge(as.data.frame(resMut), 
                      txMapping, by = 'row.names')
head(resMutMapped)

ix = which.min(resMutMapped$padj)
resMutMapped <- resMutMapped[order(resMutMapped$padj),]
head(resMutMapped)

plotCounts(ddsTxi, 
           gene = "NRRL3_02325", 
           intgroup = c("mutation"))

write.table(as.data.frame(resMutMapped), 
            file = "Mutation_kexB_vs_WT.csv",
            row.names = FALSE,
            sep = ";")

## Figure stress Genes ---------------------------------------------------------

stressGenesNames <- c("achI","agsA","agsB","agsE","chsA","chsB","chsD","chsF",
                      "csmB","ddpA","ergB","gfaA","gfaB","acuH","gelA","gelF",
                      "pcmA","dfgC","bgxA","fatD","gnaA","IplB","wscB","vcxA",
                      "CrzA","MsnA","RlmA","SrbA")
stressGenes <- c("NRRL3_06777","NRRL3_07454","NRRL3_04200","NRRL3_00248",
                 "NRRL3_04653","NRRL3_00331","NRRL3_00179","NRRL3_02932",
                 "NRRL3_06067","NRRL3_09713","NRRL3_01919","NRRL3_10711",
                 "NRRL3_08347","NRRL3_08549","NRRL3_06317","NRRL3_06843",
                 "NRRL3_10970","NRRL3_00897","NRRL3_02657","NRRL3_08145",
                 "NRRL3_03114","NRRL3_05203","NRRL3_04545","NRRL3_01896",
                 "NRRL3_10641","NRRL3_07909","NRRL3_05280","NRRL3_08408")

countData <- lapply(stressGenes, 
                    function(x) plotCounts(ddsTxi, 
                                           x, 
                                           c("mutation"),
                                           returnData = TRUE))
for(i in 1:length(countData)) countData[[i]]$gene <- stressGenesNames[i]

countData <- do.call(rbind, countData)

library('ggplot2')
ggplot(countData, 
       aes(x = mutation, 
           y = count, 
           colour = mutation,
           fill = mutation)) + 
  # scale_y_log10() +
  geom_point(size = 2) +
  # geom_dotplot(binaxis="y", 
  #              stackdir="center",
  #              dotsize = 3) +
  # stat_summary(aes(y = count,
  #                  group = mutation,
  #                  colour = mutation), 
  #              fun.y=mean,
  #              geom="line",
  #              size = 1) +
  facet_wrap(~gene, ncol = 7, scales = "free_y") +
  ggtitle(label =  'CWI stress genes',
          subtitle = 'Significance: adjusted pvalue < 0.05; Dots represent normalized read counts per library; Lines connect the means in each strain.') +
  xlab('') +
  ylab('Gene Read Counts') + 
  theme_minimal() +
  theme(legend.position = "bottom") 

## Gene Ontology Enrichment ----------------------------------------------------
library('dplyr')
library('tibble')
## Update GO terms from jgi.
goTermsFile <- 'goMapping.tab'
goMapping <- read.delim(goTermsFile, na.strings = "", stringsAsFactors = FALSE)
## We have to format the gene names
goMapping <- goMapping %>%
  mutate(GENEID = paste0('NRRL3_',
                         formatC(goMapping$X.proteinId,
                                 width = 5,
                                 format = "d",
                                 flag = "0")),
         GO = goAcc) %>%
  select(GENEID, GO)
#   
# # Last check
goMapping <- goMapping[complete.cases(goMapping),]

## We would normally download geneLengthData Needed for GOseq from the latest 
## Ensembl Release, but we add them from the txMapping file

geneLengthData <- txMapping$Length
names(geneLengthData) <- rownames(txMapping)
## Also, we have length info for all genes, so no need to filter.

## Get the significant genes
significantGenes <- as.data.frame(resMut) %>%
  rownames_to_column('gene') %>%
  filter(padj < 0.05) %>%
  column_to_rownames('gene')

genes <- as.integer(names(geneLengthData) %in% rownames(significantGenes))
names(genes) <- names(geneLengthData)
table(genes)

## run GoSeq
library('goseq')

# Fitting the Probability Weighting Function (PWF) to check if DE genes are
# biased by length

## Important: as we are supplying the Mapping in gene2cat, the species is not
## necessary (e.g. "NRRL3" has no function, but must be passed to work)

pwf=nullp(genes, "NRRL3", "ensGene", bias.data = geneLengthData)

GO = goseq(pwf, "NRRL3", "ensGene", gene2cat = goMapping)
GO$over_represented_padj <- p.adjust(GO$over_represented_pvalue, method="BH")
goResults <- subset(GO, over_represented_pvalue < 0.05)
rownames(goResults) <- goResults$category

GO$under_represented_padj <- p.adjust(GO$under_represented_pvalue, method="BH")
unEnrichedGo <- subset(GO, under_represented_pvalue<.05)
rownames(unEnrichedGo) <- unEnrichedGo$category

## write the enrichment results
write.table(as.data.frame(goResults), 
            file = "Mutation_kexB_vs_WT.enriched.goseq.GO.csv",
            row.names = FALSE,
            sep = ";")
write.table(as.data.frame(unEnrichedGo), 
            file = "Mutation_kexB_vs_WT.unenriched.goseq.GO.csv",
            row.names = FALSE,
            sep = ";")

## Gene Ontology Enrichment with Up & Down -------------------------------------

## Upregulated
## Get the significant genes
significantGenes <- as.data.frame(resMut) %>%
  rownames_to_column('gene') %>%
  filter(padj < 0.05) %>%
  filter(log2FoldChange < 0) %>% 
  column_to_rownames('gene')

genes <- as.integer(names(geneLengthData) %in% rownames(significantGenes))
names(genes) <- names(geneLengthData)
table(genes)

## run GoSeq
library('goseq')

pwf=nullp(genes, "NRRL3", "ensGene", bias.data = geneLengthData)

GO = goseq(pwf, "NRRL3", "ensGene", gene2cat = goMapping)
GO$over_represented_padj <- p.adjust(GO$over_represented_pvalue, method="BH")
goResults <- subset(GO, over_represented_pvalue < 0.05)
rownames(goResults) <- goResults$category

GO$under_represented_padj <- p.adjust(GO$under_represented_pvalue, method="BH")
unEnrichedGo <- subset(GO, under_represented_pvalue<.05)
rownames(unEnrichedGo) <- unEnrichedGo$category

## write the enrichment results
write.table(as.data.frame(goResults), 
            file = "Mutation_kexB_vs_WT.downregulated.enriched.goseq.GO.csv",
            row.names = FALSE,
            sep = ";")
write.table(as.data.frame(unEnrichedGo), 
            file = "Mutation_kexB_vs_WT.downregulated.unenriched.goseq.GO.csv",
            row.names = FALSE,
            sep = ";")


