list.of.packages <- c("DESeq2", "dplyr", "plyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(DESeq2)
library(dplyr)
library(plyr)
#library(tximport) used to see if DESeq2 will round counts from Tximport

#setwd("/home/coyote/Files_for_Oscar/Task09")
## there are 2 genes that are duplicates so "row.names = 1" fails
count.mat <- read.csv("mayo.path_aging.con.salmon.gene.counts.txt", sep = "\t", 
                        header = T)

## some gene names were converted to dates - probably got or store in excel
## map back using file from https://www.biostars.org/p/183018/
mangled.genes <- read.table("ExcelMangledGenes.txt", skip = 3, header = F)
count.mat$X <- mapvalues(count.mat$X, from = mangled.genes$V2, to = mangled.genes$V1)


## which genes are duplicated?
dups <- table(count.mat$X)[table(count.mat$X) > 1]
## MARCH1 and MARCH2 are in there twice so rename second of each MARCH*.1
## https://stat.ethz.ch/R-manual/R-patched/library/base/html/make.unique.html
rownames(count.mat) <- make.unique(count.mat$X)
count.mat <- count.mat[, -1]


pheno.dat <- read.csv("mayo.path_aging.con.phenotype.csv")

## sanity colnames = phenodata
all(colnames(count.mat) == pheno.dat$UID)
## col.names comes in with X* because first column header is blank
colnames(count.mat) <- sub("X", "", colnames(count.mat))
all(colnames(count.mat) == pheno.dat$UID)


## DESeq2 only takes integer counts this was just to confirm 
## https://support.bioconductor.org/p/133326/
## design for controlling for gender and age while testing Diagnosis condition
## comparisons will be based on the alphabetical order of the levels (Control is reference)
dseq.dat <- DESeqDataSetFromMatrix(round(count.mat), pheno.dat, design = ~ Diagnosis)


#########################################################
## Run DE analysis
#########################################################

## Control for Sex
## the controlling variable is the one preceded by "~", 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory
## measure the effect of the condition, controlling for batch differences = ~ batch + condition
## https://support.bioconductor.org/p/63766/ - controlling for Age & Sex comes before Phenotype in formula!
## I actually found out order does not matter. https://support.bioconductor.org/p/126713/
dseq.gender <- DESeqDataSetFromMatrix(round(count.mat),
                                      pheno.dat,
                                      design = ~ Sex + Diagnosis)
dds.gender <- DESeq(dseq.gender, betaPrior=FALSE)
res.gender <- results(dds.gender, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.gender.order.by.fc <- res.gender[order(abs(res.gender$log2FoldChange), decreasing = T), ]

p.gender <- plotCounts(dseq.gender, "POP1", intgroup = c("Diagnosis", "Sex"),
                       returnData = T, normalized = F, transform = F, pc = 1)
levels(p.gender$Sex) <- c("Female", "Male")
# ggplot(p.gender, aes(x=Diagnosis, y=log2(count), group = 1)) +
#   geom_point(position=position_jitter(w=0.1,h=0)) + 
#   stat_summary(fun = "mean", geom = "line", color = "red") +
#   facet_wrap(~ Sex) +
#   ylab("log2(counts + 1)") +
#   ggtitle("POP1") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text = element_text(size=11),
#         strip.text.x = element_text(size = 12),
#         #axis.title.x = element_text(size = 12),
#         axis.title.y = element_text(size = 12))


## Control for Age
pheno.dat2 <- pheno.dat
## make 'Age' continuous so remove "_or_above" from "90_or_above"
pheno.dat2$AgeAtDeath <- as.numeric(gsub("_.*", "", pheno.dat2$AgeAtDeath))
dseq.age <- DESeqDataSetFromMatrix(round(count.mat),
                                   pheno.dat2, 
                                   design = ~ AgeAtDeath + Diagnosis)
dds.age <- DESeq(dseq.age)
res.age <- results(dds.age, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.age.order.by.fc <-  res.age[order(abs(res.age$log2FoldChange), decreasing = T), ]


## Control for Sex and Age
dseq.gender.age <- DESeqDataSetFromMatrix(round(count.mat),
                                          pheno.dat2,
                                          design = ~ Sex + AgeAtDeath + Diagnosis)
dds.gender.age <- DESeq(dseq.gender.age)
resultsNames(dds.gender.age)
res.gender.age <- results(dds.gender.age, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.gender.age.order.by.fc <- res.gender.age[order(abs(res.gender.age$log2FoldChange), decreasing = T), ]


## write all results to files for faster reading for app
write.table(count.mat, "counts_cleaned.txt", sep = "\t", row.names = T)
write.table(res.gender.order.by.fc, "controlled_gender_res.txt", sep = "\t", row.names = T)
write.table(res.age.order.by.fc, "controlled_age_res.txt", sep = "\t", row.names = T)
write.table(res.gender.age.order.by.fc, "controlled_gender_age_res.txt", sep = "\t", row.names = T)


##############################################################################################
## tested below that DESeq2 rounded counts from salmon or really any data using Tximport
# library(readr)
# dir <- system.file("extdata", package = "tximportData") 
# samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
# files <- file.path(dir, "salmon", samples$run, "quant.sf.gz")
# tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))
# txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# 
# sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
# rownames(sampleTable) <- colnames(txi$counts)
# test <- DESeqDataSetFromTximport(txi, sampleTable, ~condition)
# all(head(counts(test, normalized=F)) == round(head(txi$counts)))

