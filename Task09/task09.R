library(DESeq2)
library(dplyr)
library(plyr)
#library(tximport) used to see if DESeq2 will round counts from 

setwd("/home/coyote/Files_for_Oscar/Task09")
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


## DESeq2 only takes integer counts
## https://support.bioconductor.org/p/133326/
## design for controlling for gender and age while testing Diagnosis condition
## comparisons will be based on the alphabetical order of the levels (Control is reference)
dseq.dat <- DESeqDataSetFromMatrix(round(count.mat), pheno.dat, design = ~ Diagnosis)
dds <- DESeq(dseq.dat)
res <- results(dds)
res

rownames(res)[which.min(res$padj)]
d1 <- plotCounts(dseq.dat, "SPCS2P4", intgroup = "Diagnosis",
                returnData = T, normalized = F)
res[which.min(res$padj), ]
groups <- grepl("Pathologic", pheno.dat$Diagnosis)
table(pheno.dat$Diagnosis)
mean(counts(dseq.dat)["SPCS2P4", groups])
mean(counts(dseq.dat)["SPCS2P4", -groups])
-1/(29.60976/53.78333)
library(MKmisc)
pairwise.fc(x=counts(dseq.dat)["SPCS2P4",], g=groups, log = F)

ggplot(d1, aes(x=Diagnosis, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

#####################################################################################
## Control for Sex
## the controlling variable is the one preceded by "~", 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory
## measure the effect of the condition, controlling for batch differences = ~ batch + condition
## https://support.bioconductor.org/p/63766/ - controlling for Age & Sex comes before Phenotype in formula!
## I actually found out order does not matter. https://support.bioconductor.org/p/126713/
dseq.sex <- DESeqDataSetFromMatrix(round(count.mat),
                                   pheno.dat2, 
                                   design = ~ Sex + Diagnosis)
dds.sex <- DESeq(dseq.sex)
res.sex <- results(dds.sex, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.sex.order.by.fc <-  res.sex[order(abs(res.sex$log2FoldChange), decreasing = T), ]

## Control for Age
pheno.dat2 <- pheno.dat
## make 'Age' continuous so remove "_or_above" from "90_or_above"
pheno.dat2$AgeAtDeath <- as.numeric(gsub("_.*", "", pheno.dat2$AgeAtDeath))
## https://support.bioconductor.org/p/126713/
dseq.age <- DESeqDataSetFromMatrix(round(count.mat),
                                   pheno.dat2, 
                                   design = ~ AgeAtDeath + Diagnosis)
dds.age <- DESeq(dseq.age)
res.age <- results(dds.age, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.age.order.by.fc <-  res.age[order(abs(res.age$log2FoldChange), decreasing = T), ]

p.age <- plotCounts(dseq.age, "AACS", intgroup = c("Diagnosis", "AgeAtDeath"),
                 returnData = T, normalized = F, transform = F, pc = 1)
ggplot(p.age) +
  ## ticks are not log2 - https://support.bioconductor.org/p/105938/
  geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
  stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
              formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
  ylab("log2(counts + 1)")


## Control for Sex and Age
dseq.gender.age <- DESeqDataSetFromMatrix(round(count.mat),
                                          pheno.dat2,
                                          design = ~ Sex + AgeAtDeath + Diagnosis)
dds.gender.age <- DESeq(dseq.gender.age)
resultsNames(dds.gender.age)
res.gender.age <- results(dds.gender.age, name = "Diagnosis_Pathologic.Aging_vs_Control")
res.gender.age.order.by.fc <- res.gender.age[order(abs(res.gender.age$log2FoldChange), decreasing = T), ]

p.gender.age <- plotCounts(dseq.gender.age, "AACS", intgroup = c("Diagnosis", "Sex", "AgeAtDeath"),
                           returnData = T, normalized = F, transform = F, pc = 1)
levels(p.gender.age$Sex) <- c("Female", "Male")

ggplot(p.gender.age) +
  ## ticks are not log2 - https://support.bioconductor.org/p/105938/
  geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
  stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
              formula = y ~ x + I(x^2), size = 1, alpha=0.25) +
  facet_wrap(~ Sex) +
  ylab("log2(counts + 1)") +
  ggtitle("AACS") +
  theme(plot.title = element_text(hjust = 0.5))



##############################################################################################





dds2 <- DESeq(dseq.dat2)
resultsNames(dds2)
res4 <- results(dds2, name = "DiagnosisControl.AgeAtDeath")
res5 <- results(dds2, name = "DiagnosisPathologic.Aging.AgeAtDeath")
res4.ord <- res4[order(abs(res4$log2FoldChange), decreasing = T), ]
res5.ord <- res5[order(abs(res5$log2FoldChange), decreasing = T), ]
res6 <- results(dds2, contrast = list("DiagnosisControl.AgeAtDeath", "DiagnosisPathologic.Aging.AgeAtDeath"))
res6.ord <- res6[order(abs(res6$log2FoldChange), decreasing = T), ]
res2 <- results(dds2)
res2
res2.ord <- res2[order(abs(res2$log2FoldChange), decreasing = T), ]
res2.ord$log2FoldChange
count.mat["A4GALT", ]
counts(dseq.dat2)["AACS", ]
## Normalized counts plus a pseudocount of 0.5 are shown by default
d2 <- plotCounts(dseq.dat2, "AACS", intgroup = c("Diagnosis", "AgeAtDeath"),
                returnData = T, normalized = F, transform = F, pc = 1)
mean(counts(dseq.dat3)["DAO", groups])
mean(counts(dseq.dat)["DAO", -groups])
res3.ord <- res3[order(abs(res3$log2FoldChange), decreasing = T), ]
ggplot(d2) +
  ## https://support.bioconductor.org/p/105938/
  geom_point(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis)) +
  stat_smooth(aes(x=AgeAtDeath, y=log2(count), color=Diagnosis), method = "lm",
              formula = y ~ x + I(x^2), size = 1, alpha=0.25)
table(pheno.dat2[,c("Diagnosis", "AgeAtDeath")])
#####################################################################################
## Control for Sex
## the controlling variable is the one preceded by "~", 
## https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory
## measure the effect of the condition, controlling for batch differences = ~ batch + condition
## https://support.bioconductor.org/p/63766/ - controlling for Age & Aex comes before Phenotype in formula!
## I actually found out order does not matter.
dseq.dat3 <- DESeqDataSetFromMatrix(round(count.mat), 
                                    pheno.dat2, design = ~ Diagnosis + Sex + Diagnosis:Sex)
dseq.dat4 <- DESeqDataSetFromMatrix(round(count.mat), 
                                    pheno.dat2, design = ~ Diagnosis + Diagnosis:Sex)
dseq.dat5 <- DESeqDataSetFromMatrix(round(count.mat), 
                                    pheno.dat2, design = ~ Sex + Diagnosis + Sex:Diagnosis)
dseq.dat6 <- DESeqDataSetFromMatrix(round(count.mat), 
                                    pheno.dat2, design = ~ Sex + Diagnosis)
dseq.dat7 <- DESeqDataSetFromMatrix(round(count.mat), 
                                    pheno.dat2, design = ~ Diagnosis + Sex)

dds11 <- DESeq(dseq.dat6)
dds12 <- DESeq(dseq.dat7)
resultsNames(dds11)
resultsNames(dds12)
results(dds11, name="Sex_M_vs_F")
results(dds12, name="Sex_M_vs_F")
results(dds11, name="Diagnosis_Pathologic.Aging_vs_Control")
results(dds12, name="Diagnosis_Pathologic.Aging_vs_Control")
#dds3 <- DESeq(dseq.dat3)
#dds4 <- DESeq(dseq.dat4)
dds5 <- DESeq(dseq.dat5)
res9 <- results(dds5, name="SexM.DiagnosisPathologic.Aging")
results(dds5, name="Diagnosis_Pathologic.Aging_vs_Control")
results(dds3, name = "Diagnosis_Pathologic.Aging_vs_Control")
res10 <- results(dds3, name = "DiagnosisPathologic.Aging.SexM")
round(res9$log2FoldChange,2) == round(res10$log2FoldChange,2)
resultsNames(dds3)
resultsNames(dds4)
resultsNames(dds5)
res3 <- results(dds3)


res7 <- results(dds3, contrast = list("Sex_M_vs_F", 
                                      "DiagnosisPathologic.Aging.SexM"))
results(dds3, contrast = c("Sex", "M", "F"))
results(dds3, name = "Sex_M_vs_F")

res8 <- results(dds4, contrast = list("DiagnosisControl.SexM", 
                                      "DiagnosisPathologic.Aging.SexM"))
results(dds4, contrast = list("DiagnosisPathologic.Aging.SexM", 
                              "DiagnosisControl.SexM"))
res7.ord <- res7[order(abs(res7$log2FoldChange), decreasing = T), ]
res8.ord <- res8[order(abs(res8$log2FoldChange), decreasing = T), ]
library(ggplot2)
res3[which.max(abs(res3$log2FoldChange)), ]
d <- plotCounts(dseq.dat3, "XAGE1E", intgroup = c("Diagnosis", "Sex"),
                returnData = T)
mean(counts(dseq.dat3)["DAO", groups])
mean(counts(dseq.dat)["DAO", -groups])
res3.ord <- res3[order(abs(res3$log2FoldChange), decreasing = T), ]
ggplot(d, aes(x=Diagnosis, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400)) + 
  facet_wrap(~ Sex)

p.gender <- plotCounts(dseq.sex, "DAO", intgroup = c("Diagnosis", "Sex"),
                       returnData = T, normalized = F, transform = F, pc = 1)
levels(p.gender$Sex) <- c("Female", "Male")
ggplot(p.gender, aes(x=Diagnosis, y=log2(count), group = 1)) +
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  stat_summary(fun = "mean", geom = "line", color = "red") +
  facet_wrap(~ Sex) +
  ylab("log2(counts + 1)") +
  ggtitle("DAO") +
  theme(plot.title = element_text(hjust = 0.5))

## tested below that DESeq2 rounded counts from salmon or any data
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

