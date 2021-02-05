
getwd()
setwd("C:/Users/gottams/Desktop/Mexican_Allfiles")
setwd("C:/Users/gottams/Desktop")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

chooseCRANmirror()

options(repos = c(CRAN = "http://cran.rstudio.com"))

options(repos = c(CRAN = "http://cran.rstudio.com"))
install.packages("dplyr")

install.packages("ggplot2", repos = "http://cran.rstudio.com")

BiocManager::install("qqman")

library("qqman")

BiocManager::install("GWASTools")
library("GWASTools")

BiocManager::install("ggplot2")
library("ggplot2")

BiocManager::install("RColorBrewer")

BiocManager::install("scatterplot3d")

BiocManager::install("genFun")

BiocManager::install("sqldf")
library(sqldf)

BiocManager::install("survcomp")
library("survcomp")

BiocManager::install("metap")
library("metap")

library(dplyr)
library(tidyr)

memory.limit()

memory.limit(10000000000000)

BiocManager::install("QCEWAS")
library("QCEWAS")

BiocManager::install("forestplot")

BiocManager::install("lattice")

BiocManager::install("manhattanly")

BiocManager::install("rsnps")
library("rsnps")

BiocManager::install("CMplot")
library("CMplot")

BiocManager::install("plotrix")
library("plotrix")

install.packages("remotes")
remotes::install_github("drveera/ggman")

install.packages("devtools")
library(devtools)

devtools::install_github("drveera/ggman")

install.packages("gridgraphics")


###########################################################################

MainInfile <- read.delim("SUM.txt",header = TRUE, sep = "\t")

snp346 <- read.delim("346snplist.txt", header = TRUE)

sapply(snp346, class)

x_fac <- as.factor(snp346$X100k...)

y_fac <- as.factor(snp346$X100k....1)

x_fac <- gsub(",","",x_fac)

y_fac <- gsub(",","",y_fac)

x_fac_to_num <- as.numeric(as.character(x_fac))

y_fac_to_num <- as.numeric(as.character(y_fac))

x <- sqldf("select * from MainInfile where BP = x_fac_to_num & BP = y_fac_to_num")



###############################################################
AA_GWAS <- read.table("aagwas pc1-3.plink",header = T)

Mex1 <- read.table("chrX2.plink",header = T)

filter77 <- Mex[Mex$BP  == 153275890,]
filter78 <- Mex[Mex$BP  == 153278829,]
filter79 <- Mex[Mex$BP  == 153284192,]


filter80 <- Mex1[Mex1$BP  == 153275890,]
filter81 <- Mex1[Mex1$BP  == 153278829,]
filter82 <- Mex1[Mex1$BP  == 153284192,]

y1 <- write.csv(filter77,"filter77.csv", row.names = TRUE)

y2 <- write.csv(filter78,"filter78.csv", row.names = TRUE)

y3 <- write.csv(filter79,"filter79.csv", row.names = TRUE)

y4 <- write.csv(filter80,"filter80.csv", row.names = TRUE)

y5 <- write.csv(filter81,"filter81.csv", row.names = TRUE)

y6 <- write.csv(filter82,"filter82.csv", row.names = TRUE)

GWAS <- read.table("filter78.plink",header = T)

GWAS2 <- read.delim("filter79.plink", header = T)

GWAS40 <- read.delim("PASSLUPUS.tsv",header = T)

head(GWAS40)

GWAS44 <- GWAS40[,c(1,3,6)]

head(GWAS44)

GWAS48 <- GWAS44[GWAS44$chrom >= 1 & GWAS44$chrom <= 22,]



GWAS46 <- write.table(GWAS44,"file5.txt")

head(GWAS2)

GWAS3 <- GWAS2[GWAS2$CHR == 12 & GWAS2$BP == 12870695,]

head(GWAS)

GWAS1 <- GWAS[GWAS$CHR == 12 & GWAS$BP == 12870695,]

filter555 <- write.csv(GWAS1,"GWAS1.csv", row.names = TRUE)

GWAS4 <- read.table("chrX12_Rsq30.assoc",header = T)




GWAS77 <- read.table("Spain_Results.txt.gz",header = T)

#############################################################################

Agwashrs <- read.table("AAgwasHRSAdditional.plink",header = T)

GSA5plate <- read.table("GSA_5_plates.plink",header = TRUE )

GSA12plate <- read.table("MexGSAHG19.assoc",header = TRUE)

GSA13plate <- read.table("GSA_5_plates.assoc",header = TRUE)

GSA14plate <- read.table("GSA_5plates_Mexican_geno10_hwe0001_mind20_maf005_assoc.assoc",header = TRUE)

GSA20plate <- read.table("Celi_MexicanGSA5plate.assoc",header = TRUE)

GSA21plate <- read.table("GSA5Plate_noOBIqc.assoc",header = TRUE)

GSA22plate <- read.table("GSA5Platemain.assoc",header = TRUE)

GSA23plate <- read.table("GSAnoobi_1.assoc",header = TRUE)

GSplate <- read.table("GS.assoc",header = TRUE)

GS1plate <- read.table("GS1.assoc",header = TRUE)

GSA5platenoobi <- read.table("GSA5plate.plink",header = TRUE)

chrX <- read.table("chrX5.txt",header = TRUE)

chrX1 <- read.table("chrX2.plink",header = TRUE)

MexGSA12plate <- read.table("MainFile_filtered.assoc",header = TRUE )

MexGSA12plate_Mex <- read.table("MainFile_Mexican.assoc",header = TRUE)

MexGSA <- read.table("MainFile_Mexican.assoc",header = TRUE)

YukGSA <- read.table("Yuk.assoc",header = TRUE)

PSSGSA <- read.table("PSS_1.assoc",header = TRUE)

MexGSA1 <- read.table("New_2.assoc",header = TRUE)

Mexx <- read.table("GSA12plateMexAutosome.assoc",header = TRUE)

filter991 <- Mexx[Mexx$F_A != 0 & Mexx$F_U != 0,]

head(filter991)

YukkGSA <- read.table("Yuk_filtered_1.assoc",header = TRUE)

MexGSA_5p <- read.table("GSA5plate_noobi_qc.assoc",header = TRUE)

qq(MexGSA_5p$P)

PSSGSA <- read.table("New_PSS_1.assoc",header = TRUE)

Mex <- read.table("NewMexfil.assoc",header = TRUE)

Yuk <- read.table("NewYukfil.assoc",header = TRUE)

PSS <- read.table("NewPSSfil.assoc",header = TRUE)

Mex_11 <- read.table("Mex_001.txt",header = TRUE)

NewMex <- read.table("GSA12plate_auto.assoc",header = TRUE)

NewMexx <- read.table("GSA12plate_Mexican_qc1.assoc",header = TRUE)

NewYukk <- read.table("GSA12plate_Yukatan_qc1.assoc",header = TRUE)

NewPSSS <- read.table("GSA12plate_PSS_qc1.assoc",header = TRUE)

NewMain <- read.table("MainFile.txt",header = TRUE)



qq(Mex_11$P)

Yuk_11 <- read.table("Yuk_001.txt",header = TRUE)

MetaSE <- read.table("InversebasedMexicanMeta.TBL",header = TRUE)

MetaSample <- read.table("SamplebasedMexicanMeta.TBL",header = TRUE)

Finalreport <- read.table("finalreport.txt",header = TRUE)

final <- read.table("final2.txt",header = TRUE, sep = "/t")

Newfile <- read.table("Newfile_filtered.assoc",header = TRUE)

NewMexFile <- read.table("Newfile_filtered_Mexican.assoc",header = TRUE)

NewYukFile <- read.table("Newfile_filtered_Yukatan.assoc",header = TRUE)

NewPSSFile <- read.table("Newfile_filtered_PSS.assoc",header = TRUE)

Newobi <- read.table("obi_filtered.assoc",header = TRUE)

NewGSA <- read.table("GSAA.assoc",header = TRUE)

NewGSAHG19 <- read.table("NewGSA.assoc",header = TRUE)

NewGSAA <- read.table("GSA_NoYukatan.assoc",header = TRUE)

NewGSS <- read.table("NewMexNoYuk.txt",header = TRUE)

NewMexx <- read.table("GSA_NoYukatan_001_Mex.assoc",header = TRUE)

NewPS <- read.table("GSA_NoYukatan_001_PSS.assoc",header = TRUE)

New <- read.table("GSAMerge.assoc",header = TRUE)

NewHG <- read.table("GSA_5_12_merge_HG19_controls_a.assoc",header = TRUE)

NewOBI <- read.table("GSAObi.assoc",header = TRUE)

New1 <- read.table("Controls1.assoc",header = TRUE)

New2 <- read.table("NewGSA12.assoc",header = TRUE)

New4 <- read.table("GSA5PlateNew.assoc",header = TRUE)

plink <- read.table("plink.assoc",header = TRUE)

plink1 <- read.table("NewGSA12platefinal.assoc",header = TRUE)

plink2 <- read.table("Twoplates.assoc",header = TRUE)

plink3 <- read.table("Twoplates_Controls_1.assoc",header = TRUE)

plink4 <- read.table("GSA12plateMexUnfiltered.assoc",header = TRUE)

plink5 <- read.table("filteredGSAPC1-3.txt",header = TRUE)

plink6 <- read.table("UnfilteredGSAPC1-3.txt",header = TRUE)

p5 <- read.table("plink5.txt",header = TRUE)

p6 <- read.table("NewGSA12plate_Mex.assoc",header = TRUE)

p7 <- read.table("TwoplateControls_MexNew.assoc",header = TRUE)

p8 <- read.table("GSA5.assoc",header = TRUE)

p9 <- read.table("GSA12.assoc",header = TRUE)

p10 <- read.table("Samplebased.txt",header = TRUE)

p11 <- read.table("Inversebased.txt",header = TRUE)

p12 <- read.table("ColumbiaAffy.txt",header = TRUE)

p13 <- read.table("3setsMetaAnalysis.TBL",header = TRUE)

p14 <- read.table("chrx_Rsq30.assoc",header = TRUE)

p15 <- read.table("plink.assoc",header = TRUE)

p16 <- read.table("SamplebasedChrX.txt",header = TRUE)

p17 <- read.table("chrX2.plink",header = TRUE)

p18 <- read.table("GSA12plateChrX_filtered.assoc",header = TRUE)

p20 <- read.table("Final.txt",header = TRUE)

p22 <- read.table("AllChrSet2.txt",header = TRUE)

p23 <- read.table("SampleMeta.TBL",header = TRUE)

p24 <- read.table("SampleMetaNew.TBL",header = TRUE)

p25 <- read.table("InverseMetaNew.TBL",header = TRUE)

p26 <- read.table("SampleMeta.txt",header = TRUE)

p31 <- read.table("InverseMeta.txt",header = TRUE)

p33 <- read.table("SampleNew.txt",header = TRUE)

p34 <- read.table("InverseNew.txt",header = TRUE)

p35 <- read.table("Set1New.txt",header = TRUE)

p36 <- read.table("Set2New.txt",header = TRUE)

p27 <- read.table("Set1Mexican.txt",header = TRUE)

p28 <- read.table("Set2Mexican.txt",header = TRUE)

p36 <- read.table("GSA5plate.plink",header = TRUE)

Set2chrX <- read.table("GSASet2chrXMex_hg19New.assoc",header = TRUE)

Set2chrXunfilter <- read.table("GSASet2chrXMex_hg19.assoc",header = TRUE)


#############################################################################################

#Extracting the snp name and Chr number for the snps in the meta analysis result of Columbian and Set1 Mexican data from the original data

p44 <- read.table("Set1andColumbian.TBL",header = TRUE)

p48 <- read.table("Set1Impute.txt",header = TRUE)

p49 <- read.table("ColumImputed.txt",header = TRUE)
p48_49 <- rbind.data.frame(p48[,c("CHR", "BP", "SNP")], p49[,c("CHR", "BP", "SNP")])
p48_49$SNP <- gsub("GSA-", "", p48_49$SNP)

p44 <- cbind.data.frame(p44, p48_49[match(p44$MarkerName, p48_49$SNP),])

######################################################################################################################

#Extracting the snp name and Chr number for the snps in the meta analysis result of Columbian Imputed data and Set1 Imputed data

p74 <- read.table("Samplebased_.TBL",header = TRUE)

p75 <- read.table("GSA5PlateImputed.txt",header = TRUE)

p76 <- read.table("ColumbianSLEImputed.txt",header = TRUE)

head(p74)

p62_63 <- rbind.data.frame(p62[,c("CHR", "BP", "SNP")], p63[,c("CHR", "BP", "SNP")], p63[,c("CHR", "BP", "SNP")] )

p62_63_1 <- rbind.data.frame(p62[,c("SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")], p63[,c("SNP", "A1", "F_A" , "F_U", "A2", "info", "P", "OR" )])

rownames(p62) <- p62$SNP
rownames(p63) <- p63$SNP

p62_1 <- p62[filter211$MarkerName, c("SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p63_1 <- p63[filter211$MarkerName, ,c("SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]





head(p74)

qqman::manhattan(p74, chr="CHR", bp="BP", snp="MarkerName", p="P.value", ylim =c(0,20))


####################################################################################

#Finding the common snps between Set1 Mexican data and Columbian Data

p60 <- length(intersect(p48$SNP, p49$SNP))

p99 <- length(intersect(filter0081$SNP, filter0082$SNP))

p103 <- length(union(filter0081$SNP, filter0082$SNP))

p100 <- length(intersect(filter0082$SNP, filter0083$SNP))

p101 <- length(intersect(filter0081$SNP, filter0083$SNP))

p102 <- length(union(union(filter0081$SNP,filter0082$SNP),filter0083$SNP))

#####################################################################################

#Preparing the files for the Meta analysis of Set1 Mexican Data and Columbian Imputed Data


p61 <- read.table("SLEcolAffy.plink",header = TRUE)


head(p61)

filter990 <- p61[p61$F_A > 0.01 & p61$F_U > 0.01 & p61$info > 0.70,]

SLEMetaColumbia <- p61[,c("SNP","A1", "A2", "OR", "F_A", "P", "SE")]
SLEMetaColumbia$OR <- log(SLEMetaColumbia$OR)
colnames(SLEMetaColumbia)[4] <- "EFFECT"

p70 <- p63[,c("SNP","A1", "A2", "OR", "F_A", "P", "SE")]
p70$OR <- log(p70$OR)
colnames(p70)[4] <- "EFFECT"
p70$N <- 1682
p70 <- p70[,c(1:3,8,4:7)]

p71 <- p62[,c("SNP","A1", "A2", "OR", "F_A", "P", "SE")]
p71$OR <- log(p71$OR)
colnames(p71)[4] <- "EFFECT"
p71$N <- 465
p71 <- p71[,c(1:3,8,4:7)]


p93 <- p92[,c("SNP", "A1", "A2", "OR", "F_A", "P", "SE")]

head(p90)


p83 <- p81[p81$F_A > 0.01 & p81$F_U > 0.01 & p81$info > 0.70,]


head(p70)

head(p71)


head(SLEMetaColumbia)

head(p61)


head(p36)


p62 <- p36[p36$F_A > 0.01 & p36$F_U > 0.01 & p36$info > 0.70,]

p63 <- p61[p61$F_A > 0.01 & p61$F_U > 0.01 & p61$info > 0.70,]

head(p62)

head(p63)

head(p63_1)

#################################################################################################

p008 <- read.table("GSA5plate.plink", header = TRUE)

GWAS11 <- filter0081[filter0081$CHR == 16 & filter0081$BP == 31276811,]

GWAS12 <- filter0082[filter0082$CHR == 16 & filter0082$BP == 31276811,]

GWAS13 <- filter0083[filter0083$CHR == 16 & filter0083$BP == 31276811,]

GWAS14 <- filter0081[filter0081$CHR == 7 & filter0081$BP == 128579666,]

GWAS15 <- filter0082[filter0082$CHR == 7 & filter0082$BP == 128579666,]

GWAS16 <- filter0083[filter0083$CHR == 7 & filter0083$BP == 128579666,]


p80 <- read.table("GSA12plateclean.plink", header = TRUE)

GWAS17 <- p80[p80$CHR == 16 & p80$BP == 31276811,]

p88 <- read.table("Set1Set2Common.TBL",header = TRUE)

Set13 <- read.table("Set1Set3.TBL",header = TRUE)

Set23 <- read.table("Set2Set3.TBL",header = TRUE)

p1_1 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP")], filter0082[,c("CHR", "BP", "SNP")])

p88 <- cbind.data.frame(p88, p1_1[match(p88$MarkerName, p1_1$SNP),])

head(p88)

head(p88)

p100_1 <- filter0081[filter0082$SNP, c("CHR" ,"BP ", "SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p011 <- rbind.data.frame(filter0081[,c("CHR","SNP", "BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")])


p02_03 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP")], filter0082[,c("CHR", "BP", "SNP")], filter0083[,c("CHR", "BP", "SNP")]) 

head(Set13)

head(Set23)

p02_2 <- cbind.data.frame(p02_2, p02_03[match(p02_2$MarkerName, p02_03$SNP),])

p02_5 <- rbind.data.frame(p88[,c("SNP","P.value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")])

p02_5 <- rbind.data.frame(p88[,c("SNP","P.value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")])

p02_5 <- rbind.data.frame(p88[,c("SNP","P.value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")])

p01_8 <- rbind.data.frame(Set13[,c("SNP","P.value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")])

p01_9 <- rbind.data.frame(Set23[,c("SNP","P.value", "Direction","HetISq","HetChiSq","HetDf","HetPVal")])

p02_6 <- cbind.data.frame(filter980, p02_5[match(filter980$MarkerName, p02_5$SNP),])

head(p02_6)

p02_11 <- cbind.data.frame(filter980, p01_8[match(filter980$MarkerName, p01_8$SNP),])


p02_12 <- cbind.data.frame(filter980, p01_9[match(filter980$MarkerName, p01_9$SNP),])

head(p02_2)



head(filter0081)

merge(filter0081,filter0082, by="SNP")

p81 <- read.table("GSA12plateclean.plink",header = TRUE )

p_cases <- read.table("TwoplateCases1.assoc",header = TRUE)



p_casesuniq <- read.table("TwoplateCasesUniq_1.assoc",header = TRUE)

p_controlsuniq <- read.table("TwoPlateControls_1.assoc",header = TRUE)

p_controlsuniq1 <- read.table("ControlsUniq.txt",header = TRUE)

p_controlsuniq <- na.omit(p_controlsuniq)

p_casesuniq <- na.omit(p_casesuniq)

p_controls <- read.table("TwoplatesControl1.assoc",header = TRUE)

TwoPlateMerge <- read.table("TwoPlateMerge_QC.assoc",header = TRUE)

filter2119 <- TwoPlateMerge[TwoPlateMerge$P < 0.00001 & TwoPlateMerge$P <=1,]

TwoPlateMerge_noATGC <- read.table("TwoPlateMerge_QC_noATGC.assoc",header = TRUE)

TwoPlateControl <- read.table("TwoPlateControlsQC.assoc",header = TRUE)

TwoPlateCase <- read.table("TwoplateCasesQC.assoc",header = TRUE)

Set1 <- read.table("GSA5PlateAutsomeNoOBI.assoc", header = TRUE)

Set2 <- read.table("GSA12plate_hg19_Mexautosome.assoc", header = TRUE)

gen <- read.table("1000GMXLMerged.txt",header = TRUE)

Columbian <- read.table("ColumbianSLE.txt",header = TRUE)

ColumbianMeta <- read.table("ColumbianMeta.txt",header = TRUE)

pc1_20 <- read.table("TwoPlateMerge_QC_PC1-20_1.assoc.logistic",header = TRUE) 

pc1_3 <- read.table("TwoPlateMerge_QC_PC1-3_1.assoc.logistic",header = TRUE)

pc1_5 <- read.table("TwoPlateMerge_QC_PC1-5_1.assoc.logistic",header = TRUE)

pc1_10 <- read.table("TwoPlateMerge_QC_PC1-10_1.assoc.logistic",header = TRUE)

plate5chrX <- read.table("GSA5platechrX.txt",header = TRUE)

plate12chrX <- read.table("GSA12platechrX.txt",header = TRUE)

platemerge <- read.table("chrXMerge2plate.txt",header = TRUE)

plate5 <- read.table("5platebeforeimp.assoc",header = TRUE)

plate12 <- read.table("12platebeforeimp.assoc",header = TRUE)

platemerge1 <- read.table("Mergeplatebeforeimp.assoc",header = TRUE)

chrX5plate <- read.table("chrXGSA5plate.plink",header = TRUE)

chrX12plate <- read.table("chrXSet2.plink",header = TRUE)

chrXSampleMeta <- read.table("Set1Set2chrXmaf0005Sample.TBL",header = TRUE)

chrXInverseMeta <- read.table("Set1Set2chrXmaf0005Inverse.TBL",header = TRUE)

chrXMetaSample <- read.table("chrXSampleMeta.txt",header = TRUE)

chrXInverseSample <- read.table("chrXInverseMeta.txt",header = TRUE)

chrX1000g <- read.table("1000gMXLMergedchrX2.assoc",header = TRUE)

head(chrXSampleMeta)

filter01 <- chrX12plate[chrX12plate$F_A >= 0.005 & chrX12plate$F_U >= 0.005 & chrX12plate$F_A <= 1 & chrX12plate$F_U <= 1 & chrX12plate$info > 0.30, ]

filter02 <- chrX5plate[chrX5plate$F_A >= 0.005 & chrX5plate$F_U >= 0.005 & chrX5plate$F_A <= 1 & chrX5plate$F_U <= 1 & chrX5plate$info > 0.30, ]

filter03 <- chrX12plate[chrX12plate$F_A >= 0.005 & chrX12plate$F_U >= 0.005 & chrX12plate$F_A <= 1 & chrX12plate$F_U <= 1 & chrX12plate$info > 0.70, ]

filter04 <- chrX5plate[chrX5plate$F_A >= 0.005 & chrX5plate$F_U >= 0.005 & chrX5plate$F_A <= 1 & chrX5plate$F_U <= 1 & chrX5plate$info > 0.70, ]

head(chrX5plate)

pc1_20$cordinate <- paste(pc1_20$CHR, pc1_20$BP, sep = ":")

pc1_3$cordinate <- paste(pc1_3$CHR, pc1_3$BP, sep = ":")

pc1_5$cordinate <- paste(pc1_5$CHR, pc1_5$BP, sep = ":")

pc1_10$cordinate <- paste(pc1_10$CHR, pc1_10$BP, sep = ":")

head(pc1_10)

p888 <- TwoPlateMerge_noATGC[match(TwoPlateMerge$CHR.BP,TwoPlateMerge_noATGC$CHR.BP ),]


head(pc1_20)

p008_3 <- rbind.data.frame(TwoPlateMerge_noATGC[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_4 <- rbind.data.frame(TwoPlateControl[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_5 <- rbind.data.frame(TwoPlateCase[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_6 <- rbind.data.frame(Set1[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_7 <- rbind.data.frame(Set2[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_8 <- rbind.data.frame(gen[,c("CHR.BP","SNP", "A1","A2","F_U")])

p008_9 <- rbind.data.frame(Columbian[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_10 <- rbind.data.frame(p61[,c("cordinate","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_11 <- rbind.data.frame(ColumbianMeta[,c("CHR.BP","SNP", "A1", "F_A", "F_U", "A2", "P", "OR")])

p008_12 <- rbind.data.frame(pc1_3[,c("cordinate","SNP", "A1", "P", "OR")])

p008_13 <- rbind.data.frame(pc1_5[,c("cordinate","SNP", "A1", "P", "OR")])

p008_14 <- rbind.data.frame(pc1_10[,c("cordinate","SNP", "A1", "P", "OR")])

p008_15 <- rbind.data.frame(pc1_20[,c("cordinate","SNP", "A1", "P", "OR")])

p008_16 <- rbind.data.frame(plate5chrX[,c("CHR.BP","A1", "A2", "F_A", "F_U","P","OR","Rsq")])

p008_17 <- rbind.data.frame(plate12chrX[,c("CHR.BP","A1", "A2", "F_A", "F_U","P","OR","Rsq")])

p008_18 <- rbind.data.frame(chrX5plate[,c("CHR", "BP", "SNP")], chrX12plate[,c("CHR", "BP", "SNP")])

p008_19 <- rbind.data.frame(chrX5plate[,c("CHR","BP","SNP","A1","A2","F_A","F_U","P","OR")])

p008_20 <- rbind.data.frame(chrX12plate[,c("CHR","BP","SNP","A1","A2","F_A","F_U","P","OR")])

p008_21 <- rbind.data.frame(chrXInverseSample[,c("CHR","BP","SNP","Allele1","Allele2","P.value","Direction","Effect")])

p008_22 <- rbind.data.frame(chrX1000g[,c("CHR.BP","A1","A2","F_U")])

head(chrXInverseSample)

p97_1 <- cbind.data.frame(TwoPlateMerge, p008_3[match(TwoPlateMerge$CHR.BP, p008_3$CHR.BP),])

p97_2 <- cbind.data.frame(TwoPlateMerge, p008_4[match(TwoPlateMerge$CHR.BP, p008_4$CHR.BP),])

p97_3 <- cbind.data.frame(TwoPlateMerge, p008_5[match(TwoPlateMerge$CHR.BP, p008_5$CHR.BP),])

p97_4 <- cbind.data.frame(TwoPlateMerge, p008_6[match(TwoPlateMerge$CHR.BP, p008_6$CHR.BP),])

p97_5 <- cbind.data.frame(TwoPlateMerge, p008_7[match(TwoPlateMerge$CHR.BP, p008_7$CHR.BP),])

p97_6 <- cbind.data.frame(TwoPlateMerge, p008_8[match(TwoPlateMerge$CHR.BP, p008_8$CHR.BP),])

p97_7 <- cbind.data.frame(TwoPlateMerge, p008_9[match(TwoPlateMerge$CHR.BP, p008_9$CHR.BP),])

p97_8 <- cbind.data.frame(TwoPlateMerge, p008_10[match(TwoPlateMerge$CHR.BP, p008_10$cordinate),])

p97_9 <- cbind.data.frame(TwoPlateMerge, p008_11[match(TwoPlateMerge$CHR.BP, p008_11$CHR.BP),])

p97_10 <- cbind.data.frame(TwoPlateMerge, p008_12[match(TwoPlateMerge$CHR.BP, p008_12$cordinate),])

p97_11 <- cbind.data.frame(TwoPlateMerge, p008_13[match(TwoPlateMerge$CHR.BP, p008_13$cordinate),])

p97_12 <- cbind.data.frame(TwoPlateMerge, p008_14[match(TwoPlateMerge$CHR.BP, p008_14$cordinate),])

p97_13 <- cbind.data.frame(TwoPlateMerge, p008_15[match(TwoPlateMerge$CHR.BP, p008_15$cordinate),])

p97_14 <- cbind.data.frame(platemerge, p008_16[match(platemerge$CHR.BP, p008_16$CHR.BP),])

p97_15 <- cbind.data.frame(platemerge, p008_17[match(platemerge$CHR.BP, p008_17$CHR.BP),])

p97_16 <- cbind.data.frame(chrXSampleMeta, p008_18[match(chrXSampleMeta$MarkerName, p008_18$SNP),])

p97_17 <- cbind.data.frame(chrXInverseMeta, p008_18[match(chrXInverseMeta$MarkerName, p008_18$SNP),])

p97_18 <- cbind.data.frame(chrXMetaSample, p008_19[match(chrXMetaSample$SNP, p008_19$SNP),])

p97_19 <- cbind.data.frame(chrXMetaSample, p008_20[match(chrXMetaSample$SNP, p008_20$SNP),])

p97_20 <- cbind.data.frame(chrXMetaSample, p008_21[match(chrXMetaSample$SNP, p008_21$SNP),])

p97_21 <- cbind.data.frame(chrXMetaSample, p008_22[match(chrXMetaSample$CHR.BP, p008_22$CHR.BP),])

head(p97_8)

head(chrXMetaSample)

filt318 <- write.table(p97_1,"MainFile_NoATGCComp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt319 <- write.table(p97_2,"MainFile_ControlsComp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt320 <- write.table(p97_3,"MainFile_CasesComp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt321 <- write.table(filter2119,"SigSNPSMainFile.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt322 <- write.table(p97_4,"MainFile_Set1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt323 <- write.table(p97_5,"MainFile_Set2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt324 <- write.table(p97_6,"MainFile_1000G.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt325 <- write.table(p97_7,"MainFile_Columbian.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt327 <- write.table(p97_9,"MainFile_ColumbianMeta.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt328 <- write.table(p97_10,"MainFile_PC1-3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt329 <- write.table(p97_11,"MainFile_PC1-5.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt330 <- write.table(p97_12,"MainFile_PC1-10.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt331 <- write.table(p97_13,"MainFile_PC1-20.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt332 <- write.table(p97_14,"MainFILE_Plate5.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt333 <- write.table(p97_15,"MainFile_Plate12.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt326 <- write.table(p97_8,"MainFile_ColumbianImputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt336 <- write.table(p97_18,"MetaSample_Set1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt337 <- write.table(p97_19,"MetaSample_Set2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt338 <- write.table(p97_20,"MetaSample_Inverse.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt334 <- write.table(p97_16,"chrXSampleMeta.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt335 <- write.table(p97_21,"chrX1000G_MetaSample.txt", row.names = FALSE, quote = FALSE, sep = "\t")




head(p61)

p87 <- p61[p61$P < 0.00001 & p61$P <= 1,]

p89 <- na.omit(p87)


p1_1 <- length(intersect(p_cases$SNP, p_controls$SNP))

p_controls_1 <- na.omit(p_controls)



p_casecontrolmerged <- read.table("Twoplatesmerge.assoc",header = TRUE)

p_1 <- na.omit(p_casecontrolmerged)



filter2111 <- p_cases[p_cases$P < 0.05 & p_cases$P <=1,]

filter2112 <- p_cases[p_cases$P < 0.001 & p_cases$P <=1,]

filter2113 <- p_cases[p_cases$F_A > 0.01 & p_cases$F_U > 0.01,]

filter2114 <- p_controls[p_controls$F_A > 0.01 & p_controls$F_U > 0.01,]

filtercontrols <- p_controlsuniq1[p_controlsuniq1$F_A > 0.01 & p_controlsuniq1$F_U > 0.01,]

filtercases <- p_casesuniq[p_casesuniq$F_A > 0.01 & p_casesuniq$F_U > 0.01,]

filter2115 <- p_cases[p_cases$F_A == 0,]

filter2115_1 <- na.omit(filter2115)

filter2116 <- p_cases[p_cases$F_U == 0,]

filter2117 <- p_controls[p_controls$F_A == 0,]

filter2117_1 <- na.omit(filter2117)

filter2118 <- p_controls[p_controls$F_U == 0,]

filter2118_1 <- na.omit(filter2118)

head(filter2117)

filter2116 <- p_controls[p_controls$F_A > 0.01 & p_controls$F_U > 0.01,]




head(p81)


df %>% separate(p81, 
                c("CHR1", "BP1","A11", "A12"))



filter991 <- p81[ p81$CHR == 16 & p81$BP == 31276811,]

filter995 <- p81[ p81$CHR == 12 & p81$BP == 12870695,]

filter996 <- p008[ p008$CHR == 12 & p008$BP == 12870695,]

filter997 <- p61[ p61$CHR == 12 & p61$BP == 12870695,]

filter0081 <- p008[p008$F_A > 0.01 & p008$F_U > 0.01 & p008$info > 0.70,]

filter997 <- filter0081[filter0081$CHR == 3 & filter0081$BP == 58183636,]

filter998 <- filter0082[filter0082$CHR == 3 & filter0082$BP == 58183636,]

fil1 <- filter0081[filter0081$CHR == 22,]

filter0082 <- p81[p81$F_A > 0.01 & p81$F_U > 0.01 & p81$info > 0.70,]

fil2 <- filter0082[filter0082$CHR == 22,]

filter0083 <- p61[p61$F_A > 0.01 & p61$F_U > 0.01 & p61$info > 0.70,]

fil3 <- filter0083[filter0083$CHR == 22,]

p96 <- read.table("Set1Set2Sample.TBL", header = TRUE)

p01 <- rbind.data.frame(filter0081[,c("CHR","SNP", "BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")])

p02 <- rbind.data.frame(filter0082[,c("CHR","SNP", "BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")])

p03 <- rbind.data.frame(filter0083[,c("CHR","SNP", "BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")])

p09 <- rbind.data.frame()

head(p01)

p97_1 <- cbind.data.frame(filter980, p01[match(filter980$MarkerName, p01$SNP),])

p97_2 <- cbind.data.frame(filter980, p02[match(filter980$MarkerName, p02$SNP),])

p97_4 <- cbind.data.frame(filter980, p03[match(filter980$MarkerName, p03$SNP),])

head(p97_4)

p97_5 <- na.omit(p97_4)

p97_7 <- na.omit(p97_6)

head(p97_7)

p97_7_1 <- p97_7[c("CHR","BP","SNP","A1","F_A","F_U","A2","info","P","OR","SE","L95","U95","EFFECT","HWE_ctrls","CHR_1","SNP_1","BP_1","A1_1","F_A_1","F_U1","A2_1","info_1","P_1","OR_1")]

head(p97_7)

SLEMetaColumbia <- p61[,c("SNP","A1", "A2", "OR", "F_A", "P", "SE")]
  
  
p97_6 <- cbind.data.frame(filter0081, p02[match(filter0081$SNP, p02$SNP),])

p99_1 <- na.omit(p97_6)

p99_2 <- p99_1[,c(1:15)]

head(p99_2)

p99_6 <- p99_2[,c(3,4,7,9,10,5)]

head(p99_6)

p99_6$OR <- log(p99_6$OR)
colnames(p99_6)[5] <- "EFFECT"
p99_6$N <- 465
head(p99_6)
p99_6 <- p99_6[,c(1,2,4,3,6,7,5)]

p99_6 <- p99_6[,c("SNP","A1","A2","N","EFFECT","F_A","P")]



head(p99_6)

head(p99_4)


p99_3 <- p99_1[,c(16:25)]

p99_4 <- p99_3[,c(2,4,7,5,9,10)]

head(p99_4)

p99_4$OR <- log(p99_4$OR)
colnames(p99_4)[6] <- "EFFECT"
p99_4$N <- 630
p99_4 <- p99_4[,c(1,2,3,7,6,4,5)]

head(p99_4)


p001_2 <- cbind.data.frame(p99_3, p02[match(p99_3$SNP, p03$SNP),])

p001_3 <- na.omit(p001_2)

head(p001_3)

head(p001_2)

head(p99_3)

head(p97_2)

head(p97_5)

p97_3 <- cbind.data.frame(filter889, p03[match(filter889$MarkerName, p03$SNP),])

p97 <- read.table("Set1Set2Set3Sample.TBL", header = TRUE)

p02_2 <- read.table("Set1Set2Set3Common.TBL", header = TRUE)

p09 <- read.table("Set1Set3")


head(p02_2)

p02_2 <- p02_2[order(p02_2$CHR, p02_2$BP),]

filter980 <- p02_2[p02_2$P.value < 0.00001 & p02_2$P.value <=1,]

head(filter980)

#p0081_0082 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP")], filter0082[,c("CHR", "BP", "SNP"), filter0082[,c("CHR", "BP", "SNP"), filter0083[,c("CHR", "BP", "SNP") ])

#p0081_0082_0083 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP"), filter0082[,c("CHR", "BP", "SNP"), filter0083[,c("CHR", "BP", "SNP"

p081_082_083 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP")], filter0082[,c("CHR", "BP", "SNP")], filter0083[,c("CHR", "BP", "SNP")] )









Set13 <- cbind.data.frame(Set13, p081_082_083[match(Set13$MarkerName, p081_082_083$SNP),])

Set23 <- cbind.data.frame(Set23, p081_082_083[match(Set23$MarkerName, p081_082_083$SNP),])


head(Set13)

Set13 <- Set13[order(Set13$CHR, Set13$BP),]

Set23 <- Set23[order(Set23$CHR, Set23$BP),]


p97 <- cbind.data.frame(p97, p081_082_083[match(p97$MarkerName, p081_082_083$SNP),])

p081_082 <- rbind.data.frame(filter0081[,c("CHR", "BP", "SNP")], filter0082[,c("CHR", "BP", "SNP")])


p96_0 <- cbind.data.frame(filter0083, p99_6[match(filter0083$SNP, p99_6$SNP),])

head(p96_0)

p96_1 <- na.omit(p96_0)

head(p96_1)

p07 <- p96_1[,c(1:15)]

p07_1 <- p07[,c(3,4,7,14,5,9)]

p07_1$N <- 1682

head(p07_1)

p07_2 <- p07_1[,c(1,2,3,7,4,5,6)]

head(p07_2)

head(p07)





head(p99_4)

head(p99_6)


head(p97)

head(p96)

p96 <- p96[order(p96$CHR, p96$BP),]

p97 <- p97[order(p97$CHR, p97$BP),]

p88 <- p88[order(p88$CHR, p88$BP),]

head(p88)

filter888 <- p96[p96$P.value < 0.00001 & p96$P.value <=1,]

filter900 <- p88[p88$P.value < 0.00001 & p88$P.value <=1,]

head(filter888)

filter889 <- p97[p97$P.value < 0.00001 & p97$P.value <=1,]

p61_1 <- filter0083[filter889$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p8888 <- filter0081[filter900$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p97_0 <- cbind.data.frame(filter900, p008_3[match(filter900$MarkerName, p008_3$SNP),])

p97_10 <- cbind.data.frame(filter900, p008_4[match(filter900$MarkerName, p008_4$SNP),])

head(p97_3)

p61_1 <- filter0083[filter889$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p81_1 <- filter0082[filter889$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p008_1 <- filter0081[filter889$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p81_2 <- filter0082[filter888$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p008_2 <- filter0081[filter888$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p8888 <- filter0081[filter900$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p8889 <- filter0082[filter900$MarkerName, c("CHR","BP","SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")]


filt292 <- write.table(p61_1,"ColImp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt293 <- write.table(p81_1,"12plateimp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt294 <- write.table(p008_1,"5plateimp.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt295 <- write.table(p81_2,"5plateimp_1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt296 <- write.table(p008_2,"12plateimp_1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt297 <- write.table(filter888,"SigSNPS2Sets.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt298 <- write.table(filter889,"SigSNPS3Sets.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt299 <- write.table(p97_1,"1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt300 <- write.table(p97_2,"2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt301 <- write.table(p97_3,"3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt302 <- write.table(p02_6,"Set1Set2Columns.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt303 <- write.table(p02_11,"Set1Set3Columns.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt304 <- write.table(p02_12,"Set2Set3Columns.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt305 <- write.table(GWAS11,"Snp1Set1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt306 <- write.table(GWAS12,"Snp1Set2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt307 <- write.table(GWAS13,"Snp1Set3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt308 <- write.table(GWAS14,"Snp2Set1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt309 <- write.table(GWAS15,"Snp2Set2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt310 <- write.table(GWAS16,"Snp2Set3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt311 <- write.table(fil1,"Set1_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt312 <- write.table(fil2,"Set2_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt313 <- write.table(fil3,"Set3_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt314 <- write.table(filter995,"rs34330Set2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt315 <- write.table(filter996,"rs34330Set1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt316 <- write.table(filter997,"rs34330Set3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt317 <- write.table(filter2118,"F_U0.txt", row.names = FALSE, quote = FALSE, sep = "\t")









head(p61_1)

head(filter889)

head(p92)

p62_63 <- rbind.data.frame(p62[,c("CHR", "BP", "SNP")], p63[,c("CHR", "BP", "SNP")])

p008_3 <- rbind.data.frame(filter0081[,c("SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")])

p008_4 <- rbind.data.frame(filter0082[,c("SNP", "A1", "F_A", "F_U", "A2", "info", "P", "OR")])

rownames(p62) <- p62$SNP
rownames(p63) <- p63$SNP

p92_1 <- p92[filter888$MarkerName, c("SNP","CHR" ,"BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")]

p61_11 <- p61[filter990$MarkerName, ,c("SNP", "CHR" , "BP" , "A1", "F_A", "F_U", "A2", "info", "P", "OR")]










head(chrXMetaSample)











head(p48)

head(p49)

head(p44)

head(p36)

##################################################################################################

#Extracting the Set1 GSA 5 Plate before imputation snps from the Set1 Imputed Data


p27$coordinate <- paste(p27$CHR, p27$BP, sep = ":")
p36$coordinate <- paste(p36$CHR, p36$BP, sep = ":")
p61$cordinate <- paste(p61$CHR, p61$BP, sep = ":")

head(p61)

p38 <- p36[match(p27$coordinate, p36$coordinate),]

p39 <- na.omit(p38)

p40 <- read.table("Columbian.txt",header = TRUE)

head(p40)

p41 <- na.omit(p40)

head(p41)


head(p38)

p29 <- unique(rbind.data.frame(p27[,c("CHR", "BP", "SNP")], p28[,c("CHR", "BP", "SNP")]))

p30 <- cbind.data.frame(p26, p29[match(p26$MarkerName, p29$SNP),])

p32 <- cbind.data.frame(p31, p29[match(p31$MarkerName, p29$SNP),])

write.table()


############################################################################################

filter211 <- p74[p74$P.value < 0.00001 & p74$P.value <=1,]






###############################################################################################################



filter10 <- NewMexFile[NewMexFile$F_A >= 0.01 & NewMexFile$F_U >= 0.01 & NewMexFile$F_A <= 1 & NewMexFile$F_U <= 1,]

filter11 <- NewYukFile[NewYukFile$F_A >= 0.01 & NewYukFile$F_U >= 0.01 & NewYukFile$F_A <= 1 & NewYukFile$F_U <= 1,]

filter12 <- NewPSSFile[NewPSSFile$F_A >= 0.01 & NewPSSFile$F_U >= 0.01 & NewPSSFile$F_A <= 1 & NewPSSFile$F_U <= 1,]

filter13 <- Newfile[Newfile$F_A >= 0.01 & Newfile$F_U >= 0.01 & Newfile$F_A <= 1 & Newfile$F_U <= 1,]

filter14 <- NewGSA[NewGSA$F_A >= 0.01 & NewGSA$F_U >= 0.01 & NewGSA$F_A <= 1 & NewGSA$F_U <= 1,]

filter15 <- NewGSAHG19[NewGSAHG19$F_A >= 0.01 & NewGSAHG19$F_U >= 0.01 & NewGSAHG19$F_A <= 1 & NewGSAHG19$F_U <= 1,]

filter16 <- NewGSAA[NewGSAA$F_A >= 0.01 & NewGSAA$F_U >= 0.01 & NewGSAA$F_A <= 1 & NewGSAA$F_U <= 1,]

filter17 <- New[New$F_A >= 0.01 & New$F_U >= 0.01 & New$F_A <= 1 & New$F_U <= 1,]

filter111 <- New2[New2$F_A >= 0.01 & New2$F_U >= 0.01 & New2$F_A <= 1 & New2$F_U <= 1,]

filter112 <- New4[New4$F_A >= 0.01 & New4$F_U >= 0.01 & New4$F_A <= 1 & New4$F_U <= 1,]

filter114 <- plink3[plink3$F_A >= 0.01 & plink3$F_U >= 0.01 & plink3$F_A <= 1 & plink3$F_U <= 1,]

filter115 <- p6[p6$F_A >= 0.01 & p6$F_U >= 0.01 & p6$F_A <= 1 & p6$F_U <= 1,]

filter116 <- p7[p7$F_A >= 0.01 & p7$F_U >= 0.01 & p7$F_A <= 1 & p7$F_U <= 1,]

filter117 <- p8[p8$F_A >= 0.01 & p8$F_U >= 0.01 & p8$F_A <= 1 & p8$F_U <= 1,]

filter118 <- p9[p9$F_A >= 0.01 & p9$F_U >= 0.01 & p9$F_A <= 1 & p9$F_U <= 1,]

filter119 <- p12[p12$F_A >= 0.01 & p12$F_U >= 0.01 & p12$F_A <= 1 & p12$F_U <= 1,]

filter120 <- p15[p15$F_A >= 0.01 & p15$F_U >= 0.01 & p15$F_A <= 1 & p15$F_U <= 1,]

filter121 <- p16[p16$F_A >= 0.01 & p16$F_U >= 0.01 & p16$F_A <= 1 & p16$F_U <= 1,]

filter122 <- p17[p17$F_A >= 0.01 & p17$F_U >= 0.01 & p17$F_A <= 1 & p17$F_U <= 1,]

filter123 <- p22[p22$A1 != 0 & p22$A2 != 0,]

x20 <- na.omit(filter119)



filter18 <- NewHG[NewHG$F_A >= 0.01 & NewHG$F_U >= 0.01 & NewHG$F_A <= 1 & NewHG$F_U <= 1,]

filter19 <- NewOBI[NewOBI$F_A >= 0.01 & NewOBI$F_U >= 0.0 &NewHG$F_A <= 1 & NewHG$F_U <= 1,]

filter29 <- plink1[plink1$F_A >= 0.01 & plink1$F_U >= 0.0 & plink1$F_A <= 1 & plink1$F_U <= 1,]


filter94 <- Mex[Mex$F_A >= 0.01 & Mex$F_U >= 0.01 & Mex$F_A <= 1 & Mex$F_U <= 1,]

filter95 <- Yuk[Yuk$F_A >= 0.01 & Yuk$F_U >= 0.01 & Yuk$F_A <= 1 & Yuk$F_U <= 1,]

filter96 <- Yuk[Yuk$F_A >= 0.01 & Yuk$F_U >= 0.01 & Yuk$F_A <= 1 & Yuk$F_U <= 1,]

filter97 <- PSS[PSS$F_A >= 0.01 & PSS$F_U >= 0.01 & PSS$F_A <= 1 & PSS$F_U <= 1,]

filter98 <- p14[p14$F_A >= 0.01 & p14$F_U >= 0.01 & p14$F_A <= 1 & p14$F_U <= 1,]



filter04 <- Newfile[Newfile$F_A >= 0.01 & Newfile$F_U >= 0.01 & Newfile$F_A <= 1 & Newfile$F_U <= 1 & Newfile$P <= 0.05 & Newfile$P <= 1,]

filter01 <- p23[p23$F_A >= 0.01 & p23$F_U >= 0.01 & p23$F_A <= 1 & p23$F_U <= 1 & p23$P.value <= 0.05 & p23$P.value <= 1,]

filter08 <- p10[p10$P <= 0.05 & p10$P <= 1,]

filter008 <- p24[p24$P.value <= 0.05 & p24$P.value <= 1,]

filter009 <- p25[p25$P.value <= 0.05 & p25$P.value <= 1,]

filter67 <- p13[p13$P.value <= 0.05 & p13$P.value <= 1,]

filter05 <- NewMexFile[NewMexFile$F_A >= 0.01 & NewMexFile$F_U >= 0.01 & NewMexFile$F_A <= 1 & NewMexFile$F_U <= 1 & NewMexFile$P <= 0.05 & NewMexFile$P <= 1,]

filter05 <- NewMexFile[NewMexFile$F_A >= 0.01 & NewMexFile$F_U >= 0.01 & NewMexFile$F_A <= 1 & NewMexFile$F_U <= 1 & NewMexFile$P <= 0.05 & NewMexFile$P <= 1,]

filter10 <- New1[New1$F_A >= 0.01 & New1$F_U >= 0.01 & New1$F_A <= 1 & New1$F_U <= 1 & New1$P <= 0.05 & New1$P <= 1,]

filter07 <- NewPSSFile[NewPSSFile$F_A >= 0.01 & NewPSSFile$F_U >= 0.01 & NewPSSFile$F_A <= 1 & NewPSSFile$F_U <= 1 & NewPSSFile$P <= 0.05 & NewPSSFile$P <= 1,]

filter08 <- Newobi[Newobi$F_A >= 0.01 & Newobi$F_U >= 0.01 & Newobi$F_A <= 1 & Newobi$F_U <= 1 & Newobi$P <= 0.05 & Newobi$P <= 1,]

filter99 <- MexGSA_5p[MexGSA_5p$F_A >= 0.01 & MexGSA_5p$F_U >= 0.01 & MexGSA_5p$F_A <= 1 & MexGSA_5p$F_U <= 1,]

filter01 <- MexGSA_5p[MexGSA_5p$F_A >= 0.01 & MexGSA_5p$F_U >= 0.01 & MexGSA_5p$F_A <= 1 & MexGSA_5p$F_U <= 1 & MexGSA_5p$P <= 0.05 & MexGSA_5p$P <= 1,]

filter001 <- NewMex[NewMex$F_A >= 0.01 & NewMex$F_U >= 0.01 & NewMex$F_A <= 1 & NewMex$F_U <= 1 & NewMex$P <= 0.05 & NewMex$P <= 1,]

filter002 <- NewMexx[NewMexx$F_A >= 0.01 & NewMexx$F_U >= 0.01 & NewMexx$F_A <= 1 & NewMexx$F_U <= 1,]


filter003 <- NewYukk[NewYukk$F_A >= 0.01 & NewYukk$F_U >= 0.01 & NewYukk$F_A <= 1 & NewYukk$F_U <= 1,]

filter004 <- NewPSSS[NewPSSS$F_A >= 0.01 & NewPSSS$F_U >= 0.01 & NewPSSS$F_A <= 1 & NewPSSS$F_U <= 1,]




















filter02 <- Mex_11[Mex_11$F_A >= 0.01 & Mex_11$F_U >= 0.01 & Mex_11$F_A <= 1 & Mex_11$F_U <= 1 & Mex_11$P <= 0.05 & Mex_11$P <= 1,]

filt234 <- write.table(filter94,"np6.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt235 <- write.table(filter95,"np7.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt237 <- write.table(filter97,"np8.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt238 <- write.table(filter99,"np9.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt239 <- write.table(filter01,"np10.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt240 <- write.table(filter02,"np11.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt241 <- write.table(filter001,"np12.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt242 <- write.table(filter002,"np13.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt243 <- write.table(filter003,"np14.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt244 <- write.table(filter004,"np15.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt245 <- write.table(filter04,"np16.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt246 <- write.table(filter04,"np17.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt247 <- write.table(filter05,"np18.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt248 <- write.table(filter06,"np19.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt249 <- write.table(filter07,"np20.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt250 <- write.table(filter08,"np21.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt251 <- write.table(filter10,"np22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt252 <- write.table(filter11,"np23.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt253 <- write.table(filter12,"np24.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt254 <- write.table(filter14,"np25.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt255 <- write.table(filter15,"np26.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt256 <- write.table(filter16,"np27.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt257 <- write.table(filter17,"np28.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt258 <- write.table(filter111,"np29.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt259 <- write.table(filter112,"np30.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt260 <- write.table(filter29,"np31.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt261 <- write.table(filter114,"np32.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt262 <- write.table(filter115,"np33.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt263 <- write.table(filter116,"np34.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt264 <- write.table(filter117,"np35.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt265 <- write.table(filter118,"np36.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt266 <- write.table(filter08,"np37.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt268 <- write.table(x20,"np38.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt269 <- write.table(filter67,"np39.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt270 <- write.table(filter123,"np40.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt271 <- write.table(filter008,"np42.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt272 <- write.table(p30,"SampleNew.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt273 <- write.table(p32,"InverseNew.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt274 <- write.table(p3,"InverseNew.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt275 <- write.table(p38,"Set1Imputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt276 <- write.table(p39,"Set1ImputedNoNA.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt277 <- write.table(p41,"ColumbianImputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt278 <- write.table(p71,"GSA5PlateImputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt279 <- write.table(p70,"ColumbianSLEImputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt280 <- write.table(filter2111,"Set1SigSNPS.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt281 <- write.table(filter2112,"ColumbianSigSNPS.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt282 <- write.table(p93,"GSA12plateImputed.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt283 <- write.table(p99_6,"Set1_MexicanCommon.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt284 <- write.table(p99_4,"Set2_MexicanCommon.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt285 <- write.table(p97_0,"Set1_MexicanCommonmeta.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt286 <- write.table(p97_10,"Set2_MexicanCommonmeta.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt287 <- write.table(filter900,"MetaResult.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt288 <- write.table(p07_2,"Set3_Columbianmeta.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt289 <- write.table(p02_6,"metaNew.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt290 <- write.table(p97_1,"Set1New.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt291 <- write.table(p97_2,"Set2New.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt292 <- write.table(p97_4,"ColumbianNew.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt293 <- write.table(fil1,"Set1_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt294 <- write.table(fil2,"Set2_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt295 <- write.table(fil3,"Set3_chr22.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt296 <- write.table(GWAS11,"Set1_1CHR.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt297 <- write.table(GWAS12,"Set2_1CHR.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt298 <- write.table(GWAS13,"Set3_1CHR.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt299 <- write.table(filter2112,"CasesSet1SigSNPS.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt300 <- write.table(p_controlsuniq,"ControlsUniq.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt301 <- write.table(p_casesuniq,"CasesUniq.txt", row.names = FALSE, quote = FALSE, sep = "\t")


x14 <- na.omit(chrX1)
head(x14)

x13 <- na.omit(chrX)

head(chrX)

head(GSplate)


head(GSA5platenoobi)

head(GSA12plate)

mm2 <- manhattanPlot(MexGSA_5p$P
                     , MexGSA_5p$CHR, 
                     ylim = c(0,15), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm3 <- manhattanPlot(Yuk$P
                     , Yuk$CHR, 
                     ylim = c(0,20), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm4 <- manhattanPlot(PSS$P
                     , PSS$CHR, 
                     ylim = c(0,50), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm5 <- manhattanPlot(Yuk_11$P
                     , Yuk_11$CHR, 
                     ylim = c(0,20), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm8 <- manhattanPlot(NewGSS$P
                     , NewGSS$CHR, 
                     ylim = c(0,25), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm9 <- manhattanPlot(NewMexx$P
                     , NewMexx$CHR, 
                     ylim = c(0,25), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm10 <- manhattanPlot(NewPS$P
                     , NewPS$CHR, 
                     ylim = c(0,25), trunc.lines = TRUE,
                     signif = (5e-8), thinThreshold= NULL)

mm11 <- manhattanPlot(filter18$P
                      ,filter18$CHR, 
                      ylim = c(0,70), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm12 <- manhattanPlot(New1$P
                      ,New1$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm13 <- manhattanPlot(plink$P
                      ,plink$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm14 <- manhattanPlot(plink1$P
                      ,plink1$CHR, 
                      ylim = c(0,20), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)



mm15 <- manhattanPlot(plink2$P
                      ,plink2$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm16 <- manhattanPlot(p5$P
                      ,p5$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm17 <- manhattanPlot(filter115$P
                      ,filter115$CHR, 
                      ylim = c(0,20), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm18 <- manhattanPlot(filter117$P
                      ,filter117$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm19 <- manhattanPlot(filter118$P
                      ,filter118$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm22 <- manhattanPlot(filter121$P
                      ,filter121$CHR, 
                      ylim = c(0,15), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm24 <- manhattanPlot(filter991$P
                      ,filter991$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)



mm20 <- manhattanPlot(p10$P
                      ,p10$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm21 <- manhattanPlot(p11$P
                      ,p11$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm22 <- manhattanPlot(p34$P
                      ,p34$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm23 <- manhattanPlot(p35$P
                      ,p35$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm24 <- manhattanPlot(p36$P
                      ,p36$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm25 <- manhattanPlot(p40$P
                      ,p40$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm26 <- manhattanPlot(p44$P.value
                      ,p44$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm27 <- manhattanPlot(p39$P
                      ,p39$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm28 <- manhattanPlot(p80$P
                      ,p80$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm29 <- manhattanPlot(p96$P
                      ,p96$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm30 <- manhattanPlot(filter0083$P
                      ,filter0083$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm31 <- manhattanPlot(Set13$P
                      ,Set13$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm32 <- manhattanPlot(Set23$P
                      ,Set23$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm33 <- manhattanPlot(p_cases$P
                      ,p_cases$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm34 <- manhattanPlot(p_controls$P
                      ,p_controls$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm35 <- manhattanPlot(p_casecontrolmerged$P
                      ,p_casecontrolmerged$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm36 <- manhattanPlot(filter2113$P
                      ,filter2113$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm37 <- manhattanPlot(p_casesuniq$P
                      ,p_casesuniq$CHR, 
                      ylim = c(0,20), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm38 <- manhattanPlot(p_controlsuniq1$P
                      ,p_controlsuniq1$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm39 <- manhattanPlot(p_controls_1$P
                      ,p_controls_1$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm40 <- manhattanPlot(filtercontrols$P
                      ,filtercontrols$CHR, 
                      ylim = c(0,20), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm41 <- manhattanPlot(filtercases$P
                      ,filtercases$CHR, 
                      ylim = c(0,20), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm42 <- manhattanPlot(TwoPlateMerge$P
                      ,TwoPlateMerge$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)


mm43 <- manhattanPlot(TwoPlateMerge_noATGC$P
                      ,TwoPlateMerge_noATGC$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm44 <- manhattanPlot(TwoPlateControl$P
                      ,TwoPlateControl$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)


mm45 <- manhattanPlot(TwoPlateCase$P
                      ,TwoPlateCase$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm46 <- manhattanPlot(ColumbianMeta$P
                      ,ColumbianMeta$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm47 <- manhattanPlot(Columbian$P
                      ,Columbian$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)


mm48 <- manhattanPlot(filter0083$P
                      ,filter0083$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm49 <- manhattanPlot(pc1_20$P
                      ,pc1_20$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm50 <- manhattanPlot(pc1_5$P
                      ,pc1_5$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm51 <- manhattanPlot(pc1_10$P
                      ,pc1_10$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm52 <- manhattanPlot(filter0082$P
                      ,filter0082$CHR, 
                      ylim = c(0,30), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm53 <- manhattanPlot(plink4$P
                      ,plink4$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)

mm54 <- manhattanPlot(GSANew$P
                      ,GSANew$CHR, 
                      ylim = c(0,25), trunc.lines = TRUE,
                      signif = (5e-8), thinThreshold= NULL)


qq(GSANew$P)


qq(p_casesuniq$P)

qq(p_controlsuniq1$P)

qq(p_cases$P)

qq(p_controls$P)

qq(p_casecontrolmerged$P)

qq(TwoPlateCase$P)

qq(TwoPlateMerge$P)

qq(TwoPlateMerge_noATGC$P)

qq(pc1_3$P)

qq(pc1_5$P)

qq(pc1_10$P)

qq(pc1_20$P)

qq(filter991$P)

qq(TwoPlateMerge_noATGC$P)

qq(plink4$P)

head(p97)



p50 <- na.omit(p44)






qq(p10$P)

qq(plink$P)

qq(plink1$P)

qq(p5$P)

qq(filter115$P)

qq(filter116$P)

qq(filter117$P)

qq(filter118$P)

head(plate5chrX)

qq(PSS$P)

qq(Mex_11$P)

qq(NewGSS$P)

qq(NewMexx$P)

qq(Yuk_11$P)

qq(NewPS$P)


library(qqman)
qqman::manhattan(filter121, chr="CHR", bp="BP", snp="MarkerName", p="P", ylim =c(0,20))
# # library(plotrix)
# axis.break(axis=2, breakpos=15.2, pos=NULL, bgcol="white", breakcol="black", style="zigzag", brw=0.02)
qqman::manhattan(p96, chr="CHR", bp="BP", snp="MarkerName", p="P.value", ylim =c(0,20))

qqman::manhattan(plate5chrX, chr="CHR", bp="POS", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(plate12chrX, chr="CHR", bp="POS", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(platemerge, chr="CHR", bp="POS", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(plate5, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(plate12, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(platemerge1, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(chrX5plate, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(chrX12plate, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(filter01, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,20))

qqman::manhattan(chrXMetaSample, chr="CHR", bp="BP", snp="SNP", p="P.value", ylim =c(0,20))

qqman::manhattan(chrXInverseSample, chr="CHR", bp="BP", snp="SNP", p="P.value", ylim =c(0,20))

qqman::manhattan(Set2chrX, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,10))

                
                
qqman::manhattan(Set2chrXunfilter, chr="CHR", bp="BP", snp="SNP", p="P", ylim =c(0,10))

qq(filte)

qq(Set2chrX$P)

qq(Set2chrXunfilter$P)


qq(filter120$P)


filter94 <- GSA5platenoobi[GSA5platenoobi$P < 0.0001 & GSA5platenoobi$P <= 1,]

filter91 <- GSA5platenoobi[GSA5platenoobi$CHR == 5 & GSA5platenoobi$BP == 159879978,]

filter92 <- GSA5platenoobi[GSA5platenoobi$CHR == 5 & GSA5platenoobi$BP == 159894847,]

filter93 <- GSA5platenoobi[GSA5platenoobi$CHR == 5 & GSA5platenoobi$BP == 159912418,]

filt230 <- write.table(filter91,"np1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt231 <- write.table(filter92,"np2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt232 <- write.table(filter93,"np3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt233 <- write.table(filter01,"chrXSet2_maf0005.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt234 <- write.table(filter02,"chrXSet1_maf005.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt235 <- write.table(filter991,"NewMexican.txt", row.names = FALSE, quote = FALSE, sep = "\t")

x10 <- na.omit(filter90)

x11 <- na.omit(filter94)

filt229 <- write.table(x11,"GSAnoObi10-4.txt", row.names = FALSE, quote = FALSE, sep = "\t")

qq(GSplate$P)




P_lambda(GSA5plate$P)










head(GSA5plate)

filter40 <- GSA5plate[GSA5plate$F_U > 0.01 & GSA5plate$info > 0.70,]

filter41 <- GSA5plate[GSA5plate$F_U > 0.03 & GSA5plate$info > 0.70,]

filter42 <- GSA5plate[GSA5plate$CHR == 5 & GSA5plate$BP == 159879978,]

filter43 <- GSA5plate[GSA5plate$CHR == 5 & GSA5plate$BP == 159894847,]

filter44 <- GSA5plate[GSA5plate$CHR == 5 & GSA5plate$BP == 159912418,]

filter45 <- GSA5plate[GSA5plate$CHR == 5 & GSA5plate$BP >= 159779978 & GSA5plate$BP <= 159979978,]


filt222 <- write.table(filter42,"snp1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt223 <- write.table(filter43,"snp2.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt224 <- write.table(filter44,"snp3.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt225 <- write.table(filter45,"snp4.txt", row.names = FALSE, quote = FALSE, sep = "\t")



Sletimvyse <- read.table("slegwasMetaCH1_22ld.txt")

eurset1set2 <- read.delim("EURSET1SET2.TBL",header = T)

eur1 <- read.table("EURSET1SET2.TBL",header = T)

eur2 <- read.table("SLEEUR.TBL",header = T)

head(eurset1set2)

head(eur2)

filter100 <- eur1[eur1$Direction == '+??',]

filter101 <- eur1[eur1$Direction == '-??',]

filter102 <- eur1[eur1$Direction == '--?',]

filter103 <- eur1[eur1$Direction == '-+?',]

filter104 <- eur1[eur1$Direction == '+-?',]

filter105 <- eur1[eur1$Direction == '++?',]

filter106 <- eur1[eur1$Direction == '?--',]

filter107 <- eur1[eur1$Direction == '?-+',]

filter108 <- eur1[eur1$Direction == '?++',]

filter109 <- eur1[eur1$MarkerName == 'rs17162854',]

filter110 <- eur2[eur2$MarkerName == 'rs7643092',]



head(GWAS77)

colnames(GWAS77) <- c("CHR", "SNP", "BP", "A1", "A2", "OR", "OR_lower", "OR_upper", "P")

head(GWAS67)

GSA_All <- read.table("GSA12plate_All.hwe",header = TRUE)

head(GSA_All)

filter75 <- GSA_All[GSA_All$TEST == AFF,]

filter29 <- GWAS77[GWAS77$CHR == 21, ]
filter30 <- GWAS77[GWAS77$CHR == 20, ]
filter31 <- GWAS77[GWAS77$CHR == 19, ]
filter32 <- GWAS77[GWAS77$CHR == 18, ]
filter33 <- GWAS77[GWAS77$CHR == 17, ]
filter34 <- GWAS77[GWAS77$CHR == 16, ]
filter35 <- GWAS77[GWAS77$CHR == 15, ]
filter36 <- GWAS77[GWAS77$CHR == 14, ]
filter37 <- GWAS77[GWAS77$CHR == 13, ]
filter38 <- GWAS77[GWAS77$CHR == 12, ]
filter39 <- GWAS77[GWAS77$CHR == 11, ]
filter40 <- GWAS77[GWAS77$CHR == 10, ]
filter41 <- GWAS77[GWAS77$CHR == 9, ]
filter42 <- GWAS77[GWAS77$CHR == 7, ]
filter43 <- GWAS77[GWAS77$CHR == 6, ]
filter44 <- GWAS77[GWAS77$CHR == 5, ]
filter45 <- GWAS77[GWAS77$CHR == 4, ]
filter46 <- GWAS77[GWAS77$CHR == 3, ]
filter47 <- GWAS77[GWAS77$CHR == 2, ]
filter48 <- GWAS77[GWAS77$CHR == 1, ]


filter49 <- GWAS77[GWAS77$CHR >= 23 & GWAS77$CHR <= 25, ]


f7 <- write.table(GWAS77,"F9.txt", quote = FALSE)


AA_GWAS4 <- read.delim("META2.txt",header = T)

head(AA_GWAS4)

colnames(AA_GWAS4) <- c("CHR", "BP", "SNP", "A2", "F_A", "F_U", "A1", "info", "P", "OR", "SE", "L95", "U95", "EFFECT", "HWE_ctrls")

head(AA_GWAS4)

GWAS4 <- read.csv("summary.result.observatory_dataset.csv",header = T)



head(GWAS4)



GWAS5 <- read.table("AAgwasdata2.assoc",header = T)

GWAS6 <- read.table("AAgwasdata1qc.assoc",header = T)

AGWAS77 <- read.delim("ChineseGWASmeta.txt", header = T)

GWAS7 <- read.table("AAgwasdata2qc.assoc",header = T)

GWAS8 <- read.table("AAgwasAdd.assoc",header = T)

GWAS9 <- read.table("Height4.QC.gz",header = T)

spain <- read.table("Spain_Results.txt", header = T)

meta <- read.table("SUMMARYFILE.txt.gz", header = T)

lupus <- read.delim("PASSLUPUS.tsv", header = T)

GSANew <- read.table("GSA12plateMexAutosomeNew.assoc", header = T)

head(meta)

head(lupus)

BMI <- read.table("body_BMIz.sumstats.gz", header = T)

BMIPAPER <- read.table("BMIText.txt", header = T)

Genes94 <- read.table("94Genes.txt", header = T)

Sczpaper2 <- read.table("SchizopreniasummarywithchrX.txt", header = T)

head(Sczpaper)

colnames(Sczpaper) <- c("CHR", "SNP", "A1", "A2", "BP", "INFO", "OR", "SE", "P", "Ngt")

filter81 <- Sczpaper[Sczpaper$CHR == chr1]

head(BMIPAPER)

pcaplot <- ggplot(BMIPAPER, aes(colour = Tissue)) + geom_point() 
pcaplot


data(diamonds)
p.dia <- ggplot(data = diamonds, mapping = aes(x = clarity))

p <- p.dia + layer(geom = "bar", mapping = aes(fill = cut))

GSAAllPop <- read.table("Set2AllPopPCA.txt",header = TRUE)

GSANewPop <- read.table("GSANewfilteredPCA1-3.txt",header = TRUE)

df <- read.table("GSAPCA.txt",header = TRUE)

df1 <- read.table("GSSPCA.txt",header = TRUE)

df2 <- read.table("MergedDataPCA.txt",header = TRUE)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(9, "Set2"))(nb.cols)
mycolors <- colorRampPalette(c("black", "red"))(9)

# Create a ggplot with 18 colors 
# Use scale_fill_manual
ggplot(GSAAllPop) + 
  geom_point(aes(x=PC1, y=PC2, col = POP)) +
  scale_color_manual(values = c("brown", "green", "blue", "darkgrey", "yellow", "orange", "black", "red", "purple")) +
  theme_minimal() +
  theme(legend.position = "right")


ggplot(GSAAllPop) + 
  geom_point(aes(x=PC2, y=PC3, col = POP)) +
  scale_color_manual(values = c("brown", "green", "blue", "darkgrey", "yellow", "orange", "black", "red", "purple")) +
  theme_minimal() +
  theme(legend.position = "right")


ggplot(GSAAllPop) + 
  geom_point(aes(x=PC1, y=PC3, col = POP)) +
  scale_color_manual(values = c("brown", "green", "blue", "darkgrey", "yellow", "orange", "black", "red", "purple")) +
  theme_minimal() +
  theme(legend.position = "right")


















pp <- ggplot() +
  geom_jitter(data=BMIPAPER, aes(x=x, y=BMIPAPER$X.log10p, color=BMIPAPER$Tissue)) 

pp
pp + coord_cartesian(x = c(0,250))

ggplot(BMIPAPER,as.data.frame(list(x = c(0,50,100,150,200), y = c(0,1,2,3,4,5,6,7,8))),
       aes(x = x, y = y, color = BMIPAPER$Tissue)) + geom_point()

pcaplot4 <- ggplot(GSANewPop, aes(x=PC1, y=PC2, colour = POP)) + geom_point() 
pcaplot4

pcaplot8 <- ggplot(GSANewPop, aes(x=PC2, y=PC3, colour = POP)) + geom_point() 
pcaplot8

pcaplot9 <- ggplot(GSANewPop, aes(x=PC1, y=PC3, colour = POP)) + geom_point() 
pcaplot9

pcaplot10 <- ggplot(plink6, aes(x=plink6$PC1, y=plink6$PC2, colour = plink6$POP)) + geom_point() 
pcaplot10

pcaplot11 <- ggplot(plink6, aes(x=plink6$PC2, y=plink6$PC3, colour = plink6$POP)) + geom_point() 
pcaplot11

pcaplot12 <- ggplot(plink6, aes(x=plink6$PC1, y=plink6$PC3, colour = plink6$POP)) + geom_point() 
pcaplot12





pcaplot9 <- ggplot(df2, aes(x=df2$PC1, y=df2$PC3, colour = df2$POP)) + geom_point() 
pcaplot9

pcaplot10 <- ggplot(df1, aes(x=df1$PC1, y=df1$PC3, colour = df1$POP)) + geom_point()
pcaplot10



pcaplot5 <- ggplot(BMIPAPER, aes(x=, y=aa_gwas_maf_data$PC3, colour = aa_gwas_maf_data$Ph)) + geom_point() 
pcaplot5

ggplot(as.data.frame(list(x = c(0, 200,100), y = c(7500000,10000000,2000000))), 
       aes(x = x, y = y)) +
  geom_point()


filter9 <- spain[spain$CHR >= 17 & spain$CHR <=22,]

filter28 <- meta[meta$CHR >= 12 & meta$CHR <= 22,]

filter29 <- meta[meta$CHR = 1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,]

filter43 <- meta[meta$CHR == 6 & meta$BP >= 26000026 & meta$BP <= 33999992,]

meta1 <- read.table("filead.txt", header = T)

meta2 <- read.table("fileae.txt", header = T)

filter40 <- meta1[meta1$CHR == 4:5,]

head(filter40)

##############################################################################

set1 <- read.delim("GWASset1.txt", header = T)

set2 <- read.delim("GWASset2.txt", header = T)

set3 <- read.delim("ChineseGWASmeta.txt", header = T)

SET4 <- read.delim("MetaAnalysis.txt", header = T)




filter41 <- meta2[meta2$CHR >= 7 & meta2$CHR <= 8,]



head(filter27)


head(filter5)

GWAS17 <- read.table("pgc.bip.full.2012-04.txt", header = T)

GWAS18 <- read.table("pgc.scz.full.2012-04.txt", header = T)

GWAS21 <- read.table("bipHeight.noambig.gz", header = T)

head(GWAS21)

x5 <- na.omit(GWAS17)

x6 <- na.omit(GWAS18)

colnames(GWAS21) <- c("SNP", "CHR", "BP", "A1", "A2", "OR", "SE", "P", "INFO", "Ngt", "CEUaf")

memory.limit()
memory.limit(10000000000000)

head(x3)

data1 <- as.data.frame(cbind(x3$CHR,x3$BP,x3$SNP,x3$P));

data1 <- as.data.frame(cbind(CHR,BP,SNP,P))

#vek <- as.numeric(gsub(data1$V1.pattern="chr",replacement=""))%%2)

PP1 <- P_lambda(x11$P)


head(filter23)

x1 <- na.omit(AA_GWAS)

x2 <- na.omit(GWAS1)

x3 <- na.omit(AA_GWAS4)

AA_GWAS2 <- read.table("allpop.txt",header = T)
gg <- combine.test(AA_GWAS2$P,weight = W,method = "fisher", na.rm = TRUE)


GG2 <- pchisq(AA_GWAS2$P, AA_GWAS2, ncp = 0, lower.tail = TRUE, log.p = FALSE)


require(sqldf)

CHR2_1 <- sqldf("select * from AA_GWAS3 where CHR = 2")

head(CHR2_1)

head(x1)

head(x2)

filter22 <- x2[x2$F_U > 0.01, ]

filter30 <- x1[x1$F_U > 0.01 & x1$info > 0.70,]

filter24 <- x2[x2$F_U > 0.03,]

filter31 <- x2[x2$F_U >= 0.005 & x2$F_A >= 0.005 & x2$P < 0.00001,]


filter33 <- x3[x3$F_U > 0.01 & x3$info > 0.70,]

filter32 <- x3[x3$F_U > 0.03 & x3$info > 0.70,]

filter20 <- x2[x2$info > 0.70 & x2$info <= 1,]

filter21 <- x1[x1$info > 0.70 & x1$info <= 1,]

filter26 <- x3[x3$P < 0.00001 & x3$P <=1,]

filter35 <- x2[x2$P < 0.00000001 & x2$P <=1,]

filter36 <- x3[x3$P < 0.00000001 & x3$P <=1,]

filter2 <- GWAS1[GWAS1$P < 0.00001 & GWAS1$P <=1,]


filter1 <- x4[x4$CHR == 22 & x4$BP >= 18885242 & x4$BP <= 39796097,]

filter58 <- x3[x3$P < 0.001 & x3$info > 0.70,]

filter4 <- GWAS[GWAS$F_U > 0.005,]

filt <- x3[x3$P < 0.05 & x3$info > 0.70,]

filt1 <- filt[filt$CHR >=1 & filt$CHR <= 4,]

filt2 <- filt[filt$CHR >=5 & filt$CHR <=10,]

filt3 <- filt[filt$CHR >=11 & filt$CHR <=16,]

filt4 <- filt[filt$CHR >=17 & filt$CHR <=22,]


filt2 <- write.table(filter41,"Meta78.txt", row.names = FALSE, quote = FALSE, sep = "\t")

filt6 <- write.csv(filt2,"dat2.csv", row.names = TRUE)

filt7 <- write.csv(filt3,"dat3.csv", row.names = TRUE)

filt8 <- write.csv(filt4,"dat4.csv", row.names = TRUE)

#######################################################################


##############################################


filter46 <- x2[x2$BP == "131803537",]

filter47 <- x3[x3$BP == "31276811",]

filter49 <- x2[x2$BP == "74208362",]

filter51 <- x3[x3$BP == "74208362",]

filter52 <- x3[x3$BP == "74217856",]

filter43 <- x3[x3$SNP == "rs17622517",]

filter44 <- x3[x3$SNP == "rs17622517",]


filter11 <- write.csv(filter5,"B3.csv", row.names = TRUE)

filter10 <- write.csv(filter2,"B2.csv", row.names = TRUE)

filter40 <- write.csv(filter26,"data12.csv", row.names = TRUE)

filter54 <- write.csv(filter51,"data14.csv", row.names = TRUE)

filter55 <- write.csv(filter50,"data15.csv", row.names = TRUE)

filter60 <- write.csv(filt,"data25.csv", row.names = TRUE)

head(filter5)

head(filter6)


filter3 <- CHR2_1[CHR2_1$BP > 73000554 & CHR2_1$BP < 75999874,]

filter4 <- write.csv(filter31,"t13.csv", row.names = TRUE)

filter8 <- write.csv(filter7,"data2.csv", row.names = TRUE)

filter10 <- write.csv(filter9,"data3.csv", row.names = TRUE)

first_position <- "73813530"
second_positon <- "74735302"


CCCHR2_TET3 <- sqldf("select * from CHR2 where BP in (74820862,76628006)")

head()

manhattan(AA_GWAS, chr="CHR", bp="BP", snp = "SNP", p = "P", genomewideline = 5e-12, suggestiveline = 5e-8)

manhattanPlot(AA_GWAS, chr="CHR", bp="BP", snp = "SNP", p = "P", signif=-log10(5e-8), logp = FALSE, ylim = c(0,20) )

manhattanPlot(AA_GWAS, chromosome = "AA_GWAS$CHR", bp="BP", snp = "SNP", p = "AA_GWAS$P", signif=5e-8, ylim = c(20) )

MainIn <- read.delim("HCC.txt",header = TRUE, sep = "\t")

MainIn1 <- read.delim("HC_GWAS_MetaSam.txt",header = TRUE, sep = "\t")

MainIn4 <- read.delim("file8.txt",header = TRUE, sep = "\t")

MainIn6 <- read.delim("file99.txt",header = TRUE, sep = "\t")

MainIn10 <- read.delim("hj5.annot",header = TRUE, sep = "\t")

MainIn14 <- read.delim("MetaCircle.txt",header = TRUE, sep = "\t")

MainIn15 <- read.delim("MetaCircle3.txt",header = TRUE, sep = "\t")

MainIn20 <- read.delim("MetaCircle10.txt",header = TRUE, sep = "\t")

MainIn21 <- read.delim("MetaCircle10.txt",header = TRUE, sep = "\t")

MainIn23 <- read.delim("MetaCircle14.txt",header = TRUE, sep = "\t")

CO2.data = read.csv2("MetaCircle14.txt", sep=";", dec = ",", stringsAsFactors=FALSE)

class(CO2.data$BP)

MainIn11 <- MainIn10[order(MainIn10$CHR, MainIn10$BP),]

SNPS7 <- read.table("SNP6.txt",header = TRUE)

SNP8 <- 

S61 <- read.table("SP.txt",header = TRUE)

length(S61)

S61


S41 <- read.table("S6.txt",header = TRUE,sep = "\t" )

SP15 <- list(MainIn15$SNP[MainIn15$trait1<1e-6])

genes <- paste("GENE", 1:44, sep="_")

gene1 <- read.table("Gene1.txt",header = TRUE,sep = "\t")

genes


head(SNPS7)

head(MainIn11)

filter333 <- MainIn5[MainIn5$P.value < 0.05 & MainIn5$P.value <=1,]

filter334 <- MainIn7[MainIn7$P.value < 0.05 & MainIn7$P.value <=1,]

filter335 <- MainIn5[MainIn5$P.value <= 0.001 & MainIn5$P.value <=1,]

filter336 <- MainIn7[MainIn7$P.value <= 0.001 & MainIn7$P.value <=1,]

fi23 <- write.table(fill3,"HCMetasample1.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi25 <- write.table(filter333,"HCSample.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi26 <- write.table(filter334,"HCSE.txt", row.names = FALSE, quote = FALSE, sep = "\t")



head(MainIn4)

MainIn5 <- MainIn4[order(MainIn4$CHR, MainIn4$BP),]

MainIn7 <- MainIn6[order(MainIn6$CHR, MainIn6$BP),]

SNPss <- list(SNPS7(SNPS7$SNP[SNPS7$Tier == "Tier1",],SNPS7(SNPS7$SNP[SNPS7$Tier == "Tier2",],)))

head(MainIn7)

source("https://raw.githubusercontent.com/YinLiLin/CMplot/master/R/CMplot.r")

Main8 <- CMplot(MainIn15, plot.type="m",col=c("grey30","grey60"),highlight=S6,
       highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=genes,      
       highlight.text.col=c("red","blue","green"),threshold=c(1e-05,1e-08),threshold.lty=2,   
       amplify=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)

Main10 <- CMplot(MainIn15, plot.type="m",col=c("grey60","grey30"),highlight=S61,
                highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=gene1,      
                highlight.text.col=c("red","blue","green"),threshold=c(1e-05,1e-08),threshold.lty=2) 

Main17 <- CMplot(MainIn15, plot.type="m",col=c("grey60","grey30"),highlight=S61,
                 highlight.col=c("red","blue","green"),highlight.cex=1,highlight.pch=c(15:17), highlight.text=gene1, ylab.pos=3, xticks.pos=1,   
                 highlight.text.col=c("red","blue","green"),threshold=c(1e-05,1e-08),threshold.lty=2) 




plot(mm2, main = "Axis break test")                
                
axis.break(2, 2.9, style = "zigzag")    

class(MainIn20$CHR)
                
MainIn22$CHR <- as.numeric(MainIn20$CHR)
MainIn22$BP <- as.numeric(MainIn20$BP)

class(MainIn22$CHR)


                
Main9 <- CMplot(MainIn15,type = "p", plot.type="c",chr.labels=paste("Chr",c(1:23),sep=""),highlight=S6,outward = FALSE,
                highlight.col=c("red","blue","green"),cir.legend.col="black",cir.chr.h=1,chr.den.col="black",r=1,cir.legend=TRUE,highlight.cex=1,highlight.pch=c(15:17), highlight.text=genes,      
                highlight.text.col=c("red","blue","green"),threshold=c(1e-05,1e-08),threshold.lty=2,   
                amplify=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=14,height=6)


Main2 <- CMplot(MainIn15,type="p",plot.type="c",chr.labels=paste("Chr",c(1:23),sep=""),threshold=c(1e-08),r=1,cir.legend=TRUE,highlight = SNPS3,highlight.col = "green",highlight.cex = 1,highlight.pch = c(15,17),
       outward=FALSE,cir.legend.col="black",cir.chr.h=1,threshold.col=c("red"),chr.den.col="black",file="jpg",
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

Main3 <- CMplot(MainIn15,type="p",plot.type="c",chr.labels=paste("Chr",c(1:23),sep=""),threshold=c(1e-08),r=1,cir.legend=TRUE,highlight = S61,highlight.col = "green",highlight.text=genes,highlight.text.col="red",highlight.cex = 1,highlight.pch = 19,
                outward=FALSE,cir.legend.col="black",cir.chr.h=1,threshold.col="red",chr.den.col="black",file="jpg",bin.size = 1e6,
                memo="",dpi=300,file.output=TRUE,axis.break(2, 2.9, style = "zigzag"),verbose=TRUE,width=14,height=14)

Main4 <- CMplot(MainIn15,type="p",plot.type="m",chr.labels=paste("Chr",c(1:23),sep=""),threshold=c(1e-08),highlight = SP15,highlight.col = "green",highlight.text=SP15,highlight.text.col="red",highlight.cex = 1,highlight.pch = c(15,17),amplify = TRUE,
                outward=FALSE,threshold.col=c("red"),chr.den.col="black",file="jpg",bin.size = 1e6,signal.col = "green",signal.cex = 1,
                memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


class(MainIn22$CHR)

SeaAnnot <- read.delim("SeaAnnot.txt",header = TRUE, sep = "\t")

head(SeaAnnot)

col_names = c("position","rsID","functionGVS","geneList")

SeaAnnot1 <- SeaAnnot[, col_names]

fi266 <- write.table(SeaAnnot1,"SeaAnnot1.txt", row.names = FALSE, quote = FALSE, sep = "\t")


head(MainIn7)

########################################################################################################################

##################################################################################################

HC_GWAS2 <- read.delim("HCGWAS_2_Sample.TBL",header = TRUE, sep = "\t")

HC_GWAS2SE <- read.delim("HC_GWAS_2_SE.TBL",header = TRUE, sep = "\t")

HC_GWAS3 <- read.delim("HCGWAS_3_Sample.TBL",header = TRUE, sep = "\t")

HC_GWAS3SE <- read.delim("HCGWAS_3_SE.TBL",header = TRUE, sep = "\t")

HC_GWAS4 <- read.delim("HCGWAS_4_Sample.TBL",header = TRUE, sep = "\t")

HC_GWAS4SE <- read.delim("HC_GWAS_4_SE.TBL",header = TRUE, sep = "\t")

#10-3
fi1 <- HC_GWAS2[HC_GWAS2$P <= 0.001 & HC_GWAS2$P <=1,]

fi2 <- HC_GWAS2SE[HC_GWAS2SE$P <= 0.001 & HC_GWAS2SE$P <=1,]

fi3 <- HC_GWAS3[HC_GWAS3$P <= 0.001 & HC_GWAS3$P <=1,]

fi4 <- HC_GWAS3SE[HC_GWAS3SE$P <= 0.001 & HC_GWAS3SE$P <=1,]

fi5 <- HC_GWAS4[HC_GWAS4$P <= 0.001 & HC_GWAS4$P <=1,]

fi6 <- HC_GWAS4SE[HC_GWAS4SE$P <= 0.001 & HC_GWAS4SE$P <=1,]


fi7 <- write.table(fi1,"HC2Sample.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi8 <- write.table(fi2,"HC2SE.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi9 <- write.table(fi3,"HC3Sample.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi10 <- write.table(fi4,"HC3SE.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi11 <- write.table(fi5,"HC4Sample.txt", row.names = FALSE, quote = FALSE, sep = "\t")

fi12 <- write.table(fi6,"HC4SE.txt", row.names = FALSE, quote = FALSE, sep = "\t")


##############################################################################

## make the plot 
ggplot(MainIn23, aes(x=MainIn23$CHR, y=-log10(MainIn23$P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(MainIn23$CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#81D3EB", "#0C234B"), each = 1, len = 23)) +
  
  # custom X axis:
  scale_x_continuous( label = MainIn23$CHR) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 25, 35, 50)) +     # remove space between plot area and x axis
  
  # add genome-wide sig and sugg lines
  #geom_hline(yintercept = -log10(sig), color = "gray") +
  #geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray") 
  
  # Add highlighted points
  #geom_point(data=subset(don2, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Add label using ggrepel to avoid overlapping
 # geom_label_repel( data=subset(don2, is_annotate=="yes"), aes(label= topgenes$name2), size=3,
                    #segment.color = "transparent", nudge_x = -2) +
  
  #expand_limits(y = 10) + 
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(x = "", y = "") 



ggsave("manhattan.png", 
       plot = last_plot(), 
       device = "png", 
       width = 12, height = 6, units = "in")

######################################################################################################################

## load packages
library(readr)
library(ggrepel)
library(tidyverse)
library(data.table)


## load data desktop

dat <- read.delim("MetaCircle10.txt",header = TRUE, sep = "\t")

## load gene refernce data


gene_result = gene_result %>% 
  mutate(chrom = str_replace(gene_result$chrom, "chr", "")) %>% 
  select(-name)


## define significance 
sig = 5e-8 # significant threshold line
sugg = 1e-5 # suggestive threshold line




don <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP, -desc(P)) %>%
  mutate( BPcum=BP+tot) %>% 
  distinct(SNP, .keep_all = T)


## REMOVE P FILTER BEFORE PUBLISHING
don2 = don %>% filter(CHR <= 22 & CHR > 0,
                      P < .1)

print("CONGRATULATIONS! Data was cleaned")

## for testing only comment out for plots 

## define genes for annotating
## just top snp
topsnp = don2[which.min(don2$P),]

## top SNP per chromosome
topsnps = don2 %>%  group_by(CHR) %>% slice(which.min(P)) %>% filter(P < sugg)

## this gives just TOP Gene
topgene = gene_result %>% filter(topsnp$CHR == gene_result$chrom & 
                                   topsnp$BP >  gene_result$txStart &
                                   gene_result$txEnd > topsnp$BP)


## top genes per chromosome
topgene_out_fn = function(topsnps){
  gene_result %>%     filter(topsnps$CHR == gene_result$chrom &
                               topsnps$BP >  gene_result$txStart - 500000 &                                     gene_result$txEnd + 500000 > topsnps$BP)
  
}

topgene_in_fn = function(topsnps){
  gene_result %>%
    filter(topsnps$CHR == gene_result$chrom &
             topsnps$BP >  gene_result$txStart &
             gene_result$txEnd  > topsnps$BP)
  
}

topgenes_out = by(topsnps, topsnps$CHR, function(topsnps) topgene_out_fn(topsnps))
topgenes_in = by(topsnps, topsnps$CHR, function(topsnps) topgene_in_fn(topsnps))


topgenes_out = do.call(rbind.data.frame, topgenes_out)
topgenes_in = do.call(rbind.data.frame, topgenes_in)

## remove dups? 
is_duplicate_out <- lapply(X = topgenes_out, FUN = duplicated, incomparables = FALSE)


drop_idx_out <- which(is_duplicate_out$name)
topgenes_out = topgenes_out[-drop_idx_out, ]


is_duplicate_in <- lapply(X = topgenes_in, FUN = duplicated, incomparables = FALSE)


drop_idx_in<- which(is_duplicate_in$name)
topgenes_in = topgenes_in[-drop_idx_in, ]





topgenes = list()

for (i in 1:nrow(topsnps)) {
  if(topsnps$CHR[i] %in% topgenes_in$chrom){
    topgenes[[i]] = topgenes_in[topgenes_in$chrom == topsnps$CHR[i],]}
  else{
    targetgenes =  topgenes_out[topgenes_out$chrom == topsnps$CHR[i],]
    
    ##Get the distance to each target gene
    distances<-list()
    for(j in 1:nrow(targetgenes)){
      if(targetgenes$exonStarts[j] > topsnps$BP[i]){
        distances[[j]] <-targetgenes$exonStarts[j]-topsnps$BP[i]
      }else{
        distances[[j]] <- topsnps$BP[i]-targetgenes$exonEnds[j]
      }
    }
    
    topgenes[[i]] =targetgenes[which.min(unlist(distances)),]
  }
  print(i) 
}

topgenes = do.call(rbind,topgenes)


## define snps of interest
snpsOfInterest = don2 %>% filter(topsnp$BP -50000 <  don2$BP & don2$BP < topsnp$BP +50000 & 
                                   don2$CHR == topsnp$CHR)

don2 = don2 %>% 
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(BP %in% snpsOfInterest$BP, "yes", "no")) %>% 
  mutate( is_annotate=ifelse(SNP %in% topsnps$SNP, "yes", "no")) 


### prepare x axis 
## chromosome name 
axisdf = don2 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


## make the plot 
ggplot(don2, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#81D3EB", "#0C234B"), each = 1, len = 22)) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 2, 4, 6, 8)) +     # remove space between plot area and x axis
  
  # add genome-wide sig and sugg lines
  geom_hline(yintercept = -log10(sig), color = "gray") +
  geom_hline(yintercept = -log10(sugg), linetype="dashed", color = "gray") +
  
  # Add highlighted points
  geom_point(data=subset(don2, is_highlight=="yes"), color="#EF4056", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don2, is_annotate=="yes"), aes(label= topgenes$name2), size=3,
                    segment.color = "transparent", nudge_x = -2) +
  
  expand_limits(y = 10) + 
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  ) +
  labs(x = "", y = "") 



ggsave("manhattan.png", 
       plot = last_plot(), 
       device = "png", 
       width = 12, height = 6, units = "in")











#########################################################################################


#filtermeta <- MainIn5[MainIn5 <= ]

head(MainIn15)

head(MainIn)

class(MainIn15$trait1)
class(MainIn15$chr)

MainIn15$chr <- as.numeric(MainIn15$chr)

MainIn15$trait1 <- -log10(MainIn15$trait1)

MainIn15

mm2 <- manhattanPlot(MainIn15$trait1
            , MainIn15$chr, 
              ylim = c(0,250), trunc.lines = TRUE,
              signif = (5e-8), thinThreshold= NULL)

axis(2, at=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,100,150,200))

plot(plot10.png, main = "Axis break test")

install.packages("GWASinspector")
library("GWASinspector")

manhattan.plot(dataset = MainIn15, chr = 'chr', pvalue = 'trait1', position = 'pos',
               plot.title = 'Manhattan plot', plot.subtitle = 'This data is fabricated!',
               fileName = png , p.threshold = '0.5')


install.packages("EBImage")

install.packages("png")
library("png")

install.packages('magick')
library("magick")

img <- magick::image_read('plot10.png')


install.packages("tidyverse")
library("tidyverse")


############################################################################################################################

sig.dat <- MainIn15 %>% 
  subset(MainIn15$trait1 < 0.05)
notsig.dat <- MainIn15 %>% 
  subset(MainIn15$trait1 >= 0.05) %>%
  slice(sample(nrow(.), nrow(.) / 5))
gwas.dat <- rbind(sig.dat,notsig.dat)

nCHR <- length(unique(gwas.dat$chr))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()
for (i in unique(gwas.dat$chr)){
  nbp[i] <- max(gwas.dat[gwas.dat$chr == i,]$pos)
  gwas.dat[gwas.dat$chr == i,"BPcum"] <- gwas.dat[gwas.dat$chr == i,"BP"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>% 
  group_by(gwas.dat$chr) %>% 
  
ylim <- abs(floor(log10(min(gwas.dat$trait1)))) + 2 
sig <- 5e-8




############################################################################################################################
plot(img)

plot(img, main = "Axis break test")
class(img)
# put a break at the default axis and position
axis.break()
axis.break(2, 0+15, breakcol="black", style="zigzag")
axis.break(2, 15, style = "zigzag")

plot(3:10, main = "Axis break test")
# put a break at the default axis and position
axis.break()
axis.break(2, 2.9, style = "zigzag")
twogrp <- c(rnorm(10) + 4, rnorm(10) + 20)
gap.plot(twogrp,gap = c(8,16), xlab = "Index", ylab = "Group values",
         main = "Two separated groups with gap axis break",
         col = c(rep(2, 10), rep(3, 10)), ytics = c(3, 5, 18, 20))
legend(12, 6, c("Low group", "High group"), pch = 1, col = 2:3)


qq(MainIn$P)

head(GWAS1)

library(manhattanly)

manhattanly(x3, snp = "SNP", gene = "GENE")


qq(GWAS$P)

qq(aa_$P, main = "Q-Q plot of GWAS p-values", xlim = c(0, 7), ylim = c(0, 
                                                                               12), pch = 18, col = "blue4", cex = 1.5, las = 1)

manhattan(AA_GWAS, annotatePval = 0.01)

hwe <- read.table("AA_GWAS_chr6.hwe", header = T, stringsAsFactors = F)

str(hwe)

nrow(hwe)

hweall <- hwe[which(hwe$TEST == "ALL"), ]

hwecases <- hwe[which(hwe$TEST == "AFF"), ]

###################################################

eigenvec_table <- read.table('AA_GWAS.eigenvec')

plot(data=eigenvec_table, V3~V4)


aa_gwas_data <- read.delim("GSA1604PCA.txt")


library(RColorBrewer)

#Defining the number of colors 

nb.cols <- 6

mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)




pcaplot <- ggplot(aa_gwas_data, aes(x=aa_gwas_data$PC1, y=aa_gwas_data$PC2, colour = aa_gwas_data$POP)) + geom_point() + scale_fill_manual(values = mycolors)
pcaplot

pcaplot + scale_color_manual(name="aa_gwas_data$SPOP",
                             labels = c("CHB","JPT","CHS","CDS","KHV","CEU","TSI","FIN","GBR","IBS","YRI","LWK","GWD","MSL","ESN","ASW","ACB","MXL","PUR","CLM","PEL","GIH","PJL","BEB","STU","ITU"),
                             values = c("CHB"="red",
                                         "JPT"="red",
                                          "CHS"="red",
                                          "CDS"="red",
                                          "KHV"="red"),
                                           "CEU"="blue",
                                            "TSI"="blue",
                                             "FIN"="blue",
                                             "GBR"="blue",
                                              "IBS"="blue",
                                              "YRI"="green",
                                              "LWK"="green",
                                              "GWD"="green",
                                               "MSL"="green",
                                               "ESN"="green",
                                               "ASW"="green",
                                               "ACB"="green",
                                              "MXL"="orange",
                                               "PUR"="orange",
                                               "CLM"="orange",
                                                "PEL"="orange",
                                                 "GIH"="brown",
                                                  "PJL"="black",
                                                  "BEB"="grey",
                                                  "STU"="pink",
                                                  "ITU"="violet")
                                           

pcaplot + expand_limits(x = c(-0.08,0.08), y = c(-0.08,0.08))

pcaplot1 <- ggplot(aa_gwas_data, aes(x=aa_gwas_data$PC2, y=aa_gwas_data$PC3, colour = aa_gwas_data$POP)) + geom_point() + scale_fill_manual(values = mycolors)
pcaplot1

pcaplot2 <- ggplot(aa_gwas_data, aes(x=aa_gwas_data$PC1, y=aa_gwas_data$PC3, colour = aa_gwas_data$POP)) + geom_point() + scale_fill_manual(values = mycolors)
pcaplot2


aa_gwas_maf_data <- read.delim("AA_GWAS_MAF_pca.txt")
aa_gwas_maf_data

pcaplot3 <- ggplot(aa_gwas_maf_data, aes(x=aa_gwas_maf_data$PC1, y=aa_gwas_maf_data$PC2, colour = aa_gwas_maf_data$Ph)) + geom_point() 
pcaplot3

pcaplot4 <- ggplot(aa_gwas_maf_data, aes(x=aa_gwas_maf_data$PC2, y=aa_gwas_maf_data$PC3, colour = aa_gwas_maf_data$Ph)) + geom_point() 
pcaplot4






aa_gwas_controls_data <- read.delim("gwas_pca_controls.txt")

pcplot1 <- ggplot(aa_gwas_controls_data, aes(x=aa_gwas_controls_data$PC1, y=aa_gwas_controls_data$PC2, colour = aa_gwas_controls_data$Ph)) + geom_point() 
pcplot1
pcaplot + expand_limits(x = c(-0.08,0.08), y = c(-0.08,0.08))

pcplot2 <- ggplot(aa_gwas_controls_data, aes(x=aa_gwas_controls_data$PC2, y=aa_gwas_controls_data$PC3, colour = aa_gwas_controls_data$Ph)) + geom_point()
pcplot2

pcplot3 <- ggplot(aa_gwas_controls_data, aes(x=aa_gwas_controls_data$PC1, y=aa_gwas_controls_data$PC3, colour = aa_gwas_controls_data$Ph)) + geom_point()
pcplot3
#############################################################



########################################################

options(scipen=100, digits=3)


eigen <- data.frame(read.table("AAGWAS5ADDMerge1.eigenvec", header = FALSE, skip = 0, sep = " " ))
row.names(eigen) <- eigen[,2]
eigen <- eigen[,3:ncol(eigen)]


eigen2 <- read.delim("MDS1-20.txt")



eigen1 <- princomp(eigen, cor = TRUE)

summary(eigen1)

summary(eigen)

Pov <- eigen$sdev^2/sum(eigen$sdev^2)

#Determine propertion of variance of each component

proportionofvariance <- ((apply(eigen, 1, sd)^2) / (sum(apply(eigen, 1, sd)^2)))*100

proportionofvariance

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table("1000gfullpop.eigenvec", header = FALSE, skip=0, sep = " ")
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("Principal Component ", c(1:10), sep = "")

# read in the PED data
PED <- read.table("20130606_g1k.ped", header = TRUE, skip = 0, sep = "\t")
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID == rownames(eigenvec)) == TRUE


# set colours
require("RColorBrewer")
display.brewer.all()

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

PED$Population1 <- factor(PED$Population, levels=c(
  "CDX","CHB","CHS","JPT","KHV",
  "BEB","GIH","ITU","PJL","STU"
))



table(PED$Population)
table(PED$Population1)

PED$Population

nlevels(PED$Population("ACB"));



col <- colorRampPalette(c(
  "purple","pink","red","brown","orange","green","yellow",
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "grey","grey","grey","grey","grey",
  "royalblue","royalblue","royalblue","royalblue","royalblue",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

col2 <- colorRampPalette(c(
   "blue","red","black","yellow","grey",
   "purple","pink","orange","brown","green"
))(length(unique(PED$Population1)))[factor(PED$Population1)]

col1 <- colorRampPalette(c(
  "yellow","yellow","yellow",
  "forestgreen","forestgreen",
  "grey","grey",
  "royalblue","royalblue",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)


par(mar=c(3,3,3,3), cex=1.2, cex.main=2, cex.axis=1.5, cex.lab=1.5, mfrow=c(1,2))

plot(project.pca[,1], project.pca[,2],
     type="n",
     main="A",
     adj=0.5,
     xlab="PC1",
     ylab="PC2",
     font=1,
     font.lab=1)
points(project.pca[,1], project.pca[,2], col=col2, pch=20, cex=1.5)
legend("bottomright",
       bty="n",
       cex=1,
       title="",
       c("cH-Dai","Han-Ch","SouthHan","Japan","Vietnam","Bengali","Gujrati","Telugu","Punjabi","Srilanka"),
       fill=c("blue","red","black","yellow","grey","purple","pink","orange","brown","green"))

plot(project.pca[,1], project.pca[,3],
     type="n",
     main="B",
     adj=1,
     xlab="PC1",
     ylab="PC3",
     font=1,
     font.lab=2)
points(project.pca[,1], project.pca[,3], col=col2, pch=20, cex=1.5)
legend("bottomright",
       bty="n",
       cex=1,
       title="",
       c("cH-Dai","Han-Ch","SouthHan","Japan","Vietnam","Bengali","Gujrati","Telugu","Punjabi","Srilanka"),
       fill=c("blue","red","black","yellow","grey","purple","pink","orange","brown","green"))

plot(project.pca[,2], project.pca[,3],
     type="n",
     main="C",
     adj=1,
     xlab="PC2",
     ylab="PC3",
     font=1,
     font.lab=2)
points(project.pca[,2], project.pca[,3], col=col2, pch=20, cex=1.5)
legend("bottomright",
       bty="n",
       cex=1,
       title="",
       c("cH-Dai","Han-Ch","SouthHan","Japan","Vietnam","Bengali","Gujrati","Telugu","Punjabi","Srilanka"),
       fill=c("blue","red","black","yellow","grey","purple","pink","orange","brown","green"))
###################################################

library(scatterplot3d)

# read in the vcf2eigen.pca.evec file
df <- read.table("1000 G SUB POP.txt")

# give column names
names(df) <- c("FID", "IID", "PC1", "PC2", "PC3","POP","SPOP")

# plot
with(df, scatterplot3d(PC1, PC2, PC3, color = as.numeric(POP), pch=10, main="PCA 1000 Genomes Exome Data (3D)")) 

# add legend
legend("topleft", pch=19, col=c("black", "green", "red", "blue"), legend=c("AFR", "EUR", "SAS", "EAS"))


# make plot
with(df, scatterplot3d(PC1, PC2, PC3, color = as.numeric(SPOP), pch=as.numeric(SPOP), main="PCA 1000 Genomes Exome Data (3D)", y.margin.add=1)) 

# add legend using x, y coordinates to position legend
legend(3.5, 4, pch=as.numeric(popIDs), col=as.numeric(popIDs), legend=popIDs)



###########################################################################################################

options(scipen=100, digits=3)

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table("1000gfullpop.eigenvec", header = FALSE, skip=0, sep = " ")
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste("Principal Component ", c(1:10), sep = "")

# read in the PED data
PED <- read.table("20130606_g1k.ped", header = TRUE, skip = 0, sep = "\t")
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID == rownames(eigenvec)) == TRUE


# set colours
require("RColorBrewer")

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

PED$Population1 <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "BEB","GIH","ITU","PJL","STU"))

PED$Population2 <- factor(PED$Population, levels=c(
  
  "CLM","MXL","PEL","PUR",
  "BEB","GIH","ITU","PJL","STU"))



col <- colorRampPalette(c(
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

col1 <- colorRampPalette(c(
  "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)


par(mar=c(2,2,2,2), cex=1.0, cex.main=2, cex.axis=1.75, cex.lab=1.75, mfrow=c(1,2))

plot(project.pca[,1], project.pca[,2],
     type="n",
     main="A",
     adj=0.5,
     xlab="First component",
     ylab="Second component",
     font=2,
     font.lab=2)
points(project.pca[,1], project.pca[,2], col=col, pch=20, cex=2.25)
legend("bottomright",
       bty="n",
       cex=1.0,
       title="",
       c("American","South Asian"),
       fill=c("forestgreen","black"))

plot(project.pca[,1], project.pca[,3],
     type="n",
     main="B",
     adj=0.5,
     xlab="First component",
     ylab="Third component",
     font=2,
     font.lab=2)
points(project.pca[,1], project.pca[,3], col=col, pch=20, cex=2.25)
legend("bottomright",
       bty="n",
       cex=1.0,
       title="",
       c("African","South Asian"),
       fill=c("yellow","black"))
plot(project.pca[,2], project.pca[,3],
     type="n",
     main="B",
     adj=0.5,
     xlab="First component",
     ylab="Third component",
     font=2,
     font.lab=2)
points(project.pca[,2], project.pca[,3], col=col, pch=20, cex=2.25)
legend("bottomright",
       bty="n",
       cex=1.0,
       title="",
       c("African","Hispanic","East-Asian","Caucasian","South Asian"),
       fill=c("yellow","forestgreen","grey","royalblue","black"))



############################################################

d <- read.delim("mds2")

plot(d$C1, d$C2, pch=20, cex=2, col = d$SOL+1)

mdsplot1 <- ggplot(d, aes(d$C1, y=d$C2, colour = d$POP)) + geom_point() 
mdsplot1

mdsplot2 <- ggplot(d, aes(d$C2, y=d$C3, colour = d$POP)) + geom_point() 
mdsplot2

mdsplot3 <- ggplot(d, aes(d$C1, y=d$C3, colour = d$POP)) + geom_point() 
mdsplot3

###########################################################################

library(forestplot)

cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA, 0.86,0.85,0.66,0.67,0.79,0.84,0.72,0.79), 
    lower = c(NA, 0.77,0.68,0.47,0.45,0.71,0.67,0.63,0.75),
    upper = c(NA, 0.96,1.07,0.93,1.00,0.87,1.07,0.81,0.84)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

tabletext<-cbind(
  c("Study", "KRIC", "HANIC", "MCIC", "Colombia","AAgwasHRS", "EAVyse", "HCVyse", 
  "Meta-Pval
I^2=-5.6,HetP=0.27,
FisherP=2.37E-13"),
  c("N(case/con)", "1710/3168","490/493","285/287","316/1366","1486/2786","4036/6959","1659/3398","9982/18457"),
  c("P-value", "6.42E-03","1.49E-01","1.36E-02","3.58E-02","6.64E-06","1.68E-01","3.23E-08","8.757E-16"),
  c("OR", "0.86","0.85","0.66","0.67","0.79","0.84","0.72",
    "0.79
    (0.75-0.84)"))

library(lattice)
trellis.device(device="windows", height = 25, width = 40, color=TRUE)



forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE,graphwidth = "auto", lineheight = "auto",
           is.summary=c(TRUE,rep(FALSE,7),TRUE),
           clip=c(0.45,1.10), 
           xlog=TRUE, 
           #txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.5)),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title="rs6705628 (TET3)" )


png("forestplot.png", width = 800, height = 600)
plot(...)
dev.off()


##############################################################################

library(forestplot)

cochrane_from_rmeta <- 
  structure(list(
    mean  = c(NA, 0.87,0.87,0.71,0.71,0.81,0.86,0.86,0.75), 
    lower = c(NA, 0.77,0.68,0.47,0.45,0.72,0.67,0.77,0.65),
    upper = c(NA, 0.96,1.07,0.94,1.00,0.90,1.07,0.96,0.84)),
    .Names = c("mean", "lower", "upper"), 
    row.names = c(NA, -8L), 
    class = "data.frame")

tabletext<-cbind(
  c("Study", "KRIC", "HANIC", "MCIC", "Columbia","AAgwasHRS", "EAVyse", "HCVyse", "Meta-Pval
I^2=2.6,HetP=0.22,
FisherP=1.46E-12"),c("N(case/con)", "1710/3168","490/493","285/287","316/1366","1486/2786","4036/6959","1659/3398","9982/18457"),c("P-value", "6.87E-03","1.47E-01","1.73E-02","3.43E-02","9.90E-05","1.62E-01","1.42E-08","9.08E-15"),c("OR", "0.86","0.84","0.66","0.66","0.81","0.86","0.84",
    "0.75
    (0.65-0.84)"))

library(lattice)
trellis.device(device="windows", height = 25, width = 40, color=TRUE)



forestplot(tabletext, 
           cochrane_from_rmeta,new_page = TRUE, graphwidth = "auto", lineheight = "auto",
           is.summary=c(TRUE,rep(FALSE,7),TRUE),
           clip=c(0.45,1.10), 
           xlog=TRUE, 
           #txt_gp = fpTxtGp(label = gpar(fontfamily = "", cex=0.5)),
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           title="rs6546883 (TET3)" )

png("forestplot102.png", width = 169, height = 100, res = 300)
print(pp)
dev.off


###################################################################

BiocManager::install("haploR", dependencies = TRUE)
library(haploR)
x <- queryHaploreg(file = "CHR1.txt", package = "haploR")

x1 <- queryHaploreg(file ="chr23.txt")

x6 <- queryHaploreg(file ="2CHR.txt")

x7 <- queryHaploreg(file ="3CHR.txt")

x8 <- queryHaploreg(file ="45CHR.txt")

x9 <- queryHaploreg(file ="6CHR.txt")

x10 <- queryHaploreg(file ="711CHR.txt")

x11 <- queryHaploreg(file ="1219CHR.txt")

x12 <- queryHaploreg(file ="6-1CHR.txt")

x13 <- queryHaploreg(file ="6-2CHR.txt")

x14 <- queryHaploreg(file ="6-3CHR.txt")

y1 <- write.csv(x1,"da.csv", row.names = TRUE)

y2 <- write.csv(x7,"da2.csv", row.names = TRUE)

y3 <- write.csv(x8,"da3.csv", row.names = TRUE)

y4 <- write.csv(x9,"da4.csv", row.names = TRUE)

y5 <- write.csv(x10,"da5.csv", row.names = TRUE)

y6 <- write.csv(x11,"da6.csv", row.names = TRUE)


######################################################

snptest <- read.delim("summaryfile.txt")

snptest3 <- read.delim("summaryfile1.txt")

snptest4 <- snptest3[snptest3$R2 > 0.2 & snptest3$IMD > 500000,]

snptest2 <- write.csv(snptest1,"data22.csv", row.names = TRUE)

snptest5 <- write.csv(snptest4,"data23.csv", row.names = TRUE)


##########################################################

BiocManager::install("RCurl")
library(RCurl)


myfile <- getURL('https://pubs.broadinstitute.org/mammals/haploreg/haploreg_v2.php',ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

class(myfile)

my.dat <- read.csv(textConnection(myfile), header = T)

head(my.dat)
