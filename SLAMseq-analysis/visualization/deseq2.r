rm(list=ls())
getwd()
setwd("/public/home/liunGroup/liangyan/project/zzZA/20250115_lzh-SLAMseq/grandslam")

library("ggplot2")

# 读取1.GrandSlam.tsv.gz
df <- read.table("1.GrandSlam.tsv.gz", header = T, sep = "\t")
# 提取包含conversions的列
conv <- df[,grep("Conversions", colnames(df))]
coverage <- df[,grep("[1-2].Coverage", colnames(df))]

conv_rate <- apply(conv,sum,MARGIN=2)/apply(coverage,sum,MARGIN=2)
rate = data.frame(conv_rate)
# rownames作为单独一列
rate$sample = gsub("^mapping\\.|\\.Conversions$", "", rownames(rate))
colnames(rate) <- c("Conversions","Sample")

# 去除灰色背景
ggplot(rate)+
  geom_bar(aes(x=Sample,y=Conversions),stat = "identity")+
  theme_bw()
ggsave("1.GrandSlam_conversions.pdf")

# DE analysis
library(grandR)
library(ggplot2)
library(patchwork)
library(dplyr)
library(DESeq2)
library(openxlsx)

dat <- ReadGRAND("1.GrandSlam.tsv",design=c("Condition","Replicate"))
dat$coldata["group"] <- c(
    rep("Rapa",4),
    rep("AP21967",4),
    rep(c("CD34-1d-0nM","CD34-1d-2nM","CD34-1d-8nM"),each=2),
    rep(c("HUDEP2-1d-0nM","HUDEP2-1d-2nM","HUDEP2-1d-10nM"),each=2)
    ) %>% factor()
# dat = LFC(dat,mode="total",normalize="total",contrasts=GetContrasts(dat,contrast=c("Condition")))
dat$data[["new"]] = dat$data[["count"]] * dat$data[["ntr"]]



de = function(columns,dat){
    # columns=c("C2.1","C2.2","C1.1","C1.2")
    df.total = dat$data$count[,columns]
    df.new = dat$data$new[,columns]
    coldata = data.frame(
        row.names = colnames(df.total),
        condition = c(rep("control",2),rep("treatment",2))
    )

    dds1 <- DESeqDataSetFromMatrix(countData = df.total, colData = coldata, design= ~condition)
    dds1 <- DESeq(dds1)
    res.total <- results(dds1, contrast = c("condition", "treatment", "control"))

    dds2 <- DESeqDataSetFromMatrix(countData = round(df.new), colData = coldata, design= ~condition)
    sizeFactors(dds2) <- dds1$sizeFactor
    dds2 <- estimateDispersions(dds2)
    dds2 <- nbinomWaldTest(dds2)
    res.new <- results(dds2, contrast = c("condition", "treatment", "control"))

    res.total = merge(data.frame(res.total),dat$gene.info[,1:2],by.x="row.names",by.y="Gene",alll=False)
    res.new = merge(data.frame(res.new),dat$gene.info[,1:2],by.x="row.names",by.y="Gene",alll=False)
    res.total <- res.total[order(res.total$padj, res.total$log2FoldChange, decreasing = c(FALSE, TRUE)),]
    res.new <- res.new[order(res.new$padj, res.new$log2FoldChange, decreasing = c(FALSE, TRUE)),]
    res = list("total" = res.total,"new" = res.new)
    return(res)
}

C1.C2 = de(c("C2.1","C2.2","C1.1","C1.2"),dat)
C3.C4 = de(c("C4.1","C4.2","C3.1","C3.2"),dat)

O2.O1 = de(c("O1.1","O1.2","O2.1","O2.2"),dat)
O3.O1 = de(c("O1.1","O1.2","O3.1","O3.2"),dat)

O5.O4 = de(c("O4.1","O4.2","O5.1","O5.2"),dat)
O6.O4 = de(c("O4.1","O4.2","O6.1","O6.2"),dat)

write.xlsx(C1.C2, file = "C1.C2.xlsx")
write.xlsx(C3.C4, file = "C3.C4.xlsx")
write.xlsx(O2.O1, file = "O2.O1.xlsx")
write.xlsx(O3.O1, file = "O3.O1.xlsx")
write.xlsx(O5.O4, file = "O5.O4.xlsx")
write.xlsx(O6.O4, file = "O6.O4.xlsx")