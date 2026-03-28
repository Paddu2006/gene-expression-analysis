
packages <- c("ggplot2", "RColorBrewer")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

set.seed(42)
n <- 200
genes <- paste0("Gene_", sprintf("%03d", 1:n))
ctrl <- matrix(rnbinom(n*3, mu=200, size=10), nrow=n, ncol=3)
trt  <- matrix(rnbinom(n*3, mu=200, size=10), nrow=n, ncol=3)
for(i in 1:20)  trt[i,] <- ctrl[i,] * 8
for(i in 21:40) trt[i,] <- ctrl[i,] * 0.1
mat <- cbind(ctrl, trt)
colnames(mat) <- c("C1","C2","C3","T1","T2","T3")
rownames(mat) <- genes
lmat <- log2(mat+1)
fc <- rowMeans(lmat[,4:6]) - rowMeans(lmat[,1:3])
pv <- apply(lmat,1,function(r) tryCatch(t.test(r[1:3],r[4:6])$p.value,error=function(e)1))
pa <- p.adjust(pv,"BH")
res <- data.frame(gene=genes,log2FC=round(fc,4),padj=round(pa,6))
res$status <- "NS"
res$status[pa<0.05 & fc>1]   <- "Up"
res$status[pa<0.05 & fc< -1] <- "Down"
cat("Up:", sum(res$status=="Up"),
    "Down:", sum(res$status=="Down"),
    "NS:", sum(res$status=="NS"), "\n")
library(ggplot2)
res$nlp <- -log10(pa+1e-10)
p <- ggplot(res,aes(log2FC,nlp,color=status))+
  geom_point(size=2,alpha=0.7)+
  scale_color_manual(values=c(Up="#E91E63",Down="#2196F3",NS="gray"))+
  geom_vline(xintercept=c(-1,1),linetype="dashed")+
  labs(title="Volcano Plot",x="Log2FC",y="-Log10 Padj")+
  theme_minimal()
ggsave("volcano_plot.png",p,width=8,height=6,dpi=150)
cat("Saved: volcano_plot.png\n")
pca <- prcomp(t(lmat),scale.=TRUE)
pdf <- data.frame(
  PC1=pca$x[,1],PC2=pca$x[,2],
  Sample=colnames(lmat),
  Condition=c("Control","Control","Control","Treated","Treated","Treated"))
p2 <- ggplot(pdf,aes(PC1,PC2,color=Condition,label=Sample))+
  geom_point(size=4)+geom_text(vjust=-1,size=3)+
  scale_color_manual(values=c(Control="#4CAF50",Treated="#E91E63"))+
  labs(title="PCA Plot")+theme_minimal()
ggsave("pca_plot.png",p2,width=8,height=6,dpi=150)
cat("Saved: pca_plot.png\n")
write.csv(res[order(res$padj),],"results.csv",row.names=FALSE)
cat("Saved: results.csv\n")
cat("ANALYSIS COMPLETE!\n")
