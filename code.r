# Do not forget to set working directory to source file location
##Step 1. Preparing the data
if(!require(readxl)){
  install.packages("readxl")
  library(readxl)
}
dat <- data.frame(read_excel("data.xlsx", sheet=2))
head(dat)

##creating_the_venn_diagram
dat_Pulp_dif <- dat[which(rowMeans(!is.na(dat[,c(3,7)])) > 0.99), ]
dat_sv_dif <- dat[which(rowMeans(!is.na(dat[,c(4,8)])) > 0.99), ]
dat_Pulp_contr <- dat[which(rowMeans(!is.na(dat[,c(5,10)])) > 0.99), ]
dat_sv_contr <- dat[which(rowMeans(!is.na(dat[,c(6,9)])) > 0.99), ]


if(!require(VennDiagram)){
  install.packages("VennDiagram")
  library(VennDiagram)
}

if(!require(RColorBrewer)){
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

myCol1 <- brewer.pal(4, "Pastel2")

#this_code_will_save_the_venn_diagram_to_the_WD
venn.diagram(
  x = list(dat_Pulp_dif$Accession, dat_sv_dif$Accession, dat_Pulp_contr$Accession, dat_sv_contr$Accession),
  category.names = c("Pulp_dif" , "sv_dif" , "Pulp_contr", "sv_contr"),
  filename = '#venn_diagramm.png',
  resolution = 600,
  fill = myCol1,
  output=TRUE
)

if(!require(gplots)){
  install.packages("gplots")
  library(gplots)
}

#lists of group-specific names
v.table1 <- venn(list(dat_Pulp_dif$Name, dat_sv_dif$Name, dat_Pulp_contr$Name, dat_sv_contr$Name))
print(v.table1)
v.table2 <- venn(list(dat_Pulp_dif$Accession, dat_sv_dif$Accession, dat_Pulp_contr$Accession, dat_sv_contr$Accession))
print(v.table2)

##Step 2. preparing_the_data_for_quantitative_analysis
dat1 <- dat[which(rowMeans(!is.na(dat)) > 0.85), ]
rownames(dat1) <- dat1[,1]
dat1 <- dat1[,-1]
dat1 <- dat1[,-1]
str(dat1)
mean(complete.cases(dat1))
colSums(is.na(dat1))

#imputation_of_missed_values_by_knn
if(!require(impute)){
  install.packages("impute")
  library(impute)
}

tdat <- t(dat1)
dat_knn1 <- impute.knn(tdat, k = 5)
dat_knn <- t(dat_knn1$data)
colSums(is.na(dat_knn))
mean(complete.cases(dat_knn))

#opening_legend_with_the_factors
fact <- data.frame(read_excel("legend.xlsx"))

#rename samples
colnames(dat_knn) <- fact$Sample
head(dat_knn)

#factor_matrix
rownames(fact) <- fact[,1]
fact <- fact[,-1]

fact$Cell_type <- as.factor(fact$Cell_type)
fact$Donor <- as.factor(fact$Donor)
fact$Differentiation <- as.factor(fact$Differentiation)
fact$Differentiation
fact$comb <- as.factor(fact$comb)
str(fact)

#looking_to_the_data
pal <- brewer.pal(n = 9, name = "Set1")
cols <- pal[fact$Differentiation]
boxplot(dat_knn, outline = FALSE, col = cols, main = "Raw data")
legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
colSums(dat_knn)
#log transformation
dat_log <- log2(dat_knn+1)
head(dat_log)
boxplot(dat_log, outline = FALSE, col = cols, main = "Log-transformed data")
legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
#Quantile normalization
if(!require(limma)){
  install.packages("limma")
  library(limma)
}
dat_norm <- normalizeQuantiles(dat_log)
head(dat_norm)
boxplot(dat_norm, col = cols, main = "Normalized data")
legend("topright", levels(fact$Differentiation), fill = pal, bty = "n", xpd = T)
mean(complete.cases(dat_norm))

#MA-plot
maplot <- function(X1, X2, pch = 21, main = "MA-plot", xlab = "Average log-expression", ylab = "Expression log-ratio", lpars = list(col = "blue", lwd = 2), ...){
  X <- (rowMeans(X2) + rowMeans(X1)) / 2
  Y <- rowMeans(X2) - rowMeans(X1)
  scatter.smooth(x = X, y = Y,
                 main = main, pch = pch,
                 xlab = xlab, ylab = ylab,
                 lpars = lpars, ...)
  abline(h = c(-1, 0, 1), lty = c(2, 1, 2))
}

maplot(dat_log[, c(1,3,5,8)], dat_log[, c(2,4,6,7)], main = "Log-expression data")
maplot(dat_norm[, c(1,3,5,8)], dat_norm[, c(2,4,6,7)], main = "Normalized data")


## Step 3. Clusterizations
t_dat_norm <- t(dat_norm)
# nMDS
if(!require(vegan)){
  install.packages("vegan")
  library(vegan)
}
dat_ord <- metaMDS(t_dat_norm,
                   distance = "euclidean",
                   autotransform = FALSE)
dat_ord$stress

ord_ggplo <- data.frame(fact, scores(dat_ord, display = "sites"))
head(ord_ggplo, 2)

ord_ggplo_sp <- data.frame(scores(dat_ord, display = "species"))
ord_ggplo_sp$Species <- rownames(ord_ggplo_sp)
head(ord_ggplo_sp, 2)

if(!require(ggplot2)){
  install.packages("vegan")
  library(ggplot2)
}
#tiff('nMDS.tiff', units="in", width=12, height=8, res=600, compression = 'lzw')
ggplot() +
  geom_point(data = ord_ggplo, 
             aes(x = NMDS1, y = NMDS2, colour = Cell_type, 
                 shape = Differentiation), size = 4)
#dev.off()


#PCA
if(!require(mixOmics)){
  install.packages("mixOmics")
  library(mixOmics)
}


dat_pca <- pca(t(dat_norm), ncomp = 8, center = TRUE)
dat_pca
plot(dat_pca)
dat_pca <- pca(t(dat_norm), ncomp = 4, center = TRUE)

#tiff('PCA_cell_type_2D.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')
plotIndiv(dat_pca, comp = c(1, 3), ind.names = F, 
          group = fact$Cell_type, legend = TRUE, ellipse = TRUE,
          title = 'PCA')
#dev.off()

if(!require(scatterplot3d)){
  install.packages("scatterplot3d")
  library(scatterplot3d)
}
colors <- c("#E69F00", "#56B4E9") 
colors <- colors[as.numeric(fact$Cell_type)] 
#tiff('PCA_cell_typeEd.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')
s3d <-scatterplot3d(dat_pca$X[,1], dat_pca$X[, 2],dat_pca$X[, 3],xlab="PC1(36%)",ylab="PC2(17%)", zlab="PC3(15%)", pch = 16, color=colors)
legend("right", legend = levels(fact$Cell_type),
       col =  c("#E69F00", "#56B4E9"), pch = 16)
#dev.off()


#PLSDA_cell_type
ordination.optimum.splsda <- splsda(t_dat_norm, fact$Cell_type, ncomp = 3, keepX = c(15,15,15))
selectVar(ordination.optimum.splsda, comp=1)
selectVar(ordination.optimum.splsda, comp=2)
selectVar(ordination.optimum.splsda, comp=3)

#tiff('PLSDA_cell.tiff', units="in", width=12, height=8, res=600, compression = 'lzw')
layout(matrix(c(1, 2, 3, 3, 3, 3), 2, 3))
plotLoadings(ordination.optimum.splsda, comp = 1, size.name = 1, size.title = 1.2, title = "Loadings\n on 1st component", contrib = "max", legend = FALSE, col.ties="black", ndisplay = 15)
plotLoadings(ordination.optimum.splsda, comp = 2, size.name = 1, size.title = 1.2, title = "Loadings\n on 2nd component", contrib = "max",ndisplay = 15,  legend = FALSE, col.ties="black")
plotIndiv(ordination.optimum.splsda, ind.names = F, ellipse = T, style = "graphics", abline = TRUE, cex = 1.5, size.axis = 1.2, size.xlabel = 1.5, size.ylabel = 1.5, title = "PLS-DA ordination of differnt cell types", size.title = 1.5, legend=TRUE)
#dev.off()
layout(1,1)



## Step 4. Limma_differential_protein_expression
dim(dat_norm)
dim(fact)
head(dat_norm)

fact$Cell_type
X <- model.matrix(~ fact$Cell_type)
X

fit <- lmFit(dat_norm, design = X, method = "robust", maxit = 10000)

# Empirical Bayes statistics
efit <- eBayes(fit)

# Dif_expr_table
topTable(efit, coef = 2)
numGenes <- nrow(dat_norm)
full_list_efit <- topTable(efit, number = numGenes)
#write.csv(full_list_efit,'Dif_expr.csv')
head(full_list_efit)

#Vulcanoplot representation of differential expression
if(!require(EnhancedVolcano)){
  install.packages("EnhancedVolcano")
  library(EnhancedVolcano)
}

#tiff('Vulcano.tiff', units="in", width=11, height=8, res=600, compression = 'lzw')
EnhancedVolcano(full_list_efit,
                lab = rownames(full_list_efit),
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1.5,
                xlim = c(-5, 10),
                ylim = c(0, 2),
                title ="DPSCs versus PDLSCs",
                colAlpha = 1)
#dev.off()

