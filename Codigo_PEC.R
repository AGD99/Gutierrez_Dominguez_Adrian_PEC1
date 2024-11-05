library(SummarizedExperiment)
#Creación de objeto SummarisedExperiment
data_values<-DataValues_S013[,-1]
data_values
assay_data<-as.matrix(data_values[, -c(1:5)])
t_assay_data<-t(assay_data)
sample_info<-data.frame(
  SUBJECTS = data_values$SUBJECTS,
  SURGERY = data_values$SURGERY,
  AGE = data_values$AGE,
  GENDER = data_values$GENDER,
  Group = data_values$Group
)
se <- SummarizedExperiment(assays = list(counts = t_assay_data), colData = sample_info)
metadata<-DataInfo_S013[,-1]
metadata(se)<-metadata
se
#Exploración de datos, Histogramas
GLU<-c("GLU_T0", "GLU_T2", "GLU_T4", "GLU_T5")

par(mfrow = c(3, 2))
for (x in GLU) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}
PESO <- c("PESO_T0", "PESO_T2", "PESO_T4", "PESO_T5")
COL <- c("COL_T0", "COL_T2", "COL_T4", "COL_T5")
LDL <- c("LDL_T0", "LDL_T2", "LDL_T4", "LDL_T5")
HDL <- c("HDL_T0", "HDL_T2", "HDL_T4", "HDL_T5")
LEU <- c("Leu_T0", "Leu_T2", "Leu_T4", "Leu_T5")
PRO <- c("Pro_T0", "Pro_T2", "Pro_T4", "Pro_T5")
for (x in PESO) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}

for (x in COL) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}
for (x in LDL) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}

for (x in HDL) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}
for (x in LEU) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}
for (x in PRO) {
  values <- assay(se)[x, ]
  hist(values, main = paste("Histograma de", x), xlab = x, breaks = 10)
}

#PCA
library(ggplot2)

datos_pca <- assay(se)
datos_pca_clean <- as.data.frame(datos_pca)
datos_pca_clean <- na.omit(datos_pca_clean)
pca<- scale(datos_pca_clean)
pca_result <- prcomp(pca)

pca_result$sdev
pca_result$rotation[,1]
head(pca_result$x)


par(mfrow = c(1, 1)
pca_scores <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
  geom_point(color = "blue", size = 2) +
  labs(title = "PCA primeros dos componentes", x = "PC1", y = "PC2") + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


save(se,file="ObjetoSE.Rda")
