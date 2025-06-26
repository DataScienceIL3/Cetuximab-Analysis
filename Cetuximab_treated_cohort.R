# 1.Cargar los paquetes necesarios
library(Seurat)
library(SingleCellExperiment)
library(Matrix.utils)
library(data.table)
library(base)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(DGEobj.utils)
library(plotly)
library(pheatmap)
library(corrplot)
library(viridis)
library(corto)
library(DESeq2)
library(SeuratObject)
library(openxlsx)
library(dplyr)
library(EPIC)
library(CMSclassifier)
library(janitor)
library(ggpubr)
library(ggplot2)
library(survival)
library(survminer)
library(GEOquery)
library(biomaRt)
library('org.Hs.eg.db')

# 2.Descargar los datos de expresión de GEO 
library(GEOquery)
data <- getGEO("gse5851", GSEMatrix=T)
expr_gse5851<-exprs(data[[1]])

#2. Leer i preparar los datos clínicos

clinical<-read.csv("Cetu treated patients_Clinical info.csv",
                   header=T,
                   sep=";")

#Sólo hasta la fila 80 porque es de los pacientes que hay información
clinical<-clinical[1:80,]

#Asignamos como rownames los identificadores de las muestras, para que encajen con la expresión
row.names(clinical)<-clinical$GSE5851

#3. Mapear las sondas de Affymetrix a los nombres de los genes

affyprobe<-row.names(expr_gse5851)

#Conectar con biomaRT para obtener la equivalencia entre sondas y nombres de genes
mart<-useEnsembl(biomart="ensembl",
                 dataset="hsapiens_gene_ensembl")

probetogene<-getBM(attributes=c("hgnc_symbol","affy_hg_u133a_2"),
                   filters="affy_hg_u133a_2",
                   values=affyprobe,
                   mart=mart)

#Eliminar sondas duplicadas
probetogene<-probetogene[!duplicated(probetogene$affy_hg_u133a_2),]

#Asignamos nombres de las sondas como rownames
row.names(probetogene)<-probetogene$affy_hg_u133a_2

#4. Fusionar los datos de expresión con el nombre de los genes

expr<-merge(probetogene,expr_gse5851,by="row.names")

#5. Colapsar datos de múltiples sondas por gen

data_gene_collapse<-aggregate(expr,by=list(expr$hgnc_symbol),FUN="max")

data_gene_collapse <- data_gene_collapse[-which(data_gene_collapse$hgnc_symbol == ""), ]
row.names(data_gene_collapse)<-data_gene_collapse$hgnc_symbol
data_gene_collapse<-data_gene_collapse[,-c(1:4)]

#5. Crear una variable de respuesta simplificada (R/NR/UTD)

clinical$DCG<-ifelse(clinical$Best.Clinical.Response.Assessment %in% c("SD","CR","PR"),"R",
                     ifelse(clinical$Best.Clinical.Response.Assessment=="UTD","UTD","NR"))

#6. Filtrar por pacientes KRAS WT (son los que pueden responder a cetuximab)

clinical_wt<-subset(clinical,KRAS.Mutation=="WT")

#7. Crear una tabla de respuesta

library(janitor)
tabyl(clinical_wt,DCG,Best.Clinical.Response.Assessment)



#####AHORA HACEMOS UN ANÁLISIS DE GENES DIFERENCIALMENTE EXPRESADOS (dge) ENTRE RESPONDERS Y NO RESPONDERS

library(limma)

#1. Preparación de los datos

#Descartamos los pacientes UTD (sin info de respuesta)
clinical_clean <- subset(clinical_wt, DCG %in% c("R", "NR"))

#Passem a logaritmo los valores de expresión
expr_log2 <- log2(data_gene_collapse + 1)

#Submatriu de expresión con las muestras seleccionadas
expr_matrix <- expr_log2[, clinical_clean$GSE5851]

#Comprovar que las columnas de la matriz y las filas del metadata coincidan
all(colnames(expr_matrix) == clinical_clean$GSE5851)  # tiene que salir TRUE


#2. Crea el diseño del modelo y ajusta el análisis con limma
#Define la condición de respuesta
group <- factor(clinical_clean$DCG, levels = c("NR", "R"))  # 'NR' como referencia

#Crea el modelo del diseño
design <- model.matrix(~ group)

#Ajusta el modelo con limma
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

#Extrae los resultados
res <- topTable(fit, coef = "groupR", number = Inf)

#3. Explora los resultados

#Genes diferenciales con FDR < 0.05 y loFC > 1.5
sig_genes <- res[which(res$P.Value < 0.05 & abs(res$logFC) > 1.5), ]

#Ver los primeros 10 genes diferenciales (OPCIONAL)
head(sig_genes, 10)


####REPRESENTACIÓN DE LOS DATOS
#1. VOLCANO PLOT
library(ggplot2)
library(dplyr)
library(ggrepel)

res$threshold <- with(res, P.Value < 0.05 & abs(logFC) > 1.5)

res <- res %>%
  mutate(Significance = ifelse(P.Value < 0.05 & abs(logFC) > 1.5,
                               "Significant",
                               "Not significant"))

##Volcano plot sin nombre de los genes
ggplot(res, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme(legend.position = "bottom")

##Volcano plot con nombre de los genes
sig_genes <- res %>% filter(Significance == "Significant")

ggplot(res, aes(x = logFC, y = -log10(P.Value), color = Significance)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("Not significant" = "grey", "Significant" = "red")) +
  geom_text_repel(data = sig_genes, aes(label = rownames(sig_genes)),
                  size = 3, max.overlaps = 15) +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)",
       color = "Significance") +
  theme(legend.position = "bottom")


#2. HEATMAP
library(pheatmap)

# Filtra genes significativos (usamos el p-valor) ojo!!! Hemos cambiado a 1.5 el logFC!!!

sig_genes <- res[res$P.Value < 0.05 & abs(res$logFC) > 1.5, ]
genes_sig <- rownames(sig_genes)

expr_sig <- expr_log2[genes_sig, clinical_clean$GSE5851]


pheatmap(expr_sig,
         scale = "row",
         show_rownames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "Heatmap - Gens diferencials")





######EXTRACCIÓN DE SIGNATURAS REPRESENTATIVAS DE CADA GRUPO
#1. Filtra los genes significativos
sig_genes <- res[res$P.Val < 0.05 & abs(res$logFC) > 1, ]

#2. Signatura responders (upregulated in responders)
genes_responder <- rownames(sig_genes[sig_genes$logFC > 0, ])

#3. Signatura no responders (upregulated in no responders)
genes_nonresponder <- rownames(sig_genes[sig_genes$logFC < 0, ])

#4. Ordenamos los genes per logFC decreciente
genes_responder <- rownames(sig_genes[sig_genes$logFC > 0, ])
genes_responder <- genes_responder[order(sig_genes[genes_responder, "logFC"], decreasing = TRUE)]

genes_nonresponder <- rownames(sig_genes[sig_genes$logFC < 0, ])
genes_nonresponder <- genes_nonresponder[order(sig_genes[genes_nonresponder, "logFC"], decreasing = FALSE)]  # més negatiu = més diferencial

#5. Descargar las signaturas
write.table(genes_responder, "signature_responders.CSV", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(genes_nonresponder, "signature_nonresponders.CSV", row.names = FALSE, col.names = FALSE, quote = FALSE)


#######VIAS Y FUNCIONES BIOLOGICAS RELACIONADAS CON LAS SIGNATURAS???
###Enrichment analysis (GO, KEGG, Hallmarks)
if (!requireNamespace("clusterProfiler", quietly=TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly=TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("enrichplot", quietly=TRUE)) BiocManager::install("enrichplot")
if (!requireNamespace("msigdbr", quietly=TRUE)) install.packages("msigdbr")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(msigdbr)

#1. Converteix de noms de gens a ENTREZID:
entrez_responder <- bitr(genes_responder, fromType = "SYMBOL", 
                         toType = "ENTREZID", OrgDb = org.Hs.eg.db)

entrez_nonresponder <- bitr(genes_nonresponder, fromType = "SYMBOL", 
                         toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#2.  Enrichment GO
ego <- enrichGO(gene = entrez_responder$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP", # canvia a "MF" o "CC" si vols
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05,
                readable = TRUE)

barplot(ego, showCategory = 15)


#3. GSEA clàssic (amb tots els gens ordenats)

gene_list <- res$logFC
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

#GSEA GO
gse_go <- gseGO(geneList = gene_list,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                keyType = "SYMBOL",
                minGSSize = 10,
                maxGSSize = 500,
                pvalueCutoff = 0.05,
                verbose = FALSE)

gseaplot2(gse_go, geneSetID = 1)  # pots canviar l'índex





























