setwd('./htseq_counts')

library(edgeR)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

files_Control <- c("Control-1_counts.txt", "Control-2_counts.txt", "Control-3_counts.txt")
sample_names_Control <- gsub("_counts.txt", "", files_Control)

files_DrugA <- c("DrugA-1_counts.txt", "DrugA-2_counts.txt", "DrugA-3_counts.txt")
sample_names_DrugA <- gsub("_counts.txt", "", files_DrugA)

files_DrugB <- c("DrugB-1_counts.txt", "DrugB-2_counts.txt", "DrugB-3_counts.txt")
sample_names_DrugB <- gsub("_counts.txt", "", files_DrugB)

files_DrugC <- c("DrugC-1_counts.txt", "DrugC-2_counts.txt", "DrugC-3_counts.txt")
sample_names_DrugC <- gsub("_counts.txt", "", files_DrugC)

# 1. DEGs
## DrugA vs Control ## (Other comparisons are the same as this code.)

files_Control_DrugA <- c("Control-1_counts.txt", "Control-2_counts.txt", "Control-3_counts.txt", 
                "DrugA-1_counts.txt", "DrugA-2_counts.txt", "DrugA-3_counts.txt")
sample_names_Control_DrugA <- gsub("_counts.txt", "", files_Control_DrugA)

counts_list_Control_DrugA <- lapply(files_Control_DrugA, function(file) {
  df <- read.delim(file, header = FALSE, stringsAsFactors = FALSE)
  df <- df[!grepl("^__", df$V1), ]  
  rownames(df) <- df$V1         
  return(df[, 3, drop = FALSE])  
})

counts_matrix_Control_DrugA <- do.call(cbind, counts_list_Control_DrugA)
colnames(counts_matrix_Control_DrugA) <- sample_names_Control_DrugA

group_Control_DrugA <- factor(c(rep("Control", 3), rep("DrugA", 3)), levels = c("Control", "DrugA"))

dge_Control_DrugA <- DGEList(counts = counts_matrix_Control_DrugA, group = group_Control_DrugA)
dge_Control_DrugA <- dge_Control_DrugA[filterByExpr(dge_Control_DrugA), , keep.lib.sizes = FALSE]
dge_Control_DrugA <- calcNormFactors(dge_Control_DrugA)
design_Control_DrugA <- model.matrix(~ group_Control_DrugA)
dge_Control_DrugA <- estimateDisp(dge_Control_DrugA, design_Control_DrugA)
fit_Control_DrugA <- glmQLFit(dge_Control_DrugA, design_Control_DrugA)
qlf_Control_DrugA <- glmQLFTest(fit_Control_DrugA, coef = 2)
res_Control_DrugA <- topTags(qlf_Control_DrugA, n = Inf)$table

res_Control_DrugA$gene <- rownames(res_Control_DrugA)
res_Control_DrugA$threshold <- as.factor(abs(res_Control_DrugA$logFC) > 1 & res_Control_DrugA$FDR < 0.05)

ggplot(res_Control_DrugA, aes(x = logFC, y = -log10(FDR), color = threshold)) +
  geom_point(alpha = 0.6) +
  theme_minimal(base_size = 14) +
  # geom_text_repel(data = head(res[order(res$FDR), ], 10),
  #                 aes(label = gene), size = 3) +
  scale_color_manual(values = c("grey", "red"), guide = "none") +
  labs(title = "Control vs DrugA",
       x = "log2 FC",
       y = "-log10(FDR)") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme(plot.title = element_text(hjust = 0.5))

deg_genes_Control_DrugA <- rownames(res_Control_DrugA[res_Control_DrugA$threshold == TRUE, ])

logCPM_Control_DrugA <- cpm(dge_Control_DrugA, log = TRUE)

annotation_col_Control_DrugA <- data.frame(Group = group_Control_DrugA)
rownames(annotation_col_Control_DrugA) <- colnames(logCPM_Control_DrugA)

deg_logCPM_Control_DrugA <- logCPM_Control_DrugA[deg_genes_Control_DrugA, ]

deg_logCPM_Control_DrugA <- deg_logCPM_Control_DrugA[, order(group_Control_DrugA)]
annotation_col_Control_DrugA <- annotation_col_Control_DrugA[order(group_Control_DrugA), , drop = FALSE]

pheatmap(deg_logCPM_Control_DrugA,
         annotation_col = annotation_col_Control_DrugA,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         scale = "row",
         fontsize_col = 10,
         annotation_names_col = FALSE,
         border_color = NA)

# 2. GO Analysis
## DrugA vs Control ## (Other comparisons are the same as this code.)

up_gene_entrez_Control_DrugA <- bitr(up_gene_list_Control_DrugA,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

down_gene_entrez_Control_DrugA <- bitr(down_gene_list_Control_DrugA,
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Hs.eg.db)

up_gene_go_bp_Control_DrugA <- enrichGO(gene = up_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

up_gene_go_mf_Control_DrugA <- enrichGO(gene = up_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

up_gene_go_cc_Control_DrugA <- enrichGO(gene = up_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "CC",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

down_gene_go_bp_Control_DrugA <- enrichGO(gene = down_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "BP",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

down_gene_go_mf_Control_DrugA <- enrichGO(gene = down_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "MF",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

down_gene_go_cc_Control_DrugA <- enrichGO(gene = down_gene_entrez_Control_DrugA$ENTREZID,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "CC",
                  pAdjustMethod= "BH",
                  pvalueCutoff = 0.05,
                  readable     = TRUE)

up_gene_go_bp_Control_DrugA@result$ONTOLOGY <- "Biological Process"
up_gene_go_mf_Control_DrugA@result$ONTOLOGY <- "Molecular Function"
up_gene_go_cc_Control_DrugA@result$ONTOLOGY <- "Cellular Component"

down_gene_go_bp_Control_DrugA@result$ONTOLOGY <- "Biological Process"
down_gene_go_mf_Control_DrugA@result$ONTOLOGY <- "Molecular Function"
down_gene_go_cc_Control_DrugA@result$ONTOLOGY <- "Cellular Component"

up_gene_go_all_Control_DrugA <- bind_rows(up_gene_go_bp_Control_DrugA@result, up_gene_go_mf_Control_DrugA@result, up_gene_go_cc_Control_DrugA@result)
down_gene_go_all_Control_DrugA <- bind_rows(down_gene_go_bp_Control_DrugA@result, down_gene_go_mf_Control_DrugA@result, down_gene_go_cc_Control_DrugA@result)

up_gene_go_top_Control_DrugA <- up_gene_go_all_Control_DrugA %>%
  group_by(ONTOLOGY) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

down_gene_go_top_Control_DrugA <- down_gene_go_all_Control_DrugA %>%
  group_by(ONTOLOGY) %>%
  slice_min(p.adjust, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

ggplot(up_gene_go_top_Control_DrugA, aes(x = Description, y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width=0.8, color='black') +
  coord_flip(clip='off') +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 20, 10, 70),
  ) +
  scale_fill_manual(values = c("Biological Process" = "#E41A1C", "Molecular Function" = "#4DAF4A", "Cellular Component" = "#377EB8")) +
  labs(title = "GO Fuctional Analysis(Up gene)",
       x = NULL,
       y = "-log10 Adjusted P-value",
       fill = "Ontology") +
  theme(plot.title = element_text(hjust = 0.5))

ggplot(down_gene_go_top_Control_DrugA, aes(x = Description, y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity", width=0.8, color='black') +
  coord_flip(clip='off') +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(size = 12, color = "black", face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 20, 10, 100),
  ) +
  scale_fill_manual(values = c("Biological Process" = "#E41A1C", "Molecular Function" = "#4DAF4A", "Cellular Component" = "#377EB8")) +
  labs(title = "GO Fuctional Analysis(Down gene)",
       x = NULL,
       y = "-log10 Adjusted P-value",
       fill = "Ontology") +
  theme(plot.title = element_text(hjust = 0.5))
