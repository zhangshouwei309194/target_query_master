## 用于进行GTEx数据可视化
library(jsonlite)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
odir <- args[2]

## function
boxplot_func <- function(data_df,value_cols,category_cols,type_col,order_level,custom_colors,title,xlabel,ylabel) {
    # 确保分组列是因子，并按 order_level 排序（即使某些分组不存在于数据中）
    data_df[[category_cols]] <- factor(
        data_df[[category_cols]],
        levels = order_level
    )
    # 提取数据并去除缺失值（保留因子顺序）
    subdf <- data_df[, c(value_cols, category_cols,type_col)] %>% 
        na.omit() %>%
        mutate(across(all_of(category_cols), ~ factor(., levels = order_level)))  # 保持因子顺序
    
    p_Boxplot <- ggplot(subdf, 
              aes(x = !!sym(category_cols), 
                  y = !!sym(value_cols), 
                  fill = !!sym(type_col))) # 使用原始列名进行填充
    # 箱线图
    p_Boxplot <- p_Boxplot + geom_boxplot(
      width = 0.5,
      alpha = 0.7,
      outliers = TRUE,
      color = "black")
    
    p_Boxplot <- p_Boxplot + scale_fill_manual(values = custom_colors)
    p_Boxplot <- p_Boxplot + theme_classic() + theme(
            plot.title = element_text(hjust = 0,vjust= 0.2, face = "bold", size = 14),
            plot.title.position = "panel",
            axis.title = element_text(size = 16),  # 轴标签字体大小
            axis.text = element_text(size = 12),   # 刻度标签字体大小
            axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1,face="bold",size=12,color="black"),
            legend.text = element_text(size = 16),  # 图例字体大小
            legend.position = "none",
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 1.5, 
                                        linetype = "solid"))
    p_Boxplot <- p_Boxplot + labs(x = xlabel,y = ylabel,title=title) 
    return (p_Boxplot)
}

baplot_func <- function(data_df,value_cols,category_cols,order_level,custom_colors,title,xlabel,ylabel) {
    # 确保分组列是因子，并按 order_level 排序（即使某些分组不存在于数据中）
    data_df[[category_cols]] <- factor(
        data_df[[category_cols]],
        levels = order_level
    )
    # 提取数据并去除缺失值（保留因子顺序）
    subdf <- data_df[, c(value_cols, category_cols)] %>% 
        na.omit() %>%
        mutate(across(all_of(category_cols), ~ factor(., levels = order_level)))  # 保持因子顺序
    
    p_Barplot <- ggplot(subdf, 
              aes(x = !!sym(category_cols), 
                  y = !!sym(value_cols)))
    # 箱线图
    p_Barplot <- p_Barplot + geom_bar(stat = "identity", width = 0.4,fill=custom_colors,color="black")
    p_Barplot <- p_Barplot + theme_classic() + theme(
            plot.title = element_text(hjust = 0,vjust= 0.2, face = "bold", size = 14),
            plot.title.position = "panel",
            axis.title = element_text(size = 16),  # 轴标签字体大小
            axis.text = element_text(size = 12),   # 刻度标签字体大小
            axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1,face="bold",size=12,color="black"),
            legend.text = element_text(size = 16),  # 图例字体大小
            legend.position = "none",
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 1.5, 
                                        linetype = "solid"))
    p_Barplot <- p_Barplot + labs(x = xlabel,y = ylabel,title=title) 
    return (p_Barplot)
}

order_level_func <- function(data_df,category_col) {
    order_level <- data_df %>%
        group_by(.data[[category_col]]) %>%
        summarise(median_value = median(expression_value, na.rm = TRUE)) %>%
        arrange(desc(median_value)) %>%
        pull(.data[[category_col]])
    return (order_level)
}

## load json data
gtex_gene_data <- fromJSON(input_file)

sample_level_expre_df <- gtex_gene_data$sample_level_expression$samples

sample_level_expre_df <- sample_level_expre_df %>%
  group_by(tissue_type) %>%
  mutate(
    tissue_type_count = n(),  # 计算每种组织类型的样本数
    tissue_type_label = paste0(tissue_type, "(n=", tissue_type_count, ")")  # 创建标签
  ) %>%
  group_by(tissue_subtype) %>%
  mutate(
    tissue_subtype_count = n(),  # 计算每种组织子类型的样本数
    tissue_subtype_label = paste0(tissue_subtype, "(n=", tissue_subtype_count, ")")  # 创建标签
  ) %>%
  ungroup()

sample_level_expre_df$tissue_type_count <- NULL
sample_level_expre_df$tissue_subtype_count <- NULL

sample_level_expre_df$SAMPLETYPE = 'Normal Tissue'

sample_level_expre_df$expression_value_log2 <- log2(sample_level_expre_df$expression_value+1)

tissue_level_expre_df <- gtex_gene_data$tissue_level_expression$tissue_data

tissue_level_expre_df$expression_value_log2 <- log2(tissue_level_expre_df$expression_value+1)

tau_score <- gtex_gene_data$tissue_specificity_analysis$tau_score
specificity_type <- gtex_gene_data$tissue_specificity_analysis$specificity_type
specific_tissues <- gtex_gene_data$tissue_specificity_analysis$specific_tissues

gene <- gtex_gene_data$gene_information$gene_name

## visual expression in boxplot
order_level <- order_level_func(sample_level_expre_df,"tissue_type_label")
custom_colors <- c('Normal Tissue'='#00aa36')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- boxplot_func(sample_level_expre_df,"expression_value",'tissue_type_label','SAMPLETYPE',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"GTEx.expression.tpm.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

order_level <- order_level_func(sample_level_expre_df,"tissue_type_label")
custom_colors <- c('Normal Tissue'='#00aa36')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- boxplot_func(sample_level_expre_df,"expression_value_log2",'tissue_type_label','SAMPLETYPE',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"GTEx.expression.tpm.log.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

## visual expression in borplot
order_level <- order_level_func(tissue_level_expre_df,"tissue_name")
custom_colors = "#00aa36"
title=paste(gene," median expression of tissue type",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- baplot_func(tissue_level_expre_df,"expression_value",'tissue_name',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"GTEx.expression.median.tpm.bytissue.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

order_level <- order_level_func(tissue_level_expre_df,"tissue_name")
custom_colors = "#00aa36"
title=paste(gene," median expression of tissue type",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- baplot_func(tissue_level_expre_df,"expression_value_log2",'tissue_name',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"GTEx.expression.median.tpm.log.bytissue.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)