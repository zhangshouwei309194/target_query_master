library(jsonlite)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(tidyr)
library(rstatix)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
odir <- args[2]

parse_json_with_nan <- function(file_path) {
  json_content <- readLines(file_path, warn = FALSE)
  json_content <- gsub("NaN", "null", json_content)
  # 同时处理其他可能的非标准数值表示
  json_content <- gsub("Infinity", "null", json_content)
  json_content <- gsub("-Infinity", "null", json_content)
  return(jsonlite::fromJSON(json_content))
}

## 数据清洗和预处理
depmap_gene_data <- parse_json_with_nan(input_file)
sample_level_expre_df <- depmap_gene_data$sample_level_expression$samples
sample_level_expre_df <- sample_level_expre_df[sample_level_expre_df$ModelType!="Organoid",]
sample_level_expre_df <- sample_level_expre_df[sample_level_expre_df$OncotreePrimaryDisease!="Non-Cancerous",]
sample_level_expre_df <- subset(sample_level_expre_df, !is.na(expression_value))
sample_level_expre_df$expression_value_log2 <- log2(sample_level_expre_df$expression_value+1)

### 保留OncotreeLineage中cell line number大于10的Lineage用于后续可视化
sample_level_expre_filtered_bylineage_df <- sample_level_expre_df %>%
  group_by(OncotreeLineage) %>%
  filter(n() >= 10) %>%
  ungroup()

sample_level_expre_filtered_bylineage_df <- sample_level_expre_filtered_bylineage_df %>%
  group_by(OncotreeLineage) %>%
  mutate(
    OncotreeLineage_count = n(),  # 计算每种组织类型的样本数
    OncotreeLineage_label = paste0(OncotreeLineage, "(n=", OncotreeLineage_count, ")")  # 创建标签
  ) %>%
  ungroup()

### OncotreePrimaryDisease中cell line number大于10的primary disease用于后续可视化
sample_level_expre_filtered_bydisease_df <- sample_level_expre_df %>%
  group_by(OncotreePrimaryDisease) %>%
  filter(n() >= 10) %>%
  ungroup()

sample_level_expre_filtered_bydisease_df <- sample_level_expre_filtered_bydisease_df %>%
  group_by(OncotreePrimaryDisease) %>%
  mutate(
    OncotreePrimaryDisease_count = n(),  # 计算每种组织类型的样本数
    OncotreePrimaryDisease_label = paste0(OncotreePrimaryDisease, "(n=", OncotreePrimaryDisease_count, ")")  # 创建标签
  ) %>%
  ungroup()

## 可视化不同OncotreeLineage的基因表达分布，同时可视化不同OncotreePrimaryDisease的基因表达分布
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
            legend.title = element_blank(),
            legend.position = "none",
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 1.5, 
                                        linetype = "solid"))
    p_Boxplot <- p_Boxplot + labs(x = xlabel,y = ylabel,title=title) 
    return (p_Boxplot)
}

order_level_func <- function(data_df,category_col) {
    order_level <- data_df %>%
        group_by(.data[[category_col]]) %>%
        summarise(median_value = median(expression_value, na.rm = TRUE)) %>%
        arrange(median_value) %>%
        pull(.data[[category_col]])
    return (order_level)
}

gene <- depmap_gene_data$gene_information$gene_name

options(repr.plot.width=20, repr.plot.height=10)

order_level <- order_level_func(sample_level_expre_filtered_bylineage_df,"OncotreeLineage_label")

custom_colors <- c('Cell Line'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_filtered_bylineage_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- boxplot_func(sample_level_expre_filtered_bylineage_df,"expression_value",'OncotreeLineage_label','ModelType',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreeLineage.samplesgt10.expression.tpm.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

custom_colors <- c('Cell Line'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_filtered_bylineage_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- boxplot_func(sample_level_expre_filtered_bylineage_df,"expression_value_log2",'OncotreeLineage_label','ModelType',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreeLineage.samplesgt10.expression.tpm.log.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

order_level <- order_level_func(sample_level_expre_filtered_bydisease_df,"OncotreePrimaryDisease_label")

custom_colors <- c('Cell Line'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_filtered_bydisease_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- boxplot_func(sample_level_expre_filtered_bydisease_df,"expression_value",'OncotreePrimaryDisease_label','ModelType',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreePrimaryDisease.samplesgt10.expression.tpm.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

custom_colors <- c('Cell Line'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_filtered_bydisease_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- boxplot_func(sample_level_expre_filtered_bydisease_df,"expression_value_log2",'OncotreePrimaryDisease_label','ModelType',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreePrimaryDisease.samplesgt10.expression.tpm.log.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

sample_level_expre_primaryormetastasis_df <- subset(sample_level_expre_df, !is.na(PrimaryOrMetastasis))
sample_level_expre_primaryormetastasis_df <- sample_level_expre_primaryormetastasis_df[sample_level_expre_primaryormetastasis_df$PrimaryOrMetastasis %in% c('Metastatic','Primary'),]
sample_level_expre_primaryormetastasis_filtered_df <- sample_level_expre_primaryormetastasis_df %>%
  group_by(OncotreePrimaryDisease) %>%
  filter(n() >= 10) %>%
  ungroup()

sample_level_expre_primaryormetastasis_filtered_df <- sample_level_expre_primaryormetastasis_filtered_df %>%
  group_by(OncotreePrimaryDisease) %>%
  mutate(
    OncotreePrimaryDisease_count = n(),  # 计算每种组织类型的样本数
    OncotreePrimaryDisease_label = paste0(OncotreePrimaryDisease, "(n=", OncotreePrimaryDisease_count, ")")  # 创建标签
  ) %>%
  ungroup()

# 主函数：生成基因表达箱线图（无填充版本）
generate_gene_expression_plot <- function(data,value_cols,custom_colors,title,xlabel,ylabel) {
  
  # 1. 数据预处理
    processed_data <- data %>%
      mutate(
        # 正确引用列名：使用!!sym()或直接使用列名
        !!value_cols := as.numeric(!!sym(value_cols))  # 动态列名引用
      ) %>%
      filter(!is.na(!!sym(value_cols)))  # 正确过滤缺失值

  # 2. 计算每个癌症类型的肿瘤样本中位值用于排序 (原发灶：Primary)
    tumor_medians <- processed_data %>%
      filter(PrimaryOrMetastasis == "Primary") %>%
      group_by(OncotreePrimaryDisease) %>%
      summarise(median_value = median(!!sym(value_cols), na.rm = TRUE)) %>%
      arrange(median_value)
  
  # 3. 确定癌症类型的顺序
  cancer_order <- tumor_medians$OncotreePrimaryDisease
  
  # 应用排序
  processed_data <- processed_data %>%
    mutate(OncotreePrimaryDisease = factor(OncotreePrimaryDisease, levels = cancer_order))
  
  # 4. 进行Wilcoxon检验（仅对同时有肿瘤和正常样本的癌症类型）
    stat_test_results <- data.frame()
    
    # 获取所有癌症类型
    cancer_types <- unique(processed_data$OncotreePrimaryDisease)
    
    for (cancer in cancer_types) {
      # 过滤特定癌症类型的数据
      cancer_data <- processed_data %>% 
        filter(OncotreePrimaryDisease == cancer)
      
      # 检查是否同时存在肿瘤和正常样本
      sample_types_present <- unique(cancer_data$PrimaryOrMetastasis)
      has_both_types <- all(c("Primary", "Metastatic") %in% sample_types_present)
      
      if (has_both_types) {
        # 构建公式字符串（只对value_cols使用变量）
        formula_str <- paste(value_cols, "~ PrimaryOrMetastasis")
        
        # 进行Wilcoxon检验
        stat_test <- tryCatch({
          cancer_data %>%
            wilcox_test(as.formula(formula_str)) %>%
            add_xy_position(x = "OncotreePrimaryDisease", dodge = 0.8) %>%
            mutate(OncotreePrimaryDisease = cancer)  # 添加癌症类型列
        }, error = function(e) {
          # 如果检验失败，返回空数据框
          data.frame()
        })
        
        if (nrow(stat_test) > 0) {
          stat_test_results <- bind_rows(stat_test_results, stat_test)
        }
      }
    }
  
  # 5. 创建显著性标记
  if (nrow(stat_test_results) > 0) {
    stat_test_results <- stat_test_results %>%
      mutate(
        significance = case_when(
          p < 0.001 ~ "***",
          p < 0.01 ~ "**",
          p < 0.05 ~ "*",
          TRUE ~ "NS"
        )
      )
  }
    
    legend_labels <- c(
    "Metastatic" = "Metastatic",
    "Primary" = "Primary"
    )

  # 6. 修改：生成新的x轴标签，同时显示Primary和Metastatic的样本量
label_data <- processed_data %>%
  group_by(OncotreePrimaryDisease, PrimaryOrMetastasis) %>%
  summarise(count = n(), .groups = 'drop') %>%
  # 使用pivot_wider将Primary和Metastatic展开为列
  pivot_wider(
    names_from = PrimaryOrMetastasis,
    values_from = count,
    values_fill = 0
  ) %>%
  # 确保存在Primary和Metastatic列（如果数据中缺少某些类型）
  {
    if (!"Primary" %in% colnames(.)) {
      mutate(., Primary = 0)
    } else .
  } %>%
  {
    if (!"Metastatic" %in% colnames(.)) {
      mutate(., Metastatic = 0)
    } else .
  } %>%
  # 生成标签
  mutate(
    label = case_when(
      Primary > 0 & Metastatic > 0 ~ 
        paste0(OncotreePrimaryDisease, " (Primary:", Primary, ", Metastatic:", Metastatic, ")"),
      Primary > 0 ~ 
        paste0(OncotreePrimaryDisease, " (Primary:", Primary, ")"),
      Metastatic > 0 ~ 
        paste0(OncotreePrimaryDisease, " (Metastatic:", Metastatic, ")"),
      TRUE ~ OncotreePrimaryDisease  # 如果没有样本，只显示疾病名称
    )
  )
  
  # 7. 绘图（无填充版本）
  p <- ggplot(processed_data, aes(x = OncotreePrimaryDisease, y = !!sym(value_cols))) +
    # 箱线图 - 无填充版本
    geom_boxplot(
      aes(color = PrimaryOrMetastasis),  # 使用颜色而非填充
      fill = "transparent",       # 设置填充为透明
      outlier.shape = NA, 
      alpha = 0.7,
      size = 0.8                  # 调整边框粗细
    ) +
    # 散点图（jitter points）
    geom_jitter(
      aes(color = PrimaryOrMetastasis), 
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = 1.5, 
      alpha = 0.6
    ) +
    # 自定义颜色
    scale_color_manual(
      values = custom_colors,
      labels = legend_labels
    ) +
    # x轴标签
    scale_x_discrete(
      labels = setNames(label_data$label, label_data$OncotreePrimaryDisease)
    ) +
    # 添加显著性标记
    {if (nrow(stat_test_results) > 0)
      geom_text(
        data = stat_test_results,
        aes(x = OncotreePrimaryDisease, 
            y = max(processed_data[[value_cols]], na.rm = TRUE) * 1.1, 
            label = significance),
        inherit.aes = FALSE, 
        size = 5, 
        fontface = "bold"
      )
    }+
    # 图表美化
    labs(
      title = title,
      x = xlabel,
      y = ylabel
    ) +
    theme_classic() + theme(
            plot.title = element_text(hjust = 0,vjust= 0.2, face = "bold", size = 14),
            plot.title.position = "panel",
            axis.title = element_text(size = 16),  # 轴标签字体大小
            axis.text = element_text(size = 12),   # 刻度标签字体大小
            axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1,face="bold",size=12,color="black"),
            legend.text = element_text(size = 16),  # 图例字体大小
            legend.title = element_blank(),
            legend.position = "top",
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 1.5, 
                                        linetype = "solid"))
  
  return(p)
}

options(repr.plot.width=20, repr.plot.height=12)

custom_colors <- c('Metastatic'='#ffc125','Primary'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_primaryormetastasis_filtered_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- generate_gene_expression_plot(sample_level_expre_primaryormetastasis_filtered_df,"expression_value",custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreePrimaryDisease.samplesgt10.expression.tpm.wilcoxon.test.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=12,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 12, dpi = 300)

custom_colors <- c('Metastatic'='#ffc125','Primary'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_primaryormetastasis_filtered_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- generate_gene_expression_plot(sample_level_expre_primaryormetastasis_filtered_df,"expression_value_log2",custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"DepMap.OncotreePrimaryDisease.samplesgt10.expression.tpm.log.wilcoxon.test.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=12,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 12, dpi = 300)














