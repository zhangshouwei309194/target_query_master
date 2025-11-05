library(jsonlite)
library(dplyr)  # 用于数据操作
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

tcga_gene_data <- parse_json_with_nan(input_file)

sample_level_expre_df <- tcga_gene_data$sample_level_expression$samples

sample_level_expre_df <- subset(sample_level_expre_df, !is.na(expression_value))

sample_level_expre_df <- sample_level_expre_df[sample_level_expre_df$sample_type %in% c('Primary Tumor','Solid Tissue Normal'),]

generate_disease_stats <- function(data) {
  # 数据清洗
  clean_data <- data %>%
    mutate(sample_type = case_when(
      sample_type == "Solid Tissue Normal" ~ "matched normal",
      sample_type == "Primary Tumor" ~ "Primary Tumor",
      TRUE ~ sample_type  # 保留其他样本类型
    ))
  
  # 统计
  stats <- clean_data %>%
    group_by(primary_disease, cancer_type, sample_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    arrange(primary_disease, sample_type)
  
  # 格式化输出并创建标签
  stats_formatted <- stats %>%
    mutate(
      primary_disease_label = ifelse(
        sample_type == "matched normal",
        paste0(primary_disease, "(", cancer_type, " ", sample_type, ",", count, ")"),
        paste0(primary_disease, "(", cancer_type, ",", count, ")")
      )
    )
  
  # 将标签信息合并回原始数据
  result_data <- clean_data %>%
    left_join(stats_formatted %>% 
                select(primary_disease, cancer_type, sample_type, primary_disease_label),
              by = c("primary_disease", "cancer_type", "sample_type"))
  
  return(result_data)
}

sample_level_expre_df <- generate_disease_stats(sample_level_expre_df)
sample_level_expre_df <- sample_level_expre_df %>%
    mutate(sample_type = case_when(
      sample_type == "matched normal" ~ "Solid Tissue Normal",
      sample_type == "Primary Tumor" ~ "Primary Tumor",
      TRUE ~ sample_type
    ))

sample_level_expre_df$expression_value_log2 <- log2(sample_level_expre_df$expression_value+1)
sample_type_level_expre_df <- tcga_gene_data$cancer_sample_type_level_expression$cancer_sample_type_data
sample_type_level_expre_df$expression_value_log2 <- log2(sample_type_level_expre_df$expression_value+1)

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
            legend.position = "top",
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

## style1: Primary Tumor和Solid Tumor Normal基因表达分布图 （图例按照表达量从低到高，左边为Normal,右边为Primary）
gene <- tcga_gene_data$gene_information$gene_name
sample_level_expre_primary_df <- sample_level_expre_df[sample_level_expre_df$sample_type=="Primary Tumor",]
sample_level_expre_normal_df <- sample_level_expre_df[sample_level_expre_df$sample_type=="Solid Tissue Normal",]
order_level1 <- order_level_func(sample_level_expre_primary_df,"primary_disease_label")
order_level2 <- order_level_func(sample_level_expre_normal_df,"primary_disease_label")
order_level <- c(order_level2,order_level1)

custom_colors <- c('Solid Tissue Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- boxplot_func(sample_level_expre_df,"expression_value",'primary_disease_label','sample_type',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.sampletype.expression.tpm.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

custom_colors <- c('Solid Tissue Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- boxplot_func(sample_level_expre_df,"expression_value_log2",'primary_disease_label','sample_type',order_level,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.sampletype.expression.tpm.log.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

## style2: Primary Tumor和Solid Tumor Normal基因表达差异分析和表达分布图 （图例按照Primary表达量从低到高）
# 生成基因表达箱线图（无填充版本）
generate_gene_expression_plot <- function(data, value_cols, custom_colors, title, xlabel, ylabel) {
  # 1. 数据预处理
  processed_data <- data %>%
    mutate(
      # 将Solid Tissue Normal转换为Matched Normal
      sample_type = ifelse(sample_type == "Solid Tissue Normal", "Matched Normal", sample_type),
      # 正确引用列名
      !!value_cols := as.numeric(!!sym(value_cols))
    ) %>%
    filter(!is.na(!!sym(value_cols)))
  
  # 2. 计算每个癌症类型的肿瘤样本中位值用于排序
  tumor_medians <- processed_data %>%
    filter(sample_type == "Primary Tumor") %>%
    group_by(cancer_type) %>%
    summarise(median_value = median(!!sym(value_cols), na.rm = TRUE)) %>%
    arrange(median_value)
  
  # 3. 确定癌症类型的顺序
  cancer_order <- tumor_medians$cancer_type
  processed_data <- processed_data %>%
    mutate(cancer_type = factor(cancer_type, levels = cancer_order))
  
  # 4. 进行Wilcoxon检验（仅对同时有肿瘤和正常样本的癌症类型）
  stat_test_results <- data.frame()
  
  # 获取所有癌症类型
  cancer_types <- unique(processed_data$cancer_type)
  
  for (cancer in cancer_types) {
    # 过滤特定癌症类型的数据
    cancer_data <- processed_data %>% 
      filter(cancer_type == cancer)
    
    # 检查是否同时存在肿瘤和正常样本
    sample_types_present <- unique(cancer_data$sample_type)
    has_both_types <- all(c("Primary Tumor", "Matched Normal") %in% sample_types_present)
    
    if (has_both_types) {
      # 构建公式字符串
      formula_str <- paste(value_cols, "~ sample_type")
      
      # 进行Wilcoxon检验
      stat_test <- tryCatch({
        cancer_data %>%
          wilcox_test(as.formula(formula_str)) %>%
          add_xy_position(x = "cancer_type", dodge = 0.8) %>%
          mutate(cancer_type = cancer)
      }, error = function(e) {
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
  
  # 6. 修改：生成新的x轴标签，同时显示Tumor和Normal的样本量
  label_data <- processed_data %>%
    group_by(cancer_type, primary_disease, sample_type) %>%
    summarise(count = n(), .groups = 'drop') %>%
    # 分别统计肿瘤和正常样本的数量
    group_by(cancer_type, primary_disease) %>%
    summarise(
      tumor_count = ifelse("Primary Tumor" %in% sample_type, 
                          sum(count[sample_type == "Primary Tumor"]), 0),
      normal_count = ifelse("Matched Normal" %in% sample_type, 
                           sum(count[sample_type == "Matched Normal"]), 0),
      .groups = 'drop'
    ) %>%
    # 生成标签：如果同时有肿瘤和正常样本，显示两者数量
    mutate(
      label = ifelse(
        tumor_count > 0 & normal_count > 0,
        # 同时存在肿瘤和正常样本：显示T和N的数量
        paste0(primary_disease, " (", cancer_type, ", T:", tumor_count, ", N:", normal_count, ")"),
        # 只有肿瘤样本：只显示肿瘤数量
        ifelse(
          tumor_count > 0,
          paste0(primary_disease, " (", cancer_type, ", T:", tumor_count, ")"),
          # 只有正常样本：只显示正常样本数量（这种情况较少）
          paste0(primary_disease, " (", cancer_type, ", N:", normal_count, ")")
        )
      )
    )
  
  # 图例标签
  legend_labels <- c(
    "Matched Normal" = "Solid Tissue Normal",
    "Primary Tumor" = "Primary Tumor"
  )  
  
  # 7. 绘图（无填充版本）
  p <- ggplot(processed_data, aes(x = cancer_type, y = !!sym(value_cols))) +
    # 箱线图 - 无填充版本
    geom_boxplot(
      aes(color = sample_type),
      fill = "transparent",
      outlier.shape = NA, 
      alpha = 0.7,
      size = 0.8
    ) +
    # 散点图（jitter points）
    geom_jitter(
      aes(color = sample_type), 
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
      size = 1.5, 
      alpha = 0.6
    ) +
    # 自定义颜色
    scale_color_manual(
      values = custom_colors,
      labels = legend_labels
    ) +
    # 修改：使用新的x轴标签
    scale_x_discrete(
      labels = setNames(label_data$label, label_data$cancer_type)
    ) +
    # 添加显著性标记
    {if (nrow(stat_test_results) > 0)
      geom_text(
        data = stat_test_results,
        aes(x = cancer_type, 
            y = max(processed_data[[value_cols]], na.rm = TRUE) * 1.1, 
            label = significance),
        inherit.aes = FALSE, 
        size = 5, 
        fontface = "bold"
      )
    } +
    # 图表美化
    labs(
      title = title,
      x = xlabel,
      y = ylabel,
      color = "Sample Type"
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0, vjust = 0.2, face = "bold", size = 14),
      plot.title.position = "panel",
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, 
                                face = "bold", size = 12, color = "black"),
      legend.text = element_text(size = 16),
      legend.title = element_blank(),
      legend.position = "top",
      panel.border = element_rect(colour = "black", fill = NA, size = 1.5, 
                                 linetype = "solid")
    )
  
  return(p)
}

custom_colors <- c('Matched Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- generate_gene_expression_plot(sample_level_expre_df,"expression_value",custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.sampletype.expression.tpm.wilcoxon.test.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=12,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 12, dpi = 300)

custom_colors <- c('Matched Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," expression (n=",dim(sample_level_expre_df)[1],")",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- generate_gene_expression_plot(sample_level_expre_df,"expression_value_log2",custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.sampletype.expression.tpm.log.wilcoxon.test.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=12,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 12, dpi = 300)

## style3: 按照Primary Tumor和Solid Tumor Normal计算中位表达后分布图 （图例按照Primary表达量从低到高）
barplot_func <- function(data_df,value_cols,category_cols,type_col,order_level,custom_colors,title,xlabel,ylabel) {
    # 确保分组列是因子，并按 order_level 排序（即使某些分组不存在于数据中）
    data_df[[category_cols]] <- factor(
        data_df[[category_cols]],
        levels = order_level
    )
    # 提取数据并去除缺失值（保留因子顺序）
    subdf <- data_df[, c(value_cols, category_cols,type_col)] %>% 
        na.omit() %>%
        mutate(across(all_of(category_cols), ~ factor(., levels = order_level)))  # 保持因子顺序
    
    p_Barplot <- ggplot(subdf, 
              aes(x = !!sym(category_cols), 
                  y = !!sym(value_cols),
                 fill = !!sym(type_col)))
    # 箱线图
    p_Barplot <- p_Barplot + geom_bar(stat = "identity",  position = position_dodge(width = 0.7),width = 0.6,color="black") + scale_fill_manual(values = custom_colors)
    p_Barplot <- p_Barplot + theme_classic() + theme(
            plot.title = element_text(hjust = 0,vjust= 0.2, face = "bold", size = 14),
            plot.title.position = "panel",
            axis.title = element_text(size = 16),  # 轴标签字体大小
            axis.text = element_text(size = 12),   # 刻度标签字体大小
            axis.text.x = element_text(angle=90,vjust=0.5,hjust = 1,face="bold",size=12,color="black"),
            legend.text = element_text(size = 16),  # 图例字体大小
            legend.title = element_blank(),
            legend.key.size = unit(1.5, "lines"),
            legend.spacing.x = unit(0.5, "cm"),
            legend.position = "top",
            panel.border = element_rect(colour = "black", 
                                        fill = NA, 
                                        size = 1.5, 
                                        linetype = "solid"))
    p_Barplot <- p_Barplot + labs(x = xlabel,y = ylabel,title=title) 
    return (p_Barplot)
}

tumor_medians <- sample_type_level_expre_df %>%
  filter(sample_type == "Primary Tumor") %>%
  group_by(cancer_type) %>%
  summarise(median_value = median(expression_value, na.rm = TRUE)) %>%
  arrange(median_value)

# 获取排序后的 cancer_type 顺序
cancer_order <- tumor_medians$cancer_type

# 创建 disease_label，并根据样本类型调整格式
sample_type_level_expre_df <- sample_type_level_expre_df %>%
  mutate(
    # 生成 disease_label
    primary_disease_label = case_when(
      sample_type == "Solid Tissue Normal" ~ 
        paste0(primary_disease, " (", cancer_type, " matched normal)"),
      TRUE ~ paste0(primary_disease, " (", cancer_type, ")")
    ),
    
    # 将 cancer_type 转换为因子，按照中位值排序
    cancer_type = factor(cancer_type, levels = cancer_order)
  )

# 创建正确的 disease_label 顺序
# 1. 获取每个癌症类型的样本类型情况
cancer_type_info <- sample_type_level_expre_df %>%
  group_by(cancer_type) %>%
  summarise(
    has_normal = any(sample_type == "Solid Tissue Normal"),
    .groups = 'drop'
  ) %>%
  arrange(cancer_type)  # 按照中位值排序

# 2. 构建正确的 disease_label 顺序
disease_label_order <- c()
for (cancer in cancer_order) {
  cancer_info <- cancer_type_info %>% filter(cancer_type == cancer)
  
  if (cancer_info$has_normal) {
    # 如果存在正常样本，先添加正常样本的标签
    normal_label <- paste0(
      unique(sample_type_level_expre_df$primary_disease[sample_type_level_expre_df$cancer_type == cancer]),
      " (", cancer, " matched normal)"
    )
    disease_label_order <- c(disease_label_order, normal_label)
  }
  
  # 然后添加肿瘤样本的标签
  tumor_label <- paste0(
    unique(sample_type_level_expre_df$primary_disease[sample_type_level_expre_df$cancer_type == cancer]),
    " (", cancer, ")"
  )
  disease_label_order <- c(disease_label_order, tumor_label)
}

# 3. 将 disease_label 转换为因子，按照正确的顺序
sample_type_level_expre_df <- sample_type_level_expre_df %>%
  mutate(
    primary_disease_label = factor(primary_disease_label, levels = disease_label_order)
  )

# 4. 按照正确的顺序排序数据
sample_type_level_expre_df <- sample_type_level_expre_df %>%
  arrange(primary_disease_label)

# 获取并返回 disease_label 的顺序
disease_label_order <- levels(sample_type_level_expre_df$primary_disease_label)

custom_colors <- c('Solid Tissue Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," median expression of different cancer types",sep='')
xlabel = ""
ylabel="Transcripts Per Million (TPM)"
expression_plot <- barplot_func(sample_type_level_expre_df,"expression_value",'primary_disease_label',"sample_type",disease_label_order,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.expression.median.tpm.bycancertype.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)

custom_colors <- c('Solid Tissue Normal'='#00aa36','Primary Tumor'='#df0028')
title=paste(gene," median expression of different cancer types",sep='')
xlabel = ""
ylabel="Expression (log2(TPM+1))"
expression_plot <- barplot_func(sample_type_level_expre_df,"expression_value_log2",'primary_disease_label',"sample_type",disease_label_order,custom_colors,title,xlabel,ylabel)
expression_plot_file <- paste(gene,"TCGA.expression.median.tpm.log.bycancertype.boxplot.pdf",sep=".")
expression_plot_file <- file.path(odir,expression_plot_file)
pdf(expression_plot_file,height=10,width=20)
print(expression_plot)
dev.off()
expression_plot_file <- gsub("pdf","png",expression_plot_file)
ggsave(expression_plot_file, plot = expression_plot, width = 20, height = 10, dpi = 300)


















