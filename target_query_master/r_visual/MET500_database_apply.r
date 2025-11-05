## 用于进行MET500数据可视化
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

parse_json_with_nan <- function(file_path) {
	  json_content <- readLines(file_path, warn = FALSE)
  json_content <- gsub("NaN", "null", json_content)
    # 同时处理其他可能的非标准数值表示
    json_content <- gsub("Infinity", "null", json_content)
    json_content <- gsub("-Infinity", "null", json_content)
      return(jsonlite::fromJSON(json_content))
}

met500_gene_data <- parse_json_with_nan(input_file)

sample_level_expre_df <- met500_gene_data$sample_level_expression$samples

sample_level_expre_df <- sample_level_expre_df %>%
  # 首先按tissue分组，计算每种组织的总样本数
  group_by(tissue) %>%
  mutate(
    tissue_count = n(),  # 计算每种组织类型的样本数
    tissue_label = paste0(tissue, "(n=", tissue_count, ")")  # 创建组织标签
  ) %>%
  # 然后按tissue和biopsy_tissue分组，计算每种活检组织的样本数
  group_by(tissue, biopsy_tissue) %>%
  mutate(
    biopsy_tissue_count = n(),  # 计算在每种组织下各活检组织的样本数
    biopsy_tissue_label = paste0(biopsy_tissue, "(n=", biopsy_tissue_count, ")")  # 创建活检组织标签
  ) %>%
  ungroup() 

sample_level_expre_df$tissue_count <- NULL
sample_level_expre_df$biopsy_tissue_count <- NULL

sample_level_expre_df$expression_value_log2 <- log2(sample_level_expre_df$expression_value+1)

boxplot_notest_func <- function(score_cols, category_cols, gene_expre_df, order_level,custom_colors,title,xlabel,ylabel) {
    # 确保分组列是因子，并按 order_level 排序（即使某些分组不存在于数据中）
    gene_expre_df[[category_cols]] <- factor(
        gene_expre_df[[category_cols]],
        levels = order_level
    )
    
    # 提取数据并去除缺失值（保留因子顺序）
    subdf <- gene_expre_df[, c(score_cols, category_cols,"sample_id")] %>% 
        na.omit() %>%
        mutate(across(all_of(category_cols), ~ factor(., levels = order_level)))  # 保持因子顺序

    sample_counts <- subdf %>%
      distinct(across(all_of(c(category_cols, "sample_id")))) %>%
      count(across(all_of(category_cols)), name = "n") %>%
      mutate(label = paste0(!!sym(category_cols), "\n(n=", n, ")"))
    
    # 创建一个命名向量，用于scale_x_discrete的labels映射
    label_mapping <- setNames(sample_counts$label, sample_counts[[category_cols]])    

    p_Boxplot <- ggplot(subdf, 
              aes(x = !!sym(category_cols), 
                  y = !!sym(score_cols), 
                  fill = !!sym(category_cols))) # 使用原始列名进行填充
    # 箱线图
    p_Boxplot <- p_Boxplot + geom_boxplot(
      width = 0.5,
      alpha = 0.7,
      outlier.shape = NA,
      color = "black")
    
    # 抖动点
    p_Boxplot <- p_Boxplot + geom_jitter(
      width = 0.3,
      height = 0,
      size = 1.5,
      alpha = 0.4,
      shape = 21,
      color = "black",
      fill = "black")
    p_Boxplot <- p_Boxplot + scale_fill_manual(values = custom_colors)
    p_Boxplot <- p_Boxplot + scale_x_discrete(labels = label_mapping)
    p_Boxplot <- p_Boxplot + theme_classic() + theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.title = element_text(size = 16,face = "bold"),  # 轴标签字体大小
            axis.text = element_text(size = 12),   # 刻度标签字体大小
            axis.text.x = element_text(angle=15,vjust=0.5,face="bold",color="black"),
            legend.text = element_text(size = 16),  # 图例字体大小
            legend.position = "none") 
    p_Boxplot <- p_Boxplot + labs(x = xlabel,y = ylabel,title=title)
    return (p_Boxplot)
}

order_level_func <- function(data_df,category_col) {
    order_level <- data_df %>%
        group_by(.data[[category_col]]) %>%
        summarise(median_value = median(expression_value, na.rm = TRUE)) %>%
        arrange(desc(median_value)) %>%
        pull(.data[[category_col]])
    return (order_level)
}

# 定义函数：根据条形数量设置图形尺寸
set_plot_dimensions <- function(n_bars) {
  # 输入验证
  if (!is.numeric(n_bars) || n_bars <= 0) {
    stop("条形数量必须是大于0的数值")
  }
  
  # 基于条形数量设置不同的尺寸比例
  if (n_bars <= 5) {
    # 条形数量较少时使用较小尺寸
    plot_width <- 8
    plot_height <- 6
    ratio_type <- "小尺寸比例"
  } else if (n_bars <= 15) {
    # 中等数量条形使用适中尺寸
    plot_width <- 12
    plot_height <- 8
    ratio_type <- "中等尺寸比例"
  } else {
    # 条形数量多时使用较大尺寸
    plot_width <- 16
    plot_height <- 10
    ratio_type <- "大尺寸比例"
  }
  
  # 计算长宽比
  aspect_ratio <- plot_width / plot_height
  
  # 设置图形选项
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)
  
  # 返回信息
  return(list(
    width = plot_width,
    height = plot_height,
    aspect_ratio = round(aspect_ratio, 2),
    ratio_type = ratio_type,
    n_bars = n_bars
  ))
}

gene <- met500_gene_data$gene_information$gene_name

tissue_vector = unique(sample_level_expre_df$tissue)

tissue_vector <- na.omit(tissue_vector)

for (tissue in tissue_vector ) {
    tissue_sample_level_expre_df <- sample_level_expre_df[sample_level_expre_df$tissue==tissue,]
    order_level <- order_level_func(tissue_sample_level_expre_df,"biopsy_tissue")
    order_level <- order_level[!order_level %in% c(NA)]
    num_color <- length(order_level)
    value_color <- colorRampPalette(c("#08306B","#DEEBF7"))(num_color)
    names(value_color) <- order_level
    cancer <- paste(tissue,"cancer",sep=" ")
    ylabel <- paste(gene," Expression (FPKM)")
    plot_setting <- set_plot_dimensions(num_color)
    width <- plot_setting$width
    height <- plot_setting$height
    p <- boxplot_notest_func('expression_value', 'biopsy_tissue',tissue_sample_level_expre_df, order_level,value_color,cancer,'',ylabel)
    expression_plot_file <- paste(gene,tissue,"MET500.cancer.metastasis.expression.fpkm.boxplot.pdf",sep=".")
    expression_plot_file <- file.path(odir,expression_plot_file)
    pdf(expression_plot_file,height=height,width=width)
    print(p)
    dev.off()
    expression_plot_file <- gsub("pdf","png",expression_plot_file)
    ggsave(expression_plot_file, plot = p, width = width, height = height, dpi = 300)

    ## expression value calculatated by log2(x+1) 
    ylabel <- paste(gene," Expression-log2(FPKM+1)")
    p <- boxplot_notest_func('expression_value_log2', 'biopsy_tissue',tissue_sample_level_expre_df, order_level,value_color,cancer,'',ylabel)
    expression_plot_file <- paste(gene,tissue,"MET500.cancer.metastasis.expression.fpkm.log.boxplot.pdf",sep=".")
    expression_plot_file <- file.path(odir,expression_plot_file)
    pdf(expression_plot_file,height=height,width=width)
    print(p)
    dev.off()
    expression_plot_file <- gsub("pdf","png",expression_plot_file)
    ggsave(expression_plot_file, plot = p, width = width, height = height, dpi = 300)        
}
