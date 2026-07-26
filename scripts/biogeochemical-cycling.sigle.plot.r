#!/usr/bin/env Rscript

# 加载包
library(tidyverse)
library(pheatmap)
library(grid)
library(gridExtra)

# ============================================================
# 配置：四个循环
# ============================================================

cycles <- list(
  Nitrogen = list(
    file = "N.txt",
    col2 = "passway",
    heatmap_colors = c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"),
    pathway_colors = c(
      "Assimilatory nitrate reduction" = "#85C1E9",
      "Dissimilatory nitrate reduction" = "#7DCEA0",
      "Denitrification" = "#F1948A",
      "Organic degradation and synthesis" = "#F0B27A",
      "Others" = "#BDC3C7"
    ),
    cycle_color = "#3182BD"
  ),
  Sulfur = list(
    file = "S.txt",
    col2 = "passway",
    heatmap_colors = c("#F7FCF5", "#A6DBA0", "#4DAF4A", "#006837"),
    pathway_colors = c(
      "Organic sulphur transformation" = "#E64B35",
      "Assimilatory sulphate reduction" = "#4DBBD5",
      "Dissimilatory sulphur reduction and oxidation" = "#00A087",
      "Linkages between inorganic and organic sulphur transformation" = "#7E6148",
      "Sulphur oxidation" = "#3C5488"
    ),
    cycle_color = "#4DAF4A"
  ),
  Methane = list(
    file = "CH4.txt",
    col2 = "passway",
    heatmap_colors = c("#FFF7EC", "#FEE8C8", "#FDD49E", "#FDBB84", "#FC8D59", "#D7301F"),
    pathway_colors = c(
      "Aceticlastic methanogenesis" = "#D95F02",
      "Anaerobic oxidation of methane (AOM)" = "#7570B3",
      "Hydrogenotrophic methanogenesis" = "#E7298A",
      "Oxidation of formaldehyde" = "#66A61E",
      "Oxidation of formate" = "#E6AB02",
      "RuMP cycle" = "#A6761D",
      "Serine cycle" = "#666666"
    ),
    cycle_color = "#FC8D59"
  ),
  Carbon = list(
    file = "Ccyc.txt",
    col2 = "passway",
    heatmap_colors = c("#F3F0FF", "#D9D2E9", "#B4A7D6", "#8E7CC3", "#674EA7", "#351C75"),
    pathway_colors = c(
      "Carbon Release" = "#E64B35",
      "Organic Biosynthesis" = "#4DBBD5",
      "Organic Degration" = "#00A087",
      "Organic Transformation" = "#7E6148",
      "Transport" = "#3C5488"
    ),
    cycle_color = "#8E7CC3"
  )
)

# ============================================================
# 处理所有循环，合并到一个矩阵（按基因展示）
# ============================================================

all_mats <- list()
all_gene_names <- c()
all_pathway_annotations <- c()  # 每个基因对应的通路
all_cycle_labels <- c()
all_pathway_colors_combined <- c()
cycle_color_map <- c()

# 存储第一个文件的物种顺序
strain_order <- NULL

# 用于重名基因处理
gene_counter <- 1

for (cycle_name in names(cycles)) {
  cat("\n========== 处理", cycle_name, "==========\n")
  
  config <- cycles[[cycle_name]]
  
  df <- read.delim(config$file, header = TRUE, sep = "\t", 
                   stringsAsFactors = FALSE, check.names = FALSE)
  
  colnames(df)[2] <- "passway"
  
  # 名称统一
  df$passway[df$passway == "Anaerobic oxidation of methane (AOM)；Central methanogenic pathway"] <- "Anaerobic oxidation of methane (AOM)"
  df$passway[df$passway == "Hydrogenotrophic methanogenesis；Anaerobic oxidation of methane (AOM)"] <- "Hydrogenotrophic methanogenesis"
  
  df$passway <- gsub("\n", " ", df$passway)
  df$passway <- gsub("  +", " ", df$passway)
  df$passway <- trimws(df$passway)
  
  df <- df[!is.na(df$Gene), ]
  df <- df[df$Gene != "", ]
  
  strain_cols <- colnames(df)[3:ncol(df)]
  
  # 用第一个文件的物种顺序
  if (is.null(strain_order)) {
    strain_order <- strain_cols
  }
  
  # ========== 改动：不汇总，每个基因单独作为一行 ==========
  # 直接用原始数据，基因列作为行名
  mat <- df %>%
    select(all_of(strain_cols)) %>%
    as.matrix()
  mode(mat) <- "numeric"
  
  # 处理基因名，确保唯一性
  gene_names <- df$Gene
  # 如果有重复基因名，添加后缀
  duplicated_genes <- duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE)
  if (any(duplicated_genes)) {
    warning(paste("在", cycle_name, "中发现重复基因名，将添加后缀以区分"))
    gene_names[duplicated_genes] <- paste0(gene_names[duplicated_genes], "_", 
                                           seq_len(sum(duplicated_genes)))
  }
  
  rownames(mat) <- gene_names
  colnames(mat) <- strain_cols
  mat[is.na(mat)] <- 0
  
  all_mats[[cycle_name]] <- mat
  
  # 记录基因名、对应的通路和循环
  all_gene_names <- c(all_gene_names, gene_names)
  all_pathway_annotations <- c(all_pathway_annotations, df$passway)
  all_cycle_labels <- c(all_cycle_labels, rep(cycle_name, length(gene_names)))
  
  # 通路颜色
  pathway_colors <- config$pathway_colors
  for (p in df$passway) {
    if (!(p %in% names(pathway_colors))) {
      pathway_colors[p] <- "#BEBEBE"
    }
  }
  all_pathway_colors_combined <- c(all_pathway_colors_combined, pathway_colors[df$passway])
  
  cycle_color_map[cycle_name] <- config$cycle_color
  
  cat("基因数:", nrow(mat), "\n")
}

cat("\n物种顺序（按第一个文件）:", paste(strain_order, collapse = ", "), "\n")
cat("总基因数:", length(all_gene_names), "\n")

# ============================================================
# 合并所有矩阵（按基因展示）
# ============================================================

# 构建合并矩阵：行=基因，列=菌株
combined_mat <- matrix(0, nrow = length(all_gene_names), ncol = length(strain_order))
rownames(combined_mat) <- all_gene_names
colnames(combined_mat) <- strain_order

# 填充矩阵
row_idx <- 1
for (name in names(all_mats)) {
  mat <- all_mats[[name]]
  for (i in seq_len(nrow(mat))) {
    for (j in seq_len(ncol(mat))) {
      strain <- colnames(mat)[j]
      if (strain %in% strain_order) {
        combined_mat[row_idx, strain] <- mat[i, j]
      }
    }
    row_idx <- row_idx + 1
  }
}

cat("\n合并矩阵:", nrow(combined_mat), "基因 ×", ncol(combined_mat), "菌株\n")

# ============================================================
# 行注释（基因所属通路 + 所属循环）
# ============================================================

# 创建双重行注释
annotation_row <- data.frame(
  Pathway = all_pathway_annotations,
  Cycle = all_cycle_labels,
  row.names = all_gene_names,
  stringsAsFactors = FALSE
)

# ============================================================
# 行分隔线位置（按循环分隔）
# ============================================================

gaps <- c()
current_row <- 0
for (name in names(cycles)) {
  n_genes <- sum(all_cycle_labels == name)
  current_row <- current_row + n_genes
  if (current_row < nrow(combined_mat)) {
    gaps <- c(gaps, current_row)
  }
}

# ============================================================
# 图片尺寸（因为基因数可能很多，调整尺寸）
# ============================================================

n_rows <- nrow(combined_mat)
n_cols <- ncol(combined_mat)

# 根据基因数调整行高和图片尺寸
height <- max(20, n_rows * 0.3)   # 每个基因0.3英寸
width <- max(22, n_cols * 0.5)

cell_width <- 13
cell_height <- 10
row_font <- ifelse(n_rows > 100, 5, 7)   # 基因多时减小字体
col_font <- 7
number_fontsize <- ifelse(n_rows > 100, 5, 6.5)

cat("\n图片尺寸:", round(width, 1), "x", round(height, 1), "inches\n")

# ============================================================
# 图例
# ============================================================

max_val <- max(combined_mat)

if (max_val <= 10) {
  legend_step <- 2
} else if (max_val <= 20) {
  legend_step <- 5
} else {
  legend_step <- 10
}

# ============================================================
# 绘制热图 - 行=基因，列=菌株
# ============================================================

# 为通路创建颜色映射
unique_pathways <- unique(all_pathway_annotations)
pathway_color_map <- all_pathway_colors_combined[match(unique_pathways, all_pathway_annotations)]
names(pathway_color_map) <- unique_pathways

pdf("heatmap_all_cycles.pdf", width = width, height = height)

pheatmap(combined_mat,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_colors = list(
           Pathway = pathway_color_map,
           Cycle = cycle_color_map
         ),
         color = colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))(100),
         main = "Gene Copy Number in Talaromyces (by Gene)",
         fontsize = 10,
         fontsize_row = row_font,
         fontsize_col = col_font,
         angle_col = 90,
         border_color = "white",
         cellwidth = cell_width,
         cellheight = cell_height,
         gaps_row = gaps,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_legend = TRUE,
         annotation_names_row = TRUE,
         legend = TRUE,
         breaks = seq(0, max_val, length.out = 100),
         legend_breaks = seq(0, max_val, by = legend_step),
         legend_labels = as.character(seq(0, max_val, by = legend_step)),
         display_numbers = TRUE,
         number_color = "#333333",
         number_format = "%.0f",
         fontsize_number = number_fontsize
)

dev.off()

cat("\n✅ PDF已保存为 heatmap_all_cycles.pdf\n")

# ============================================================
# PNG
# ============================================================

png("heatmap_all_cycles.png", width = width * 150, height = height * 150, res = 150)

pheatmap(combined_mat,
         scale = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_row = annotation_row,
         annotation_colors = list(
           Pathway = pathway_color_map,
           Cycle = cycle_color_map
         ),
         color = colorRampPalette(c("#F7FBFF", "#DEEBF7", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C"))(100),
         main = "Gene Copy Number in Talaromyces (by Gene)",
         fontsize = 10,
         fontsize_row = row_font,
         fontsize_col = col_font,
         angle_col = 45,
         border_color = "white",
         cellwidth = cell_width,
         cellheight = cell_height,
         gaps_row = gaps,
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_legend = TRUE,
         annotation_names_row = TRUE,
         legend = TRUE,
         legend_breaks = seq(0, max_val, by = legend_step),
         legend_labels = as.character(seq(0, max_val, by = legend_step)),
         display_numbers = TRUE,
         number_color = "#333333",
         number_format = "%.0f",
         fontsize_number = number_fontsize
)

dev.off()

cat("✅ PNG已保存为 heatmap_all_cycles.png\n")

# ============================================================
# 通路图例
# ============================================================

pdf("pathway_legend_all_cycles.pdf", width = 14, height = 10)

legend_list <- list()

for (cycle_name in names(cycles)) {
  config <- cycles[[cycle_name]]
  pw_names <- names(config$pathway_colors)
  pw_colors <- config$pathway_colors
  
  legend_df <- data.frame(
    Pathway = factor(pw_names, levels = pw_names),
    Color = pw_colors,
    stringsAsFactors = FALSE
  )
  
  p <- ggplot(legend_df, aes(x = Pathway, y = 1, fill = Pathway)) +
    geom_tile(color = "black", linewidth = 0.5, height = 0.6) +
    scale_fill_manual(values = pw_colors) +
    labs(title = cycle_name) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
    coord_flip()
  
  legend_list[[cycle_name]] <- p
}

grid.arrange(grobs = legend_list, ncol = 1)
dev.off()

cat("✅ 通路图例已保存为 pathway_legend_all_cycles.pdf\n")

cat("\n========== 处理完成 ==========\n")
