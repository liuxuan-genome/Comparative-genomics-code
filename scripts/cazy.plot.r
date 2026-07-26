library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
# 读取数据
data <- read.table("cazy.special.gene.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# 保持物种原始顺序
original_species_order <- colnames(data)[-1]

# 处理通路信息
pathway_info <- data[, "passway", drop = FALSE]
gene_data <- data[, original_species_order, drop = FALSE]

# 转换数据格式为长格式
gene_data$Gene <- rownames(gene_data)
pathway_info$Gene <- rownames(pathway_info)

gene_data_long <- melt(gene_data, id.vars = "Gene", 
                       variable.name = "Species", 
                       value.name = "Copy_Number")

gene_data_long <- merge(gene_data_long, pathway_info, by = "Gene")

# 确保Species因子水平严格按照原始顺序
gene_data_long$Species <- factor(gene_data_long$Species, levels = original_species_order)

# 反转基因顺序使第一个基因在顶部
gene_data_long$Gene <- factor(gene_data_long$Gene, levels = rev(rownames(gene_data)))

# **只保留拷贝数大于0的数据**
gene_data_filtered <- gene_data_long %>%
  filter(Copy_Number > 0)

# 创建形状类型列，用于图例
gene_data_filtered$Shape_Type <- ifelse(gene_data_filtered$Species == "Talaromyces_barcinensis", 
                                        "★ T. barcinensis", "● Other species")

# 分离目标菌株和其他菌株的数据
target_data <- gene_data_filtered[gene_data_filtered$Shape_Type == "★ T. barcinensis", ]
other_data <- gene_data_filtered[gene_data_filtered$Shape_Type == "● Other species", ]

# 创建CAZy类别颜色方案（更鲜明的颜色）
cazy_colors <- c(
  "CBM" = "#E41A1C",  # 红色
  "GH" = "#377EB8",   # 蓝色
  "AA" = "#4DAF4A"    # 绿色
)

# 最终版本：大小=拷贝数，颜色=CAZy类别
p_final <- ggplot() +
  # 其他物种用圆圈
  geom_point(data = other_data,
             aes(x = Species, y = Gene, 
                 size = Copy_Number, 
                 color = passway,
                 shape = Shape_Type),
             alpha = 0.8) +
  # 目标物种用菱形（实心，无边框）
  geom_point(data = target_data,
             aes(x = Species, y = Gene, 
                 size = Copy_Number, 
                 color = passway,
                 shape = Shape_Type),
             alpha = 1) +
  # 形状映射：圆圈和菱形
  scale_shape_manual(values = c("★ T. barcinensis" = 18, "● Other species" = 16),
                     name = "Species") +
  # 气泡大小映射拷贝数（最小值3确保拷贝数1的点清晰可见）
  scale_size_continuous(range = c(3, 12), 
                       name = "Copy Number",
                       breaks = c(1, 2, 4, 6, 8, 10, 15, 20, 25)) +
  # 颜色映射CAZy类别
  scale_color_manual(values = cazy_colors,
                    name = "CAZy Class") +
  # 分面展示
  facet_grid(passway ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 8),
    strip.text.y = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title = "CAZymes Gene Copy Number Distribution Across Talaromyces Species",
    subtitle = "◆ Talaromyces barcinensis  |  ● Other species  |  Color = CAZy class  |  Size = Copy number",
    x = "Talaromyces Species",
    y = "Gene Families"
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 6)),
    color = guide_legend(order = 2, override.aes = list(size = 5)),
    size = guide_legend(order = 3)
  )

# 显示图形
print(p_final)

# 保存图形
ggsave("CAZymes_faceted_final_improved.pdf", p_final, width = 18, height = 12, dpi = 300)
ggsave("CAZymes_faceted_final_improved.png", p_final, width = 18, height = 12, dpi = 300)

# 统计输出
cat("\n========================================\n")
cat("Talaromyces barcinensis Gene Copy Numbers\n")
cat("========================================\n")
barcinensis_data <- target_data[order(target_data$Gene), ]
print(barcinensis_data[, c("Gene", "passway", "Copy_Number")], row.names = FALSE)

cat("\nSummary for T. barcinensis:\n")
cat("  Genes detected:", nrow(barcinensis_data), "\n")
cat("  Total copies:", sum(barcinensis_data$Copy_Number), "\n")
cat("  Mean copies per gene:", round(mean(barcinensis_data$Copy_Number), 2), "\n")
cat("  Range:", min(barcinensis_data$Copy_Number), "-", max(barcinensis_data$Copy_Number), "\n")

# 按CAZy类别统计T. barcinensis
cat("\nBy CAZy class:\n")
barcinensis_by_class <- barcinensis_data %>%
  group_by(passway) %>%
  summarise(
    genes = n(),
    total_copies = sum(Copy_Number),
    mean_copies = round(mean(Copy_Number), 1),
    .groups = 'drop'
  )
print(barcinensis_by_class)

# 确认物种顺序
cat("\n========================================\n")
cat("Species order on x-axis (left to right):\n")
cat("========================================\n")
for(i in 1:length(original_species_order)) {
  marker <- ifelse(original_species_order[i] == "Talaromyces_barcinensis", 
                   " ← TARGET", "")
  cat(sprintf("%2d. %s%s\n", i, original_species_order[i], marker))
}
