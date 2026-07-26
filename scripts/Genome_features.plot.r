# ============================================================
# 脚本名: bubble_plot_final.R
# 功能: 生成 Figure 3 风格气泡图（数字在气泡右上方，不同指标不同颜色）
# ============================================================

# 加载包
library(tidyverse)
library(RColorBrewer)

# ========== 1. 读取数据 ==========
file_path <- "species_data.txt"  # 请修改为你的文件路径

df <- read.delim(file_path, 
                 sep = "\t", 
                 header = TRUE, 
                 fileEncoding = "UTF-8-BOM")

print("读取到的列名：")
print(colnames(df))

# ========== 2. 数据清洗和重命名 ==========
plot_data <- df %>%
  select(Species,
         Contigs = Number.of.contigs,
         Size_Mb = Size..Mb.,
         GC = GC.content....,
         N50_Mb = N50..Mb.,
         Gene_len_kb = Average.gene.length..kb.,
         Exons_per_gene = Average.exons.per.gene,
         Exon_len_kb = Average.exon.length..kb.,
         Introns_per_gene = Average.number.of.introns.per.gene,
         Intron_len_bp = Average.intron.length..bp.)

# 检查数据
print("数据前3行：")
print(head(plot_data, 3))

# ========== 3. 转换成长格式 ==========
long_data <- plot_data %>%
  pivot_longer(cols = -Species, 
               names_to = "Metric", 
               values_to = "Value")

# 处理缺失值
if(any(is.na(long_data$Value))) {
  warning("数据中存在缺失值，将用该指标的中位数填充")
  long_data <- long_data %>%
    group_by(Metric) %>%
    mutate(Value = ifelse(is.na(Value), median(Value, na.rm = TRUE), Value)) %>%
    ungroup()
}

# ========== 4. 添加显示文本（格式化数字）==========
long_data <- long_data %>%
  mutate(DisplayText = case_when(
    # 整数类型
    Metric == "Contigs" ~ as.character(round(Value, 0)),
    Metric == "Exons_per_gene" ~ sprintf("%.2f", Value),
    Metric == "Introns_per_gene" ~ sprintf("%.2f", Value),
    # 百分比（保留2位小数）
    Metric == "GC" ~ sprintf("%.2f", Value),
    # 大小类（保留2位小数）
    Metric == "Size_Mb" ~ sprintf("%.2f", Value),
    Metric == "N50_Mb" ~ sprintf("%.2f", Value),
    # 长度类（kb 单位，保留2位小数）
    Metric == "Gene_len_kb" ~ sprintf("%.2f", Value),
    Metric == "Exon_len_kb" ~ sprintf("%.2f", Value),
    # 长度类（bp 单位，保留2位小数）
    Metric == "Intron_len_bp" ~ sprintf("%.2f", Value),
    # 默认
    TRUE ~ sprintf("%.1f", Value)
  ))

# ========== 5. 设置菌株顺序 ==========
# 按照数据原始顺序（第一列的顺序）
species_order <- plot_data$Species
long_data$Species <- factor(long_data$Species, levels = rev(species_order))

# ========== 6. 设置指标顺序和标签 ==========
metric_order <- c("Contigs", "Size_Mb", "GC", "N50_Mb", 
                  "Gene_len_kb", "Exons_per_gene", "Exon_len_kb", 
                  "Introns_per_gene", "Intron_len_bp")

metric_labels <- c("Number of contigs", "Size (Mb)", "GC content (%)", "N50 (Mb)",
                   "Average gene length (kb)", "Average exons per gene", 
                   "Average exon length (kb)", "Average number of introns per gene", 
                   "Average intron length (bp)")

long_data$Metric <- factor(long_data$Metric, 
                           levels = metric_order,
                           labels = metric_labels)

# ========== 7. 为每个指标分配颜色 ==========
colors_manual <- c(
  "Number of contigs" = "#E41A1C",
  "Size (Mb)" = "#377EB8",
  "GC content (%)" = "#4DAF4A",
  "N50 (Mb)" = "#984EA3",
  "Average gene length (kb)" = "#FF7F00",
  "Average exons per gene" = "#FFFF33",
  "Average exon length (kb)" = "#A65628",
  "Average number of introns per gene" = "#F781BF",
  "Average intron length (bp)" = "#999999"
)

# ========== 8. 创建文本数据（用于显示数字）==========
text_data <- long_data %>%
  mutate(
    x_numeric = as.numeric(Metric),
    text_x = x_numeric + 0.12,
    text_y = as.numeric(Species) + 0.12
  )

text_data$Species <- factor(text_data$Species, levels = rev(species_order))

# ========== 9. 绘制气泡图 ==========
p <- ggplot() +
  
  # 气泡层
  geom_point(data = long_data, 
             aes(x = Metric, y = Species, 
                 size = Value, 
                 fill = Metric),
             shape = 21, 
             color = "black", 
             stroke = 0.2) +
  
  # 数字层
  geom_text(data = text_data, 
            aes(x = text_x, y = text_y, label = DisplayText),
            color = "black", 
            size = 2.5,
            hjust = 0,
            vjust = 0,
            show.legend = FALSE) +
  
  # 气泡大小范围
  scale_size_continuous(range = c(2, 12), 
                        name = "Value",
                        breaks = pretty(range(long_data$Value, na.rm = TRUE), 4)) +
  
  # 颜色填充
  scale_fill_manual(name = "Genomic features", 
                    values = colors_manual,
                    labels = metric_labels) +
  
  # 主题美化
  theme_minimal(base_size = 10) +
  
  theme(
    axis.text.x = element_text(angle = 45, 
                               hjust = 1, 
                               size = 9,
                               face = "bold"),
    axis.text.y = element_text(size = 7, hjust = 1),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  
  labs(title = "Genomic features of Talaromyces and related strains")

# ========== 10. 显示和保存 ==========
print(p)

n_species <- length(unique(long_data$Species))
plot_height <- max(10, n_species * 0.35)

ggsave("bubble_plot_colored.png", 
       plot = p, 
       width = 16,
       height = plot_height, 
       dpi = 300, 
       limitsize = FALSE)

ggsave("bubble_plot_colored.pdf", 
       plot = p, 
       width = 16, 
       height = plot_height, 
       limitsize = FALSE)

cat("\n========== 绘图完成 ==========\n")
cat(sprintf("物种数量: %d\n", n_species))
cat(sprintf("图片尺寸: 16 x %.1f 英寸\n", plot_height))
cat("输出文件:\n")
cat("  - bubble_plot_colored.png\n")
cat("  - bubble_plot_colored.pdf\n")
cat("===============================\n")
