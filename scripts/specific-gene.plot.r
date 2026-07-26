# ============================================================
# 气泡图：功能蛋白分组分布图（X轴和Y轴都随机排序）
# ============================================================

library(ggplot2)
library(dplyr)

# 读取数据
data <- read.table("data.txt", header = TRUE, sep = "\t", quote = "")

# --- 数据处理：统计每个功能蛋白在各分组中的分布情况 ---
bubble_data <- data %>%
  group_by(Groups, Description) %>%
  summarise(Frequency = n(), .groups = 'drop') %>%
  rename(Functional_Group = Groups,
         Functional_Protein = Description,
         Occurrence = Frequency)

# --- X轴（功能蛋白）随机打乱顺序 ---
set.seed(123)

all_proteins <- unique(bubble_data$Functional_Protein)
random_order_x <- sample(all_proteins)
bubble_data$Functional_Protein <- factor(bubble_data$Functional_Protein, 
                                          levels = random_order_x)

# --- Y轴（功能分组）随机打乱顺序 ---
all_groups <- unique(bubble_data$Functional_Group)
random_order_y <- sample(all_groups)
bubble_data$Functional_Group <- factor(bubble_data$Functional_Group, 
                                        levels = random_order_y)

# --- 颜色定义 ---
n_groups <- length(unique(bubble_data$Functional_Group))
group_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
)[1:n_groups]

names(group_colors) <- levels(bubble_data$Functional_Group)

# --- 绘制气泡图（X轴和Y轴都随机）---
p <- ggplot(bubble_data, aes(x = Functional_Protein, 
                             y = Functional_Group, 
                             size = Occurrence)) +
  geom_point(aes(color = Functional_Group), alpha = 0.85) +
  
  # 调整气泡大小范围
  scale_size_continuous(
    range = c(2, 3.5),
    name = "出现次数",
    breaks = c(1, 2, 3, 4),
    labels = c("1", "2", "3", "4+")
  ) +
  
  # 颜色
  scale_color_manual(values = group_colors, name = "功能分组") +
  
  # 紧凑布局
  scale_x_discrete(expand = c(0.05, 0.05)) +
  scale_y_discrete(expand = c(0.05, 0.05)) +
  
  # 主题
  theme_bw() +
  theme(
    # X轴文字45度倾斜
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, 
                               size = 8, face = "bold"),
    axis.text.y = element_text(size = 9, hjust = 1),
    axis.title = element_text(size = 10, face = "bold"),
    panel.border = element_rect(linewidth = 0.5),
    panel.grid.major = element_line(color = "gray95", linewidth = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "right",
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.5, "cm"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "gray30")
  ) +
  
  labs(
    title = "功能蛋白在各分组中的分布",
    subtitle = "气泡大小表示出现次数（X轴和Y轴均随机排序）",
    x = "功能蛋白",
    y = "功能分组"
  ) +
  
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5)))

# 显示图形
print(p)

# 保存图形
ggsave("function_protein_bubble_random_both.png", 
       plot = p, width = 14, height = 8, dpi = 300, bg = "white")
ggsave("function_protein_bubble_random_both.pdf", 
       plot = p, width = 14, height = 8, bg = "white")

# 输出随机顺序
cat("\n当前 X 轴功能蛋白的随机顺序（从左到右）：\n")
cat(paste(levels(bubble_data$Functional_Protein), collapse = " -> "))
cat("\n\n当前 Y 轴功能分组的随机顺序（从上到下）：\n")
cat(paste(levels(bubble_data$Functional_Group), collapse = " -> "))
cat("\n")
