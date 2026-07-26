# ============================================================
# 气泡图：功能蛋白分组分布图（调整气泡大小版）
# 打乱 Y 轴功能蛋白的顺序，避免对角线分布
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

# --- 关键修改：将 Y 轴（功能蛋白）的顺序随机打乱 ---
set.seed(123)  # 固定随机种子，方便复现

# 获取所有功能蛋白并打乱顺序
all_proteins <- unique(bubble_data$Functional_Protein)
random_order <- sample(all_proteins)
bubble_data$Functional_Protein <- factor(bubble_data$Functional_Protein, 
                                          levels = random_order)

# X 轴（功能分组）按总数排序
group_order <- bubble_data %>%
  group_by(Functional_Group) %>%
  summarise(Total = sum(Occurrence), .groups = 'drop') %>%
  arrange(desc(Total)) %>%
  pull(Functional_Group)

bubble_data$Functional_Group <- factor(bubble_data$Functional_Group, 
                                        levels = group_order)

# --- 颜色定义 ---
n_groups <- length(unique(bubble_data$Functional_Group))
group_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494"
)[1:n_groups]

names(group_colors) <- levels(bubble_data$Functional_Group)

# --- 绘制气泡图（调整气泡大小范围）---
p <- ggplot(bubble_data, aes(x = Functional_Group, 
                             y = Functional_Protein, 
                             size = Occurrence)) +
  geom_point(aes(color = Functional_Group), alpha = 0.85) +
  
  # 调整气泡大小：缩小范围，让1和2的差距不那么明显
  scale_size_continuous(
    range = c(2, 4),  # 从原来的 c(2,8) 改为 c(2,5)，差距变小
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, 
                               size = 9, face = "bold"
			       ),
    axis.text.y = element_text(size = 7, hjust = 1),
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
    subtitle = "气泡大小表示该功能蛋白在分组中出现的次数（Y轴随机排序）",
    x = "功能分组",
    y = "功能蛋白"
  ) +
  
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 2.5)))

# 显示图形
print(p)

# 保存图形
ggsave("function_protein_bubble_random.png", 
       plot = p, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("function_protein_bubble_random.pdf", 
       plot = p, width = 14, height = 8, bg = "white")
