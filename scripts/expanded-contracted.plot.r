# ============================================================
# 分组条形图 - 按新指定顺序，扩张在上，收缩在下
# ============================================================

rm(list = ls())

library(ggplot2)
library(dplyr)

# 直接创建数据
data <- data.frame(
  group = c(
    "expanded", "expanded", "expanded", "expanded", "expanded",
    "expanded", "expanded", "expanded", "expanded", "expanded",
    "expanded",
    "contracted", "contracted", "contracted", "contracted", "contracted",
    "contracted", "contracted", "contracted", "contracted", "contracted"
  ),
  Class = c(
    "DNA repair and recombination",
    "DNA repair and recombination",
    "Transposon and retrotransposon-related elements",
    "Transposon and retrotransposon-related elements",
    "Transposon and retrotransposon-related elements",
    "Carbohydrate Metabolism",
    "Carbohydrate Metabolism",
    "Signal transduction",
    "Transports",
    "Transports",
    "Transports",
    "Transporters",
    "Transporters",
    "Transporters",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic",
    "Secondary metabolite biosynthetic"
  ),
  Description = c(
    "ATP-dependent DNA helicase",
    "Bifunctional 3'-5' exonuclease/ATP-dependent helicase",
    "Retrovirus-related Pol polyprotein from transposon",
    "Transposon Ty3-I Gag-Pol polyprotein",
    "Transposon Tf2-9 polyprotein",
    "Endo-1,4-beta-xylanase",
    "glucosidase",
    "Serine/threonine-protein kinase",
    "Lactose permease",
    "Maltose permease",
    "ABC multidrug transporter",
    "MFS-type transporter",
    "ABC-type transporter",
    "Uncharacterized transporter",
    "Nonribosomal peptide synthetase",
    "Polyketide synthase",
    "Cytochrome P450 monooygenase",
    "Methyltransferase",
    "FAD-dependent monooxygenase",
    "FAD-binding monooxygenase",
    "ABC multidrug transporter"
  ),
  number = c(13, 4, 28, 6, 14, 6, 7, 7, 12, 6, 3, 28, 4, 4, 4, 13, 15, 9, 3, 3, 9),
  stringsAsFactors = FALSE
)

# 为收缩的数量添加负值（扩张为正，收缩为负）
data_plot <- data %>%
  mutate(number_plot = ifelse(group == "expanded", number, -number))

# ============================================================
# 你指定的视觉顺序（从上到下）
# ============================================================
visual_order_from_top_to_bottom <- c(
  "ATP-dependent DNA helicase",
  "Bifunctional 3'-5' exonuclease/ATP-dependent helicase",
  "Retrovirus-related Pol polyprotein from transposon",
  "Transposon Ty3-I Gag-Pol polyprotein",
  "Transposon Tf2-9 polyprotein",
  "Endo-1,4-beta-xylanase",
  "glucosidase",
  "Serine/threonine-protein kinase",
  "Lactose permease",
  "Maltose permease",
  "ABC multidrug transporter",
  "MFS-type transporter",
  "ABC-type transporter",
  "Uncharacterized transporter",
  "Nonribosomal peptide synthetase",
  "Polyketide synthase",
  "Cytochrome P450 monooygenase",
  "Methyltransferase",
  "FAD-dependent monooxygenase",
  "FAD-binding monooxygenase"
)

# 关键：因为 coord_flip() 会反转显示顺序
# 所以因子顺序 = 反转后的视觉顺序
factor_order <- rev(visual_order_from_top_to_bottom)

cat("===== 顺序说明 =====\n")
cat("视觉顺序（从上到下，共", length(visual_order_from_top_to_bottom), "个）：\n")
for(i in seq_along(visual_order_from_top_to_bottom)) {
  if(i <= 10) {
    type <- "扩张"
  } else if(i == 11) {
    type <- "中间（同时存在）"
  } else {
    type <- "收缩"
  }
  cat("  ", i, ". ", visual_order_from_top_to_bottom[i], " [", type, "]\n")
}

cat("\n因子顺序（第一个显示在最底部，最后一个显示在最顶部）：\n")
for(i in seq_along(factor_order)) {
  cat("  ", i, ". ", factor_order[i], "\n")
}

# 应用顺序
data_plot$Description <- factor(data_plot$Description, levels = factor_order)

# 设置Class的顺序
class_order <- c(
  "DNA repair and recombination",
  "Transposon and retrotransposon-related elements",
  "Carbohydrate Metabolism",
  "Signal transduction",
  "Transports",
  "Transporters",
  "Secondary metabolite biosynthetic"
)
data_plot$Class <- factor(data_plot$Class, levels = class_order)

# ============================================================
# 配色：扩张 #D480B2，收缩 #569BBF
# ============================================================
group_colors <- c(
  "expanded" = "#D480B2",
  "contracted" = "#569BBF"
)

# 获取数量的最大值
max_num <- max(abs(data_plot$number_plot))
x_limit <- max_num * 1.2

# ============================================================
# 绘图
# ============================================================
p <- ggplot(data_plot, aes(x = Description, y = number_plot, fill = group)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(aes(label = number,
                hjust = ifelse(number_plot > 0, -0.2, 1.2)),
            size = 3.2) +
  scale_fill_manual(values = group_colors, name = "") +
  coord_flip() +
  scale_y_continuous(
    name = "Number of genes/proteins",
    limits = c(-x_limit, x_limit),
    breaks = seq(-max_num, max_num, by = 5),
    labels = function(x) abs(x)
  ) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 8, hjust = 1, color = "black"),
    axis.text.x = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11, face = "plain"),
    legend.position = "top",
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 12, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "gray40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "gray85", linewidth = 0.2),
    panel.border = element_rect(color = "black", linewidth = 0.5),
    plot.background = element_rect(fill = "white"),
    panel.background = element_rect(fill = "white")
  ) +
  labs(
    title = "Gene family expansion and contraction",
    subtitle = "Expanded (above, #D480B2) | Contracted (below, #569BBF)",
    x = "",
    y = "Number of genes/proteins"
  )

# 显示并保存
print(p)

ggsave("diverging_bar_final.png", plot = p, width = 12, height = 10, dpi = 300, bg = "white")
ggsave("diverging_bar_final.pdf", plot = p, width = 12, height = 10, bg = "white")

cat("\n===== 绘图完成！=====\n")
cat("文件保存为: diverging_bar_final.png\n")
