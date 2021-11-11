#fig3 other data
#calculate different tissue percent
t <- table(seurat.other$cell_label, seurat.other$sample_type)
t <- as.matrix(t)
t <- as.data.frame(t)
library(reshape2)
t <- dcast(t, Var1~Var2)
rownames(t) <- t$Var1
t <- t[,-1]

t3 <- apply(t, 1, function(x) x/t2 *100)
t3 <- t(t3)
t4 <- apply(t3, 1, function(x) x/mean(x))
t4 <- t(t4)
t4 <- melt(t4)

t$erichment <- t4$value

flag <- grepl("DC-|Macro-|Mast",rownames(t4)) & !grepl("SCG",rownames(t4))

t4_s <- t4[flag,]
t <- t[flag,]

t <- data.frame(cell_label = rownames(t), t)
t_melt <- melt(t)


t4_s <- data.frame(cell_label = rownames(t4_s), t4_s)
t4_melt <- melt(t4_s)

t_melt$percent <- t4_melt$value

ggplot(t_melt, aes(variable, factor(cell_label, levels = rev(color_map_s$cell_label)))) + geom_point(aes(color = percent, size = value)) +
  scale_colour_gradient(low = 'white', high = "#e76f51") + theme_classic() + theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(hjust = 1, angle = 90),
    legend.position = "left"
  ) + scale_y_discrete(position = "right")

