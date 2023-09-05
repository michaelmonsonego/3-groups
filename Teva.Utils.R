library(ggpubr)
library(ggrepel)
colors <- c("#CFBAE1", "#5C7D9D","#A4DEF9")


interaction_barplot <- function(interactions, r, l, cols = colors){
  interactions$dataset <- factor(interactions$dataset, levels = c('control', 'treatment'))
  
  aa <- interactions %>% filter(receptor == r) %>% filter(ligand == l)
  aa <- aa %>% group_by(dataset) %>% dplyr::summarise(prob = sum(prob))
  
  if(nrow(aa) == 0){
    return ("Invalid")
  }
  
  p <- ggplot(aa, aes(x=dataset, y= prob, fill = dataset))+
    geom_col(width = 0.8) +
    ggtitle(paste0(l, " - ", r, " interaction")) +
    xlab("") +
    ylab("Interaction Probability") + 
    scale_fill_manual(name = "Sample", labels = c("Saline", "\u03b1PD1(RMP1-30)-IL2"), values = cols, limits = levels(interactions$dataset)) +
    scale_x_discrete(labels = c("Saline", "\u03b1PD1(RMP1-30)-IL2")) +
    theme_classic() +
    theme(legend.position = "none", 
          axis.text.x= element_text(size = 14), 
          axis.title.y = element_text(size = 14))
  
  return(p)
}



bar_plot <- function(file_name, object, choice, dodged = F, intercept = F, geom = "bar"){
  # choice = 0 for count & normalized counts
  # choice = 1 for normalized counts
  # intercept for vline=50%
  
  frque <- table(Idents(object), object$Sample)
  frque <- as.data.frame(frque)
  frque <- frque %>% 
    dplyr::rename(Cluster = Var1, Sample = Var2)
  
  frque <- ddply(frque, .(Sample), transform, percentperstatus = Freq/sum(Freq))
  frque <- ddply(frque, .(Sample), transform, percentperstatus = Freq/sum(Freq))
  frque <- ddply(frque, .(Cluster), transform, percentpercluster = Freq/sum(Freq))
  frque <- ddply(frque, .(Cluster), transform, normalizedppstatus = percentperstatus/sum(percentperstatus)*100)
  frque <- ddply(frque, .(Sample), transform, normalizedpcluster = percentpercluster/sum(percentpercluster))
  
  frque$Cluster <- fct_reorder(frque$Cluster, frque$Freq, .desc = T)
  
  df.median = frque %>% 
    group_by(Cluster) %>% 
    mutate(ymedian = median(normalizedppstatus))
  
  
  p1 <- ggplot(frque, aes(x= Cluster, y= normalizedppstatus, fill = Sample)) +
    theme_classic() +
    scale_fill_manual(values =colors) +
    # geom_bar(stat = "identity", width = 0.7) + 
    ggtitle("Normalized Cluster Counts") +
    ylab("Normalized Proportions in %") + xlab(NULL) +
    theme(
      axis.text.x = element_text(size = 11, color="black",angle=330, hjust = 0),
      axis.text.y = element_text(size = 11, color="black"),
      axis.title.x = element_text(size = 13, color="black"),
      axis.title.y = element_text(size = 13, color="black"),
      legend.text = element_text(size = 13, color="black"),
      title = element_text(size = 13),
      legend.justification = "top",
      panel.border = element_blank()) 
  
    # geom_hline(yintercept = 50, linetype='dotted') +
    # annotate("text", x = levels(object@active.ident)[length(levels(object@active.ident))], y = 50, label = "50%",
    # vjust = -0.5)
  
  if(geom == "bar"){
    if(dodged == T){
      p1 <- p1 + geom_bar(stat = "identity", width = 0.7, position = "dodge")
    }
    else{
      p1 <- p1 + geom_bar(stat = "identity", width = 0.7) 
    }
  }
  
  else if (geom == "point"){
    p1 <- p1 + geom_jitter(size = 10, width = 0.2, aes(color = Sample)) +
      geom_errorbar(data=df.median, aes(Cluster, ymax = ymedian, ymin = ymedian),
                    size=0.5, inherit.aes = F, width = 0.8) +
      scale_color_manual(values =colors) +
      geom_vline(xintercept = seq(1.5,length(levels(frque$Cluster))-0.5),linetype='dashed',color='gray')
  }
  else{ return("geom can be bar/point")}

  
  if(intercept == T){
    p1 <- p1 + geom_hline(yintercept = 50, linetype='dotted') +
      annotate("text", x = levels(object@active.ident)[length(levels(object@active.ident))], y = 50, label = "50%",
      vjust = -0.5)
  }
  
  
  p2<-ggplot(frque, aes(x= Cluster, y= Freq, fill = Sample)) +
    theme_classic() +
    scale_fill_manual(values =colors) +
    geom_bar(stat = "identity", width = 0.7) +
    ggtitle("Cluster Counts") + ylab("Number of cells") + xlab(NULL) +
    theme(
      axis.text.x = element_text(size = 11, color="black",angle=330, hjust = 0),
      axis.text.y = element_text(size = 11, color="black"),
      axis.title.x = element_text(size = 13, color="black"),
      axis.title.y = element_text(size = 13, color="black"),
      title = element_text(size = 13),
      legend.position = "none",
      panel.border = element_blank())
  
  
  if(choice == 0){
    p3 <- grid.arrange(p2,p1, nrow= 1, widths=c(0.4, 0.6))
    ggsave(file_name, p3,width = 11, height = 5, dpi = 300)
    # ggsave(file_name, p1,width = 7, height = 6)
  }
  else if (choice == 1){
    ggsave(file_name, p1,width = 7, height = 6, dpi = 300)
  }
}

# example 
# bar_plot("C:\\Users\\tomer\\Desktop\\Master\\Projects\\a.png", merged, 1)

wrap_text <- function(x, chars = 20) {
  stringr::str_wrap(x, chars)
}

Vln_Plot <- function(object, feature, my_comparisons, idents = object@active.ident, facet = 'up', box = F, stat = T, assay = 'RNA'){
  df <- object@assays[[assay]]@data %>% as.data.frame() 
  max <- df %>% filter(row.names(df) %in% c(feature)) %>% max()
  
  if(length(levels(object$Sample)) == 2){
    colors <- colors[c(1,3)]
    max.y = max *1.15
  }
  else{
    max.y = max *1.3
  }
  
  
  # print(levels(object$Sample))
  
  
  p <- VlnPlot(object, features = feature, pt.size = 0, group.by = 'Sample', 
               cols = colors, y.max = max.y, assay = assay) +
    xlab(NULL)
  if(facet == 'up'){
    p <- p+ facet_wrap(idents, nrow = 1) 
  }
  else{
    p <- p+ facet_wrap2(idents, nrow = 1, strip.position = "bottom", labeller = as_labeller(wrap_text), 
                        strip = strip_vanilla(clip = "off")) +
      theme(strip.text.x = element_text(angle=330, hjust = 0), 
            strip.placement = "outside",
            strip.background = element_blank())
  }
  
  if(box == T){
    p <- p + geom_boxplot(outlier.shape = NA, width = .4, alpha = .3, show.legend = FALSE, coef = 0, 
                          position = position_dodge(width = .9))
  }
  
  if(stat == T){
    p <- p + stat_compare_means(aes(group = 'Sample'), comparisons = my_comparisons, label = "p.signif", 
                                method = "wilcox.test", size = 6, bracket.size = 0.8, hide.ns = F) +
      theme(axis.text.x = element_blank())
  }
  
  if(assay == 'sig'){
    p <- p +   scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, by = 0.25)) +
      ylab("Signature Score") +
      theme(legend.position="none")  
  }
  
  return (p)
}

# example
# comparisons <- list(levels(CD3subset$Sample))
# Vln_Plot(CD3subset, feature = 'Havcr2', my_comparisons = comparisons, idents = CD3subset$Annotation)


# Sys.setlocale(category = "LC_ALL", locale = "English")




Volcano_plot <- function(table, highlight_genes = c(), top_n = 10, pajd_thr = 0.05, LFC_cutoff = 0.25, title, subtitle = "Atteunukine vs Saline"){
  table <- table %>%
    mutate(geneLabel = ifelse(p_val_adj <= pajd_thr, gene, NA))
  
  
  table <- table %>%
    mutate(diffexpressed = 
             ifelse(p_val_adj <= 0.05 & avg_log2FC >= LFC_cutoff, "UP", 
                    ifelse(p_val_adj <= 0.05 & avg_log2FC <= -LFC_cutoff, "DOWN", "NO")))
  
  
  table$highlighted <- NA
  table$highlighted[which(table$gene %in% highlight_genes)] <- table$gene[which(table$gene %in% highlight_genes)]
  
  for(gene in highlight_genes){
    if(!gene %in% table$gene){
      print(paste0(gene, " was not found in table"))
    }
  }
  

  p <- ggplot(data=table, aes(x = avg_log2FC, y=-log10(p_val_adj), 
                                  col = diffexpressed, label =  geneLabel)) + 
    theme_classic() +
    geom_point() +
    scale_color_manual(values=c("DOWN" = "#1d5cf0", "NO" = "black", "UP" = "#c70202"), na.translate = F) +
    geom_vline(xintercept = LFC_cutoff, linetype = "dashed") +
    geom_vline(xintercept = -LFC_cutoff, linetype = "dashed") +
    geom_hline(yintercept = -log10(pajd_thr), linetype = "dashed") +
    xlab("Avg_log2FC") + 
    ylab("-Log10(padj)") +
    ggtitle(title) +
    labs(subtitle = subtitle) +
    theme(axis.line = element_line(),
          legend.position="none", 
          plot.title = element_text(hjust = 0),
          axis.title = element_text(size = 10)) 
  
  if(length(highlight_genes) >0){
    p <- p + geom_label_repel(aes(label = highlighted), max.overlaps = Inf, box.padding = 0.2, fill = "white",color = "black", 
                     arrow = arrow(length = unit(0.015, "npc")),  show.legend=F)
  }
  else{
    p <- p + geom_label_repel(aes(label = geneLabel), max.overlaps = Inf, box.padding = 0.2, fill = "white",color = "black", 
                     arrow = arrow(length = unit(0.015, "npc")),  show.legend=F)
  }
  

    
    
    return(p)
}

  
  
  
  
  
  