######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/3/2020            #######################################################################
######################################################################################################

#This code is to generate FACS-Patchseq side by side barplot and layer distribution

Build_layer_barplot_perlayer_data <- function(facs.anno,
                                              patch.anno,
                                              dend,
                                              dendcluster_ids,
                                              ypad = 0.05, 
                                              xpad = 0.1){
  
  facs.anno <- as.data.frame(facs.anno)
  patch.anno <- as.data.frame(patch.anno)
  keep_layers <- c("L1", "L2/3", "L4", "L5", "L6")
  
  patch.anno$layer_label <- patch.anno$Revisited_layer_label
  patch.anno$layer_id <- patch.anno$Revisited_layer_id
  patch.anno$layer_color <- patch.anno$Revisited_layer_color
  
  filtered.facs.anno <- facs.anno %>%
    filter(dendcluster_id %in% dendcluster_ids) 
  
  filtered.facs.anno <- filtered.facs.anno %>%
    filter(layer_label %in% keep_layers) 
  
  filtered.patch.anno <- patch.anno %>%
    filter(dendcluster_id %in% dendcluster_ids)
  
  #Layer range is the same for both patch and facs
  facs.layer_ranges <- data.frame(layer_label = rev(keep_layers),
                                  ymin = seq(1, 3, by=0.5) -1 + ypad,
                                  ymax = seq(1, 3, by=0.5) -0.5 - ypad) %>% mutate(ymid = (ymin + ymax)/2) 
  
  
  facs.cluster_ranges <- filtered.facs.anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = 1:n() - 1 + xpad,
           xmax = 1:n()     - xpad,
           xmid = 1:n() - 0.5)
  
  facs.layer_rects <- data.frame()
  i = 1
  for (l in unique(keep_layers)) {
    print(l)
    tmp <- facs.cluster_ranges %>% select(cluster_label,xmin, xmax) %>% 
      mutate(xmin = facs.cluster_ranges$xmin , 
             xmax = facs.cluster_ranges$xmax ,
             layer_label = l) %>% 
      left_join(facs.layer_ranges)
    tmp$fill <- rep(c("#E0E0E0", "#FFFFFF"), dim(facs.cluster_ranges)[1])[1:dim(facs.cluster_ranges)[1]]
    facs.layer_rects <- rbind.data.frame(facs.layer_rects, tmp)
    i <- i+1
  }
  
    # Layer rectangles
  facs.plot_data <- filtered.facs.anno %>%
    select(cluster_id, cluster_label, cluster_color, layer_id, layer_label, layer_color) %>%
    left_join(facs.layer_rects[,c("cluster_label", "xmin", "xmax", "layer_label", "ymin", "ymax")]) %>% 
    group_by(cluster_id, layer_id, layer_label) %>%
    mutate(ly_n = n()) %>%
    ungroup() %>%
    group_by(cluster_id) %>%
    arrange(layer_id) %>%
    mutate(cluster_n = n(),
           ly_frac = ly_n/cluster_n) %>%
    unique() %>%
    arrange(cluster_id, layer_id) %>%
    group_by(cluster_id) %>%
    mutate(xmin_bar = xmin,
           xmax_bar = xmin + 0.4,
           ymax_bar = ymin + (ymax - ymin) * ly_frac,
           ymin_bar = ymin)
  
  patch.plot_data <- filtered.patch.anno %>%
    select(cluster_id, cluster_label, cluster_color, layer_id, layer_label, layer_color) %>%
    left_join(facs.layer_rects[,c("cluster_label", "layer_label", "xmin", "xmax", "ymin", "ymax")]) %>% 
    group_by(cluster_id, layer_id, layer_label) %>%
    mutate(ly_n = n()) %>%
    ungroup() %>%
    group_by(cluster_id) %>%
    arrange(layer_id) %>%
    mutate(cluster_n = n(),
           ly_frac = ly_n/cluster_n) %>%
    unique() %>%
    arrange(layer_id) %>%
    mutate(ly_cum_frac = cumsum(ly_frac)) %>%
    ungroup() %>%
    arrange(cluster_id, layer_id) %>%
    group_by(cluster_id) %>%
    mutate(xmin_bar = xmin + 0.4,
           xmax_bar = xmax,
           ymax_bar = ymin + (ymax - ymin) * ly_frac,
           ymin_bar = ymin)
  
  prune_dend_labels <- labels(dend)[!labels(dend) %in% facs.plot_data$cluster_label]
  library(scrattch.hicat)
  filtered_dend <- prune_dend(dend = dend, rm.labels = prune_dend_labels)
  
  dend_seg <- as.ggdend(filtered_dend)$segments %>%
    mutate(y = (y/max(y))*3 + max(facs.layer_rects$ymax) + ypad,
           yend = (yend/max(yend))*3 + max(facs.layer_rects$ymax) + ypad,
           x = x - 0.5,
           xend = xend - 0.5)
 
  shrink_value =  0.15 
  flat_plot <- ggplot() +
    geom_segment(data = dend_seg,
                 aes(x = x, xend = xend,
                     y = y *shrink_value + -1, yend = (yend*shrink_value) + -1,
                     color = "Black"),
                 lineend = "square") +
    geom_point(data = patch.plot_data,
               aes(x =  (xmax + xmin)/2 , 
                   y = -0.75,
                   size = cluster_n,
                   fill = cluster_color),
               color= patch.plot_data$cluster_color,
               pch = 21) + 
    geom_text(data = patch.plot_data,
              aes(x = (xmax + xmin)/2 + 0.2,
                  y = -1.5,
                  label = cluster_n,
                  color = "Black"),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 3) +
    geom_text(data = facs.plot_data,
              aes(x = (xmax + xmin)/2 -0.2,
                  y = -1.5,
                  label = cluster_n,
                  color = "Black"),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 3) +
    geom_rect(data = facs.layer_rects,
              aes(xmin = xmin ,
                  xmax = xmax ,
                  ymin = ymin,
                  ymax = ymax,
                  fill = fill)) +
    geom_rect(data = facs.plot_data,
              aes(xmin = xmin_bar ,
                  xmax = xmax_bar-0.05 ,
                  ymin = ymin_bar,
                  ymax = ymax_bar,
                  fill = facs.plot_data$cluster_color)) +
    geom_rect(data = patch.plot_data,
              aes(xmin = xmin_bar ,
                  xmax = xmax_bar -0.05,
                  ymin = ymin_bar,
                  ymax = ymax_bar,
                  fill = patch.plot_data$cluster_color)) +
    geom_rect(data = facs.cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2 ,
                  ymax = 5,
                  ymin = 2.6),
              fill = "#E0E0E0")+
    geom_text(data = facs.cluster_ranges,
              aes(x = xmid,
                  y = 2.6 +  ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 3) +
    geom_rect(data = facs.cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2 ,
                  ymax = 5.1 - ypad,
                  ymin = 5.1 - 2 * ypad),
              fill = facs.cluster_ranges$cluster_color) +
    geom_rect(data = facs.cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2 ,
                  ymax = 2.65 - ypad,
                  ymin = 2.65 - 2 * ypad),
              fill = facs.cluster_ranges$cluster_color) +
    scale_size_area(max_size = 5,
                    breaks = c(1,10,50,100,200)) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_y_continuous(expand = c(0.25,0)) +
    theme_void()
  flat_plot
}

group_violin_FACS_patch_plot  <- function (genes = c("Hspa8", "Snap25", "Gad2", "Vip"), group_by = "final", 
                                           clusters = 1:10,  sort = F, logscale = F, 
                                           showcounts = T, rotatecounts = F, fontsize = 7, labelheight = 25, 
                                           max_width = 10, REF.anno.feather = REF.anno.feather, REF.data.feather = REF.data.feather, 
                                           map.anno.feather = map.anno.feather, map.data.feather = map.data.feather,dend = dend ) 
{
  library(dplyr)
  library(ggplot2)
  
  genes <- rev(genes)
  
  tmp <- get_violin_data(anno.feather = REF.anno.feather, data.feather = REF.data.feather, 
                         FACS.or.patch = "FACS", dend = dend, genes = genes, logscale = logscale, group_by = group_by, clusters = clusters)
  REF.data <- tmp$data
  REF.max_vals <- tmp$max_vals
  tmp <- get_violin_data(anno.feather = map.anno.feather, data.feather = map.data.feather, 
                         FACS.or.patch = "Patch", dend = dend, genes = genes, logscale = logscale, group_by=group_by, clusters= clusters)
  map.data <- tmp$data
  map.max_vals <- tmp$max_vals
  
  genes <- tmp$genes
  ngenes <- length(genes)
  nclust <- length(union(unique(REF.data$plot_id), unique(map.data$plot_id)))
  
  header_labels <- build_header_labels(data = REF.data, grouping = "plot", 
                                       ymin = ngenes + 1,
                                       label_height = labelheight, 
                                       label_type = "simple")
  
  max_val <- cbind(REF.max_vals, map.max_vals)
  max_val <- apply(max_val, 1, max)
  max_val <- log10(max_val)
  max_val <- format(round(max_val, 2), nsmall = 2)

  max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
                                 0.5, label = max_val)

  max_header <- data.frame(x = (nclust + 0.5) * 1.01, y = ngenes + 
                             1, label = "Max value")
  max_width <- nclust * (max_width/100)/(1 - max_width/100)
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  
  cluster_data <- REF.data %>% group_by(plot_label, plot_color, 
                                        plot_id) %>% summarise(cn = n()) %>% as.data.frame(stringsAsFactors = F) %>% 
    arrange(plot_id) %>% mutate(labely = ngenes + label_y_size * 
                                  0.05, cny = max(header_labels$ymax) - 0.1 * label_y_size, 
                                xpos = plot_id)
  
  vertical_background_rect <- header_labels %>% mutate(ymin = 1,
                                                       ymax = length(genes) + 1,
                                                       xmin = seq(0.5,nclust*4,2)[1:nclust],
                                                       xmax = seq(2.5,nclust*4,2)[1:nclust], 
                                                       color = rep(c("#E0E0E0", "#FFFFFF"), dim(header_labels)[1])[1:dim(header_labels)[1]])
  

  p <- ggplot() +
    scale_fill_identity() +
    scale_y_continuous("",
                       breaks = 1:length(genes) + 0.45,
                       labels = genes,
                       expand = c(0,0)) +
    scale_x_continuous("", expand = c(0, 0, 0, 5)) +
    theme_classic(fontsize) +
    theme(axis.text = element_text(size = rel(1)),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(face = "italic", color = "#000000"),
          legend.position = "none",
          axis.line = element_line(size = 0.1)) +
    geom_rect(data = vertical_background_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = color)) +
    geom_rect(data=header_labels, 
              aes(xmin = xmin, 
                  xmax = xmax + 1,
                  ymin = min(vertical_background_rect$ymax) , 
                  ymax = min(vertical_background_rect$ymax) + 0.5,
                  fill = color))
  
 p 
 
  for (i in 1:length(genes)) {
    p <- p + geom_violin(data = REF.data, 
                         aes_string(x = "xpos", 
                                    y = genes[i],
                                    group = "xpos",
                                    fill = "cluster_color"),
                         #fill = "#0000CC",
                         scale = "width",
                         width = 0.9, 
                         size = 0.05,
                         adjust = 2) +
      stat_summary(data = REF.data, aes_string(x = "xpos", 
                                               y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1) +
      geom_violin(data = map.data, 
                  aes_string(x = "xpos", 
                             y = genes[i],
                             group = "xpos",
                             fill = "cluster_color"),
                  #fill = "#FF0000",
                  scale = "width",
                  width = 0.9, 
                  size = 0.05,
                  adjust = 2) +
      stat_summary(data = map.data, aes_string(x = "xpos", 
                                               y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1) 

  }

   p <- p + geom_text(data = max_labels, 
                     aes(x = 2*x  , 
                         y = y, 
                         label = as.character(label)),
                     color = "Black",
                     hjust = 0, 
                     vjust = 0.65, 
                     size = pt2mm(fontsize), 
                     parse = TRUE) 
  
  if (showcounts) {
    if (rotatecounts) {
      p <- p + geom_text(data = cluster_data, aes(y = cny,
                                                  x = xpos, label = cn), angle = 90, vjust = 0.35,
                         hjust = 1, size = pt2mm(fontsize))
    }
    else {
      p <- p + geom_text(data = cluster_data, aes(y = cny,
                                                  x = xpos, label = cn), size = pt2mm(fontsize))
    }
  }
  p
}

get_violin_data <- function(anno.feather = anno.feather, data.feather= data.feather, 
                            FACS.or.patch= "FACS", dend=dend, genes = genes, logscale = logscale, group_by= group_by, clusters = clusters){
  
  anno <- anno.feather %>% dplyr::mutate_if(is.factor, as.character)
  if (FACS.or.patch == "FACS") {
    anno <- anno[anno$cluster_label %in% labels(dend) & anno$region_label == "VISp",]
  }
  data <- get_feather_data(anno = anno, data = data.feather, genes = genes, group_by = group_by,
                           group_ids =  clusters) 
  
  genes <- sub("-", ".", genes)
  genes <- genes[genes %in% names(data)]
  if (FACS.or.patch == "FACS") {
    data <- data %>% select(-xpos) %>% mutate(xpos = plot_id + (plot_id -1 ))
  } else {
    data <- data %>% select(-xpos) %>% mutate(xpos = plot_id  * 2)
  }
  
  genes[grepl("^[0-9]", genes)] <- paste0("X", genes[grepl("^[0-9]", genes)])
  names(data)[grepl("^[0-9]", genes)] <- paste0("X", names(data)[grepl("^[0-9]", 
                                                                       genes)])
  
  max_vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% 
    unlist()
  data[genes] <- data[genes] + runif(nrow(data), 0, 1e-05)
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if (logscale) {
      data[gene] <- log10(data[gene] + 1)/log10(gene_max + 1) * 0.9 + i
    }
    else {
      data[gene] <- data[gene]/gene_max * 0.9 + i
    }
  }
  return(list(data = data,max_vals= max_vals, genes = genes))
}


get_feather_data <- function (anno, data, genes, group_by, group_ids)#, dend=NULL) 
{
  id_cols <- names(anno)[grepl("_id$", names(anno)) & names(anno) != 
                           "sample_id"]
  anno[id_cols] <- lapply(anno[id_cols], as.numeric)
  data_names <- names(data)
  if (sum(genes %in% data_names) != length(genes)) {
    not_found <- genes[!toupper(genes) %in% toupper(data_names)]
    warning(paste(paste0(not_found, collapse = ", "), "not found in feather data!"))
    genes <- data_names[toupper(data_names) %in% toupper(genes)]
  }
  data_cols <- which(data_names %in% c("sample_id", genes))
  gene_data <- data[, data_cols]
  colnames(gene_data) <- gsub("-", ".", colnames(gene_data))
  genes <- gsub("-", ".", genes)
  all_anno <- anno %>% 
    dplyr::rename_(plot_id = paste0(group_by,"_id"), plot_label = paste0(group_by, "_label"), 
                   plot_color = paste0(group_by,"_color"))
  cluster_order <- data.frame(group_ids = group_ids) %>% dplyr::mutate(cluster_x = 1:n())
  data <- dplyr::left_join(all_anno, gene_data, by = "sample_id") %>% 
    dplyr::filter(plot_id %in% group_ids) %>% 
    dplyr::left_join(cluster_order, by = c(plot_id = "group_ids")) %>%
    dplyr::arrange(cluster_x) %>% 
    dplyr::mutate(xpos = 1:n()) %>% 
    dplyr::select(-plot_id) %>% 
    dplyr::rename_(plot_id = "cluster_x")
  return(data)
}

