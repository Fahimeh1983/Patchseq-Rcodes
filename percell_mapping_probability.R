######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/3/2020            #######################################################################
######################################################################################################

# plot percell mapping probabilities 

library(tibble)
library(feather)
library(dplyr)
library(matrixStats)
library(scrattch.hicat)
library(cowplot)
library(feather)
library(dendextend)
library(scrattch.vis)
library(scrattch.io)
library("cowplot")
library(ggbeeswarm)
library(reshape2)
.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")


######################################################################################################
### Input path: ######################################################################################
######################################################################################################
work.dir <- "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019"
memb.path <- paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/",
                    "mouse_patchseq_VISp_20200318_collapsed40_cpm/")

locked.data <- paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/",
                      "revised_manuscript_specimen_ids.txt")

facs.dend.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/",
                        "star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData")

facs.anno.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/anno.feather")

batch_date="20200318_BT014-RSC-252"

######################################################################################################
### Read input: ######################################################################################
######################################################################################################

FACS_dend <- readRDS(facs.dend.path)

patch.anno <- read_feather(paste0(memb.path, "/anno.feather"))
locked.spec.id <- read.csv(locked.data, header = FALSE)$V1
locked.sample.id <- patch.anno[patch.anno$spec_id_label %in% locked.spec.id, "sample_id"]$sample_id

#Read facs.anno and kepp only neuronal part
facs.anno <- read_feather(facs.anno.path)
facs.anno <- facs.anno[facs.anno$class_label %in% c("GABAergic", "Glutamatergic"),]

# Read memb and some cosmetic changes
memb <- read.csv(paste0(memb.path, "/mapping.memb.with.bp.40.csv"))
memb = memb[memb$X %in% locked.sample.id, ]
rownames(memb) <- memb$X
memb <- memb[,setdiff(colnames(memb), "X")]
temp <- gsub("2.3", "2/3", colnames(memb), fixed=TRUE)
temp <- gsub(".", " ", temp, fixed=TRUE)
colnames(memb) <- temp
cls <- temp[(temp %in% labels(FACS_dend))]
memb <- memb[, cls]

exc_types = unique(facs.anno[facs.anno$class_label == "Glutamatergic", "cluster_label"]$cluster_label)
Inh_types = unique(facs.anno[facs.anno$class_label == "GABAergic", "cluster_label"]$cluster_label)
Excitatory_cell_types <- rowSums(memb[,colnames(memb)[colnames(memb) %in% exc_types]])

memb <- cbind(memb[,colnames(memb)[colnames(memb) %in% Inh_types]], Excitatory_cell_types)
dim(memb)
######################################################################################################
### creat a long version of memb: ####################################################################
######################################################################################################


memb <- rownames_to_column(memb)
df <- melt(memb)
df <- df[df$value > 0,] # Keep only the values more than zero
colnames(df) <- c("sample_id", "cluster_label", "mapping_prob")

# The final tree call for each cell is in topleaf label column so we need to add that to our data frame
Topleaf_label = as.data.frame(patch.anno[,c("sample_id", "topLeaf_label")])

# Also we need the dendcluster id for plotting purposes, that should be taken from FACS data
cl.df <- as.data.frame(unique(facs.anno[,c("dendcluster_id", "cluster_label")]))
colnames(cl.df) <- c("dendcluster_id", "topLeaf_label")
cl.df

# Now we add topleaf label and dendcluster id to the data frame
df <- left_join(df, Topleaf_label)
df <- left_join(df, cl.df, by= "topLeaf_label")
df <- df[order(df$dendcluster_id),] #order it by the dendcluster_id

######################################################################################################
### core, i1, i2, i3 and poorQ cells:#################################################################
######################################################################################################

patch.anno <- patch.anno[patch.anno$sample_id %in% locked.sample.id, ]
Core <- patch.anno[patch.anno$Tree_call_label %in% c("Core"), "sample_id"]$sample_id
I1 <- patch.anno[patch.anno$Tree_call_label %in% c("I1"), "sample_id"]$sample_id
I2 <- patch.anno[patch.anno$Tree_call_label %in% c("I2"), "sample_id"]$sample_id
I3 <- patch.anno[patch.anno$Tree_call_label %in% c("I3"), "sample_id"]$sample_id
PoorQ <- patch.anno[patch.anno$Tree_call_label %in% c("PoorQ"), "sample_id"]$sample_id

#load("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/memb.df.rda")

highly_con <- df[df$sample_id %in% c(Core, I1), ]
moderate_con <- df[df$sample_id %in% c(I2,I3), ]
incon <- df[df$sample_id %in% PoorQ, ]

highly_con$sample_id <- as.numeric(factor(highly_con$sample_id))
moderate_con$sample_id <- as.numeric(factor(moderate_con$sample_id))
incon$sample_id <- as.numeric(factor(incon$sample_id))

# For each individual df, we need to sort the cells based on their topleaf label to have a
# diagonal shape for plotting purposes

sort_cells <- function (which.df) {
  
  tmp_sample_id <- c()
  tmp_id <- c()
  track_sample_id <- c()
  count_i = 1
  
  for (t in labels(FACS_dend)) {
    tmpdf = which.df[which.df$topLeaf_label == t, ]
    tmpdf = tmpdf[order(tmpdf$mapping_prob, decreasing = TRUE),]
    for (i in tmpdf$sample_id){
      track_sample_id <- c(track_sample_id, i)
      if (i %in% tmp_sample_id){
        idx = match(i, tmp_sample_id)
        tmp_id <- c(tmp_id, idx)
      }else{
        tmp_sample_id <- c(tmp_sample_id, i)
        idx = count_i
        tmp_id <- c(tmp_id, idx)
        count_i = count_i + 1
      }
    }
  }
  df_reordered <- cbind.data.frame(track_sample_id, tmp_id)
  colnames(df_reordered) <- c("sample_id", "id")
  which.df <- unique(left_join(which.df, df_reordered))
  
  return(which.df)
}

highly_con <- sort_cells(highly_con)
moderate_con <- sort_cells(moderate_con)
incon <- sort_cells(incon)

moderate_con$id <- moderate_con$id + length(unique(highly_con$id)) 
incon$id <- incon$id + length(unique(moderate_con$id)) +  length(unique(highly_con$id)) 

all <- rbind.data.frame(highly_con, moderate_con, incon)



ggplot2::ggplot(all, ggplot2::aes(x = cluster_label, y = id)) +
  ggplot2::geom_point(ggplot2::aes(size = sqrt(mapping_prob), color = mapping_prob))+
  ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.1, hjust = 1, angle = 90, size = 7)) +
  ggplot2::scale_color_gradient(low = "yellow", high = "darkblue") + scale_size_area(max_size = 3,
                                                                                   breaks = c(0.25, 0.5, 0.75, 1))+
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000)) + ylab("Cell ids") +xlab("Cell types")


p <- ggplot(all, aes(cluster_label, id)) + 
  geom_tile(aes(fill = mapping_prob)) + 
  scale_fill_gradient(low = " orange", high = "darkblue")+
  ggplot2::theme(axis.text.x = ggplot2::element_text(vjust = 0.1, hjust = 1, angle = 90, size = 12)) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 17)) +
  scale_y_continuous(breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000)) + 
  ylab("Cell ids") +xlab("Cell types")

p

save_plot(paste0(work.dir, "/Figures/memb_heatmap_June25.pdf"), p, base_width = 20, base_height = 30, dpi = 700)

                                                                                          
