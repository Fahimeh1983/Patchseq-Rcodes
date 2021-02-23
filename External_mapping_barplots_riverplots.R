######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 2/20/2021             ######################################################################
######################################################################################################

# Plotting (barplots, riverplots) M1 tolias and Callaway data

.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")
library(dplyr)
library(ggplot2)
library(cowplot)
library(feather)
library(dendextend)
library(scrattch.vis)
library(scrattch.io)
library(tibble)
library("cowplot")
library(ggbeeswarm)
library(matrixStats)
library(reshape2)
library(scrattch.hicat)
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/",
              "Manuscript_patchseq_2019/Rcodes/Formating_patchanno_helper_functions.R"))
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/",
              "Manuscript_patchseq_2019/Rcodes/Stat_helper_functions.R"))

######################################################################################################
### Input path: ######################################################################################
######################################################################################################

facs.anno.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/anno.feather")

facs.data.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/data.feather")

facs.dend.V1.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/",
                           "star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData")

facs.dend.V1.ALM.labeled.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                                       "patch_seq/star/mouse_patchseq_script_repository/V1_ALM/",
                                       "dend.labeled.rda")

facs.dend.V1.ALM.notlabeled.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                                          "patch_seq/star/mouse_patchseq_script_repository/V1_ALM/",
                                          "neuron_only_dend_V1_ALM_markers_attached.rda")

patch.seq.data.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/",
                             "Mouse/patchseq/R_Object/")

select.markers.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                             "Taxonomies/AIT2.3.1/select.markers.rda")

work.dir = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb",
                  "/EXC_patchseq_paper_2020/")

locked.data.file <- "locked_Feb3.xlsx"

Tolias.data.path = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Tolias_M1_patch/"


Tolias.patch.m1.file = "m1_patchseq_cpm.csv"

Callaway.data.file = "GSE133230_rawCount_945filtered_cells.tsv"

batch_date="20201216_BT014-RSC-263"

######################################################################################################
### Reading ref data, V1+ALM taxonomy and V1 taxonomy ################################################
######################################################################################################
### TODO: the following two lines should be modified if the number of FACS cells has changed 

anno <- read_feather(facs.anno.path)
#checking 
table(anno$class_label)
#removing low quality and non neuronal cells as V1+ALM taxonomy does not have those types
#if you want to include those in the mapping, then we need to remake the tree 
anno <- anno[anno$class_label %in% c("GABAergic", "Glutamatergic"),]

#In total 22439 cells remain
dim(anno)

cl.df <- unique(as.data.frame(anno[,c("cluster_id", "cluster_label", 
                                      "cluster_color", "subclass_id", 
                                      "subclass_label", "class_id", 
                                      "class_label", "cl", "dendcluster_id")]))

dim(cl.df)
cl.df <- cl.df[order(cl.df$dendcluster_id),]

#cl <- set_names(anno[, "cluster_label"]$cluster_label, anno[, "sample_id"]$sample_id)
#cl <- factor(cl)

######################################################################################################
### Loading the patchseq data ########################################################################
######################################################################################################

tmp<-load(paste0(patch.seq.data.path,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(patch.seq.data.path,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),
                      which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),] #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))
patch.seq.norm = query.dat.norm
dim(patch.seq.norm)

######################################################################################################
### river plots Patchseq: ############################################################################
######################################################################################################

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R")

#Read mapping results:
patch.results <- read.csv(paste0(work.dir,"/patch_results_on_v1_and_v1_alm_taxonomies.csv"))
patch.results <- patch.results %>% column_to_rownames("X")

#First lets compare how patchseq data map to V1 vs V1+ALM taxonomy
patch.v1.alm.results <- patch.results %>%
  rownames_to_column("sample_id") %>% 
  select(c("sample_id", "v1_alm_first_cl", "v1_alm_call")) %>%
  rename(cluster_label = v1_alm_first_cl)

patch.v1.alm.results <- merge(cl.df[, c("dendcluster_id", "cluster_color", "cluster_label", "class_label")], patch.v1.alm.results, by = "cluster_label")
colnames(patch.v1.alm.results) <- c("v1_alm_cluster_label", "v1_alm_dendcluster_id", "v1_alm_cluster_color", "v1_alm_class_label", "sample_id", "v1_alm_call" )

patch.v1.results <- patch.results %>% 
  rownames_to_column("sample_id") %>% 
  select(c("sample_id", "v1_first_cl", "v1_call")) %>%
  rename(cluster_label = v1_first_cl)

patch.v1.results <- merge(cl.df[, c("dendcluster_id", "cluster_color", "cluster_label", "class_label")], patch.v1.results, by = "cluster_label")
colnames(patch.v1.results) <- c("v1_cluster_label", "v1_dendcluster_id", "v1_cluster_color", "v1_class_label", "sample_id", "v1_call" )

agg.patch <- merge(patch.v1.alm.results, patch.v1.results, by = "sample_id")
dim(agg.patch)

#library("readxl")
locked.data <- read_excel(paste0(work.dir, locked.data.file))
Specimen_id <- locked.data$`Cell Specimen Id`
locked_ids <- samp.dat[samp.dat$cell_id %in% Specimen_id, "patched_cell_container"]
sum(locked_ids %in% agg.patch$sample_id)


river_plot(agg.patch %>% 
             rename(cluster_id = v1_dendcluster_id, 
                    cluster_label = v1_cluster_label,
                    cluster_color = v1_cluster_color,
                    map_cluster_id = v1_alm_dendcluster_id, 
                    map_cluster_label = v1_alm_cluster_label, 
                    map_cluster_color = v1_alm_cluster_color), 
           min.cells = 0)

sum(as.character(patch.results$v1_alm_first_cl) == as.character(patch.results$v1_first_cl))
table(patch.results[patch.results$v1_alm_call %in% c("Core", "I1", "I2", "I3"), "v1_alm_first_cl"])

######################################################################################################
### river plots Callaway: ############################################################################
######################################################################################################

source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Plot_utils.R")

#Read mapping results:
Callaway.results <- read.csv(paste0(work.dir,"/Callaway_results_on_v1_and_v1_alm_taxonomies.csv"))
Callaway.results <- Callaway.results %>% column_to_rownames("X")

#First lets compare how patchseq data map to V1 vs V1+ALM taxonomy
Callaway.v1.alm.results <- Callaway.results %>%
  rownames_to_column("sample_id") %>% 
  select(c("sample_id", "v1_alm_first_cl", "v1_alm_call")) %>%
  rename(cluster_label = v1_alm_first_cl)

Callaway.v1.alm.results <- merge(cl.df[, c("dendcluster_id", "cluster_color", "cluster_label", "class_label")], Callaway.v1.alm.results, by = "cluster_label")
colnames(Callaway.v1.alm.results) <- c("v1_alm_cluster_label", "v1_alm_dendcluster_id", "v1_alm_cluster_color", "v1_alm_class_label", "sample_id", "v1_alm_call" )

Callaway.v1.results <- Callaway.results %>% 
  rownames_to_column("sample_id") %>% 
  select(c("sample_id", "v1_first_cl", "v1_call")) %>%
  rename(cluster_label = v1_first_cl)

Callaway.v1.results <- merge(cl.df[, c("dendcluster_id", "cluster_color", "cluster_label", "class_label")], Callaway.v1.results, by = "cluster_label")
colnames(Callaway.v1.results) <- c("v1_cluster_label", "v1_dendcluster_id", "v1_cluster_color", "v1_class_label", "sample_id", "v1_call" )

agg.Callaway <- merge(Callaway.v1.alm.results, Callaway.v1.results, by = "sample_id")
dim(agg.Callaway)

river_plot(agg.Callaway[agg.Callaway$v1_alm_call %in% c("Core", "I1", "I2") & agg.Callaway$v1_call %in% c("Core", "I1", "I2"),] %>% 
             rename(cluster_id = v1_dendcluster_id, 
                    cluster_label = v1_cluster_label,
                    cluster_color = v1_cluster_color,
                    map_cluster_id = v1_alm_dendcluster_id, 
                    map_cluster_label = v1_alm_cluster_label, 
                    map_cluster_color = v1_alm_cluster_color), 
           min.cells = 0)

######################################################################################################
### Bar plots:   #####################################################################################
######################################################################################################

locked.data <- read_excel(paste0(work.dir, "/locked_Feb3.xlsx"))
locked_Specimen_id <- locked.data$`Cell Specimen Id`
samp.dat %>% select("cell_id", "patched_cell_container") %>% filter("cell_id" %in% locked_Specimen_id)
locked_sample_id <- samp.dat[samp.dat$cell_id %in% locked_Specimen_id, "patched_cell_container"]

locked.patch.results <- patch.results[locked_sample_id,]

M1.results <- read.csv(paste0(work.dir,"/M1_results_on_v1_and_v1_alm_taxonomies.csv"))
M1.results <- M1.results %>% column_to_rownames("X")

Callaway.results <- read.csv(paste0(work.dir,"/Callaway_results_on_v1_and_v1_alm_taxonomies.csv"))
Callaway.results <- Callaway.results %>% column_to_rownames("X")


# Compare the exc types for M1 and patchseq data on V1_ALM taxonomy
exc_types <- cl.df[cl.df$class_label == "Glutamatergic", "cluster_label"]
counts.M1.locked.patchseq <- rbind(table(droplevels(locked.patch.results[locked.patch.results$v1_alm_first_cl %in% exc_types &
                                                                           locked.patch.results$v1_alm_call %in% c("Core", "I1", "I2", "I3"), c("v1_alm_first_cl")])), 
                                   table(droplevels(M1.results[M1.results$v1_alm_first_cl %in% exc_types &
                                                                 M1.results$v1_alm_call %in% c("Core", "I1", "I2", "I3"), c("v1_alm_first_cl")])))
rownames(counts.M1.locked.patchseq) <- c("Locked_patchseq_data", "M1_Tolias")
ordered_types <- exc_types[exc_types %in% colnames(counts.M1.locked.patchseq)]

par(mar=c(14, 4.1, 4.1, 2.1), xpd=TRUE)
barplot(counts.M1.locked.patchseq[, ordered_types], main="Comparing M1(red) and locked_patchseq data(blue) on V1_ALM taxonomy (Core, I1, I2 and I3 cells)",
        col=c("darkblue","red"), 
        legend = rownames(counts.M1.locked.patchseq), 
        beside=TRUE, 
        legend.text=TRUE,
        las=2,
        args.legend = list(x = "topleft", bty="n"))


# Compare the exc types for M1 and patchseq data on V1 taxonomy
exc_types <- cl.df[cl.df$class_label == "Glutamatergic", "cluster_label"]
counts.M1.locked.patchseq <- rbind(table(droplevels(locked.patch.results[locked.patch.results$v1_first_cl %in% exc_types &
                                                                           locked.patch.results$v1_call %in% c("Core", "I1", "I2", "I3"), c("v1_first_cl")])), 
                                   table(droplevels(M1.results[M1.results$v1_first_cl %in% exc_types &
                                                                 M1.results$v1_call %in% c("Core", "I1", "I2", "I3"), c("v1_first_cl")])))
rownames(counts.M1.locked.patchseq) <- c("Locked_patchseq_data", "M1_Tolias")
ordered_types <- exc_types[exc_types %in% colnames(counts.M1.locked.patchseq)]

par(mar=c(14, 4.1, 4.1, 2.1), xpd=TRUE)
barplot(counts.M1.locked.patchseq[, ordered_types], main="Comparing M1(red) and locked_patchseq data(blue) on V1 taxonomy (Core, I1, I2 and I3 cells)",
        col=c("darkblue","red"), 
        legend = rownames(counts.M1.locked.patchseq), 
        beside=TRUE, 
        legend.text=TRUE,
        las=2,
        args.legend = list(x = "topleft", bty="n"))

# Compare the exc types for Callaway and patchseq data on V1_ALM taxonomy
exc_types <- cl.df[cl.df$class_label == "Glutamatergic", "cluster_label"]
counts.Callaway.locked.patchseq <- rbind(table(droplevels(locked.patch.results[locked.patch.results$v1_alm_first_cl %in% exc_types &
                                                                           locked.patch.results$v1_alm_call %in% c("Core", "I1", "I2", "I3"), c("v1_alm_first_cl")])), 
                                   table(droplevels(Callaway.results[Callaway.results$v1_alm_first_cl %in% exc_types &
                                                                       Callaway.results$v1_alm_call %in% c("Core", "I1", "I2", "I3"), c("v1_alm_first_cl")])))
rownames(counts.Callaway.locked.patchseq) <- c("Locked_patchseq_data", "M1_Tolias")
ordered_types <- exc_types[exc_types %in% colnames(counts.Callaway.locked.patchseq)]

par(mar=c(14, 4.1, 4.1, 2.1), xpd=TRUE)
barplot(counts.Callaway.locked.patchseq[, ordered_types], main="Comparing Callaway(red) and locked_patchseq data(blue) on V1_ALM taxonomy (Core, I1, I2 and I3 cells)",
        col=c("darkblue","red"), 
        legend = rownames(counts.Callaway.locked.patchseq), 
        beside=TRUE, 
        legend.text=TRUE,
        las=2,
        args.legend = list(x = "topleft", bty="n"))


# Compare the exc types for Callaway and patchseq data on V1_ALM taxonomy
exc_types <- cl.df[cl.df$class_label == "Glutamatergic", "cluster_label"]
counts.Callaway.locked.patchseq <- rbind(table(droplevels(locked.patch.results[locked.patch.results$v1_first_cl %in% exc_types &
                                                                                 locked.patch.results$v1_call %in% c("Core", "I1", "I2", "I3"), c("v1_first_cl")])), 
                                         table(droplevels(Callaway.results[Callaway.results$v1_first_cl %in% exc_types &
                                                                             Callaway.results$v1_call %in% c("Core", "I1", "I2", "I3"), c("v1_first_cl")])))
rownames(counts.Callaway.locked.patchseq) <- c("Locked_patchseq_data", "M1_Tolias")
ordered_types <- exc_types[exc_types %in% colnames(counts.Callaway.locked.patchseq)]

par(mar=c(14, 4.1, 4.1, 2.1), xpd=TRUE)
barplot(counts.Callaway.locked.patchseq[, ordered_types], main="Comparing Callaway(red) and locked_patchseq data(blue) on V1 taxonomy (Core, I1, I2 and I3 cells)",
        col=c("darkblue","red"), 
        legend = rownames(counts.Callaway.locked.patchseq), 
        beside=TRUE, 
        legend.text=TRUE,
        las=2,
        args.legend = list(x = "topleft", bty="n"))





