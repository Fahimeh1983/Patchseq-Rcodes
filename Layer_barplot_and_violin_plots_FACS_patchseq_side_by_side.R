######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 10/3/2020            #######################################################################
######################################################################################################

# This code is to plot the layer distribution comparison of patchseq exc cells with respect to facs 
# Also to plot the gene expression comparison between the patchseq and facs

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
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Patchseq-Rcodes",
              "/Formating_patchanno_helper_functions.R"))
source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Patchseq-Rcodes/",
              "layer_jitter_barplot_dotplot.R"))


######################################################################################################
### Input path: ######################################################################################
######################################################################################################

facs.anno.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/anno.feather")
facs.data.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/",
                        "mouse_V1_ALM_20180520/data.feather")
facs.dend.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/",
                        "star/mouse_patchseq_VISp_20180626_collapsed40_cpm/dend.RData")
patchseq.path = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq",
                       "/star/mouse_patchseq_VISp_current")
select.markers.path = paste0("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                             "Taxonomies/AIT2.3.1/select.markers.rda")
work.dir = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb",
                  "/EXC_patchseq_paper_2020/")

query.dir = paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/",
                   "Mouse/patchseq/R_Object/")


locked.data = "exc_mouse_data_brl_rd_withEphysFeatures.csv"
#batch_date="20200318_BT014-RSC-252"
layer.label.data = "topview_locations.csv"

######################################################################################################
### Prepare facs and patchseq input with no filtering: ###############################################
######################################################################################################

FACS_dend <- readRDS(facs.dend.path)

plot(FACS_dend)

out <- Build_ref_dendclusterid_from_facs(facs.anno.path = facs.anno.path, 
                                         facs.dend.path = facs.dend.path,
                                         region = "VISp")
ref_dendclusterid_from_facs <- out$REF
facs.anno <- out$facs.anno
rm(out)

patchseq.anno <- Assign_clusterlabel_clusterid_for_patchseq(patchseq.path = patchseq.path, 
                                                            facs.dend.path = facs.dend.path, 
                                                            ref_cluster_label_color = ref_dendclusterid_from_facs)

dim(patchseq.anno)
Locked_spec_id <- as.character(read.csv(paste0(work.dir, locked.data), header = TRUE)[["Cell.Specimen.Id"]])
Locked_sample_id <- patchseq.anno[patchseq.anno$spec_id_label %in% Locked_spec_id, "sample_id"]
patchseq.anno <- patchseq.anno[patchseq.anno$sample_id %in% Locked_sample_id,]
dim(patchseq.anno)
#ci1 <- patchseq.anno[patchseq.anno$Tree_call_label %in% c("Core", "I1"), "spec_id_label"]

layer.labels <- read.csv(paste0(work.dir, layer.label.data)) %>% select("specimen_id", "lims_struct")
colnames(layer.labels) <- c("spec_id_label", "lims_struct")
layer.labels$lims_struct <- as.character(layer.labels$lims_struct)
patchseq.anno <- merge(patchseq.anno, layer.labels)
patchseq.anno$structure_label <- patchseq.anno$lims_struct


######################################################################################################
### WE  APPLY FILTER HERE: ###########################################################################
######################################################################################################

facs.anno <- facs.anno[facs.anno$class_label == "Glutamatergic",]
dim(facs.anno)
dim(patchseq.anno)

patchseq.anno <- patchseq.anno %>% filter(Tree_call_label %in% c("Core", "I1"))

# Here we modify the layer information for patchseq data. Meaning that if some cells are comoing from
# multiple different regions in L4 ... we call all of them L4 and so on ... 
# We are combining all the cells that bel
# Modify patchseq layer
patchseq.anno <- Modify_layer_label_GABAcells(patchseq.anno)
unique(patchseq.anno[,c("structure_label", "Revisited_layer_label")])
dim(patchseq.anno)

# Select common cres for layer comparison
table(patchseq.anno$cre_label)

common_cres <- unique(patchseq.anno$cre_label)[
  unique(patchseq.anno$cre_label) %in% unique(facs.anno$cre_label)]

facs.plot.data <- facs.anno %>% 
  filter(cre_label %in% common_cres) 

patch.plot.data <- patchseq.anno %>% 
  filter(cre_label %in% common_cres)

dim(facs.plot.data)
dim(patch.plot.data)

# remove small clusters
#tmp <- table(patch.plot.data$topLeaf_label)
#large_clusters <- names(tmp)[tmp>=5]
#patch.plot.data <- patch.plot.data[patch.plot.data$topLeaf_label %in% large_clusters,]
#dim(patch.plot.data)

patchseq.fig1.cells <- patch.plot.data[,"sample_id"]
#facs.fig1.cells <- facs.plot.data[,"sample_id"]$sample_id

#Total number of facs cells in the plot
facs.fig1.cells <- facs.plot.data[facs.plot.data$cluster_label %in% unique(patchseq.anno[patchseq.anno$sample_id %in% patchseq.fig1.cells, "cluster_label"]),]$sample_id
rownames(facs.anno) <- facs.anno$sample_id

facs.plot.data <- facs.anno[facs.fig1.cells, ]

######################################################################################################
### Plotting layers: #################################################################################
######################################################################################################

dendcluster_ids <- 
  sort(unique(patch.plot.data[,c("cluster_label", "dendcluster_id")])[,2])

#sst_pvalb_dendcluster_ids <- 
#  sort(unique(patch.plot.data[patch.plot.data$subclass_label %in% 
#                                c("Sst", "Pvalb"), "dendcluster_id"]))

#Lamp5_vip_sncg_dendcluster_ids <- 
#  sort(unique(patch.plot.data[patch.plot.data$subclass_label %in% 
#                                c("Lamp5", "Sncg", "Vip", "Serpinf1"), "dendcluster_id"]))

l6_color = unique(facs.plot.data[facs.plot.data$layer_label == "L6", "layer_color"])
l6_id = unique(facs.plot.data[facs.plot.data$layer_label == "L6", "layer_id"])


p1 <- Build_layer_barplot_perlayer_data(as.data.frame(facs.plot.data %>% 
                                                        mutate(layer_label = ifelse(layer_label == "L6b", "L6", layer_label)) %>%
                                                        mutate(layer_color = ifelse(layer_label == "L6", l6_color$layer_color, layer_color)) %>%
                                                        mutate(layer_id = ifelse(layer_label == "L6", l6_id$layer_id, layer_id))),
                                        as.data.frame(patch.plot.data %>% 
                                                        mutate(Revisited_layer_label = case_when(
                                                          Revisited_layer_label == "VISp1" ~ "L1",
                                                          Revisited_layer_label == "VISp2/3" ~ "L2/3",
                                                          Revisited_layer_label == "VISp4" ~ "L4",
                                                          Revisited_layer_label == "VISp5" ~ "L5",
                                                          Revisited_layer_label == "VISp6" ~ "L6"))),
                                        FACS_dend,
                                        dendcluster_ids)

ggsave(paste0(work.dir, "/Figures/layer_barplot_Oct13.pdf"),
       p1, 
       width = 16, height = 12, 
       useDingbats = F, dpi = 700)


######################################################################################################
### Plotting violins: ################################################################################
######################################################################################################

REF.anno.feather <- feather::read_feather(facs.anno.path) 
REF.anno.feather <- REF.anno.feather[REF.anno.feather$sample_id %in% facs.fig1.cells,]

REF.data.feather <- feather::feather(facs.data.path)
REF.data.feather <- REF.data.feather[REF.data.feather$sample_id %in% facs.fig1.cells,]
REF.data.feather <- REF.data.feather[REF.data.feather$sample_id %in% REF.anno.feather$sample_id,]

data_file <- paste0(patchseq.path, "/data.feather")
map.data.feather <- feather::feather(data_file)
map.data.feather <- map.data.feather[map.data.feather$sample_id %in% patchseq.fig1.cells,]


exc_genes <- c("Slc30a3", "Cux2", "Rorb", "Deptor", "Scnn1a", "Rspo1", "Hsd11b1",
                     "Batf3", "Oprk1", "Osr1", "Car3", "Fam84b", "Chrna6", "Pvalb", "Pappa2",
                     "Foxp2", "Slc17a8", "Trhr", "Tshz2", "Rapgef3", "trh", "Gpr139",
                     "Nxph4", "Rprm", "Crym")


p3<- group_violin_FACS_patch_plot(group_by = "dendcluster", 
                                  clusters = dendcluster_ids,
                                  genes = exc_genes,
                                  logscale = TRUE,
                                  labelheight = 2,
                                  max_width = 10,
                                  fontsize = 10,
                                  showcounts = FALSE, 
                                  REF.anno.feather = REF.anno.feather,
                                  REF.data.feather = REF.data.feather,
                                  map.anno.feather = patch.plot.data[
                                    patch.plot.data$sample_id %in% patchseq.fig1.cells,],
                                  map.data.feather = map.data.feather,
                                  dend = FACS_dend)

ggsave(paste0(work.dir, "/Figures/Violin_gene_expressions_Oct13.pdf"),
       p3, 
       width = 16, height = 12, 
       useDingbats = F, dpi = 700)
