######################################################################################################
### Author: Fahimeh Baftizadeh #######################################################################
### Date: 2/20/2021            #######################################################################
######################################################################################################

# this code is to map our Patchseq, Tolias M1 data and Callaway on to V1(AIT2.3.1) and V1+ALM taxonomies

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


load(select.markers.path)

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###load reference data and tree
load(facs.dend.V1.ALM.labeled.path)
plot(dend.labeled)
#These are all the cell types in the V1+ALM tree
labels(dend.labeled)

#Reading FACS data which are the facs cpm values and anno file which has all the metadata
data <- read_feather(facs.data.path)
anno <- read_feather(facs.anno.path)
#checking 
table(anno$class_label)
#removing low quality and non neuronal cells as V1+ALM taxonomy does not have those types
#if you want to include those in the mapping, then we need to remake the tree 
anno <- anno[anno$class_label %in% c("GABAergic", "Glutamatergic"),]
data <- data[data$sample_id %in% anno$sample_id,]
data <- data %>% column_to_rownames("sample_id")
norm.dat <- log2(data + 1)
norm.dat <- t(norm.dat)

#In total 22439 cells remain
dim(anno)
dim(norm.dat)

cl.df <- unique(as.data.frame(anno[,c("cluster_id", "cluster_label", 
                                      "cluster_color", "subclass_id", 
                                      "subclass_label", "class_id", 
                                      "class_label", "cl", "dendcluster_id")]))

dim(cl.df)
cl.df <- cl.df[order(cl.df$dendcluster_id),]

cl <- set_names(anno[, "cluster_label"]$cluster_label, anno[, "sample_id"]$sample_id)
cl <- factor(cl)

facs.cells <- names(cl)

#loading the v1 taxonomy
dend.v1 <- readRDS(facs.dend.V1.path)
plot(dend.v1)
labels(dend.v1)
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
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]
patch.seq.norm = query.dat.norm
dim(patch.seq.norm)

######################################################################################################
### Reading M1 data            #######################################################################
######################################################################################################

introns <-read.csv(paste0(Tolias.data.path, "m1_patchseq_intron_counts.csv"))
M1.tolias.cpm <- read.csv(paste0(Tolias.data.path, Tolias.patch.m1.file))

rownames(introns) <- introns$X
introns <- introns[,setdiff(colnames(introns), "X")]
dim(introns)

rownames(M1.tolias.cpm) <- rownames(introns)
M1.tolias.cpm <- M1.tolias.cpm[,setdiff(colnames(M1.tolias.cpm), "X")]
dim(M1.tolias.cpm)

M1.tolias.logcpm <- log2(M1.tolias.cpm + 1)
rm(M1.tolias.cpm)
dim(M1.tolias.logcpm)

######################################################################################################
### Reading Callaway data            #################################################################
######################################################################################################

callaway.counts <- read.table(file=paste0(work.dir, "GSE133230_rawCount_945filtered_cells.tsv"), 
                              sep = '\t', 
                              header = TRUE)

# We need to convert the counts to logcpm to be able to map them on FACS
dim(callaway.counts)
callaway.genes <- callaway.counts$Symbol
callaway.cells <- setdiff(colnames(callaway.counts), c("GeneID", "Symbol"))
callaway.log.cpm <- logCPM(as.matrix(callaway.counts[,callaway.cells]))

callaway.log.cpm %>% as.data.frame(callaway.log.cpm)
rownames(callaway.log.cpm) <- callaway.genes
dim(callaway.log.cpm)

######################################################################################################
### Check if genes are the same in V1 and M1 data     ################################################
######################################################################################################


Total_number_M1_genes <- length(rownames(M1.tolias.logcpm))
Total_number_M1_genes
#42466

Total_number_of_V1_genes <- length(rownames(norm.dat))
Total_number_of_V1_genes
#45768

Total_number_of_callaway_genes <- length(callaway.genes)
Total_number_of_callaway_genes
#48440

#common genes for M1 and V1
M1.V1.common.genes <- rownames(M1.tolias.logcpm)[rownames(M1.tolias.logcpm) %in% rownames(norm.dat)]
#29158

#common genes for callaway and V1
Callaway.V1.common.genes <- tolower(rownames(callaway.log.cpm))[tolower(rownames(callaway.log.cpm)) %in% tolower(rownames(norm.dat))]
#31003
#Callaway gene names are all uppercase, So I am going to change them to how we call them in V1 dataset
Callaway.V1.common.genes.correctnames <- rownames(norm.dat)[match(Callaway.V1.common.genes, tolower(rownames(norm.dat)))]
sum(tolower(Callaway.V1.common.genes.correctnames) == Callaway.V1.common.genes)
#Now I change all their gene names to lower case 
rownames(callaway.log.cpm) <- tolower(rownames(callaway.log.cpm))


# subsetting facs, patchseq, Callaway and M1 data for common genes only
M1.tolias.logcpm <- M1.tolias.logcpm[M1.V1.common.genes,]
dim(M1.tolias.logcpm)
callaway.log.cpm <- callaway.log.cpm[Callaway.V1.common.genes,]
dim(callaway.log.cpm)
## Here I change the names to the gene names in our dataset
rownames(callaway.log.cpm) <- Callaway.V1.common.genes.correctnames

M1.query.dat.cells <- colnames(M1.tolias.logcpm)
patch.query.cells <- colnames(patch.seq.norm)
Callaway.cells <- colnames(callaway.log.cpm)


# We need to compute cluster med for each gene to be used in mapping later
cl.med = get_cl_means(norm.dat, cl)
dim(cl.med)
######################################################################################################
### Mapping data using V1+ALM Tree ###################################################################
######################################################################################################

set.seed(1983)

M1.v1_alm_taxonomy.memb = map_dend_membership(dend.labeled, cl, cl.med[M1.V1.common.genes,], norm.dat[M1.V1.common.genes,], M1.tolias.logcpm, M1.query.dat.cells, 
                                              bs.num=100, p=0.7, low.th=0.15)
dim(M1.v1_alm_taxonomy.memb)

Callaway.v1_alm_taxonomy.memb = map_dend_membership(dend.labeled, cl, cl.med[Callaway.V1.common.genes.correctnames, ], norm.dat[Callaway.V1.common.genes.correctnames,], 
                                                    callaway.log.cpm, callaway.cells, bs.num=100, p=0.7, low.th=0.15)
dim(Callaway.v1_alm_taxonomy.memb)

patch.v1_alm_taxonomy.memb = map_dend_membership(dend.labeled, cl, cl.med, norm.dat, patch.seq.norm, patch.query.cells, 
                                   bs.num=100, p=0.7, low.th=0.15)
dim(patch.v1_alm_taxonomy.memb)

######################################################################################################
### Mapping data using V1 Tree #######################################################################
######################################################################################################
#dend.v1 is the v1 taxonomy and dend.labeled is the v1+ALM taxonomy. First lets check if all the 
#v1 cell types exist in v1+ALM cell types and there is no change in the names 
labels(dend.v1) %in% labels(dend.labeled)

#Now we should find the cells that are in v1 taxonomy
v1.cells <- anno[anno$cluster_label %in% labels(dend.v1) & anno$region_label=="VISp","sample_id"]$sample_id
#And compute the cl.med for the v1 cells only
v1.cl <- droplevels(cl[v1.cells])
v1.cl.med = get_cl_means(norm.dat[,v1.cells], v1.cl)
dim(v1.cl.med)

M1.v1_taxonomy.memb = map_dend_membership(dend.v1, v1.cl, v1.cl.med[M1.V1.common.genes,], norm.dat[M1.V1.common.genes, v1.cells], 
                                          M1.tolias.logcpm, M1.query.dat.cells, bs.num=100, p=0.7, low.th=0.15)
dim(M1.v1_taxonomy.memb)

Callaway.v1_taxonomy.memb = map_dend_membership(dend.v1, v1.cl, v1.cl.med[Callaway.V1.common.genes.correctnames, ], norm.dat[Callaway.V1.common.genes.correctnames, v1.cells], 
                                                callaway.log.cpm, callaway.cells, bs.num=100, p=0.7, low.th=0.15)
dim(Callaway.v1_taxonomy.memb)

patch.v1_taxonomy.memb = map_dend_membership(dend.v1, v1.cl, v1.cl.med, v1.norm.dat, patch.seq.norm, patch.query.cells, 
                                                 bs.num=100, p=0.7, low.th=0.15)
dim(patch.v1_taxonomy.memb)

######################################################################################################
### KLdiv for v1_alm taxonomy ########################################################################
######################################################################################################
script_rep="//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"
source(file.path(script_rep,"patchseq/Mapping_helper_functions.R"))

# File Tree_mapping_prob is being loaded, this is a 117 by 117 matrix which shows what is the 
#probability of each cell type being confused with other cell types. This matrix was computed
#using mapping_ref_on_ref.R code, which is basically mapping FACS onto FACS data. The facs data used
#was 22439 neuron only cells
load(paste0(work.dir, "/V1_ALM_REF_PROB_MAT.rda"))
dim(Tree_mapping_probability)
v1_alm.select.cl = labels(dend.labeled)

M1.v1_alm_taxonomy.KLdiv = compute_KLdiv(select.cl = v1_alm.select.cl,
                                         select.cells = M1.query.dat.cells, 
                                         mapping_probability = Tree_mapping_probability, 
                                         memb = M1.v1_alm_taxonomy.memb)
dim(M1.v1_alm_taxonomy.KLdiv)
#1329 117


# Callaway data did not map to a lot of types from ALM
missing.cls <- v1_alm.select.cl[!v1_alm.select.cl %in% colnames(Callaway.v1_alm_taxonomy.memb)]
# in order to run compute KL_div I add these types with mapping values equal to zero to the 
# membership matrix
Callaway.v1_alm_taxonomy.memb <- as.data.frame.matrix(Callaway.v1_alm_taxonomy.memb)
#now you have to add the missing cols, here only one type was missing
Callaway.v1_alm_taxonomy.memb$"Pvalb Akr1c18 Ntf3" <- 0.0

Callaway.v1_alm_taxonomy.KLdiv = compute_KLdiv(select.cl = v1_alm.select.cl,
                                         select.cells = callaway.cells, 
                                         mapping_probability = Tree_mapping_probability, 
                                         memb = Callaway.v1_alm_taxonomy.memb)
dim(Callaway.v1_alm_taxonomy.KLdiv)

patch.v1_alm_taxonomy.KLdiv = compute_KLdiv(select.cl = v1_alm.select.cl,
                              select.cells = patch.query.cells, 
                              mapping_probability = Tree_mapping_probability, 
                              memb = patch.v1_alm_taxonomy.memb)
#11278 117
dim(patch.v1_alm_taxonomy.KLdiv)

######################################################################################################
### KLdiv for v1 taxonomy ############################################################################
######################################################################################################
path.to.ref.prob.mat <- paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/",
                               "Taxonomies/AIT2.3.1/REF_mapping_probability.rda")
load(path.to.ref.prob.mat)
dim(Tree_mapping_probability)

v1.select.cl = labels(dend.v1)
M1.v1_taxonomy.KLdiv = compute_KLdiv(select.cl = v1.select.cl,
                                     select.cells = M1.query.dat.cells, 
                                     mapping_probability = Tree_mapping_probability, 
                                     memb = M1.v1_taxonomy.memb)
#1329 93
dim(M1.v1_taxonomy.KLdiv)

# Callaway data did not map to come typs
missing.cls <- v1.select.cl[!v1.select.cl %in% colnames(Callaway.v1_taxonomy.memb)]
# in order to run compute KL_div I add these types with mapping values equal to zero to the 
# membership matrix
Callaway.v1_taxonomy.memb <- as.data.frame.matrix(Callaway.v1_taxonomy.memb)
#now you have to add the missing cols, here only one type was missing
Callaway.v1_taxonomy.memb$"Vip Lmo1 Fam159b" <- 0.0
Callaway.v1_taxonomy.memb$"Vip Lmo1 Myl1" <- 0.0
Callaway.v1_taxonomy.memb$"Vip Gpc3 Slc18a3" <- 0.0
Callaway.v1_taxonomy.memb$"Pvalb Akr1c18 Ntf3" <- 0.0 

Callaway.v1_taxonomy.KLdiv = compute_KLdiv(select.cl = v1.select.cl,
                                     select.cells = Callaway.cells, 
                                     mapping_probability = Tree_mapping_probability, 
                                     memb = Callaway.v1_taxonomy.memb)
#945 93
dim(Callaway.v1_taxonomy.KLdiv)


patch.v1_taxonomy.KLdiv = compute_KLdiv(select.cl = v1.select.cl,
                                        select.cells = patch.query.cells, 
                                        mapping_probability = Tree_mapping_probability, 
                                        memb = patch.v1_taxonomy.memb)

#11278 93
dim(patch.v1_taxonomy.KLdiv)

######################################################################################################
### Correlation for v1_alm taxonomy ##################################################################
######################################################################################################
M1.select.markers <- select.markers[select.markers %in% rownames(M1.tolias.logcpm)]
Callaway.select.markers <- select.markers[select.markers %in% rownames(callaway.log.cpm)]

M1.facs.cor.v1_alm_taxonomy <- Compute_correlation_mat(markers = M1.select.markers, cells = M1.query.dat.cells,
                                                       query.dat.norm = M1.tolias.logcpm, train.cl.med = cl.med)

Callaway.facs.cor.v1_alm_taxonomy <- Compute_correlation_mat(markers = Callaway.select.markers, cells = callaway.cells,
                                                       query.dat.norm = callaway.log.cpm, train.cl.med = cl.med)

patch.facs.cor.v1_alm_taxonomy <- Compute_correlation_mat(markers = select.markers, cells = patch.query.cells,
                                                          query.dat.norm = patch.seq.norm, train.cl.med = cl.med)

######################################################################################################
### Correlation for v1 taxonomy ######################################################################
######################################################################################################

M1.facs.cor.v1_taxonomy <- Compute_correlation_mat(markers = M1.select.markers, cells = M1.query.dat.cells,
                                                   query.dat.norm = M1.tolias.logcpm, train.cl.med = v1.cl.med)

Callaway.facs.cor.v1_taxonomy <- Compute_correlation_mat(markers = Callaway.select.markers, cells = callaway.cells,
                                                             query.dat.norm = callaway.log.cpm, train.cl.med = v1.cl.med)

patch.facs.cor.v1_taxonomy <- Compute_correlation_mat(markers = select.markers, cells = patch.query.cells,
                                                      query.dat.norm = patch.seq.norm, train.cl.med = v1.cl.med)

######################################################################################################
### Aggreagte results ################################################################################
######################################################################################################

M1.v1_alm_taxonomy.memb <- M1.v1_alm_taxonomy.memb[M1.query.dat.cells, v1_alm.select.cl]
dim(M1.v1_alm_taxonomy.memb)
M1.v1_taxonomy.memb <- M1.v1_taxonomy.memb[M1.query.dat.cells, v1.select.cl]
dim(M1.v1_taxonomy.memb)

Callaway.v1_alm_taxonomy.memb <- Callaway.v1_alm_taxonomy.memb[callaway.cells, v1_alm.select.cl]
dim(Callaway.v1_alm_taxonomy.memb)
Callaway.v1_taxonomy.memb <- Callaway.v1_taxonomy.memb[callaway.cells, v1.select.cl]
dim(Callaway.v1_taxonomy.memb)

patch.v1_alm_taxonomy.memb <- patch.v1_alm_taxonomy.memb[patch.query.cells, v1_alm.select.cl]
dim(patch.v1_alm_taxonomy.memb)
patch.v1_taxonomy.memb <- patch.v1_taxonomy.memb[patch.query.cells, v1.select.cl]
dim(patch.v1_taxonomy.memb)

v1_alm_3_cls <- Get_3_best_cl(M1.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_cls) <- c("v1_alm_first_cl", "v1_alm_second_cl", "v1_alm_third_cl")

v1_3_cls <- Get_3_best_cl(M1.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_cls) <- c("v1_first_cl", "v1_second_cl", "v1_third_cl")

v1_alm_3_bts <- Get_3_best_bt(M1.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_bts) <- c("v1_alm_first_bt", "v1_alm_second_bt", "v1_alm_third_bt")

v1_3_bts <- Get_3_best_bt(M1.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_bts) <- c("v1_first_bt", "v1_second_bt", "v1_third_bt")

v1_alm_3_KL <- Get_3_best_KL(memb = M1.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, KLdiv = M1.v1_alm_taxonomy.KLdiv)
colnames(v1_alm_3_KL) <-c("v1_alm_first_KL", "v1_alm_second_KL", "v1_alm_third_KL")

v1_3_KL <- Get_3_best_KL(memb = M1.v1_taxonomy.memb, ref.cl = v1.select.cl, KLdiv = M1.v1_taxonomy.KLdiv)
colnames(v1_3_KL) <-c("v1_first_KL", "v1_second_KL", "v1_third_KL")

v1_alm_3_cor <- Get_3_best_cor(memb = M1.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, cor = M1.facs.cor.v1_alm_taxonomy)
colnames(v1_alm_3_cor) <- c("v1_alm_first_cor", "v1_alm_second_cor", "v1_alm_third_cor")

v1_3_cor <- Get_3_best_cor(memb = M1.v1_taxonomy.memb, ref.cl = v1.select.cl, cor = M1.facs.cor.v1_taxonomy)
colnames(v1_3_cor) <- c("v1_first_cor", "v1_second_cor", "v1_third_cor")

cells <- M1.query.dat.cells
results <- cbind.data.frame(v1_alm_3_cls[cells,],
                            v1_3_cls[cells,],
                            v1_alm_3_bts[cells,],
                            v1_3_bts[cells,],
                            v1_alm_3_KL[cells,],
                            v1_3_KL[cells,],
                            v1_alm_3_cor[cells,],
                            v1_3_cor[cells,])


Original_cols <- colnames(results)

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(v1_alm_not_finall_call = ifelse(v1_alm_first_cor > 0.5  & v1_alm_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_alm_call = case_when(v1_alm_not_finall_call == "Good" & v1_alm_first_bt >= 0.9 ~ "Core",
                               v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                 v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                 v1_alm_first_bt / v1_alm_second_bt >= 2 ~ "I1", 
                               v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                 v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                 v1_alm_first_bt / v1_alm_second_bt < 2 ~ "I2",
                               v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                 v1_alm_first_bt + v1_alm_second_bt < 0.7 ~ "I3",
                               v1_alm_not_finall_call == "PoorQ" ~ "PoorQ",
                               TRUE ~ "Other")) %>%
  mutate(v1_not_finall_call = ifelse(v1_first_cor > 0.5  & v1_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_call = case_when(v1_not_finall_call == "Good" & v1_first_bt >= 0.9 ~ "Core",
                                 v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                                   v1_first_bt + v1_second_bt >= 0.7 &
                                   v1_first_bt / v1_second_bt >= 2 ~ "I1", 
                                 v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                                   v1_first_bt + v1_second_bt >= 0.7 &
                                   v1_first_bt / v1_second_bt < 2 ~ "I2",
                                 v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                                   v1_first_bt + v1_second_bt < 0.7 ~ "I3",
                                 v1_not_finall_call == "PoorQ" ~ "PoorQ",
                                 TRUE ~ "Other")) %>%
  
  column_to_rownames("id") 

results <- results[,c(Original_cols, "v1_alm_call", "v1_call")]

M1.results <- results

write.csv(M1.results, paste0(work.dir, "/M1_results_on_v1_and_v1_alm_taxonomies.csv"))


v1_alm_3_cls <- Get_3_best_cl(Callaway.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_cls) <- c("v1_alm_first_cl", "v1_alm_second_cl", "v1_alm_third_cl")

v1_3_cls <- Get_3_best_cl(Callaway.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_cls) <- c("v1_first_cl", "v1_second_cl", "v1_third_cl")

v1_alm_3_bts <- Get_3_best_bt(Callaway.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_bts) <- c("v1_alm_first_bt", "v1_alm_second_bt", "v1_alm_third_bt")

v1_3_bts <- Get_3_best_bt(Callaway.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_bts) <- c("v1_first_bt", "v1_second_bt", "v1_third_bt")

v1_alm_3_KL <- Get_3_best_KL(memb = Callaway.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, KLdiv = Callaway.v1_alm_taxonomy.KLdiv)
colnames(v1_alm_3_KL) <-c("v1_alm_first_KL", "v1_alm_second_KL", "v1_alm_third_KL")

v1_3_KL <- Get_3_best_KL(memb = Callaway.v1_taxonomy.memb, ref.cl = v1.select.cl, KLdiv = Callaway.v1_taxonomy.KLdiv)
colnames(v1_3_KL) <-c("v1_first_KL", "v1_second_KL", "v1_third_KL")

v1_alm_3_cor <- Get_3_best_cor(memb = Callaway.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, cor = Callaway.facs.cor.v1_alm_taxonomy)
colnames(v1_alm_3_cor) <- c("v1_alm_first_cor", "v1_alm_second_cor", "v1_alm_third_cor")

v1_3_cor <- Get_3_best_cor(memb = Callaway.v1_taxonomy.memb, ref.cl = v1.select.cl, cor = Callaway.facs.cor.v1_taxonomy)
colnames(v1_3_cor) <- c("v1_first_cor", "v1_second_cor", "v1_third_cor")

cells <- callaway.cells
results <- cbind.data.frame(v1_alm_3_cls[cells,],
                            v1_3_cls[cells,],
                            v1_alm_3_bts[cells,],
                            v1_3_bts[cells,],
                            v1_alm_3_KL[cells,],
                            v1_3_KL[cells,],
                            v1_alm_3_cor[cells,],
                            v1_3_cor[cells,])


Original_cols <- colnames(results)

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(v1_alm_not_finall_call = ifelse(v1_alm_first_cor > 0.5  & v1_alm_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_alm_call = case_when(v1_alm_not_finall_call == "Good" & v1_alm_first_bt >= 0.9 ~ "Core",
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                   v1_alm_first_bt / v1_alm_second_bt >= 2 ~ "I1", 
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                   v1_alm_first_bt / v1_alm_second_bt < 2 ~ "I2",
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt < 0.7 ~ "I3",
                                 v1_alm_not_finall_call == "PoorQ" ~ "PoorQ",
                                 TRUE ~ "Other")) %>%
  mutate(v1_not_finall_call = ifelse(v1_first_cor > 0.5  & v1_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_call = case_when(v1_not_finall_call == "Good" & v1_first_bt >= 0.9 ~ "Core",
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt >= 0.7 &
                               v1_first_bt / v1_second_bt >= 2 ~ "I1", 
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt >= 0.7 &
                               v1_first_bt / v1_second_bt < 2 ~ "I2",
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt < 0.7 ~ "I3",
                             v1_not_finall_call == "PoorQ" ~ "PoorQ",
                             TRUE ~ "Other")) %>%
  
  column_to_rownames("id") 

results <- results[,c(Original_cols, "v1_alm_call", "v1_call")]

Callaway.results <- results

write.csv(Callaway.results,  paste0(work.dir,"/Callaway_results_on_v1_and_v1_alm_taxonomies.csv"))


v1_alm_3_cls <- Get_3_best_cl(patch.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_cls) <- c("v1_alm_first_cl", "v1_alm_second_cl", "v1_alm_third_cl")

v1_3_cls <- Get_3_best_cl(patch.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_cls) <- c("v1_first_cl", "v1_second_cl", "v1_third_cl")

v1_alm_3_bts <- Get_3_best_bt(patch.v1_alm_taxonomy.memb, v1_alm.select.cl)
colnames(v1_alm_3_bts) <- c("v1_alm_first_bt", "v1_alm_second_bt", "v1_alm_third_bt")

v1_3_bts <- Get_3_best_bt(patch.v1_taxonomy.memb, v1.select.cl)
colnames(v1_3_bts) <- c("v1_first_bt", "v1_second_bt", "v1_third_bt")

v1_alm_3_KL <- Get_3_best_KL(memb = patch.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, KLdiv = patch.v1_alm_taxonomy.KLdiv)
colnames(v1_alm_3_KL) <-c("v1_alm_first_KL", "v1_alm_second_KL", "v1_alm_third_KL")

v1_3_KL <- Get_3_best_KL(memb = patch.v1_taxonomy.memb, ref.cl = v1.select.cl, KLdiv = patch.v1_taxonomy.KLdiv)
colnames(v1_3_KL) <-c("v1_first_KL", "v1_second_KL", "v1_third_KL")

v1_alm_3_cor <- Get_3_best_cor(memb = patch.v1_alm_taxonomy.memb, ref.cl = v1_alm.select.cl, cor = patch.facs.cor.v1_alm_taxonomy)
colnames(v1_alm_3_cor) <- c("v1_alm_first_cor", "v1_alm_second_cor", "v1_alm_third_cor")

v1_3_cor <- Get_3_best_cor(memb = patch.v1_taxonomy.memb, ref.cl = v1.select.cl, cor = patch.facs.cor.v1_taxonomy)
colnames(v1_3_cor) <- c("v1_first_cor", "v1_second_cor", "v1_third_cor")

cells <- patch.query.cells
results <- cbind.data.frame(v1_alm_3_cls[cells,],
                            v1_3_cls[cells,],
                            v1_alm_3_bts[cells,],
                            v1_3_bts[cells,],
                            v1_alm_3_KL[cells,],
                            v1_3_KL[cells,],
                            v1_alm_3_cor[cells,],
                            v1_3_cor[cells,])


Original_cols <- colnames(results)

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(v1_alm_not_finall_call = ifelse(v1_alm_first_cor > 0.5  & v1_alm_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_alm_call = case_when(v1_alm_not_finall_call == "Good" & v1_alm_first_bt >= 0.9 ~ "Core",
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                   v1_alm_first_bt / v1_alm_second_bt >= 2 ~ "I1", 
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt >= 0.7 &
                                   v1_alm_first_bt / v1_alm_second_bt < 2 ~ "I2",
                                 v1_alm_not_finall_call == "Good" & v1_alm_first_bt < 0.9 &
                                   v1_alm_first_bt + v1_alm_second_bt < 0.7 ~ "I3",
                                 v1_alm_not_finall_call == "PoorQ" ~ "PoorQ",
                                 TRUE ~ "Other")) %>%
  mutate(v1_not_finall_call = ifelse(v1_first_cor > 0.5  & v1_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(v1_call = case_when(v1_not_finall_call == "Good" & v1_first_bt >= 0.9 ~ "Core",
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt >= 0.7 &
                               v1_first_bt / v1_second_bt >= 2 ~ "I1", 
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt >= 0.7 &
                               v1_first_bt / v1_second_bt < 2 ~ "I2",
                             v1_not_finall_call == "Good" & v1_first_bt < 0.9 &
                               v1_first_bt + v1_second_bt < 0.7 ~ "I3",
                             v1_not_finall_call == "PoorQ" ~ "PoorQ",
                             TRUE ~ "Other")) %>%
  
  column_to_rownames("id") 

results <- results[,c(Original_cols, "v1_alm_call", "v1_call")]

patch.results <- results
dim(patch.results)

write.csv(patch.results,  paste0(work.dir,"/patch_results_on_v1_and_v1_alm_taxonomies.csv"))

#Leaf_node_cells <- rownames(M1.mapping.df)[M1.mapping.df$resolution.index==1]
#Internal_node_cells <- rownames(M1.mapping.df)[(M1.mapping.df$resolution.index < 1 & M1.mapping.df$resolution.index > 0.7)]
#PoorQ_cells <- setdiff(rownames(M1.mapping.df), c(Leaf_node_cells, Internal_node_cells))

#results[Leaf_node_cells, "Old_call"] <- c("Leaf_node")
#results[Internal_node_cells, "Old_call"] <- c("Internal_node")
#results[PoorQ_cells, "Old_call"] <- c("PoorQ")

#ggplot(results, aes(Tree_first_bt , fill = Tree_call)) + 
#  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
#  xlab("first cluster call confidence") + ylab("Density")

#ggplot(results, aes(Tree_second_bt  , fill = Tree_call)) + 
#  geom_density(alpha = 0.3) + xlim(c(0,1)) +
#  xlab("second cluster call confidence") + ylab("Density")

#ggplot(results, aes(Tree_first_KL , fill = Tree_call)) + 
#  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
#  xlab("first cluster call confidence") + ylab("Density")

#ggplot(results, aes(Tree_second_KL  , fill = Tree_call)) + 
#  geom_density(alpha = 0.3) + xlim(c(0,1)) +
#  xlab("second cluster call confidence") + ylab("Density")


######################################################################################################
### Plot              ################################################################################
######################################################################################################
facs.anno <- anno
ref.dendcluster.ids <- as.data.frame(unique(facs.anno[,c("dendcluster_id", "cluster_label", "cluster_color", "cluster_id")]))
rownames(ref.dendcluster.ids) <- ref.dendcluster.ids$cluster_label
ref.dendcluster.ids <- ref.dendcluster.ids[labels(dend.labeled), ]
rownames(ref.dendcluster.ids) <- 1:dim(ref.dendcluster.ids)[1]

source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/",
              "Manuscript_patchseq_2019/Rcodes/layer_jitter_barplot_dotplot.R"))


M1.highlyconsistent.anno <-  patch.results %>% rownames_to_column()
M1.highlyconsistent.anno["cluster_label"] <- M1.highlyconsistent.anno$Tree_first_cl
M1.highlyconsistent.anno <- merge(M1.highlyconsistent.anno, ref.dendcluster.ids, on="cluster_label")

dim(M1.highlyconsistent.anno)
#Fake columns
M1.highlyconsistent.anno["Revisited_layer_label"] <- "L6"
M1.highlyconsistent.anno["Revisited_layer_id"] <- 10
M1.highlyconsistent.anno["Revisited_layer_color"] <- "#FFFFFF"
l6_color = unique(facs.anno[facs.anno$layer_label == "L6", "layer_color"])
l6_id = unique(facs.anno[facs.anno$layer_label == "L6", "layer_id"])

source(paste0("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/",
              "Manuscript_patchseq_2019/Rcodes/layer_jitter_barplot_dotplot.R"))

p1 <- Build_layer_barplot_perlayer_data(as.data.frame(facs.anno %>% 
                                                        mutate(layer_label = ifelse(layer_label == "L6b", "L6", layer_label))),
                                        as.data.frame(M1.highlyconsistent.anno),
                                        dend.labeled,
                                        ref.dendcluster.ids$dendcluster_id)
p1
######################################################################################################
### plot per cell mapping prob: ####################################################################
######################################################################################################

M1.Tree.memb <- rownames_to_column(as.data.frame.matrix(M1.Tree.memb))
df <- melt(M1.Tree.memb)
df <- df[df$value > 0,] # Keep only the values more than zero
colnames(df) <- c("sample_id", "cluster_label", "mapping_prob")

# The final tree call for each cell is in topleaf label column so we need to add that to our data frame
results <- results %>% rownames_to_column(var="sample_id")
Topleaf_label <- as.data.frame(results[,c("sample_id", "Tree_first_cl")])

# Also we need the dendcluster id for plotting purposes, that should be taken from FACS data
cl.df <- as.data.frame(unique(facs.anno[,c("dendcluster_id", "cluster_label")]))
colnames(cl.df) <- c("dendcluster_id", "Tree_first_cl")
cl.df

# Now we add topleaf label and dendcluster id to the data frame
df <- left_join(df, Topleaf_label)
df <- left_join(df, cl.df, by= "Tree_first_cl")
df <- df[order(df$dendcluster_id),] #order it by the dendcluster_id

######################################################################################################
### core, i1, i2, i3 and poorQ cells:#################################################################
######################################################################################################

Core <- results[results$Tree_call %in% c("Core"), "sample_id"]
I1 <- results[results$Tree_call %in% c("I1"), "sample_id"]
I2 <- results[results$Tree_call %in% c("I2"), "sample_id"]
I3 <- results[results$Tree_call %in% c("I3"), "sample_id"]
PoorQ <- results[results$Tree_call %in% c("PoorQ"), "sample_id"]

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
  
  for (t in labels(dend)) {
    tmpdf = which.df[which.df$Tree_first_cl == t, ]
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

#save_plot(paste0(work.dir, "/Figures/memb_heatmap_June25.pdf"), p, base_width = 20, base_height = 30, dpi = 700)


