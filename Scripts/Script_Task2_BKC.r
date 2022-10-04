library(Giotto)
library(mclust)

results_folder = '/gr/'

python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path)
data_path = 'Task2'

my_visium = createGiottoVisiumObject(
  visium_dir = data_path,
  expr_data = "filter",
  gene_column_index = 1,
  h5_visium_path = 'Task2/Train_filtered_feature_bc_matrix.h5',
  h5_gene_ids = "symbols",
  h5_tissue_positions_path = 'Task2/Train_tissue_positions_list.csv',
  h5_image_png_path = 'Task2/Train_tissue_hires_image.png',
  h5_json_scalefactors_path = 'Task2/Train_json.json',
  png_name = 'tissue_lowres_image.png',
  do_manual_adj = FALSE,
  xmax_adj = 0,
  xmin_adj = 0,
  ymax_adj = 0,
  ymin_adj = 0,
  instructions = instrs,
  cores = NA,
  verbose = TRUE
)

pDataDT(my_visium)

showGiottoImageNames(my_visium)

spatPlot(gobject = my_visium, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7)
metadata = pDataDT(my_visium)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
my_visium = subsetGiotto(my_visium, cell_ids = in_tissue_barcodes)

my_visium <- filterGiotto(gobject = my_visium,
                             expression_threshold = 1,
                             feat_det_in_min_cells = 50,
                             min_det_feats_per_cell = 1000,
                           expression_values = c('raw'),
                            verbose = T)

my_visium <- normalizeGiotto(gobject = my_visium, scalefactor = 6000, verbose = T)

my_visium <- addStatistics(gobject = my_visium)


#Dimension Reduction ---------------------------------------------------
my_visium <- calculateHVF(gobject = my_visium)
my_visium <- runPCA(gobject = my_visium)
my_visium <- runUMAP(my_visium, dimensions_to_use = 1:10)
my_visium <- runtSNE(my_visium, dimensions_to_use = 1:10)
#Dimension Reduction ---------------------------------------------------


#Selecting the best #cluster using kmeans-------------------------------
Data = pDataDT(gobject = my_visium)
df = Data
ds <- read.csv("Task2/Train_metadata.csv")
Id <- c()
Cluster <- c()
dc <- data.frame(Id, Cluster)
i = 1
for(idx in 1: nrow(ds))
{
  if(ds[idx,'Id'] == df[i,'cell_ID']){
    i <- i + 1
    dc[nrow(dc) + 1,'Id'] <- ds[idx,'Id'] 
    dc[nrow(dc),'Cluster'] <- ds[idx,'Cluster']
  }
}
best = -1
my_k = 0
for (val in 5: 10)
{
  my_visium <- doKmeans(
    gobject=my_visium,
    feat_type = NULL,
    spat_unit = NULL,
    expression_values = "normalized",
    feats_to_use = NULL,
    genes_to_use = NULL,
    dim_reduction_to_use = "pca",
    dim_reduction_name = "pca",
    dimensions_to_use = 1:10,
    distance_method = "original",
    centers = val,
    iter_max = 100,
    nstart = 1000,
    algorithm = "Hartigan-Wong",
    name = "kmeans",
    return_gobject = TRUE,
    set_seed = TRUE,
    seed_number = 1234
  )
  ##my_visium <- doKmeans(gobject = my_visium, dimensions_to_use = 1:10, k = val)
  Data = pDataDT(gobject = my_visium)
  df = Data
  Id <- c()
  kmeans <- c()
  dq <- data.frame(Id, kmeans)
  for(idx in 1: nrow(df))
  {
      dq[nrow(dq) + 1,'kmeans'] <- df[idx,'kmeans']
  }
  dx = unlist(dq['kmeans'])
  dy = unlist(dc['Cluster'])
  
  Ari = adjustedRandIndex(dx,dy)  
  if(Ari > best){
    best = Ari
    my_k = val
  }
}
#Selecting the best #cluster using kmeans-------------------------------

#Selecting the top 100 genes for performing HMRF------------------------
my_visium <- createSpatialGrid(gobject = my_visium,
                                   sdimx_stepsize = 400,
                                   sdimy_stepsize = 400,
                                   minimum_padding = 0)
my_visium = createSpatialNetwork(gobject = my_visium, minimum_k = 0)
kmtest = binSpect(my_visium)

ext_spatial_genes = kmtest[1:500]$feats
spat_cor_netw_DT = detectSpatialCorFeats(my_visium,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network',
                                         subset_feats = ext_spatial_genes)

spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = my_k)

netw_ranks = rankSpatialCorGroups(my_visium, spatCorObject = spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                  save_param = c(save_name = '22-z2-rank_correlated_groups',
                                                 base_height = 3, base_width = 5))
top_netw_spat_cluster = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus',
                                            selected_clusters = my_k, show_top_feats = 1)

# 5. create metagene enrichment score for clusters
cluster_genes_DT = showSpatialCorFeats(spat_cor_netw_DT, use_clus_name = 'spat_netw_clus', show_top_feats = 1)
cluster_genes = cluster_genes_DT$clus; names(cluster_genes) = cluster_genes_DT$feat_ID

my_visium = createMetafeats(my_visium, feat_clusters = cluster_genes, name = 'cluster_metagene')

my_visium = createSpatialNetwork(gobject = my_visium, minimum_k = 2, name = 'Delaunay_full')
my_spatial_genes <- kmtest[1:100]$feats
#Selecting the top 100 genes for performing HMRF------------------------


# HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF -----------------------------------------------------------------------------------------------------------------------

hmrf_folder = paste0(results_folder,'/','HMRF_results/')
if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

HMRF_spatial_genes = doHMRF(gobject = my_visium,
                            expression_values = 'scaled',
                            spatial_network_name = 'Delaunay_full',
                            spatial_genes = my_spatial_genes,
                            k = my_k,
                            betas = c(0, 1, 2),
                            output_folder = paste0(hmrf_folder, '/', 'Spatial_genes/SG_topgenes_k7_scaled'))


my_visium = addHMRF(gobject = my_visium,
                        HMRFoutput = HMRF_spatial_genes,
                        k = my_k, betas_to_add = c(0, 2),
                        hmrf_name = 'HMRF')
spatPlot(gobject = my_visium, cell_color = 'HMRF_k7_b.0', point_size = 5)
combineMetadata(my_visium)
# HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF HMRF -----------------------------------------------------------------------------------------------------------------------


viewer_folder = paste0(results_folder, '/', 'mouse_my_visium_viewer')
ansData = pDataDT(gobject = my_visium)
df = ansData
df <- subset(df, select = c(1,9))
colnames(df) <- c('Id','Expected')
write.csv(df, file = "./result.csv", row.names = FALSE)


