library(Giotto)

results_folder = '/T1_bkc/'
python_path = NULL
if(is.null(python_path)) {
  installGiottoEnvironment()
}
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = python_path)

data_path = 'BKC_Task1'
my_visium = createGiottoVisiumObject(
  visium_dir = data_path,
  expr_data = "filter",
  gene_column_index = 1,
  h5_visium_path = 'BKC_Task1/Train_filtered_feature_bc_matrix.h5',
  h5_gene_ids = "symbols",
  h5_tissue_positions_path = 'BKC_Task1/Train_tissue_positions_list.csv',
  h5_image_png_path = 'BKC_Task1/Train_tissue_hires_image.png',
  h5_json_scalefactors_path = 'BKC_Task1/Train_json.json',
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

spatPlot2D(gobject = my_visium, show_image = T, point_alpha = 0.7)
spatPlot2D(gobject = my_visium, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_feats', color_as_factor = F)
my_visium <- calculateHVF(gobject = my_visium)
my_visium <- runPCA(gobject = my_visium)
my_visium <- runUMAP(my_visium, dimensions_to_use = 1:10)
my_visium <- runtSNE(my_visium, dimensions_to_use = 1:10)
my_visium <- createNearestNetwork(gobject = my_visium, dimensions_to_use = 1:10, k = 15)
my_visium <- doLeidenCluster(gobject = my_visium, resolution = 0.4, n_iterations = 1000)
plotUMAP(gobject = my_visium, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)


combineMetadata(my_visium)
ansData = pDataDT(gobject = my_visium)
df = ansData
df <- subset(df, select = c(1,8))
colnames(df) <- c('Id','Expected')
write.csv(df, file = "./result.csv", row.names = FALSE)