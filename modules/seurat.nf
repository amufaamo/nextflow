process SEURAT_PROCESS {
    tag "$sample_id"
    container 'satijalab/seurat:4.3.0'
    publishDir "${params.outdir}/seurat/${sample_id}", mode: 'copy'

    input:
    tuple val(sample_id), path(pipseeker_dir)

    output:
    path "${sample_id}.rds", emit: rds
    path "*.pdf", emit: plots

    script:
    def sensitivity = params.sc_sensitivity
    
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(dplyr)
    library(patchwork)
    library(ggplot2)

    data_dir <- paste0("${pipseeker_dir}/filtered_matrix/sensitivity_", "${sensitivity}")
    
    if (!dir.exists(data_dir)) {
        stop("Error: Sensitivity directory not found.")
    }

    counts_data <- Read10X(data.dir = data_dir)
    seu_obj <- CreateSeuratObject(counts = counts_data, project = "${sample_id}", min.cells = 3, min.features = 200)

    seu_obj <- seu_obj %>%
        NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
        FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
        ScaleData() %>%
        RunPCA() %>%
        FindNeighbors(dims = 1:20) %>%
        FindClusters(resolution = 0.5) %>%
        RunUMAP(dims = 1:20)

    pdf("${sample_id}_umap.pdf")
    print(DimPlot(seu_obj, reduction = "umap", label = TRUE) + NoLegend())
    dev.off()

    saveRDS(seu_obj, file = "${sample_id}.rds")
    """
}
