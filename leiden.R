library(leiden)
library(igraph)
graph_object <- graph_from_adjacency_matrix(seurat.obj.ki@graphs$RNA_nn, mode = "directed")
adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)
partition <- leiden(adjacency_matrix, resolution_parameter = 1)
seurat.obj.ki$leiden.res.1 <- partition

partition <- leiden(adjacency_matrix, n_iterations = 4)

partition <- leiden(adjacency_matrix, partition_type = "SignificanceVertexPartition")
DimPlot(seurat.obj.im, reduction = "umap", label = T,group.by = "leiden.res.1")
seurat.obj.im$leiden.res.1.6 <- partition
partition <- leiden(adjacency_matrix, resolution_parameter = 1.5)
seurat.obj.im$leiden.res.1.5 <- partition
partition <- leiden(adjacency_matrix, resolution_parameter = 1.2)
seurat.obj.im$leiden.res.1.2 <- partition

partition <- leiden(adjacency_matrix, resolution_parameter = 1.3)
seurat.obj.im$leiden.res.1.3 <- partition

graph_object <- graph_from_adjacency_matrix(seurat.obj.im@graphs$RNA_snn, weighted = T, mode = "directed")
adjacency_matrix <- igraph::as_adjacency_matrix(graph_object)
partition <- leiden(adjacency_matrix)
table(partition)
seurat.obj.im$leiden.weight.1 <- partition
seurat.obj.im$leiden.iter.4 <- partition
