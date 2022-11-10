# RNASeq
RNASeq Analysis for Double Stimulated Macrophages

Project Overview: 
Macrophages are immune system cells that respond to threats in the body. 
We conducted experiments to elucidate the genetic mechanisms behind this immune response, 
with the goal of idenitfying specific genes and pathways responsible for macrophage immune response.

Code Overview: 
1. Apply RPKM and TPM normalization to raw gene expression values
2. Add pseudocounts to normalized data
3. Create 4 dataframes of log2FC values of control and primary stimulus condition for each set of normalized data (with pseudocounts added)
4. Select for induced genes by applying a cutoff for log2FC values
5. Subset induced genes from normalized dataset with pseudocounts
6. Select for expressed genes by applying various cutoff values to the induced set of genes
7. Cluster data using kmeans clustering
8. Plot data in heatmaps to visualize clustering patterns
9. Plot line graphs of potentiated genes
