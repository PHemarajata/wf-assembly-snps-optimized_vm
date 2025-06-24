# Example command to run the pipeline with PopPUNK clusters:
nextflow run main.nf \
  -profile gcp_vm \
  --input . \
  --poppunk_clusters poppunk_prepared/filtered_poppunk_clusters.csv \
  --min_cluster_size 3 \
  --outdir results_clustered \
  --snp_package parsnp