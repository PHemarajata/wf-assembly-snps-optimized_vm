process PARSE_POPPUNK_CLUSTERS {
    tag "PopPUNK cluster parsing"
    label "process_single"

    conda "conda-forge::python=3.9 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.3.3' :
        'quay.io/biocontainers/pandas:1.3.3' }"

    input:
    path poppunk_assignments
    path input_files

    output:
    path "cluster_*.txt", emit: cluster_files
    path "cluster_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3
    
    import pandas as pd
    import os
    import glob
    from pathlib import Path
    
    # Read PopPUNK assignments
    df = pd.read_csv('${poppunk_assignments}')
    
    # Get list of input files from the current directory
    input_files = []
    # Look for all files in current directory that were staged by Nextflow
    for file_path in glob.glob('*'):
        if os.path.isfile(file_path) and file_path != '${poppunk_assignments}':
            # Check if it's a FASTA file
            for ext in ['.fasta', '.fas', '.fna', '.fsa', '.fa', '.fasta.gz', '.fas.gz', '.fna.gz', '.fsa.gz', '.fa.gz']:
                if file_path.endswith(ext):
                    input_files.append(file_path)
                    break
    
    # Create mapping of sample names to files
    sample_to_file = {}
    for file_path in input_files:
        file = os.path.basename(file_path)
        # Extract sample name from filename (remove extensions)
        sample_name = Path(file).stem
        if sample_name.endswith('.fasta'):
            sample_name = sample_name[:-6]
        elif sample_name.endswith('.fas'):
            sample_name = sample_name[:-4]
        elif sample_name.endswith('.fna'):
            sample_name = sample_name[:-4]
        elif sample_name.endswith('.fsa'):
            sample_name = sample_name[:-4]
        elif sample_name.endswith('.fa'):
            sample_name = sample_name[:-3]
        
        # Try different matching strategies
        # 1. Exact match
        sample_to_file[sample_name] = file_path
        # 2. Remove common suffixes
        clean_name = sample_name.replace('-SPAdes', '').replace('_contigs', '')
        sample_to_file[clean_name] = file_path
        # 3. Extract base sample ID (everything before first underscore or dash)
        base_name = sample_name.split('_')[0].split('-')[0]
        if len(base_name) > 3:  # Only if meaningful length
            sample_to_file[base_name] = file_path
    
    # Group by cluster
    clusters = df.groupby('Cluster')
    
    cluster_summary = []
    
    for cluster_id, group in clusters:
        cluster_files = []
        cluster_filename = f"cluster_{cluster_id}.txt"
        
        with open(cluster_filename, 'w') as f:
            for _, row in group.iterrows():
                sample_name = row['Taxon']  # Adjust column name if different
                if sample_name in sample_to_file:
                    file_path = sample_to_file[sample_name]
                    f.write(f"{file_path}\\n")
                    cluster_files.append(file_path)
        
        cluster_summary.append({
            'Cluster': cluster_id,
            'Sample_Count': len(cluster_files),
            'Samples': ';'.join([Path(f).stem for f in cluster_files])
        })
    
    # Write summary
    summary_df = pd.DataFrame(cluster_summary)
    summary_df.to_csv('cluster_summary.tsv', sep='\\t', index=False)
    
    print(f"Created {len(clusters)} cluster files")
    for _, row in summary_df.iterrows():
        print(f"Cluster {row['Cluster']}: {row['Sample_Count']} samples")
    
    # Write versions
    with open('versions.yml', 'w') as f:
        f.write('"${task.process}":\\n')
        f.write('    python: "3.9"\\n')
        f.write('    pandas: "1.3.3"\\n')
    """
}