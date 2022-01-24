Steps for Snakemake:
--------------------

1. activate conda environment where the snakemake program is downloaded
    
2. add SRR to the config.yaml file under `SRR`
  
  2a. add the merged SRR/bioproject accession under `Biosample`

3. from the Snakemake dir, run `snakemake --cores [N] --use-conda --conda-frontend conda --verbose -p [local_rule]`
    #The local rules to run are in the Snakemake/workflow/Snakefile document, run them sequentially:
    
    # generate_fastq
    # generate_bam
    # scExecute


  3b. After generating fastqs, trim and merge the files
      and move trimmed fastq back into the fastq dir with the same naming convention (important for the script!)
  


  
