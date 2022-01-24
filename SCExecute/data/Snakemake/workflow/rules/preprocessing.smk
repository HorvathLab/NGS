rule fasterq:
    output:
        fastq_1= "fastqs/{accession}_1.fastq", 
        fastq_2= "fastqs/{accession}_2.fastq"
    conda:
        "../envs/fastqgen.yaml"
    params:
        extra="--skip-technical"
#    log:
#        "2> log/fasterq.log"
    wrapper:
        "0.78.0/bio/sra-tools/fasterq-dump"

## Add in rule for combing SRR's into Biosample, for now must do manually

## Add in rule for trimming {sample}_1.fastq, for now must do manually
rule STARsolo:
    input:
        genome= config["Genome_dir"], 
        GTF= config["GTF_file"], 
        fastq_1_biosam= "fastqs/{sample}_1.fastq", 
        fastq_2_biosam= "fastqs/{sample}_2.fastq", 
        whitelist= config["STAR"]["Whitelist"], 
#        UMI_len= config["STAR"]["UMI_length"], 
#        overhang_len= config["STAR"]["Overhang"]
    output:
        og_bam= "fastqs/{sample}_wasp_Aligned.sortedByCoord.out.bam", 
        solo_out_barcode= "fastqs/{sample}_wasp_Solo.out/Gene/filtered/barcodes.tsv"
    threads: 2
    conda:
        "../envs/STARsolomapping.yaml"
    shell:
        "for i in {config[Biosample]}; do STAR --genomeDir {input.genome} --sjdbGTFfile {input.GTF} --sjdbOverhang 150 --limitOutSJcollapsed 2000000 --twopassMode Basic --readFilesIn fastqs/${{i}}_2.fastq fastqs/${{i}}_1.fastq --soloType CB_UMI_Simple --soloUMIlen 12 --soloCBwhitelist {input.whitelist} --outSAMattributes NH HI AS nM NM MD CB UB CR UR GX GN sS sQ sM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix fastqs/${{i}}_wasp_ --soloFeatures Gene GeneFull SJ Velocyto --runThreadN {threads} ; done"
        

rule samtoolsprelude:
    input:
        og_bam= "fastqs/{sample}_wasp_Aligned.sortedByCoord.out.bam"
    output: 
        new_bam= "fastqs/{sample}_wasp_Aligned.sortedByCoord_vW_filt.bam", 
        bam_bai= "fastqs/{sample}_wasp_Aligned.sortedByCoord_vW_filt.bam.bai"
    conda:
        "../envs/samtoolsprelude.yaml"
    shell:
        "for i in {config[Biosample]}; do cd {config[Path2fastqs]}; samtools view -H {input.og_bam} > header_${{i}}'.sam' ; samtools view {input.og_bam} | awk -F '\t' '$0 !~ /vW:i:[2-4]/{{print $0}}' > temp_aa.sam ; cat header_${{i}}'.sam' temp_aa.sam > ${{i}}'_wasp_vW_filt.sam' ; samtools sort -T temp_aa -o {output.new_bam} ${{i}}'_wasp_vW_filt.sam' ; samtools index {output.new_bam} > {output.bam_bai}; done"


rule scExecute_gatk:
    input:
        new_bam= expand("fastqs/{sample}_wasp_Aligned.sortedByCoord_vW_filt.bam", sample=config["Biosample"])
    output: 
        gatk_filt_dir= directory("fastqs/gatk_output")
    threads: 2
    conda:
        "../envs/gatk.yaml"
    shell:
        "ml scExecute; mkdir -p {config[Path2logs]}; mkdir -p {config[Path2gatk]}; mkdir -p {config[Path2samples]}; for i in {config[Biosample]}; do cd {config[Path2fastqs]}; scExecute -t {threads} -i -r ${{i}}_wasp_Aligned.sortedByCoord_vW_filt.bam -b ${{i}}_wasp_Solo.out/Gene/filtered/barcodes.tsv -C '/data/lab/Snakemake/Snakemake/workflow/scripts/GATK_scExec_111821.sh {{}}'; done; mv *_gatk_filt.vcf {config[Path2gatk]}; rm *_split_rg.bam; rm *_gatk.vcf"

barcodes= glob_wildcards("fastqs/gatk_output/{sample}_wasp_Aligned.sortedByCoord_vW_filt.{barcode}_gatk_filt.vcf").barcode 
print(barcodes)

rule scExecute_strelka:
    input:
        new_bam= expand("fastqs/{sample}_wasp_Aligned.sortedByCoord_vW_filt.bam", sample=config["Biosample"])
    output:
        strelka_filt_dir= directory("fastqs/strelka")
    threads: 2
    conda:
        "../envs/strelka.yaml"
    shell:
        "ml scExecute; mkdir -p {config[Path2strelka]}; for i in {config[Biosample]}; do cd {config[Path2fastqs]}; scExecute -t {threads} -i -r ${{i}}_wasp_Aligned.sortedByCoord_vW_filt.bam -b ${{i}}_wasp_Solo.out/Gene/filtered/barcodes.tsv -C '/data/lab/Snakemake/Snakemake/workflow/scripts/strelka_scExec_111821.sh Germline {{}}'; done; mv *_strelka_filt.vcf {output.strelka_filt_dir}; mv log_* {config[Path2logs]}; mv _${{i}}_* {config[Path2strelka]}; rm *_strelka.vcf.gz; rm *_strelka.vcf.gz.csi"
        
                
#################################

