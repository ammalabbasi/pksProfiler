process mapReads {

    scratch true
    label 'mapBothReads'
    publishDir("${params.mapped_reads_dir}", mode: 'copy')
    conda "${params.minimap2_env}"

    input:
    tuple val(sampleID), path(r1_fastq), path(r2_fastq)

    output:
    tuple val(sampleID),
        path("${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"),
        path("${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz")

    script:
    """
	set -euo pipefail

    out1="${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"
    out2="${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz"

    # Skip condition
    if [[ -f "${params.mapped_reads_dir}/\$out1" && -f "${params.mapped_reads_dir}/\$out2" ]]; then
        echo "Skipping mapReads: Found \$out1, \$out2 in publishDir"

        for f in "\$out1" "\$out2"; do
            if [[ ! -f "\$f" ]]; then
                ln -s "${params.mapped_reads_dir}/\$f" . 2>/dev/null || cp "${params.mapped_reads_dir}/\$f" .
            fi
        done
        exit 0
    fi

    # ========= Actual execution =========
    echo "Running on R1"
    # Run minimap2 on hg38 reference
    minimap2 -2 -ax sr -t 16 ${params.hg38_db} ${r1_fastq} -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz

    # Run minimap2 on t2t_phix reference
    minimap2 -2 -ax sr -t 16 ${params.t2t_phix_db} ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R1.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz

    
    echo "Running on R2"
    # Run minimap2 on hg38 reference
    minimap2 -2 -ax sr -t 16 ${params.hg38_db} ${r2_fastq} -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz

    # Run minimap2 on t2t_phix reference
    minimap2 -2 -ax sr -t 16 ${params.t2t_phix_db} ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.fastq.gz -a | samtools fastq -@ 16 -f 4 -F 256 | gzip > ${sampleID}.R2.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz

    """
}
