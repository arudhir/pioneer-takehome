rule fq_to_fa:
    input:
        fastq = '{outdir}/inserts/{sample}.fastq'
    output:
        fasta = '{outdir}/inserts/{sample}.fasta'
    run:
        shell(
            'seqtk seq -a {input.fastq} > {output.fasta}'
        )

rule get_genomic_inserts:
    input:
        fasta = rules.fq_to_fa.output.fasta,
        ref = config['ref_genome']
    output:
        bam = '{outdir}/bam/{sample}.bam',
        bai = '{outdir}/bam/{sample}.bam.bai',
        temp_sam = temp('{outdir}/bam/{sample}.temp.sam')  # Temporary file that will be auto-deleted
    params:
        preset = config.get("minimap2_preset", "asm20")
    log:
        '{outdir}/logs/align/{sample}.log'
    threads: config.get("threads", {}).get("align", 8)
    run:
        shell(
            """
            # Run minimap2 and output SAM
            minimap2 -ax {params.preset} -t {threads} {input.ref} {input.fasta} > {output.temp_sam} 2> {log}
            
            # Convert to sorted BAM
            samtools view -bS {output.temp_sam} | samtools sort -@ {threads} -o {output.bam}
            
            # Index the BAM
            samtools index {output.bam}
            """
        )

rule feature_overlap:
    input:
        bam = rules.get_genomic_inserts.output.bam,
        bai = rules.get_genomic_inserts.output.bai,
        gff = config['gff']
    output:
        bed = '{outdir}/overlaps/{sample}.overlaps.bed'
    log:
        '{outdir}/logs/overlaps/{sample}.log'
    run:
        shell(
            'bedtools intersect -b {input.gff} -a {input.bam} -bed -wo > {output.bed} 2> {log}'
        )