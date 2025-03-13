rule qc_reads:
    input:
        fastq = '{outdir}/raw/{sample}.fastq'
    output:
        json = '{outdir}/raw_nanoq/{sample}.json'
    log:
        '{outdir}/logs/raw_nanoq/{sample}.log'
    run:
        shell(
            'nanoq -i {input.fastq} --json --report {output.json} > {log} 2>&1'
        )

rule extract_inserts:
    input:
        fastq = '{outdir}/raw/{sample}.fastq',
        json = rules.qc_reads.output.json
    output:
        inserts = '{outdir}/inserts/{sample}.fastq'
    params:
        barcode5 = config['barcode5'],
        barcode3 = config['barcode3'],
        error_rate = config.get("cutadapt_error_rate", "0.1"),
        min_length = config.get("cutadapt_min_length", "100")
    log:
        '{outdir}/logs/extract/{sample}_cutadapt.log'
    run:
        shell(
            """
            cutadapt \
                -g "{params.barcode5}...{params.barcode3}" \
                -e {params.error_rate} \
                -m {params.min_length} \
                --info-file={output.inserts}.info \
                -o {output.inserts} \
                {input.fastq} > {log} 2>&1
            """
        )

rule qc_inserts:
    input:
        fastq = rules.extract_inserts.output.inserts
    output:
        json = '{outdir}/inserts_nanoq/{sample}.json'
    log:
        '{outdir}/logs/inserts_nanoq/{sample}.log'
    run:
        shell(
            'nanoq -i {input.fastq} --json --report {output.json} > {log} 2>&1'
        )