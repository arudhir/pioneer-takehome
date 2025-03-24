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
        homology1 = config['homology1'],
        homology2 = config['homology2'],
        error_rate = config.get("cutadapt_error_rate", "0.1"),
        min_length = config.get("cutadapt_min_length", "100")
    log:
        '{outdir}/logs/extract/{sample}_cutadapt.log'
    run:
        shell(
            """
            cutadapt \
                -g "{params.homology1}...{params.homology2}" \
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

rule find_barcodes:
    input:
        fastq = '{outdir}/inserts/{sample}.fastq'
    output:
        matches = '{outdir}/barcodes/{sample}.tsv'
    log:
        '{outdir}/logs/barcodes/{sample}.log'
    params:
        barcode = config['barcode']
    shell:
        r"""
        set -euo pipefail

        # header: read_id, strand, start, end, matched_sequence
        echo -e "read_id\tstrand\tstart\tend\tmatched_seq" > {output.matches}

        seqkit locate -i --degenerate -p "{params.barcode}" {input.fastq} | \
        awk -F'\t' 'NR > 1 {{print $1 "\t" $3 "\t" $4 "\t" $5 "\t" $7}}' >> {output.matches} 2> {log}
        """
