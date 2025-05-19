#!/usr/bin/env nextflow

process msi_sensor {
    tag "${sample_id}"
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/trimmed/", mode: params.publish_dir_mode, failOnError: true
    container "${params.alignment_image_uri}:${params.alignment_image_version}"
    cpus params.trimming_cores
    memory "${params.trimming_mem}"

    input:
    tuple val(sample_id), path(fastqs)
    path(adapter_fa)
    
    output:
    val(sample_id), emit: sample_id
    path('trimmed/*.filtered.fastq.gz'), emit: filter_reads
    
    script:
    adapter_str = params.trim_illumina_adapters ? "--adapter_fa ${adapter_fa}" : ""
    trim_str = params.trim_illumina_adapters ? "--trim_illumina_adapters" : ""
    nochop_str = params.nochop ? "--skip_porechop" : ""
    """
    filter_wrapper.py \
        --threads ${params.trimming_cores} \
        --platform ${params.platform} \
    """
}
