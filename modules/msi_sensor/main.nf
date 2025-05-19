#!/usr/bin/env nextflow

process msi_sensor {
    tag "${sample_id}"
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.results_path}/msi_sensor/", mode: params.publish_dir_mode, failOnError: true
    container "${params.cite_seq_image_uri}:${params.cite_seq_image_version}"
    cpus params.msi_sensor_cores
    memory "${params.msi_sensor_mem}"

    input:
    tuple val(sample_id), path(data)
    
    output:
    val(sample_id), emit: sample_id
    path('msi_sensor/*.msisensor_input.csv'), emit: msi_input
    path('msi_sensor/*.msisensor_output.csv'), emit: outfile
    
    script:
    """
    msi_sensor.py \
        --data ${params.trimming_cores} \
        --outdir ${params.results_path} \
        --sample_id ${params.sample_id} \
        --sample_col ${params.sample_col} \
        --min_counts ${params.min_counts} \
        --min_cells ${params.min_cells} \
        --target_sum ${params.target_sum} \
        --impute_method ${params.impute_method} \
        --model ${params.model}
    """
}
