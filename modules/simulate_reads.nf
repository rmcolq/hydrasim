#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_csv = ""
params.dataset_csv = ""
params.reference_sample_size = 1
params.dataset_sample_size = 1
params.dataset_index = "public_database_accession"
params.dataset_coverage = "10k"
params.reference_dir = "store/references"

params.simulator = "badread" // Default simulator
params.input_fasta = "genome.fasta" // Input reference genome
params.output_dir = "./output" // Output directory
params.reads_count = 100000 // Number of reads to simulate

process subset_reference_accessions {
    input:
    path reference_csv
    val sample_size

    output:
    path "references.csv"

    script:
    """
    python $launchDir/bin/subset_accessions.py ${reference_csv} "references.csv" category_id ${sample_size}
    """
}

process subset_dataset_accessions {
    input:
    path dataset_csv
    val sample_size

    output:
    path "dataset.csv"

    script:
    """
    python $launchDir/bin/subset_accessions.py ${dataset_csv} "dataset.csv" platform ${sample_size}
    """
}

process download_reference_fasta {
    maxForks 1
    storeDir "${params.reference_dir}/${category}"

    input:
    tuple val(accession), val(category)

    output:
    path "${accession}_genomic.fna"

    script:
    """
    datasets download genome accession ${accession}
    sleep 10
    unzip -o ncbi_dataset.zip
    sleep 5
    mv ncbi_dataset/data/*/*_genomic.fna ${accession}_genomic.fna
    """
}


process download_dataset_accession {
    input:
    tuple val(accession), val(platform)
    val coverage

    output:
    path "*.fq.gz"

    script:
    paired_flag = ""
    if ("${platform}" == "illumina"){
        paired_flag = " --paired"
    }
    """
    python $launchDir/bin/download_sra_accession.py --accession ${accession} --coverage ${coverage} ${paired_flag}
    """
}

process simulate_reads {
    tag "${params.simulator}"
    
    input:
    path fasta_file from params.input_fasta
    
    output:
    path "${params.output_dir}/simulated_reads.fasta" into simulated_reads
    
    script:
    if (params.simulator == "badread") {
        """
        badread simulate --reference ${fasta_file} --quantity ${params.reads_count}x --output ${params.output_dir}/simulated_reads.fasta
        """
    } else if (params.simulator == "nanosim") {
        """
        NanoSim-2.6.0/src/simulator.py fasta -i ${fasta_file} -o ${params.output_dir}/simulated_reads -n ${params.reads_count}
        """
    } else if (params.simulator == "wgsim") {
        """
        wgsim -N ${params.reads_count} -1 150 -2 150 -r 0.001 -R 0.15 ${fasta_file} ${params.output_dir}/simulated_reads_1.fq ${params.output_dir}/simulated_reads_2.fq
        """
    } else {
        error "Unsupported simulator: ${params.simulator}"
    }
}

workflow {
    main:
        reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
        subset_reference_accessions(reference_csv, params.reference_sample_size)
        subset_reference_accessions.out.splitCsv(header: true).map { row -> tuple("${row.accession}","${row.category_id}") }.set{ reference_tuples }
        download_reference_fasta(reference_tuples.unique())

        dataset_csv = file(params.dataset_csv, type: "file", checkIfExists:true)
        subset_dataset_accessions(dataset_csv, params.dataset_sample_size)
        subset_dataset_accessions.out.splitCsv(header: true).map { row -> tuple("${row.public_database_accession}", "${row.platform}") }.set{ dataset_tuples }
        dataset_tuples.view()
        download_dataset_accession(dataset_tuples, params.dataset_coverage)
}
