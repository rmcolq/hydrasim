#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_csv = ""
params.dataset_csv = ""
params.reference_sample_size = 1
params.dataset_sample_size = 1
params.dataset_index = "public_database_accession"
params.dataset_coverage = "10k"
params.reference_dir = "store/references"
params.dataset_dir = "store/datasets"

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
    storeDir "${params.dataset_dir}/"
    input:
    tuple val(accession), val(platform)

    output:
    path("${accession}/*.fastq.gz")
    tuple val(accession), val(platform), path("${accession}/*_pass_1.fastq.gz"), path("${accession}/*_pass_2.fastq.gz"), optional: true, emit:paired
    tuple val(accession), val(platform), path("${accession}/*_pass.fastq.gz"), val(null), optional:true, emit: unpaired

    script:
    """
    prefetch ${accession}
    fastq-dump --outdir ${accession} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip */*.sra
    """
}

process downsample_dataset_accession {
    input:
    tuple val(accession), val(platform), path(raw_reads1)
    val coverage

    output:
    tuple val(accession), val(platform), path("${accession}.subsampled.fastq.gz")

    script:
    if ("${platform}" == "illumina"){
    """
    rasusa reads --bases ${coverage} ${raw_reads1} -o "${accession}.subsampled.fastq.gz"
    """
    } else {
    """
    rasusa reads --bases ${coverage} ${raw_reads1} ${raw_reads2} -o "${accession}_1.subsampled.fastq.gz" -o "${accession}_2.subsampled.fastq.gz"
    """
    }
}

process downsample_dataset_accession_paired {
    input:
    tuple val(accession), val(platform), path(raw_reads1), path(raw_reads2)
    val coverage

    output:
    tuple val(accession), val(platform), path("${accession}_1.subsampled.fastq.gz"), path("${accession}_2.subsampled.fastq.gz"), optional: true

    script:
    if ("${platform}" == "illumina"){
    """
    rasusa reads --bases ${coverage} ${raw_reads1} -o "${accession}.subsampled.fastq.gz"
    """
    } else {
    """
    rasusa reads --bases ${coverage} ${raw_reads1} ${raw_reads2} -o "${accession}_1.subsampled.fastq.gz" -o "${accession}_2.subsampled.fastq.gz"
    """
    }
}


workflow {
    main:
        //reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
        //subset_reference_accessions(reference_csv, params.reference_sample_size)
        //subset_reference_accessions.out.splitCsv(header: true).map { row -> tuple("${row.accession}","${row.category_id}") }.set{ reference_tuples }
        //download_reference_fasta(reference_tuples.unique())

        dataset_csv = file(params.dataset_csv, type: "file", checkIfExists:true)
        subset_dataset_accessions(dataset_csv, params.dataset_sample_size)
        subset_dataset_accessions.out.splitCsv(header: true).unique().map{row -> ["${row.public_database_accession}","${row.platform}","${row.human_filtered_reads_1}","${row.human_filtered_reads_2}"]}.set{ dataset_accessions }

        dataset_accessions.branch { accession, platform, reads1, reads2 ->
                    download: reads1 == ""
                        return tuple(accession, platform)
                    other: true
                        return tuple(accession, platform, reads1, reads2)
                    }.set { result }
        download_dataset_accession(result.download)
        download_dataset_accession.out.paired.concat(download_dataset_accession.out.unpaired, result.other).set{to_downsample}


        downsample_dataset_accession(to_downsample, "${params.dataset_coverage}")


}
