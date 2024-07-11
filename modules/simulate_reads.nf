#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_csv = ""
params.dataset_csv = ""
params.reference_sample_size = 1
params.dataset_sample_size = 3
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
    tuple val(accession), val(platform), path("${accession}/*_pass.fastq.gz")

    script:
    """
    prefetch ${accession}
    fastq-dump --outdir ${accession} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip */*.sra
    """
}

process download_dataset_accession_paired {
    storeDir "${params.dataset_dir}/"
    input:
    tuple val(accession), val(platform)

    output:
    tuple val(accession), val(platform), path("${accession}/*_pass_1.fastq.gz"), path("${accession}/*_pass_2.fastq.gz")

    script:
    """
    prefetch ${accession}
    fastq-dump --outdir ${accession} --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip */*.sra
    """
}

process downsample_dataset_accession {
    input:
    tuple val(accession), val(platform), val(index), path(raw_reads1)
    val coverage

    output:
    tuple val(accession), val(platform), val(index), path("${accession}.subsampled.fastq.gz")

    script:
    """
    rasusa reads --bases ${coverage} ${raw_reads1} -o "${accession}.subsampled.fastq.gz"
    """
}

process downsample_dataset_accession_paired {
    input:
    tuple val(accession), val(platform), val(index), path(raw_reads1), path(raw_reads2)
    val coverage

    output:
    tuple val(accession), val(platform), val(index), path("${accession}_1.subsampled.fastq.gz"), path("${accession}_2.subsampled.fastq.gz")

    script:
    """
    rasusa reads --bases ${coverage} ${raw_reads1} ${raw_reads2} -o "${accession}_1.subsampled.fastq.gz" -o "${accession}_2.subsampled.fastq.gz"
    """
}

workflow download_and_downsample {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(accession, platform, index, reads1)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession.out, by: 0).map{ accession, platform, index, platform1, reads1 -> [accession, platform, index, reads1]}.set{downloaded}
        result.other.concat(downloaded).set{ to_downsample }
        downsample_dataset_accession(to_downsample, "${params.dataset_coverage}")
    emit:
        downsample_dataset_accession.out
}

workflow download_and_downsample_paired {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1, reads2 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(accession, platform, index, reads1, reads2)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession_paired(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession_paired.out, by: 0).map{ accession, platform, index, platform1, reads1, reads2 -> [accession, platform, index, reads1, reads2]}.set{downloaded}
        result.other.concat(downloaded).set{ to_downsample }
        downsample_dataset_accession_paired(to_downsample, "${params.dataset_coverage}")
    emit:
        downsample_dataset_accession_paired.out
}

workflow get_base_datasets {
    main:
        dataset_csv = file(params.dataset_csv, type: "file", checkIfExists:true)
        subset_dataset_accessions(dataset_csv, params.dataset_sample_size)
        subset_dataset_accessions.out.splitCsv(header: true).map{row -> ["${row.public_database_accession}","${row.platform}","${row.index}","${row.human_filtered_reads_1}","${row.human_filtered_reads_2}"]}.set{ dataset_accessions }

        dataset_accessions.branch { accession, platform, index, reads1, reads2 ->
            paired: platform == "illumina"
                return tuple(accession, platform, index, reads1, reads2)
            unpaired: true
                return tuple(accession, platform, index, reads1)
            }.set { by_platform }

        download_and_downsample_paired(by_platform.paired)
        download_and_downsample(by_platform.unpaired)
        download_and_downsample.out.concat(download_and_downsample_paired.out).set{ datasets }
        //datasets.view()
    emit:
        datasets
}

workflow get_references {
    main:
        reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
        subset_reference_accessions(reference_csv, params.reference_sample_size)
        subset_reference_accessions.out.splitCsv(header: true).map { row -> tuple("${row.accession}","${row.category_id}", "${row.index}") }.set{ reference_accessions }
        reference_accessions.tap{ to_download }

        download_reference_fasta(to_download.map{ accession, category, index -> [accession, category] }.unique())
        reference_accessions.combine(download_reference_fasta.out, by: 0).map{ accession, category, index, category1, fasta -> [accession, category, index, fasta]}.set{downloaded}
        //downloaded.view()
    emit:
        downloaded
}

workflow {
    main:
        reference_fasta = get_references()
        reference_fasta.view()
        base_datasets = get_base_datasets()
        base_datasets.view()
}
