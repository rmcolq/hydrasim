#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_csv = ""
params.reference_dir = "store/references"

params.dataset_csv = ""
params.dataset_index = "public_database_accession"
params.dataset_dir = "store/datasets"

params.read_types = ["illumina","ont"]
params.coverages = [0.1, 1, 10, 100]
params.num_iterations = 1
params.output_dir = "output" // Output directory

def scriptDir = workflow.projectDir

process subset_reference_accessions {
    input:
    path scriptDir
    path reference_csv
    val sample_size

    output:
    path "references.csv"

    script:
    """
    python ${scriptDir}/../bin/subset_accessions.py ${reference_csv} "references.csv" category_id ${sample_size} 
    """
}

process subset_dataset_accessions {
    input:
    path scriptDir
    path dataset_csv
    val sample_size

    output:
    path "dataset.csv"

    script:
    """
    python ${scriptDir}/../bin/subset_accessions.py ${dataset_csv} "dataset.csv" platform ${sample_size} 
    """
}

process download_reference_fasta {
    maxForks 1
    storeDir "${params.reference_dir}/${category}"

    input:
    tuple val(accession), val(category)

    output:
    tuple val(accession), val(category), path("${accession}_genomic.fna")

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

workflow get_base_fastq {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(index, accession, platform, reads1)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession.out, by: 0).map{ accession, platform, index, platform1, reads1 -> [index, accession, platform, reads1]}.set{downloaded}
        result.other.concat(downloaded).set{ base_fastq }
    emit:
        base_fastq
}

workflow get_base_fastq_paired {
    take:
        dataset_row
    main:
        dataset_row.branch { accession, platform, index, reads1, reads2 ->
            download: reads1 == ""
                return tuple(accession, platform, index)
            other: true
                return tuple(index, accession, platform, reads1, reads2)
            }.set { result }
        result.download.tap{ to_download }
        download_dataset_accession_paired(to_download.map{ accession, platform, index -> [accession, platform] }.unique())
        result.download.combine(download_dataset_accession_paired.out, by: 0).map{ accession, platform, index, platform1, reads1, reads2 -> [index, accession, platform, reads1, reads2]}.set{downloaded}
        result.other.concat(downloaded).set{ base_fastq }
    emit:
        base_fastq
}

workflow get_reference_fastas {
    main:
        reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
        subset_reference_accessions(scriptDir, reference_csv, params.num_iterations)
        subset_reference_accessions.out.splitCsv(header: true).map { row -> tuple("${row.accession}","${row.category_id}", "${row.index}") }.set{ reference_accessions }

        reference_accessions.tap{ to_download }
        download_reference_fasta(to_download.map{ accession, category, index -> [accession, category] }.unique())
        reference_accessions.combine(download_reference_fasta.out, by: 0).map{ accession, category, index, category1, fasta -> [index, accession, category, fasta]}.set{downloaded}
    emit:
        downloaded
}

workflow get_base_datasets {
    main:
        dataset_csv = file(params.dataset_csv, type: "file", checkIfExists:true)
        subset_dataset_accessions(scriptDir, dataset_csv, params.num_iterations)
        subset_dataset_accessions.out.splitCsv(header: true).map{row -> ["${row.public_database_accession}","${row.platform}","${row.index}","${row.human_filtered_reads_1}","${row.human_filtered_reads_2}"]}.set{ dataset_accessions }

        dataset_accessions.branch { accession, platform, index, reads1, reads2 ->
            paired: platform == "illumina"
                return tuple(accession, platform, index, reads1, reads2)
            unpaired: true
                return tuple(accession, platform, index, reads1)
            }.set { by_platform }

        get_base_fastq_paired(by_platform.paired)
        get_base_fastq(by_platform.unpaired)

    emit:
        paired = get_base_fastq_paired.out
        unpaired = get_base_fastq.out
}

workflow generate_recipes {
    main:
        get_reference_fastas()
        coverages = channel.from(params.coverages)
        get_reference_fastas.out.combine(coverages).set{ references }

        get_base_datasets()
        references.combine(get_base_datasets.out.paired, by: 0).set{ paired_recipes }
        references.combine(get_base_datasets.out.unpaired, by: 0).set{ unpaired_recipes }
        paired_recipes.view()
        unpaired_recipes.view()
     emit:
        paired = paired_recipes
        unpaired = unpaired_recipes
}
