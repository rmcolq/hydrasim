#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.reference_csv = ""
params.reference_sample_size = 1
params.reference_dir = "store/references"

params.dataset_csv = ""
params.dataset_sample_size = 3
params.dataset_index = "public_database_accession"
params.dataset_coverage = "10k"
params.dataset_dir = "store/datasets"

params.simulators = ["badread","wgsim"]
lookup = ["illumina": "wgsim", "ont": "badread"]
params.coverages = [0.1, 1, 10, 100]
params.output_dir = "output" // Output directory

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

process simulate_reads {
    input:
        tuple val(simulator), val(coverage), val(accession), val(category), val(index), path(fasta)

    output:
        tuple val(accession), val(category), val(index), val(read_type), val(coverage), path("${accession}_${index}_${simulator}_${coverage}x_reads*.fq.gz")

    script:
    read_type = "other"
    if ("${simulator}" == "badread") {
        read_type = "ont"
        """
        badread simulate --reference ${fasta} --quantity ${coverage}x > ${accession}_${index}_${simulator}_${coverage}x_reads.fq
        gzip ${accession}_${index}_${simulator}_${coverage}x_reads.fq
        """
    } else if ("${simulator}" == "wgsim") {
        read_type = "illumina"
        """
        genome_length=\$(cat ${fasta} | wc -c)
        num_reads=\$(echo "\$genome_length*${coverage}/300" | bc)
        echo "\$genome_length \$num_reads"
        wgsim -N \$num_reads -1 150 -2 150 -r 0.001 -R 0.15 ${fasta} ${accession}_${index}_${simulator}_${coverage}x_reads_1.fq ${accession}_${index}_${simulator}_${coverage}x_reads_2.fq
        gzip ${accession}_${index}_${simulator}_${coverage}x_reads_1.fq ${accession}_${index}_${simulator}_${coverage}x_reads_2.fq
        """
    } else {
        error "Unsupported simulator: ${simulator}"
    }
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

workflow get_reference_fastas {
    main:
        reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
        subset_reference_accessions(reference_csv, params.reference_sample_size)
        subset_reference_accessions.out.splitCsv(header: true).map { row -> tuple("${row.accession}","${row.category_id}", "${row.index}") }.set{ reference_accessions }
        reference_accessions.tap{ to_download }

        download_reference_fasta(to_download.map{ accession, category, index -> [accession, category] }.unique())
        reference_accessions.combine(download_reference_fasta.out, by: 0).map{ accession, category, index, category1, fasta -> [accession, category, index, fasta]}.set{downloaded}
    emit:
        downloaded
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

workflow simulate_reference_reads {
    take:
        references
    main:
        simulators = channel.from(params.simulators)
        coverages = channel.from(params.coverages)
        simulators.combine(coverages).combine(references).set{ combinations }
        combinations.view()
        simulate_reads(combinations)
}

workflow {
    main:
        get_reference_fastas()
        //base_datasets = get_base_datasets()
        //base_datasets.view()
        simulate_reference_reads(get_reference_fastas.out)
}
