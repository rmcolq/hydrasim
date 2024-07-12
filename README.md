# hydrasim
Hybrid metagenome simulator

User inputs:
1. References CSV containing `accession` and `category_id` columns
2. Datasets CSV, indexed by params.dataset_index and containing a `platform` column. If `human_filtered_reads_1` and `human_filtered_reads_1` columns are provided should use the reads at those paths instead of downloading them
3. Other parameters specified at the top of the modules files

Expected outputs:
The output directory will be populated with folders labelled by `category_id`. Within each folder there will be simulated dataset files named `${dataset_accession}_${ref_accession}_${simulator}_${ref_coverage}x_${index}[_R1,_R2,].fq.gz` where the index numbers the `params.num_iterations` simulations for each ref, dataset, simulator and coverage.

Expected behaviour:
Note that this pipeline makes use of a store directory, where datasets and reference fasta are stored so they can be called on several times without needing to download again
Also note that in places maxForks has been used where having more than 1 fork led to errors either due to clashes over reading a file by rasusa, or errors when downloading if too many simultaneous calls.
