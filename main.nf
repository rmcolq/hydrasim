include { generate_recipes } from './subworkflows/generate_recipes'
include { simulate_datasets } from './subworkflows/simulate_datasets'

process checkPath {    
    script:
    """
    echo "Current PATH: \$PATH" >temp.txt
    """
}


workflow {
    unique_id = "${params.unique_id}"

    // // Set Unique ID
    // if (unique_id == "null") {
    //     unique_id = "${params.timestamp}.out"
    // }
    // println "Unique ID: ${unique_id}"

    // Check input files
    if (params.reference_csv) {
        reference_csv = path(params.reference_csv, type: "file", checkIfExists:true)
    } else {
        exit 1, "Reference CSV be provided -- aborting"
    }
    if (params.dataset_csv) {
        dataset_csv = path(params.dataset_csv, type: "file", checkIfExists:true)
    } else {
        exit 1, "Dataset CSV be provided -- aborting"
    }

    println reference_csv
    println dataset_csv
    // Check input files and run 
    // simulate_datasets(unique_id)

}