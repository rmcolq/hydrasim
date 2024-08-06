// Check input files


    if (params.reference_csv) {
            reference_csv = file(params.reference_csv, type: "file", checkIfExists:true)
    } else {
        exit 1, "Reference CSV be provided -- aborting"
    }
    if (params.dataset_csv) {
            reference_csv = file(params.dataset_csv, type: "file", checkIfExists:true)
    } else {
        exit 1, "Dataset CSV be provided -- aborting"
    }