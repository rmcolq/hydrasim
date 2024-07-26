include { get_params_and_versions } from '../modules/get_params_and_versions'

workflow ingest {
    take:
        unique_id
    main:
        get_params_and_versions(unique_id)

}