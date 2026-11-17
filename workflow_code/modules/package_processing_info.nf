process PACKAGE_PROCESSING_INFO {
    executor 'local'
    cpus 1
    memory '256 MB'
    cache false
    errorStrategy 'terminate'

    input:
        val(records)
        val(suffix)

    output:
        path("processing_info${suffix}.zip"), emit: zip

    exec:
        WorkflowProvenance.packageRecords(records, suffix, task.workDir.resolve("processing_info${suffix}.zip"))
}
