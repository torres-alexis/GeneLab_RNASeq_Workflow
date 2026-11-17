process CAPTURE_PROCESSING_RECORDS {
    executor 'local'
    cpus 1
    memory '256 MB'
    cache false
    errorStrategy 'terminate'

    input:
        val(completed_outputs)
        val(output_dir)
        val(suffix)
        val(run)
        val(samples)
        val(processed_directory)

    exec:
        WorkflowProvenance.save(output_dir, suffix, run, samples, file(processed_directory))
}
