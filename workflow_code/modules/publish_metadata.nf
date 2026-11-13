process PUBLISH_METADATA {
    publishDir path: { "${publishdir}/Metadata" },
        mode: params.publish_dir_mode,
        pattern: "metadata/*",
        saveAs: { filename -> filename.replaceFirst('^metadata/', '') }

    input:
        val(publishdir)
        path(runsheet)
        path(isa_archive)

    output:
        path("metadata/*")

    script:
        def copy_isa = isa_archive ? "cp -L '${isa_archive}' metadata/" : ''
        """
        mkdir -p metadata
        cp -L '${runsheet}' metadata/
        ${copy_isa}
        """
}
