import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.time.OffsetDateTime
import java.util.regex.Pattern
import java.util.zip.ZipEntry
import java.util.zip.ZipOutputStream
import nextflow.Global
import nextflow.file.FileHelper
import nextflow.trace.config.TraceConfig

class WorkflowProvenance {
    static TraceConfig traceConfig() {
        def config = new TraceConfig(Global.session.config.trace as Map ?: [:])
        if (!config.enabled)
            throw new IllegalStateException('Processing-record capture requires trace.enabled = true')
        def fields = config.fields.collect { it.split(':', 2)[0].trim() }
        if (!fields.containsAll(['name', 'status', 'workdir']))
            throw new IllegalStateException('Processing-record capture requires trace.fields to include name, status and workdir')
        if (fields.any { it in ['script', 'env'] })
            throw new IllegalStateException('Processing-record capture does not support multiline script or env trace fields')
        if (!config.sep || config.sep.contains('\n') || config.sep.contains('\r'))
            throw new IllegalStateException('Processing-record capture requires a nonempty, single-line trace separator')
        return config
    }

    static Path traceFile() {
        def session = Global.session
        def config = traceConfig()
        def directory = config.hasProperty('directory') ? config.directory : null
        return directory
            ? session.outputDir.resolve(directory).resolve(config.file).normalize()
            : FileHelper.asPath(config.file)
    }

    static synchronized void begin(Path output, Map run) {
        def trace = traceFile()
        if (Files.exists(trace) && Files.lines(trace).withCloseable { lines -> lines.skip(1).anyMatch { !it.trim().isEmpty() } })
            throw new IllegalStateException('Execution trace already contains task records; use a new trace filename or enable trace.overwrite')
        def root = output.resolve('processing_info')
        def archives = root.resolve('runs')
        Files.createDirectories(archives)
        def current = root.resolve('workflow_run.json')
        def archive = archives.resolve(run.run_name as String)
        Files.createDirectory(archive)
        Files.writeString(archive.resolve('workflow_run.json'), JsonOutput.prettyPrint(JsonOutput.toJson(run + [capture_status: 'pending'])) + '\n')
        Files.copy(archive.resolve('workflow_run.json'), current, StandardCopyOption.REPLACE_EXISTING)
    }

    static void requireCurrent(Path root, Map run) {
        def current = new JsonSlurper().parseText(Files.readString(root.resolve('workflow_run.json')))
        if (current.run_name != run.run_name || current.session_id != run.session_id || current.capture_status != 'pending')
            throw new IllegalStateException('Processing records belong to another run; use a separate output directory')
    }

    static synchronized void finish(Path output, String runName, String sessionId, boolean success) {
        def root = output.resolve('processing_info')
        def current = root.resolve('workflow_run.json')
        if (!Files.isRegularFile(current)) return
        def run = new JsonSlurper().parseText(Files.readString(current))
        if (run.run_name != runName || run.session_id != sessionId) return
        def missingCapture = success && run.capture_status != 'captured'
        run.capture_status = success && !missingCapture ? 'complete' : 'failed'
        if (success && !missingCapture) run.completed = OffsetDateTime.now().toString()
        def content = JsonOutput.prettyPrint(JsonOutput.toJson(run)) + '\n'
        def archived = root.resolve('runs').resolve(runName).resolve('workflow_run.json')
        if (success || Files.isRegularFile(archived)) Files.writeString(archived, content)
        Files.writeString(current, content)
        if (missingCapture) throw new IllegalStateException('Processing-record capture did not complete')
    }

    static synchronized void save(Path output, String suffix, Map run, Path samples, Path processed) {
        def session = Global.session
        def root = output.resolve('processing_info')
        requireCurrent(root, run)
        def stats
        // Process outputs can close before Nextflow records the last task completion.
        for (int attempt = 0; attempt < 100; attempt++) {
            stats = session.statsObserver.stats
            if (stats.pendingCount == 0 && stats.submittedCount == 0 && stats.runningCount == 1) break
            Thread.sleep(100)
        }
        if (stats.pendingCount != 0 || stats.submittedCount != 0 || stats.runningCount != 1)
            throw new IllegalStateException('Analysis tasks have not finished before processing-record capture')
        int expected = stats.succeededCount + stats.cachedCount + stats.failedCount
        Path trace = traceFile()
        def separator = Pattern.quote(traceConfig().sep)
        List<String> header = []
        def rows = []
        boolean traceReady = false
        // The trace writer is asynchronous; wait for every completed/cached task.
        for (int attempt = 0; attempt < 100; attempt++) {
            def text = Files.exists(trace) ? Files.readString(trace) : ''
            def lines = text.endsWith('\n') ? text.readLines() : []
            if (lines.size() >= expected + 1) {
                header = lines.first().split(separator, -1).toList()
                if (!header.containsAll(['name', 'status', 'workdir']))
                    throw new IllegalStateException('Execution trace must include name, status and workdir')
                rows = lines.drop(1).findAll { it }.collect { it.split(separator, -1) }
                if (rows.size() == expected && rows.every { it.size() == header.size() }) {
                    traceReady = true
                    break
                }
            }
            Thread.sleep(100)
        }
        if (!traceReady)
            throw new IllegalStateException("Incomplete or malformed execution trace: ${trace}")

        def records = rows.collect { fields ->
            header.withIndex().collectEntries { key, i -> [(key): fields[i]] }
        }
        def scripts = records.collect { record ->
            if (!(record.status in ['COMPLETED', 'CACHED', 'FAILED']))
                throw new IllegalStateException("Unexpected execution trace status for ${record.name}: ${record.status}")
            def script = record.workdir && record.workdir != '-' ? FileHelper.asPath(record.workdir).resolve('.command.sh') : null
            def content
            if (script && Files.isRegularFile(script))
                content = Files.readString(script)
            else if (record.status == 'FAILED')
                content = '[shell script unavailable for failed task]\n'
            else
                throw new IllegalStateException("Missing shell script for ${record.name}: ${script ?: record.workdir}")
            "${record.name}\t${record.status}\n${content}"
        }.join('\n')
        def archive = root.resolve('runs').resolve(run.run_name as String)
        def captured = run + [capture_status: 'captured', output_root: output.relativize(processed).toString()]
        def files = [
            ("nextflow_run_command${suffix}.txt".toString()): "${run.command}\n",
            'samples.txt': Files.readString(samples),
            ("nextflow_processing_info${suffix}.txt".toString()): scripts,
            'workflow_run.json': JsonOutput.prettyPrint(JsonOutput.toJson(captured)) + '\n'
        ]
        files.each { name, content ->
            Files.writeString(archive.resolve(name), content as String)
        }
        requireCurrent(root, run)
        files.each { name, content ->
            Files.copy(archive.resolve(name), root.resolve(name), StandardCopyOption.REPLACE_EXISTING)
        }
    }

    static Map postRecord(Path output, String accession, List aliases = []) {
        def current = output.resolve('processing_info/workflow_run.json')
        if (!Files.isRegularFile(current)) return null
        def run = new JsonSlurper().parseText(Files.readString(current))
        def requested = accession?.toUpperCase()
        def requested_accessions = ([requested] + aliases).findAll { it }
            .collectMany { it.split(',').collect { id -> id.trim().toUpperCase() } }
        def matches = requested
            ? requested_accessions.any { it in [run.accession?.toUpperCase(), run.output_root?.toUpperCase()] }
            : !run.accession
        if (!matches) return null
        if (run.containsKey('capture_status') && run.capture_status != 'complete')
            throw new IllegalStateException("Processing records are ${run.capture_status}; complete the main workflow before post-processing")
        return run as Map
    }

    static void packageRecords(Path directory, String suffix, Path archive) {
        def names = ['samples.txt', "nextflow_processing_info${suffix}.txt", "nextflow_run_command${suffix}.txt"]
        if (Files.isRegularFile(directory.resolve('workflow_run.json'))) {
            def run = new JsonSlurper().parseText(Files.readString(directory.resolve('workflow_run.json')))
            if (run.containsKey('capture_status') && run.capture_status != 'complete')
                throw new IllegalStateException('Only completed processing records can be packaged')
            names.add('workflow_run.json')
        }
        new ZipOutputStream(Files.newOutputStream(archive)).withCloseable { zip ->
            names.each { name ->
                zip.putNextEntry(new ZipEntry("processing_info/${name}"))
                Files.copy(directory.resolve(name.toString()), zip)
                zip.closeEntry()
            }
        }
    }
}
