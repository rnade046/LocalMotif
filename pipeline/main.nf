nextflow.enable.dsl = 2
// input = java properties file
//params.properties = file("localMotifs.properties")
//params.localmotif_home = "src/main/java/"
//params.build_dir = "${workDir}/build/classes"

process compile {
    
    publishDir "${workDir}", mode: 'move'

    input:
    path srcdir 
    path jars

    output:
    path 'classes', emit: compiled_dir

    script:
    """
    mkdir -p "classes"
 
    # Build classpath
    CP=\$(echo ${jars.collect{ it.getName() }.join(' ')} | tr ' ' ':')
    
    # Compile all Java files
    javac -d classes -cp "\$CP" \$(find 'java/' -name '*.java')
    """
}

workflow  {

    Channel.fromPath('../src/main/java', type: 'dir').set { CH_SRCDIR }
    Channel.fromPath('../lib/*.jar').collect().set { CH_JARS }

    compile(CH_SRCDIR, CH_JARS)
}