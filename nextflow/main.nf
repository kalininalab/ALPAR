
nextflow.enable.dsl=2

/*
 * Automated Learning Pipeline for Antimicrobial Resistance (ALPAR)
 * Nextflow Implementation
 */

params.input = "input_folder"
params.output = "results"
params.reference = "reference.gbk"
params.bakta_db = null
params.threads = 4
params.ram = 16
params.annotation_tool = "bakta" // bakta or prokka

// Define channels
input_dir_ch = Channel.fromPath(params.input, type: 'dir')
ref_ch = Channel.fromPath(params.reference)
bakta_db_ch = params.bakta_db ? Channel.fromPath(params.bakta_db, type: 'dir') : Channel.empty()

process INITIALIZE {
    executor 'local'
    
    input:
    path input_dir
    
    output:
    path "strains.txt", emit: strains_list
    path "random_names.txt", emit: names_dict
    path "phenotype_table.tsv", emit: phenotype_table
    path "strains_info.csv", emit: strains_info
    
    script:
    """
    alpar create_phenotype_table -i $input_dir -o . --overwrite
    
    # Generate mapping for Nextflow
    python3 -c "
    import os
    names = {}
    with open('random_names.txt') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) == 2:
                names[parts[0]] = parts[1]
    
    with open('strains.txt') as f, open('strains_info.csv', 'w') as out:
        for line in f:
            path = line.strip()
            orig_name = os.path.splitext(os.path.basename(path))[0]
            if orig_name in names:
                random_name = names[orig_name]
                out.write(f'{random_name},{path}\\n')
    "
    """
}

process RUN_SNIPPY {
    tag "$strain_id"
    conda 'bioconda::snippy=4.6.0'
    
    input:
    tuple val(strain_id), path(fasta)
    path reference
    
    output:
    tuple val(strain_id), path("${strain_id}"), emit: out_dir
    
    script:
    """
    snippy --outdir ${strain_id} --ref $reference --ctgs $fasta --prefix $strain_id --cpus ${params.threads} --ram ${params.ram}
    """
}

process RUN_BAKTA {
    tag "$strain_id"
    conda 'bioconda::bakta=1.9.3'
    
    input:
    tuple val(strain_id), path(fasta)
    path bakta_db
    
    output:
    tuple val(strain_id), path("${strain_id}_bakta"), emit: out_dir
    
    script:
    def db_cmd = bakta_db ? "--db ${bakta_db}" : ""
    """
    bakta $db_cmd --output ${strain_id}_bakta --prefix $strain_id --threads ${params.threads} $fasta
    """
}

process MERGE_RESULTS {
    publishDir "${params.output}", mode: 'copy'
    
    input:
    path "snippy_all/*"
    path "bakta_all/*"
    path "random_names.txt"
    path "strains.txt"
    path "phenotype_table.tsv"
    path reference
    
    output:
    path "binary_mutation_table.tsv"
    path "mutations_annotations.tsv"
    path "binary_mutation_table_with_gene_presence_absence.tsv", optional: true
    
    script:
    """
    mkdir -p snippy bakta
    # Move collected results into expected structure
    find snippy_all -maxdepth 1 -type d -not -name 'snippy_all' -exec mv {} snippy/ \\;
    find bakta_all -maxdepth 1 -type d -not -name 'bakta_all' -exec mv {} bakta/ \\;
    
    # We need a dummy status.txt to skip tool runs in alpar
    mkdir -p temp
    echo "2" > temp/status.txt
    
    # Run alpar merge logic (using create_binary_tables in checkpoint mode)
    # We use a wrapper or modified alpar to avoid tool checks
    # For now, we assume alpar is available in the env
    export PYTHONPATH=\$PYTHONPATH:${projectDir}
    python3 -m sr_amr create_binary_tables -i . -o . --reference $reference --checkpoint --threads ${params.threads} --temp temp
    """
}

workflow {
    INITIALIZE(input_dir_ch)
    
    strains_ch = INITIALIZE.out.strains_info
        .splitCsv()
        .map { row -> tuple(row[0], file(row[1])) }
        
    RUN_SNIPPY(strains_ch, ref_ch)
    
    if (params.annotation_tool == 'bakta') {
        RUN_BAKTA(strains_ch, bakta_db_ch.first())
        bakta_out = RUN_BAKTA.out.out_dir.map{ it[1] }.collect()
    } else {
        // Handle prokka or skip
        bakta_out = Channel.empty().collect()
    }
    
    MERGE_RESULTS(
        RUN_SNIPPY.out.out_dir.map{ it[1] }.collect(),
        bakta_out,
        INITIALIZE.out.names_dict,
        INITIALIZE.out.strains_list,
        INITIALIZE.out.phenotype_table,
        ref_ch
    )
}
