
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
params.no_ml = false
params.fast = false // if true, skip PanACoTA
params.ml_algorithms = "rf,svm,gb,histgb,lr,xgb"

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
    if os.path.exists('random_names.txt'):
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
    
    input:
    tuple val(strain_id), path(fasta)
    path bakta_db
    
    output:
    tuple val(strain_id), path("${strain_id}"), emit: out_dir
    
    script:
    def db_cmd = bakta_db ? "--db ${bakta_db}" : ""
    """
    bakta $db_cmd --output ${strain_id} --prefix $strain_id --threads ${params.threads} $fasta
    """
}

process MERGE_RESULTS {
    publishDir "${params.output}", mode: 'copy', pattern: "*.tsv"
    publishDir "${params.output}", mode: 'copy', pattern: "mutations_annotations.tsv"
    
    input:
    path "snippy_all/*"
    path "bakta_all/*"
    path "random_names.txt"
    path "strains.txt"
    path "phenotype_table.tsv"
    path reference
    
    output:
    path "binary_mutation_table.tsv", emit: mutation_table
    path "mutations_annotations.tsv", emit: annotations
    path "binary_mutation_table_with_gene_presence_absence.tsv", optional: true, emit: full_table
    path "random_names.txt", emit: names_dict
    path "strains.txt", emit: strains_list
    path "phenotype_table.tsv", emit: phenotype_table
    
    script:
    def qc_cmd = params.run_qc ? "--run_qc" : ""
    """
    mkdir -p snippy ${params.annotation_tool}
    # Move collected results into expected structure for alpar
    cp -r snippy_all/* snippy/
    cp -r bakta_all/* ${params.annotation_tool}/
    
    mkdir -p temp
    echo "2" > temp/status.txt
    
    alpar create_binary_tables -i strains.txt -o . --reference $reference --checkpoint \
        --threads ${params.threads} --temp temp --annotation_tool ${params.annotation_tool} \
        --gene_presence_absence_analysis_tool ${params.gpa_tool} $qc_cmd
    """
}

process THRESHOLDING {
    publishDir "${params.output}/binary_table_threshold", mode: 'copy'
    
    input:
    path mutation_table
    
    output:
    path "*.tsv", emit: thresholded_table
    
    script:
    """
    alpar binary_table_threshold -i $mutation_table -o .
    """
}

process PANACOTA {
    publishDir "${params.output}/panacota", mode: 'copy'
    
    input:
    path strains_list
    path names_dict
    
    output:
    path "phylogenetic_tree.newick", emit: tree
    
    script:
    """
    alpar panacota -i $strains_list -o . --threads ${params.threads} --random_names_dict $names_dict
    """
}

process MASH_TREE {
    publishDir "${params.output}/phylogeny", mode: 'copy'
    
    input:
    path strains_list
    path names_dict
    
    output:
    path "phylogenetic_tree.tree", emit: tree
    
    script:
    """
    alpar phylogenetic_tree -i $strains_list -o . --random_names_dict $names_dict
    """
}

process PRPS {
    publishDir "${params.output}/prps", mode: 'copy'
    
    input:
    path binary_table
    path tree
    
    output:
    path "prps_score.tsv", emit: prps_score
    
    script:
    """
    alpar prps -i $binary_table -o . -t $tree --threads ${params.threads}
    """
}

process GWAS {
    publishDir "${params.output}/gwas", mode: 'copy'
    
    input:
    path binary_table
    path phenotype_table
    path tree
    
    output:
    path "gwas_results/*", emit: results
    path "similarity_matrix.tsv"
    path "genotype_matrix.tsv"
    
    script:
    """
    alpar gwas -i $binary_table -p $phenotype_table -t $tree --threads ${params.threads} -o .
    """
}

process ML {
    publishDir "${params.output}/ml_results/${algorithm}/${antibiotic}", mode: 'copy'
    tag "${algorithm}_${antibiotic}"
    
    input:
    path binary_table
    path phenotype_table
    path prps_score
    path annotations
    path strains_list
    val antibiotic
    val algorithm
    
    output:
    path "*"
    
    script:
    def datasail = params.fast ? "" : "--sail $strains_list"
    """
    alpar ml -i $binary_table -o . -p $phenotype_table -a $antibiotic \
        --prps $prps_score --annotation $annotations \
        --ml_algorithm $algorithm --threads ${params.threads} \
        --save_model --parameter_optimization --feature_importance_analysis $datasail
    """
}

workflow {
    INITIALIZE(input_dir_ch)
    
    strains_ch = INITIALIZE.out.strains_info
        .splitCsv()
        .map { row -> tuple(row[0], file(row[1])) }
        
    RUN_SNIPPY(strains_ch, ref_ch)
    
    if (params.annotation_tool == 'bakta') {
        RUN_BAKTA(strains_ch, bakta_db_ch.first().ifEmpty([]))
        annot_out = RUN_BAKTA.out.out_dir.map{ it[1] }.collect()
    } else {
        // Placeholder for prokka if needed
        annot_out = Channel.empty().collect()
    }
    
    MERGE_RESULTS(
        RUN_SNIPPY.out.out_dir.map{ it[1] }.collect(),
        annot_out,
        INITIALIZE.out.names_dict,
        INITIALIZE.out.strains_list,
        INITIALIZE.out.phenotype_table,
        ref_ch
    )
    
    THRESHOLDING(MERGE_RESULTS.out.full_table.ifEmpty(MERGE_RESULTS.out.mutation_table))
    
    if (!params.fast) {
        PANACOTA(MERGE_RESULTS.out.strains_list, MERGE_RESULTS.out.names_dict)
        tree_ch = PANACOTA.out.tree
    } else {
        MASH_TREE(MERGE_RESULTS.out.strains_list, MERGE_RESULTS.out.names_dict)
        tree_ch = MASH_TREE.out.tree
    }
    
    PRPS(THRESHOLDING.out.thresholded_table, tree_ch)
    
    GWAS(THRESHOLDING.out.thresholded_table, MERGE_RESULTS.out.phenotype_table, tree_ch)
    
    if (!params.no_ml) {
        // Extract antibiotics from phenotype table
        antibiotics_ch = MERGE_RESULTS.out.phenotype_table
            .splitCsv(sep: '\t', limit: 1)
            .flatMap { it.drop(1) } // drop Strain column
            .filter { !it.startsWith('.') }
            
        algorithms_ch = Channel.fromList(params.ml_algorithms.split(',').toList())
        
        ML(
            THRESHOLDING.out.thresholded_table,
            MERGE_RESULTS.out.phenotype_table,
            PRPS.out.prps_score,
            MERGE_RESULTS.out.annotations,
            MERGE_RESULTS.out.strains_list,
            antibiotics_ch,
            algorithms_ch
        )
    }
}
