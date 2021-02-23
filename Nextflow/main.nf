#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/proline
========================================================================================
 nf-core/compomics_workflow Analysis Pipeline.
#### Homepage / Documentation
TODO https://github.com/nf-core/compomics
----------------------------------------------------------------------------------------
*/


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    Usage:

    The typical command for running the pipeline is as follows:
    nextflow run main.nf --raws '*.raw' --fasta '*.fasta' --experiment_design 'test.txt'  -profile docker

    For testing purposes:
    nextflow run main.nf  -profile docker,test 

    Mandatory arguments:
      --raws                            Path to input data (must be surrounded with quotes)
      --fasta                Fasta file for database search
      -profile                          Configuration profile to use. Can use multiple (comma separated)
                                        Available: standard, conda, docker, singularity, awsbatch, test

    Mass Spectrometry Search:
      --peptide_min_length              Minimum peptide length for filtering
      --peptide_max_length              Maximum peptide length for filtering
      --precursor_mass_tolerance        Mass tolerance of precursor mass (ppm)
      --fragment_mass_tolerance         Mass tolerance of fragment mass bin (Da)
      --fions                          Forward ions for spectral matching
      --rions                          Reverse ions for spectral matching
//      --fdr_threshold                   Threshold for FDR filtering
//      --fdr_level                       Level of FDR calculation ('peptide-level-fdrs', 'psm-level-fdrs', 'protein-level-fdrs')
      --enzyme                          Enzymatic cleavage (e.g. 'Trypsin', see SearchGUI documentation)
      --miscleavages            Number of allowed miscleavages
//      --number_mods                     Maximum number of modifications of PSMs
      --fixed_mods                      Fixed modifications ('Carbamidomethyl of C', see SearchGUI modifications)
      --variable_mods                   Variable modifications ('Oxidation of M', see SearchGUI modifications)
     --min_charge                       Minimal precursor charge 
     --max_charge                       Maximal precursor charge 
      --skip_decoy_generation           Use a fasta databse that already includes decoy sequences
//      --quantification_fdr              Assess and assign ids matched between runs with an additional quantification FDR
      --run_xtandem                     SearchGui runs xtandem database search
      --run_msgf                        SearchGui runs msgf+ database search
      --run_comet                       SearchGui runs comet database search
      --run_ms_amanda                   SearchGui runs msamanda database search
      --run_myrimatch                   SearchGui runs myrimatch database search
      
    Options for flashLFQ:
      --mbr				true/false, apply match between runs
      --max_rt_alignment_shift          FlashLFQ parameter: Maximum MBR window (minutes) - The retention-time error allowed in match-between runs. 
      --experiment_design                text-file containing 2 columns: first with mzDB file names and second with names for experimental conditions

    Other options:
      --outdir                          The output directory where the results will be saved
      --email                           Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                             Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

    AWSBatch options:
      --awsqueue                        The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                       The AWS Region for your AWS Batch job to run on
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Validate inputs
params.raws = params.raws ?: { log.error "No read data provided. Make sure you have used the '--raws' option."; exit 1 }()
params.fasta = params.fasta ?: { log.error "No fasta file provided. Make sure you have used the '--fasta' option."; exit 1 }()
params.outdir = params.outdir ?: { log.warn "No output directory provided. Will put the results into './results'"; return "./results" }()


/*
 * Define the default parameters
 */

//MS params
params.peptide_min_length = 8
params.peptide_max_length = 12
params.fragment_mass_tolerance = 0.5
params.precursor_mass_tolerance = 30
params.fions = "b"
params.rions = "y"
//params.fdr_threshold = 0.01
//params.fdr_level = 'peptide-level-fdrs'
//fdr_level = (params.fdr_level == 'psm-level-fdrs') ? '' : '-'+params.fdr_level
//params.number_mods = 3

params.min_charge = 2
params.max_charge = 3

params.enzyme = 'Trypsin'
params.miscleavages = 1
params.fixed_mods = 'Carbamidomethylation of C'
params.variable_mods = 'Oxidation of M'
params.run_xtandem = 0
params.run_msgf = 0
params.run_comet = 0
params.run_ms_amanda = 0
params.run_myrimatch = 1
if (params.run_xtandem == 0  && params.run_msgf == 0 && params.run_comet == 0 && params.run_ms_amanda == 0 && params.run_myrimatch == 0) {
           log.error "No database engine defined. Make sure you have set one of the --run_searchengine options to 1 (searchengine can be xtandem, msgf, comet, ms_amanda, myrimatch)."; exit 1 
}

params.skip_decoy_generation = false
if (params.skip_decoy_generation) {
log.warn "Be aware: skipping decoy generation will prevent generating variants and subset FDR refinement"
log.warn "Decoys have to be named with DECOY_ as prefix in your fasta database"
}

params.experiment_design = "none"
// FlashLFQ parameters
params.mbr = "true"
params.max_rt_alignment_shift = 2.5

/*params.quantification_fdr = false
if (params.quantification_fdr) {
   log.warn "Quantification FDR enabled"
}
*/

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Configurable variables
params.name = false
params.email = false
params.plaintext_email = false


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input raw files
 */
Channel
        .fromPath( params.raws ).into {input_raw; input_raw2; input_raw3; input_raw4}

/*
 * Create a channel for fasta file
 */
  Channel
        .fromPath( params.fasta ).into {input_fasta; input_fasta2; input_fasta3}
  
        
/* 
 * Create a channel for proline experimental design file
 */
input_exp_design =  Channel.fromPath(params.experiment_design)
if (params.experiment_design == "none") {
    log.warn "No experimental design! All raw files will be considered being from the one and the same experimental condition."
} else if(!(file(params.experiment_design).exists())) {
        log.error "File with experimental design does not exit"; exit 1
    
}


/*
 * STEP 1 - convert raw files to mzML
 */
process convert_raw_mzml {
    publishDir "${params.outdir}"
    input:
      file rawfile from input_raw
    
    output:
     file "${rawfile.baseName}.mzML" into (mzmls, mzmls2, mzmls3, mzmls4)

    script:
     """
     ThermoRawFileParser.sh -i ${rawfile} -o ./  -f 1
     """

}

/*
 * STEP 2 - create  decoy database
 */
process create_decoy_database {
    publishDir "${params.outdir}"
    input:
      file fasta from input_fasta

    output:
      file "${fasta.baseName}_concatenated_target_decoy.fasta" into (fasta_with_decoy, fasta_with_decoy2)

    when:
      !params.skip_decoy_generation

    script:
     """
     searchgui eu.isas.searchgui.cmd.FastaCLI -in ${fasta} -decoy
     """    
}    



/*
 * STEP 3 - create  searchgui parameter file
 */
process create_searchgui_paramfile {
    publishDir "${params.outdir}"
    input:

      

    output:
      file "searchgui.par" into (searchgui_param, searchgui_param2)

    script:
     """
         searchgui eu.isas.searchgui.cmd.IdentificationParametersCLI -prec_tol ${params.precursor_mass_tolerance} -prec_ppm 1 \\
         -frag_tol ${params.fragment_mass_tolerance} -frag_ppm 0 -enzyme "${params.enzyme}" -mc ${params.miscleavages}  \\
             -fixed_mods "${params.fixed_mods}" -variable_mods "${params.variable_mods}" -min_charge ${params.min_charge} -max_charge ${params.max_charge} \\
         -fi "${params.fions}" -ri "${params.rions}" -import_peptide_length_min ${params.peptide_min_length} \\
         -import_peptide_length_max ${params.peptide_max_length} -xtandem_quick_acetyl 0 -xtandem_quick_pyro 0 \\
          -out searchgui.par 
         """    
} 



/*
 * STEP 4 - run database search
 */
process run_searchgui_search{
    publishDir "${params.outdir}"
    input:
      each file(mzmlfile) from mzmls
      file paramfile from searchgui_param
      file fasta_decoy from fasta_with_decoy.ifEmpty(input_fasta2)      

    output:
     tuple file("${mzmlfile.baseName}.zip"), file(mzmlfile) into (searchgui_out)

    script:
     """
	# needed for Myrimatch, see https://github.com/compomics/searchgui/issues/245
	LANG=/usr/lib/locale/en_US
	export LC_ALL=C; unset LANGUAGE
        mkdir tmp
        mkdir log       
        searchgui eu.isas.searchgui.cmd.PathSettingsCLI -temp_folder ./tmp -log ./log
         searchgui eu.isas.searchgui.cmd.SearchCLI -spectrum_files ./  -output_folder ./ -fasta_file "./${fasta_decoy}"  -id_params "./${paramfile}" -threads ${task.cpus} \\
         -xtandem ${params.run_xtandem} -msgf ${params.run_msgf} -comet ${params.run_comet} -ms_amanda ${params.run_ms_amanda} -myrimatch ${params.run_myrimatch}
         mv searchgui_out.zip ${mzmlfile.baseName}.zip

         """    
}

/*
 * STEP 5 - fdr, ... by PeptideShaker
 */
process run_peptideshaker {
    publishDir "${params.outdir}"
    input:
      tuple file(search_out), file(mzmlfile) from searchgui_out
      each file(fasta_decoy) from fasta_with_decoy2.ifEmpty(input_fasta3)      

    output:
      tuple file("${mzmlfile.baseName}.psdb"), file(mzmlfile) into (peptideshaker_file, peptideshaker_file2)

    script:
    mem = " ${task.memory}"
    mem = mem.replaceAll(" ","")
    mem = mem.replaceAll("B","")
     """
        mkdir tmp
        mkdir log    
        unzip ${search_out} searchgui.par
        peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.PeptideShakerCLI -spectrum_files "./${mzmlfile}"  -identification_files "./${search_out}"  -id_params searchgui.par \\
        -fasta_file "${fasta_decoy}" -reference "${params.name}" -out "./${mzmlfile.baseName}.psdb" -threads ${task.cpus} -Xmx${mem}
         """    
}

/*
 * STEP 6 - get full PSM report from peptideshaker results (tsv)
 */
// We might need to add the fasta file as input, as being outside the work folder
process get_peptideshaker_tsv {
    publishDir "${params.outdir}"
    input:
        tuple file(pepshaker), file(mzmlfile) from peptideshaker_file
	

    output:
      file "${mzmlfile.baseName}.txt"  into (peptideshaker_tsv_file)
      file "${mzmlfile.baseName}_filtered.txt"  into (peptideshaker_tsv_file_filtered)

    script:
     """
        peptide-shaker eu.isas.peptideshaker.cmd.PathSettingsCLI  -temp_folder ./tmp -log ./log
        peptide-shaker eu.isas.peptideshaker.cmd.ReportCLI -in "./${pepshaker}" -out_reports "./" -reports "3,4"
	mv "${params.name}_Default_PSM_Report_with_non-validated_matches.txt" "${mzmlfile.baseName}.txt"
	mv "${params.name}_Default_PSM_Report.txt" "${mzmlfile.baseName}_filtered.txt"
         """    
}

/*
 * STEP 7 - run flashLFQ for quantification
*/
process flashLFQ_all {
    echo true
    publishDir "${params.outdir}"

    input:
	file peptideshaker_out from peptideshaker_tsv_file_filtered.collect()
        file mzmlfiles from mzmls2.collect()

    output:
//        stdout stdout_channel
	file "QuantifiedPeaks.tsv" into flashlfq_peaks
	file "QuantifiedPeptides.tsv" into flashlfq_peptides
	file "QuantifiedProteins.tsv" into flashlfq_proteins


    script:
        """

       first_line=""
       for file in *.txt		
       	   do
	   echo \$file
	   tail -n +2 "\$file" >> tlfq_ident.tabular
	   first_line=\$(head -n1 "\$file")
        done
        echo "\$first_line" | cat - tlfq_ident.tabular > lfq_ident.tabular
        FlashLFQ --idt "lfq_ident.tabular" --rep "./" --out ./ --mbr ${params.mbr} --mrt ${params.max_rt_alignment_shift} --ppm ${params.precursor_mass_tolerance} --thr ${task.cpus}
        """
}

/*
 * STEP 9 - run MSqRob for stats
*/
process run_msqrob {
    publishDir "${params.outdir}"

    input:
      file exp_design from input_exp_design
      file rawfiles from input_raw4.collect()
      file quant_tab from flashlfq_peptides
       
    output:
      file "MSqRobOut.csv"  into msqrob_prot_out

    script:
   // no file provided
    expdesign_text = "run\tgenotype\tbiorep"
    if (exp_design.getName() == "none") {
        if(rawfiles[1] != null) {
                for( int i=0; i<rawfiles.size(); i++ ) {
                   biorep = i+1
                   expdesign_text += "\n${rawfiles[i].getBaseName()}\tMain\tA${biorep}"
                }
        } else {
          expdesign_text += "\n${rawfiles.getBaseName()}\tMain\tA1"
        }
    }

    """
    echo "${expdesign_text}" > none
    cp "${exp_design}" exp_design.tsv
    mv "${quant_tab}" q_input.txt
    Rscript $baseDir/runMSqRob.R
    """

 }




/*
 * STEP 9 - run PolySTest for stats

process run_polystest {
    publishDir "${params.outdir}"

    input:
      file exp_design from input_exp_design
      file moff_res from mbr_output
       
    output:
      file "polystest_prot_res.csv"  into polystest_prot_out
      file "polystest_pep_res.csv"  into polystest_pep_out

    script:
    """
    convertFromProline.R "${exp_design}" "${moff_res}"
    sed -i "s/threads: 2/threads: ${task.cpus}/g" pep_param.yml
    sed -i "s/threads: 2/threads: ${task.cpus}/g" prot_param.yml
    runPolySTestCLI.R pep_param.yml
    runPolySTestCLI.R prot_param.yml    
    """

 }
*/

workflow.onComplete {
    log.info ( workflow.success ? "\nDone! Open the files in the following folder --> $params.outdir\n" : "Oops .. something went wrong" )
}

       
