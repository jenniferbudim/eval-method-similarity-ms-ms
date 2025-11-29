#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_LS = "$moduleDir/bin"

// Input data
params.input_csv = "/home/cc/nf_output/ML_ready_subset_positive/selected_summary.csv"
//"https://fileserver.wanglab.science/p_ml_cleanup/to_upload/cleaning_outputs/ML_ready_subset_positive/selected_summary.csv"
params.input_mgf = "/home/cc/nf_output/ML_ready_subset_positive/selected_spectra.mgf"
//"https://fileserver.wanglab.science/p_ml_cleanup/to_upload/cleaning_outputs/ML_ready_subset_positive/selected_spectra.mgf"

// Which task subset to use
params.subset = "Structural_Similarity_Prediction"
params.split_type = "sample_structures_smart_inchikey"//"structure_smart"  // 'basic_sampling_scheme', 'structure_smart', 'random' (random spectra), or 'structure' (random inchi14)
// params.split_type = "structure"
params.test_set_num = 4891

// Batch Global Generation Parameters
params.num_epochs = 150
params.split_size = 10  // How many epochs each worker has, num_epochs should be divisible by this number
params.batch_size = 32
params.num_turns = 2    // Number of times each inchikey shows up per epoch
params.save_dir = "./nf_output/${params.subset}/${params.split_type}/"

// Which mass analyzers to include in filtered data.
params.mass_analyzer_lst = "" // Default is "", which is all mass analyzers. Otherwise, a semicolon delimited list of mass analyzers to include

// Generate test sets
params.generate_test_set = true

// Parallelism Settings
params.parallelism = 24

// Submodule for batch generation
include { sampleSeveralEpochs as unfiltered_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_biased_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_biased_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_biased_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_biased_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_biased_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_biased_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_biased_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_biased_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_biased_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_biased_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_biased_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_biased_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'


include { exhaustivePairEnumeration as unfiltered_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'
include { exhaustivePairEnumeration as filtered_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'
include { exhaustivePairEnumeration as filtered_strict_ce_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'


process generate_subset {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path cleaned_csv
  path cleaned_mgf

  output:
  path "summary/*", emit: output_summary
  path "spectra/*.parquet", emit: output_parquet, optional: true
  path "spectra/*.mgf", emit: output_mgf
  path "json_outputs/*.json", emit: output_json, optional: true
  path "util/*", optional: true

  """
  python3 $TOOL_FOLDER_LS/subset_generator.py "$params.subset"
  """
}

process calculate_pairwise_similarity {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path metadata_csv

  output:
  path 'pairwise_similarities.csv', emit: pairwise_similarities

  """
  python3 $TOOL_FOLDER_LS/calculate_pairwise_similarity.py --input_metadata $metadata_csv --output_filename pairwise_similarities.csv
  """
}

/*
 * Generates a random subset of params.test_set_num  data points to use as a test set
 * These data points should be consistent across splitting methodology
 */
process generate_test_set {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path csv_file
  path pairwise_similarities
  path mgf_file

  output:
  path 'test_rows.csv',  emit: test_rows_csv
  path 'test_rows.mgf',  emit: test_rows_mgf
  path 'test_rows.json', emit: test_rows_json
  path 'test_similarities.csv', emit: test_similarities_csv
  path 'train_rows.csv', emit: train_rows_csv
  path 'train_rows.mgf', emit: train_rows_mgf
  path 'train_rows.json', emit: train_rows_json
  path 'train_similarities.csv', emit: train_similarities_csv
  path 'val_rows.csv', emit: val_rows_csv
  path 'val_rows.mgf', emit: val_rows_mgf
  path 'val_rows.json', emit: val_rows_json
  path 'val_similarities.csv', emit: val_similarities_csv
  path 'train_test_similarities.csv', emit: train_test_similarities_csv

  """
  python3 $TOOL_FOLDER_LS/calc_test.py \
          --input_csv "$csv_file" \
          --input_similarities "$pairwise_similarities" \
          --input_mgf "$mgf_file" \
          --num_test_points $params.test_set_num \
          --sampling_strategy "$params.split_type" \
          --threshold 0.5 \
          --debug
  """
}

process plot_data {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'

  cache false

  input:
  path test_rows_csv
  path train_test_similarities_csv
  path test_similarities_csv

  output:
  path '*.png'

  """
  python3 $TOOL_FOLDER_LS/plot_data.py \
          --test_csv "$test_rows_csv" \
          --train_test_similarities "$train_test_similarities_csv" \
          --test_similarities "$test_similarities_csv"
  """
}

// Outputs the results (avoids seding data to output if proceses are included in another script)
process output_handler {
  publishDir "./nf_output", mode: 'copy'

  input:
  path train_rows_csv
  path test_rows_csv
  path train_rows_mgf
  path test_rows_mgf

  output: 
  path "*", includeInputs: true

  """
  echo "Outputting data to nf_output"
  """
}

workflow {
  csv_file  = Channel.fromPath(params.input_csv)
  mgf_file  = Channel.fromPath(params.input_mgf)

  generate_subset(
    csv_file, 
    mgf_file
  )

  calculate_pairwise_similarity(
    generate_subset.out.output_summary
  )

  generate_test_set(
    generate_subset.out.output_summary, 
    calculate_pairwise_similarity.out.pairwise_similarities, 
    generate_subset.out.output_mgf
  )  

  // Generate Plot(s)
  plot_data(
    generate_test_set.out.test_rows_csv, 
    generate_test_set.out.train_test_similarities_csv, 
    generate_test_set.out.test_similarities_csv
  )
}