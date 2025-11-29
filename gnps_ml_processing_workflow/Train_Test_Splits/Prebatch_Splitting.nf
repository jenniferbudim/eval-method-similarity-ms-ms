#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_LS = "$moduleDir/bin"

// Which task subset to use
params.subset = "Structural_Similarity_Prediction"
params.split_type = "sample_structures_smart_inchikey"//"structure_smart"  // 'basic_sampling_scheme', 'structure_smart', 'random' (random spectra), or 'structure' (random inchi14)
// params.split_type = "structure"
params.test_set_num = 4891

// Input data from Dataset_Splitting.nf
params.train_csv   = "nf_output/train_rows.csv"
params.val_csv     = "nf_output/val_rows.csv"
params.test_csv    = "nf_output/test_rows.csv"

params.train_mgf   = "nf_output/train_rows.mgf"
params.val_mgf     = "nf_output/val_rows.mgf"
params.test_mgf    = "nf_output/test_rows.mgf"

params.train_sim   = "nf_output/train_similarities.csv"
params.val_sim     = "nf_output/val_similarities.csv"
params.test_sim    = "nf_output/test_similarities.csv"
params.train_test_sim = "nf_output/train_test_similarities.csv"

// default params  
params.lowest_structural_threshold = 0.4
params.structural_similarity_fingerprint = "morgan"
params.thresholds = 0.5

// Batch Global Generation Parameters
params.num_epochs = 150
params.split_size = 10  // How many epochs each worker has, num_epochs should be divisible by this number
params.batch_size = 20
params.num_turns = 2  // Number of times each inchikey shows up per epoch
params.save_dir = "./nf_output/${params.subset}/${params.split_type}/"

// Which mass analyzers to include in filtered data.
params.mass_analyzer_lst = "" // Default is "", which is all mass analyzers. Otherwise, a semicolon delimited list of mass analyzers to include

// Prebatching Parameters
params.prebatch = true
params.unfiltered = false
params.filtered = true
params.include_biased = false
params.include_strict_ce = false

// Generate test sets
params.generate_test_set = true

// Parallelism Settings
params.parallelism = 3

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

process structural_similarity_calculation {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'
  
  cache true

  input:
  path train_rows_csv
  path test_rows_csv

  output:
  path 'output.csv', emit: structural_similarities

  """
  python3 $TOOL_FOLDER_LS/structural_similarity.py \
          --train_csv "$train_rows_csv" \
          --test_csv "$test_rows_csv" \
          --similarity_threshold $params.lowest_structural_threshold \
          --fingerprint "$params.structural_similarity_fingerprint"
  """
}

process split_data {
  publishDir './nf_output', mode: 'copy'
  cache false

  // Gather inputs avoiding name collisions
  input:
  path train_rows_csv, stageAs: 'train_rows.csv'
  path test_rows_csv, stageAs: 'test_rows.csv'
  path train_rows_mgf, stageAs: 'train_rows.mgf'
  path test_rows_mgf, stageAs: 'test_rows.mgf'
  path spectral_similarities
  path structural_similarities

  output:
  path 'train_rows_*.csv'
  path 'train_rows_*.mgf'
  path 'test_rows_*.csv'
  path 'test_rows_*.mgf'
  
  """
  python3 $TOOL_FOLDER_LS/split_data.py \
          --input_train_csv $train_rows_csv \
          --input_train_mgf $train_rows_mgf \
          --input_test_csv $test_rows_csv \
          --input_test_mgf $test_rows_mgf \
          --spectral_similarities "$spectral_similarities" \
          --strucutral_similarities "$structural_similarities" \
          --similarity_thresholds $params.thresholds  # No quotes around this param is necessary
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
  train_rows_ch = Channel.fromPath(params.train_csv)
  val_rows_ch   = Channel.fromPath(params.val_csv)
  test_rows_ch  = Channel.fromPath(params.test_csv)

  train_mgf_ch  = Channel.fromPath(params.train_mgf)
  val_mgf_ch    = Channel.fromPath(params.val_mgf)
  test_mgf_ch   = Channel.fromPath(params.test_mgf)

  train_sim_ch  = Channel.fromPath(params.train_sim)
  val_sim_ch    = Channel.fromPath(params.val_sim)
  test_sim_ch   = Channel.fromPath(params.test_sim)
  train_test_sim_ch = Channel.fromPath(params.train_test_sim)

  // To split the epochs
  epoch_ch = Channel.of(1..params.num_epochs/params.split_size)
  val_epoch_ch = Channel.of(1,)
  
  // Unfiltered
  if (params.prebatch && params.unfiltered) {
    unfiltered_train_sampler(
                        epoch_ch,
                        train_rows_ch, 
                        train_sim_ch, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "standard",  // mode
                                "",  // mass_analyzer_lst
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                false, // force_one_epoch
                                "train_unfiltered.hdf5", // output_name
                            ))
    unfiltered_train_assembler(unfiltered_train_sampler.out.prebatched_files.collect(), 
                              unfiltered_train_sampler.out.name_file.take(1))

    unfiltered_val_sampler(
                            val_epoch_ch,
                            val_rows_ch, 
                            val_sim_ch, 
                            tuple(  params.batch_size,
                                    params.num_turns,
                                    "standard",  // mode
                                    "",   // mass_analyzer_lst
                                    10,  // num_bins
                                    "", // exponential_bins
                                    false, // strict_collision_energy
                                    true, // force_one_epoch
                                    "val_unfiltered.hdf5", // output_name
                                ))
    unfiltered_val_assembler(unfiltered_val_sampler.out.prebatched_files.collect(), 
                              unfiltered_val_sampler.out.name_file.take(1))

  // Unfiltered Biased
    if (params.include_biased) {
      unfiltered_biased_train_sampler(
                          epoch_ch,
                          train_rows_ch, 
                          train_sim_ch, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "standard",  // filter
                                  "",  // mass_analyzer_lst
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  false, // force_one_epoch
                                  "train_unfiltered_biased.hdf5", // output_name
                              ))
      unfiltered_biased_train_assembler(unfiltered_biased_train_sampler.out.prebatched_files.collect(), 
                                unfiltered_biased_train_sampler.out.name_file.take(1))

      unfiltered_biased_val_sampler(
                          val_epoch_ch,
                          val_rows_ch,
                          val_sim_ch, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "standard",  // mode
                                  "",  // mass_analyzer_lst
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  true, // force_one_epoch
                                  "val_unfiltered_biased.hdf5", // output_name
                              ))
      unfiltered_biased_val_assembler(unfiltered_biased_val_sampler.out.prebatched_files.collect(),
                              unfiltered_biased_val_sampler.out.name_file.take(1))
    }
  }

  // Filtered
  if (params.prebatch && params.filtered) {
    filtered_train_sampler(
                        epoch_ch,
                        train_rows_ch, 
                        train_sim_ch, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "filter",  // mode
                                params.mass_analyzer_lst,
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                false, // force_one_epoch
                                "train_filtered.hdf5", // output_name
                            ))
    filtered_train_assembler(filtered_train_sampler.out.prebatched_files.collect(), 
                              filtered_train_sampler.out.name_file.take(1))

    filtered_val_sampler(
                        val_epoch_ch,
                        val_rows_ch,
                        val_sim_ch, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "filter",  // mode
                                params.mass_analyzer_lst,
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                true, // force_one_epoch
                                "val_filtered.hdf5", // output_name
                            ))
    filtered_val_assembler(filtered_val_sampler.out.prebatched_files.collect(),
                            filtered_val_sampler.out.name_file.take(1))
    // Filtered Biased
    if (params.include_biased) {
      filtered_biased_train_sampler(
                      epoch_ch,
                      train_rows_ch, 
                      train_sim_ch, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              20,  // num_bins
                              0.3, // exponential_bins
                              false, // strict_collision_energy
                              false, // force_one_epoch
                              "train_filtered_biased.hdf5", // output_name
                          ))
      filtered_biased_train_assembler(filtered_biased_train_sampler.out.prebatched_files.collect(), 
                                filtered_biased_train_sampler.out.name_file.take(1))  

      filtered_biased_val_sampler(
                          val_epoch_ch,
                          val_rows_ch,
                          val_sim_ch, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "filter",  // mode
                                  params.mass_analyzer_lst,
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  true, // force_one_epoch
                                  "val_filtered_biased.hdf5", // output_name
                              ))

      filtered_biased_val_assembler(filtered_biased_val_sampler.out.prebatched_files.collect(),
                            filtered_biased_val_sampler.out.name_file.take(1))

    }

    // Filtered Strict CE
    if (params.include_strict_ce) {
      filtered_strict_ce_train_sampler(
                      epoch_ch,
                      train_rows_ch, 
                      train_sim_ch, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              10,  // num_bins
                              "", // exponential_bins
                              true, // strict_collision_energy
                              false, // force_one_epoch
                              "train_filtered_strict_collision_energy.hdf5", // output_name
                          ))
      filtered_strict_ce_train_assembler(filtered_strict_ce_train_sampler.out.prebatched_files.collect(), 
                            filtered_strict_ce_train_sampler.out.name_file.take(1))

      filtered_strict_ce_val_sampler(
                      val_epoch_ch,
                      val_rows_ch,
                      val_sim_ch, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              10,  // num_bins
                              "", // exponential_bins
                              true, // strict_collision_energy
                              true, // force_one_epoch
                              "val_filtered_strict_collision_energy.hdf5", // output_name
                          ))

      filtered_strict_ce_val_assembler(filtered_strict_ce_val_sampler.out.prebatched_files.collect(),
                          filtered_strict_ce_val_sampler.out.name_file.take(1))
    }
  }

  if (params.generate_test_set) {
    // Unfiltered Test
    /*
    unfiltered_test_enumerator(test_rows_ch, 
                              test_sim_ch,
                              train_test_sim_ch,
                              tuple( false, // mode
                                      "", // mass_analyzer_lst
                                      false, // strict_collision_energy
                                      "unfiltered_test.parquet"
                              )
                              )
    */

    // Filtered Test
    filtered_test_enumerator(test_rows_ch, 
                            test_sim_ch,
                            train_test_sim_ch,
                            tuple( true, // mode
                                    params.mass_analyzer_lst, // mass_analyzer_lst
                                    false, // strict_collision_energy
                                    "filtered_test.parquet"
                            )
                            )

    if (params.include_strict_ce) {
      // Filtered Test (Strict CE)
      filtered_strict_ce_test_enumerator(test_rows_ch, 
                              test_sim_ch,
                              train_test_sim_ch,
                              tuple( true, // mode
                                      params.mass_analyzer_lst, // mass_analyzer_lst
                                      true, // strict_collision_energy
                                      "filtered_strict_collision_energy_test.parquet"
                              )
                              )
    }
  }
}