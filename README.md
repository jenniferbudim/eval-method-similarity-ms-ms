
# Reproduced an Evaluation Methodology For Machine Learning-Based Tandem Mass Spectra Similarity Prediction

> **Michael Strobel, Alberto Gil-de-la-Fuente, Mohammad Reza Zare Shahneh, Yasin El Abiead, Roman Bushuiev, Anton Bushuiev, Tomáš Pluskal, Mingxun Wang**  
> *BioMed Central Bioinformatics, 2025*  
> [Paper](https://rdcu.be/eQkWY) | [Video](https://example.com) | Original Source Code ([Preprocessing Pipeline](https://github.com/Wang-Bioinformatics-Lab/gnps_ml_processing_workflow)) ([Benchmark](https://github.com/Wang-Bioinformatics-Lab/Baselines_For_Benchmark))

---

# Overview

This repository contains a reproduction attempt for the experiments from the BMC Bioinformatics 2025 publication. The project focuses on creating a benchmark for the basis of development for future machine learning approaches in MS/MS similarity.

The Preprocessing Pipeline is used to harmonize the dataset's metadata (data is taken from GNPS Libraries and MassBank EU) and discard entries that doesn't follow the rules. The data is then split and prebatched for the Train and Test process. 

The Benchmarking Process is done by Training the MS2DeepScore model with the All-Pairs dataset (proposed method) and testing the all_pairs model using filtered data to minimize overhead. 

---
## Preprocessing Pipeline and Baseline For Benchmark

- **Processor Infrastructure:** Chameleon Cloud CPU (Intel(R) Xeon(R) Gold 6126 CPU @ 2.60GHz) at CHI@UC
- **Platform:** https://chameleoncloud.org/
- **Source Data:** ([GNPS](https://gnps.ucsd.edu/ProteoSAFe/libraries.jsp)) ([MassBank EU](https://massbank.eu/MassBank/))
- **Library Depedencies:**
    * Python
    * Conda (preferably Mamba)
    * Nextflow
      * JDK 17


## Running the Workflows

1. Preprocessing GNPS and MassBank EU Data

```bash
  cd ~/gnps_ml_processing_workflow/GNPS_ML_Processing
  mamba env create --file ./bin/conda_env.yml --prefix ./bin/gnps_ml_processing_env2/
  mamba env create --file ./bin/gnps_ml_processing_matchms.yml --prefix ./bin/gnps_ml_processing_matchms_env/
  nohup nextflow run nf_workflow.nf --subset=Structural_Similarity_Prediction -bg &
```

2. Split and Prebatch Dataset

```bash
  cd ~/gnps_ml_processing_workflow/Train_Test_Splits
  nextflow run Dataset_Splitting.nf
```
    
3. Training MS2DeepScore Model with All Pairs Dataset

```bash
  cd ~/Baselines_For_Benchmark
  ./build_environment.sh
  cd ~/scripts/ms2deepscore/train/fresh_dataset/
  ./train_all_pairs.sh
```

4. Testing Model

```bash
  cd ./scripts/ms2deepscore/train/fresh_dataset/
  ./all_pairs_to_all_filters.sh
```
## Acknowledgements

 - Computing Resources: Chameleon Cloud (CHI@UC) for providing CPU infrastructure
 - Original Authors: For making the training code publicly available
 - BMC Bioinformatics 2025: For publishing the foundational research

