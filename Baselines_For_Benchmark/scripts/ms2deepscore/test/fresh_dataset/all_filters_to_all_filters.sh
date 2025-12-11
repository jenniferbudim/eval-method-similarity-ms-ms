#!/usr/bin/env bash
eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"

FOLDER="fresh_dataset"
DATASET="filtered_test"  

BIN_DIR="/home/cc/eval-method-similarity-ms-ms/Baselines_For_Benchmark/src/ms2deepscore"
DATA_DIR="/home/cc/eval-method-similarity-ms-ms/Baselines_For_Benchmark/data/$FOLDER/"
MODEL_PATH="./$FOLDER/model_04_12_2025_00_33_15.hdf5"

conda activate "/home/cc/eval-method-similarity-ms-ms/Baselines_For_Benchmark/src/shared/conda_env"

echo "Performing Inference for MS2DeepScore with a basic split, filtered by pairs"

cd $BIN_DIR

python3 test_presampled.py \
    --test_path $DATA_DIR"/processed/data/ALL_GNPS_positive_test_split.pickle" \
    --presampled_pairs_path $DATA_DIR"/raw/test_pairs/filtered_test.parquet/" \
    --save_dir "./$FOLDER" \
    --n_jobs 8 \
    --model_path $MODEL_PATH

DIR_NAME=$(dirname "$MODEL_PATH")
BASE_NAME_WITHOUT_EXT=$(basename "$MODEL_PATH" .hdf5)

python3 ~/eval-method-similarity-ms-ms/Baselines_For_Benchmark/src/shared/metrics_and_plotting.py \
    --prediction_path "$DIR_NAME/$BASE_NAME_WITHOUT_EXT/filtered_test/dask_output.parquet" \
    --save_dir "$DIR_NAME/$BASE_NAME_WITHOUT_EXT/filtered_test/" \
    --n_jobs 12