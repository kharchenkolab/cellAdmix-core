#!/bin/bash
# Thinning ladder on the pancreas pseudo-replicates: for each spread of the
# gene-wise log-efficiency deviation applied to the right half, cut the
# bundles, run the separate, consensus, and joint arms with bridge scoring,
# and compare them. One ladder step at a time; results land under
# analysis/multisample/results/ladder_<sigma>/.
#
# Usage: run_ladder.sh <threads> <sigma> [<sigma> ...]
set -u
THREADS=$1; shift
WT=/home/pkharchenko/cellAdmix/wt-dev_multisample
D=/home/pkharchenko/cellAdmix/datasets
RUNS=/home/pkharchenko/cellAdmix/multisample_runs
R=/usr/local/R/R-4.4.3/bin/Rscript
PY=$WT/.venv-sdata/bin/python
SRC=/home/pkharchenko/baysor/rewrite/Baysor/examples/Xenium_pancreas_membrane_377/data
ANN=/home/pkharchenko/cellAdmix/cellAdmix-core/examples/xenium_pancreas_membrane_377_full/annotations/annotation.csv.gz
mkdir -p $RUNS/logs
for SIGMA in "$@"; do
  TAG=ladder_${SIGMA}
  DATA=$D/pancreas_${TAG}
  OUT=$RUNS/pancreas_${TAG}
  LOG=$RUNS/logs/${TAG}.log
  echo "$(date +%H:%M:%S) sigma=$SIGMA start" | tee -a $LOG
  $PY -I $WT/analysis/multisample/make_pseudo_replicates.py $SRC $DATA --sigma $SIGMA --seed 1 --annotation $ANN >> $LOG 2>&1
  for S in left right; do
    $R $WT/analysis/multisample/run_separate.R $S $DATA/$S $DATA/${S}_annotation.csv.gz $OUT/$S $THREADS >> $LOG 2>&1 &
  done
  wait
  SPEC="left=$DATA/left=$DATA/left_annotation.csv.gz=$OUT/left right=$DATA/right=$DATA/right_annotation.csv.gz=$OUT/right"
  $R $WT/analysis/multisample/run_collection.R $OUT/consensus $THREADS consensus $SPEC >> $LOG 2>&1 &
  $PY -I $WT/analysis/multisample/make_joint_table.py $D/joint/pancreas_${TAG} 50000 left=$DATA/left=$DATA/left_annotation.csv.gz right=$DATA/right=$DATA/right_annotation.csv.gz >> $LOG 2>&1 \
    && $PY -I -c "import pyarrow.parquet as pq; pq.read_table('$D/joint/pancreas_${TAG}_molecules.parquet').to_pandas().to_csv('$D/joint/pancreas_${TAG}_molecules.csv.gz', index=False)" >> $LOG 2>&1 \
    && $R $WT/analysis/multisample/run_joint.R $D/joint/pancreas_${TAG}_molecules.csv.gz $D/joint/pancreas_${TAG}_annotation.csv.gz $OUT/joint $THREADS >> $LOG 2>&1 &
  wait
  $R $WT/analysis/multisample/compare.R pancreas_${TAG} $RUNS $WT/analysis/multisample/results/${TAG} left=$DATA/left_annotation.csv.gz right=$DATA/right_annotation.csv.gz >> $LOG 2>&1
  echo "$(date +%H:%M:%S) sigma=$SIGMA done" | tee -a $LOG
done
