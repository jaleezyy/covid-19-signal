
# export
export LMAT_DIR=/workspace/tsangkk2/LMAT-1.2.6/runtime_inputs

# run lmat
bash  /workspace/tsangkk2/LMAT-1.2.6/bin/run_rl.sh \
--db_file=/workspace/tsangkk2/LMAT-1.2.6/data/kML+Human.v4-14.20.g10.db \
--query_file=assembly.tiled.fasta \
--odir=output \
--overwrite --verbose \
--threads=20
