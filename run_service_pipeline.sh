#!/usr/bin/env bash
# I am running pipeline.py to test breaking the main functions of each
# .py file in the dvc pipeline to take/give objects when possible and avoid
# read/write intermediate files
src/motif_mining/pipeline.py \
	data/zhong2021/expression.csv \
  data/zhong2021/meta.csv \
	data/string/9606.protein.physical.links.detailed.v12.0.txt.gz \
  data/string/string_v12_idmap.csv \
	data/mining \
	test_pipeline_exec \
	--value npx \
	--var_id symbol

