#!/bin/bash
docker run --rm \
  -v "$PWD/data-raw":/project/data-raw \
  -v "$PWD/R":/project/R \
  -v "$PWD/data":/project/data \
  -w /project \
  b2a_res_single_fiber \
  /opt/conda/bin/conda run -n snakemake_env --no-capture-output snakemake --cores 10