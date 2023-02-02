#!/bin/bash

# change directory
cd "gitlab/uASPIre_UTR_CDS"

# execute scripts
for i in {1..5}
do
  printf "Executing script ${i} ...\n"
  Rscript "./Figure_0${i}.R"
done
