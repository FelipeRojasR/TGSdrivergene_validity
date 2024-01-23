#!/bin/bash

args=("$@")
echo "First->" ${args[0]}
echo "Second->" ${args[1]}

# inside the for loop you can put the argument and iterate over every argument in the folder where the loop is directing.
mutations=${args[0]}
var2=${args[1]}

cd PATH/2020plus/
snakemake -s Snakefile pretrained_predict -p --cores 20 \
	--config mutations="${args[0]}" output_dir="${args[1]}" trained_classifier="PATH/2020plus/2020plus_100k.Rdata"
