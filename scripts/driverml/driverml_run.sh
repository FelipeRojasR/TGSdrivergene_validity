#!/bin/bash

args=("$@")
echo "First->" ${args[0]}
echo "Second->" ${args[1]}
#echo "Third->" ${args[2]}

# inside the for loop you can put the argument and iterate over every argument in the folder where the loop is directing.
input=${args[0]} ## PATH FOR INPUT DATA
output=${args[1]} ## FILE OUTPUT NAME
#summary=${args[2]} ## SUMMARY FILE NAME

cd PATH/DriverML-master/
bash PATH/DriverML-master/run.driverml.sh -w PATH/DriverML-master -i "${args[0]}" -f PATH/hg19.fa -r PATH/hg19.fa -m 10 -o "${args[1]}"
