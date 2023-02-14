#!/bin/bash
input="uniform_samples_1e_8_1e_3.txt"
export PATH=$PATH:/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch
echo $PATH
while IFS= read -r line
do
   filename="sfs_out_selection_{$line}.txt"
   #echo $filename
  ./GOFish_57_epoch $line 113000 $filename
done < "$input"