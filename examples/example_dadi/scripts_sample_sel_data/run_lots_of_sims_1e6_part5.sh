#!/bin/bash
input="scripts_sample_sel_data/uniform_samples_1e_6_1e_4_part5.txt"
if ! [ -x "$(command -v sinfo)" ]; then
  export PATH=$PATH:/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch
  echo 'on local computer'
  echo $PATH
else
  export PATH=$PATH:/project2/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch
  echo 'on cluster'
  echo $PATH
fi

while IFS= read -r line
do
  # make sure negative selection
  if [[ ${line:0:1} == "-" ]]
  then
   filename="sfs_out_selection_{$line}.txt" 
   #echo "Output file name is $filename"
   #echo "Already negative selection"
  ./GOFish_57_epoch $line 113770 $filename
  else
    line="-$line"
    filename="sfs_out_selection_{$line}.txt" 
    #echo "Output file name is $filename"
    ./GOFish_57_epoch $line 113770 $filename
  fi
done < "$input"