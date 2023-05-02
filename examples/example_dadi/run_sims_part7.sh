#!/bin/bash
#input="scripts_sample_sel_data/uniform_samples_5e_3_015_part2.txt"
if ! [ -x "$(command -v sinfo)" ]; then
  export PATH=$PATH:/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch
  echo 'on local computer'
  echo $PATH
else
  nameofhost = hostname
  if [[$nameofhost =~ "midway2"]]
  then
    export PATH=$PATH:/project2/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch
    echo 'on midway2 cluster'
    echo $PATH
  else
    export PATH=$PATH:/project/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch
    echo 'on midway3 cluster'
    echo $PATH
    fi
fi

files=$( ls $1/*_7*.txt ) #Add () to convert output to array
for i in $files; do
    input=$i
    echo "reading file $input"

while IFS= read -r line
do
  # make sure negative selection
  if [[ ${line:0:1} == "-" ]]
  then
   filename="sfs_out_selection_{$line}.txt"
   echo "Output file name is $filename"
   echo "Already negative selection"
  ./GOFish_57_epoch $line 111710 $filename
  else
    line="-$line"
    filename="sfs_out_selection_{$line}.txt"
    #echo "Output fileclear name is $filename"
    ./GOFish_57_epoch $line 111710 $filename
  fi
done < "$input"
done
