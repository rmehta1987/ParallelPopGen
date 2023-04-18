#!/bin/bash
if ! [ -x "$(command -v sinfo)" ]; then
  #export PATH=$PATH:/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch
  acmd='/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch'
  echo 'on local computer'
  echo $PATH
else
  nameofhost = hostname
  if [[$nameofhost =~ "midway2"]]
  then
    #export PATH=$PATH:/project2/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch
    acmd='/project2/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch'
    echo 'on midway2 cluster'
    echo $PATH
  else
    #export PATH=$PATH:/project/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch
    acmd='/project/jjberg/mehta5/ParallelPopGen/examples/example_dadi/GOFish_57_epoch'
    echo 'on midway3 cluster'
    echo $PATH
    fi
fi

for file in more_samples_part2_less_than_point_9/*.txt; do
while read -r line
do
  # make sure negative selection
  if [[ ${line:0:1} == "-" ]]
  then
   filename="sfs_out_selection_{$line}.txt" 
   #echo "Output file name is $filename"
   #echo "Already negative selection"
  $acmd $line 111710 $filename
  else
    line="-$line"
    filename="sfs_out_selection_{$line}.txt" 
    #echo "Output file name is $filename"
    $acmd $line 111710 $filename
  fi
  done < "$file"
  echo 'Finished File'
  echo $file
done 
