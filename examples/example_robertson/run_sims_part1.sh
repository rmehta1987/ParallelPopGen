#!/bin/bash
#input="needtoknow.txt"
if ! [ -x "$(command -v sinfo)" ]; then
  export PATH=$PATH:/home/rahul/PopGen/ParallelPopGen-0.3.2/examples/example_dadi/GOFish_57_epoch
  echo 'on local computer'
  echo $PATH
else
  nameofhost=hostname
  if [$nameofhost =~ "midway2"]
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

files=$( ls $1/*_1*.txt ) #Add () to convert output to array
for i in $files; do
    input=$i
    echo "reading file $input"

while IFS= read -r line
do
  line="$line"
  filename="sfs_out_selection_{$line}.txt"
  #echo "Output fileclear name is $filename"
  ./GOFISH_uk_robertson $line 360388 $filename
done < "$input"

done
