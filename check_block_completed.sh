#! /usr/bin/bash

# uso {} [folders to check]

folder=$1
prefixout=$2

for i in $folder 
do 
  for j in $(ls $folder"/"$prefixout".out")
  do
    a=($(tail -n2 $j |head -n1))
    #echo $j" "${a[0]}" "${a[1]} 
    if [ ${a[0]}' '${a[1]} != 'JOB DONE.' ]
         then  
          echo '#### ERRORE '$j 
    fi
  done 
done







