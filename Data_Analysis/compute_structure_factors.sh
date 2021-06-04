#!/bin/bash

sourcepath=$1

files=($sourcepath/State*.dat)

filename=${files[0]}

[[ ${filename} =~ State_(.*_epsAB=[[:digit:]]*.[[:digit:]]*)_ ]] && extracted_file_info=${BASH_REMATCH[1]}


echo ${files[@]} | ./structure_factor $extracted_file_info $2

