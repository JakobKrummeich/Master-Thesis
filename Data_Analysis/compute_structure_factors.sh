#!/bin/bash

sourcepath=$1

files=($sourcepath/State*.dat)

echo ${files[@]} | ./structure_factor $2
