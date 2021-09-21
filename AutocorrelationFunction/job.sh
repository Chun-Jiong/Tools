#!/bin/bash

if [[  $# -ne 1 ]]; then
	echo "enter the name of file."
	exit 1
fi

file=$1
if ! [ -f $file ]; then
	echo "No file: $file"
	exit 2
fi

tn=$(wc -l $file| awk '{print $1}')

for i in $(seq 1000 1000 $tn); do
	echo "$i $(./a.out $file $i)"
done

exit 0
