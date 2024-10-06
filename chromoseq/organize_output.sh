#!/bin/bash

for sample in `cat samples.txt`
do
	cd output_chromoseq/

	matching_files=$(ls | grep "$sample")

	if [ -n "$matching_files" ]; then
		mkdir -p $sample
		
		for file in $matching_files; do
            mv "$file" "$sample/"
        done		
	fi

	cd ../
done



