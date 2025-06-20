#!/bin/bash

nextflow -log $(pwd)/.logs/nextflow.log \
run $(find $(pwd) -maxdepth 2 -type f -name "main.nf") \
-resume \
-work-dir $(pwd)/.work/ \
-with-report $(pwd)/.logs/run-report.html \
-with-trace $(pwd)/.logs/trace.txt \
-entry PREPROCESS \
--sample_table $(find $(pwd) -type f -name "sample_table.csv") \
--input_type "fastq" \
--genome "GRCh38" \
--star \
--screen