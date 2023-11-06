#!/bin/bash

nextflow run "./RNAseq_Nextflow/main.nf" -with-report -with-trace --genome 'GRCm38RMSKsalmon' --multiqc_report_title 'MultiQC Report Title' --sample_table "./RNAseq_Nextflow/sample_table_input.csv" -resume
