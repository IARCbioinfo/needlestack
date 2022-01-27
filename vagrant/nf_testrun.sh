#!/bin/sh

data_dir=/vagrant/data_test

nextflow run iarcbioinfo/needlestack -with-docker \
                 --bed $data_dir/BED/TP53_all.bed \
                 --input_bams $data_dir/BAM/BAM_multiple/ \
                 --ref $data_dir/REF/17.fasta \
                 --output_vcf all_variants.vcf
