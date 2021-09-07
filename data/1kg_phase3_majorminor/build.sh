#!/bin/bash

BUILD_VCFS=1
if [[ $BUILD_VCFS == 1 ]]
then

    # download full 1000 genomes phase 3 site vcf (hg19)
    wget -N \
    ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502//ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz

fi

BUILD_TSVS=1
if [[ $BUILD_TSVS == 1 ]]
then

    # make output dirs
    mkdir -p hg18/ hg19/ hg38/

    # MAF thresholds for file generation
    THRESHES=(0.001 0)

    # for each MAF threshold .. 
    for maf_thresh in ${THRESHES[@]}
    do  
    
        # build tsv for hg19 
        gunzip -c ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz \
        | ./1kg_vcf_to_table.py stdin ${maf_thresh} hg19/1000_genomes_phase3.hg19.maf_gt_${maf_thresh}
    
    done

fi

exit
