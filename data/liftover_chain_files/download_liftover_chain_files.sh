#!/bin/bash

### chain files to download :

## hg38 
#
#  hg38 -> hg19
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

## hg19
#
# hg19 -> hg38
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# hg19 -> hg18
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg18.over.chain.gz

## hg18
#
#hg18 -> hg38 
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
# hg18 -> hg19
wget -N http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz

exit
