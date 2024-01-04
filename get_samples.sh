#!/usr/bin/env sh

# Orig
#parallel fasterq-dump {} ::: SRR020080 SRR020135

# Study: High-Throughput Methylome Analysis Reveals Differential Methylation During for
#        Early and Late Onset Preeclampsia for Mothers and Their Newborn Children
# Name: GSM7848911


mkdir -p samples
parallel fasterq-dump -O samples {} ::: SRR26440518 SRR26440528
