#!/bin/sh

# Run this in results to prepare the data for R processing of density data.
zcat complete/G.tsv.gz | tail -q -n +6 | grep All | gzip > G-all.tsv.gz
zcat complete/G.tsv.gz | tail -q -n +6 | grep -v All | gzip > G-classes.tsv.gz

zcat complete/*.tsv.gz | tail -q -n +6 | grep All | gzip > all-all.tsv.gz
(parallel 'zcat {} | tail -q -n +6 | grep -v All' ::: complete/*.tsv.gz)| gzip > all-classes.tsv.gz
