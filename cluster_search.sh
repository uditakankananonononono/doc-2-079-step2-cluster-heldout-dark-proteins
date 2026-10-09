#!/bin/bash
# DOC-2-079 step 2: sequence clustering and remote-homology search with MMseqs2 (version string recorded in run log).
set -e; M=${MMSEQS:-mmseqs}
for A in armA armB; do
  $M easy-cluster $A.fasta ${A}_clu tmp_$A --min-seq-id 0.4 -c 0.8 --cov-mode 1 -v 1
done
$M easy-search armB.fasta armB.fasta armB_hits.m8 tmp_search -s 7.5 -e 1e-3 --min-seq-id 0.2 --max-seqs 2000 --format-output query,target,fident,evalue -v 1
echo CLUSTER_SEARCH_DONE
