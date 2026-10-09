# DOC-2-079 step 2 RESULTS - review pending (independent gate has not cleared; no claim is final)

Labels (set mechanically by analysis_s2.py): **Arm A: NO-DEFICIT-ABOVE-MDG. Arm B: INCONCLUSIVE-POWER.** G1 controls pass in both arms.

## Arm A: same Pfam families, leave-cluster-out (40% identity) nearest neighbour
| item | value |
|---|---|
| sample | 990 proteins (328 dark, 662 characterized), 35 families, 531 clusters at 40% identity (expected about 47 families; fewer after exact-sequence dedup) |
| evaluable proteins (>= 1 same-family protein in another cluster) | dark 328, characterized 662 |
| ESM-2 1-NN family accuracy | dark 0.924, characterized 0.944 |
| composition baseline C (descriptive) | dark 0.366, characterized 0.295 |
| G1 permuted-family hit rate | dark 0.018, char 0.023; threshold 3/35 = 0.086; pass |
| gap (characterized minus dark) | +0.020, 95% family-bootstrap CI [-0.009, +0.053]; CI half-width 0.031 (<= 0.05); SE 0.016; MDG (80% power) = 0.045 |
| dark E minus C (reported only) | +0.558, CI [0.466, 0.660] |

Reading: with redundancy removed and neighbours forced into other 40% clusters, ESM-2 still recovers the Pfam family for 92% of dark and 94% of characterized proteins. No deficit larger than the MDG (about 0.045) was detectable for dark proteins; nothing is said below that. The point gap (+0.020) is positive but its CI includes 0.

## Arm B: no-Pfam dark vs characterized, truth = alignment-defined remote homologs
| item | value |
|---|---|
| sample | 3,000 proteins (1,500 dark no-Pfam, 1,500 characterized), 2,367 clusters |
| evaluable proteins (>= 1 remote homolog at 20-40% identity, e <= 1e-3, inside the sample) | **dark no-Pfam 41, characterized 901** |
| ESM-2 1-NN remote-homolog hit rate | dark no-Pfam 0.512 (21 of 41), characterized 0.657 |
| composition baseline C (descriptive) | dark no-Pfam 0.171, characterized 0.072 |
| G1 random eligible neighbour hit rate | dark 0.0, char 0.0 (threshold < 0.10); pass |
| gap | +0.145, CI [-0.014, +0.313]; CI half-width 0.164 (> 0.05); MDG = 0.234 |

Reading: **INCONCLUSIVE-POWER.** Only 41 of 1,500 no-Pfam dark proteins have an alignment-detectable remote homolog inside the sample, so the dark arm has too few evaluable proteins. The point estimate (dark 0.51 vs characterized 0.66) is a larger gap than arm A but its CI spans from -0.014 to +0.313; no deficit is claimed and no absence of deficit is claimed. This arm does not answer the no-Pfam question.

## Limits (carried from PROTOCOL.md)
- The gate-required numbers (families at >= 5 each, about 47) were expectations, not guarantees. Realised: 35 families in arm A and 41 evaluable dark proteins in arm B; an INCONCLUSIVE-POWER result is acceptable under the protocol.
- "Characterized" = EC-annotated only. Arm B characterized proteins may have Pfam while dark ones do not, so arm B compares annotation status together with Pfam coverage, not annotation status alone. Remote homology truth is alignment-defined at 20-40% identity with e <= 1e-3; it misses homologs below detectability and counts only pairs inside the sample, so arm B accuracies are relative, not absolute.
- Clusters are 40% identity; proteins at 30-40% identity can still be related, so leave-cluster-out is not leave-superfamily-out. MMseqs2 clustering and search depend on the version and parameters recorded in the run log (version string 564f40d8857f4eca4e1dfe100c67c155b1933e70).
- One model size, one pooling, single seed, CPU only. No simulated data in results; the smoke test of analysis_s2.py was on synthetic arms only, no record saved.

## Disclosures
- Run once. analysis_s2.py md5 fda894cdd907a2d67cf7572eb75949cc equals the lock-1 file; tag_tree_check.txt records the lock-1 tag tree (6 files). run_log.txt: UTC 17:37:29-17:37:34, exit 0, MMseqs2 564f40d8..., numpy 2.2.6, pandas 2.3.3, torch 2.14.1+cpu, transformers 5.19.0, python 3.10.12. input_md5.txt holds md5s of all inputs, embeddings, cluster/hit files and scripts. No amendments.
- The build_s2.py output line was lost to a shell timeout; counts above were recomputed from armA.tsv and armB.tsv (the files used by the run).
- Large files (armA/armB embeddings, FASTA, cluster and hit files, raw TSVs) are in the Drive folder.
