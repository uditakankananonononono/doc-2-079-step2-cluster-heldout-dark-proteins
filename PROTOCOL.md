# DOC-2-079 step 2: with redundancy removed, are "uncharacterized" proteins still organised by ESM-2 as well as characterized ones, and do no-Pfam dark proteins recover remote homologs as well? (frozen protocol, lock-1)

Written 2026-10-09 IST and committed BEFORE any step-2 data was downloaded, sampled, embedded or clustered. Follows the step-1 gate verdict (repo doc-2-079-dark-protein-family-organisation, step 1 label ORGANISED; no deficit detectable, ceiling, 4 errors in 340). The ledger has no spec text for DOC-2-079; this is a design response to the step-1 gate findings, not the full atlas.

## What step 1 did and why step 2 differs
Step 1 was at ceiling (98% / 99%), 17 families, 31 of 340 exact duplicates, and 93% of proteins had a neighbour with cosine > 0.95, so it could not tell dark from characterized. Step 2 changes the design in the five ways the gate required: (a) sequence clusters at 40% identity (MMseqs2) with leave-cluster-out nearest-neighbour scoring and exact duplicates removed; (b) relaxed family filter (>= 5 dark and >= 5 characterized); (c) a separate arm for dark proteins with NO Pfam, whose truth is not Pfam; (d) a pre-registered minimum detectable gap (MDG) rule that can return INCONCLUSIVE; (e) the dark vs characterized comparison stays within the same families in arm A. Not a re-ask of the earlier unit, which could not answer the question; it is its registered redesign.
Prior art: embeddings find remote homologs (e.g. embedding-based homology search) is known; no novelty claimed.

## Data (frozen)
- Frozen inputs: dark_raw.tsv from step 1 (md5 a9c57ef14d325e9f7bfc1dea650dabca, UniProt 2026_03), uniprot_raw.tsv from DOC-2-072 (md5 aa086b6e99d5e84c44789b1a85bc71f1; EC-annotated = characterized), and nopfam_raw.tsv fetched by acquire_nopfam.py (reviewed, name contains "uncharacterized", no Pfam, no EC, 100-400 aa; release and md5 recorded). Drop X/B/Z/U/O residues. Exact duplicate sequences removed.
- Arm A (same Pfam families): single-Pfam proteins; families with >= 5 dark and >= 5 characterized after dedup; up to 20 per group per family (seed 801).
- Arm B (no-Pfam dark): 1,500 dark no-Pfam proteins vs 1,500 characterized EC proteins (any Pfam; seed 802), dedup across both.
- Embeddings: ESM-2 35M mean-pooled, embed.py identical to DOC-2-072's (md5 8437eda4de9a3e4c05a98684e6cb9623), run on armA.tsv and armB.tsv separately.
- Clusters: MMseqs2 easy-cluster at --min-seq-id 0.4 -c 0.8 --cov-mode 1 per arm (version string recorded in the run log). Remote homology truth (arm B): MMseqs2 easy-search all-vs-all within arm B, -s 7.5, -e 1e-3, --min-seq-id 0.2; a remote homolog pair has identity in [0.20, 0.40) and e <= 1e-3; pairs with identity >= 0.40 are "close" and excluded.

## Method
- Arm A: for each protein, cosine 1-NN over all proteins in OTHER 40% clusters (same-cluster and self masked). Hit = same Pfam family. Only "evaluable" proteins (at least one same-family protein in another cluster) are scored. Accuracy by group (dark, characterized); composition baseline C (20 AA freq + log length, z-scored) descriptive.
- Arm B: for each protein, cosine 1-NN over proteins that are not close homologs, not in its own 40% cluster, not itself. Hit = that neighbour is a remote homolog of the query (alignment-defined truth, independent of Pfam and of embeddings). Evaluable = proteins with at least one remote homolog.
- Bootstrap: arm A resamples families; arm B resamples 40% clusters; 2,000 resamples, seed 12345, percentile 95%.

## Gates and labels (per arm, mechanical)
- G1 (control): arm A, permuted family labels give 1-NN hit rate < 3 / n_families for both groups. Arm B, a random eligible neighbour gives hit rate < 0.10 for both groups. Failure = INVALID, protocol stop.
- Gap = accuracy(characterized) - accuracy(dark or dark no-Pfam), E features. SE = bootstrap SD; CI half-width = half the 95% CI width; MDG (80% power, two-sided 5%) = 2.8 x SE.
- Label: INVALID if G1 fails; INCONCLUSIVE-POWER if the CI half-width of the gap exceeds 0.05; DARK-IS-DARKER if gap >= 0.05 and CI lower bound > 0; otherwise NO-DEFICIT-ABOVE-MDG (meaning: no deficit larger than the reported MDG was detectable; nothing is said below it).
- Reported only: arm A dark E minus C (embedding adds beyond composition), evaluable counts, cluster counts.

## Limits stated up front
- The gate-required numbers (families at >= 5 each: about 47; clusters; evaluable counts) are expectations, not guarantees; realised counts will be reported against them and an INCONCLUSIVE-POWER result is acceptable.
- "Characterized" = EC-annotated only; arm B characterized proteins may have Pfam while dark ones do not, so arm B compares annotation status together with Pfam coverage, not annotation status alone. Remote homology truth is alignment-defined at 20-40% identity with e <= 1e-3; it misses homologs below detectability and counts only pairs inside the sample, so arm B accuracies are relative, not absolute.
- Clusters are 40% identity; proteins at 30-40% identity can still be related, so leave-cluster-out is not leave-superfamily-out. MMseqs2 clustering and search depend on the version and parameters recorded in the run log.
- One model size, one pooling, single seed, CPU only. No simulated data in results; analysis_s2.py smoke test was on synthetic arms only, no record saved. Any crash fix is a dated AMENDMENT-N.md committed before outcomes exist.
