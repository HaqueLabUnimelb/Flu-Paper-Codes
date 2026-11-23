Day 4 PR8 spatial randomness test

Samples: 2_2 and 2_3, tested separately.
PR8 expression: summed raw RNA UMI across detected PR8 features; thresholds use PR8_total_umi >= t.
Spatial graph: symmetrized 6-nearest-neighbour graph.
Null: PR8 expression values are exchangeable among spots within each of 10 rank-based total-RNA UMI strata.
Permutations: 1999; random seed: 20260803.
Continuous omnibus test: upper-tailed Moran's I on log1p(PR8 total UMI). Two-sample Bonferroni P values are reported.
Threshold test: upper-tailed positive-positive join count. Candidate thresholds must leave at least 100 positive and negative spots.
Threshold multiplicity: permutation maxT family-wise correction within each sample, then Bonferroni correction across the two samples.
Interpretation: the first significant threshold is an operational detection boundary, not a biological infection cutoff or a causal transition.
Fisher's exact test is not used because it ignores spatial adjacency; it remains suitable for a prespecified cell-type-by-positive 2x2 enrichment test.
