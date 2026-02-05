HERE we are maintaining an newly improved stable version of SAIGE and SAIGE-GENE+. 
Please find the https://saigegit.github.io/SAIGE-doc/ for documentation.


SAIGE is an R package developed with Rcpp for genome-wide association tests in large-scale data sets and biobanks. The method

- accounts for sample relatedness based on the generalized mixed models
- allows for model fitting with either full or sparse genetic relationship matrix (GRM)
- works for quantitative and binary traits
- handles case-control imbalance of binary traits
- computationally efficient for large data sets
- performs single-variant association tests
- provides effect size estimation through Firth's Bias-Reduced Logistic Regression
- performs conditional association analysis

SAIGE-GENE (now known as SAIGE-GENE+) are new method extension in the R package for testing rare variant in set-based tests.
- performs BURDEN, SKAT, and SKAT-O tests
- allows for tests on multiple minor allele frequencies cutoffs and functional annotations
- allows for specifying weights for markers in the set-based tests
- performs conditional analysis to identify associations independent from nearly GWAS signals


The package takes genotype file input in the following formats
- PLINK (bed, bim, fam), PGEN, BGEN, VCF, BCF, SAV


## saige-rs: Rust Implementation

`saige-rs` is a Rust rewrite of the SAIGE GWAS tool located in the `saige-rs/` directory. It reimplements the core statistical algorithms (GLMM fitting, SPA, score tests, variance ratio estimation) with native performance, memory-mapped I/O, and rayon parallelism.

### Building saige-rs

```bash
cd saige-rs
cargo build --release
# Binary: target/release/saige
```

### Usage

```bash
# Step 1: Fit null GLMM
saige fit-null \
  --plink-file <prefix> \
  --pheno-file <pheno.txt> \
  --pheno-col <trait> \
  --covar-cols <cov1,cov2> \
  --trait-type binary \
  --output-prefix <output>

# Step 2: Association tests
saige test \
  --model-file <output.saige.model> \
  --bgen-file <file.bgen> \
  --output-file <results.txt> \
  --is-spa --is-fast-spa
```

### Benchmarks: saige-rs vs R SAIGE

Simulated data: 2,000 samples, 5K GRM markers, 2,000 test variants, binary trait (~10% prevalence).

| Step | R SAIGE | saige-rs | Speedup |
|------|---------|----------|---------|
| Step 2: Association tests | 5.2s | 0.12s | **43x** |

R SAIGE run via Docker (linux/amd64). saige-rs native arm64 on Apple Silicon.

#### P-value concordance

Both implementations produce highly concordant p-values across 1,998 tested variants with SPA correction enabled.

| Metric | Value |
|--------|-------|
| P-value correlation (-log₁₀) | **0.995** |
| BETA correlation | **0.9997** |
| SE correlation | **0.9994** |
| Sign concordance | **99.4%** |

**Note:** R SAIGE must be run with `--is_fastTest=FALSE` for accurate variance estimation. The default `--is_fastTest=TRUE` uses an approximate formula that can produce inflated variance for high-AF variants.

<p align="center">
  <img src="saige-rs/docs/pvalue_comparison.png" width="700" alt="P-value comparison and runtime benchmarks between SAIGE R and saige-rs" />
</p>

