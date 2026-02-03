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

### Benchmarks: saige-rs vs SAIGE R

Benchmarked on the included test data (1,000 samples, 100 markers, covariates x1 + x2).

**Environment:**
- **saige-rs**: Native arm64 release binary on Apple Silicon (M-series Mac)
- **SAIGE R v1.5.1**: Docker container (linux/amd64 via Rosetta 2 emulation)
- Both single-threaded unless noted

#### Step 1: Fit Null GLMM

| Trait | SAIGE R | saige-rs | Speedup |
|---|---|---|---|
| Binary | 60.8s | 1.15s | **53x** |
| Quantitative | 64.3s | 0.34s | **189x** |
| Quantitative (4 threads) | -- | 0.16s | **402x** |

The R wall times include Docker startup (~2s), R package loading (~5s), genotype I/O (~30s for 10K-sample PLINK file), and GLMM fitting. saige-rs times include the full pipeline end-to-end.

Adjusting for Rosetta 2 overhead (~1.5-2x), saige-rs is conservatively **25-95x faster** than native R on the same hardware.

#### Step 2: Single-Variant Association Tests (70 variants tested)

| Format | SAIGE R | saige-rs | Speedup |
|---|---|---|---|
| PLINK | 20.2s | 0.025s | **808x** |
| BGEN | -- | 0.022s | -- |

The R Step 2 wall time includes Docker startup, R loading, genotype reading, and testing 70 variants. saige-rs times are end-to-end including I/O.

#### Notes

- R times include Rosetta 2 emulation overhead since the SAIGE R package requires a Linux x86_64 environment (building natively on macOS arm64 requires pixi and plink-ng). Native Linux R performance would be faster.
- saige-rs multi-threading (via rayon) shows near-linear scaling for GRM construction and trace estimation.
- Both implementations tested the same 70 variants and produced concordant results.
- Test dataset: `extdata/input/genotype_100markers` (1,000 phenotyped samples out of 10,000 in the genotype file).

