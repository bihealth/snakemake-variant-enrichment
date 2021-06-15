# snakemake-variant-enrichment

A Snakemake workflow for germline variant burden analyses.

## Input

- **case**
    - a multi-VCF file with variants from one or multiple case groups
    - a TSV file mapping each sample to one group
- **control** -- you can make your pick either between
    - gnomAD population frequencies (hg37 and hg38 data sources available, exomes on hg38 limited to lift-over at the moment)
    - the IGSR data set of 2504 unrelated individuals (hg38 only)

## Manual

In the future we will provide pre-built background data.
However, for now you will have to build this yourself.

### Step 1: Build Background Data

```bash
cd background
mamba env create -f environment.yaml
conda activate snakemake-variant-enrichment
```
