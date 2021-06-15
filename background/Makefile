.PHONY: default
default:

.PHONY: black
black:
	snakefmt -l 100 .
	black -l 100 .

.PHONY: snakemake-run
snakemake:
	snakemake -j1 -p --use-conda --conda-frontend=mamba

.PHONY: snakemake-create-envs-only
snakemake-create-envs-only:
	snakemake -j1 -p --use-conda --conda-create-envs-only
