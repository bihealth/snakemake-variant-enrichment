from snakemake.shell import shell

shell(
    r"""
set -x
set -euo pipefail

name=$(basename {snakemake.output.tsv})

trap "rm -f {snakemake.output.tsv} {snakemake.output.tbi}" ERR

curl https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz.tbi \
> {snakemake.output.tbi}

curl https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1/coverage/exomes/gnomad.exomes.coverage.summary.tsv.bgz \
> {snakemake.output.tsv}
""")
