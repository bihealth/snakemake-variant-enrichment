from snakemake.shell import shell

shell(
    r"""
set -x
set -euo pipefail

name=$(basename {snakemake.output.vcf})

trap "rm -f {snakemake.output.vcf} {snakemake.output.tbi}" ERR

curl https://azureopendatastorage.blob.core.windows.net/gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{snakemake.wildcards.chrom}.vcf.bgz \
> {snakemake.output.vcf}

curl https://azureopendatastorage.blob.core.windows.net/gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.{snakemake.wildcards.chrom}.vcf.bgz.tbi \
> {snakemake.output.tbi}
""")
