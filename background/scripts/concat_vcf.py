from snakemake.shell import shell

shell(
    r"""
set -euo pipefail

name=$(basename {snakemake.output.vcf})

trap "rm -f {snakemake.output.vcf} {snakemake.output.tbi}" ERR

bcftools concat \
    -o {snakemake.output.vcf} \
    -O z \
    {snakemake.input.vcf}

tabix -f {snakemake.output.vcf}
"""
)
