"""Pre-filter VCF.

Only include records that have a LOW, MEDIUM, or HIGH effect and those that
have an overall gnomAD-exomes frequency below 1%.  The rationale is to save
subsequent running time.
"""

from snakemake.shell import shell


shell(
    r"""
set -x
set -euo pipefail

bcftools view \
    -i '((ANN ~ "MEDIUM") | (ANN ~ "HIGH") | (ANN ~ "LOW")) & ((GNOMAD_EXOMES_AF = ".") | (GNOMAD_EXOMES_AF <= 0.01))' \
    -O z \
    -o {snakemake.output.vcf} \
    {snakemake.input.vcf}

tabix -f {snakemake.output.vcf}
""")
