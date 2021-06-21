"""Strip gnomAD

Remove all annotations from gnomAD to speed up subsequent processing.
"""

from snakemake.shell import shell


shell(
    r"""
set -x
set -euo pipefail


zcat {snakemake.input.vcf} \
| awk -F $'\t' \
    'BEGIN {{ OFS=FS }}
     (/^#/) {{ print }}
     (!/^#/) {{ $8 = "."; print }}' \
| bgzip -c \
> {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}
""")
