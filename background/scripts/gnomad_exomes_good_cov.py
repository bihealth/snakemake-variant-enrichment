from snakemake.shell import shell

shell(
    r"""
set -x
set -euo pipefail

zcat {snakemake.input.tsv} \
| tail -n +2 \
| awk -F $'\t' 'BEGIN {{ OFS = FS; }} ($7 > 0.9) {{ print $1, $2 - 1, $2 }}' \
| bedtools merge -i stdin \
> {snakemake.output.bed}

bgzip -c {snakemake.output.bed} > {snakemake.output.bed_gz}
tabix -f -p bed {snakemake.output.bed_gz}
""")
