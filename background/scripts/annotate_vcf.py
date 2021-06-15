from snakemake.shell import shell

shell(
    r"""
set -x
set -euo pipefail

name=$(basename {snakemake.output.vcf})
db={snakemake.input.db}

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR {snakemake.output.vcf} {snakemake.output.tbi}" ERR

jannovar annotate-vcf \
    -d {snakemake.input.ser} \
    -i {snakemake.input.vcf} \
    -o $TMPDIR/tmp.vcf.gz
tabix -f $TMPDIR/tmp.vcf.gz

jannovar vardb-annotate \
    --genome-build {snakemake.params.genome_release} \
    --database-file ${{db%.h2.db}} \
    --input-vcf $TMPDIR/tmp.vcf.gz \
    --output-vcf {snakemake.output.vcf} \
    --table-names {snakemake.params.table_name}

tabix -f {snakemake.output.vcf}
""")
