from snakemake.shell import shell

shell(
    r"""
set -x
set -euo pipefail

output={snakemake.output.db}

jannovar vardb-import \
    --genome-build "{snakemake.params.genome_release}" \
    --database-file "${{output%.h2.hdb}}" \
    --vcf-files "$(echo {snakemake.input.vcf} | tr ' ' ',')" \
    --table-name "{snakemake.params.table_name}" \
    --db-name "gnomAD-exomes" \
    --db-version "2.2.1" \
    --default-prefix "GNOMAD_EXOMES_" \
    --vcf-info-fields "AF,AC,AN" \
    --truncate-table
"""
)
