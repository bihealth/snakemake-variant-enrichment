from snakemake.shell import shell

BASE_URL = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

shell(
    r"""
set -x
set -euo pipefail

name=$(basename {snakemake.output.vcf})
export TMPDIR=$(mktemp -d)

trap "rm -rf $TMPDIR {snakemake.output.vcf} {snakemake.output.tbi}" ERR

curl {BASE_URL}/$name \
> $TMPDIR/$(basename {snakemake.output.vcf})

curl {BASE_URL}/$name.tbi \
> $TMPDIR/$(basename {snakemake.output.tbi})

bcftools view -e '(VT="SV")' \
    -o {snakemake.output.vcf} \
    -O z \
    $TMPDIR/$(basename {snakemake.output.vcf})

tabix -f {snakemake.output.vcf}
"""
)
