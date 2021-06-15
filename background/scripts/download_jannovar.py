from snakemake.shell import shell

BASE_URL = "https://zenodo.org/record/4916051/files"

shell(
    r"""
set -x
set -euo pipefail

name=$(basename {snakemake.output.ser})

trap "rm -rf {snakemake.output.ser} {snakemake.output.ser_md5}" ERR

curl {BASE_URL}/$name?download=1 >{snakemake.output.ser}
curl {BASE_URL}/$name.md5?download=1 >{snakemake.output.ser_md5}

cd $(dirname {snakemake.output.ser})
md5sum --check $(basename {snakemake.output.ser_md5})
"""
)
