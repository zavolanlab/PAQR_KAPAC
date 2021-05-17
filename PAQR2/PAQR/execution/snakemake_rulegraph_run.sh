#!/usr/bin/env bash

# Create the Rule Graph for the workflow

cleanup () {
    rc=$?
    rm -rf ../images/rulegraph.svg
    cd "$user_dir"
    echo "Exit status: $rc"
}
trap cleanup SIGINT

set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
set -x  # facilitates debugging by printing out executed commands

user_dir=$PWD
pipeline_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd "$pipeline_dir"

snakemake \
    --snakefile="../Snakefile" \
    --configfile="../configs/config.yml" \
    --printshellcmds \
    --dryrun \
    --verbose \
    --rulegraph \
    | dot -Tsvg \
    > ../images/rulegraph.svg
