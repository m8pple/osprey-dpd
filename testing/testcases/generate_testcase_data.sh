#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Enable job control
set -m

dpd_executable=$1


if [[ $# -eq 1 ]] ; then
    for i in dmpci.* ; do 
        $0 $@ $i
    done
    exit 0
fi

dmpci_file=$2

function on_ctrl_c() {
    >&2 echo "Received ctrl-c. Killing $generate_pid and children. Will now summarise data, press ctrl-c again to cancel."
    kill %+
    trap - INT
}

../tools/run_replicates.py $dpd_executable $dmpci_file &
generate_pid=$!
trap on_ctrl_c INT

set +e
wait $generate_pid
set -e

>&2 echo "Extracting $dmpci_file"
../tools/extract_dmpcas_distributions.py $dmpci_file
>&2 echo "Generating graphs"
../tools/convert_dmpcas_distributions_to_graphs.py $dmpci_file
