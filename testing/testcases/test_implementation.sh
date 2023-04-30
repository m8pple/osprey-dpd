#!/bin/bash
set -euo pipefail
IFS=$'\n\t'

# Enable job control
set -m

SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

function usage()
{
    >&2 echo "test_implementation.sh path_to_dpd dmpci_file? --working_dir=dir"
    >&2 echo ""
    >&2 echo "   path_to_dpd : the executable to test"
    >&2 echo "   dmpci_file? : a specific dmpci file to execute"
    >&2 echo "                 If no dmpci file is specified, test all in current directory."
    >&2 echo ""
    >&2 echo "   --working_dir : Where to put working files and outputs."
    >&2 echo "                   Will use 'mktemp -d' by default."
    exit 1
}


dpd_executable=""
dmpci_file=""
working_dir=""
replicate_count=4

while [[ $# -gt 0 ]]; do
  case $1 in
    --working_dir=*)
        working_dir="${1#*=}"
        shift
        ;;
    -*|--*)
        echo "Unknown option $1"
        exit 1
        ;;
    *)
        if [[ "$dpd_executable" == "" ]]; then
            >&2 echo "Set dpd_executable=$1"
            dpd_executable="$1"
            shift
        elif [[ "$dmpci_file" == "" ]]; then 
            >&2 echo "Set dmpci_file=$1"
            dmpci_file="$1"
            shift
        else
            echo "Unexpected parameter $1"
            exit 1
        fi
        ;;
  esac
done

if [[ "$working_dir" == "" ]] ; then
    working_dir=$(mktemp -d)
    >&2 echo "Temporay working dir $working_dir"
fi

>&2 echo "dpd=$dpd_executable, dmpci=$dmpci_file"

if [[ "$dmpci_file" == "" ]] ; then
    RES=0
    >&2 echo "Recurse"
    for i in dmpci.* ; do 
        >&2 echo "Recurse $i"
        echo "$0 --working_dir=$working_dir $dpd_executable $i" 
        set +e
        $0 --working_dir=$working_dir $dpd_executable $i 
        RES=$((RES+1))
        set -e
    done
    exit $RES
fi

dmpci_name=$(basename $dmpci_file)
dmpci_name=${dmpci_name/dmpci./}
dmpci_hash=$(${SCRIPTPATH}/../tools/calculate_dmpci_file_hash.py "$dmpci_file")

>&2 echo "Running ${replicate_count} replicates"
../tools/run_replicates.py $dpd_executable $dmpci_file ${replicate_count} $working_dir

>&2 echo "Extracting $dmpci_file"
../tools/extract_dmpcas_distributions.py $dmpci_file $working_dir $working_dir

>&2 echo "Generating graphs"
../tools/convert_dmpcas_distributions_to_graphs.py $working_dir/$dmpci_name-$dmpci_hash.h5

>&2 echo "Testing"
../tools/evaluate_series.py $working_dir $dmpci_file | tee $working_dir/${dmpci_name}-${dmpci_hash}.csv
