#!/bin/bash
source `dirname $0`/../../env.sh

# parse experiment information
source `dirname $0`/parse_exp_info.sh

# variables
bin_dir="$WORK_DIR/bin"
lib_dir="$WORK_DIR/lib"

echo "[INFO] geneset analysis processing"
