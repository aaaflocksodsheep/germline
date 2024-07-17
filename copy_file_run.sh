#!/bin/bash
IFS='\n'
RUN_ID="240110_NDX550773_RUO_test_AH7LFNBGXV"
SAMPLE_LIST="/data/analysis/service/hd/240110_NDX550773_RUO_test_AH7LFNBGXV_samplelist.txt"

while read -r SAMPLE_NAME; do
  #echo ${SAMPLE_NAME}
  python  /data/script/lims/Hereditary_Disease/copy_file.py "$RUN_ID" "$SAMPLE_NAME" &
done < "$SAMPLE_LIST"

wait

